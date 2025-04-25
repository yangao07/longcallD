#include <ctype.h>
#include "bam_utils.h"
#include "utils.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "sdust.h"

extern int LONGCALLD_VERBOSE;

char test_read_name[1024];
char test_chr[1024];

read_var_profile_t *init_read_var_profile(int n_reads, int n_total_vars) {
    read_var_profile_t *p = (read_var_profile_t*)malloc(n_reads * sizeof(read_var_profile_t));
    for (int i = 0; i < n_reads; ++i) {
        p[i].read_id = i;
        p[i].start_var_idx = -1; p[i].end_var_idx = -2;
        p[i].alleles = (int*)malloc(n_total_vars * sizeof(int));
        for (int j = 0; j < n_total_vars; ++j) p[i].alleles[j] = -1; // init as unused
    }
    return p;
}

read_var_profile_t *init_read_var_profile_with_ids(int n_reads, int *read_ids, int n_total_vars) {
    read_var_profile_t *p = (read_var_profile_t*)malloc(n_reads * sizeof(read_var_profile_t));
    for (int i = 0; i < n_reads; ++i) {
        p[i].read_id = read_ids[i];
        p[i].start_var_idx = -1; p[i].end_var_idx = -2;
        p[i].alleles = (int*)malloc(n_total_vars * sizeof(int));
        for (int j = 0; j < n_total_vars; ++j) p[i].alleles[j] = -1; // init as unused
    }
    return p;
}

void free_read_var_profile(read_var_profile_t *p, int n_reads) {
    for (int i = 0; i < n_reads; ++i) {
        if (p[i].alleles) free(p[i].alleles);
    }
    free(p);
}

int has_equal_X_in_bam_cigar(bam1_t *read) {
    // Get the CIGAR string for the read
    const uint32_t *cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    if (n_cigar <= 0) return 0;
    // Check if '=' or 'X' is in the CIGAR operations
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            return 1;
        } else if (op == BAM_CMATCH) {
            return 0;
        }
    }
    return 0;
}

int has_cs_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "cs");
    return s != NULL;
}

int has_MD_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "MD");
    return s != NULL;
}

int get_aux_int_from_bam(bam1_t *b, const char *tag) {
    uint8_t *s = bam_aux_get(b, tag);
    if (s == NULL) return -1;
    return bam_aux2i(s);
}

char *get_aux_str_from_bam(bam1_t *b, const char *tag) {
    uint8_t *s = bam_aux_get(b, tag);
    if (s == NULL) return NULL;
    return (char*)s;
}

void print_digar(digar_t *digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tis_low_qual\n");
    for (int i = 0; i < digar->n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d\n", digar->digars[i].pos, BAM_CIGAR_STR[digar->digars[i].type], digar->digars[i].len, digar->digars[i].qi, digar->digars[i].is_low_qual);
    }
}

typedef struct {
    hts_pos_t *pos; int *lens;
    int *counts, *is_dense;
    int front, rear, count;
    int max_s, win; // max_s: max sub/gaps in win
    int size_m;
} xid_queue_t;

xid_queue_t *init_xid_queue(int max_sites, int max_s, int win) {
    xid_queue_t *q = (xid_queue_t*)malloc(sizeof(xid_queue_t));
    q->pos = (hts_pos_t*)malloc(max_sites * sizeof(hts_pos_t));
    q->counts = (int*)malloc(max_sites * sizeof(int));
    q->lens = (int*)malloc(max_sites * sizeof(int));
    q->is_dense = (int*)calloc(max_sites, sizeof(int));
    q->size_m = max_sites;
    q->front = 0; q->rear = -1; q->count = 0;
    q->max_s = max_s; q->win = win;
    return q;
}

void realloc_xid_queue(xid_queue_t *q) {
    if (q->rear+1 >= q->size_m) {
        fprintf(stderr, "realloc xid_queue: %s %s\n", test_chr, test_read_name);
        q->size_m *= 2;
        q->pos = (hts_pos_t*)realloc(q->pos, q->size_m * sizeof(hts_pos_t));
        q->counts = (int*)realloc(q->counts, q->size_m * sizeof(int));
        q->lens = (int*)realloc(q->lens, q->size_m * sizeof(int));
        q->is_dense = (int*)realloc(q->is_dense, q->size_m * sizeof(int));
        for (int i = q->rear+1; i < q->size_m; ++i) q->is_dense[i] = 0;
    }
}

void free_xid_queue(xid_queue_t *q) {
    free(q->is_dense);
    free(q->pos); free(q->counts); 
    free(q->lens);
    free(q);
}

void push_xid_size_queue_win(xid_queue_t *q, hts_pos_t pos, int len, int count,
                             cgranges_t *cr, hts_pos_t *cr_cur_start, hts_pos_t *cr_cur_end,
                             int *cr_q_start, int *cr_q_end) {
    realloc_xid_queue(q);
    q->pos[++q->rear] = pos; // left-most pos of the region
    q->lens[q->rear] = len; q->counts[q->rear] = count;
    q->count += count;

    hts_pos_t noisy_start = -1, noisy_end = -1;
    while (q->pos[q->front]+q->lens[q->front]-1 <= pos - q->win) {
        q->count -= q->counts[q->front]; // q->lens[q->front];
        q->front++;
    }
    // fprintf(stderr, "pos: %" PRIi64 ", len: %d, count: %d, front: %d, rear: %d, count: %d\n", pos, len, count, q->front, q->rear, q->count);
    if (count > 0) {
        if (q->count > q->max_s) {
            for (int j = q->front; j <= q->rear; j++)
                q->is_dense[j] = 1; // mark all dense regions
            noisy_start = q->pos[q->front]; noisy_end = q->pos[q->rear]+q->lens[q->rear];
            if (*cr_cur_start == -1) {
                *cr_cur_start = noisy_start; *cr_cur_end = noisy_end;
                *cr_q_start = q->front; *cr_q_end = q->rear;
            } else {
                if (noisy_start <= *cr_cur_end) { // merge
                    *cr_cur_end = noisy_end; // *cr_qi stays the same
                    *cr_q_end = q->rear; // update end
                } else {
                    // add the previous region
                    int var_size = 0;
                    for (int i = *cr_q_start; i <= *cr_q_end; i++) var_size += q->counts[i];
                    var_size = var_size > (*cr_cur_end - *cr_cur_start + 1) ? var_size : (*cr_cur_end - *cr_cur_start + 1);
                    cr_add(cr, "cr", *cr_cur_start-1, *cr_cur_end, var_size);
                    *cr_cur_start = noisy_start; *cr_cur_end = noisy_end;
                    *cr_q_start = q->front; *cr_q_end = q->rear;
                }
            }
            // fprintf(stderr, "noisy region: %" PRId64 "-%" PRId64 "\n", noisy_start, noisy_end);
        }
    }
}

int get_var_start(cand_var_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start) {
    int i;
    for (i = cur_site_i; i < n_total_pos; ++i) {
        if (var_sites[i].pos >= start) return i;
    }
    return i;
}

int get_var_site_start(var_site_t *var_sites, int cur_site_i, int n_total_pos, hts_pos_t start) {
    int i;
    for (i = cur_site_i; i < n_total_pos; ++i) {
        if (var_sites[i].pos >= start) return i;
    }
    return i;
}

// only for mismatch and insertion
uint8_t *make_alt_seq(digar1_t *digar, const uint8_t *bseq) {
    if (digar->type == BAM_CDEL) return NULL;
    uint8_t *alt_seq = (uint8_t*)malloc((digar->len) * sizeof(uint8_t));
    for (int i = 0; i < digar->len; ++i) {
        alt_seq[i] = seq_nt16_int[bam_seqi(bseq, digar->qi+i)];
    }
    return alt_seq;
}

void update_var_site_with_allele(cand_var_t *cand_var, int is_low_qual, uint8_t strand, int allele_i) {
    if (is_low_qual) { // low qual base quality, allele_i considered as -1
        cand_var->low_qual_cov++;
        return;
    }
    cand_var->total_cov++;
    cand_var->alle_covs[allele_i] += 1; // ref/alt
    cand_var->strand_to_alle_covs[strand][allele_i] += 1; // strand-wise: 0:forward/1:reverse -> allele_i -> read count
}

void update_read_var_profile_with_allele(int var_i, int allele_i, read_var_profile_t *read_var_profile) {
    // even allele_i is -1, we still set it XXX
    if (read_var_profile->start_var_idx == -1) read_var_profile->start_var_idx = var_i;
    read_var_profile->end_var_idx = var_i;
    int _var_i = var_i - read_var_profile->start_var_idx;
    read_var_profile->alleles[_var_i] = allele_i; // ref:0, alt:1~n, or -1: not ref/alt
}


int get_digar_ave_qual(digar1_t *digar, const uint8_t *qual) {
    if (digar->qi < 0) return 0;
    int q_start, q_end;
    if (digar->type == BAM_CDEL) {
        if (digar->qi == 0) {
            q_start = 0;
            q_end = 0;
        } else {
            q_start = digar->qi - 1;
            q_end = digar->qi;
        }
    } else {
        q_start = digar->qi;
        q_end = digar->qi + digar->len - 1;
    }
    int ave_qual = 0, qlen = q_end - q_start + 1;
    for (int i = q_start; i <= q_end; ++i) {
        ave_qual += qual[i];
    }
    return ave_qual / qlen;
}

// cand_vars:
// X: 1 -> M SNPs
// I: 0 -> M INSs
// D: M -> 1 DEL
// XXX no minor_alt_allele here, will be used in read_var_profile
// here, minor_alt_allele will be simply treated as non-alt, i.e., ref allele
int update_cand_vars_from_digar(const call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_var_sites, var_site_t *var_sites, int start_i, cand_var_t *cand_vars) {
    int cur_start_i, site_i, digar_i = 0;
    // hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    hts_pos_t pos_start = digar->beg, pos_end = digar->end;
    digar1_t *digars = digar->digars; uint8_t strand = digar->is_rev;

    cur_start_i = site_i = get_var_site_start(var_sites, start_i, n_var_sites, pos_start);
    int tid = chunk->tid; // read->core.tid;
    // compare sorted var_sites and digars
    // both var_sites and digars are sorted by pos/type/ref_len
    for (; site_i < n_var_sites && digar_i < digar->n_digar; ) {
        if (digars[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar->digars+digar_i);
        int ave_qual = get_digar_ave_qual(digar->digars+digar_i, digar->qual);
        int ret = comp_var_site(var_sites+site_i, &digar_var_site);
        if (ret < 0) { // var_site < digar_var_site
            // update_var_site_with_allele(cand_vars+site_i, is_in_noisy_reg(cand_vars[site_i].pos, digar->noisy_regs), strand, 0);
            update_var_site_with_allele(cand_vars+site_i, 0, strand, 0);
            site_i++;
        } else if (ret == 0) { // merge together, update/add alt allele
            update_var_site_with_allele(cand_vars+site_i, (digar->digars[digar_i].is_low_qual) || (ave_qual < opt->min_bq), strand, 1);
            site_i++;
        } else { // ret > 0, var_site > digar_var_site
            digar_i++;
        }
    }
    for (; site_i < n_var_sites; ++site_i) {
        if (var_sites[site_i].pos > pos_end) break;
        // update_var_site_with_allele(cand_vars+site_i, is_in_noisy_reg(cand_vars[site_i].pos, digar->noisy_regs), strand, 0);
        update_var_site_with_allele(cand_vars+site_i, 0, strand, 0);
    }
    return cur_start_i;
}

int get_alt_allele_idx(cand_var_t *cand_var, digar1_t *digar, const uint8_t *bseq) {
    if (digar->is_low_qual) return -1; // low qual base quality
    else return 1; // alt allele
    // if (cand_var->var_type == BAM_CDEL) return 1; // only 0/1 for deletion
    // uint8_t *alt_seq = make_alt_seq(digar, bseq);
    // int alt_len = digar->len;
    // int allele_i = -1;
    // for (int i = 0; i < cand_var->n_uniq_alles-1; ++i) {
    //     if (cand_var->alt_lens[i] == alt_len && memcmp(cand_var->alt_seqs[i], alt_seq, alt_len) == 0) {
    //         allele_i = i+1;
    //         break;
    //     }
    // }
    // free(alt_seq);
    // return allele_i;
    // if (allele_i == -1) {
    //     // _err_error_exit("Alt allele not found: %d\n", allele_i);
    //     allele_i = cand_var->n_uniq_alles;
    // }
    // return allele_i;
}

int update_read_var_profile_from_digar(const call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, cand_var_t *cand_vars, int start_var_i, read_var_profile_t *read_var_profile) {
    int cur_start_i, var_i, digar_i = 0, allele_i=-1;
    // hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    hts_pos_t pos_start = digar->beg, pos_end = digar->end;
    digar1_t *digar1 = digar->digars;
    cur_start_i = var_i = get_var_start(cand_vars, start_var_i, n_cand_vars, pos_start);
    int tid = chunk->tid; //read->core.tid;
    // compare sorted cand_vars and digars
    // both cand_vars and digars are sorted by pos/type/ref_len
    for (; var_i < n_cand_vars && digar_i < digar->n_digar; ) {
        if (digar1[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t var_site0 = make_var_site_from_cand_var(cand_vars+var_i);
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar1+digar_i);
        int ave_qual = get_digar_ave_qual(digar1+digar_i, digar->qual);
        // if (ave_qual < opt->min_bq) {
            // digar_i++; continue;
        // }
        int is_ovlp, ret;
        ret = comp_ovlp_var_site(&var_site0, &digar_var_site, &is_ovlp);
        if (is_ovlp == 0) { // no overlapping
            if (ret < 0) { // var_site < digar_var_site: ref_allele
                update_read_var_profile_with_allele(var_i, 0, read_var_profile); //is_in_noisy_reg(var_site0.pos, digar->noisy_regs));
                var_i++;
            } else if (ret > 0) { // var_site > digar_var_site
                digar_i++;
            } else { // var_site == digar_var_site
                _err_error("Unexpected case: is_ovlp == 0 && var_site == digar_var_site, %d-%c-%d-%d\n", cand_vars[var_i].pos, BAM_CIGAR_STR[cand_vars[var_i].var_type], digar_var_site.pos, BAM_CIGAR_STR[digar_var_site.var_type]);
                var_i++; digar_i++;
            }
        } else { // overlap
            if (ret == 0) { // exact the same with var_site
                allele_i = get_alt_allele_idx(cand_vars+var_i, digar->digars+digar_i, digar->bseq);
                if (ave_qual < opt->min_bq) allele_i = -1;
                update_read_var_profile_with_allele(var_i, allele_i, read_var_profile); // is_in_noisy_reg(var_site0.pos, digar->noisy_regs));
                var_i++;
            } else {
                update_read_var_profile_with_allele(var_i, -1, read_var_profile); // is_in_noisy_reg(var_site0.pos, digar->noisy_regs));
                var_i++;
            }
        }
    }
    for (; var_i < n_cand_vars; ++var_i) {
        if (cand_vars[var_i].pos > pos_end) break;
        if (is_in_noisy_reg(cand_vars[var_i].pos, digar->noisy_regs)) continue;
        update_read_var_profile_with_allele(var_i, 0, read_var_profile); // is_in_noisy_reg(cand_vars[var_i].pos, digar->noisy_regs));
    }
    return cur_start_i;
}

void set_digar(digar1_t *digar, hts_pos_t pos, int type, int len, int qi, int is_low_qual, uint8_t *alt_seq) {
    digar->pos = pos; digar->type = type; digar->len = len; digar->qi = qi; digar->alt_seq = alt_seq; digar->is_low_qual = is_low_qual;
}

void push_digar1(digar_t *digar, digar1_t d) {
    _uni_realloc(digar->digars, digar->n_digar, digar->m_digar, digar1_t);
    digar->digars[digar->n_digar].pos = d.pos; digar->digars[digar->n_digar].type = d.type;
    digar->digars[digar->n_digar].len = d.len; digar->digars[digar->n_digar].alt_seq = d.alt_seq;
    digar->digars[digar->n_digar].qi = d.qi; digar->digars[digar->n_digar].is_low_qual = d.is_low_qual;
    digar->n_digar++;
}

void print_digar1(digar1_t *digars, int n_digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d\n", digars[i].pos, BAM_CIGAR_STR[digars[i].type], digars[i].len, digars[i].qi, digars[i].is_low_qual);
    }
}

int collect_noisy_region_len(cgranges_t *noisy_reg) {
    int len = 0;
    for (int i = 0; i < noisy_reg->n_r; ++i) {
        hts_pos_t noisy_reg_beg = cr_start(noisy_reg, i), noisy_reg_end = cr_end(noisy_reg, i);
        len += (noisy_reg_end - noisy_reg_beg + 1);
    }
    return len;
}

// if [sa_pos, sa_end] overlaps with [primary_pos, primary_end] by over 90% of the length, return 1, else return 0
int check_ont_palindrome(hts_pos_t primary_pos, hts_pos_t primary_end, hts_pos_t sa_pos, hts_pos_t sa_end) {
    hts_pos_t primary_len = primary_end - primary_pos + 1, sa_len = sa_end - sa_pos + 1;
    hts_pos_t overlap_len = 0;
    if (sa_pos <= primary_pos) {
        if (sa_end >= primary_end) overlap_len = primary_len;
        else if (sa_end >= primary_pos) overlap_len = sa_end - primary_pos + 1;
    } else if (sa_pos <= primary_end) {
        if (sa_end >= primary_end) overlap_len = primary_end - sa_pos + 1;
        else overlap_len = sa_len;
    }
    if (overlap_len >= primary_len * 0.9) return 1;
    return 0;
}

// return 
//   1: if clipping is potentially noisy
//   0: if clipping is supplementary alignment (palindrome for ont)
int check_noisy_end_clip(const call_var_opt_t *opt, bam1_t *read) {
    if (opt->is_ont)  { // check SA tag
        hts_pos_t primary_pos = read->core.pos+1, primary_end = bam_endpos(read);
        uint8_t *sa_tag = bam_aux_get(read, "SA");  // Get SA tag
        if (sa_tag) { // collect mapping chrom, pos, strand, end pos, MAPQ, NM
            char *sa_value = bam_aux2Z(sa_tag);
            char *sa_entry = strtok(sa_value, ";");
            while (sa_entry) {
                char rname[100], cigar[100], strand;
                int pos, mapq, flag;
                
                sscanf(sa_entry, "%[^,],%d,%c,%[^,],%d,%d", rname, &pos, &strand, cigar, &mapq, &flag);

                // Estimate end position from CIGAR
                int32_t sa_end_pos = pos;
                char *cptr = cigar;
                while (*cptr) {
                    int len = strtol(cptr, &cptr, 10);
                    if (*cptr == 'M' || *cptr == 'D' || *cptr == '=' || *cptr == 'X') {
                        sa_end_pos += len;
                    }
                    cptr++;
                }
                sa_entry = strtok(NULL, ";");
                // if any SA entry is a palindrome, then it is not a noisy region
                if (check_ont_palindrome(primary_pos, primary_end, pos, sa_end_pos)) {
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "palindrome: %s %" PRIi64 "-%" PRIi64 ", sa: %s %d-%d\n", bam_get_qname(read), primary_pos, primary_end, rname, pos, sa_end_pos);
                    return 0;
                }
            }
        }
    }
    return 1;
}

// is_low_qual: low base quality
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;
    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    uint32_t qlen = read->core.l_qseq;
    digar->bseq = (uint8_t*)malloc((qlen+1)/2 * sizeof(uint8_t));
    for (int i = 0; i < (qlen+1)/2; ++i) digar->bseq[i] = bam_get_seq(read)[i];
    digar->qual = (uint8_t*)malloc(qlen * sizeof(uint8_t));
    for (int i = 0; i < qlen; ++i) digar->qual[i] = bam_get_qual(read)[i];
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win); // for noisy region
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                if (digar->qual[qi] >= opt->min_bq) { // skip low-quality bases for noisy regions
                    push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                    set_digar(_digars+_n_digar, pos, op, 1, qi, 0, x_seq);
                } else set_digar(_digars+_n_digar, pos, op, 1, qi, 1, x_seq);
                _n_digar++; // push_xid_queue(q, pos, 1, 1);
                n_total_cand_vars++;
                pos++; qi++;
            }
        } else if (op == BAM_CEQUAL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            // push_xid_queue_win(q, pos, len, 0, digar->noisy_regs, &noisy_start, &noisy_end);
            pos += len; qi += len;
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, op, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, op, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL); // clipping
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (check_noisy_end_clip(opt, read)) {
                    if (i == 0) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                    } else {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                    }
                    n_total_cand_vars++;
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) { // XXX no action for N op
            pos += len;
        } else if (op == BAM_CMATCH) {
            _err_error_exit("CIGAR operation 'M' is not expected in EQX CIGAR: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) { // XXX not set noisy-region for digars
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    // if total noisy region length > mapped length * 0.8, skip the read
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "SkipRead: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) { // we may double-count reads with multiple nearby noisy regions 
                                                        // downside is minor, as we just consider an easy region as noisy, variant calling result should not be affected
                                                        // read1:  ---- [noisy] ---- [noisy] ----
                                                        // read2:  --------- [  noisy  ] --------
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end)) {
                // set cr_label as the largest insertion in noisy region
                // the insertion size is used to determine the merge_distance:
                // ins_size: merge_dis
                //   <= 100: 100
                //    < 500: ins_size
                //   >= 500: 500
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            }
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "DIGAR1: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

// XXX TODO make it consistent with collect_digar_from_eqx_cigar
int collect_digar_from_cs_tag(bam_chunk_t *chunk, bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    uint8_t *cs_tag = bam_aux_get(read, "cs");
    if (cs_tag == NULL) { _err_error_exit("cs tag not found in the BAM file: %s", bam_get_qname(read)); }
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    if (n_cigar <= 0) { _err_error_exit("CIGAR not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    digar->bseq = (uint8_t*)malloc((qlen+1)/2 * sizeof(uint8_t));
    for (int i = 0; i < (qlen+1)/2; ++i) digar->bseq[i] = bam_get_seq(read)[i];
    digar->qual = (uint8_t*)malloc(qlen * sizeof(uint8_t));
    for (int i = 0; i < qlen; ++i) digar->qual[i] = bam_get_qual(read)[i];
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *cs = bam_aux2Z(cs_tag); int cs_i = 0, cs_len = strlen(cs);
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    // left-end clipping
    if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[0]) == BAM_CHARD_CLIP) {
        int len = bam_cigar_oplen(cigar[0]);
        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
        set_digar(_digars+_n_digar, pos, bam_cigar_op(cigar[0]), len, qi, 0, NULL); // clipping
        _n_digar++; // push_xid_queue(q, pos, 0, 0);
        if (len > end_clip_reg) {
            if (check_noisy_end_clip(opt, read)) {
                if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                // if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                n_total_cand_vars++;
            }
        }
        if (bam_cigar_op(cigar[0]) == BAM_CSOFT_CLIP) qi += len;
    }
    while (*cs) {
        if (*cs == ':') { // identical seq length
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            int len = strtol(&cs[1], &cs, 10);
            set_digar(_digars+_n_digar, pos, BAM_CEQUAL, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            pos += len; qi += len;
        } else if (*cs == '=') { // identical sequence
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            set_digar(_digars+_n_digar, pos, BAM_CEQUAL, len, qi, 0, NULL);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            pos += len; qi += len;
        } else if (*cs == '*') { // mismatch: 1bp
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
            x_seq[0] = nst_nt4_table[(int)*(cs+2)];
            if (digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
            } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
            _n_digar++; pos++; qi++; cs += 3;
            n_total_cand_vars++;
        } else if (*cs == '+') { // insertion
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int i = 0; i < len; ++i) ins_seq[i] = nst_nt4_table[(int)*(cs-len+i)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi+_i] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq);
            _n_digar++; qi += len;
            n_total_cand_vars++;
        } else if (*cs == '-') { // deletion
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            cs++; int len = 0;
            while (isalpha(*cs)) {
                len++; cs++;
            }
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; pos += len;
            n_total_cand_vars++;
        } else if (*cs == '~') { // intron/splice
            cs++;
            while (isalpha(*cs) || isdigit(*cs)) cs++;
        } else {
            _err_error_exit("Unexpected character in cs tag: %c", *cs);
        }
    }
    
    // right-end clipping
    if (bam_cigar_op(cigar[n_cigar-1]) == BAM_CSOFT_CLIP || bam_cigar_op(cigar[n_cigar-1]) == BAM_CHARD_CLIP) {
        int len = bam_cigar_oplen(cigar[n_cigar-1]);
        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
        set_digar(_digars+_n_digar, pos, bam_cigar_op(cigar[n_cigar-1]), len, qi, 0, NULL); // clipping
        _n_digar++; // push_xid_queue(q, pos, 0, 0);
        if (len > end_clip_reg) {
            if (check_noisy_end_clip(opt, read)) {
                if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                n_total_cand_vars++;
            }
        }
        if (bam_cigar_op(cigar[n_cigar-1]) == BAM_CSOFT_CLIP) qi += len;
    }

    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "Skip: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "DIGAR2: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

int collect_digar_from_MD_tag(bam_chunk_t *chunk, bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    uint8_t *s = bam_aux_get(read, "MD");
    if (s == NULL) { _err_error_exit("MD tag not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    digar->bseq = (uint8_t*)malloc((qlen+1)/2 * sizeof(uint8_t));
    for (int i = 0; i < (qlen+1)/2; ++i) digar->bseq[i] = bam_get_seq(read)[i];
    digar->qual = (uint8_t*)malloc(qlen * sizeof(uint8_t));
    for (int i = 0; i < qlen; ++i) digar->qual[i] = bam_get_qual(read)[i];
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *md = bam_aux2Z(s); int md_i = 0;
    // printf("MD: %s\n", md);
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    int last_eq_len = 0;
    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        // fprintf(stderr, "%c %d, md_i: %d, MD: %s\n", BAM_CIGAR_STR[op], len, md_i, md);
        if (op == BAM_CMATCH) {
            // check if base is equal or mismatch based on MD tag
            int m_len = len, eq_len;
            while (1) { // MD
                if (last_eq_len > 0) {
                    if (last_eq_len >= m_len) {
                        eq_len = m_len;
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos, 0, 0);
                        pos += eq_len; qi += eq_len;
                        last_eq_len -= m_len; m_len = 0;
                    } else { // last_eq_len < m_len
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, last_eq_len, qi, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos, 0, 0);
                        pos += last_eq_len; qi += last_eq_len;
                        m_len -= last_eq_len; md_i = 0; // start of the next run
                        last_eq_len = 0;
                    }
                } else if (isdigit(md[md_i])) { // =
                    eq_len = strtol(&md[md_i], &md, 10);
                    if (eq_len > m_len) { //_err_error_exit("MD and CIGAR do not match: %s (%d : %d)", bam_get_qname(read), eq_len, m_len);
                        last_eq_len = eq_len - m_len;
                        eq_len = m_len;
                    } else if (eq_len == 0) {
                        md_i=0; continue;
                    }
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0, NULL); // qual unset XXX
                    _n_digar++; // push_xid_queue(q, pos, 0, 0);
                    pos += eq_len; qi += eq_len;
                    m_len -= eq_len; md_i = 0; // start of the next run
                } else if (isalpha(md[md_i])) { // X
                    // if (qual[qi] >= min_bq) {
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                    x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                    if (digar->qual[qi] >= opt->min_bq) {
                        push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
                    _n_digar++; // push_xid_queue(q, pos, 1, 1);
                    n_total_cand_vars++;
                    pos++; qi++; m_len -= 1; 
                    if (md[md_i+1] == '\0' || md[md_i+1] != '0') md_i++;
                    else md_i += 2; // skip the 0 after X
                } else _err_error_exit("MD and CIGAR do not match: %s (%s)", bam_get_qname(read), bam_aux2Z(s));
                if (m_len <= 0) break;
            }
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
            // MD
            md_i++;
            while (md[md_i] && isalpha(md[md_i])) md_i++;
            if (md[md_i] == '0') md_i++; // skip 0 after D
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0, NULL); // clipping
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (check_noisy_end_clip(opt, read)) {
                    if (i == 0) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                    } else {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                    }
                    n_total_cand_vars++;
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        // fprintf(stderr, "Skip: %s %d-noisy %d-var %d-ref_span\n", bam_get_qname(read), total_noisy_reg_len, n_total_cand_vars, mapped_len);
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "DIGAR2: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

// XXX TODO make it consistent with collect_digar_from_eqx_cigar
int collect_digar_from_ref_seq(bam_chunk_t *chunk, bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read); digar->is_rev = bam_is_rev(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    uint32_t qlen = read->core.l_qseq;
    digar->bseq = (uint8_t*)malloc((qlen+1)/2 * sizeof(uint8_t));
    for (int i = 0; i < (qlen+1)/2; ++i) digar->bseq[i] = bam_get_seq(read)[i];
    digar->qual = (uint8_t*)malloc(qlen);
    for (int i = 0; i < qlen; ++i) digar->qual[i] = bam_get_qual(read)[i];
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = chunk->whole_ref_len;
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    strcpy(test_read_name, bam_get_qname(read));
    strcpy(test_chr, chunk->tname);

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CDIFF || op == BAM_CEQUAL) {
            int eq_len = 0;
            for (int j = 0; j < len; ++j) {
                // Get the reference base
                if (pos <= ref_beg || pos > ref_end) {
                    if (LONGCALLD_VERBOSE >= 2) {
                        fprintf(stderr, "pos: %" PRIi64 " (%" PRIi64 "-%" PRIi64 ")\t", pos, ref_beg, ref_end);
                        fprintf(stderr, "Read exceed reference region sequence: %s", bam_get_qname(read));
                    }
                    pos++; qi++; continue;
                }
                int ref_base = nst_nt4_table[(int)ref_seq[pos-ref_beg]];
                int read_base = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                if (ref_base != read_base) {
                    if (eq_len > 0) {
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0, NULL);
                        _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                        eq_len = 0;
                    }
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                    x_seq[0] = read_base; //seq_nt16_int[bam_seqi(digar->bseq, qi)];
                    if (digar->qual[qi] >= opt->min_bq) {
                        push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
                    _n_digar++;
                    n_total_cand_vars++;
                } else eq_len++;
                pos++; qi++;
            }
            if (eq_len > 0) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0, NULL);
                _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                eq_len = 0;
            }
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0, NULL);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1, NULL);
            _n_digar++; // push_xid_queue(q, pos, len, 1);
            n_total_cand_vars++;
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            uint8_t *ins_seq = (uint8_t*)malloc(len * sizeof(uint8_t));
            for (int j = 0; j < len; ++j) ins_seq[j] = seq_nt16_int[bam_seqi(digar->bseq, qi+j)];
            // set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi); // insertion
            int is_low_qual = 1;
            for (int _i = 0; _i < len; ++_i) {
                if (digar->qual[qi] >= opt->min_bq) {
                    is_low_qual = 0; break;
                }
            }
            if (!is_low_qual) push_xid_size_queue_win(q, pos, 0, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq);
            _n_digar++; //push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CSOFT_CLIP, len, qi, 0, NULL); // clipping 
            _n_digar++; //push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (check_noisy_end_clip(opt, read)) {
                    if (i == 0) {
                        if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                    } else {
                        if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                    }
                    n_total_cand_vars++;
                }
            }
            if (op == BAM_CSOFT_CLIP) qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } // else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            // _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        // }
    }
    for (int i = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
        // fprintf(stderr, "var_size: %d, noisy_start: %" PRIi64 ", noisy_end: %" PRIi64 "\n", var_size, noisy_start, noisy_end);
        cr_add(digar->noisy_regs, "cr", noisy_start-1, noisy_end, var_size);
    }
    cr_index(digar->noisy_regs);
    int skip = 0;
    int total_noisy_reg_len = collect_noisy_region_len(digar->noisy_regs);
    int mapped_len = digar->end - digar->beg + 1;
    if (total_noisy_reg_len > mapped_len * max_noisy_frac_per_read || n_total_cand_vars > mapped_len * max_var_ratio_per_read) {
        skip = 1;
    } else {
        for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "NoisyRead %s %s:%d-%d %d\n", bam_get_qname(read), chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
            if (is_overlap_reg(cr_start(digar->noisy_regs, i)+1, cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "DIGAR3: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
            if (digar->noisy_regs->n_r > 0) fprintf(stderr, "Noisy-%s\n", bam_get_qname(read));
            for (int i = 0; i < digar->noisy_regs->n_r; ++i) {
                fprintf(stderr, "Noisy: %s:%d-%d\n", chunk->tname, cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i));
            }
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

int bam_chunk_init0(bam_chunk_t *chunk, int n_reads) {
    // input
    chunk->n_reads = 0; chunk->m_reads = n_reads;
    chunk->n_up_ovlp_reads = 0; chunk->up_ovlp_read_i = (int*)malloc(n_reads * sizeof(int));
    chunk->n_down_ovlp_reads = 0; chunk->down_ovlp_read_i = (int*)malloc(n_reads * sizeof(int));
    chunk->reads = (bam1_t**)malloc(n_reads * sizeof(bam1_t*));
    chunk->ref_seq = NULL;
    chunk->low_comp_cr = NULL;
    // intermediate
    chunk->is_skipped = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->digars = (digar_t*)calloc(n_reads, sizeof(digar_t));
    for (int i = 0; i < n_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
    }
    // noisy regions
    chunk->chunk_noisy_regs = NULL; chunk->noisy_reg_to_reads = NULL; chunk->noisy_reg_to_n_reads = NULL;
    // variant
    chunk->n_cand_vars = 0; chunk->cand_vars = NULL; chunk->var_i_to_cate = NULL;
    chunk->read_var_profile = NULL; chunk->read_var_cr = NULL;
    // output:
    chunk->flip_hap = 0;
    chunk->haps = (int*)calloc(n_reads, sizeof(int)); // read-wise
    chunk->phase_sets = (hts_pos_t*)calloc(n_reads, sizeof(hts_pos_t)); // read-wise
    return 0;
}

int bam_chunk_realloc(bam_chunk_t *chunk) {
    int m_reads = chunk->m_reads * 2;
    chunk->reads = (bam1_t**)realloc(chunk->reads, m_reads * sizeof(bam1_t*));
    chunk->up_ovlp_read_i = (int*)realloc(chunk->up_ovlp_read_i, m_reads * sizeof(int));
    chunk->down_ovlp_read_i = (int*)realloc(chunk->down_ovlp_read_i, m_reads * sizeof(int));
    chunk->is_skipped = (uint8_t*)realloc(chunk->is_skipped, m_reads * sizeof(uint8_t));
    chunk->digars = (digar_t*)realloc(chunk->digars, m_reads * sizeof(digar_t));
    chunk->haps = (int*)realloc(chunk->haps, m_reads * sizeof(int));
    chunk->phase_sets = (hts_pos_t*)realloc(chunk->phase_sets, m_reads * sizeof(hts_pos_t));
    for (int i = chunk->m_reads; i < m_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
        chunk->is_skipped[i] = 0; chunk->haps[i] = 0; chunk->phase_sets[i] = 0;
    }
    chunk->m_reads = m_reads;
    return 0;
}

void bam_chunk_free(bam_chunk_t *chunk) {
    if (chunk->ref_seq != NULL) free(chunk->ref_seq);
    if (chunk->low_comp_cr != NULL) cr_destroy(chunk->low_comp_cr);
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) {
            for (int j = 0; j < chunk->digars[i].n_digar; j++) {
                if (chunk->digars[i].digars[j].type == BAM_CINS || chunk->digars[i].digars[j].type == BAM_CDIFF) 
                    free(chunk->digars[i].digars[j].alt_seq);
            }
            free(chunk->digars[i].bseq); free(chunk->digars[i].qual);
            free(chunk->digars[i].digars);
            cr_destroy(chunk->digars[i].noisy_regs);
        }
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < chunk->m_reads; i++) {
            bam_destroy1(chunk->reads[i]);
        }
        free(chunk->reads);
    }
    if (chunk->noisy_reg_to_reads != NULL) {
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; i++) {
            free(chunk->noisy_reg_to_reads[i]);
        } free(chunk->noisy_reg_to_reads);
    }
    if (chunk->noisy_reg_to_n_reads != NULL) free(chunk->noisy_reg_to_n_reads);
    if (chunk->chunk_noisy_regs != NULL) cr_destroy(chunk->chunk_noisy_regs);
    if (chunk->cand_vars != NULL) free_cand_vars(chunk->cand_vars, chunk->n_cand_vars);
    if (chunk->var_i_to_cate != NULL) free(chunk->var_i_to_cate);

    if (chunk->read_var_profile != NULL) free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    if (chunk->read_var_cr != NULL) cr_destroy(chunk->read_var_cr);
    free(chunk->digars); free(chunk->is_skipped); free(chunk->haps); free(chunk->phase_sets);
    if (chunk->up_ovlp_read_i != NULL) free(chunk->up_ovlp_read_i);
    if (chunk->down_ovlp_read_i != NULL) free(chunk->down_ovlp_read_i);
}

void bam_chunk_post_free(bam_chunk_t *chunk) {
    if (chunk->ref_seq != NULL) free(chunk->ref_seq);
    if (chunk->cand_vars != NULL) free_cand_vars(chunk->cand_vars, chunk->n_cand_vars);
    if (chunk->var_i_to_cate != NULL) free(chunk->var_i_to_cate);
    free(chunk->is_skipped); free(chunk->haps); free(chunk->phase_sets);
    if (chunk->up_ovlp_read_i != NULL) free(chunk->up_ovlp_read_i);
    if (chunk->down_ovlp_read_i != NULL) free(chunk->down_ovlp_read_i);
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < chunk->m_reads; i++) {
            bam_destroy1(chunk->reads[i]);
        }
        free(chunk->reads);
    }
}

// free variables that will not be used in stitch & make variants
void bam_chunk_mid_free(bam_chunk_t *chunk) {
    if (chunk->low_comp_cr != NULL) cr_destroy(chunk->low_comp_cr);
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) {
            for (int j = 0; j < chunk->digars[i].n_digar; j++) {
                if (chunk->digars[i].digars[j].type == BAM_CINS || chunk->digars[i].digars[j].type == BAM_CDIFF) 
                    free(chunk->digars[i].digars[j].alt_seq);
            }
            free(chunk->digars[i].bseq); free(chunk->digars[i].qual);
            free(chunk->digars[i].digars);
            cr_destroy(chunk->digars[i].noisy_regs);
        }
    }
    free(chunk->digars); 
    if (chunk->chunk_noisy_regs != NULL) cr_destroy(chunk->chunk_noisy_regs);
    if (chunk->noisy_reg_to_reads != NULL) {
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; i++) {
            free(chunk->noisy_reg_to_reads[i]);
        } free(chunk->noisy_reg_to_reads);
    }
    if (chunk->noisy_reg_to_n_reads != NULL) free(chunk->noisy_reg_to_n_reads);
    if (chunk->read_var_profile != NULL) free_read_var_profile(chunk->read_var_profile, chunk->n_reads);
    if (chunk->read_var_cr != NULL) cr_destroy(chunk->read_var_cr);
}

void bam_chunks_mid_free(bam_chunk_t *chunks, int n_chunks) {
    for (int i = 0; i < n_chunks; i++) {
        bam_chunk_mid_free(chunks+i);
    }
}

void bam_chunks_post_free(bam_chunk_t *chunks, int n_chunks) {
    for (int i = 0; i < n_chunks; i++) {
        bam_chunk_post_free(chunks+i);
    }
    free(chunks);
}

void get_bam_chunk_reg_ref_seq0(faidx_t *fai, bam_chunk_t *chunk) {
    assert(fai != NULL); // faidx_t *fai = fai_load(ref_fasta);
    int ref_seq_len = faidx_seq_len(fai, chunk->tname);
    int flank_len = 50000; // [reg_beg-flank_len, reg_end+flank_len]
    hts_pos_t ref_beg = MAX_OF_TWO(flank_len, chunk->reg_beg-1)-flank_len;
    hts_pos_t ref_end = MIN_OF_TWO(ref_seq_len-flank_len-1, chunk->reg_end-1)+flank_len;
    chunk->whole_ref_len = ref_seq_len;

    int len;
    chunk->ref_seq = faidx_fetch_seq(fai, chunk->tname, ref_beg, ref_end, &len); // ref_beg & ref_end: 0-based
    chunk->ref_beg = ref_beg+1;   // 1-based
    chunk->ref_end = ref_beg+len; // 1-based

    // collect low-complexity regions
    chunk->low_comp_cr = cr_init();
    uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
    r = sdust(0, (uint8_t*)chunk->ref_seq-chunk->ref_beg+chunk->reg_beg, chunk->reg_end - chunk->reg_beg+1, T, W, &n);
    for (int i = 0; i < n; ++i) {
        // fprintf(stderr, "low_comp_cr: %s %d-%d\n", chunk->tname, chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1);
        cr_add(chunk->low_comp_cr, "cr", chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1, 0);
    }
    cr_index(chunk->low_comp_cr); free(r);
}

static int is_ovlp_with_prev_region(const struct call_var_pl_t *pl, bam_chunk_t *chunk, bam1_t *read) {
    int tid = chunk->tid, reg_chunk_i = chunk->reg_chunk_i, reg_i = chunk->reg_i;
    if (reg_i <= 0) return 0;
    int pre_reg_chunk_i = reg_chunk_i, pre_reg_i = reg_i-1;
    int pre_tid = pl->reg_chunks[pre_reg_chunk_i].reg_tids[pre_reg_i];
    if (tid != pre_tid) return 0;
    hts_pos_t pre_reg_beg = pl->reg_chunks[pre_reg_chunk_i].reg_begs[pre_reg_i];
    hts_pos_t pre_reg_end = pl->reg_chunks[pre_reg_chunk_i].reg_ends[pre_reg_i];
    hts_pos_t read_beg = read->core.pos+1, read_end = bam_endpos(read);
    if (read_end < pre_reg_beg || read_beg > pre_reg_end) return 0;
    else {
        // fprintf(stderr, "ovlp_pre: %s %" PRIi64 "-%" PRIi64 " %" PRId64 "-%" PRId64 "\n", bam_get_qname(read), read_beg, read_end, chunk->reg_beg, chunk->reg_end);
        return 1;
    }
}

static int is_ovlp_with_next_region(const struct call_var_pl_t *pl, bam_chunk_t *chunk, bam1_t *read) {
    int tid = chunk->tid, reg_chunk_i = chunk->reg_chunk_i, reg_i = chunk->reg_i;
    if (reg_i >= pl->reg_chunks[reg_chunk_i].n_regions-1) return 0;
    int next_reg_chunk_i = reg_chunk_i, next_reg_i = reg_i+1;
    int next_tid = pl->reg_chunks[next_reg_chunk_i].reg_tids[next_reg_i];
    if (tid != next_tid) return 0;
    hts_pos_t next_reg_beg = pl->reg_chunks[next_reg_chunk_i].reg_begs[next_reg_i];
    hts_pos_t next_reg_end = pl->reg_chunks[next_reg_chunk_i].reg_ends[next_reg_i];
    hts_pos_t read_beg = read->core.pos+1, read_end = bam_endpos(read);
    if (read_end < next_reg_beg || read_beg > next_reg_end) return 0;
    else {
        // fprintf(stderr, "ovlp_next: %s %" PRIi64 "-%" PRIi64 " %" PRId64 "-%" PRId64 "\n", bam_get_qname(read), read_beg, read_end, chunk->reg_beg, chunk->reg_end);
        return 1;
    }
}

// load ref_seq/read in reg_chunks[reg_chunk_i]->tid/beg/end[reg_i] to chunks
int collect_ref_seq_bam_main(const struct call_var_pl_t *pl, struct call_var_io_aux_t *io_aux, int reg_chunk_i, int reg_i, bam_chunk_t *chunk) {
    assert(reg_chunk_i < pl->n_reg_chunks); assert(reg_i < pl->reg_chunks[reg_chunk_i].n_regions);
    int tid = pl->reg_chunks[reg_chunk_i].reg_tids[reg_i];
    hts_pos_t reg_beg = pl->reg_chunks[reg_chunk_i].reg_begs[reg_i];
    hts_pos_t reg_end = pl->reg_chunks[reg_chunk_i].reg_ends[reg_i];

    // create iterator for the region
    hts_itr_t *iter = sam_itr_queryi(io_aux->idx, tid, reg_beg-1, reg_end); // (reg_beg-1, reg_end] ==> [reg_beg, reg_end]
    if (iter == NULL) {
        _err_warning("Failed to create iterator for region: %s:%" PRId64 "-%" PRId64 "\n", io_aux->header->target_name[tid], reg_beg, reg_end);
        return -1;
    } else {
        call_var_opt_t *opt = pl->opt;
        // load bam records
        samFile *in_bam = io_aux->bam; bam_hdr_t *header = io_aux->header;
        int min_mapq = opt->min_mq;
        bam_chunk_init0(chunk, 4096);
        chunk->reg_chunk_i = reg_chunk_i; chunk->reg_i = reg_i;
        chunk->tid = tid; chunk->tname = header->target_name[tid]; chunk->reg_beg = reg_beg; chunk->reg_end = reg_end;
        while (sam_itr_next(in_bam, iter, chunk->reads[chunk->n_reads]) >= 0) {
            if (chunk->reads[chunk->n_reads]->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue; // BAM_FSUPPLEMENTARY
            if (chunk->reads[chunk->n_reads]->core.qual < min_mapq) continue;
            // fprintf(stderr, "CHUNK-READ: %s %d-%d\n", bam_get_qname(chunk->reads[chunk->n_reads]), chunk->reg_beg, chunk->reg_end);
            // check if read is overlapping with previous region
            if (is_ovlp_with_prev_region(pl, chunk, chunk->reads[chunk->n_reads])) {
                chunk->up_ovlp_read_i[chunk->n_up_ovlp_reads++] = chunk->n_reads;
            }
            if (is_ovlp_with_next_region(pl, chunk, chunk->reads[chunk->n_reads])) {
                chunk->down_ovlp_read_i[chunk->n_down_ovlp_reads++] = chunk->n_reads;
            }
            // check if read is overlapping with next region
            chunk->n_reads++;
            if (chunk->n_reads == chunk->m_reads) bam_chunk_realloc(chunk);
        }
        bam_itr_destroy(iter);
        if (chunk->n_reads <= 0) return 0;
        // load ref seq
        get_bam_chunk_reg_ref_seq0(io_aux->fai, chunk);
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "CHUNK: tname: %s, tid: %d, beg: %" PRId64 ", end: %" PRId64 ", n_reads: %d\n", chunk->tname, chunk->tid, chunk->reg_beg, chunk->reg_end, chunk->n_reads);
        }
        return chunk->n_reads;
    }

}

// open & read bam, then output phased bam
int write_read_to_bam(bam_chunk_t *chunk, const struct call_var_opt_t *opt) {
    int n_out_reads = 0;
    int min_mapq = opt->min_mq;
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end; int tid = chunk->tid;
    samFile *in_bam = sam_open(opt->in_bam_fn, "r"); if (in_bam == NULL) _err_error_exit("Failed to open alignment file \'%s\'\n", opt->in_bam_fn);
    const htsFormat *fmt = hts_get_format(in_bam); if (fmt == NULL) _err_error_exit("Failed to get format of alignment file \'%s\'\n", opt->in_bam_fn);
    if (fmt->format != bam && fmt->format != cram) {
        sam_close(in_bam); _err_error_exit("Input file must be BAM or CRAM format.\n");
    }
    hts_tpool *p = hts_tpool_init(MAX_OF_TWO(1, opt->pl_threads));
    htsThreadPool thread_pool = {p, 0};
    if (hts_set_thread_pool(in_bam, &thread_pool) != 0) _err_error_exit("Failed to set thread pool.\n");
    bam_hdr_t *header = sam_hdr_read(in_bam); if (header == NULL) _err_error_exit("Failed to read header \'%s\'\n", opt->in_bam_fn);
    faidx_t *fai = fai_load3(opt->ref_fa_fn, opt->ref_fa_fai_fn, NULL, FAI_CREATE);
    if (fai == NULL) _err_error_exit("Failed to load/build reference fasta index \'%s\'\n", opt->ref_fa_fn);
    if (fmt->format == cram) {
        if (hts_set_fai_filename(in_bam, opt->ref_fa_fn) != 0) {
            fai_destroy(fai); hts_tpool_destroy(p); sam_hdr_destroy(header); sam_close(in_bam);
            _err_error_exit("Failed to set reference file for CRAM decoding: %s %s\n", opt->in_bam_fn, opt->ref_fa_fn);
        }
    }
    hts_idx_t *idx = sam_index_load(in_bam, opt->in_bam_fn); if (idx == NULL) _err_error_exit("Failed to load index \'%s\'\n", opt->in_bam_fn);
    hts_itr_t *iter = sam_itr_queryi(idx, tid, reg_beg-1, reg_end); // (reg_beg-1, reg_end] ==> [reg_beg, reg_end]
    if (iter == NULL) {
        _err_warning("Failed to create iterator for region: %s:%" PRId64 "-%" PRId64 ". Skipping.\n", header->target_name[tid], reg_beg, reg_end);
    } else {
        // load bam records
        bam1_t *read = bam_init1();
        int read_i = 0;
        int min_mapq = opt->min_mq;
        while (sam_itr_next(in_bam, iter, read) >= 0) {
            if (read->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue; // BAM_FSUPPLEMENTARY
            if (read->core.qual < min_mapq) continue;
            // check if read is overlapping with previous region
            if (read_i >= chunk->n_up_ovlp_reads) {
                // Add HP tag
                int hap = chunk->haps[read_i];
                if (hap != 0) {
                    // check if HP tag exists
                    uint8_t *hp = bam_aux_get(read, "HP");
                    if (hp != NULL) {
                        if (hap != bam_aux2i(hp)) {
                            bam_aux_del(read, hp);
                            bam_aux_append(read, "HP", 'i', 4, (uint8_t *)&(hap));
                        }
                    } else {
                        bam_aux_append(read, "HP", 'i', 4, (uint8_t *)&(hap));
                    }
                }
                // Add PS tag
                hts_pos_t ps = chunk->phase_sets[read_i];
                if (ps != -1) {
                    // check if PS tag exists
                    uint8_t *ps_tag = bam_aux_get(read, "PS");
                    if (ps_tag != NULL) {
                        if (ps != bam_aux2i(ps_tag))
                        {
                            bam_aux_del(read, ps_tag);
                            bam_aux_append(read, "PS", 'i', 4, (uint8_t *)&(ps));
                        }
                    } else {
                        bam_aux_append(read, "PS", 'i', 4, (uint8_t *)&(ps));
                    }
                }
                // write to bam
                if (sam_write1(opt->out_bam, header, read) < 0) _err_error_exit("Failed to write BAM record. %s:%" PRIi64 " %s\n", chunk->tname, read->core.pos+1, bam_get_qname(read));
                n_out_reads++;
            }
            read_i++;
        }
        bam_itr_destroy(iter); bam_destroy1(read);
    }
    fai_destroy(fai); hts_idx_destroy(idx);
    sam_hdr_destroy(header); sam_close(in_bam); hts_tpool_destroy(p); 
    return n_out_reads;
}

// check if multiple RG/SM tag are the same, if not same, output warning message
char *extract_sample_name_from_bam_header(bam_hdr_t *header) {
    int n_rg = sam_hdr_count_lines(header, "RG");
    char *sample_name = NULL;
    for (int i = 0; i < n_rg; ++i) {
        kstring_t ks = KS_INITIALIZE;
        if (sam_hdr_find_tag_pos(header, "RG", i, "SM", &ks) == 0) {
            if (i == 0) {
                sample_name = ks.s;
            } else {
                if (strcmp(sample_name, ks.s) != 0) {
                    free(ks.s);
                    _err_warning("Multiple RG/SM found in the BAM/CRAM header, using the first one: %s\n", sample_name);
                    return sample_name;
                } free(ks.s);
            }
        }
    }
    if (sample_name != NULL) _err_info("Sample name from BAM/CRAM header: %s\n", sample_name);
    return sample_name;
}
