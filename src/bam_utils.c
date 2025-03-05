#include <ctype.h>
#include "bam_utils.h"
#include "utils.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "sdust.h"

extern int LONGCALLD_VERBOSE;

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
} xid_queue_t;

xid_queue_t *init_xid_queue(int max_sites, int max_s, int win) {
    xid_queue_t *q = (xid_queue_t*)malloc(sizeof(xid_queue_t));
    q->pos = (hts_pos_t*)malloc(max_sites * sizeof(hts_pos_t));
    q->counts = (int*)malloc(max_sites * sizeof(int));
    q->lens = (int*)malloc(max_sites * sizeof(int));
    q->is_dense = (int*)calloc(max_sites, sizeof(int));
    q->front = 0; q->rear = -1; q->count = 0;
    q->max_s = max_s; q->win = win;
    return q;
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
    q->pos[++q->rear] = pos; // left-most pos of the region
    q->lens[q->rear] = len; q->counts[q->rear] = count;
    q->count += count;

    hts_pos_t noisy_start = -1, noisy_end = -1;
    while (q->pos[q->front]+q->lens[q->front]-1 <= pos - q->win) {
        q->count -= q->counts[q->front]; // q->lens[q->front];
        q->front++;
    }
    // fprintf(stderr, "pos: %ld, len: %d, count: %d, front: %d, rear: %d, count: %d\n", pos, len, count, q->front, q->rear, q->count);
    if (count > 0) {
        if (q->count > q->max_s) {
            int max_ins_len = 0;
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

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD) {
    int n_to_check_reads = 10;
    // if any of n_to_check_reads reads do NOT have eqx or MD tag
    // set tag has NO
    bam1_t *b = bam_init1();
    int r, n_checked_reads = 0;
    *has_eqx = 1; *has_MD = 1;
    while ((r = sam_read1(in_bam, header, b)) >= 0) {
        if (has_equal_X_in_bam_cigar(b) == 0) *has_eqx = 0;
        if (has_MD_in_bam(b) == 0) *has_MD = 0;
        n_checked_reads++;
        if ((*has_eqx == 0 && *has_MD == 0) || n_checked_reads >= n_to_check_reads)
            break;
    }
    bam_destroy1(b);
    // rewind the file
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

// cand_vars:
// X: 1 -> M SNPs
// I: 0 -> M INSs
// D: M -> 1 DEL
// XXX no minor_alt_allele here, will be used in read_var_profile
// here, minor_alt_allele will be simply treated as non-alt, i.e., ref allele
int update_cand_vars_from_digar(digar_t *digar, bam1_t *read, int n_var_sites, var_site_t *var_sites, int start_i, cand_var_t *cand_vars) {
    int cur_start_i, site_i, digar_i = 0;
    hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    digar1_t *digars = digar->digars; uint8_t strand = bam_is_rev(read);

    cur_start_i = site_i = get_var_site_start(var_sites, start_i, n_var_sites, pos_start);
    int tid = read->core.tid;
    // compare sorted var_sites and digars
    // both var_sites and digars are sorted by pos/type/ref_len
    for (; site_i < n_var_sites && digar_i < digar->n_digar; ) {
        if (digars[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar->digars+digar_i);
        int ret = comp_var_site(var_sites+site_i, &digar_var_site);
        if (ret < 0) { // var_site < digar_var_site
            // update_var_site_with_allele(cand_vars+site_i, is_in_noisy_reg(cand_vars[site_i].pos, digar->noisy_regs), strand, 0);
            update_var_site_with_allele(cand_vars+site_i, 0, strand, 0);
            site_i++;
        } else if (ret == 0) { // merge together, update/add alt allele
            update_var_site_with_allele(cand_vars+site_i, digar->digars[digar_i].is_low_qual, strand, 1);
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

int update_read_var_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_vars, cand_var_t *cand_vars, int start_var_i, read_var_profile_t *read_var_profile) {
    int cur_start_i, var_i, digar_i = 0, allele_i=-1;
    hts_pos_t pos_start = read->core.pos+1, pos_end = bam_endpos(read);
    digar1_t *digar1 = digar->digars;
    cur_start_i = var_i = get_var_start(cand_vars, start_var_i, n_cand_vars, pos_start);
    int tid = read->core.tid;
    // compare sorted cand_vars and digars
    // both cand_vars and digars are sorted by pos/type/ref_len
    for (; var_i < n_cand_vars && digar_i < digar->n_digar; ) {
        if (digar1[digar_i].type == BAM_CEQUAL) {
            digar_i++; continue;
        }
        var_site_t var_site0 = make_var_site_from_cand_var(cand_vars+var_i);
        // XXX should always update read_var_profile
        // if (is_in_noisy_reg(var_site0.pos, digar->noisy_regs)) {
        //     var_i++; continue;
        // }
        var_site_t digar_var_site = make_var_site_from_digar(tid, digar1+digar_i);
        int is_ovlp, ret;
        ret = comp_ovlp_var_site(&var_site0, &digar_var_site, &is_ovlp);
        if (is_ovlp == 0) { // no overlapping
            if (ret < 0) { // var_site < digar_var_site: ref_allele
                update_read_var_profile_with_allele(var_i, 0, read_var_profile); //is_in_noisy_reg(var_site0.pos, digar->noisy_regs));
                var_i++;
            } else if (ret > 0) { // var_site > digar_var_site
                digar_i++;
            } else { // var_site == digar_var_site
                _err_error_exit("Unexpected case: is_ovlp == 0 && var_site == digar_var_site, %d-%c-%d-%d\n", cand_vars[var_i].pos, BAM_CIGAR_STR[cand_vars[var_i].var_type], digar_var_site.pos, BAM_CIGAR_STR[digar_var_site.var_type]);
                // if (var_site0.var_type != BAM_CINS || digar_var_site.var_type != BAM_CINS) {
                    // _err_error_exit("Unexpected case: is_ovlp == 0 && var_site == digar_var_site, %d-%c-%d-%d\n", cand_vars[var_i].pos, BAM_CIGAR_STR[cand_vars[var_i].var_type], digar_var_site.pos, BAM_CIGAR_STR[digar_var_site.var_type]);
                // }
                // allele_i = get_alt_allele_idx(cand_vars+var_i, digar->digars+digar_i, digar->bseq);
                // update_read_var_profile_with_allele(var_i, allele_i, read_var_profile); // is_in_noisy_reg(var_site0.pos, digar->noisy_regs));
                // var_i++;
            }
        } else { // overlap
            if (ret == 0) { // exact the same with var_site
                allele_i = get_alt_allele_idx(cand_vars+var_i, digar->digars+digar_i, digar->bseq);
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

void set_digar(digar1_t *digar, hts_pos_t pos, int type, int len, int qi, int is_low_qual) {
    // fprintf(stderr, "Set digar: %ld, %d%c, %d\n", pos, len, BAM_CIGAR_STR[type], qi);
    // uint8_t d_qual; if (qi > 0) d_qual = (qual[qi] + qual[qi-1]) / 2; else d_qual = qual[qi];
    digar->pos = pos; digar->type = type; digar->len = len; digar->qi = qi; digar->alt_seq = NULL; digar->is_low_qual = is_low_qual;
}

void set_seq_digar(digar1_t *digar, hts_pos_t pos, int type, int len, int qi, int is_low_qual, uint8_t *alt_seq) {
    digar->pos = pos; digar->type = type; digar->len = len; digar->qi = qi; digar->is_low_qual = is_low_qual; digar->alt_seq = alt_seq;
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
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "palindrome: %s %ld-%ld, sa: %s %d-%d\n", bam_get_qname(read), primary_pos, primary_end, rname, pos, sa_end_pos);
                    return 0;
                }
            }
        }
    }
    return 1;
}

// is_low_qual: low base quality
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar) {
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    hts_pos_t reg_beg = chunk->active_reg_beg, reg_end = chunk->active_reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    // int dens_reg_flank_win = opt->dens_reg_flank_win, indel_flank_win = opt->indel_flank_win;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;
    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    digar->beg = pos; digar->end = bam_endpos(read);
    digar->bseq = bam_get_seq(read); digar->qual = bam_get_qual(read);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = faidx_seq_len(pl->fai, chunk->tname);
    xid_queue_t *q = init_xid_queue(rlen, max_s, win); // for noisy region
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                // XXX qual[qi] == 255 : base quality not available
                // if (qual[qi] >= min_bq) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                if (digar->qual[qi] >= opt->min_bq) { // skip low-quality bases for noisy regions
                    push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                    set_seq_digar(_digars+_n_digar, pos, op, 1, qi, 0, x_seq);
                } else set_seq_digar(_digars+_n_digar, pos, op, 1, qi, 1, x_seq);
                _n_digar++; // push_xid_queue(q, pos, 1, 1);
                n_total_cand_vars++;
                pos++; qi++;
            }
        } else if (op == BAM_CEQUAL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0);
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            // push_xid_queue_win(q, pos, len, 0, digar->noisy_regs, &noisy_start, &noisy_end);
            pos += len; qi += len;
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, op, len, qi, 0);
            } else set_digar(_digars+_n_digar, pos, op, len, qi, 1);
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
            set_seq_digar(_digars+_n_digar, pos, op, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0); // clipping
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
            qi += len;
        } else if (op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, 0); // clipping
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (i == 0) {
                    // the very start of the contig
                    if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                } else {
                    if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                }
            }
            // push_xid_queue_win(q, pos, 0, 0, digar->noisy_regs, &noisy_start, &noisy_end);
        } else if (op == BAM_CREF_SKIP) { // XXX no action for N op
            pos += len;
        } else if (op == BAM_CMATCH) {
            _err_error_exit("CIGAR operation 'M' is not expected in EQX CIGAR: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0, j = 0; i < _n_digar; ++i) { // XXX not set noisy-region for digars
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
        var_size = var_size > (noisy_end - noisy_start + 1) ? var_size : (noisy_end - noisy_start + 1);
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
            if (is_overlap_reg(cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), reg_beg, reg_end)) {
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
int collect_digar_from_cs_tag(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar, cgranges_t *chunk_noisy_regs) {
    return 0;
}

int collect_digar_from_MD_tag(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar) {
    hts_pos_t reg_beg = chunk->active_reg_beg, reg_end = chunk->active_reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    uint8_t *s = bam_aux_get(read, "MD");
    if (s == NULL) { _err_error_exit("MD tag not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    digar->bseq = bam_get_seq(read); digar->qual = bam_get_qual(read);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *md = bam_aux2Z(s); int md_i = 0;
    // printf("MD: %s\n", md);
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = faidx_seq_len(pl->fai, chunk->tname);
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
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0);
                        _n_digar++; // push_xid_queue(q, pos, 0, 0);
                        pos += eq_len; qi += eq_len;
                        last_eq_len -= m_len; m_len = 0;
                    } else { // last_eq_len < m_len
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, last_eq_len, qi, 0);
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
                    set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, 0); // qual unset XXX
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
                        set_seq_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_seq_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
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
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1);
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
            set_seq_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CSOFT_CLIP, len, qi, 0); // clipping
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
            qi += len;
        } else if (op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0); // insertion
            _n_digar++; // push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (i == 0) {
                    // the very start of the contig
                    if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                } else {
                    if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                }
            }
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0, j = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
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
            if (is_overlap_reg(cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), cr_label(digar->noisy_regs, i));
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
int collect_digar_from_ref_seq(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar) {
    return 0;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    hts_pos_t reg_beg = chunk->active_reg_beg, reg_end = chunk->active_reg_end; cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t pos = read->core.pos+1, qi = 0;
    digar->beg = pos; digar->end = bam_endpos(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int max_s = opt->noisy_reg_max_xgaps, win = opt->noisy_reg_slide_win;
    double max_noisy_frac_per_read = opt->max_noisy_frac_per_read, max_var_ratio_per_read = opt->max_var_ratio_per_read;
    int end_clip_reg = opt->end_clip_reg, end_clip_reg_flank_win = opt->end_clip_reg_flank_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    digar->noisy_regs = cr_init();
    digar->bseq = bam_get_seq(read); digar->qual = bam_get_qual(read);
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar); int tlen = faidx_seq_len(pl->fai, chunk->tname);
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);
    hts_pos_t noisy_start = -1, noisy_end = -1; int cr_q_start = -1, cr_q_end = -1;
    int n_total_cand_vars = 0;

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH ) {
            int eq_len = 0;
            for (int j = 0; j < len; ++j) {
                // Get the reference base
                if (pos <= ref_beg || pos > ref_end) _err_error_exit("Read exceed reference region sequence: %s", bam_get_qname(read));
                char ref_base = ref_seq[pos-ref_beg];
                // char ref_base = get_ref_base_from_cr(ref_seq, chunk_cr, pos);
                // Get the read base
                char read_base = seq_nt16_str[bam_seqi(digar->bseq, qi)];
                if (ref_base != read_base) {
                    if (eq_len > 0) {
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0);
                        _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                        eq_len = 0;
                    }
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    uint8_t *x_seq = (uint8_t*)malloc(sizeof(uint8_t));
                    x_seq[0] = seq_nt16_int[bam_seqi(digar->bseq, qi)];
                    if (digar->qual[qi] >= opt->min_bq) {
                        push_xid_size_queue_win(q, pos, 1, 1, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                        set_seq_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 0, x_seq);
                    } else set_seq_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, 1, x_seq);
                    _n_digar++;
                    n_total_cand_vars++;
                } else eq_len++;
                pos++; qi++;
            }
            if (eq_len > 0) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, 0);
                _n_digar++; // push_xid_queue(q, pos-eq_len, 0, 0);
                eq_len = 0;
            }
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            if ((qi == 0 || digar->qual[qi-1] >= opt->min_bq) && digar->qual[qi] >= opt->min_bq) {
                push_xid_size_queue_win(q, pos, len, len, digar->noisy_regs, &noisy_start, &noisy_end, &cr_q_start, &cr_q_end);
                set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 0);
            } else set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, 1);
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
            set_seq_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, is_low_qual, ins_seq);
            _n_digar++; //push_xid_queue(q, pos, 0, 1);
            n_total_cand_vars++;
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CSOFT_CLIP, len, qi, 0); // clipping 
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
            qi += len;
        } else if (op == BAM_CHARD_CLIP) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CHARD_CLIP, len, qi, 0); // clipping 
            _n_digar++; //push_xid_queue(q, pos, 0, 0);
            if (len > end_clip_reg) {
                if (i == 0) {
                    // the very start of the contig
                    if (pos > 1) cr_add(digar->noisy_regs, "cr", pos-1, pos+end_clip_reg_flank_win, 0); // left end
                } else {
                    if (pos < tlen) cr_add(digar->noisy_regs, "cr", pos-1-end_clip_reg_flank_win, pos, 0); // right end
                }
            }
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s\n", bam_get_qname(read));
        }
    }
    for (int i = 0, j = 0; i < _n_digar; ++i) {
        push_digar1(digar, _digars[i]);
    }
    if (noisy_start != -1) {
        int var_size = 0;
        for (int i = cr_q_start; i <= cr_q_end; ++i) var_size += q->counts[i];
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
            if (is_overlap_reg(cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), reg_beg, reg_end))
                cr_add(chunk_noisy_regs, "cr", cr_start(digar->noisy_regs, i), cr_end(digar->noisy_regs, i), 1);
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "DIGAR3: %s\n", bam_get_qname(read));
            print_digar(digar, stderr);
        }
    }
    free(_digars); free_xid_queue(q);
    return skip == 1 ? -1 : 0;
}

void copy_bam_chunk(bam_chunk_t *from_chunk, bam_chunk_t *to_chunk, int *chunk_read_i, int n_reads, hts_pos_t *chunk_read_beg) {
    to_chunk->up_ovlp_read_i = (int*)malloc(n_reads * sizeof(int));
    for (int i = 0; i < n_reads; i++) {
        int from_read_i = chunk_read_i[i];
        if (bam_copy1(to_chunk->reads[to_chunk->n_reads + i], from_chunk->reads[from_read_i]) == NULL)
            _err_error_exit("Failed to copy BAM record: %s (%d,%d)", bam_get_qname(from_chunk->reads[from_read_i]), from_read_i, n_reads);
        if (*chunk_read_beg == -1) *chunk_read_beg = from_chunk->reads[from_read_i]->core.pos;
        // to_chunk->reads[to_chunk->n_reads + i] =  from_chunk->reads[read_i];
        if (from_chunk->is_ovlp[from_read_i]) {
            // fprintf(stderr, "ovlp-read: %s\n", bam_get_qname(from_chunk->reads[from_read_i]));
            to_chunk->up_ovlp_read_i[to_chunk->n_up_ovlp_reads++] = from_read_i;
            to_chunk->is_ovlp[to_chunk->n_reads + i] = 1;
        }
        to_chunk->is_skipped[to_chunk->n_reads + i] = 0;
        to_chunk->haps[to_chunk->n_reads + i] = 0;
        to_chunk->PS[to_chunk->n_reads + i] = -1;
    }
    to_chunk->n_reads += n_reads;
}

void copy_bam_chunk0(bam_chunk_t *from_chunk, bam_chunk_t *to_chunk) {
    to_chunk->tid = from_chunk->tid; to_chunk->active_reg_beg = from_chunk->active_reg_beg; to_chunk->active_reg_end = from_chunk->active_reg_end;
    to_chunk->n_reads = to_chunk->m_reads = from_chunk->n_reads;
    to_chunk->reads = (bam1_t**)malloc(from_chunk->n_reads * sizeof(bam1_t*));
    to_chunk->is_ovlp = (uint8_t*)calloc(from_chunk->n_reads, sizeof(uint8_t));
    to_chunk->haps = (int*)malloc(from_chunk->n_reads * sizeof(int));
    to_chunk->PS = (hts_pos_t*)malloc(from_chunk->n_reads * sizeof(hts_pos_t));
    to_chunk->is_skipped = (uint8_t*)malloc(from_chunk->n_reads * sizeof(uint8_t));
    for (int read_i = 0; read_i < from_chunk->n_reads; read_i++) {
        to_chunk->reads[read_i] = bam_init1();
        if (bam_copy1(to_chunk->reads[read_i], from_chunk->reads[read_i]) == NULL)
            _err_error_exit("Failed to copy BAM record: %s (%d,%d)", bam_get_qname(from_chunk->reads[read_i]), read_i, from_chunk->n_reads);
        to_chunk->is_ovlp[read_i] = from_chunk->is_ovlp[read_i];
        to_chunk->haps[read_i] = from_chunk->haps[read_i];
        to_chunk->PS[read_i] = from_chunk->PS[read_i];
        to_chunk->is_skipped[read_i] = from_chunk->is_skipped[read_i];
        // fprintf(stderr, "copy_bam_chunk0: ovlp-read: %s\n", bam_get_qname(from_chunk->reads[read_i]));
        // fprintf(stderr, "hap: %d, PS: %ld\n", to_chunk->haps[read_i], to_chunk->PS[read_i]);
    }
}

int bam_chunk_init(bam_chunk_t *chunk, int n_reads, bam_chunk_t *last_chunk, int *last_chunk_read_i, int n_last_chunk_reads, hts_pos_t *chunk_read_beg) {
    if (n_last_chunk_reads > n_reads) n_reads = n_last_chunk_reads;

    // input
    chunk->n_reads = 0; chunk->m_reads = n_reads;
    chunk->n_up_ovlp_reads = 0; chunk->up_ovlp_read_i = NULL;
    chunk->reads = (bam1_t**)malloc(n_reads * sizeof(bam1_t*));
    chunk->reg_cr = NULL; chunk->low_comp_cr = NULL;
    // intermediate
    chunk->is_skipped = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->is_ovlp = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->digars = (digar_t*)calloc(n_reads, sizeof(digar_t));
    for (int i = 0; i < n_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
    }
    // noisy regions
    chunk->chunk_noisy_regs = NULL; chunk->noisy_reg_to_reads = NULL; chunk->noisy_reg_to_n_reads = NULL;
    // variant
    chunk->n_cand_vars = 0; chunk->cand_vars = NULL;
    chunk->var_i_to_cate = NULL;
    chunk->read_var_profile = NULL; chunk->read_var_cr = NULL;
    // output
    chunk->haps = (int*)calloc(n_reads, sizeof(int));
    chunk->PS = (hts_pos_t*)malloc(n_reads * sizeof(hts_pos_t));
    for (int i = 0; i < n_reads; i++) chunk->PS[i] = -1;
    chunk->flip_hap = 0; chunk->pre_chunk_reg_end = -1; chunk->flip_pre_chunk_PS = -1; chunk->flip_cur_chunk_PS = -1;

    if (n_last_chunk_reads > 0) {
        copy_bam_chunk(last_chunk, chunk, last_chunk_read_i, n_last_chunk_reads, chunk_read_beg);
        return chunk->reads[0]->core.tid;
    } else return -1;
}

int bam_chunk_realloc(bam_chunk_t *chunk) {
    int m_reads = chunk->m_reads * 2;
    chunk->reads = (bam1_t**)realloc(chunk->reads, m_reads * sizeof(bam1_t*));
    chunk->is_skipped = (uint8_t*)realloc(chunk->is_skipped, m_reads * sizeof(uint8_t));
    chunk->is_ovlp = (uint8_t*)realloc(chunk->is_ovlp, m_reads * sizeof(uint8_t));
    chunk->digars = (digar_t*)realloc(chunk->digars, m_reads * sizeof(digar_t));
    chunk->haps = (int*)realloc(chunk->haps, m_reads * sizeof(int));
    chunk->PS = (hts_pos_t*)realloc(chunk->PS, m_reads * sizeof(hts_pos_t));
    for (int i = chunk->m_reads; i < m_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
        chunk->is_skipped[i] = 0; chunk->haps[i] = 0; chunk->PS[i] = 0;
        chunk->is_ovlp[i] = 0;
    }
    chunk->m_reads = m_reads;
    return 0;
}

void bam_ovlp_chunk_free(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->m_reads; i++) {
        bam_destroy1(chunk->reads[i]);
    }
    free(chunk->reads); free(chunk->is_ovlp);
    free(chunk->haps); free(chunk->PS); free(chunk->is_skipped);
}

void bam_chunk_free(bam_chunk_t *chunk) {
    if (chunk->need_free_ref_seq) free(chunk->ref_seq);
    if (chunk->reg_cr != NULL) cr_destroy(chunk->reg_cr);
    if (chunk->low_comp_cr != NULL) cr_destroy(chunk->low_comp_cr);
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) {
            for (int j = 0; j < chunk->digars[i].n_digar; j++) {
                if (chunk->digars[i].digars[j].type == BAM_CINS || chunk->digars[i].digars[j].type == BAM_CDIFF) 
                    free(chunk->digars[i].digars[j].alt_seq);
            }
            free(chunk->digars[i].digars);
            cr_destroy(chunk->digars[i].noisy_regs);
        }
    }
    for (int i = 0; i < chunk->m_reads; i++) {
        bam_destroy1(chunk->reads[i]);
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
    free(chunk->reads); free(chunk->digars);
    free(chunk->is_skipped); free(chunk->is_ovlp); free(chunk->haps); free(chunk->PS);
    if (chunk->up_ovlp_read_i != NULL) free(chunk->up_ovlp_read_i);
}

void bam_chunks_free(bam_chunk_t *chunks, int n_chunks) {
    for (int i = 0; i < n_chunks; i++) {
        bam_chunk_free(chunks+i);
    }
    free(chunks);
}

char *get_region_seq(ref_reg_seq_t *r, const char *tname, hts_pos_t beg, hts_pos_t end) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(r->reg_cr, tname, beg, end, &ovlp_b, &max_b);
    if (ovlp_n <= 0) {
        _err_error_exit("Region not found in the reference: %s:%d-%d", tname, beg, end);
    } else if (ovlp_n > 1) {
        _err_error_exit("Multiple regions found in the reference: %s:%d-%d", tname, beg, end);
    }
    int reg_i = cr_label(r->reg_cr, ovlp_b[0]);
    // if (beg < r->reg_seq[reg_i].beg || end > r->reg_seq[reg_i].end)
        // _err_error_exit("Region not found in the reference: %s:%d-%d", tname, beg, end);
    free(ovlp_b);
    return r->reg_seq[reg_i].seq.s + (beg - r->reg_seq[reg_i].beg);
}

void get_bam_chunk_reg_ref_seq(faidx_t *fai, ref_reg_seq_t *ref_reg_seq, bam_chunk_t *chunk) {
    cgranges_t *reg_cr = chunk->reg_cr;
    chunk->low_comp_cr = cr_init();
    assert(reg_cr->n_r > 0);
    if (reg_cr->n_r == 1) {
        char *tname = chunk->tname; hts_pos_t ref_beg = cr_start(reg_cr, 0), ref_end = cr_end(reg_cr, 0);
        // get the reference sequence
        int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
        ovlp_n = cr_overlap(ref_reg_seq->reg_cr, tname, ref_beg-1, ref_end, &ovlp_b, &max_b);
        assert(ovlp_n == 1);
        int ref_i = cr_label(ref_reg_seq->reg_cr, ovlp_b[0]); free(ovlp_b);

        chunk->ref_beg = ref_reg_seq->reg_seq[ref_i].beg;
        chunk->ref_end = ref_reg_seq->reg_seq[ref_i].end;
        chunk->need_free_ref_seq = 0;
        chunk->ref_seq = ref_reg_seq->reg_seq[ref_i].seq.s;
    } else { // if reg_cr->n_r > 1: allocate a new memory for the reference sequence, containing gaps between regions
        assert(fai != NULL); // faidx_t *fai = fai_load(ref_fasta);
        chunk->need_free_ref_seq = 1;
        hts_pos_t ref_beg = cr_start(reg_cr, 0), ref_end = cr_end(reg_cr, 0);
        for (int64_t i = 1; i < reg_cr->n_r; i++) {
            if (cr_start(reg_cr, i) < ref_beg) ref_beg = cr_start(reg_cr, i);
            if (cr_end(reg_cr, i) > ref_end) ref_end = cr_end(reg_cr, i);
        }
        int len;
        hts_pos_t _ref_beg = MAX_OF_TWO(0, ref_beg - 1000);
        hts_pos_t _ref_end = ref_end + 1000;
        chunk->ref_seq = faidx_fetch_seq(fai, chunk->tname, _ref_beg, _ref_end, &len);
        chunk->ref_beg = _ref_beg+1;
        chunk->ref_end = _ref_beg+len;
    }
    // collect low-complexity regions
    uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
    r = sdust(0, (uint8_t*)chunk->ref_seq-chunk->ref_beg+chunk->active_reg_beg, chunk->active_reg_end - chunk->active_reg_beg+1, T, W, &n);
    for (int i = 0; i < n; ++i) {
        cr_add(chunk->low_comp_cr, "cr", chunk->active_reg_beg+(int)(r[i]>>32)-1, chunk->active_reg_beg+(int)r[i]-1, 0);
    }
    cr_index(chunk->low_comp_cr); free(r);
}

void get_bam_chunk_reg_cr(cgranges_t *ref_seq_reg_cr, bam_chunk_t *chunk, hts_pos_t chunk_active_reg_beg, hts_pos_t chunk_active_reg_end) {
    chunk->reg_cr = cr_init();
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(ref_seq_reg_cr, chunk->tname, chunk_active_reg_beg-1, chunk_active_reg_end, &ovlp_b, &max_b);
    assert(ovlp_n > 0);
    hts_pos_t _beg = chunk_active_reg_end, _end = chunk_active_reg_beg-1;
    for (ovlp_i = 0; ovlp_i < ovlp_n; ovlp_i++) {
        // fprintf(stderr, "reg_cr, ovlp: %s:%d-%d\n", chunk->tname, cr_start(ref_seq_reg_cr, ovlp_b[ovlp_i]), cr_end(ref_seq_reg_cr, ovlp_b[ovlp_i]));
        cr_add(chunk->reg_cr, chunk->tname, cr_start(ref_seq_reg_cr, ovlp_b[ovlp_i]), cr_end(ref_seq_reg_cr, ovlp_b[ovlp_i]), cr_label(ref_seq_reg_cr, ovlp_b[ovlp_i]));
        if (cr_start(ref_seq_reg_cr, ovlp_b[ovlp_i]) < _beg) _beg = cr_start(ref_seq_reg_cr, ovlp_b[ovlp_i]);
        if (cr_end(ref_seq_reg_cr, ovlp_b[ovlp_i]) > _end) _end = cr_end(ref_seq_reg_cr, ovlp_b[ovlp_i]);
    }
    // XXX TODO: variant calling region: [reg_beg, reg_end]
    //           non-overlapping with nearby chunks, so that any variants will only be called in one chunk/region
    //           take care of the very first & last chunk
    chunk->active_reg_beg = MAX_OF_TWO(_beg+1, chunk_active_reg_beg); chunk->active_reg_end = MIN_OF_TWO(_end, chunk_active_reg_end);
    // fprintf(stderr, "CHUNK_REG: %ld %ld %ld %ld\n", _beg+1, chunk_active_reg_beg, _end, chunk_active_reg_end);
    free(ovlp_b);
    cr_index(chunk->reg_cr);
}

// collect reads from a BAM file, region: [start, start+max_reg_len_per_chunk]
// *n_ovlp_reads: number of overlapping reads between previous chunk and current chunk, update it to "overlapping reads between current chunk and next chunk" after collecting reads
// int collect_bam_chunk(samFile *in_bam, bam_hdr_t *header, hts_itr_t *iter, int use_iter, int max_reg_len_per_chunk, int **ovlp_read_i, int *n_ovlp_reads, hts_pos_t *last_reg_end, bam_chunk_t *chunk) {

// collect reads from a BAM file, region: [start, start+max_reg_len_per_chunk]
// **last_chunk_read_i: indices of reads from the last chunk that will be collected in current chunk
// *n_last_chunk_reads: number of reads in the last chunk that will be collected in current chunk
int collect_bam_chunk(call_var_pl_t *pl, bam_chunk_t *chunks, int chunk_i) {
    bam_chunk_t *chunk = chunks + chunk_i;
    // int **last_chunk_read_i, int *n_last_chunk_reads, hts_pos_t *cur_active_reg_beg, 
    samFile *in_bam = pl->bam; bam_hdr_t *header = pl->header; hts_itr_t *iter = pl->iter; 
    int use_iter = pl->use_iter; int max_reg_len_per_chunk = pl->max_reg_len_per_chunk, ovlp_reg_len = pl->ovlp_region_len;
    int r, min_mq = pl->opt->min_mq, has_eqx_cigar=0, has_MD=0; // has_cs=0; XXX
    int reg_tid=-1, tid0;
    hts_pos_t beg0=-1, end0=-1; // read-wise, end: end of all reads in current chunk, beg0/end0: start/end of the current read
    hts_pos_t active_reg_beg = -1, active_reg_end = -1, next_active_reg_beg = -1;
    hts_pos_t chunk_read_reg_beg=-1; // beg of the current chunk, based on all reads in the chunk

    // if *n_last_chunk_reads > 0, *cur_reg_beg was set during collecting the last chunk; tid will be set as the tid of last_chunk_reads
    // else if tid==-1, and *cur_reg_beg needs to be set as the start of the current region
    // *cur_active_reg_beg: -1 means cur_active_reg_beg is not set yet, set as the very first read in the chunk, including overlapping ones

    bam_chunk_t *last_chunk = pl->last_chunk; // last chunk
    if (chunk_i > 0) last_chunk = chunk-1;
    reg_tid = bam_chunk_init(chunk, 4096, last_chunk, pl->last_chunk_read_i, pl->n_last_chunk_reads, &chunk_read_reg_beg);
    // fprintf(stderr, "chunk_beg: %ld\n", chunk_read_reg_beg);
    if (pl->cur_active_reg_beg != -1) {
        // assert(reg_tid != -1);
        if (chunk_read_reg_beg != -1) pl->cur_active_reg_beg = MAX_OF_TWO(pl->cur_active_reg_beg, chunk_read_reg_beg);
        active_reg_beg = pl->cur_active_reg_beg;
        active_reg_end = active_reg_beg + max_reg_len_per_chunk + ovlp_reg_len - 1;
        // fprintf(stderr, "active_beg1: %ld active_end1: %ld\n", active_reg_beg, active_reg_end);
        next_active_reg_beg = active_reg_end - ovlp_reg_len * 2 + 1;
    } else {
        if (chunk_read_reg_beg != -1) {
            active_reg_beg = chunk_read_reg_beg;
            active_reg_end = active_reg_beg + max_reg_len_per_chunk + ovlp_reg_len - 1;
            next_active_reg_beg = active_reg_end - ovlp_reg_len * 2 + 1;
            // fprintf(stderr, "active_beg2: %ld active_end2: %ld\n", active_reg_beg, active_reg_end);
        } // if (chunk_read_reg_beg == -1): wait unit read the first read
    }
    if (pl->last_chunk_read_i != NULL) free(pl->last_chunk_read_i);
    int m_last_chunk_reads = chunk->m_reads; // push to last_chunk_read_i if end >= reg_end or tid != cur_tid, break if beg > reg_end or tid != cur_tid
    int _n_last_chunk_reads = 0;             // for tid != cur_id, only one read will be collect for next chunk (*n_last_chunk_reads=1)
    pl->last_chunk_read_i = (int*)malloc(m_last_chunk_reads * sizeof(int));
    while (1) {
        if (chunk->n_reads == chunk->m_reads) bam_chunk_realloc(chunk);
        if (use_iter) r = sam_itr_next(in_bam, iter, chunk->reads[chunk->n_reads]);
        else r = sam_read1(in_bam, header, chunk->reads[chunk->n_reads]);
        if (r < 0) break;
        // fprintf(stderr, "%s: %ld\n", bam_get_qname(chunk->reads[chunk->n_reads]), chunk->reads[chunk->n_reads]->core.pos);
        // skip unmapped/secondary/supplementary, low MAPQ alignments
        if (chunk->reads[chunk->n_reads]->core.flag & (BAM_FUNMAP | BAM_FSECONDARY | BAM_FSUPPLEMENTARY)) continue; // BAM_FSUPPLEMENTARY
        if (chunk->reads[chunk->n_reads]->core.qual < min_mq) continue;

        tid0 = chunk->reads[chunk->n_reads]->core.tid;
        beg0 = chunk->reads[chunk->n_reads]->core.pos + 1;
        end0 = bam_endpos(chunk->reads[chunk->n_reads]);

        if (reg_tid == -1) { // first read in the first chunk
            reg_tid = tid0;
            if (active_reg_beg == -1) {
                active_reg_beg = beg0;
                active_reg_end = active_reg_beg + max_reg_len_per_chunk + ovlp_reg_len - 1;
                next_active_reg_beg = active_reg_end - 2 * ovlp_reg_len + 1;
            }
            // if (end0 >= active_reg_end) { // extreamly long reads that span the whole region
            if (end0 >= next_active_reg_beg) {
                chunk->is_ovlp[chunk->n_reads] = 1; // overlap between current and next chunk
                if (_n_last_chunk_reads == m_last_chunk_reads) {
                    m_last_chunk_reads *= 2;
                    pl->last_chunk_read_i = (int*)realloc(pl->last_chunk_read_i, m_last_chunk_reads * sizeof(int));
                }
                pl->last_chunk_read_i[_n_last_chunk_reads++] = chunk->n_reads;
            }
            // fprintf(stderr, "active_beg3: %ld active_end3: %ld\n", active_reg_beg, active_reg_end);
        } else if (reg_tid != tid0) { // diff chr: set next_reg_beg, n_last_chunk_reads, last_chunk_read_i
            // reg_end = -1; // next chunk is diff chr, so set current reg_end as -1
            next_active_reg_beg = beg0;
            _n_last_chunk_reads = 0;
            pl->last_chunk_read_i[_n_last_chunk_reads++] = chunk->n_reads; // only keep one read for the next chunk
            chunk->is_ovlp[chunk->n_reads] = 0; // no overlap with the current chunk
            chunk->is_skipped[chunk->n_reads] = 1; // skip this read
            ++chunk->n_reads; // add the current read to the current chunk
            break;
        } else { // same chromosome
            // reg_beg/reg_end is already set
            assert(active_reg_beg != -1); assert(active_reg_end != -1);
            // if (reg_end == -1) reg_end = beg0 + max_reg_len_per_chunk - 1;
            // if (end0 >= active_reg_end) {
            if (end0 >= next_active_reg_beg) {
                chunk->is_ovlp[chunk->n_reads] = 1; // overlap between current and next chunk
                if (_n_last_chunk_reads == m_last_chunk_reads) {
                    m_last_chunk_reads *= 2;
                    pl->last_chunk_read_i = (int*)realloc(pl->last_chunk_read_i, m_last_chunk_reads * sizeof(int));
                }
                pl->last_chunk_read_i[_n_last_chunk_reads++] = chunk->n_reads;
            }
            if (beg0 > active_reg_end) { // out of the current region, first read in the next chunk
                // fprintf(stderr, "beg0: %ld %s, active_reg_end: %ld\n", beg0, bam_get_qname(chunk->reads[chunk->n_reads]), active_reg_end);
                chunk->is_ovlp[chunk->n_reads] = 0; // no overlap with the current chunk
                chunk->is_skipped[chunk->n_reads] = 1; // skip this read
                ++chunk->n_reads;
                break;
            }
        }
        // check eqx cigar or MD tag
        if (has_eqx_cigar == 0) has_eqx_cigar = has_equal_X_in_bam_cigar(chunk->reads[chunk->n_reads]);
        if (has_MD == 0) has_MD = has_MD_in_bam(chunk->reads[chunk->n_reads]);
        ++chunk->n_reads;
    }
    if (chunk->n_reads > pl->n_last_chunk_reads || _n_last_chunk_reads > 0) {
        chunk->tid = reg_tid; chunk->tname = header->target_name[reg_tid];
        get_bam_chunk_reg_cr(pl->ref_reg_seq->reg_cr, chunk, active_reg_beg, active_reg_end);
        get_bam_chunk_reg_ref_seq(pl->fai, pl->ref_reg_seq, chunk);
        chunk->bam_has_eqx_cigar = has_eqx_cigar; chunk->bam_has_md_tag = has_MD;
        // for the next chunk
        pl->cur_active_reg_beg = next_active_reg_beg;
        pl->n_last_chunk_reads = _n_last_chunk_reads;
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "CHUNK: tname: %s, tid: %d, beg: %" PRId64 ", end: %" PRId64 ", n_reads: %d, n_ovlp: %d, next_act_beg: %ld\n", chunk->tname, reg_tid, chunk->active_reg_beg, chunk->active_reg_end, chunk->n_reads, _n_last_chunk_reads, next_active_reg_beg);
        }
        // print reads in current chunk
        // if (LONGCALLD_VERBOSE >= 2) {
        //     for (int i = 0; i < chunk->n_reads; i++) {
        //         if (chunk->is_skipped[i] == 1) {
        //             fprintf(stderr, "Cur-CHUNK-SKIP: %s %ld-%ld\n", bam_get_qname(chunk->reads[i]), chunk->reads[i]->core.pos+1, bam_endpos(chunk->reads[i]));
        //         } else fprintf(stderr, "Cur-CHUNK-READ: %s %ld-%ld\n", bam_get_qname(chunk->reads[i]), chunk->reads[i]->core.pos+1, bam_endpos(chunk->reads[i]));
        //     }
        // }
    } else {
        pl->n_last_chunk_reads = 0;
        chunk->n_reads = 0;
        bam_chunk_free(chunk);
    }
    return r; // chunk->n_reads;
}

char *extract_sample_name_from_bam_header(bam_hdr_t *header) {
    int n_rg = sam_hdr_count_lines(header, "RG");
    for (int i = 0; i < n_rg; ++i) {
        kstring_t ks = KS_INITIALIZE;
        if (sam_hdr_find_tag_pos(header, "RG", i, "SM", &ks) == 0) {
            if (i < n_rg-1) _err_warning("Multiple RG/SM found in the BAM header, using the first one: %s\n", ks.s);
            else _err_info("Sample name extracted from the BAM header: %s\n", ks.s);
            return ks.s;
        }
    }
    return NULL;
}
