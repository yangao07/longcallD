#include <ctype.h>
#include "bam_utils.h"
#include "utils.h"
#include "call_var.h"

extern int LONGCALLD_VERBOSE;
read_snp_profile_t *init_read_snp_profile(int n_reads, int n_total_snps) {
    read_snp_profile_t *p = (read_snp_profile_t*)malloc(n_reads * sizeof(read_snp_profile_t));
    for (int i = 0; i < n_reads; ++i) {
        p[i].read_id = i;
        p[i].start_snp_idx = -1; p[i].end_snp_idx = -2;
        p[i].snp_is_used = (uint8_t*)calloc(n_total_snps, sizeof(uint8_t));
        p[i].snp_bases = (uint8_t*)malloc(n_total_snps * sizeof(uint8_t));
        p[i].snp_qual = (uint8_t*)malloc(n_total_snps * sizeof(uint8_t));
    }
    return p;
}

void free_read_snp_profile(read_snp_profile_t *p, int n_reads) {
    for (int i = 0; i < n_reads; ++i) {
        if (p[i].snp_is_used) free(p[i].snp_is_used);
        if (p[i].snp_bases) free(p[i].snp_bases);
        if (p[i].snp_qual) free(p[i].snp_qual);
    }
    free(p);
}

int has_equal_X_in_bam_cigar(bam1_t *read) {
    // Get the CIGAR string for the read
    const uint32_t *cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    // Check if '=' or 'X' is in the CIGAR operations
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            return 1;
        } else if (op == BAM_CMATCH) {
            return 0;
        }
    }
    return 2;
}

int has_MD_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "MD");
    return s != NULL;
}

void print_digar(digar_t *digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < digar->n_digar; ++i) {
        fprintf(fp, "%ld\t%c\t%d\t%d\t%c\t%d\n", digar->digars[i].pos, BAM_CIGAR_STR[digar->digars[i].type], digar->digars[i].len, digar->digars[i].qi, LONGCALLD_BAM_BASE_STR[digar->digars[i].base], digar->digars[i].is_low_qual);
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
    free(q->is_dense); free(q->pos); free(q->counts); 
    free(q->lens); 
    free(q);
}

// pos/len: region pos and length
// count: number of sub/gaps in the region
//        equal: 0
//        sub: 1 (for each sub)
//        ins/del: 1 (entire gap, ignore the length)
void push_xid_queue(xid_queue_t *q, hts_pos_t pos, int len, int count) {
    q->pos[++q->rear] = pos; // left-most pos of the region
    q->lens[q->rear] = len; q->counts[q->rear] = count;
    q->count += count;

    while (q->pos[q->front]+q->lens[q->front]-1 <= pos - q->win) {
        q->count -= q->counts[q->front]; // q->lens[q->front];
        q->front++;
    }
    // fprintf(stderr, "pos: %ld, len: %d, count: %d, front: %d, rear: %d, count: %d\n", pos, len, count, q->front, q->rear, q->count);
    if (count > 0) {
        if (q->count > q->max_s) {
            for (int j = q->front; j <= q->rear; j++) {
                q->is_dense[j] = 1; // mark all dense regions
            }
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

int get_x_site_index(hts_pos_t *x_sites, int cur_site_i, int n_total_pos, hts_pos_t pos) {
    int site_i = -1;
    for (int i = cur_site_i; i < n_total_pos; i++) {
        if (x_sites[i] == pos) {
            site_i = i;
            break;
        } else if (x_sites[i] > pos) {
            break;
        }
    }
    if (site_i == -1) {
        if (LONGCALLD_VERBOSE >= 2)
        _err_func_printf("Mismatch position is ignored (no high-qual bases): %ld\n", pos);
    }
    return site_i;
}

// int check_snp_site_index(hts_pos_t *mis_pos, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos, int *hit) {
//     int mis_i = -1; *hit = 0;
//     for (int i = cur_mis_i; i < n_total_mis_pos; i++) {
//         if (mis_pos[i] >= pos) {
//             mis_i = i;
//             if (mis_pos[i] == pos) *hit = 1;
//             break;
//         }
//     }
//     if (mis_i == -1) mis_i = n_total_mis_pos;
//         // _err_error_exit("Mismatch position not found: %ld", pos);
//     return mis_i;
// }

int get_snp_start(cand_snp_t *snp_sites, int cur_site_i, int n_total_pos, hts_pos_t start) {
    int i;
    for (i = cur_site_i; i < n_total_pos; ++i) {
        if (snp_sites[i].pos >= start) return i;
    }
    return i;
}

int get_x_site_start(hts_pos_t *pos, int cur_site_i, int n_total_pos, hts_pos_t start) {
    int i;
    for (i = cur_site_i; i < n_total_pos; ++i) {
        if (pos[i] >= start) return i;
    }
    return i;
}

// for high-qual sites: update both depth and base
// for low-qual sites: only update depth
void merge_site1(cand_snp_t *snp_site, uint8_t base, uint8_t is_low_qual) {
    if (is_low_qual) {
        snp_site->n_low_depth++;
        return;
    }
    snp_site->n_depth++;
    int base_i = snp_site->n_uniq_bases, exist=0;
    for (int i = 0; i < snp_site->n_uniq_bases; ++i) {
        if (snp_site->bases[i] == base) {
            exist = 1; base_i = i;
        }
    }
    if (exist == 0) {
        snp_site->bases[base_i] = base;
        snp_site->base_covs[base_i] = 1;
        snp_site->base_to_i[base] = base_i;
        snp_site->n_uniq_bases += 1;
    } else {
        snp_site->base_covs[base_i] += 1;
    }
}

int merge_diff_site(cand_snp_t *cand_snps, hts_pos_t *x_sites, int cur_site_i, int n_total_x_sites, hts_pos_t pos, uint8_t base, uint8_t is_low_qual) {
    int new_site_i = get_x_site_index(x_sites, cur_site_i, n_total_x_sites, pos);
    if (new_site_i == -1) return cur_site_i;
    merge_site1(cand_snps+new_site_i, base, is_low_qual);
    return new_site_i;
}

int merge_equal_sites(cand_snp_t *snp_sites, hts_pos_t *x_sites, int cur_site_i, int n_total_pos, hts_pos_t pos_start, hts_pos_t pos_end, int qi, const uint8_t *bam_seq, uint8_t is_low_qual) {
    int new_site_i = cur_site_i;
    for (int i = cur_site_i; i < n_total_pos; i++) {
        if (x_sites[i] >= pos_start && x_sites[i] < pos_end) {
            merge_site1(snp_sites+i, LONGCALLD_BAM_REF_BASE_IDX, is_low_qual);
            if (snp_sites[i].ref_base == -1)
                snp_sites[i].ref_base = seq_nt16_int[bam_seqi(bam_seq, qi + snp_sites[i].pos - pos_start)];
            // printf("Merge ref site: %ld (%d): %d\n", snp_sites[i].pos, i, snp_sites[i].n_depth);
            new_site_i = i+1;
        } else if (x_sites[i] >= pos_end) {
            break;
        }
    }
    return new_site_i;
}

int merge_del_sites(cand_snp_t *snp_sites, hts_pos_t *x_sites, int cur_site_i, int n_total_pos, hts_pos_t pos_start, hts_pos_t pos_end, int qi) {
    int new_site_i = cur_site_i;
    for (int i = cur_site_i; i < n_total_pos; i++) {
        if (x_sites[i] >= pos_start && x_sites[i] < pos_end) {
            merge_site1(snp_sites+i, LONGCALLD_BAM_DEL_BASE_IDX, 1);
            // printf("Merge ref site: %ld (%d): %d\n", snp_sites[i].pos, i, snp_sites[i].n_depth);
            new_site_i = i+1;
        } else if (x_sites[i] >= pos_end) {
            break;
        }
    }
    return new_site_i;
}


int update_cand_snps_from_digar(digar_t *digar, bam1_t *read, int n_x_sites, hts_pos_t *x_sites, int start_i, cand_snp_t *cand_snps) {
    int cur_start_i = -1, site_i = -1;
    const uint8_t *bam_seq = bam_get_seq(read);

    for (int i = 0; i < digar->n_digar; ++i) {
        if (cur_start_i == -1) site_i = cur_start_i = get_x_site_start(x_sites, start_i, n_x_sites, digar->digars[i].pos);
        // if (digar->digars[i].is_low_qual) continue;
        if (digar->digars[i].type == BAM_CDIFF) {
            site_i = merge_diff_site(cand_snps, x_sites, site_i, n_x_sites, digar->digars[i].pos, digar->digars[i].base, digar->digars[i].is_low_qual);
        } else if (digar->digars[i].type == BAM_CEQUAL) {
            site_i = merge_equal_sites(cand_snps, x_sites, site_i, n_x_sites, digar->digars[i].pos, digar->digars[i].pos+digar->digars[i].len, digar->digars[i].qi, bam_seq, digar->digars[i].is_low_qual);
        } else if (digar->digars[i].type == BAM_CDEL) {
            site_i = merge_del_sites(cand_snps, x_sites, site_i, n_x_sites, digar->digars[i].pos, digar->digars[i].pos+digar->digars[i].len, digar->digars[i].qi);
        }
   }
    return cur_start_i;
}

int update_read_snp_profile_from_D_sites(cand_snp_t *snp_sites, int cur_site_i, int n_total_pos, hts_pos_t pos_start, hts_pos_t pos_end, 
                                         uint8_t qual, read_snp_profile_t *read_snp_profile) {
    int new_site_i = cur_site_i;
    for (int i = cur_site_i; i < n_total_pos; ++i) {
        if (snp_sites[i].pos >= pos_start && snp_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
            /* 6(D): ref allele */
            int snp_idx = i - read_snp_profile->start_snp_idx;
            read_snp_profile->snp_bases[snp_idx] = LONGCALLD_BAM_DEL_BASE_IDX; 
            read_snp_profile->snp_qual[snp_idx] = qual;
            new_site_i = i+1;
            read_snp_profile->snp_is_used[snp_idx] = 0; // D sites: not used
        } else if (snp_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_site_i;
}

int udpate_read_snp_profile_from_equal_sites(cand_snp_t *snp_sites, int cur_site_i, int n_total_pos, hts_pos_t pos_start, hts_pos_t pos_end, 
                                             const uint8_t *qual, read_snp_profile_t *read_snp_profile, uint8_t is_low_qual) {
    int new_site_i = cur_site_i;
    for (int i = cur_site_i; i < n_total_pos; ++i) {
        if (snp_sites[i].pos >= pos_start && snp_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
            int snp_idx = i - read_snp_profile->start_snp_idx;
             /* 5(.): ref allele */
            read_snp_profile->snp_bases[snp_idx] = LONGCALLD_BAM_REF_BASE_IDX; 
            read_snp_profile->snp_qual[snp_idx] = qual[snp_sites[i].pos-pos_start];
            new_site_i = i+1;
            if (!is_low_qual) read_snp_profile->snp_is_used[snp_idx] = 1;
        } else if (snp_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_site_i;
}

int update_read_snp_profile_from_diff_site(cand_snp_t *snp_sites, int cur_site_i, int n_total_pos, hts_pos_t pos_start, hts_pos_t pos_end, 
                                            int qi, const uint8_t *bam_seq, const uint8_t *qual, read_snp_profile_t *read_snp_profile, uint8_t is_low_qual) {
    int new_site_i = cur_site_i;
    for (int i = cur_site_i; i < n_total_pos; ++i) {
        if (snp_sites[i].pos >= pos_start && snp_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
            int snp_idx = i - read_snp_profile->start_snp_idx;
            read_snp_profile->snp_bases[snp_idx] = seq_nt16_int[bam_seqi(bam_seq, qi+snp_sites[i].pos-pos_start)];
            read_snp_profile->snp_qual[snp_idx] = qual[qi+snp_sites[i].pos-pos_start];
            new_site_i = i+1;
            if (!is_low_qual) read_snp_profile->snp_is_used[snp_idx] = 1;
        } else if (snp_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_site_i;
}

int update_read_snp_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_snps, cand_snp_t *cand_snps, int start_snp_i, read_snp_profile_t *read_snp_profile) {
    int cur_start_i = -1, snp_i = -1;
    const uint8_t *qual = bam_get_qual(read), *bam_seq = bam_get_seq(read);
    for (int i = 0; i < digar->n_digar; ++i) {
        if (cur_start_i == -1) snp_i = cur_start_i = get_snp_start(cand_snps, start_snp_i, n_cand_snps, digar->digars[i].pos);
        // if (digar->digars[i].is_low_qual) continue;
        if (digar->digars[i].type == BAM_CDIFF) {
            snp_i = update_read_snp_profile_from_diff_site(cand_snps, snp_i, n_cand_snps, digar->digars[i].pos, digar->digars[i].pos+1, digar->digars[i].qi, bam_seq, qual, read_snp_profile, digar->digars[i].is_low_qual);
        } else if (digar->digars[i].type == BAM_CEQUAL) {
            snp_i = udpate_read_snp_profile_from_equal_sites(cand_snps, snp_i, n_cand_snps, digar->digars[i].pos, digar->digars[i].pos+digar->digars[i].len, qual+digar->digars[i].qi, read_snp_profile, digar->digars[i].is_low_qual);
        } else if (digar->digars[i].type == BAM_CDEL) {
            snp_i = update_read_snp_profile_from_D_sites(cand_snps, snp_i, n_cand_snps, digar->digars[i].pos, digar->digars[i].pos+digar->digars[i].len, digar->digars[i].qual, read_snp_profile);
        }
    }
    return cur_start_i;
}

void set_digar(digar1_t *digar, hts_pos_t pos, int type, int len, int qi, uint8_t base, uint8_t qual) {
    // fprintf(stderr, "Set digar: %ld, %d%c, %d\n", pos, len, BAM_CIGAR_STR[type], qi);
    digar->pos = pos; digar->type = type; digar->len = len;
    digar->qi = qi; digar->base = base; digar->qual = qual;
}

void push_digar0(digar_t *digar, hts_pos_t pos, int type, int len, int qi, uint8_t base, uint8_t qual, int is_low_qual) {
    _uni_realloc(digar->digars, digar->n_digar, digar->m_digar, digar1_t);
    digar->digars[digar->n_digar].pos = pos; digar->digars[digar->n_digar].type = type; digar->digars[digar->n_digar].len = len;
    digar->digars[digar->n_digar].qi = qi; digar->digars[digar->n_digar].base = base; digar->digars[digar->n_digar].qual = qual;
    digar->digars[digar->n_digar].is_low_qual = is_low_qual;
    digar->n_digar++;
}

void push_digar1(digar_t *digar, digar1_t d) {
    _uni_realloc(digar->digars, digar->n_digar, digar->m_digar, digar1_t);
    digar->digars[digar->n_digar].pos = d.pos; digar->digars[digar->n_digar].type = d.type; digar->digars[digar->n_digar].len = d.len;
    digar->digars[digar->n_digar].qi = d.qi; digar->digars[digar->n_digar].base = d.base; digar->digars[digar->n_digar].qual = d.qual;
    digar->digars[digar->n_digar].is_low_qual = d.is_low_qual;
    digar->n_digar++;
}

void print_digar1(digar1_t *digars, int n_digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < n_digar; ++i) {
        fprintf(fp, "%ld\t%c\t%d\t%d\t%c\t%d\n", digars[i].pos, BAM_CIGAR_STR[digars[i].type], digars[i].len, digars[i].qi, LONGCALLD_BAM_BASE_STR[digars[i].base], digars[i].is_low_qual);
    }
}

void assign_down_flanking_len(digar1_t *digars, int i, int len, int *down_low_qual_len) {
    int j = i, remain_low_qual_len = len;
    while (j >= 0 && remain_low_qual_len > 0) {
        if (digars[j].len < remain_low_qual_len) {
            down_low_qual_len[j] = digars[j].len;
            remain_low_qual_len -= digars[j].len;
        } else {
            down_low_qual_len[j] = MAX_OF_TWO(remain_low_qual_len, down_low_qual_len[j]);
            remain_low_qual_len = 0;
        }
        j--;
    }
}

void assign_up_flanking_len(digar1_t *digars, int i, int n, int len, int *up_low_qual_len) {
    int j = i, remain_low_qual_len = len;
    while (j < n && remain_low_qual_len > 0) {
        if (digars[j].len < remain_low_qual_len) {
            up_low_qual_len[j] = digars[j].len;
            remain_low_qual_len -= digars[j].len;
        } else {
            up_low_qual_len[j] = MAX_OF_TWO(remain_low_qual_len, up_low_qual_len[j]);
            remain_low_qual_len = 0;
        }
        j++;
    }
}

// for each dense-region and standalone indel
// 1. dense-region: extend dense_flank_reg_size
// 2. indel: extend indel_flank_win_size
// 3. split equal regions based on flanking size
void collect_noisy_flanking_wins(digar1_t *_digars, int _n_digar, int dense_flank_reg_size, int indel_flank_win_size, int *up_low_qual_len, int *down_low_qual_len) {
    for (int i = 0; i < _n_digar; ++i) {
        // if (_digars[i].is_low_qual)
            // fprintf(stderr, "Low: %c, pos: %ld, end: %ld, len: %d\n", BAM_CIGAR_STR[_digars[i].type], _digars[i].pos, _digars[i].pos+_digars[i].len-1, _digars[i].len);
        if (_digars[i].type == BAM_CEQUAL) continue; 
        int _up_low_qual_len = 0, _down_low_qual_len = 0;
        if (i > 0 && _digars[i-1].is_low_qual == 0) { // upstream of a dense region/indel
            if (_digars[i].is_low_qual) {
                _down_low_qual_len = MAX_OF_TWO(dense_flank_reg_size, _down_low_qual_len);
            } else if (_digars[i].type == BAM_CDEL || _digars[i].type == BAM_CINS) {
                _down_low_qual_len = MAX_OF_TWO(indel_flank_win_size, _down_low_qual_len);
            }
            // assign _down_low_qual_len to i'th and previous regions
            assign_down_flanking_len(_digars, i-1, _down_low_qual_len, down_low_qual_len);
        }
        if (i < _n_digar-1 && _digars[i+1].is_low_qual == 0) { // downstream of a dense region/indel
            if (_digars[i].is_low_qual) {
                _up_low_qual_len = MAX_OF_TWO(dense_flank_reg_size, _up_low_qual_len);
            } else if (_digars[i].type == BAM_CDEL || _digars[i].type == BAM_CINS) {
                _up_low_qual_len = MAX_OF_TWO(indel_flank_win_size, _up_low_qual_len);
            }
            // assign _up_low_qual_len to i'th and following regions
            assign_up_flanking_len(_digars, i+1, _n_digar, _up_low_qual_len, up_low_qual_len);
        }
    }
    // print up and down low-qual lengths
    // for (int i = 0; i < _n_digar; ++i) {
    //     if (up_low_qual_len[i] > 0 || down_low_qual_len[i] > 0) {
    //         fprintf(stderr, "[%c] pos: %ld, end: %ld, up: %d, down: %d\n", BAM_CIGAR_STR[_digars[i].type],  _digars[i].pos, _digars[i].pos+_digars[i].len-1, up_low_qual_len[i], down_low_qual_len[i]);
    //     }
    // }
}

void post_update_digar(digar1_t *_digars, int _n_digar, const struct call_var_opt_t *opt, digar_t *digar) {
    int dense_flank_reg_size = opt->dens_reg_slide_win / (opt->dens_reg_max_sites+1);
    int indel_flank_win_size = opt->indel_flank_win_size;
    int *up_low_qual_len = (int*)calloc(_n_digar, sizeof(int));
    int *down_low_qual_len = (int*)calloc(_n_digar, sizeof(int));
    collect_noisy_flanking_wins(_digars, _n_digar, dense_flank_reg_size, indel_flank_win_size, up_low_qual_len, down_low_qual_len);
    // split equal regions based on up/down_low_qual_len
    for (int i = 0; i < _n_digar; ++i) {
        // is_low_qual: push
        if (_digars[i].is_low_qual || (up_low_qual_len[i] == 0 && down_low_qual_len[i] == 0)) {
            push_digar1(digar, _digars[i]);
        } else {
            if (_digars[i].type == BAM_CEQUAL) { // EQUAL: split if up/down_low_qual_len > 0
                if (up_low_qual_len[i] + down_low_qual_len[i] >= _digars[i].len) {
                    push_digar1(digar, _digars[i]); digar->digars[digar->n_digar-1].is_low_qual = 1;
                } else {
                    if (up_low_qual_len[i] > 0)
                        push_digar1(digar, (digar1_t){_digars[i].pos, BAM_CEQUAL, up_low_qual_len[i], _digars[i].qi, LONGCALLD_BAM_REF_BASE_IDX, 0, 1});
                    if (_digars[i].len - up_low_qual_len[i] - down_low_qual_len[i] > 0)
                        push_digar1(digar, (digar1_t){_digars[i].pos+up_low_qual_len[i], BAM_CEQUAL, _digars[i].len - up_low_qual_len[i] - down_low_qual_len[i], _digars[i].qi+up_low_qual_len[i], LONGCALLD_BAM_REF_BASE_IDX, 0, 0});
                    if (down_low_qual_len[i] > 0)
                        push_digar1(digar, (digar1_t){_digars[i].pos+_digars[i].len-down_low_qual_len[i], BAM_CEQUAL, down_low_qual_len[i], _digars[i].qi+_digars[i].len-down_low_qual_len[i], LONGCALLD_BAM_REF_BASE_IDX, 0, 1});
                }
            } else { // X/I/D: set as low_qual if up/down_low_qual_len > 0
                push_digar1(digar, _digars[i]); digar->digars[digar->n_digar-1].is_low_qual = 1;
            }
        }
    }
    free(up_low_qual_len); free(down_low_qual_len);
}

void collect_digar_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    const uint8_t *qual = bam_get_qual(read), *bam_seq = bam_get_seq(read);
    int min_bq = opt->min_bq, max_s = opt->dens_reg_max_sites, win = opt->dens_reg_slide_win;
    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar);
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);

    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                if (qual[qi] >= min_bq) {
                    // _uni_realloc(*read_xs, n_xid, m_xid, xid_t);
                    _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                    set_digar(_digars+_n_digar, pos, op, 1, qi, seq_nt16_int[bam_seqi(bam_seq, qi)], qual[qi]);
                    push_xid_queue(q, pos, 1, 1); _n_digar++;
                }
                pos++; qi++;
            }
        } else if (op == BAM_CEQUAL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, LONGCALLD_BAM_REF_BASE_IDX, 0); // unset XXX
            _n_digar++; push_xid_queue(q, pos, 0, 0);
            pos += len; qi += len;
        } else if (op == BAM_CDEL) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            int d_qual; if (qi > 0) d_qual = (qual[qi] + qual[qi-1]) / 2; else d_qual = qual[qi];
            set_digar(_digars+_n_digar, pos, op, len, qi, LONGCALLD_BAM_DEL_BASE_IDX, d_qual);
            _n_digar++; push_xid_queue(q, pos, len, 1);
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, op, len, qi, LONGCALLD_BAM_INS_BASE_IDX, 0); // insertion, unset XXX
            _n_digar++; push_xid_queue(q, pos, 0, 1);
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            qi += len;
        } else if (op == BAM_CREF_SKIP) { // XXX no action for N op
            pos += len;
        } else if (op == BAM_CMATCH) {
            _err_error_exit("CIGAR operation 'M' is not expected in EQX CIGAR: %s", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) _digars[i].is_low_qual = q->is_dense[i];
    post_update_digar(_digars, _n_digar, opt, digar);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "DIGAR: %s\n", bam_get_qname(read));
        print_digar(digar, stderr);
    }
    free(_digars); free_xid_queue(q);
}

void collect_digar_from_MD_tag(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar) {
    uint8_t *s = bam_aux_get(read, "MD");
    if (s == NULL) { _err_error_exit("MD tag not found in the BAM file: %s", bam_get_qname(read)); }
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int min_bq = opt->min_bq, max_s = opt->dens_reg_max_sites, win = opt->dens_reg_slide_win;
    const uint8_t *qual = bam_get_qual(read), *bam_seq = bam_get_seq(read);

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    char *md = bam_aux2Z(s); int md_i = 0;
    // printf("MD: %s\n", md);
    int rlen = bam_cigar2rlen(n_cigar, cigar);
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);

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
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, LONGCALLD_BAM_REF_BASE_IDX, 0); // qual unset XXX
                        _n_digar++; push_xid_queue(q, pos, 0, 0);
                        pos += eq_len; qi += eq_len;
                        last_eq_len -= m_len; m_len = 0;
                    } else { // last_eq_len < m_len
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CEQUAL, last_eq_len, qi, LONGCALLD_BAM_REF_BASE_IDX, 0); // qual unset XXX
                        _n_digar++; push_xid_queue(q, pos, 0, 0);
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
                    set_digar(_digars+_n_digar, pos, BAM_CEQUAL, eq_len, qi, LONGCALLD_BAM_REF_BASE_IDX, 0); // qual unset XXX
                    _n_digar++; push_xid_queue(q, pos, 0, 0);
                    pos += eq_len; qi += eq_len;
                    m_len -= eq_len; md_i = 0; // start of the next run
                } else if (isalpha(md[md_i])) { // X
                    if (qual[qi] >= min_bq) {
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, seq_nt16_int[bam_seqi(bam_seq, qi)], qual[qi]);
                        _n_digar++; push_xid_queue(q, pos, 1, 1);
                    }
                    pos++; qi++; m_len -= 1; 
                    if (md[md_i+1] == '\0' || md[md_i+1] != '0') md_i++;
                    else md_i += 2; // skip the 0 after X
                } else _err_error_exit("MD and CIGAR do not match: %s (%s)", bam_get_qname(read), bam_aux2Z(s));
                if (m_len <= 0) break;
            }
        } else if (op == BAM_CDEL) {
            uint8_t d_qual; if (qi > 0) d_qual = (qual[qi] + qual[qi-1]) / 2; else d_qual = qual[qi];
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, LONGCALLD_BAM_DEL_BASE_IDX, d_qual);
            _n_digar++; push_xid_queue(q, pos, len, 1);
            pos += len;
            // MD
            md_i++;
            while (md[md_i] && isalpha(md[md_i])) md_i++;
            if (md[md_i] == '0') md_i++; // skip 0 after D
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, LONGCALLD_BAM_INS_BASE_IDX, 0); // insertion, unset XXX
            _n_digar++; push_xid_queue(q, pos, 0, 1);
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) _digars[i].is_low_qual = q->is_dense[i];
    post_update_digar(_digars, _n_digar, opt, digar);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "DIGAR: %s\n", bam_get_qname(read));
        print_digar(digar, stderr);
    }
    free(_digars); free_xid_queue(q);
}

void collect_digar_from_ref_seq(bam1_t *read, const call_var_opt_t *opt, kstring_t *ref_seq, digar_t *digar) {
    hts_pos_t pos = read->core.pos+1, qi = 0;
    const uint8_t *qual = bam_get_qual(read), *bam_seq = bam_get_seq(read);
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    int min_bq = opt->min_bq, max_s = opt->dens_reg_max_sites, win = opt->dens_reg_slide_win;

    digar->n_digar = 0; digar->m_digar = 2 * n_cigar; digar->digars = (digar1_t*)malloc(n_cigar * 2 * sizeof(digar1_t));
    int _n_digar = 0, _m_digar = 2 * n_cigar; digar1_t *_digars = (digar1_t*)malloc(_m_digar * sizeof(digar1_t));
    int rlen = bam_cigar2rlen(n_cigar, cigar);
    xid_queue_t *q = init_xid_queue(rlen, max_s, win);

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH ) {
            int eq_len = 0;
            for (int j = 0; j < len; ++j) {
                if (qual[qi] >= min_bq) {
                    // Get the reference base
                    if (pos > ref_seq->l) _err_error_exit("Read exceed reference sequence length: %s", bam_get_qname(read));
                    char ref_base = ref_seq->s[pos-1];
                    // Get the read base
                    char read_base = seq_nt16_str[bam_seqi(bam_seq, qi)];
                    if (ref_base != read_base) {
                        if (eq_len > 0) {
                            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                            set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, LONGCALLD_BAM_REF_BASE_IDX, 0);
                            _n_digar++; push_xid_queue(q, pos-eq_len, 0, 0);
                            eq_len = 0;
                        }
                        _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                        set_digar(_digars+_n_digar, pos, BAM_CDIFF, 1, qi, seq_nt16_int[bam_seqi(bam_seq, qi)], qual[qi]);
                        _n_digar++; push_xid_queue(q, pos, 1, 1);
                    } else eq_len++;
                }
                pos++; qi++;
            }
            if (eq_len > 0) {
                _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
                set_digar(_digars+_n_digar, pos-eq_len, BAM_CEQUAL, eq_len, qi-eq_len, LONGCALLD_BAM_REF_BASE_IDX, 0);
                _n_digar++; push_xid_queue(q, pos-eq_len, 0, 0);
                eq_len = 0;
            }
        } else if (op == BAM_CDEL) {
            uint8_t d_qual; if (qi > 0) d_qual = (qual[qi] + qual[qi-1]) / 2; else d_qual = qual[qi];
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CDEL, len, qi, LONGCALLD_BAM_DEL_BASE_IDX, d_qual);
            _n_digar++; push_xid_queue(q, pos, len, 1);
            pos += len;
        } else if (op == BAM_CINS) {
            _uni_realloc(_digars, _n_digar, _m_digar, digar1_t);
            set_digar(_digars+_n_digar, pos, BAM_CINS, len, qi, LONGCALLD_BAM_INS_BASE_IDX, 0); // insertion, unset XXX
            _n_digar++; push_xid_queue(q, pos, 0, 1);
            qi += len;
        } else if (op == BAM_CSOFT_CLIP) {
            qi += len;
        } else if (op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            _err_error_exit("CIGAR operation '=/X' is not expected: %s", bam_get_qname(read));
        }
    }
    for (int i = 0; i < _n_digar; ++i) _digars[i].is_low_qual = q->is_dense[i];
    // if (LONGCALLD_VERBOSE >= 2) {
        // fprintf(stderr, "RAW DIGAR: %s\n", bam_get_qname(read));
        // print_digar1(_digars, _n_digar, stderr);
    // }
    post_update_digar(_digars, _n_digar, opt, digar);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "DIGAR: %s\n", bam_get_qname(read));
        print_digar(digar, stderr);
    }
    free(_digars); free_xid_queue(q);
}

void copy_bam_chunk(bam_chunk_t *from_chunk, bam_chunk_t *to_chunk, int *ovlp_read_i, int n_reads) {
    for (int i = 0; i < n_reads; i++) {
        int read_i = ovlp_read_i[i];
        // if (bam_copy1(to_chunk->reads[to_chunk->n_reads + i], from_chunk->reads[read_i]) == NULL)
            // _err_error_exit("Failed to copy BAM record: %s (%d,%d)", bam_get_qname(from_chunk->reads[read_i]), read_i, n_reads);
        to_chunk->reads[to_chunk->n_reads + i] =  from_chunk->reads[read_i];
    }
    to_chunk->n_reads += n_reads;
}

int bam_read_chunk_init(bam_chunk_t *chunk, int n_reads, int *ovlp_read_i, int n_ovlp_reads) {
    if (n_ovlp_reads > n_reads) n_reads = n_ovlp_reads;

    chunk->reads = (bam1_t**)malloc(n_reads * sizeof(bam1_t*));
    chunk->digars = (digar_t*)malloc(n_reads * sizeof(digar_t));
    for (int i = 0; i < n_reads; i++) {
        if (i >= n_ovlp_reads) chunk->reads[i] = bam_init1();
        else chunk->reads[i] = NULL;
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
    }
    chunk->is_skipped = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->haps = (int*)calloc(n_reads, sizeof(int));
    chunk->n_reads = 0; chunk->m_reads = n_reads;

    if (n_ovlp_reads > 0) {
        copy_bam_chunk(chunk-1, chunk, ovlp_read_i, n_ovlp_reads);
        chunk->n_ovlp_reads = n_ovlp_reads;
    } 
    return 0;
}

int bam_read_chunk_realloc(bam_chunk_t *chunk) {
    int m_reads = chunk->m_reads * 2;
    chunk->reads = (bam1_t**)realloc(chunk->reads, m_reads * sizeof(bam1_t*));
    chunk->digars = (digar_t*)realloc(chunk->digars, m_reads * sizeof(digar_t));
    chunk->is_skipped = (uint8_t*)realloc(chunk->is_skipped, m_reads * sizeof(uint8_t));
    chunk->haps = (int*)realloc(chunk->haps, m_reads * sizeof(int));
    for (int i = chunk->m_reads; i < m_reads; i++) {
        chunk->reads[i] = bam_init1();
        chunk->digars[i].n_digar = chunk->digars[i].m_digar = 0;
        chunk->is_skipped[i] = 0; chunk->haps[i] = 0;
    }
    chunk->m_reads = m_reads;
    return 0;
}

void bam_read_chunk_free(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->m_reads; i++) {
        if (chunk->digars[i].m_digar > 0) free(chunk->digars[i].digars);
    }
    for (int i = chunk->n_ovlp_reads; i < chunk->m_reads; i++) {
        bam_destroy1(chunk->reads[i]);
    }
    free(chunk->reads); free(chunk->digars);
    free(chunk->is_skipped); free(chunk->haps);
}


// collect reads from a BAM file, region: [start, start+max_reg_len_per_chunk]
// *n_ovlp_reads: number of overlapping reads between previous chunk and current chunk, update it to "overlapping reads between current chunk and next chunk" after collecting reads
int collect_bam_chunk(samFile *in_bam, bam_hdr_t *header, hts_itr_t *iter, int use_iter, int max_reg_len_per_chunk, int **ovlp_read_i, int *n_ovlp_reads, bam_chunk_t *chunk) {
    int r, has_eqx_cigar=0, has_MD=0;
    int tid = -1; char *tname = NULL; 
    hts_pos_t beg=-1, _beg=-1, end=-1, _end = -1, reg_end=-1;

    bam_read_chunk_init(chunk, 4096, *ovlp_read_i, *n_ovlp_reads);
    if (*ovlp_read_i != NULL) free(*ovlp_read_i);

    int m_ovlp_reads = chunk->m_reads, _n_ovlp_reads = 0; // push to ovlp_read_i if end >= reg_end, break if beg > reg_end
    *ovlp_read_i = (int*)malloc(m_ovlp_reads * sizeof(int));
    while (1) {
        if (use_iter) r = sam_itr_next(in_bam, iter, chunk->reads[chunk->n_reads]);
        else r = sam_read1(in_bam, header, chunk->reads[chunk->n_reads]);
        if (r < 0) break;
        // skip unmapped/secondary alignments
        if (chunk->reads[chunk->n_reads]->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue; // BAM_FSUPPLEMENTARY
        if (tid == -1) { // first in the chunk
            tid = chunk->reads[chunk->n_reads]->core.tid;
            tname = header->target_name[tid];
            _beg = chunk->reads[chunk->n_reads]->core.pos + 1;
            reg_end = _beg + max_reg_len_per_chunk - 1;
        } else if (tid != chunk->reads[chunk->n_reads]->core.tid) break; // different chromosome
        else { // same chromosome
            beg = chunk->reads[chunk->n_reads]->core.pos + 1;
            end = beg + bam_cigar2rlen(chunk->reads[chunk->n_reads]->core.n_cigar, bam_get_cigar(chunk->reads[chunk->n_reads])) - 1;
            if (reg_end == -1) reg_end = beg + max_reg_len_per_chunk - 1;
            if (end >= reg_end) {
                (*ovlp_read_i)[_n_ovlp_reads] = chunk->n_reads;
                if (++_n_ovlp_reads >= m_ovlp_reads) {
                    m_ovlp_reads *= 2;
                    *ovlp_read_i = (int*)realloc(*ovlp_read_i, m_ovlp_reads * sizeof(int));
                }
            }
            if (beg > reg_end) break; // out of the region
            if (end > _end) _end = end; 
        }
        // fprintf(stderr, "qname: %s, pos: %ld\n", bam_get_qname(chunk->reads[chunk->n_reads]), beg);
        // check eqx cigar or MD tag
        if (has_eqx_cigar == 0) has_eqx_cigar = has_equal_X_in_bam_cigar(chunk->reads[chunk->n_reads]);
        if (has_MD == 0) has_MD = has_MD_in_bam(chunk->reads[chunk->n_reads]);
        if (++(chunk->n_reads) >= chunk->m_reads) bam_read_chunk_realloc(chunk);
    }

    if (chunk->n_reads > *n_ovlp_reads) {
        chunk->tid = tid; chunk->tname = tname;
        chunk->beg = _beg; chunk->end = _end;
        chunk->bam_has_eqx_cigar = has_eqx_cigar; chunk->bam_has_md_tag = has_MD;
        // fprintf(stderr, "CHUNK: eqx: %d, MD: %d, n_reads: %d, n_ovlp: %d (%d)\n", has_eqx_cigar, has_MD, chunk->n_reads, *n_ovlp_reads, _n_ovlp_reads);
        *n_ovlp_reads = _n_ovlp_reads;
    } else {
        *n_ovlp_reads = 0;
        chunk->n_reads = 0;
        bam_read_chunk_free(chunk);
    }
    return chunk->n_reads;
}