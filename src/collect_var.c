#include <stdio.h>
#include <stdlib.h>
#include "collect_var.h"
#include "call_var_main.h"
#include "utils.h"
#include "bam_utils.h"
#include "align.h"
#include "seq.h"
#include "assign_hap.h"
#include "vcf_utils.h"
#include "kmedoids.h"

extern int LONGCALLD_VERBOSE;
char test_read_name[1024] = "m84039_230928_213653_s3/23396787/ccs";

// only 1 variant is stored in each cand_var_t, i.e., SNP/INS/DEL
// each read with be 
//    0: ref ,or 1: alt, or
//   -1: other alt. mostly sequencing errors, should be considered as non-informative during haplotype assignment, i.e. -1!=-1
cand_var_t *init_cand_vars_based_on_sites(int n_var_sites, var_site_t *var_sites) {
    cand_var_t *cand_vars = (cand_var_t*)malloc(n_var_sites * sizeof(cand_var_t));
    for (int i = 0; i < n_var_sites; ++i) {
        // static information
        // from var_sites
        cand_vars[i].tid = var_sites[i].tid;
        cand_vars[i].pos = var_sites[i].pos; 
        cand_vars[i].phase_set = 0; // unset
        cand_vars[i].var_type = var_sites[i].var_type;
        cand_vars[i].ref_len = var_sites[i].ref_len; cand_vars[i].ref_base = 5; // unknown

        cand_vars[i].total_cov = 0; cand_vars[i].low_qual_cov = 0; 
        cand_vars[i].n_uniq_alles = 2; // ref/alt
        cand_vars[i].alle_covs = (int*)calloc(2, sizeof(int)); // ref/alt
        cand_vars[i].strand_to_alle_covs = (int**)malloc(2 * sizeof(int*));
        for (int j = 0; j < 2; ++j) {
            cand_vars[i].strand_to_alle_covs[j] = (int*)calloc(2, sizeof(int)); // ref/alt
        }
        cand_vars[i].alt_len = var_sites[i].alt_len; // cand_vars[i].alt_seq = var_sites[i].alt_seq;
        if (cand_vars[i].var_type == BAM_CINS || cand_vars[i].var_type == BAM_CDIFF) {
            cand_vars[i].alt_seq = (uint8_t*)malloc(cand_vars[i].alt_len * sizeof(uint8_t));
            for (int j = 0; j < cand_vars[i].alt_len; ++j) cand_vars[i].alt_seq[j] = var_sites[i].alt_seq[j];
        } else cand_vars[i].alt_seq = NULL;
        // dynamic information, allocate and update during haplotype assignment
        cand_vars[i].hap_to_alle_profile = NULL; cand_vars[i].hap_to_cons_alle = NULL;
    }
    return cand_vars;
}

void free_cand_vars1(cand_var_t *cand_vars) {
    if (cand_vars->alle_covs != NULL) free(cand_vars->alle_covs); 
    if (cand_vars->strand_to_alle_covs != NULL) {
        for (int j = 0; j < 2; ++j) {
            free(cand_vars->strand_to_alle_covs[j]);
        } free(cand_vars->strand_to_alle_covs);
    }
    if (cand_vars->alt_seq != NULL) free(cand_vars->alt_seq);
    if (cand_vars->hap_to_alle_profile != NULL) {
        for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
            free(cand_vars->hap_to_alle_profile[j]);
        } free(cand_vars->hap_to_alle_profile);
    }
    if (cand_vars->hap_to_cons_alle != NULL) free(cand_vars->hap_to_cons_alle);
}


void free_cand_vars(cand_var_t *cand_vars, int m) {
    for (int i = 0; i < m; ++i) free_cand_vars1(cand_vars+i);
    free(cand_vars);
}

void free_noisy_cand_vars1(cand_var_t *cand_vars) {
    if (cand_vars->alle_covs != NULL) free(cand_vars->alle_covs); 
    if (cand_vars->strand_to_alle_covs != NULL) {
        for (int j = 0; j < 2; ++j) {
            free(cand_vars->strand_to_alle_covs[j]);
        } free(cand_vars->strand_to_alle_covs);
    }
    if (cand_vars->alt_seq != NULL) free(cand_vars->alt_seq);

    if (cand_vars->hap_to_alle_profile != NULL) {
        for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
            free(cand_vars->hap_to_alle_profile[j]);
        } free(cand_vars->hap_to_alle_profile);
    }
    if (cand_vars->hap_to_cons_alle != NULL) free(cand_vars->hap_to_cons_alle);
}

void free_noisy_cand_vars(cand_var_t *cand_vars, int m) {
    for (int i = 0; i < m; ++i) free_noisy_cand_vars1(cand_vars+i);
    free(cand_vars);
}

int collect_cand_vars(bam_chunk_t *chunk, int n_var_sites, var_site_t *var_sites) {
    chunk->cand_vars = init_cand_vars_based_on_sites(n_var_sites, var_sites);
    cand_var_t *cand_vars = chunk->cand_vars;
    // 2nd pass: update snp_sites, calculate the depth and allele frequency of each site
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        bam1_t *read = chunk->reads[i];
        if (chunk->is_skipped[i]) continue;
        start_var_i = update_cand_vars_from_digar(chunk->digars+i, read, n_var_sites, var_sites, start_var_i, cand_vars);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_var_sites; ++i) {
            fprintf(stderr, "CandVar: %" PRId64 "\t", var_sites[i].pos);
            fprintf(stderr, "Type: %c\t", BAM_CIGAR_STR[var_sites[i].var_type]);
            fprintf(stderr, "RefLen: %d\t", var_sites[i].ref_len);
            fprintf(stderr, "Depth: %d\t", cand_vars[i].total_cov);
            for (int k = 0; k < cand_vars[i].alt_len; ++k)
                fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seq[k]]);
            fprintf(stderr, ": ");
            fprintf(stderr, "%d\t", cand_vars[i].alle_covs[1]);
            fprintf(stderr, "Low-Depth: %d\n", cand_vars[i].low_qual_cov);
        }
    }
    return 0;
}

int var_is_strand_bias(cand_var_t *var) {
    // strand bias: 1) >=3 vs 0, 2) >= 3 folds
    int min_fold = 3;
    if (var->alle_covs[1] < min_fold) return 0;
    else if (var->strand_to_alle_covs[0][1] >= 2 && var->strand_to_alle_covs[1][1] >= 2) return 0;
    else if (var->strand_to_alle_covs[0][1] == 0 && var->strand_to_alle_covs[1][1] >= min_fold) return 1;
    else if (var->strand_to_alle_covs[1][1] == 0 && var->strand_to_alle_covs[0][1] >= min_fold) return 1;
    else if (var->strand_to_alle_covs[0][1] > 0 && var->strand_to_alle_covs[1][1] > 0) {
        if (var->strand_to_alle_covs[0][1] >= min_fold * var->strand_to_alle_covs[1][1]) return 1;
        else if (var->strand_to_alle_covs[1][1] >= min_fold * var->strand_to_alle_covs[0][1]) return 1;
    }
    return 0;
}

// XXX check the sequence around the variant site, not just left/right side
// filter: 
//   1) usable X/= over total non-low-qual depth >= threshold
//   2) total non-low-qual depth / total depth >= threshold

// 1-homopolymer: AAA 1*3
// 2-homopolymer: CGCGCG 2*3
// N-homopolymer: ACGACGACG N*3
int var_is_homopolymer(char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    // fprintf(stderr, "%d-%d: %s\n", ref_beg, ref_end, ref_seq);
    hts_pos_t start_pos, end_pos; // = var->pos, end_pos = var->pos;
    if (var->var_type == BAM_CDIFF) {
        start_pos = var->pos-1; end_pos = var->pos+1;
    } else if (var->var_type == BAM_CINS) {
        start_pos = var->pos-1; end_pos = var->pos;
    } else { // DEL
        start_pos = var->pos+var->ref_len-1; end_pos = var->pos;
    }
    int var_type = var->var_type;
    int is_homopolymer = 1, max_unit_len=6, n_check_copy_num = 3; // Short Tandem Repeat 1-6 bp
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i)
        ref_bases[i] = nst_nt4_table[(int) ref_seq[end_pos+i-ref_beg]];

    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (nst_nt4_table[(int) ref_seq[end_pos-ref_beg+i*rep_unit_len+j]] != ref_bases[j]) { is_homopolymer = 0; break; }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Homopolymer: %" PRId64 ", %c, unit: %d\n", var->pos, BAM_CIGAR_STR[var_type], rep_unit_len);
            break;
        }
    }
    if (is_homopolymer) return is_homopolymer;
    // reverse
    for (int i = 0; i < 6; ++i) {
        // fprintf(stderr, "%ld, %ld, %ld, %ld\n", start_pos, end_pos, ref_beg, ref_end);
        ref_bases[i] = nst_nt4_table[(int) ref_seq[start_pos-ref_beg-i]];
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (nst_nt4_table[(int) ref_seq[start_pos-ref_beg-i*rep_unit_len-j]] != ref_bases[j]) { is_homopolymer = 0; break; }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Homopolymer: %" PRId64 ", %c, unit: %d (R)\n", var->pos, BAM_CIGAR_STR[var_type], rep_unit_len);
            break;
        }
    }
    return is_homopolymer;
}

// XXX skip detection for long ins/del
int var_is_repeat_region(char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    hts_pos_t pos = var->pos; int ref_len = var->ref_len; int var_type = var->var_type;
    uint8_t *ref_bseq, *alt_bseq; int len = 0;
    int is_repeat = 1;
    if (var_type == BAM_CDEL) { // [pos, pos+del_len*N] vs [pos+del_len, pos+del_len*(N+1)]
        // nst_nt4_table['A'] -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4
        int del_len = ref_len;
        len = del_len * 3; // see if del seq is 3-fold repeat
        if (pos < ref_beg || pos+del_len+len >= ref_end) {
            fprintf(stderr, "DelLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", del_len, ref_len, pos, ref_beg, ref_end);
            return 0;
        }
        ref_bseq = get_bseq(ref_seq+pos-ref_beg, len);
        alt_bseq = get_bseq(ref_seq+pos-ref_beg+del_len, len);
        for (int i = 0; i < len; ++i) {
            if (ref_bseq[i] != alt_bseq[i]) { is_repeat = 0; break; }
        }
        free(ref_bseq); free(alt_bseq);
    } else { // if (var_type == BAM_CINS) { // [ins] * N vs [pos, pos+ ins_len * N]
        // for (int i = 0; i < var->n_uniq_alles-1; ++i) {
        int ins_len = var->alt_len;
        len = ins_len * 3; // see if ins seq is 3-fold repeat
        if (pos < ref_beg || pos+len >= ref_end) {
            fprintf(stderr, "InsLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", ins_len, ref_len, pos, ref_beg, ref_end);
            return 0;
        }
        ref_bseq = get_bseq(ref_seq+pos-ref_beg, len);
        alt_bseq = get_bseq(ref_seq+pos-ref_beg, len);
        for (int j = ins_len; j < len; ++j) alt_bseq[j] = alt_bseq[j-ins_len];
        for (int j = 0; j < ins_len; ++j) alt_bseq[j] = var->alt_seq[j];
        is_repeat = 1;
        for (int k = 0; k < len; ++k) {
            if (ref_bseq[k] != alt_bseq[k]) { is_repeat = 0; break; }
        }
        free(ref_bseq); free(alt_bseq);
        // if (is_repeat == 1) break;
        // }
    }
    if (is_repeat) {
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "RepeatRegion: %" PRId64 ", type: %c, refLen: %d\n", pos, BAM_CIGAR_STR[var_type], ref_len);
    }
    return is_repeat;
}

// 6-tier classification:
// classify all candidate variants into different categories, process each category sequentially
//   0) low total cov, skipped
//   1) clean-region germ+het/hom(SNP or indel): high cov, min_af < AF < max_af - used for phasing
//   2) noisy-region germ+het/hom: high cov, min_af < AF < max_af - used for phasing
//   4) somatic: high cov, AF < min_af - require full phasing info (alt_af < min_af)
int classify_var_cate(const call_var_opt_t *opt, char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end,
                      cand_var_t *var, int min_dp_thres, int min_alt_dp_thres, double min_somatic_af_thres,
                      double min_af_thres, double max_af_thres, double min_noisy_reg_ratio) {
    // total depth < min_dp or total alt depth < min_alt_dp: (skip)
    if (var->total_cov + var->low_qual_cov < min_dp_thres) return LONGCALLD_LOW_COV_VAR;
    // low-qual depth / total depth > min_noisy_reg_ratio: too many low-qual bases, eg., dense X/gap region
    double alt_af=0.0; int alt_dp=0;
    alt_dp = var->alle_covs[1]; alt_af = (double) alt_dp / var->total_cov;
    if (alt_dp < min_alt_dp_thres) return LONGCALLD_LOW_COV_VAR; // skip var with low alt depth 
    if (opt->is_ont && var_is_strand_bias(var)) return LONGCALLD_STRAND_BIAS_VAR; // skip var with strand bias
    // remaining categories: CLEAN_HET(SNP/INDEL), CLEAN_HOM(SNP/INDEL), NOISY_REG_VAR (others)
    if (alt_af < min_af_thres) return LONGCALLD_LOW_AF_VAR; // needs to further check
    if ((double) var->low_qual_cov / (var->low_qual_cov + var->total_cov) >= min_noisy_reg_ratio) return LONGCALLD_NOISY_REG_VAR; // too many low-qual bases, needs to further check
    if (alt_af > max_af_thres) return LONGCALLD_CAND_HOM_VAR; // unlikely germline het., likely hom or somatic, require full phasing info
    // snps & indels in homo/repeat regions
    if ((var->var_type == BAM_CINS || var->var_type == BAM_CDEL) && (var_is_homopolymer(ref_seq, ref_beg, ref_end, var) || var_is_repeat_region(ref_seq, ref_beg, ref_end, var))) return LONGCALLD_REP_HET_VAR; // require basic phasing info, MSA, provide additional phasing info
    // if (var_is_homopolymer(ref_seq, ref_beg, ref_end, var) || var_is_repeat_region(ref_seq, ref_beg, ref_end, var)) return LONGCALLD_REP_HET_VAR; // require basic phasing info, MSA, provide additional phasing info
    // not call somatic variant around homopolymer/repeat region
    // if (alt_af1 < min_af_thres) return LONGCALLD_CAND_SOMA_VAR; // XXX could be het in homopolymer region
    if (var->var_type == BAM_CDIFF) return LONGCALLD_CLEAN_HET_SNP;
    else return LONGCALLD_CLEAN_HET_INDEL;
    // return LONGCALLD_CLEAN_HET_VAR;
}

void copy_var(cand_var_t *to_var, cand_var_t *from_var) {
    to_var->pos = from_var->pos; to_var->var_type = from_var->var_type;
    to_var->total_cov = from_var->total_cov; to_var->low_qual_cov = from_var->low_qual_cov;
    to_var->ref_len = from_var->ref_len;

    // allele_covs, alt_seq, alt_len
    if (to_var->alle_covs) free(to_var->alle_covs);
    to_var->alle_covs = (int*)malloc(from_var->n_uniq_alles * sizeof(int));
    to_var->n_uniq_alles = from_var->n_uniq_alles;

    for (int j = 0; j < to_var->n_uniq_alles; ++j) {
        to_var->alle_covs[j] = from_var->alle_covs[j];
    }
    if (to_var->alt_len != from_var->alt_len) {
        to_var->alt_len = from_var->alt_len;
        free(to_var->alt_seq); to_var->alt_seq = (uint8_t*)malloc(from_var->alt_len * sizeof(uint8_t));
    }
    memcpy(to_var->alt_seq, from_var->alt_seq, from_var->alt_len);
}


// extend noisy regions by flanking length without any vars, up to noisy_reg_flank_len
void post_process_noisy_regs(bam_chunk_t *chunk, call_var_opt_t *opt) {
    cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    int n_noisy_regs = chunk_noisy_regs->n_r;
    // int max_noisy_reg_reads = opt->max_noisy_reg_reads;
    // int *noisy_reg_to_total_n_reads = (int*)calloc(chunk_noisy_regs->n_r, sizeof(int));
    
    // int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    // for (int i = 0; i < chunk->n_reads; ++i) {
    //     if (chunk->is_skipped[i]) continue;
    //     bam1_t *read = chunk->reads[i];
    //     hts_pos_t beg = chunk->digars[i].beg, end = chunk->digars[i].end;
    //     ovlp_n = cr_overlap(chunk_noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
    //     for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
    //         int r_i = ovlp_b[ovlp_i];
    //         noisy_reg_to_total_n_reads[r_i]++;
    //     }
    // }

    // int *first_left_vars = (int*)malloc(n_noisy_regs * sizeof(hts_pos_t)), *first_right_vars = (int*)malloc(n_noisy_regs * sizeof(hts_pos_t));
    // for (int i = 0; i < n_noisy_regs; ++i) {
    //     first_left_vars[i] = first_right_vars[i] = -1;
    // }
    // int var_i = 0, reg_i = 0;
    // for (; var_i < chunk->n_cand_vars && reg_i < n_noisy_regs; ) {
    //     cand_var_t *cand_var = chunk->cand_vars + var_i;
    //     hts_pos_t var_start = cand_var->pos;
    //     // if (cand_var->var_type != BAM_CDIFF) var_start--;
    //     hts_pos_t var_end = var_start + chunk->cand_vars[var_i].ref_len-1;
    //     hts_pos_t reg_start = cr_start(chunk->chunk_noisy_regs, reg_i), reg_end = cr_end(chunk->chunk_noisy_regs, reg_i);
    //     // reg_start
    //     if (var_start > reg_end) {
    //         if (first_right_vars[reg_i] == -1) {
    //             first_right_vars[reg_i] = var_i;
    //         }
    //         reg_i++;
    //     } else if (var_end < reg_start) {
    //         first_left_vars[reg_i] = var_i;
    //         var_i++;
    //     } else {
    //         _err_error_exit("Error: unexpected var/reg overlap: %s:%ld-%ld %s:%ld-%ld\n", chunk->tname, var_start, var_end, chunk->tname, reg_start, reg_end);
    //     }
    // }
    int noisy_reg_flank_len = opt->noisy_reg_flank_len, max_noisy_reg_len = opt->max_noisy_reg_len;
    cgranges_t *noisy_regs = cr_init();
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        // if (noisy_reg_to_total_n_reads[reg_i] > max_noisy_reg_reads) continue; // XXX will be skipped during collect_noisy_vars
        // int first_left_var_i = first_left_vars[reg_i], first_right_var_i = first_right_vars[reg_i];
        // int left_pos = first_left_var_i == -1 ? -1 : chunk->cand_vars[first_left_var_i].pos + chunk->cand_vars[first_left_var_i].ref_len;
        // int right_pos = first_right_var_i == -1 ? INT32_MAX : chunk->cand_vars[first_right_var_i].pos - 1;
        // int noisy_reg_start = MAX_OF_TWO(cr_start(chunk->chunk_noisy_regs, reg_i) - noisy_reg_flank_len, left_pos);
        // int noisy_reg_end = MIN_OF_TWO(cr_end(chunk->chunk_noisy_regs, reg_i) + noisy_reg_flank_len, right_pos);
        int noisy_reg_start = cr_start(chunk_noisy_regs, reg_i) - noisy_reg_flank_len;
        int noisy_reg_end = cr_end(chunk_noisy_regs, reg_i) + noisy_reg_flank_len;
        // if (noisy_reg_end - noisy_reg_start + 1 > max_noisy_reg_len) { // XXX will be skipped during collect_noisy_vars
            // if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Skipped region: %s:%d-%d %d\n", chunk->tname, noisy_reg_start, noisy_reg_end, noisy_reg_end - noisy_reg_start + 1);
        // }
        cr_add(noisy_regs, "cr", noisy_reg_start, noisy_reg_end, cr_label(chunk->chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs); cr_destroy(chunk->chunk_noisy_regs);
    chunk->chunk_noisy_regs = cr_merge(noisy_regs, 0);
    // free(first_left_vars); free(first_right_vars);
    // free(noisy_reg_to_total_n_reads); free(ovlp_b);
}

// XXX excluding all vars in the noisy regions
// update chunk_noisy_regs if any variant is overlapping with it
int classify_cand_vars(bam_chunk_t *chunk, int n_var_sites, call_var_opt_t *opt) {
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end;
    cand_var_t *cand_vars = chunk->cand_vars;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    int *var_i_to_cate = (int*)malloc(n_var_sites * sizeof(int));
    cgranges_t *var_pos_cr = cr_init(); cgranges_t *noisy_var_cr = cr_init(); // overlapping vars: DP needs to be >= min_alt_dp
    cgranges_t *low_comp_cr = chunk->low_comp_cr;
    chunk->var_i_to_cate = (int*)malloc(n_var_sites * sizeof(int));
    int min_dp = opt->min_dp, min_alt_dp = opt->min_alt_dp;
    double min_somatic_af = opt->min_somatic_af, min_af = opt->min_af, max_af = opt->max_af, min_noisy_reg_ratio = opt->min_noisy_reg_ratio;
    int var_cate = -1;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = classify_var_cate(opt, ref_seq, ref_beg, ref_end, var, min_dp, min_alt_dp, min_somatic_af, min_af, max_af, min_noisy_reg_ratio);
        var_i_to_cate[i] = var_cate;
        if (var_cate == LONGCALLD_LOW_COV_VAR) continue; // skipped
        if (opt->is_ont && var_cate == LONGCALLD_STRAND_BIAS_VAR) continue;
        if (var->var_type == BAM_CINS) cr_add(var_pos_cr, "cr", var->pos-1, var->pos, 1);
        else cr_add(var_pos_cr, "cr", var->pos-1, var->pos+var->ref_len-1, 1);
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "CandVarCate-%c: %s:%" PRId64 " %d-%c-%d %d\t", LONGCALLD_VAR_CATE_TYPE(var_cate), chunk->tname, cand_vars[i].pos, cand_vars[i].ref_len, BAM_CIGAR_STR[cand_vars[i].var_type], cand_vars[i].alt_len, cand_vars[i].total_cov);
            fprintf(stderr, "Low-Depth: %d\t", cand_vars[i].low_qual_cov);
            for (int k = 0; k < cand_vars[i].alt_len; ++k) fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seq[k]]);
            fprintf(stderr, ": %d\n", cand_vars[i].alle_covs[1]);
        }
    }
    cr_index(var_pos_cr);
    int64_t var_pos_ovlp_n, noisy_ovlp_n, *ovlp_b = 0, max_b = 0;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = var_i_to_cate[i];
        if (var_cate == LONGCALLD_LOW_COV_VAR) continue;
        // 1. var is in noisy regions
        if (chunk->chunk_noisy_regs != NULL && chunk->chunk_noisy_regs->n_r > 0) {
            int noisy_ovlp_n;
            if (var->var_type == BAM_CINS) noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos, &ovlp_b, &max_b);
            else noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos+var->ref_len-1, &ovlp_b, &max_b);
            if (noisy_ovlp_n > 0) continue; // skip all vars in noisy regions
        }
        // 2. var is overlapping with other vars in ref, i.e., two DELs with overlapping bases
        if (var->var_type == BAM_CINS) var_pos_ovlp_n = cr_overlap(var_pos_cr, "cr", var->pos-1, var->pos, &ovlp_b, &max_b);
        else var_pos_ovlp_n = cr_overlap(var_pos_cr, "cr", var->pos-1, var->pos+var->ref_len-1, &ovlp_b, &max_b);
        if (var_cate == LONGCALLD_NOISY_REG_VAR || var_cate == LONGCALLD_REP_HET_VAR || var_pos_ovlp_n > 1) {
            hts_pos_t var_start, var_end;
            if (var->pos >= reg_beg && var->pos <= reg_end) {
                if (var->var_type == BAM_CINS) {
                    var_start = var->pos; var_end = var->pos;
                } else if (var->var_type == BAM_CDEL) {
                    var_start = var->pos; var_end = var->pos+var->ref_len-1;
                } else {
                    var_start = var->pos; var_end = var->pos+var->ref_len-1;
                }
                if (low_comp_cr != NULL) {
                    int64_t *low_comp_b = 0; int64_t low_comp_n = 0, max_low_comp_n = 0;
                    low_comp_n = cr_overlap(low_comp_cr, "cr", var_start-1, var_end, &low_comp_b, &max_low_comp_n);
                    for (int64_t j = 0; j < low_comp_n; ++j) {
                        int32_t start = cr_start(low_comp_cr, low_comp_b[j]);
                        int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
                        if (start < var_start) var_start = start;
                        if (end > var_end) var_end = end;
                    }
                    free(low_comp_b);
                }
                cr_add(noisy_var_cr, "cr", var_start-1, var_end, 1); // XXX set merge_win as 1, previously: niosy_reg_flank_len
            }
            continue;
        }
        if (var_cate == LONGCALLD_LOW_AF_VAR) var_i_to_cate[i] = LONGCALLD_LOW_COV_VAR;
    }
    if (noisy_var_cr->n_r > 0) {
        cr_index(noisy_var_cr);
        cgranges_t *tmp_cr = cr_merge2(chunk->chunk_noisy_regs, noisy_var_cr, -1);
        cr_destroy(chunk->chunk_noisy_regs);
        chunk->chunk_noisy_regs = tmp_cr;
    }
    post_process_noisy_regs(chunk, opt);
    // skip all vars that overlap with noisy regions
    int cand_var_i = 0;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = var_i_to_cate[i];
        if (var_cate == LONGCALLD_LOW_COV_VAR) continue;
        if (chunk->chunk_noisy_regs != NULL && chunk->chunk_noisy_regs->n_r > 0) {
            int noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos+var->ref_len, &ovlp_b, &max_b);
            if (noisy_ovlp_n > 0) continue; // skip all vars in noisy regions
        }
        // copy i'th to cand_var_i'th
        if (i != cand_var_i) copy_var(cand_vars+cand_var_i, cand_vars+i);
        chunk->var_i_to_cate[cand_var_i++] = var_cate;
    }

    for (int i = cand_var_i; i < n_var_sites; ++i) {
        free(cand_vars[i].alle_covs);
        for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
            free(cand_vars[i].strand_to_alle_covs[j]);
        } free(cand_vars[i].strand_to_alle_covs);
        if (cand_vars[i].alt_seq != NULL) free(cand_vars[i].alt_seq);
    }
    free(var_i_to_cate); cr_destroy(var_pos_cr); free(ovlp_b); cr_destroy(noisy_var_cr);
    return(chunk->n_cand_vars = cand_var_i);
}

int collect_noisy_reg_reads1(bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t niosy_reg_end, int noisy_reg_i, int **noisy_reads) {
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    int n_noisy_reads = 0;
    (*noisy_reads) = (int*)malloc(chunk->n_reads * sizeof(int));

    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        hts_pos_t beg = chunk->digars[i].beg, end = chunk->digars[i].end;
        if (end <= noisy_reg_beg) continue;
        if (beg > niosy_reg_end) break;
        (*noisy_reads)[n_noisy_reads++] = i;
    }
    return n_noisy_reads;
}

void collect_noisy_reg_reads(bam_chunk_t *chunk) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    chunk->noisy_reg_to_n_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    chunk->noisy_reg_to_reads = (int**)malloc(noisy_regs->n_r * sizeof(int*));
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    hts_pos_t beg, end;
    
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            chunk->noisy_reg_to_n_reads[ovlp_b[ovlp_i]]++;
        }
    }
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        chunk->noisy_reg_to_reads[i] = (int*)malloc(chunk->noisy_reg_to_n_reads[i] * sizeof(int));
        chunk->noisy_reg_to_n_reads[i] = 0;
    }
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        // bam1_t *read = chunk->reads[i];
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            chunk->noisy_reg_to_reads[ovlp_b[ovlp_i]][chunk->noisy_reg_to_n_reads[ovlp_b[ovlp_i]]++] = i;
        }
    }
    // print read names for each noisy region
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (cr_label(noisy_regs, i) < 10) continue;
            fprintf(stderr, "NoisyRegReads: %s:%d-%d %d %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i) - cr_start(noisy_regs, i), cr_label(noisy_regs, i));
            for (int j = 0; j < chunk->noisy_reg_to_n_reads[i]; ++j) {
                fprintf(stderr, "Read: %s\n", bam_get_qname(chunk->reads[chunk->noisy_reg_to_reads[i][j]]));
            }
        }
    }
    free(ovlp_b);
}

void collect_digars_from_bam(bam_chunk_t *chunk, const struct call_var_pl_t *pl) {
    chunk->chunk_noisy_regs = cr_init();
    call_var_opt_t *opt = pl->opt;
    for (int i = 0; i < chunk->n_reads; ++i) {
        bam1_t *read = chunk->reads[i];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%d: qname: %s, flag: %d, pos: %" PRId64 ", end: %" PRId64 "\n", i, bam_get_qname(read), read->core.flag, read->core.pos+1, bam_endpos(read));
        if (chunk->is_skipped[i]) continue;
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        if (read->core.qual < opt->min_mq || read->core.flag & BAM_FSUPPLEMENTARY || read->core.flag & BAM_FSECONDARY) {
            chunk->is_skipped[i] = BAM_RECORD_LOW_QUAL; 
            // fprintf(stderr, "Skip read %s with low mapping quality %d\n", bam_get_qname(read), read->core.qual);
            continue;
        }
        int ret;
        if (chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            ret = collect_digar_from_eqx_cigar(chunk, read, pl, opt, chunk->digars+i);
        } else if (chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
            ret = collect_digar_from_MD_tag(chunk, read, pl, opt, chunk->digars+i);
        // XXX TODO use cs tag
        } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
            ret = collect_digar_from_ref_seq(chunk, read, pl, opt, chunk->digars+i);
        }
        if (ret < 0) chunk->is_skipped[i] = BAM_RECORD_WRONG_MAP;
    }

}

int update_digars_with_clip(bam_chunk_t *chunk, call_var_pl_t *pl) {
    return 0;
}

// XXX should be digar->pos-1 for INS/DEL
var_site_t make_var_site_from_digar(int tid, digar1_t *digar) {
    var_site_t var_site = {tid, digar->pos, digar->type, 1, digar->len, digar->alt_seq};
    if (digar->type == BAM_CINS) {
        var_site.ref_len = 0;
    } else if (digar->type == BAM_CDEL) {
        var_site.ref_len = digar->len; var_site.alt_len = 0;
    }
    return var_site;
}

var_site_t make_var_site_from_cand_var(cand_var_t *var) {
    var_site_t var_site = {var->tid, var->pos, var->var_type, var->ref_len, var->alt_len, var->alt_seq};
    return var_site;
}

// check if two var_sites overlap based on their pos & ref_len
int ovlp_var_site(var_site_t *var1, var_site_t *var2) {
    int beg1 = var1->pos; int end1 = var1->pos + var1->ref_len;
    int beg2 = var2->pos; int end2 = var2->pos + var2->ref_len;
    if (var1->ref_len == 0 && var2->ref_len == 0) { // both are INS
        if (beg1 == beg2) return 1;
        else return 0;
    } else if (var1->ref_len == 0) { // var1 is INS, var2 is SNP/DEL
        if (beg1 > beg2 && end1 < end2) return 1;
        else return 0;
    } else if (var2->ref_len == 0) { // var2 is INS, var1 is SNP/DEL
        if (beg2 > beg1 && end2 < end1) return 1;
        else return 0;
    } else { // both are SNP/DEL
        if (beg1 >= end2 || beg2 >= end1) return 0; // no overlap
        return 1; // overlap
    }
}

int ovlp_cand_vars(cand_var_t *var1, cand_var_t *var2) {
    int beg1 = var1->pos; int end1 = var1->pos + var1->ref_len;
    int beg2 = var2->pos; int end2 = var2->pos + var2->ref_len;
    if (var1->ref_len == 0 && var2->ref_len == 0) { // both are INS
        if (beg1 == beg2) return 1;
        else return 0;
    } else if (var1->ref_len == 0) { // var1 is INS, var2 is SNP/DEL
        if (beg1 > beg2 && end1 < end2) return 1;
        else return 0;
    } else if (var2->ref_len == 0) { // var2 is INS, var1 is SNP/DEL
        if (beg2 > beg1 && end2 < end1) return 1;
        else return 0;
    } else { // both are SNP/DEL
        if (beg1 >= end2 || beg2 >= end1) return 0; // no overlap
        return 1; // overlap
    }
}

int comp_ovlp_var_site(var_site_t *var1, var_site_t *var2, int *is_ovlp) {
    *is_ovlp = ovlp_var_site(var1, var2);
    return comp_var_site(var1, var2);
}

// order of variants with same pos: order by type and ref_len
int merge_var_sites(int n_total_var_sites, var_site_t **var_sites, int tid, hts_pos_t reg_beg, hts_pos_t reg_end,
                    int n_digar, digar1_t *digars) {
    int new_total_var_sites = 0;
    var_site_t *new_var_sites = (var_site_t*)malloc((n_total_var_sites + n_digar) * sizeof(var_site_t));
    int i, j;
    for (i = j = 0; i < n_total_var_sites && j < n_digar; ) {
        hts_pos_t digar_pos = digars[j].pos;
        if (reg_beg != -1 && digar_pos < reg_beg) { j++; continue; }
        if (reg_end != -1 && digar_pos > reg_end) { break; }
        // if (!is_in_reg(reg_cr, tname, digar_pos-1, digar_pos)) { j++; continue; }
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
        // if (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL) { // keep noisy-regions variants
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            int ret = comp_var_site((*var_sites)+i, &digar_var_site);
            if (ret < 0) {
                new_var_sites[new_total_var_sites++] = (*var_sites)[i];
                i++;
            } else if (ret > 0) {
                new_var_sites[new_total_var_sites++] = digar_var_site;
                j++;
            } else { // merge
                new_var_sites[new_total_var_sites++] = (*var_sites)[i];
                i++; j++;
            }
        } else j++;
    }
    for (; i < n_total_var_sites; ++i) new_var_sites[new_total_var_sites++] = (*var_sites)[i];
    for (; j < n_digar; ++j) {
        hts_pos_t digar_pos = digars[j].pos;
        if (reg_beg != -1 && digar_pos < reg_beg) { j++; continue; }
        if (reg_end != -1 && digar_pos > reg_end) { break; }
        // if (!is_in_reg(reg_cr, tname, digar_pos-1, digar_pos)) continue;
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            new_var_sites[new_total_var_sites++] = digar_var_site;
        }
    }
    free(*var_sites); *var_sites = new_var_sites;
    return new_total_var_sites;
}

// collect all candidate variant sites from digars, excluding low-quality/noisy-region ones
int collect_all_cand_var_sites(bam_chunk_t *chunk, var_site_t **var_sites) {
    int n_total_var_sites = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
       n_total_var_sites = merge_var_sites(n_total_var_sites, var_sites, 
                                           chunk->tid, chunk->reg_beg, chunk->reg_end, //->beg, chunk->end, 
                                           chunk->digars[i].n_digar, chunk->digars[i].digars);
    }
    // fprintf(stderr, "total_cand_var: %d\n", n_total_var_sites);
    // print cand_vars
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Total candidate variant sites: %d\n", n_total_var_sites);
        for (int i = 0; i < n_total_var_sites; ++i) {
            fprintf(stderr, "CandVarSite: %s:%" PRId64 "-%c\n", chunk->tname, (*var_sites)[i].pos, BAM_CIGAR_STR[(*var_sites)[i].var_type]);
        }
    }
    return n_total_var_sites;
}

int comp_cand_var(cand_var_t *var1, cand_var_t *var2) {
    var_site_t var_site1 = make_var_site_from_cand_var(var1);
    var_site_t var_site2 = make_var_site_from_cand_var(var2);
    return comp_var_site(&var_site1, &var_site2);
}

int merge_var_profile(bam_chunk_t *chunk, int n_reads_with_new_var, int n_new_vars, cand_var_t *new_vars, int *new_var_cate, read_var_profile_t *new_p) {
    if (n_new_vars <= 0) return 0;
    cand_var_t *old_vars = chunk->cand_vars;
    read_var_profile_t *old_p = chunk->read_var_profile;
    cand_var_t *merged_vars = (cand_var_t *)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(cand_var_t)); 
    read_var_profile_t *merged_p = init_read_var_profile(chunk->n_reads, chunk->n_cand_vars + n_new_vars);
    
    int *merged_var_i_to_cate = (int*)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(int));
    int old_var_i = 0, new_var_i = 0, merged_var_i = 0;
    for (; old_var_i < chunk->n_cand_vars && new_var_i < n_new_vars; ) {
        int ret = comp_cand_var(old_vars+old_var_i, new_vars+new_var_i);
        if (ret < 0) { // add old_var to merged_vars, updated read_profile
            for (int i = 0; i < chunk->n_reads; ++i) {
                if (chunk->is_skipped[i]) continue;
                read_var_profile_t *old_p1 = old_p + i;
                read_var_profile_t *merged_p1 = merged_p + i;
                if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
        } else if (ret > 0) { // add new_var to merged vars, update read_profile
            for (int i = 0; i < chunk->n_reads; ++i) {
                if (chunk->is_skipped[i]) continue;
                read_var_profile_t *new_p1 = new_p + i;
                read_var_profile_t *merged_p1 = merged_p + i;
                if (new_p1->start_var_idx > new_var_i || new_p1->end_var_idx < new_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = new_var_cate[new_var_i]; //
            merged_vars[merged_var_i++] = new_vars[new_var_i++];
        } else { // ret == 0
            // always use old_var
            for (int i = 0; i < chunk->n_reads; ++i) {
                if (chunk->is_skipped[i]) continue;
                read_var_profile_t *old_p1 = old_p + i;
                read_var_profile_t *merged_p1 = merged_p + i;
                if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
            // free new_vars
            free_noisy_cand_vars1(new_vars+new_var_i); new_var_i++;
        }
    }
    for (; old_var_i < chunk->n_cand_vars; ++old_var_i) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            if (chunk->is_skipped[i]) continue;
            read_var_profile_t *old_p1 = old_p + i;
            read_var_profile_t *merged_p1 = merged_p + i;
            if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
            update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], merged_p1);
        }
        merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
        merged_vars[merged_var_i++] = old_vars[old_var_i];
    }
    for (; new_var_i < n_new_vars; ++new_var_i) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            if (chunk->is_skipped[i]) continue;
            read_var_profile_t *new_p1 = new_p + i;
            read_var_profile_t *merged_p1 = merged_p + i;
            if (new_p1->start_var_idx > new_var_i || new_p1->end_var_idx < new_var_i) continue;
            update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], merged_p1);
        }
        merged_var_i_to_cate[merged_var_i] = new_var_cate[new_var_i];
        merged_vars[merged_var_i++] = new_vars[new_var_i];
    }
    cgranges_t *merged_read_var_cr = cr_init();
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        cr_add(merged_read_var_cr, "cr", merged_p[i].start_var_idx, merged_p[i].end_var_idx+1, i);
    }
    cr_index(merged_read_var_cr);

    if (chunk->n_cand_vars > 0) free_read_var_profile(chunk->read_var_profile, chunk->n_reads); 
    free_read_var_profile(new_p, chunk->n_reads); free(chunk->var_i_to_cate);
    free(chunk->cand_vars); free(new_vars); cr_destroy(chunk->read_var_cr);
    chunk->read_var_profile = merged_p; chunk->cand_vars = merged_vars; chunk->n_cand_vars = merged_var_i; chunk->var_i_to_cate = merged_var_i_to_cate; chunk->read_var_cr = merged_read_var_cr;

    if (LONGCALLD_VERBOSE >= 2) {
        for (int read_id = 0; read_id < chunk->n_reads; ++read_id) {
            read_var_profile_t *p1 = merged_p + read_id;
            fprintf(stderr, "MergedProfile: %s start_var_i: %d, end_var_i: %d\n", bam_get_qname(chunk->reads[read_id]), p1->start_var_idx, p1->end_var_idx);
            for (int k = 0; k <= p1->end_var_idx-p1->start_var_idx; ++k) {
                fprintf(stderr, "P\tVar: (%d) %" PRId64 "", k, merged_vars[k+p1->start_var_idx].pos);
                fprintf(stderr, " %d-%c-%d, allele: %d\n", merged_vars[k+p1->start_var_idx].ref_len, BAM_CIGAR_STR[merged_vars[k+p1->start_var_idx].var_type], merged_vars[k+p1->start_var_idx].alt_len, p1->alleles[k]);
            }
        }
    }
    return new_var_i; // n_vars
}

read_var_profile_t *collect_read_var_profile(bam_chunk_t *chunk) {
    int n_cand_vars = chunk->n_cand_vars;
    cand_var_t *cand_vars = chunk->cand_vars;
    read_var_profile_t *p = init_read_var_profile(chunk->n_reads, n_cand_vars);
    cgranges_t *read_var_cr = cr_init();
    // 3rd pass: collect read-wise SNP profiles
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        start_var_i = update_read_var_profile_from_digar(chunk->digars+i, read, n_cand_vars, cand_vars, start_var_i, p+i);
        if (p[i].start_var_idx < 0 || p[i].end_var_idx < 0) continue;
        cr_add(read_var_cr, "cr", p[i].start_var_idx, p[i].end_var_idx+1, i);
        if (LONGCALLD_VERBOSE >= 2) {
            if (p[i].start_var_idx >= 0) {
                fprintf(stderr, "Read: %s, start_var_i: %d, end_var_i: %d\n", bam_get_qname(read), p[i].start_var_idx, p[i].end_var_idx);
                for (int j = 0; j <= p[i].end_var_idx-p[i].start_var_idx; ++j) {
                    fprintf(stderr, "P\tVar: (%d) %" PRId64 "", j, cand_vars[j+p[i].start_var_idx].pos);
                    fprintf(stderr, " %d-%c-%d, allele: %d\n", cand_vars[j+p[i].start_var_idx].ref_len, BAM_CIGAR_STR[cand_vars[j+p[i].start_var_idx].var_type], cand_vars[j+p[i].start_var_idx].alt_len, p[i].alleles[j]);
                }
            }
        }
    }
    cr_index(read_var_cr); chunk->read_var_cr = read_var_cr;
    return p;
}

// sample-wise genotype quality
int cal_sample_GQ(int ref_depth, int alt_depth, double logp, double log1p, double log2, int max_gq) {
    int PL[3];
    PL[0] = (int)(-10 * (ref_depth * log1p + alt_depth * logp));
    PL[1] = (int)(10 * (ref_depth + alt_depth) * log2);
    PL[2] = (int)(-10 * (ref_depth * logp + alt_depth * log1p));
    int min_pl = INT_MAX, sec_min_pl = INT_MAX;
    for (int i = 0; i < 3; ++i) {
        if (PL[i] < min_pl) {
            sec_min_pl = min_pl;
            min_pl = PL[i];
        } else if (PL[i] < sec_min_pl) {
            sec_min_pl = PL[i];
        }
    }
    // fprintf(stderr, "%d %d %d %d %d %d %d %d\n", ref_depth, alt_depth, PL[0], PL[1], PL[2], min_pl, sec_min_pl, max_gq);
    int GQ = sec_min_pl - min_pl;
    return MIN_OF_TWO(max_gq, GQ);
}
// var-site-wise QUAL
int cal_var_QUAL1(int ref_depth, int alt_depth, double logp, double log1p, int max_qual) {
    // fprintf(stderr, "QUAL: %d %d %d\n", ref_depth, alt_depth, (int)(-10 * (ref_depth * log1p + alt_depth * logp)));
    return MIN_OF_TWO(max_qual, (int)(-10 * (ref_depth * log1p + alt_depth * logp)));
}

// XXX cand_vars: candidate variants, up to 1 alt_allele
// collect variants based on hap_to_cons_alle
// 1) filter out low-quality variants, e.g., P(var|hap,phasing)
// 2) merge variants with the overlapping pos & ref_len (e.g., ACGT -> A, ACGTCGT) XXX
// 3) extend phase blocks if possible, e.g., if ≥ 1 read supports the longer phase block
int make_variants(const call_var_opt_t *opt, bam_chunk_t *chunk, var_t **_var) {
    int n_cand_vars = chunk->n_cand_vars;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    hts_pos_t active_reg_beg = chunk->reg_beg, active_reg_end = chunk->reg_end;
    if (n_cand_vars <= 0) return 0;
    cand_var_t *cand_vars = chunk->cand_vars;
    (*_var) = (var_t*)malloc(sizeof(var_t));
    var_t *var = *_var;
    var->n = 0; var->m = n_cand_vars;
    int flip = chunk->flip_hap; int hap1_idx, hap2_idx, hom_idx=0;
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    if (flip) {
        hap1_idx = 2; hap2_idx = 1;
    } else {
        hap1_idx = 1; hap2_idx = 2;
    }
    var->vars = (var1_t*)malloc(n_cand_vars * sizeof(var1_t));
    int i = 0, is_hom, hom_alt_is_set, hom_alle, hap1_alle, hap2_alle;
    for (int cand_i = 0; cand_i < n_cand_vars; ++cand_i) {
        if (chunk->var_i_to_cate[cand_i] != LONGCALLD_CLEAN_HET_SNP &&
            chunk->var_i_to_cate[cand_i] != LONGCALLD_CLEAN_HET_INDEL &&
            chunk->var_i_to_cate[cand_i] != LONGCALLD_CAND_HOM_VAR &&
            chunk->var_i_to_cate[cand_i] != LONGCALLD_NOISY_CAND_HET_VAR &&
            chunk->var_i_to_cate[cand_i] != LONGCALLD_NOISY_CAND_HOM_VAR) continue;
        // XXX 
        hom_alle = cand_vars[cand_i].hap_to_cons_alle[hom_idx];
        hap1_alle = cand_vars[cand_i].hap_to_cons_alle[hap1_idx];
        hap2_alle = cand_vars[cand_i].hap_to_cons_alle[hap2_idx];
        is_hom = 0; hom_alt_is_set = 0;
        if (hap1_alle == -1 && hap2_alle == -1) { // homozygous
            is_hom = 1;
            hap1_alle = hap2_alle = hom_alle;
        } else if (hap1_alle == hap2_alle) is_hom = 1;
        if (hap1_alle == -1) hap1_alle = LONGCALLD_REF_ALLELE;
        if (hap2_alle == -1) hap2_alle = LONGCALLD_REF_ALLELE;

        var->vars[i].type = cand_vars[cand_i].var_type;
        var->vars[i].PS = cand_vars[cand_i].phase_set;
        var->vars[i].ref_len = cand_vars[cand_i].ref_len;
        if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) {
            var->vars[i].pos = cand_vars[cand_i].pos-1;
            var->vars[i].ref_len += 1;
        } else var->vars[i].pos = cand_vars[cand_i].pos;
        if (var->vars[i].pos < active_reg_beg || var->vars[i].pos > active_reg_end) continue;
        if (cr_overlap(chunk->reg_cr, chunk->tname, var->vars[i].pos-1, var->vars[i].pos+var->vars[i].ref_len-1, &ovlp_b, &max_b) <= 0)
            continue;
        var->vars[i].ref_bases = get_bseq(ref_seq+var->vars[i].pos-ref_beg, var->vars[i].ref_len);
        if (is_hom) {
            var->vars[i].alt_len = (int*)malloc(1 * sizeof(int));
            var->vars[i].alt_bases = (uint8_t**)malloc(1 * sizeof(uint8_t*));
        } else {
            var->vars[i].alt_len = (int*)malloc(2 * sizeof(int));
            var->vars[i].alt_bases = (uint8_t**)malloc(2 * sizeof(uint8_t*));
        }
        var->vars[i].n_alt_allele = 0;
        for (int hap=1; hap <= 2; ++hap) {
            int hap_alle = hap == 1 ? hap1_alle : hap2_alle;
            if (hap_alle != 0) { // alt allele
                if (is_hom && hom_alt_is_set) {
                    var->vars[i].GT[hap-1] = var->vars[i].n_alt_allele;
                    continue;
                }
                int alt_len = cand_vars[cand_i].alt_len; //s[hap_alle-1];
                var->vars[i].alt_bases[var->vars[i].n_alt_allele] = (uint8_t*)malloc((alt_len+1) * sizeof(uint8_t));
                if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) {
                    alt_len += 1;
                    var->vars[i].alt_bases[var->vars[i].n_alt_allele][0] = nst_nt4_table[(int)ref_seq[var->vars[i].pos-ref_beg]];
                    for (int j = 1; j < alt_len; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j] = cand_vars[cand_i].alt_seq[j-1];
                    }
                } else {
                    for (int j = 0; j < alt_len; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j] = cand_vars[cand_i].alt_seq[j];
                    }
                }
                var->vars[i].alt_len[var->vars[i].n_alt_allele] = alt_len;
                var->vars[i].GT[hap-1] = ++var->vars[i].n_alt_allele;
                if (is_hom) hom_alt_is_set = 1;
            } else var->vars[i].GT[hap-1] = 0;
        }
        var->vars[i].DP = cand_vars[cand_i].total_cov;
        for (int j = 0; j < cand_vars[cand_i].n_uniq_alles; ++j) var->vars[i].AD[j] = cand_vars[cand_i].alle_covs[j];
        var->vars[i].QUAL = cal_var_QUAL1(var->vars[i].AD[0], var->vars[i].AD[1], opt->log_p, opt->log_1p, opt->max_qual);
        var->vars[i].GQ = cal_sample_GQ(var->vars[i].AD[0], var->vars[i].AD[1], opt->log_p, opt->log_1p, opt->log_2, opt->max_gq);
        // var->vars[i].total_depth = cand_snps[i].total_cov;
        // var->vars[i].depths[0] = cand_snps[i].base_covs[cand_snps[i].base_to_i[cand_snps[i].ref_base]];
        // var->vars[i].depths[1] = cand_snps[i].total_cov - var->vars[i].depths[0];
        // var->vars[i].genotype[0] = 0; var->vars[i].genotype[1] = 0;
        i++;
    }
    free(ovlp_b);
    return(var->n = i);
}

// XXX use overlapping reads to extend phase blocks
// 1. check if haplotype of variants in bam_chunk is inconsistent with the previous bam_chunk
// 2. extend phase blocks if possible (e.g., if ≥ 1 read supports the longer phase block)
int flip_variant_hap(bam_chunk_t *prev_chunk, bam_chunk_t *cur_chunk) {
    int n_ovlp_reads = cur_chunk->n_up_ovlp_reads; if (n_ovlp_reads <= 0) return 0;
    if (prev_chunk->n_cand_vars <= 0 || cur_chunk->n_cand_vars <= 0) return 0;
    // 1) find overlapping reads that ovlp with both prev and cur variants
    int *ovlp_read_i = cur_chunk->up_ovlp_read_i; // read_i in the previous bam_chunk
    int *used_ovlp_read_i = (int*)calloc(n_ovlp_reads, sizeof(int)); int n_used_ovlp_reads = 0;
    // hts_pos_t ovlp_start_pos = -1;
    hts_pos_t ovlp_end_pos = 0;
    // int pre_var_start_i=0, pre_var_end_i=prev_chunk->n_cand_vars-1;
    int cur_var_start_i=0, cur_var_end_i=cur_chunk->n_cand_vars-1; // variants covered by overlapping reads
    // hts_pos_t pre_var_end_pos = prev_chunk->cand_vars[pre_var_end_i].pos;
    hts_pos_t cur_var_start_pos = cur_chunk->cand_vars[cur_var_start_i].pos;
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        if (prev_chunk->is_skipped[read_i]) continue;
        // hts_pos_t start_pos = prev_chunk->reads[read_i]->core.pos+1;
        hts_pos_t end_pos = bam_endpos(prev_chunk->reads[read_i]);
        // if (start_pos > pre_var_end_pos) continue;
        if (end_pos < cur_var_start_pos) continue;
        // if (ovlp_start_pos == -1) ovlp_start_pos = start_pos;
        if (end_pos > ovlp_end_pos) ovlp_end_pos = end_pos;
        used_ovlp_read_i[i] = 1; n_used_ovlp_reads++;
    }
    // collect ovlp variants from the previous bam_chunk
    // for (int i = pre_var_end_i; i >= 0; --i) {
    //     cand_var_t *var = prev_chunk->cand_vars + i;
    //    if (var->pos >= ovlp_start_pos) pre_var_start_i = i; else break;
    // }
    // collect ovlp variants from the current bam_chunk
    for (int i = cur_var_start_i; i < cur_chunk->n_cand_vars; ++i) {
        cand_var_t *var = cur_chunk->cand_vars + i;
        if (var->pos <= ovlp_end_pos)
            cur_var_end_i = i;
        else
            break;
    }
    // 2) collect read_var_profile for the overlapping reads & overlapping variants
    // 3) check if haplotype needs to be flipped
    //    var_hap == read_hap: +1
    //    var_hap != read_hap: -1
    //    var_hap == 0: 0
    int keep_hap_score = 0;
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        if (prev_chunk->is_skipped[read_i] || used_ovlp_read_i[i] == 0)
            continue;
        int read_hap = prev_chunk->haps[read_i];
        if (read_hap == 0)
            continue;
        // bam1_t *read = prev_chunk->reads[read_i];
        // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        // int pre_var_start_i0 = prev_chunk->read_var_profile[read_i].start_var_idx;
        // int pre_var_end_i0 = prev_chunk->read_var_profile[read_i].end_var_idx;
        // for (int j = pre_var_start_i; j <= pre_var_end_i; ++j) {
        // if (j < pre_var_start_i0 || j > pre_var_end_i0 || prev_chunk->read_var_profile[read_i].var_is_used[j-pre_var_start_i0] == 0) continue;
        // int allele_i = prev_chunk->read_var_profile[read_i].alleles[j-pre_var_start_i0];
        // fprintf(stderr, "\tPre-Var: %lld, %c, allele: %d\t", (long long) prev_chunk->cand_vars[j].pos, BAM_CIGAR_STR[prev_chunk->cand_vars[j].var_type], allele_i);
        // cand_var_t *var = prev_chunk->cand_vars + j;
        // if (var->hap_to_cons_alle[1] == allele_i) fprintf(stderr, "H1\n");
        // else if (var->hap_to_cons_alle[2] == allele_i) fprintf(stderr, "H2\n");
        // else fprintf(stderr, "none\n");
        // }
        int cur_var_start_i0 = cur_chunk->read_var_profile[i].start_var_idx;
        int cur_var_end_i0 = cur_chunk->read_var_profile[i].end_var_idx;
        for (int j = cur_var_start_i; j <= cur_var_end_i; ++j) {
            if (j < cur_var_start_i0 || j > cur_var_end_i0)
                continue;
            int allele_i = cur_chunk->read_var_profile[i].alleles[j - cur_var_start_i0];
            // fprintf(stderr, "\tCur-Var: %lld, %c, allele: %d\n", (long long) cur_chunk->cand_vars[j].pos, BAM_CIGAR_STR[cur_chunk->cand_vars[j].var_type], allele_i);
            cand_var_t *var = cur_chunk->cand_vars + j;
            if (var->hap_to_cons_alle[read_hap] == allele_i)
                keep_hap_score++;
            else if (var->hap_to_cons_alle[3 - read_hap] == allele_i)
                keep_hap_score--;
        }
    }
    free(used_ovlp_read_i);
    // fprintf(stderr, "pos: %" PRId64 ", keep_hap_score: %d\n", prev_chunk->ref_end, keep_hap_score);
    int flip = (keep_hap_score > 0 ? 0 : 1);
    // 4) update HP tag for the overlapping reads (in prev_bam_chunk, not cur_bam_chunk, to maintain the order in output bam file)
    // since only part of the variants (either prev or cur variants) are used to determine the HP in previous step
    // XXX currently only use HP & variant from the cur_bam_chunk
    // char test_read_name[1024] = "m84039_231005_222902_s1/234751166/ccs";
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        // bam1_t *read = prev_chunk->reads[read_i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
        // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        if (prev_chunk->haps[read_i] != 0 || cur_chunk->haps[i] == 0)
            continue;
        if (flip)
            prev_chunk->haps[read_i] = cur_chunk->haps[i] ^ 3;
        else
            prev_chunk->haps[read_i] = cur_chunk->haps[i];
    }
    // 5) update phase block for the current bam_chunk
    return flip;
}

// update phase set within bam_chunk
int extend_phase_set(bam_chunk_t *chunk) {
    if (chunk->n_cand_vars <= 1)
        return 0;
    hts_pos_t phase_set = -1;
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
    }
return 0;
}

void print_var_seqs(digar1_t *digars, uint8_t **reg_var_seqs, int reg_n_digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < reg_n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d", digars[i].pos, BAM_CIGAR_STR[digars[i].type], digars[i].len, digars[i].qi, digars[i].is_low_qual);
        if (reg_var_seqs[i]) {
            int seq_len = digars[i].len;
            if (digars[i].type == BAM_CDIFF) seq_len *= 2;
            fprintf(fp, "\t");
            for (int j = 0; j < seq_len; ++j) {
                fprintf(fp, "%c", "ACGT"[reg_var_seqs[i][j]]);
            }
        }
        fprintf(fp, "\n");
    }
}

int update_read_var_profile_with_phased_info(bam_chunk_t *chunk, int target_var_cate) {
    return 0;
}

read_var_profile_t *collect_noisy_read_var_profile(bam_chunk_t *chunk, int noisy_reg_i, int n_cand_vars, cand_var_t *cand_vars) {
    return NULL;
}

// excluding non-fully covering reads: not long enough, or clipped
// keep up to 2 medoid reads for each haplotype, in case the phasing was wrong
int *collect_medoid_noisy_reads(bam_chunk_t *chunk, int noisy_reg_i, int *n_medoids, int use_phase_info) {
    // check if alignment-based medoid selection method is needed
    // if not, use XID-profile-based medoid selection method
    return NULL;
}

uint8_t collect_non_gap_char(char *ref_seq, int ref_pos) {
    while (ref_seq[ref_pos] == '-') ref_pos--;
    if (ref_pos < 0) return 4;
    return nst_nt4_table[(int)ref_seq[ref_pos]];
}

uint8_t collect_non_gap_base(uint8_t *ref_seq, int ref_pos) {
    while (ref_seq[ref_pos] == 5) ref_pos--;
    if (ref_pos < 0) return 4;
    return ref_seq[ref_pos];
}

int make_vars_from_baln0(hts_pos_t ref_beg, uint8_t *_ref_seq, uint8_t *_query_seq, int new_msa_len, var_t **vars, int is_hom) {
    *vars = (var_t*)malloc(sizeof(var_t));
    (*vars)->vars = (var1_t*)malloc(new_msa_len * sizeof(var1_t));
    (*vars)->m = new_msa_len;
    hts_pos_t ref_pos = ref_beg;
    int n_vars = 0, i = 0;
    while (i < new_msa_len) {
        if (_ref_seq[i] == _query_seq[i]) {
            i++; ref_pos++;
            continue;
        }
        if (_ref_seq[i] != 5 && _query_seq[i] != 5) {
            // DIFF
            (*vars)->vars[n_vars].type = BAM_CDIFF;
            (*vars)->vars[n_vars].pos = ref_pos;
            (*vars)->vars[n_vars].ref_len = 1;
            (*vars)->vars[n_vars].ref_bases = (uint8_t*)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].ref_bases[0] = _ref_seq[i];
            (*vars)->vars[n_vars].n_alt_allele = 1;
            (*vars)->vars[n_vars].alt_len = (int*)malloc(1 * sizeof(int));
            (*vars)->vars[n_vars].alt_len[0] = 1;
            (*vars)->vars[n_vars].alt_bases = (uint8_t**)malloc(1 * sizeof(uint8_t*));
            (*vars)->vars[n_vars].alt_bases[0] = (uint8_t*)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].alt_bases[0][0] = _query_seq[i];
            (*vars)->vars[n_vars].QUAL = 0;
            (*vars)->vars[n_vars].PS = 0;
            i += 1; ref_pos += 1;
        } else if (_ref_seq[i] == 5) { // INS
            int gap_len = 1;
            while (i+gap_len < new_msa_len && _ref_seq[i+gap_len] == 5 && _query_seq[i+gap_len] != 5) gap_len++;
            (*vars)->vars[n_vars].type = BAM_CINS;
            (*vars)->vars[n_vars].pos = ref_pos-1;
            (*vars)->vars[n_vars].ref_len = 1;
            (*vars)->vars[n_vars].ref_bases = (uint8_t*)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].ref_bases[0] = collect_non_gap_base(_ref_seq, i-1); // _ref_seq[i-1];
            (*vars)->vars[n_vars].n_alt_allele = 1;
            (*vars)->vars[n_vars].alt_len = (int*)malloc(1 * sizeof(int));
            (*vars)->vars[n_vars].alt_len[0] = gap_len + 1;
            (*vars)->vars[n_vars].alt_bases = (uint8_t**)malloc(1 * sizeof(uint8_t*));
            (*vars)->vars[n_vars].alt_bases[0] = (uint8_t*)malloc((gap_len+1) * sizeof(uint8_t));
            (*vars)->vars[n_vars].alt_bases[0][0] = collect_non_gap_base(_query_seq, i-1); // _query_seq[i-1];
            for (int j = 1; j < gap_len+1; ++j)
                (*vars)->vars[n_vars].alt_bases[0][j] = _query_seq[i-1+j];
            (*vars)->vars[n_vars].QUAL = 0;
            (*vars)->vars[n_vars].PS = 0;
            i += gap_len;
        } else if (_query_seq[i] == 5) { // DEL
            int gap_len = 1;
            while (i+gap_len < new_msa_len && _ref_seq[i+gap_len] != 5 && _query_seq[i+gap_len] == 5) gap_len++;
            (*vars)->vars[n_vars].type = BAM_CDEL;
            (*vars)->vars[n_vars].pos = ref_pos-1;
            (*vars)->vars[n_vars].ref_len = gap_len+1;
            (*vars)->vars[n_vars].ref_bases = (uint8_t*)malloc((gap_len+1) * sizeof(uint8_t));
            (*vars)->vars[n_vars].ref_bases[0] = collect_non_gap_base(_ref_seq, i-1); // _ref_seq[i-1];
            for (int j = 1; j <= gap_len; ++j) // collect non-'-' bases
                (*vars)->vars[n_vars].ref_bases[j] = _ref_seq[i-1+j];
            (*vars)->vars[n_vars].n_alt_allele = 1;
            (*vars)->vars[n_vars].alt_len = (int*)malloc(1 * sizeof(int));
            (*vars)->vars[n_vars].alt_len[0] = 1;
            (*vars)->vars[n_vars].alt_bases = (uint8_t**)malloc(1 * sizeof(uint8_t*));
            (*vars)->vars[n_vars].alt_bases[0] = (uint8_t*)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].alt_bases[0][0] = collect_non_gap_base(_query_seq, i-1); // _query_seq[i-1];
            (*vars)->vars[n_vars].QUAL = 0;
            (*vars)->vars[n_vars].PS = 0;
            i += gap_len; ref_pos += gap_len;
        } else {
            _err_error_exit("Error: %d, %d\n", _ref_seq[i], _query_seq[i]);
        }
        if (is_hom) {
            (*vars)->vars[n_vars].GT[0] = 1;
            (*vars)->vars[n_vars].GT[1] = 1;
        } else {
            (*vars)->vars[n_vars].GT[0] = 0;
            (*vars)->vars[n_vars].GT[1] = 0;
        }
        n_vars++;
    }
    return n_vars;
}

// position with both gaps are filtered in advance
int make_cand_vars_from_baln0(int tid, hts_pos_t noisy_reg_beg, hts_pos_t active_reg_beg, hts_pos_t active_reg_end,
                             uint8_t *ref_msa_seq, uint8_t *cons_msa_seq, int msa_len, cand_var_t **cand_vars, int **var_cate) {
    *cand_vars = (cand_var_t *)malloc((msa_len + 1) * sizeof(cand_var_t));
    *var_cate = (int *)malloc((msa_len + 1) * sizeof(int));
    int n_vars = 0;
    hts_pos_t ref_pos = noisy_reg_beg; // int query_pos = 0
    int i = 0;
    while (i < msa_len) {
        if (ref_msa_seq[i] == 5 && cons_msa_seq[i] == 5) {
            i++; continue;
        }
        if (ref_msa_seq[i] == cons_msa_seq[i]) { // EQUAL
            i++; ref_pos++;
            continue;
        }
        if (ref_msa_seq[i] != 5 && cons_msa_seq[i] != 5) { // DIFF
            if (ref_pos < active_reg_beg || ref_pos > active_reg_end) {
                i ++; ref_pos++; continue;
            }
            memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
            (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos;
            (*cand_vars)[n_vars].var_type = BAM_CDIFF;
            (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
            (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
            (*cand_vars)[n_vars].ref_len = 1;
            (*cand_vars)[n_vars].ref_base = ref_msa_seq[i];
            (*cand_vars)[n_vars].alt_len = 1;
            (*cand_vars)[n_vars].alt_seq = (uint8_t*)malloc(sizeof(uint8_t)); // XXX needs to be freed
            (*cand_vars)[n_vars].alt_seq[0] = cons_msa_seq[i];
            n_vars++;
            i += 1; ref_pos += 1;
        } else if (ref_msa_seq[i] == 5) { // INS: 0->I
            int gap_len = 1;
            while (i+gap_len < msa_len && ref_msa_seq[i+gap_len] == 5 && cons_msa_seq[i+gap_len] != 5) gap_len++;
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                i += gap_len; continue;
            }
            memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
            (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos;
            (*cand_vars)[n_vars].var_type = BAM_CINS;
            (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
            (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
            (*cand_vars)[n_vars].ref_len = 0;
            (*cand_vars)[n_vars].alt_len = gap_len;
            (*cand_vars)[n_vars].alt_seq = (uint8_t*)malloc(gap_len * sizeof(uint8_t)); // XXX needs to be freed
            for (int j = 0; j < gap_len; ++j) (*cand_vars)[n_vars].alt_seq[j] = cons_msa_seq[i + j];
            n_vars++;
            i += gap_len;
        } else if (cons_msa_seq[i] == 5) { // DEL: D->0
            int gap_len = 1;
            while (i+gap_len < msa_len && ref_msa_seq[i+gap_len] != 5 && cons_msa_seq[i+gap_len] == 5) gap_len++;
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                i += gap_len; ref_pos += gap_len; continue;
            }
            memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
            (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos;
            (*cand_vars)[n_vars].var_type = BAM_CDEL;
            (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
            (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
            (*cand_vars)[n_vars].ref_len = gap_len;
            (*cand_vars)[n_vars].alt_len = 0;
            (*cand_vars)[n_vars].alt_seq = NULL;
            n_vars++;
            i += gap_len; ref_pos += gap_len;
        } else {
            _err_error_exit("Error: %d, %d\n", ref_msa_seq[i], cons_msa_seq[i]);
        }
    }
    for (int i = 0; i < n_vars; ++i) {
        (*cand_vars)[i].hap_to_alle_profile = NULL; (*cand_vars)[i].hap_to_cons_alle = NULL;
        (*var_cate)[i] = LONGCALLD_NOISY_CAND_HOM_VAR;
    }
    return n_vars;
}

int make_cand_vars_from_msa(int tid, hts_pos_t noisy_reg_beg, hts_pos_t active_reg_beg, hts_pos_t active_reg_end,
                            uint8_t *ref_msa_seq, uint8_t *cons_msa_seq, int msa_len, cand_var_t **cand_vars, int **var_cate) {
    uint8_t *_ref_msa_seq = (uint8_t *)malloc(msa_len * sizeof(uint8_t));
    uint8_t *_cons_msa_seq = (uint8_t *)malloc(msa_len * sizeof(uint8_t));
    int _msa_len = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (ref_msa_seq[i] != 5 || cons_msa_seq[i] != 5) {
            _ref_msa_seq[_msa_len] = ref_msa_seq[i];
            _cons_msa_seq[_msa_len++] = cons_msa_seq[i];
        }
    }
    int n_vars = make_cand_vars_from_baln0(tid, noisy_reg_beg, active_reg_beg, active_reg_end, _ref_msa_seq, _cons_msa_seq, _msa_len, cand_vars, var_cate);
    free(_ref_msa_seq); free(_cons_msa_seq);
    return n_vars;
}

// make cand_vars: type, ref_pos, ref_len, alt_len, alt_seq
// leave total_cov & alle_covs empty, update later
int make_cand_vars_from_aln0(int tid, hts_pos_t noisy_reg_beg, hts_pos_t active_reg_beg, hts_pos_t active_reg_end,
                             uint8_t *query_seq, int qlen, uint32_t *cigar_buf, int cigar_len,
                             cand_var_t **cand_vars) {
    *cand_vars = (cand_var_t *)malloc((qlen + 1) * sizeof(cand_var_t));
    int n_vars = 0;
    hts_pos_t ref_pos = noisy_reg_beg;
    int query_pos = 0;
    for (int i = 0; i < cigar_len; ++i) {
        int len = cigar_buf[i] >> 4;
        int op = cigar_buf[i] & 0xf;
        if (op == BAM_CEQUAL) {
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CDIFF) { // X: 1->1
            for (int j = 0; j < len; ++j) {
                if (ref_pos + j < active_reg_beg || ref_pos + j > active_reg_end)
                    continue;
                memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
                (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos + j;
                (*cand_vars)[n_vars].var_type = BAM_CDIFF;
                (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
                (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
                (*cand_vars)[n_vars].ref_len = 1;
                (*cand_vars)[n_vars].alt_len = 1;
                (*cand_vars)[n_vars].alt_seq = (uint8_t*)malloc(sizeof(uint8_t)); // XXX needs to be freed
                (*cand_vars)[n_vars].alt_seq[0] = query_seq[query_pos + j];
                n_vars++;
            }
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CINS) { // I: 0->n
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                query_pos += len; continue;
            }
            memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
            (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos;
            (*cand_vars)[n_vars].var_type = BAM_CINS;
            (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
            (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
            (*cand_vars)[n_vars].ref_len = 0;
            (*cand_vars)[n_vars].alt_len = len;
            (*cand_vars)[n_vars].alt_seq = (uint8_t*)malloc(len * sizeof(uint8_t)); // XXX needs to be freed
            for (int j = 0; j < len; ++j) (*cand_vars)[n_vars].alt_seq[j] = query_seq[query_pos + j];
            n_vars++;
            query_pos += len;
        } else if (op == BAM_CDEL) { // D: n->0
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                ref_pos += len; continue;
            }
            memset((*cand_vars) + n_vars, 0, sizeof(cand_var_t));
            (*cand_vars)[n_vars].tid = tid; (*cand_vars)[n_vars].pos = ref_pos;
            (*cand_vars)[n_vars].var_type = BAM_CDEL;
            (*cand_vars)[n_vars].n_uniq_alles = 2; // ref allele
            (*cand_vars)[n_vars].alle_covs = (int *)calloc(2, sizeof(int));
            (*cand_vars)[n_vars].ref_len = len;
            (*cand_vars)[n_vars].alt_len = 0;
            (*cand_vars)[n_vars].alt_seq = NULL;
            n_vars++;
            ref_pos += len;
        } else {
            _err_error_exit("Error: %d\n", op);
        }
    }
    return n_vars;
}

int make_vars_from_aln0(char *chunk_ref_seq, hts_pos_t chunk_ref_seq_beg, hts_pos_t chunk_ref_seq_end, hts_pos_t noisy_reg_beg,
                        hts_pos_t active_reg_beg, hts_pos_t active_reg_end,
                        uint8_t *query_seq, int qlen, uint32_t *cigar_buf, int cigar_len,
                        var_t **vars, int is_hom) {
    *vars = (var_t *)malloc(sizeof(var_t));
    (*vars)->vars = (var1_t *)malloc((qlen + 1) * sizeof(var1_t));
    (*vars)->m = (qlen + 1);
    hts_pos_t ref_pos = noisy_reg_beg; int query_pos = 0;
    int n_vars = 0, i = 0;
    for (int i = 0; i < cigar_len; ++i) {
        int len = cigar_buf[i] >> 4;
        int op = cigar_buf[i] & 0xf;
        if (op == BAM_CEQUAL) {
            ref_pos += len; query_pos += len;
        } else if (op == BAM_CDIFF) { // DIFF
            for (int j = 0; j < len; ++j) {
                if (ref_pos + j < active_reg_beg || ref_pos + j > active_reg_end) continue;
                (*vars)->vars[n_vars].type = BAM_CDIFF;
                (*vars)->vars[n_vars].pos = ref_pos + j;
                (*vars)->vars[n_vars].ref_len = 1;
                (*vars)->vars[n_vars].ref_bases = (uint8_t *)malloc(1 * sizeof(uint8_t));
                (*vars)->vars[n_vars].ref_bases[0] = nst_nt4_table[(int)chunk_ref_seq[ref_pos - chunk_ref_seq_beg + j]];
                (*vars)->vars[n_vars].n_alt_allele = 1;
                (*vars)->vars[n_vars].alt_len = (int *)malloc(1 * sizeof(int));
                (*vars)->vars[n_vars].alt_len[0] = 1;
                (*vars)->vars[n_vars].alt_bases = (uint8_t **)malloc(1 * sizeof(uint8_t *));
                (*vars)->vars[n_vars].alt_bases[0] = (uint8_t *)malloc(1 * sizeof(uint8_t));
                (*vars)->vars[n_vars].alt_bases[0][0] = query_seq[query_pos + j];
                (*vars)->vars[n_vars].QUAL = 0;
                (*vars)->vars[n_vars].PS = 0;
                if (is_hom) {
                    (*vars)->vars[n_vars].GT[0] = 1;
                    (*vars)->vars[n_vars].GT[1] = 1;
                }
                n_vars++;
            }
            ref_pos += len;
            query_pos += len;
        } else if (op == BAM_CINS) { // INS
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                query_pos += len; continue;
            }
            (*vars)->vars[n_vars].type = BAM_CINS;
            (*vars)->vars[n_vars].pos = ref_pos - 1;
            (*vars)->vars[n_vars].ref_len = 1;
            (*vars)->vars[n_vars].ref_bases = (uint8_t *)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].ref_bases[0] = nst_nt4_table[(int)chunk_ref_seq[ref_pos - 1 - chunk_ref_seq_beg]];
            (*vars)->vars[n_vars].n_alt_allele = 1;
            (*vars)->vars[n_vars].alt_len = (int *)malloc(1 * sizeof(int));
            (*vars)->vars[n_vars].alt_len[0] = len + 1;
            (*vars)->vars[n_vars].alt_bases = (uint8_t **)malloc(1 * sizeof(uint8_t *));
            (*vars)->vars[n_vars].alt_bases[0] = (uint8_t *)malloc((len + 1) * sizeof(uint8_t));
            (*vars)->vars[n_vars].alt_bases[0][0] = (*vars)->vars[n_vars].ref_bases[0];
            for (int j = 1; j <= len; ++j)
                (*vars)->vars[n_vars].alt_bases[0][j] = query_seq[query_pos + j - 1];
            (*vars)->vars[n_vars].QUAL = 0;
            (*vars)->vars[n_vars].PS = 0;
            if (is_hom) {
                (*vars)->vars[n_vars].GT[0] = 1;
                (*vars)->vars[n_vars].GT[1] = 1;
            }
            n_vars++; query_pos += len;
        } else if (op == BAM_CDEL) { // DEL
            if (ref_pos - 1 < active_reg_beg || ref_pos - 1 > active_reg_end) {
                ref_pos += len; continue;
            }
            (*vars)->vars[n_vars].type = BAM_CDEL;
            (*vars)->vars[n_vars].pos = ref_pos - 1;
            (*vars)->vars[n_vars].ref_len = len + 1;
            (*vars)->vars[n_vars].ref_bases = (uint8_t *)malloc((len + 1) * sizeof(uint8_t));
            for (int j = 0; j <= len; ++j) // collect non-'-' bases
                (*vars)->vars[n_vars].ref_bases[j] = nst_nt4_table[(int)chunk_ref_seq[ref_pos - 1 - chunk_ref_seq_beg + j]];
            (*vars)->vars[n_vars].n_alt_allele = 1;
            (*vars)->vars[n_vars].alt_len = (int *)malloc(1 * sizeof(int));
            (*vars)->vars[n_vars].alt_len[0] = 1;
            (*vars)->vars[n_vars].alt_bases = (uint8_t **)malloc(1 * sizeof(uint8_t *));
            (*vars)->vars[n_vars].alt_bases[0] = (uint8_t *)malloc(1 * sizeof(uint8_t));
            (*vars)->vars[n_vars].alt_bases[0][0] = (*vars)->vars[n_vars].ref_bases[0]; // ref base
            (*vars)->vars[n_vars].QUAL = 0;
            (*vars)->vars[n_vars].PS = 0;
            if (is_hom) {
                (*vars)->vars[n_vars].GT[0] = 1;
                (*vars)->vars[n_vars].GT[1] = 1;
            }
            n_vars++;
            ref_pos += len;
        } else {
            _err_error_exit("Error: %d\n", op);
        }
    }
    (*vars)->n = n_vars;
    return n_vars;
}

// order of variants with same pos: order by type, ref_len, alt_len
int comp_var_site(var_site_t *var1, var_site_t *var2) {
    if (var1->pos < var2->pos) return -1;
    if (var1->pos > var2->pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->alt_len < var2->alt_len) return -1;
    if (var1->alt_len > var2->alt_len) return 1;
    // fprintf(stderr, "Error: %ld, %d, %d, %d, %d\n", var1->pos, var1->var_type, var1->ref_len, var1->alt_len, 
                                                        // var2->alt_len);
    // fprintf(stderr, "1alt_seq: %d\n", var1->alt_seq[0]);
    // fprintf(stderr, "2alt_seq: %d\n", var2->alt_seq[0]);

    if (var1->var_type == BAM_CINS || var1->var_type == BAM_CDIFF) 
        return memcmp(var1->alt_seq, var2->alt_seq, var1->alt_len);
    else return 0;
}


// check if read_seq and target_seq are identical at positions [target_pos...target_pos+len-1]
// XXX reads with part of the inserted sequence ???
// for large insertion, check if similairty are >= 0.9, i.e., #(#) / #(total) >= 0.9
// read: |-------|
// cons: |  INS  |
int is_match(uint8_t *read_seq, uint8_t *target_seq, int msa_len, int target_pos, int len, float sim_thres) {
    int cur_pos = -1, n_eq = 0, n_xid = 0;
    int checked_base = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (target_seq[i] != 5) cur_pos++;
        if (cur_pos == target_pos + len) break;
        else if (cur_pos >= target_pos) {
            checked_base++;
            if (read_seq[i] == target_seq[i]) n_eq++;
            else n_xid++;
        }
        if (checked_base == len) break;
    }
    // if (len >= 100) return (n_eq >= (len * sim_thres));
    if (len >= 10) return (n_eq >= (len * sim_thres));
    else return (n_eq == len && n_xid == 0);
}

// 0: ref, 1: alt, -1: unknown
int is_match_aln_str(aln_str_t *aln_str, int target_pos, int len, float cons_sim_thres, int *full_cover) {
    int cur_pos = -1, n_eq = 0, n_xid = 0;
    int checked_base = 0;
    int cover_start = 0, cover_end = 0;
    int start_pos, end_pos;
    if (target_pos < 0) {
        start_pos = 0; end_pos = len - 1;
    } else {
        start_pos = target_pos; end_pos = target_pos + len - 1;
    }
    for (int i = 0; i < aln_str->aln_len; ++i) {
        if (aln_str->target_aln[i] != 5) cur_pos++;
        if (cur_pos == target_pos + len) break;
        if (i < aln_str->query_beg || i < aln_str->target_beg) continue;
        else if (i > aln_str->query_end || i > aln_str->target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;

        if (cur_pos >= target_pos) {
            checked_base++;
            if (aln_str->query_aln[i] == aln_str->target_aln[i]) n_eq++;
            else n_xid++;
        }
        // if (checked_base == len) break;
    }
    if (cover_start && cover_end) *full_cover = 1;
    else *full_cover = 0;
    if (len >= 10) {
        if (n_eq >= (len * cons_sim_thres)) return 1;
        else if (cover_start && cover_end) return 0;
        else return -1;
    } else {
        if (n_eq == len && n_xid == 0) return 1;
        else if (cover_start && cover_end) return 0;
        else return -1;
    }
}

// 0: ref, 1: alt, -1: unknown
int is_match_aln_str_del(aln_str_t *aln_str, int target_del_left, int target_del_right, float sim_thres, int *full_cover) {
    int n_del = 0, cur_pos = -1;
    int started_check_del = 0, n_non_del = 0;
    int cover_start = 0, cover_end = 0;
    int start_pos, end_pos;
    if (target_del_left < 0) {
        start_pos = 0; end_pos = target_del_right;
    } else {
        start_pos = target_del_left; end_pos = target_del_right;
    }
    for (int i = 0; i < aln_str->aln_len; ++i) {
        if (aln_str->target_aln[i] != 5) cur_pos++;
        if (cur_pos > target_del_right) break;

        if (i < aln_str->query_beg || i < aln_str->target_beg) continue;
        else if (i > aln_str->query_end || i > aln_str->target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;

        if (cur_pos >= target_del_left && cur_pos < target_del_right) {
            if (started_check_del == 0) {
                started_check_del = 1;
            } else {
                if (aln_str->query_aln[i] != 5) n_non_del++;
            }
        }
    }
    if (cover_start && cover_end) {
        *full_cover = 1;
        if (n_non_del == 0) return 1;
        else return 0;
    } else {
        *full_cover = 0;
        return -1;
    }
}

// XXX double-check if this is correct
int is_match_del(uint8_t *read_seq, uint8_t *target_seq, int msa_len, int target_del_left, int target_del_right, float sim_thres) {
    int n_del = 0, cur_pos = -1;
    int started_check_del = 0, n_non_del = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (target_seq[i] != 5) cur_pos++;
        if (cur_pos >= target_del_right) break;
        else if (cur_pos >= target_del_left) {
            if (started_check_del == 0) {
                started_check_del = 1;
            } else {
                if (read_seq[i] != 5) n_non_del++;
            }
        }
    }
    if (n_non_del > 0) return 0;
    else return 1;
}

int get_var_allele_i_from_cons_aln_str(aln_str_t *cons_aln_str, int var_type, int alt_pos, int alt_len, float cons_sim_thres, int *full_cover) {
    *full_cover = 0;
    if (var_type == BAM_CDIFF) {
        assert(alt_len == 1);
        return is_match_aln_str(cons_aln_str, alt_pos, 1, cons_sim_thres, full_cover);
    } else if (var_type == BAM_CINS) {
        return is_match_aln_str(cons_aln_str, alt_pos, alt_len, cons_sim_thres, full_cover);
    } else if (var_type == BAM_CDEL) { // XXX large DEL ???
        return is_match_aln_str_del(cons_aln_str, alt_pos-1, alt_pos, cons_sim_thres, full_cover);
    }
    return -1;
}

int is_cover_aln_str(aln_str_t *aln_str, int target_pos, int len) {
    int cur_pos = -1;
    int cover_start = 0, cover_end = 0;
    int start_pos, end_pos;
    if (target_pos < 0) {
        start_pos = 0; 
        end_pos = len - 1;
    } else {
        start_pos = target_pos;
        end_pos = target_pos + len - 1;
    }

    for (int i = 0; i < aln_str->aln_len; ++i) {
        if (aln_str->target_aln[i] != 5) cur_pos++;
        if (i < aln_str->query_beg || i < aln_str->target_beg) continue;
        else if (i > aln_str->query_end || i > aln_str->target_end) break;

        if (cur_pos == start_pos) cover_start = 1;
        if (cur_pos == end_pos) cover_end = 1;
        if (cover_start && cover_end) return 1;
    }
    return 0;
}

int get_full_cover_from_cons_aln_str(aln_str_t *cons_aln_str, int var_type, int alt_pos, int ref_len) {
    if (var_type == BAM_CDIFF) {
        assert(ref_len == 1);
        return is_cover_aln_str(cons_aln_str, alt_pos, 1);
    } else if (var_type == BAM_CINS) {
        return is_cover_aln_str(cons_aln_str, alt_pos, ref_len+1);
    } else if (var_type == BAM_CDEL) { // XXX large DEL ???
        return is_cover_aln_str(cons_aln_str, alt_pos-1, ref_len+1);
    }
    return 0;
}

// cons_aln_str: cons vs read
// cand_var: update total_cov & alle_covs (only count alt_allele, no ref_allele)
// p: update read_id,
void update_cand_var_profile_from_cons_aln_str(aln_str_t *cons_aln_str, hts_pos_t ref_pos_beg, cand_var_t *vars, int n_vars, read_var_profile_t *p) {
    uint8_t *_cons_seq = cons_aln_str->target_aln; uint8_t *_cons_read_seq = cons_aln_str->query_aln; int _cons_aln_len = cons_aln_str->aln_len;
    float cons_sim_thres = 0.9; // XXX
    
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "cons:\n");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_seq[i]]); fprintf(stderr, "\n");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_read_seq[i]]); fprintf(stderr, "\n");
    }
    // X: base is identical
    // I: inserted bases are identical
    // D: two flanking bases are =/X, no I/D
    // if var matches the pattern, set var_profile as 1, otherwise -1 (unused)
    int allele_i = -1;
    int delta_ref_alt = 0;
    for (int i = 0; i < n_vars; ++i) {
        int var_ref_pos = vars[i].pos - ref_pos_beg;
        int var_ref_len = vars[i].ref_len, var_alt_len = vars[i].alt_len;
        uint8_t *alt_seq = vars[i].alt_seq; int full_cover = 0;
        allele_i = get_var_allele_i_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_ref_pos-delta_ref_alt, var_alt_len, cons_sim_thres, &full_cover);
        if (full_cover) {
            vars[i].total_cov++;
            if (allele_i != -1) vars[i].alle_covs[allele_i]++;
            update_read_var_profile_with_allele(i, allele_i, p);
        }
        if (vars[i].var_type == BAM_CINS) delta_ref_alt -= var_alt_len;
        else if (vars[i].var_type == BAM_CDEL) delta_ref_alt += var_ref_len;
    }
}

// aln_strs:
// 0: ref vs cons
// 1..clu_n_seqs: cons vs read
// clu_n_seqs+1..clu_n_seqs*2: ref vs read
void update_cand_var_profile_from_cons_aln_str1(int clu_n_seqs, int *clu_read_ids, aln_str_t *clu_aln_strs, hts_pos_t ref_pos_beg,
                                           cand_var_t *noisy_vars, int n_vars, read_var_profile_t *p) {
    // uint8_t *ref_seq = msa_seqs[clu_n_seqs+1];
    for (int i = 0; i < clu_n_seqs; ++i) {
        aln_str_t *cons_aln_str = LONGCALLD_CONS_READ_ALN_STR(clu_aln_strs, i);
        int read_id = clu_read_ids[i];
        update_cand_var_profile_from_cons_aln_str(cons_aln_str, ref_pos_beg, noisy_vars, n_vars, p+read_id);
    }
}

// use cons and ref to update var profile
void update_cand_var_profile_from_cons_aln_str21(int clu_idx, aln_str_t *cons_aln_str, hts_pos_t ref_pos_beg,
                                                 cand_var_t *vars, int n_vars, int *var_from_cons_idx, read_var_profile_t *p) {
    uint8_t *_cons_seq = cons_aln_str->target_aln; uint8_t *_cons_read_seq = cons_aln_str->query_aln; int _cons_aln_len = cons_aln_str->aln_len;
    float cons_sim_thres = 0.9, ref_sim_thres = 0.9, alt_sim_thres = 1.0;

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "cons:\t");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_seq[i]]); fprintf(stderr, "\nread:\t");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_read_seq[i]]); fprintf(stderr, "\n");
    }
    int allele_i = -1, delta_ref_alt = 0;

    for (int i = 0; i < n_vars; ++i) {
        int var_ref_pos = vars[i].pos - ref_pos_beg, var_ref_len = vars[i].ref_len, var_alt_len = vars[i].alt_len;
        uint8_t *var_alt_seq = vars[i].alt_seq; int full_cover = 0;
        if (var_from_cons_idx[i] & clu_idx) { // check if read contain the var (1)
            allele_i = get_var_allele_i_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_ref_pos-delta_ref_alt, var_alt_len, cons_sim_thres, &full_cover);
        } else {
            // XXX check if read covers the var
            full_cover = get_full_cover_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_ref_pos-delta_ref_alt, var_ref_len);
            allele_i = 0;
        }
        if (full_cover) {
            vars[i].total_cov++;
            if (allele_i != -1) vars[i].alle_covs[allele_i]++;
            update_read_var_profile_with_allele(i, allele_i, p);
        }
        // update delta_ref_alt for both cons/var
        if (vars[i].var_type == BAM_CINS) {
            if ((var_from_cons_idx[i] & clu_idx) > 0) {
                delta_ref_alt -= var_alt_len;
            }
        } else if (vars[i].var_type == BAM_CDEL) {
            if ((var_from_cons_idx[i] & clu_idx) > 0) {
                delta_ref_alt += var_ref_len;
            }
        }
    }
}

// msa_seqs: clu_n_seqs + cons + cons from the other haplotype
int update_cand_var_profile_from_cons_aln_str2(bam_chunk_t *chunk, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, 
                                          hts_pos_t noisy_reg_beg, cand_var_t *hap1_vars, int *hap1_var_cate, int n_hap1_vars, cand_var_t *hap2_vars, int *hap2_var_cate, int n_hap2_vars, 
                                          cand_var_t **noisy_vars, int **noisy_var_cate, read_var_profile_t **p) {
    if (n_hap1_vars + n_hap2_vars == 0) return 0;
    int n_vars = 0;
    (*noisy_vars) = (cand_var_t *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(cand_var_t));
    (*noisy_var_cate) = (int *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(int));
    int i1 = 0, i2 = 0;
    int *var_from_cons_idx = (int *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(int));
    for (; i1 < n_hap1_vars && i2 < n_hap2_vars; ) {
        int ret = comp_cand_var(hap1_vars + i1, hap2_vars + i2);
        if (ret < 0) {
            (*noisy_var_cate)[n_vars] = LONGCALLD_NOISY_CAND_HET_VAR;
            var_from_cons_idx[n_vars] = 1;
            (*noisy_vars)[n_vars++] = hap1_vars[i1++];
        } else if (ret > 0) {
            (*noisy_var_cate)[n_vars] = LONGCALLD_NOISY_CAND_HET_VAR;
            var_from_cons_idx[n_vars] = 2;
            (*noisy_vars)[n_vars++] = hap2_vars[i2++];
        } else {
            (*noisy_var_cate)[n_vars] = LONGCALLD_NOISY_CAND_HOM_VAR;
            var_from_cons_idx[n_vars] = 3; // both
            (*noisy_vars)[n_vars++] = hap1_vars[i1++];
            free_noisy_cand_vars1(hap2_vars + i2); i2++;
        }
    } 
    for (; i1 < n_hap1_vars; ++i1) {
        (*noisy_var_cate)[n_vars] = LONGCALLD_NOISY_CAND_HET_VAR;
        var_from_cons_idx[n_vars] = 1;
        (*noisy_vars)[n_vars++] = hap1_vars[i1];
    }
    for (; i2 < n_hap2_vars; ++i2) {
        (*noisy_var_cate)[n_vars] = LONGCALLD_NOISY_CAND_HET_VAR;
        var_from_cons_idx[n_vars] = 2;
        (*noisy_vars)[n_vars++] = hap2_vars[i2];
    }
    // collect read_var_profile for each read and each var
    *p = init_read_var_profile(chunk->n_reads, n_vars);
    int p_i = 0;
    for (int i = 0; i < 2; ++i) {
        aln_str_t *clu_aln_str = aln_strs[i];
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_id = clu_read_ids[i][j];
            read_var_profile_t *p1 = *p + read_id;
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read: %s\n", bam_get_qname(chunk->reads[read_id]));
            aln_str_t *cons_aln_str = LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, j);
            update_cand_var_profile_from_cons_aln_str21(i+1, cons_aln_str, noisy_reg_beg, *noisy_vars, n_vars, var_from_cons_idx, p1);
            p_i++;
        }
    }
    free(var_from_cons_idx);
    return n_vars;
}
// XXX allow partial-aligned reads
// msa_seqs: ref + cons + clu_n_seqs
int make_vars_from_msa_cons_aln(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_reads, int *read_ids, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, hts_pos_t active_reg_beg, hts_pos_t active_reg_end,
                                int n_cons, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs,
                                cand_var_t **noisy_vars, int **noisy_var_cate, read_var_profile_t **p) {
    if (n_cons == 0) return 0;
    // char *ref_cseq = NULL, *chunk_ref_seq = chunk->ref_seq; 
    // int ref_seq_len = collect_reg_ref_cseq(chunk, noisy_reg_beg, noisy_reg_end, &ref_cseq);
    // hts_pos_t chunk_ref_seq_beg = chunk->ref_beg, chunk_ref_seq_end = chunk->ref_end;
    // uint8_t *ref_seq = NULL; int ref_seq_len = collect_reg_ref_bseq(chunk, noisy_reg_beg, noisy_reg_end, &ref_seq);


    // WFA for ref vs cons1/2
    int n_vars = 0, n_hap1_vars = 0, n_hap2_vars = 0;
    cand_var_t *hap1_vars = NULL, *hap2_vars = NULL; int *hap1_var_cate=NULL, *hap2_var_cate=NULL; read_var_profile_t *p1 = NULL, *p2 = NULL;
    // uint8_t **ref_cons_msa_seq = NULL; // size: (1+n_cons) * ref_cons_msa_len, ref/cons1/cons2
    // int ref_cons_msa_len = collect_ref_cons_msa_seqs(ref_seq, ref_seq_len, cons_seqs, cons_lens, n_cons, &ref_cons_msa_seq);
    for (int i = 0; i < n_cons; ++i) {
        uint32_t *cigar_buf = NULL;
        // int cigar_len = end2end_aln(opt, ref_cseq, ref_seq_len, cons_seqs[i], cons_lens[i], &cigar_buf);
        if (n_cons == 1) {
            aln_str_t *clu_aln_str = aln_strs[0];
            aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
            uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
            n_vars = make_cand_vars_from_msa(chunk->tid, noisy_reg_beg, active_reg_beg, active_reg_end, 
                                             ref_seq_aln, cons_seq_aln, aln_len, noisy_vars, noisy_var_cate);
            if (n_vars > 0) {
                *p = init_read_var_profile(chunk->n_reads, n_vars);
                update_cand_var_profile_from_cons_aln_str1(clu_n_seqs[0], clu_read_ids[0], clu_aln_str, noisy_reg_beg, *noisy_vars, n_vars, *p);
            } else free(*noisy_vars);
        } else {
            if (i == 0) {
                aln_str_t *clu_aln_str = aln_strs[0];
                aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
                uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
                n_hap1_vars = make_cand_vars_from_msa(chunk->tid, noisy_reg_beg, active_reg_beg, active_reg_end,
                                                      ref_seq_aln, cons_seq_aln, aln_len, &hap1_vars, &hap1_var_cate);
            } else if (i == 1) {
                aln_str_t *clu_aln_str = aln_strs[1];
                aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
                uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
                n_hap2_vars = make_cand_vars_from_msa(chunk->tid, noisy_reg_beg, active_reg_beg, active_reg_end,
                                                      ref_seq_aln, cons_seq_aln, aln_len, &hap2_vars, &hap2_var_cate);
            }
        }
        if (cigar_buf != NULL) free(cigar_buf);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_vars; ++i)
            fprintf(stderr, "HOM: %" PRId64 ", %d-%c-%d\n", (*noisy_vars)[i].pos, (*noisy_vars)[i].ref_len, BAM_CIGAR_STR[(*noisy_vars)[i].var_type], (*noisy_vars)[i].alt_len);
    }
    if (n_cons == 2) {
        n_vars = update_cand_var_profile_from_cons_aln_str2(chunk, clu_n_seqs, clu_read_ids, aln_strs, noisy_reg_beg,
                                                            hap1_vars, hap1_var_cate, n_hap1_vars, hap2_vars, hap2_var_cate, n_hap2_vars, noisy_vars, noisy_var_cate, p);
        if (LONGCALLD_VERBOSE >= 2 && n_vars > 0) {
            for (int i = 0; i < n_vars; ++i) {
                fprintf(stderr, "Var: %" PRId64 ", %d-%c-%d %d(%d,%d) %c\t", (*noisy_vars)[i].pos, (*noisy_vars)[i].ref_len, BAM_CIGAR_STR[(*noisy_vars)[i].var_type], (*noisy_vars)[i].alt_len,
                                (*noisy_vars)[i].total_cov, (*noisy_vars)[i].alle_covs[0], (*noisy_vars)[i].alle_covs[1], LONGCALLD_VAR_CATE_TYPE((*noisy_var_cate)[i]));
                for (int l = 0; l < (*noisy_vars)[i].alt_len; ++l) {
                    fprintf(stderr, "%c", "ACGTN-"[(*noisy_vars)[i].alt_seq[l]]);
                } fprintf(stderr, "\n");
            }
            for (int i = 0; i < n_cons; ++i) {
                for (int j = 0; j < clu_n_seqs[i]; ++j) {
                    int read_id = clu_read_ids[i][j];
                    read_var_profile_t *p1 = *p + read_id;
                    fprintf(stderr, "Hap%d-Read: %s start_var_i: %d, end_var_i: %d\n", i+1, bam_get_qname(chunk->reads[read_id]), p1->start_var_idx, p1->end_var_idx);
                    for (int k = 0; k <= p1->end_var_idx-p1->start_var_idx; ++k) {
                        fprintf(stderr, "P\tVar: (%d) %" PRId64 "", k, (*noisy_vars)[k+p1->start_var_idx].pos);
                        fprintf(stderr, " %d-%c-%d, allele: %d\n", (*noisy_vars)[k+p1->start_var_idx].ref_len, BAM_CIGAR_STR[(*noisy_vars)[k+p1->start_var_idx].var_type], (*noisy_vars)[k+p1->start_var_idx].alt_len, p1->alleles[k]);
                    }
                }
            }
        }
        free(hap1_vars); free(hap1_var_cate); free(hap2_vars); free(hap2_var_cate);
    }
    return n_vars;
}

// copy add_vars and old_vars to new_vars, merge them based on pos
// free unused add_vars
int merge_vars(var_t *old_vars, var_t *add_vars) {
    var1_t *new_vars = (var1_t*)malloc((old_vars->n + add_vars->n) * sizeof(var1_t));
    int old_i = 0, add_i = 0, new_i = 0;
    for (; old_i < old_vars->n && add_i < add_vars->n; ) {
        if (old_vars->vars[old_i].pos < add_vars->vars[add_i].pos) {
            new_vars[new_i++] = old_vars->vars[old_i++];
        } else if (old_vars->vars[old_i].pos > add_vars->vars[add_i].pos) {
            new_vars[new_i++] = add_vars->vars[add_i++];
        } else {
            if (old_vars->vars[old_i].ref_len == add_vars->vars[add_i].ref_len &&
                old_vars->vars[old_i].type == add_vars->vars[add_i].type && 
                old_vars->vars[old_i].n_alt_allele == add_vars->vars[add_i].n_alt_allele) {
                int same = 1;
                for (int i = 0; i < old_vars->vars[old_i].n_alt_allele-1; ++i) {
                    if (old_vars->vars[old_i].alt_len[i] != add_vars->vars[add_i].alt_len[i]) {
                        same = 0; break;
                    }
                }
                if (same) {
                    new_vars[new_i++] = old_vars->vars[old_i++];
                    var1_free(add_vars->vars+add_i);
                    add_i++;
                    continue;
                }
            } else {
                if (LONGCALLD_VERBOSE >= 2)
                    fprintf(stderr, "DIFF: old: %ld %d %c %d add: %ld %d %c %d\n", old_vars->vars[old_i].pos,
                                                                                   old_vars->vars[old_i].ref_len,
                                                                                   BAM_CIGAR_STR[old_vars->vars[old_i].type],
                                                                                   old_vars->vars[old_i].n_alt_allele,
                                                                                   add_vars->vars[add_i].pos,
                                                                                   add_vars->vars[add_i].ref_len,
                                                                                   BAM_CIGAR_STR[add_vars->vars[add_i].type],
                                                                                   add_vars->vars[add_i].n_alt_allele);
                new_vars[new_i++] = old_vars->vars[old_i++];
                new_vars[new_i++] = add_vars->vars[add_i++];
            }
        }

    }
    // fprintf(stderr, "old_n: %d, add_n: %d, new_n: %d\n", old_vars->n, add_vars->n, new_i);
    for (; old_i < old_vars->n; ) new_vars[new_i++] = old_vars->vars[old_i++];
    for (; add_i < add_vars->n; ) new_vars[new_i++] = add_vars->vars[add_i++];

    free(old_vars->vars);
    old_vars->vars = new_vars; old_vars->n = new_i;
    return new_i;
}

// XXX ignore regions without fully-spanning reads
// re-alignment of specific regions within each haplotyp, and collect candidate variants
// 1. consensus calling within each haplotype
// 2. re-align clipping reads (if exist) to the consensus sequence
// 3. re-align unphased reads (if exist) to the consensus sequence
// 4. (*optional*) re-generate consensus sequence using all reads
// 5. collect candidate variants based on MSA of ref + cons1 + cons2
// 6. update cand_var & read_var_profile
// call variant for noisy_reg_i, return 0 if no variant is called
int collect_noisy_vars1(bam_chunk_t *chunk, const call_var_opt_t *opt, int noisy_reg_i) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t noisy_reg_beg = cr_start(noisy_regs, noisy_reg_i), noisy_reg_end = cr_end(noisy_regs, noisy_reg_i);
    uint8_t *ref_seq = NULL; int ref_seq_len = collect_reg_ref_bseq(chunk, &noisy_reg_beg, &noisy_reg_end, &ref_seq);
    hts_pos_t active_reg_beg = chunk->reg_beg, active_reg_end = chunk->reg_end;
    int max_noisy_reg_len = opt->max_noisy_reg_len, max_noisy_reg_reads = opt->max_noisy_reg_reads;
    if (noisy_reg_end - noisy_reg_beg + 1 > max_noisy_reg_len) {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped long region: %s:%ld-%ld %ld (>%d)\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, max_noisy_reg_len);
        free(ref_seq);
        return 0;
    }
    int *noisy_reads; int n_noisy_reads = collect_noisy_reg_reads1(chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, &noisy_reads);
    if (n_noisy_reads > max_noisy_reg_reads) {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped deep region: %s:%ld-%ld %ld %d reads\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        free(noisy_reads); free(ref_seq);
        return 0;
    }

    double realtime0 = realtime();

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "NoisyReg: chunk_reg: %s:%ld-%ld, reg: %s:%ld-%ld %ld (%d)\n", chunk->tname, chunk->reg_beg, chunk->reg_end, 
                        chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
    }
    // MSA and consensus calling
    int n_cons = 0; int *clu_n_seqs = (int*)calloc(2, sizeof(int)); int **clu_read_ids = (int**)malloc(2 * sizeof(int*));
    // for alignment of cons vs ref, read vs cons, read vs ref
    aln_str_t **aln_strs = (aln_str_t**)malloc(2 * sizeof(aln_str_t*));
    for (int i = 0; i < 2; ++i) {
        clu_read_ids[i] = NULL;
        aln_strs[i] = (aln_str_t*)malloc((1 + n_noisy_reads) * sizeof(aln_str_t));
        for (int j = 0; j < 1+ n_noisy_reads; ++j) {
            aln_strs[i][j].target_aln = NULL; aln_strs[i][j].query_aln = NULL; aln_strs[i][j].aln_len = 0;
        }
    }

    n_cons = collect_noisy_reg_aln_strs(opt, chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, n_noisy_reads, noisy_reads, ref_seq, ref_seq_len, // both cons_seqs and msa_seqs are separated for n_cons==2
                                        clu_n_seqs, clu_read_ids, aln_strs);
    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "n_cons: %d\n", n_cons);
    free(ref_seq);

    int n_noisy_vars = 0;
    if (n_cons == 0) {
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "Skipped region: %s:%ld-%ld %ld %d reads\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        n_noisy_vars = -1;
        goto collect_noisy_vars1_end;
    }

    cand_var_t *noisy_vars = NULL; int *noisy_var_cate=NULL; read_var_profile_t *noisy_rvp = NULL;
    // make variants and update cand_vars & read_var_profile
    n_noisy_vars = make_vars_from_msa_cons_aln(opt, chunk, n_noisy_reads, noisy_reads, noisy_reg_beg, noisy_reg_end, active_reg_beg, active_reg_end, 
                                               n_cons, clu_n_seqs, clu_read_ids, aln_strs, &noisy_vars, &noisy_var_cate, &noisy_rvp);
    int n_total_reads = 0;
    for (int i = 0; i < n_cons; ++i) n_total_reads += clu_n_seqs[i];
    merge_var_profile(chunk, n_total_reads, n_noisy_vars, noisy_vars, noisy_var_cate, noisy_rvp);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "NoisyReg: chunk_reg: %s:%ld-%ld, reg: %s:%ld-%ld %ld (%d)\t", chunk->tname, chunk->reg_beg, chunk->reg_end, 
                        chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        fprintf(stderr, "Real time: %.3f sec.\n", realtime() - realtime0);
    }
    free(noisy_var_cate);
collect_noisy_vars1_end:
    for (int i = 0; i < 2; ++i) {
        if (clu_read_ids[i] != NULL) free(clu_read_ids[i]);
        for (int j = 0; j < 1 + n_noisy_reads; ++j) {
            if (aln_strs[i][j].target_aln != NULL) free(aln_strs[i][j].target_aln);
        }
        free(aln_strs[i]);
    } 
    free(clu_n_seqs); free(clu_read_ids); free(aln_strs); free(noisy_reads);
    return n_noisy_vars;
}

// XXX
// sort noisy regions by length, full-cover & phased read counts, from short to long
int *sort_noisy_regs(bam_chunk_t *chunk) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int n_noisy_regs = noisy_regs->n_r;
    int *sorted_noisy_regs = (int*)malloc(n_noisy_regs * sizeof(int));
    int *noisy_reg_lens = (int*)malloc(n_noisy_regs * sizeof(int));
    int *noisy_reg_var_sizes = (int*)malloc(n_noisy_regs * sizeof(int));
    for (int i = 0; i < n_noisy_regs; ++i) {
        noisy_reg_lens[i] = cr_end(noisy_regs, i) - cr_start(noisy_regs, i);
        noisy_reg_var_sizes[i] = cr_label(noisy_regs, i);
        sorted_noisy_regs[i] = i;
    }
    for (int i = 0; i < n_noisy_regs; ++i) {
        for (int j = i+1; j < n_noisy_regs; ++j) {
            if (noisy_reg_var_sizes[sorted_noisy_regs[i]] > noisy_reg_var_sizes[sorted_noisy_regs[j]]) {
                int tmp = sorted_noisy_regs[i]; sorted_noisy_regs[i] = sorted_noisy_regs[j]; sorted_noisy_regs[j] = tmp;
            } else if (noisy_reg_var_sizes[sorted_noisy_regs[i]] == noisy_reg_var_sizes[sorted_noisy_regs[j]]) {
                if (noisy_reg_lens[sorted_noisy_regs[i]] > noisy_reg_lens[sorted_noisy_regs[j]]) {
                    int tmp = sorted_noisy_regs[i]; sorted_noisy_regs[i] = sorted_noisy_regs[j]; sorted_noisy_regs[j] = tmp;
                }
            }
        }
    }
    free(noisy_reg_lens); free(noisy_reg_var_sizes);
    return sorted_noisy_regs;
}

// remove noisy_regions that have low coverage, i.e., either low ratio or low absolute read count
void pre_process_noisy_regs(bam_chunk_t *chunk, call_var_opt_t *opt) {
    if (chunk->chunk_noisy_regs == NULL || chunk->chunk_noisy_regs->n_r == 0) return;
    cr_index(chunk->chunk_noisy_regs);
    // fprintf(stderr, "Before merge: %s:%ld-%ld %d noisy regions\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs->n_r);
    // for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
    //     fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
    // }
    chunk->chunk_noisy_regs = cr_merge(chunk->chunk_noisy_regs, -1);
    // fprintf(stderr, "After merge: %s:%ld-%ld %d noisy regions\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs->n_r);
    // for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
    //     fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
    // }

    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    hts_pos_t beg, end;
    uint8_t *skip_noisy_reg = (uint8_t*)calloc(noisy_regs->n_r, sizeof(uint8_t));
    int *noisy_reg_to_total_n_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    int *noisy_reg_to_noisy_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int r_i = ovlp_b[ovlp_i];
            noisy_reg_to_total_n_reads[r_i]++;
            // check if the read is noisy in the region XXX
            int noisy_reg_start = cr_start(noisy_regs, r_i), noisy_reg_end = cr_end(noisy_regs, r_i);
            int64_t noisy_digar_ovlp_n, *noisy_digar_ovlp_b = 0, noisy_digar_max_b = 0;
            noisy_digar_ovlp_n = cr_overlap(chunk->digars[i].noisy_regs, "cr", noisy_reg_start, noisy_reg_end, &noisy_digar_ovlp_b, &noisy_digar_max_b);
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy_reads[r_i]++;
            }
            free(noisy_digar_ovlp_b);
        }
    }
    int min_noisy_reg_reads = opt->min_noisy_reg_reads; //, max_noisy_reg_reads = opt->max_noisy_reg_reads;
    float min_noisy_reg_ratio = opt->min_noisy_reg_ratio;
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        int n_noisy_reg_reads = noisy_reg_to_noisy_reads[i];
        if (n_noisy_reg_reads < min_noisy_reg_reads // || noisy_reg_to_total_n_reads[i] > max_noisy_reg_reads
            || (float)n_noisy_reg_reads/noisy_reg_to_total_n_reads[i] < min_noisy_reg_ratio) {
            skip_noisy_reg[i] = 1;
            if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped region: %s:%d-%d %d noisy: %d, total: %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i)-cr_start(noisy_regs, i)+1, n_noisy_reg_reads, noisy_reg_to_total_n_reads[i]);
            n_skipped++;
        }
    }
    if (n_skipped > 0) {
        cgranges_t *new_noisy_regs = cr_init();
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (skip_noisy_reg[i]) continue;
            cr_add(new_noisy_regs, "cr", cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_label(noisy_regs, i));
        }
        cr_index(new_noisy_regs); cr_destroy(noisy_regs);
        chunk->chunk_noisy_regs = new_noisy_regs;
    }
    free(ovlp_b); free(skip_noisy_reg); free(noisy_reg_to_total_n_reads); free(noisy_reg_to_noisy_reads);
}

void collect_var_main(const call_var_pl_t *pl, bam_chunk_t *chunk) {
    // first round: easy-to-call SNPs (+indels)
    // 1. collect X/I/D sites from BAM
    collect_digars_from_bam(chunk, pl);

    // 2. merge all var sites from all reads, including low-depth ones, but not including nosiy-region ones
    var_site_t *var_sites = NULL; int n_var_sites;
    n_var_sites = collect_all_cand_var_sites(chunk, &var_sites);
    if (n_var_sites > 0) {
        // 3. collect all candidate variants, not including noisy-region ones
        // collect reference and alternative alleles for all var sites
        // all cand vars, including true/false germline/somatic variants
        // XXX for noisy/repeat regions, we need to carefully pick the candidate variants, based on supporting counts & re-alignments
        collect_cand_vars(chunk, n_var_sites, var_sites); free(var_sites);
    }

    // pre-process noisy regions (>5 XIDs in 100 win), not including low-complexity regions
    pre_process_noisy_regs(chunk, pl->opt);

    // 4. identify categories for cand vars, and update chunk_noisy_regs
    //    n_cand_vars & cand_vars: only include clean region vars
    classify_cand_vars(chunk, n_var_sites, pl->opt);

    if (chunk->n_cand_vars == 0 && (chunk->chunk_noisy_regs == NULL || chunk->chunk_noisy_regs->n_r == 0))
        return; // no variant to be called
    if (chunk->n_cand_vars > 0) { // XXX what if only have homozygous vars
        // 5. collect clean region read-wise var profiles
        chunk->read_var_profile = collect_read_var_profile(chunk);

        // process candidate variants in the following order
        // within each round, use previously obtained phasing/haplotype information as intialization
        //   1. LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL
        //   2. LONGCALLD_REP_HET_VAR + LONGCALLD_DENSE_REG_VAR + ALL VARS IN NOISY_REG
        //   3. LONGCALLD_CAND_HOM_VAR
        //   4. LONGCALLD_CAND_SOMATIC_VAR
        // 6. co-phasing and variant calling using clean region SNPs + indels
        assign_hap_based_on_het_vars_kmeans(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CAND_HOM_VAR, pl->opt);
    }
    // 7. iteratively call variants in noisy regions and variant/read phasing
    if (chunk->chunk_noisy_regs != NULL && chunk->chunk_noisy_regs->n_r > 0) {
        // sort noisy regions by type, length, full-cover & phased read counts, from short to long
        int *sorted_noisy_regs = sort_noisy_regs(chunk);
        // for each noisy region, call variants and update read_var_profile, then update phasing/haplotype information
        int *noisy_reg_is_done = (int*)calloc(chunk->chunk_noisy_regs->n_r, sizeof(int));
        while (1) {
            int new_region_is_done = 0, new_var = 0;
            for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
                int noisy_reg_i = sorted_noisy_regs[i];
                if (noisy_reg_is_done[noisy_reg_i]) continue;
                int ret = collect_noisy_vars1(chunk, pl->opt, noisy_reg_i);
                if (ret >= 0) {
                    noisy_reg_is_done[noisy_reg_i] = 1; new_region_is_done = 1;
                    if (ret > 0) {
                        new_var = 1;
                        // update phasing/haplotype information every time a noisy region is processed
                        // assign_hap_based_on_het_vars(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CAND_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR, pl->opt);
                        // assign_hap_based_on_het_vars_kmeans(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CAND_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR, pl->opt);
                    } else {
                        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "No var: %s:%d-%d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, noisy_reg_i), cr_end(chunk->chunk_noisy_regs, noisy_reg_i));
                    }
                } // unable to resolve the region
            }
            if (new_var) {
                assign_hap_based_on_het_vars_kmeans(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CAND_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR, pl->opt);
            }
            if (new_region_is_done == 0) break;
        }
        free(sorted_noisy_regs); free(noisy_reg_is_done);
    }
}

// stitch ii and ii+1
void stitch_var_main(call_var_step_t *step, bam_chunk_t *chunk, var_t *var, long ii) {
    call_var_pl_t *pl = step->pl;
    call_var_opt_t *opt = pl->opt;
    // ref_seq_t *ref_seq = pl->ref_seq;
    // if (ii == 0) { // extend phase set between two adjacent bam chunks
        // bam_chunk_t *prev_bam_chunk = step->chunks+ii-1;
        // chunk->flip_hap = prev_chunk->flip_hap ^ flip_variant_hap(prev_bam_chunk, bam_chunk);
    // }

    // generate variant-related information, e.g., GT, DP, etc.
    // merge variants based on ref_pos if needed
    var_t *_vars;
    int n_vars = make_variants(opt, chunk, &_vars);
    if (n_vars > 0) {
        merge_vars(var, _vars);
        free(_vars->vars); free(_vars);
    }
}