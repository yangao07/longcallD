#include <stdio.h>
#include <stdlib.h>
#include "collect_var.h"
#include "call_var_main.h"
#include "utils.h"
#include "bam_utils.h"
#include "align.h"
#include "math_utils.h"
#include "sdust.h"
#include "seq.h"
#include "assign_hap.h"
#include "vcf_utils.h"

extern int LONGCALLD_VERBOSE;

// only 1 variant is stored in each cand_var_t, i.e., SNP/INS/DEL
// each read with be 
//    0: ref ,or 1: alt, or
//   -1: other alt. mostly sequencing errors, should be considered as non-informative during haplotype assignment, i.e. -1!=-1
cand_var_t *init_cand_vars_based_on_sites(int n_var_sites, var_site_t *var_sites) {
    cand_var_t *cand_vars = (cand_var_t*)malloc(n_var_sites * sizeof(cand_var_t));
    for (int i = 0; i < n_var_sites; ++i) {
        memset(cand_vars+i, 0, sizeof(cand_var_t));
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
        // somatic/TSD/TE information
        cand_vars[i].somatic_aux_info = NULL; // unset
        cand_vars[i].te_seq_i = -1; // unset
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
    if (cand_vars->somatic_aux_info != NULL) free_somatic_var_aux_info(cand_vars->somatic_aux_info);
    if (cand_vars->tsd_len > 0) free(cand_vars->tsd_seq);
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

int fuzzy_ovlp_var_site(var_site_t *var1, var_site_t *var2) {
    int beg1 = var1->pos; int end1 = var1->pos + var1->ref_len;
    int beg2 = var2->pos; int end2 = var2->pos + var2->ref_len;
    if (var1->var_type == BAM_CINS && var2->var_type == BAM_CINS) {
        int min_alt_len = MIN_OF_TWO(var1->alt_len, var2->alt_len);
        if (min_alt_len >= abs(beg1-beg2)) return 1;
    } else if (var1->var_type == BAM_CDEL && var2->var_type == BAM_CDEL) {
        int min_end = MIN_OF_TWO(end1, end2);
        int max_beg = MAX_OF_TWO(beg1, beg2);
        if (min_end >= max_beg) return 1;
    }
    return 0;
}

// int fuzzy_comp_seq(const call_var_opt_t *opt, uint8_t *seq1, int len1, uint8_t *seq2, int len2) {
//     int n_eq = 0, n_xid = 0;
//     wfa_heuristic_aln(seq1, len1, seq2, len2, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, &n_eq, &n_xid);
//     int max_len = MAX_OF_TWO(len1, len2);
//     if (n_eq >= max_len * 0.9) return 0;
//     else return 1;
// }

// allow cyclic permutation match
int vntr_fuzzy_comp_seq(const call_var_opt_t *opt, uint8_t *seq1, int len1, uint8_t *seq2, int len2) {
    int n_eq=0;
    uint8_t *long_seq, *short_seq; int long_len, short_len;
    if (len1 > len2) {
        long_seq = (uint8_t*)malloc(len1*2*sizeof(uint8_t));
        for (int i = 0; i < len1; ++i) {
            long_seq[i] = seq1[i];
            long_seq[i+len1] = seq1[i];
        }
        long_len = len1*2;
        short_seq = seq2; short_len = len2;
    } else {
        long_seq = (uint8_t*)malloc(len2*2*sizeof(uint8_t));
        for (int i = 0; i < len2; ++i) {
            long_seq[i] = seq2[i];
            long_seq[i+len2] = seq2[i];
        }
        long_len = len2*2;
        short_seq = seq1; short_len = len1;
    }
    uint8_t *seq1_aln_str, *seq2_aln_str; int aln_len;
    wfa_end2end_aln(long_seq, long_len, short_seq, short_len, opt->gap_aln, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, NULL, NULL, &seq1_aln_str, &seq2_aln_str, &aln_len);
    for (int i = 0; i < aln_len; ++i) {
        if (seq1_aln_str[i] != 5 && seq2_aln_str[i] == seq1_aln_str[i]) n_eq++;
    }
    free(seq1_aln_str); free(long_seq);
    if (n_eq >= short_len * 0.8) return 0;
    else return 1;
}

// comp two INS/DEL variants that fuzzy_ovlp with each other
// return 0: fuzzy match; non-0: not fuzzy match
int fuzzy_comp_var_site(const call_var_opt_t *opt, var_site_t *var1, var_site_t *var2) {
    if (var1->var_type == BAM_CDEL && var2->var_type == BAM_CDEL) {
        int min_len = MIN_OF_TWO(var1->ref_len, var2->ref_len);
        int max_len = MAX_OF_TWO(var1->ref_len, var2->ref_len);
        if (min_len >= max_len * 0.8) return 0;
    } else if (var1->var_type == BAM_CINS && var2->var_type == BAM_CINS) {
        int min_len = MIN_OF_TWO(var1->alt_len, var2->alt_len);
        int max_len = MAX_OF_TWO(var1->alt_len, var2->alt_len);
        if (min_len >= max_len * 0.8 && (vntr_fuzzy_comp_seq(opt, var1->alt_seq, var1->alt_len, var2->alt_seq, var2->alt_len) == 0)) return 0;
    }
    return exact_comp_var_site(opt, var1, var2);
}

int fuzzy_comp_ovlp_var_site(const call_var_opt_t *opt, var_site_t *var1, var_site_t *var2, int *is_ovlp) {
    if ((var1->var_type == BAM_CINS && var2->var_type == BAM_CINS && var1->alt_len >= opt->min_sv_len && var2->alt_len >= opt->min_sv_len)
     || (var1->var_type == BAM_CDEL && var2->var_type == BAM_CDEL && var1->ref_len >= opt->min_sv_len && var2->ref_len >= opt->min_sv_len)) {
        *is_ovlp = fuzzy_ovlp_var_site(var1, var2);
        if (*is_ovlp == 1) {
            int ret = fuzzy_comp_var_site(opt, var1, var2);
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "fuzzy_comp_ovlp_var_site: %" PRId64 " %d-%c-%d vs %" PRId64 " %d-%c-%d => %d\n", 
                        var1->pos, var1->ref_len, BAM_CIGAR_STR[var1->var_type], var1->alt_len,
                        var2->pos, var2->ref_len, BAM_CIGAR_STR[var2->var_type], var2->alt_len, ret);
            }
            return ret;
        } else return exact_comp_var_site(opt, var1, var2);
    }
    *is_ovlp = ovlp_var_site(var1, var2);
    return exact_comp_var_site(opt, var1, var2);
}

int fuzzy_comp_cand_var(const call_var_opt_t *opt, cand_var_t *var1, cand_var_t *var2) {
    var_site_t var_site1 = make_var_site_from_cand_var(var1);
    var_site_t var_site2 = make_var_site_from_cand_var(var2);

    if ((var_site1.var_type == BAM_CINS && var_site2.var_type == BAM_CINS && var_site1.alt_len >= opt->min_sv_len && var_site2.alt_len >= opt->min_sv_len) ||
        (var_site1.var_type == BAM_CDEL && var_site2.var_type == BAM_CDEL && var_site1.ref_len >= opt->min_sv_len && var_site2.ref_len >= opt->min_sv_len)) {
        if (fuzzy_ovlp_var_site(&var_site1, &var_site2))
            return fuzzy_comp_var_site(opt, &var_site1, &var_site2);
    } 
    return exact_comp_var_site(opt, &var_site1, &var_site2);
}

int fuzzy_comp_ins_var_low_complexity(const call_var_opt_t *opt, var_site_t *large_var, var_site_t *small_var) {
    if (large_var->var_type == BAM_CINS && small_var->var_type == BAM_CINS) {
        int large_ins_len = large_var->alt_len, small_ins_len = small_var->alt_len;
        if (large_ins_len < small_ins_len) return 0; // not low_comp_ins
        // if ((vntr_fuzzy_comp_seq(opt, large_var->alt_seq, large_var->alt_len, small_var->alt_seq, small_var->alt_len) == 0)) return 0;
        int diff_len = 0; uint8_t *diff_seq = 0;
        diff_len = wfa_collect_diff_ins_seq(opt, large_var->alt_seq, large_var->alt_len, small_var->alt_seq, small_var->alt_len, &diff_seq);
        uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
        r = sdust(0, diff_seq, diff_len, T, W, &n);
        int low_comp_len = 0;
        for (int i = 0; i < n; ++i) {
            // fprintf(stderr, "low_comp_cr: %s %d-%d\n", chunk->tname, chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1);
            low_comp_len += ((int)r[i]- (int)(r[i]>>32));
        }
        free(r); if (diff_len > 0) free(diff_seq);
        if (low_comp_len > diff_len * 0.8) return 1;
    }
    return 0;
}

// void free_noisy_cand_vars1(cand_var_t *cand_vars) {
//     if (cand_vars->alle_covs != NULL) free(cand_vars->alle_covs); 
//     if (cand_vars->strand_to_alle_covs != NULL) {
//         for (int j = 0; j < 2; ++j) {
//             free(cand_vars->strand_to_alle_covs[j]);
//         } free(cand_vars->strand_to_alle_covs);
//     }
//     if (cand_vars->alt_seq != NULL) free(cand_vars->alt_seq);
//     if (cand_vars->hap_to_alle_profile != NULL) {
//         for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
//             free(cand_vars->hap_to_alle_profile[j]);
//         } free(cand_vars->hap_to_alle_profile);
//     }
//     if (cand_vars->hap_to_cons_alle != NULL) free(cand_vars->hap_to_cons_alle);
// }

// void free_noisy_cand_vars(cand_var_t *cand_vars, int m) {
//     for (int i = 0; i < m; ++i) free_noisy_cand_vars1(cand_vars+i);
//     free(cand_vars);
// }

int collect_cand_vars(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_var_sites, var_site_t *var_sites) {
    chunk->cand_vars = init_cand_vars_based_on_sites(n_var_sites, var_sites);
    chunk->n_cand_vars = n_var_sites;
    cand_var_t *cand_vars = chunk->cand_vars;
    // 2nd pass: update snp_sites, calculate the depth and allele frequency of each site
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        start_var_i = update_cand_vars_from_digar(opt, chunk, chunk->digars+i, n_var_sites, var_sites, start_var_i, cand_vars);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Collect %d candidate variants from %d reads\n", n_var_sites, chunk->n_reads);
        for (int i = 0; i < n_var_sites; ++i) {
            fprintf(stderr, "CandVar: %" PRId64 "\t", var_sites[i].pos);
            fprintf(stderr, "Type: %c\t", BAM_CIGAR_STR[var_sites[i].var_type]);
            fprintf(stderr, "RefLen: %d\t", var_sites[i].ref_len);
            fprintf(stderr, "Depth: %d\t", cand_vars[i].total_cov);
            if (cand_vars[i].var_type == BAM_CDEL) fprintf(stderr, "DEL");
            else {
                for (int k = 0; k < cand_vars[i].alt_len; ++k)
                    fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seq[k]]);
            }
            fprintf(stderr, ": ");
            fprintf(stderr, "%d: %d\t", cand_vars[i].alt_len, cand_vars[i].alle_covs[1]);
            fprintf(stderr, "LowQual-Depth: %d\n", cand_vars[i].low_qual_cov);
        }
    }
    return 0;
}

int var_is_strand_bias(cand_var_t *var, const call_var_opt_t *opt) {
    int for_alt_cov = var->strand_to_alle_covs[0][1]; // forward strand alt coverage
    int rev_alt_cov = var->strand_to_alle_covs[1][1]; // reverse strand alt coverage
    int for_ref_cov = var->strand_to_alle_covs[0][0]; // forward strand ref coverage
    int rev_ref_cov = var->strand_to_alle_covs[1][0]; // reverse strand ref coverage
    int expected_alt_cov = (for_alt_cov + rev_alt_cov) / 2; // expected alt coverage
    if (expected_alt_cov == 0) return 0; // no alt coverage, no strand bias
    // check if strand bias is significant
    float fisher_p = fisher_exact_test(for_alt_cov, rev_alt_cov, expected_alt_cov, expected_alt_cov, opt);
    // float fisher_p = fisher_exact_test(for_alt_cov, rev_alt_cov, for_ref_cov, rev_ref_cov, opt);
    if (fisher_p < opt->strand_bias_pval) {
        // fprintf(stderr, "ONT-StrandBias: %" PRId64 " for_alt_cov=%d, rev_alt_cov=%d, for_ref_cov=%d, rev_ref_cov=%d, fisher p-value=%.5f\n",
                // var->pos, for_alt_cov, rev_alt_cov, for_ref_cov, rev_ref_cov, fisher_p);
        return 1; // significant strand bias
    } else return 0; // no significant strand bias
    // strand bias: 1) >=3 vs 0, 2) >= 3 folds
    // int min_fold = 3;
    // if (var->alle_covs[1] < min_fold) return 0;
    // else if (var->strand_to_alle_covs[0][1] >= 2 && var->strand_to_alle_covs[1][1] >= 2) return 0;
    // else if (var->strand_to_alle_covs[0][1] == 0 && var->strand_to_alle_covs[1][1] >= min_fold) return 1;
    // else if (var->strand_to_alle_covs[1][1] == 0 && var->strand_to_alle_covs[0][1] >= min_fold) return 1;
    // else if (var->strand_to_alle_covs[0][1] > 0 && var->strand_to_alle_covs[1][1] > 0) {
    //     if (var->strand_to_alle_covs[0][1] >= min_fold * var->strand_to_alle_covs[1][1]) return 1;
    //     else if (var->strand_to_alle_covs[1][1] >= min_fold * var->strand_to_alle_covs[0][1]) return 1;
    // }
    // return 0;
}

// XXX check the sequence around the variant site, not just left/right side
// filter: 
//   1) usable X/= over total non-low-qual depth >= threshold
//   2) total non-low-qual depth / total depth >= threshold

// 1-homopolymer: AAA 1*3
// 2-homopolymer: CGCGCG 2*3
// N-homopolymer: ACGACGACG N*3
int var_is_homopolymer(const call_var_opt_t *opt, char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    // fprintf(stderr, "%d-%d: %s\n", ref_beg, ref_end, ref_seq);
    hts_pos_t start_pos, end_pos; // = var->pos, end_pos = var->pos;
    int xid = opt->noisy_reg_max_xgaps; // XXX xid
    if (var->var_type == BAM_CDIFF) {
        start_pos = var->pos-1; end_pos = var->pos+1;
    } else if (var->var_type == BAM_CINS) {
        if (var->alt_len > xid) return 0; // skip checking for INDEL > 5bp
        start_pos = var->pos-1; end_pos = var->pos;
    } else { // DEL
        if (var->ref_len > xid) return 0; // skip checking for INDEL > 5bp
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
        // fprintf(stderr, "%" PRIi64 ", %" PRIi64 ", %" PRIi64 ", %" PRIi64 "\n", start_pos, end_pos, ref_beg, ref_end);
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
int var_is_repeat_region(const call_var_opt_t *opt, char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    hts_pos_t pos = var->pos; int ref_len = var->ref_len; int var_type = var->var_type;
    uint8_t *ref_bseq, *alt_bseq; int len = 0, xid = opt->noisy_reg_max_xgaps; // XXX xid
    int is_repeat = 1;
    if (var_type == BAM_CDEL) { // [pos, pos+del_len*N] vs [pos+del_len, pos+del_len*(N+1)]
        // nst_nt4_table['A'] -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4
        int del_len = ref_len;
        if (del_len > xid) return 0; // skip checking for INDEL > 5bp
        len = del_len * 3; // see if del seq is 3-fold repeat
        if (pos < ref_beg || pos+del_len+len >= ref_end) {
            // fprintf(stderr, "DelLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", del_len, ref_len, pos, ref_beg, ref_end);
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
        if (ins_len > xid) return 0; // skip checking for INDEL > 5bp
        len = ins_len * 3; // see if ins seq is 3-fold repeat
        if (pos < ref_beg || pos+len >= ref_end) {
            // fprintf(stderr, "InsLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", ins_len, ref_len, pos, ref_beg, ref_end);
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
                      cand_var_t *var, int min_dp_thres, int min_alt_dp_thres,
                      double min_af_thres, double max_af_thres, double min_noisy_reg_ratio) {
    // total depth < min_dp or total alt depth < min_alt_dp: (skip)
    if (var->total_cov + var->low_qual_cov < min_dp_thres) return LONGCALLD_LOW_COV_VAR;
    // low-qual depth / total depth > min_noisy_reg_ratio: too many low-qual bases, eg., dense X/gap region
    double alt_af=0.0; int alt_dp=0;
    alt_dp = var->alle_covs[1]; alt_af = (double) alt_dp / var->total_cov;
    if (alt_dp < min_alt_dp_thres) return LONGCALLD_LOW_COV_VAR; // skip var with low alt depth 
    if (opt->is_ont && var_is_strand_bias(var, opt)) return LONGCALLD_STRAND_BIAS_VAR; // skip var with strand bias
    // remaining categories: CLEAN_HET(SNP/INDEL), CLEAN_HOM(SNP/INDEL), NOISY_REG_VAR (others)
    if (alt_af < min_af_thres) return LONGCALLD_LOW_AF_VAR; // needs to further check
    // if ((double) var->low_qual_cov / (var->low_qual_cov + var->total_cov) >= min_noisy_reg_ratio) return LONGCALLD_NOISY_REG_VAR; // too many low-qual bases, needs to further check
    if (alt_af > max_af_thres) return LONGCALLD_CLEAN_HOM_VAR; // unlikely germline het., likely hom or somatic, require full phasing info
    // snps & indels in homo/repeat regions
    if ((var->var_type == BAM_CINS || var->var_type == BAM_CDEL) && 
        (var_is_homopolymer(opt, ref_seq, ref_beg, ref_end, var) || var_is_repeat_region(opt, ref_seq, ref_beg, ref_end, var)))
        return LONGCALLD_REP_HET_VAR; // require basic phasing info, MSA, provide additional phasing info
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
    // copy tsd-related info
    if (to_var->tsd_len > 0) free(to_var->tsd_seq);
    if (from_var->tsd_len > 0) {
        to_var->tsd_seq = (uint8_t*)malloc(from_var->tsd_len * sizeof(uint8_t));
        memcpy(to_var->tsd_seq, from_var->tsd_seq, from_var->tsd_len);
    }
    to_var->checked_tsd = from_var->checked_tsd;
    to_var->tsd_len = from_var->tsd_len; to_var->polya_len = from_var->polya_len;
    to_var->tsd_pos1 = from_var->tsd_pos1; to_var->tsd_pos2 = from_var->tsd_pos2;
    to_var->te_seq_i = from_var->te_seq_i; to_var->te_is_rev = from_var->te_is_rev;
}

void low_comp_cr_start_end(cgranges_t *low_comp_cr, int32_t start, int32_t end, int32_t *new_start, int32_t *new_end) {
    *new_start = start; *new_end = end;
    if (low_comp_cr == NULL || low_comp_cr->n_r == 0) return;
    int64_t *low_comp_b = 0; int64_t low_comp_n = 0, max_low_comp_n = 0;
    low_comp_n = cr_overlap(low_comp_cr, "cr", start-1, end, &low_comp_b, &max_low_comp_n);
    for (int64_t j = 0; j < low_comp_n; ++j) {
        int32_t start = cr_start(low_comp_cr, low_comp_b[j])+1;
        int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
        if (start < *new_start) *new_start = start;
        if (end > *new_end) *new_end = end;
    }
    if (low_comp_b) free(low_comp_b);
}

// XXX no vars spanning the boundaries of noisy regions
static void collect_noisy_reg_start_end(bam_chunk_t *chunk, const call_var_opt_t *opt, int *var_i_to_cate, int32_t *start, int32_t *end) {
    int *max_left_var_i = (int32_t*)malloc(chunk->chunk_noisy_regs->n_r * sizeof(int32_t));
    int *min_right_var_i = (int32_t*)malloc(chunk->chunk_noisy_regs->n_r * sizeof(int32_t));
    for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
        max_left_var_i[i] = -1; min_right_var_i[i] = -1;
    }
    int n_noisy_regs = chunk->chunk_noisy_regs->n_r;
    for (int reg_i = 0, var_i = 0; reg_i < n_noisy_regs && var_i < chunk->n_cand_vars;) {
        if (var_i_to_cate[var_i] == LONGCALLD_LOW_COV_VAR) {
            var_i++; continue;
        }
        cand_var_t *var = chunk->cand_vars + var_i;
        int32_t var_start = var->pos, var_end = var->pos + var->ref_len - 1;
        int32_t reg_start = cr_start(chunk->chunk_noisy_regs, reg_i)+1, reg_end = cr_end(chunk->chunk_noisy_regs, reg_i);
        if (var_start > reg_end) {
            if (min_right_var_i[reg_i] == -1) min_right_var_i[reg_i] = var_i;
            reg_i++;
        } else if (var_end < reg_start) {
            max_left_var_i[reg_i] = var_i;
            var_i++;
        } else { // var and reg overlap
            var_i++;
        }
    }
    // update noisy regions: at least flank_len bp without any vars on both ends
    int flank_len = opt->noisy_reg_flank_len;
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        if (max_left_var_i[reg_i] == -1) max_left_var_i[reg_i] = MIN_OF_TWO(chunk->n_cand_vars-1, 0);
        if (min_right_var_i[reg_i] == -1) min_right_var_i[reg_i] = MAX_OF_TWO(0, chunk->n_cand_vars-1);
        int32_t ori_reg_start = cr_start(chunk->chunk_noisy_regs, reg_i)+1, ori_reg_end = cr_end(chunk->chunk_noisy_regs, reg_i);
        // start
        int32_t cur_start, cur_end;
        cur_start = ori_reg_start-flank_len;
        // iterate from the rightmost var to the leftmost var, check if _cur_start_ meet the criteria
        for (int var_i = max_left_var_i[reg_i]; var_i >= 0; --var_i) {
            if (var_i_to_cate[var_i] == LONGCALLD_LOW_COV_VAR) continue;
            cand_var_t *var = chunk->cand_vars+var_i;
            int32_t var_start = var->pos, var_end = var->pos + var->ref_len - 1;
            if (var_end < cur_start) break;
            else cur_start = var_start-flank_len;
        }
        start[reg_i] = cur_start;
        // end
        cur_end = ori_reg_end + flank_len;
        for (int var_i = min_right_var_i[reg_i]; var_i < chunk->n_cand_vars; ++var_i) {
            if (var_i_to_cate[var_i] == LONGCALLD_LOW_COV_VAR) continue;
            cand_var_t *var = chunk->cand_vars+var_i;
            int32_t var_start = var->pos, var_end = var->pos + var->ref_len - 1;
            if (var_start > cur_end) break;
            else cur_end = var_end + flank_len;
        }
        end[reg_i] = cur_end;
    }
    free(max_left_var_i); free(min_right_var_i);
}

// only merge same type of noisy regions, i.e., DEL/INS, avoid merging large DELs with other small DELs/INSs
// remove noisy_regions that have low coverage, i.e., either low ratio or low absolute read count
void pre_process_noisy_regs(bam_chunk_t *chunk, call_var_opt_t *opt) {
    if (chunk->chunk_noisy_regs == NULL || chunk->chunk_noisy_regs->n_r == 0) return;
    cr_index(chunk->chunk_noisy_regs);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Before merge: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " noisy regions\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs->n_r);
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
    }
    chunk->chunk_noisy_regs = cr_merge(chunk->chunk_noisy_regs, -1, opt->noisy_reg_merge_dis, opt->min_sv_len);
    if (LONGCALLD_VERBOSE >= 3) {
        fprintf(stderr, "After merge: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " noisy regions\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs->n_r);
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
    }

    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    hts_pos_t beg, end;
    uint8_t *skip_noisy_reg = (uint8_t*)calloc(noisy_regs->n_r, sizeof(uint8_t));
    int *noisy_reg_to_total_n_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    int *noisy_reg_to_noisy_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    int *noisy_reg_to_label = (int*)calloc(noisy_regs->n_r, sizeof(int));
    
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int r_i = ovlp_b[ovlp_i];
            noisy_reg_to_total_n_reads[r_i]++;
            // check if the read is noisy in the region XXX
            int noisy_reg_start = cr_start(noisy_regs, r_i)+1, noisy_reg_end = cr_end(noisy_regs, r_i);
            noisy_reg_to_label[r_i] = cr_label(noisy_regs, r_i);
            int64_t noisy_digar_ovlp_n, *noisy_digar_ovlp_b = 0, noisy_digar_max_b = 0;
            noisy_digar_ovlp_n = cr_overlap(chunk->digars[i].noisy_regs, "cr", noisy_reg_start-1, noisy_reg_end, &noisy_digar_ovlp_b, &noisy_digar_max_b);
            if (noisy_digar_ovlp_n > 0) {
                noisy_reg_to_noisy_reads[r_i]++;
            }
            free(noisy_digar_ovlp_b);
        }
    }
    int min_noisy_reg_reads = opt->min_noisy_reg_reads; //, max_noisy_reg_reads = opt->max_noisy_reg_reads;
    float min_noisy_reg_ratio = opt->min_af; // opt->min_noisy_reg_ratio;
    int n_skipped = 0;
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        int n_noisy_reg_reads = noisy_reg_to_noisy_reads[i];
        // fprintf(stderr, "NoisyRegion: %s:%d-%d %d noisy: %d, total: %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i)-cr_start(noisy_regs, i)+1, n_noisy_reg_reads, noisy_reg_to_total_n_reads[i]);
        if (n_noisy_reg_reads < min_noisy_reg_reads // || noisy_reg_to_total_n_reads[i] > max_noisy_reg_reads
            || (float)n_noisy_reg_reads/noisy_reg_to_total_n_reads[i] < min_noisy_reg_ratio) {
            skip_noisy_reg[i] = 1;
            if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped region: %s:%d-%d %d noisy: %d, total: %d (%.3f < %.3f)\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i)-cr_start(noisy_regs, i)+1,
                                                        n_noisy_reg_reads, noisy_reg_to_total_n_reads[i], (float)n_noisy_reg_reads/noisy_reg_to_total_n_reads[i], min_noisy_reg_ratio);
            // potential somatic variant
            // if (opt->out_somatic)  // add to somatic_cr
            //     cr_add(somatic_cr, "cr", cr_start(noisy_regs, i), cr_end(noisy_regs, i), noisy_reg_to_label[i]);
            // if (noisy_reg_to_label[i] >= opt->min_sv_len) { 
            //     fprintf(stderr, "Skipped region: %s:%d-%d %d noisy: %d, total: %d, SV: %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i)-cr_start(noisy_regs, i)+1, 
            //                     n_noisy_reg_reads, noisy_reg_to_total_n_reads[i], noisy_reg_to_label[i]);
            // }
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
    if (LONGCALLD_VERBOSE >= 3) {
        fprintf(stderr, "After filtering: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " noisy regions\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs->n_r);
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
    }

    free(ovlp_b); free(skip_noisy_reg);
    free(noisy_reg_to_total_n_reads); free(noisy_reg_to_noisy_reads); free(noisy_reg_to_label);
}

// extend noisy regions by flanking length without any vars, up to noisy_reg_flank_len
// >=1 bp between noisy regions and flanking clean vars to avoid overlapping
void post_process_noisy_regs(bam_chunk_t *chunk, const call_var_opt_t *opt, int *var_i_to_cate) {
    cgranges_t *chunk_noisy_regs = chunk->chunk_noisy_regs;
    int n_noisy_regs = chunk_noisy_regs->n_r;
    cgranges_t *noisy_regs = cr_init();
    int32_t *noisy_reg_start = (int32_t*)malloc(n_noisy_regs * sizeof(int32_t));
    int32_t *noisy_reg_end = (int32_t*)malloc(n_noisy_regs * sizeof(int32_t));
    collect_noisy_reg_start_end(chunk, opt, var_i_to_cate, noisy_reg_start, noisy_reg_end);
    for (int reg_i = 0; reg_i < n_noisy_regs; ++reg_i) {
        cr_add(noisy_regs, "cr", noisy_reg_start[reg_i], noisy_reg_end[reg_i], cr_label(chunk->chunk_noisy_regs, reg_i));
    }
    cr_index(noisy_regs); cr_destroy(chunk->chunk_noisy_regs);
    chunk->chunk_noisy_regs = cr_merge(noisy_regs, 0, -1, -1);
    free(noisy_reg_start); free(noisy_reg_end);
}

float var_noisy_reads_ratio(bam_chunk_t *chunk, hts_pos_t var_start, hts_pos_t var_end) {
    int total_n_reads = 0;
    int noisy_reads = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        hts_pos_t beg = chunk->digars[i].beg, end = chunk->digars[i].end;
        if (beg > var_end || end < var_start) continue; // no overlap
        total_n_reads++;
        // check if there is mismatch/ins/del digars in the region
        for (int j = 0; j < chunk->digars[i].n_digar; ++j) {
            digar1_t *d1 = chunk->digars[i].digars+j;
            if (d1->type != BAM_CDIFF && d1->type != BAM_CINS && d1->type != BAM_CDEL) continue; // only check X/I/D
            hts_pos_t digar_pos = d1->pos;
            hts_pos_t digar_end = d1->pos;
            if (d1->type == BAM_CDIFF || d1->type == BAM_CDEL) digar_end += (d1->len-1);
            if (digar_end < var_start) continue;
            else if (digar_pos > var_end) break; // no overlap
            noisy_reads++; break;
        }
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "var_noisy_ratio: %s:%" PRId64 "-%" PRId64 " total_n_reads: %d, noisy_reads: %d, ratio: %.3f\n", 
                chunk->tname, var_start, var_end, total_n_reads, noisy_reads, (float) noisy_reads / (total_n_reads+0.0));
    }
    if (total_n_reads == 0) return 0.0;
    else return (float) noisy_reads / (total_n_reads+0.0);
}

// add var to the regions, extend the region based on low_comp_cr if needed
int cr_add_var_cr(const call_var_opt_t *opt, bam_chunk_t *chunk, cgranges_t *var_cr, cgranges_t *low_comp_cr, cand_var_t *var, int check_noisy_reads_ratio) {
    hts_pos_t var_start, var_end;
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
            int32_t start = cr_start(low_comp_cr, low_comp_b[j])+1;
            int32_t end = cr_end(low_comp_cr, low_comp_b[j]);
            if (start < var_start) var_start = start;
            if (end > var_end) var_end = end;
        }
        free(low_comp_b);
    }
    if (check_noisy_reads_ratio == 0 || (var_noisy_reads_ratio(chunk, var_start, var_end) >= opt->min_af))
        cr_add(var_cr, "cr", var_start-1, var_end, 1); // XXX set merge_win as 1, previously: niosy_reg_flank_len
    return 0;
}

int var_is_low_comp(bam_chunk_t *chunk, const call_var_opt_t *opt, cand_var_t *var) {
    // check if the variant is in low-complexity regions
    int64_t low_comp_n = 0, max_low_comp_n = 0;
    if (chunk->low_comp_cr == NULL || chunk->low_comp_cr->n_r == 0) return 0;
    int64_t *low_comp_b = 0; 
    low_comp_n = cr_is_contained(chunk->low_comp_cr, "cr", var->pos-1, var->pos+var->ref_len-1, &low_comp_b, &max_low_comp_n);
    free(low_comp_b);
    if (low_comp_n == 0) {
        char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
        if (var_is_homopolymer(opt, ref_seq, ref_beg, ref_end, var) || var_is_repeat_region(opt, ref_seq, ref_beg, ref_end, var)) return 1;
        return 0;
    } else return 1;
}

// XXX to do
int var_is_ovlp_others(bam_chunk_t *chunk, int var_i) {
    int ovlp = 0;
    int *var_i_to_cate = chunk->var_i_to_cate;
    for (int i = var_i - 1; i >= 0; --i) {
        int var_cate = var_i_to_cate[i];
        if (var_cate == LONGCALLD_LOW_COV_VAR || var_cate == LONGCALLD_NON_VAR || var_cate == LONGCALLD_STRAND_BIAS_VAR) continue;
        cand_var_t *var1 = chunk->cand_vars + var_i;


    }
    for (int i = var_i + 1; i < chunk->n_cand_vars; ++i) {
        int var_cate = var_i_to_cate[i];
        if (var_cate == LONGCALLD_LOW_COV_VAR || var_cate == LONGCALLD_NON_VAR || var_cate == LONGCALLD_STRAND_BIAS_VAR) continue;
    }


    return ovlp;
}

// SNV/DEL/INS: n_reads >= 2
// TE-INS:  n_reads >= 1
// no somatic small indel
int var_is_cand_somatic(bam_chunk_t *chunk, const call_var_opt_t *opt, cand_var_t *var) {
    // if (var_is_low_comp(chunk, opt, var)) return 0;
    if (var->var_type == BAM_CDIFF) {
        if (var->alle_covs[1] >= opt->min_somatic_alt_dp) return 1;
    } else { // INDEL
        int tsd_len = collect_te_info_from_var(opt, chunk, var);
        if (var->alt_len >= opt->min_sv_len || var->ref_len >= opt->min_sv_len) {
            if (var->alle_covs[1] >= opt->min_somatic_alt_dp) return 1;
            else if (var->var_type == BAM_CINS && var->alle_covs[1] >= opt->min_somatic_te_dp && tsd_len > 0) {
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "CandTE-SV: %s:%" PRId64 " %d-%c-%d %d vs %d\n", chunk->tname, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len, var->alle_covs[0], var->alle_covs[1]);
                return 1;
            }
        }
    }
    return 0;
}

// merge somatic vars in clean regions
// 1. merge exising somatic vars if they match
// 2. will collect somatic_var_read_profile in later stage, so the alle_covs may not be correct right now
void merge_clean_somatic_vars(const call_var_opt_t *opt, cand_var_t *cand_vars, int *var_i_to_cate, int n_cand_vars) {
    int n_somatic_vars = 0;
    int *somatic_var_i = (int*)malloc(n_cand_vars * sizeof(int));
    for (int i = 0; i < n_cand_vars; ++i) {
        if (var_i_to_cate[i] == LONGCALLD_CAND_SOMATIC_VAR) somatic_var_i[n_somatic_vars++] = i;
    }
    if (n_somatic_vars <= 0) {
        free(somatic_var_i); return;
    }
    // 1. merge existing somatic vars if they fuzzy match: only collect merged somatic var site
    // int *updated_somatic_var_i = (int*)calloc(n_somatic_vars, sizeof(int));
    for (int i = 0; i < n_somatic_vars-1; ++i) {
        int var_i = somatic_var_i[i];
        if (var_i_to_cate[var_i] != LONGCALLD_CAND_SOMATIC_VAR) continue;
        cand_var_t *var1 = cand_vars + var_i;
        for (int j = i+1; j < n_somatic_vars; ++j) {
            int var_j = somatic_var_i[j];
            if (var_i_to_cate[var_i] != LONGCALLD_CAND_SOMATIC_VAR) break;
            if (var_i_to_cate[var_j] != LONGCALLD_CAND_SOMATIC_VAR) continue;
            cand_var_t *var2 = cand_vars + var_j;
            // check if two vars: 1) close enough, 2) similar enough (SV INSs)
            // int ret = fuzzy_comp_cand_var(opt, cand_vars+var_i, cand_vars+var_j);
            // !!! var_i and var_j are no way to be exact match, otherwise they will be merged already !!!
            int comp = -1;
            if (var1->var_type == BAM_CINS && var2->var_type == BAM_CINS) { // two INSs
                int min_alt_len = MIN_OF_TWO(var1->alt_len, var2->alt_len);
                if (min_alt_len >= abs(var1->pos - var2->pos)) { // ovlp
                    int min_len = MIN_OF_TWO(var1->alt_len, var2->alt_len);
                    int max_len = MAX_OF_TWO(var1->alt_len, var2->alt_len);
                    if (min_len >= max_len * 0.8 && (vntr_fuzzy_comp_seq(opt, var1->alt_seq, var1->alt_len, var2->alt_seq, var2->alt_len) == 0))
                        comp = 0;
                }
            } else if (var1->var_type == BAM_CDEL && var2->var_type == BAM_CDEL) { // two DELs
                int min_end = MIN_OF_TWO(var1->pos+var1->ref_len, var2->pos+var2->ref_len);
                int max_beg = MAX_OF_TWO(var1->pos, var2->pos);
                if (min_end >= max_beg) { // ovlp
                    int min_len = MIN_OF_TWO(var1->ref_len, var2->ref_len);
                    int max_len = MAX_OF_TWO(var1->ref_len, var2->ref_len);
                    if (min_len >= max_len * 0.8) comp = 0;
                }
            }
            // if (comp >= 0) {
            //     fprintf(stderr, "SomaticVar-i: %d %" PRId64 " %d-%c-%d %d vs %d\n", var_i,
            //             cand_vars[var_i].pos, cand_vars[var_i].ref_len, BAM_CIGAR_STR[cand_vars[var_i].var_type], 
            //             cand_vars[var_i].alt_len, cand_vars[var_i].alle_covs[0], cand_vars[var_i].alle_covs[1]);
            //     fprintf(stderr, "SomaticVar-j: %d %" PRId64 " %d-%c-%d %d vs %d\n", var_j,
            //             cand_vars[var_j].pos, cand_vars[var_j].ref_len, BAM_CIGAR_STR[cand_vars[var_j].var_type], 
            //             cand_vars[var_j].alt_len, cand_vars[var_j].alle_covs[0], cand_vars[var_j].alle_covs[1]);
            // }
            if (comp == 0) { 
                // updated_somatic_var_i[i] = 1;
                if (cand_vars[var_j].alle_covs[1] > cand_vars[var_i].alle_covs[1]) { // keep var_j
                    var_i_to_cate[var_i] = LONGCALLD_NON_VAR;
                } else { // keep var_i
                    var_i_to_cate[var_j] = LONGCALLD_NON_VAR;
                }
                cand_vars[var_i].total_cov = MAX_OF_TWO(cand_vars[var_i].total_cov, cand_vars[var_j].total_cov);
            }
            // XXX now the alle_covs[1] may not be correct, need to update based on all other reads' digars
        }
    }
    for (int i = 0; i < n_somatic_vars; ++i) {
        int var_i = somatic_var_i[i];
        // if (updated_somatic_var_i[i] == 1) {
            // cand_vars[var_i].alle_covs[0] = cand_vars[var_i].total_cov - cand_vars[var_i].alle_covs[1];
        // }
        if (var_i_to_cate[var_i] == LONGCALLD_CAND_SOMATIC_VAR) {
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "AfterMerge:\tSomaticVar: % " PRId64 " %d-%c-%d %d vs %d\n", 
                        cand_vars[var_i].pos, cand_vars[var_i].ref_len, BAM_CIGAR_STR[cand_vars[var_i].var_type], cand_vars[var_i].alt_len, 
                        cand_vars[var_i].alle_covs[0], cand_vars[var_i].alle_covs[1]);
            }
        }
    }
    // free(updated_somatic_var_i);
    free(somatic_var_i); 
}

// a) identify clean region vars
// b) additional noisy regions: low-complexity regions (hps)
// c) identify clean-region candidate somatic vars (will have additional somatic vars from noisy regions)
//     1) SNV: alt_dp >= 2
//     2) TE-SV: alt_dp >= 1, with TSD & polyA
//     3) other-SV: alt_dp >= 2
// XXX excluding all vars in the noisy regions
// update chunk_noisy_regs if any variant is overlapping with it
int classify_cand_vars(bam_chunk_t *chunk, int n_var_sites, const call_var_opt_t *opt) {
    hts_pos_t reg_beg = chunk->reg_beg, reg_end = chunk->reg_end;
    cand_var_t *cand_vars = chunk->cand_vars;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    int *var_i_to_cate = (int*)malloc(n_var_sites * sizeof(int));
    cgranges_t *var_pos_cr = cr_init(); cgranges_t *noisy_var_cr = cr_init(); // overlapping vars: DP needs to be >= min_alt_dp
    cgranges_t *low_comp_cr = chunk->low_comp_cr;
    chunk->var_i_to_cate = (int*)malloc(n_var_sites * sizeof(int));
    int min_dp = opt->min_dp, min_alt_dp = opt->min_alt_dp;
    double min_af = opt->min_af, max_af = opt->max_af, min_noisy_reg_ratio = opt->min_af;//opt->min_noisy_reg_ratio;
    int var_cate = -1;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = classify_var_cate(opt, ref_seq, ref_beg, ref_end, var, min_dp, min_alt_dp, min_af, max_af, min_noisy_reg_ratio);
        var_i_to_cate[i] = var_cate;
        if (var_cate == LONGCALLD_LOW_COV_VAR) continue; 
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
        if (var_cate == LONGCALLD_NON_VAR || var_cate == LONGCALLD_STRAND_BIAS_VAR) {
            continue;
        }
        // 1. var is in noisy regions: skip
        if (chunk->chunk_noisy_regs != NULL && chunk->chunk_noisy_regs->n_r > 0) {
            if (var->var_type == BAM_CINS) noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos, &ovlp_b, &max_b);
            else noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos+var->ref_len-1, &ovlp_b, &max_b);
            if (noisy_ovlp_n > 0) {
                var_i_to_cate[i] = LONGCALLD_NON_VAR; continue; // skip all vars in noisy regions
            }
        }
        // 2. low_cov_var
        if (var_cate == LONGCALLD_LOW_COV_VAR) {
            if (opt->out_somatic && var_is_cand_somatic(chunk, opt, var)) { // clean
                var_i_to_cate[i] = LONGCALLD_CAND_SOMATIC_VAR;
            }
            continue;
        }
        // 3. var is in low-complexity regions: add to noisy regions
        if (var_cate == LONGCALLD_REP_HET_VAR) {
            hts_pos_t var_start, var_end;
            if (var->pos >= reg_beg && var->pos <= reg_end) {
                cr_add_var_cr(opt, chunk, noisy_var_cr, low_comp_cr, var, 0); // no need to check min_noisy_reg_ratio
            }
            continue;
        }
        // 4. var is overlapping with other vars in ref, i.e., two DELs with overlapping bases: add to noisy regions (should always be low-complexity regions as well)
        // if (var->pos >= reg_beg && var->pos <= reg_end) {
            // if (var_is_ovlp_others(chunk, i)) {
                // cr_add_var_cr(opt, chunk, noisy_var_cr, low_comp_cr, var, 1);
            // }
        // }
        if (var->var_type == BAM_CINS) var_pos_ovlp_n = cr_overlap(var_pos_cr, "cr", var->pos-1, var->pos, &ovlp_b, &max_b);
        else var_pos_ovlp_n = cr_overlap(var_pos_cr, "cr", var->pos-1, var->pos+var->ref_len-1, &ovlp_b, &max_b);
        if (var_pos_ovlp_n > 1) { // multiple cand vars at the same position, need to check if # noisy reads >= min_noisy_reg_ratio, if not skip
            hts_pos_t var_start, var_end;
            if (var->pos >= reg_beg && var->pos <= reg_end) {
                cr_add_var_cr(opt, chunk, noisy_var_cr, low_comp_cr, var, 1);
            }
        }
        // 5. LONGCALLD_LOW_COV_VAR or LONGCALLD_LOW_AF_VAR: potential somatic variant
        if (var_cate == LONGCALLD_LOW_AF_VAR) {
        // if (var_cate == LONGCALLD_LOW_AF_VAR || var_cate == LONGCALLD_LOW_COV_VAR) {
            if (opt->out_somatic && var_is_cand_somatic(chunk, opt, var)) {
                var_i_to_cate[i] = LONGCALLD_CAND_SOMATIC_VAR;
                continue;
            }
            // if (var_cate == LONGCALLD_LOW_COV_VAR) continue;
        }
        // remaining low_af_var: skip
        if (var_cate == LONGCALLD_LOW_AF_VAR) {
            var_i_to_cate[i] = LONGCALLD_LOW_COV_VAR;
        }
    }
    if (noisy_var_cr->n_r > 0) {
        cr_index(noisy_var_cr);
        // add additional noisy regions from low-complexity regions
        cgranges_t *tmp_cr = cr_merge2(chunk->chunk_noisy_regs, noisy_var_cr, -1, opt->noisy_reg_merge_dis, opt->min_sv_len);
        cr_destroy(chunk->chunk_noisy_regs);
        chunk->chunk_noisy_regs = tmp_cr;
    }

    // print all noisy regions
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "After classify var: %ld total noisy regions\n", chunk->chunk_noisy_regs->n_r);
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "NoisyReg: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
    }
    if (opt->out_somatic) // candidate somatic vars from clean regions
        merge_clean_somatic_vars(opt, cand_vars, var_i_to_cate, n_var_sites);
    // add flank_len to each noisy region
    post_process_noisy_regs(chunk, opt, var_i_to_cate);
    // skip all vars that fully within noisy regions, re-call them in noisy regions
    // this does not include those that are on the boundaries,
    // XXX after post_process_noisy_regs, there should be no vars on the boundaries
    int cand_var_i = 0;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = var_i_to_cate[i];
        if (var_cate == LONGCALLD_LOW_COV_VAR || var_cate == LONGCALLD_NON_VAR || var_cate == LONGCALLD_STRAND_BIAS_VAR) continue;
        if (chunk->chunk_noisy_regs != NULL && chunk->chunk_noisy_regs->n_r > 0) {
            // int noisy_ovlp_n = cr_overlap(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos+var->ref_len, &ovlp_b, &max_b);
            if (cr_is_contained(chunk->chunk_noisy_regs, "cr", var->pos-1, var->pos+var->ref_len, &ovlp_b, &max_b) > 0) {
                var_i_to_cate[i] = LONGCALLD_NON_VAR;
                continue; // skip all vars in noisy regions
            }
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
        if (cand_vars[i].tsd_len > 0) free(cand_vars[i].tsd_seq);
    }
    // print all variants
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "After classify var: %d candidate vars\n", cand_var_i);
        for (int i = 0; i < cand_var_i; ++i) {
            cand_var_t *var = cand_vars+i;
            fprintf(stderr, "CandVarCate-%c: %s:%" PRId64 " %d-%c-%d %d\t", LONGCALLD_VAR_CATE_TYPE(chunk->var_i_to_cate[i]), chunk->tname, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len, var->total_cov);
            fprintf(stderr, "Low-Depth: %d\t", var->low_qual_cov);
            for (int k = 0; k < var->alt_len; ++k) fprintf(stderr, "%c", "ACGTN"[var->alt_seq[k]]);
            fprintf(stderr, ": %d\n", var->alle_covs[1]);
        }
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

void collect_digars_from_bam(bam_chunk_t *chunk, const struct call_var_pl_t *pl) {
    chunk->chunk_noisy_regs = cr_init();
    call_var_opt_t *opt = pl->opt;
    for (int i = 0; i < chunk->n_reads; ++i) {
        bam1_t *read = chunk->reads[i];
        if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "%d: qname: %s, flag: %d, pos: %" PRId64 ", end: %" PRId64 "\n", i, bam_get_qname(read), read->core.flag, read->core.pos+1, bam_endpos(read));
        if (chunk->is_skipped[i]) continue;
        int ret;
        if (has_equal_X_in_bam_cigar(read)) {
            ret = collect_digar_from_eqx_cigar(chunk, i, opt, chunk->digars+i);
        } else if (has_cs_in_bam(read)) {
            ret = collect_digar_from_cs_tag(chunk, i, opt, chunk->digars+i);
        } else if (has_MD_in_bam(read)) {
            ret = collect_digar_from_MD_tag(chunk, i, opt, chunk->digars+i);
        } else { // no =/X in cigar and no cs/MD tag, compare bases with ref_seq
            ret = collect_digar_from_ref_seq(chunk, i, opt, chunk->digars+i);
        }
        if (ret < 0) chunk->is_skipped[i] = BAM_RECORD_WRONG_MAP;
    }
    // print chunk->qual_counts
    int n_all_quals = 256;
    int valid_quals[256], n_valid_quals = 0; // XXX MAX_QUAL = 255
    int64_t n_total_counts = 0;
    for (int i = 0; i < n_all_quals; ++i) n_total_counts += chunk->qual_counts[i];
    for (int i = 0; i < n_all_quals; ++i) {
        if (chunk->qual_counts[i] <= 0) continue; // skip 0 counts
        if (chunk->qual_counts[i] >= 0.001 * n_total_counts) valid_quals[n_valid_quals++] = i;
        // if (LONGCALLD_VERBOSE >= 0) {
            // fprintf(stderr, "qual: %d, counts: %d, ratio: %.3f%%\n", i, chunk->qual_counts[i], (float)chunk->qual_counts[i] / n_total_counts * 100);
        // }
    }
    // collect 1st quartile/median/3rd quartile
    chunk->min_qual = valid_quals[0];
    chunk->first_quar_qual = valid_quals[n_valid_quals/4];
    chunk->median_qual = valid_quals[n_valid_quals/2];
    chunk->third_quar_qual = valid_quals[n_valid_quals*3/4];
    chunk->max_qual = valid_quals[n_valid_quals-1];
    // fprintf(stderr, "n_valid: %d, min: %d, 1st quartile: %d, median: %d, 3rd quartile: %d, max: %d\n", n_valid_quals, chunk->min_qual, chunk->first_quar_qual, chunk->median_qual, chunk->third_quar_qual, chunk->max_qual);
    if (LONGCALLD_VERBOSE < 2) {
        for (int i = 0; i < chunk->m_reads; ++i) bam_destroy1(chunk->reads[i]);
        free(chunk->reads);
    }
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

int comp_ovlp_var_site(const call_var_opt_t *opt, var_site_t *var1, var_site_t *var2, int *is_ovlp) {
    *is_ovlp = ovlp_var_site(var1, var2);
    return exact_comp_var_site(opt, var1, var2);
}

// order of variants with same pos: order by type and ref_len
int merge_var_sites(const call_var_opt_t *opt, int n_total_var_sites, var_site_t **var_sites, int tid, hts_pos_t reg_beg, hts_pos_t reg_end,
                    int n_digar, digar1_t *digars) {
    int new_total_var_sites = 0;
    var_site_t *new_var_sites = (var_site_t*)malloc((n_total_var_sites + n_digar) * sizeof(var_site_t));
    int i, j;
    for (i = j = 0; i < n_total_var_sites && j < n_digar; ) {
        hts_pos_t digar_pos = digars[j].pos;
        if (reg_beg != -1 && digar_pos < reg_beg) { j++; continue; }
        if (reg_end != -1 && digar_pos > reg_end) { break; }
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
        // if (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL) { // keep noisy-regions variants
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            // allow in-exact match for large INS???
            // int ret = comp_var_site(opt, (*var_sites)+i, &digar_var_site);
            int ret = exact_comp_var_site_ins(opt, (*var_sites)+i, &digar_var_site);
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
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            new_var_sites[new_total_var_sites++] = digar_var_site;
        }
    }
    free(*var_sites); *var_sites = new_var_sites;
    return new_total_var_sites;
}

// collect all candidate variant sites from digars, excluding low-quality/noisy-region ones
// including germline hom/het and somatic ones
int collect_all_cand_var_sites(const call_var_opt_t *opt, bam_chunk_t *chunk, var_site_t **var_sites) {
    int n_total_var_sites = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
       n_total_var_sites = merge_var_sites(opt, n_total_var_sites, var_sites, 
                                           chunk->tid, chunk->reg_beg, chunk->reg_end, //->beg, chunk->end, 
                                           chunk->digars[i].n_digar, chunk->digars[i].digars);
    }
    // fprintf(stderr, "total_cand_var: %d\n", n_total_var_sites);
    // print cand_vars
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Total candidate variant sites: %d\n", n_total_var_sites);
        for (int i = 0; i < n_total_var_sites; ++i) {
            fprintf(stderr, "CandVarSite: %s:%" PRId64 " %d-%c-%d\n", chunk->tname, (*var_sites)[i].pos, (*var_sites)[i].ref_len,
                            BAM_CIGAR_STR[(*var_sites)[i].var_type], (*var_sites)[i].alt_len);
        }
    }
    return n_total_var_sites;
}

int exact_comp_cand_var(const call_var_opt_t *opt, cand_var_t *var1, cand_var_t *var2) {
    var_site_t var_site1 = make_var_site_from_cand_var(var1);
    var_site_t var_site2 = make_var_site_from_cand_var(var2);
    return exact_comp_var_site(opt, &var_site1, &var_site2);
}

int merge_var_profile(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_new_vars, cand_var_t *new_vars, int *new_var_cate, read_var_profile_t *new_p) {
    if (n_new_vars <= 0) return 0;
    cand_var_t *old_vars = chunk->cand_vars;
    read_var_profile_t *old_p = chunk->read_var_profile;
    cand_var_t *merged_vars = (cand_var_t *)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(cand_var_t)); 
    read_var_profile_t *merged_p = init_read_var_profile(chunk->n_reads, chunk->n_cand_vars + n_new_vars);
    
    int *merged_var_i_to_cate = (int*)malloc((chunk->n_cand_vars + n_new_vars) * sizeof(int));
    int old_var_i = 0, new_var_i = 0, merged_var_i = 0;
    for (; old_var_i < chunk->n_cand_vars && new_var_i < n_new_vars; ) {
        int ret = exact_comp_cand_var(opt, old_vars+old_var_i, new_vars+new_var_i);
        if (ret < 0) { // add old_var to merged_vars, updated read_profile
            for (int i = 0; i < chunk->n_reads; ++i) {
                if (chunk->is_skipped[i]) continue;
                read_var_profile_t *old_p1 = old_p + i;
                read_var_profile_t *merged_p1 = merged_p + i;
                if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
        } else if (ret > 0) { // add new_var to merged vars, update read_profile
            for (int i = 0; i < chunk->n_reads; ++i) {
                if (chunk->is_skipped[i]) continue;
                read_var_profile_t *new_p1 = new_p + i;
                read_var_profile_t *merged_p1 = merged_p + i;
                if (new_p1->start_var_idx > new_var_i || new_p1->end_var_idx < new_var_i) continue;
                update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], new_p1->alt_qi[new_var_i-new_p1->start_var_idx], merged_p1);
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
                update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
            }
            merged_var_i_to_cate[merged_var_i] = chunk->var_i_to_cate[old_var_i];
            merged_vars[merged_var_i++] = old_vars[old_var_i++];
            // free new_vars
            free_cand_vars1(new_vars+new_var_i); new_var_i++;
        }
    }
    for (; old_var_i < chunk->n_cand_vars; ++old_var_i) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            if (chunk->is_skipped[i]) continue;
            read_var_profile_t *old_p1 = old_p + i;
            read_var_profile_t *merged_p1 = merged_p + i;
            if (old_p1->start_var_idx > old_var_i || old_p1->end_var_idx < old_var_i) continue;
            update_read_var_profile_with_allele(merged_var_i, old_p1->alleles[old_var_i-old_p1->start_var_idx], old_p1->alt_qi[old_var_i-old_p1->start_var_idx], merged_p1);
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
            update_read_var_profile_with_allele(merged_var_i, new_p1->alleles[new_var_i-new_p1->start_var_idx], new_p1->alt_qi[new_var_i-new_p1->start_var_idx], merged_p1);
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
    free(chunk->cand_vars); free(new_vars); free(new_var_cate); cr_destroy(chunk->read_var_cr);
    chunk->read_var_profile = merged_p; chunk->cand_vars = merged_vars; chunk->n_cand_vars = merged_var_i; chunk->var_i_to_cate = merged_var_i_to_cate; chunk->read_var_cr = merged_read_var_cr;

    if (LONGCALLD_VERBOSE >= 2) {
        for (int read_id = 0; read_id < chunk->n_reads; ++read_id) {
            read_var_profile_t *p1 = merged_p + read_id;
            if (chunk->is_skipped[read_id]) continue;
            fprintf(stderr, "MergedProfile: %s start_var_i: %d, end_var_i: %d\n", bam_get_qname(chunk->reads[read_id]), p1->start_var_idx, p1->end_var_idx);
            for (int k = 0; k <= p1->end_var_idx-p1->start_var_idx; ++k) {
                fprintf(stderr, "P\tVar: (%d) %" PRId64 "", k, merged_vars[k+p1->start_var_idx].pos);
                fprintf(stderr, " %d-%c-%d, allele: %d\n", merged_vars[k+p1->start_var_idx].ref_len, BAM_CIGAR_STR[merged_vars[k+p1->start_var_idx].var_type], merged_vars[k+p1->start_var_idx].alt_len, p1->alleles[k]);
            }
        }
    }
    return new_var_i; // n_vars
}

read_var_profile_t *collect_read_var_profile(const call_var_opt_t *opt, bam_chunk_t *chunk) {
    int n_cand_vars = chunk->n_cand_vars;
    cand_var_t *cand_vars = chunk->cand_vars;
    int *var_i_to_cate = chunk->var_i_to_cate;
    read_var_profile_t *p = init_read_var_profile(chunk->n_reads, n_cand_vars);
    cgranges_t *read_var_cr = cr_init();
    // 3rd pass: collect read-wise SNP profiles
    if (opt->out_somatic) {
        // init allele_cov for somatic vars, as they may be incorrect in the previous pass
        for (int i = 0; i < n_cand_vars; ++i) {
            if (var_i_to_cate[i] == LONGCALLD_CAND_SOMATIC_VAR) {
                cand_vars[i].alle_covs[0] = 0; // ref allele
                cand_vars[i].alle_covs[1] = 0; // alt allele
                cand_vars[i].total_cov = 0;
            }
        }
    }
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        // update read vs all vars, include germline and somatic vars
        // for somatic vars: also update cand_vars.alle_covs/total_cov
        start_var_i = update_read_vs_all_var_profile_from_digar(opt, chunk, chunk->digars+i, n_cand_vars, cand_vars, var_i_to_cate, start_var_i, p+i);
        if (p[i].start_var_idx < 0 || p[i].end_var_idx < 0) continue;
        cr_add(read_var_cr, "cr", p[i].start_var_idx, p[i].end_var_idx+1, i);
        if (LONGCALLD_VERBOSE >= 2) {
            if (p[i].start_var_idx >= 0) {
                bam1_t *read = chunk->reads[i];
                fprintf(stderr, "Read: %s, start_var_i: %d, end_var_i: %d\n", bam_get_qname(read), p[i].start_var_idx, p[i].end_var_idx);
                for (int j = 0; j <= p[i].end_var_idx-p[i].start_var_idx; ++j) {
                    fprintf(stderr, "P\tVar: (%d) %" PRId64 " %c", j, cand_vars[j+p[i].start_var_idx].pos, LONGCALLD_VAR_CATE_TYPE(var_i_to_cate[j+p[i].start_var_idx]));
                    int alt_qi = p[i].alt_qi[j];
                    int alt_qual = -1;
                    if (alt_qi != -1) alt_qual = chunk->digars[i].qual[alt_qi];
                    fprintf(stderr, " %d-%c-%d, allele: %d, alt_pos: %d, alt_qual: %d\n", cand_vars[j+p[i].start_var_idx].ref_len, BAM_CIGAR_STR[cand_vars[j+p[i].start_var_idx].var_type], cand_vars[j+p[i].start_var_idx].alt_len, p[i].alleles[j], alt_qi, alt_qual);
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

// func: cand_var_t -> var_t
// cand_vars: candidate variants, up to 1 alt_allele
// collect variants based on hap_to_cons_alle
// 1) filter out low-quality variants, e.g., P(var|hap,phasing)
// 2) merge variants with the overlapping pos & ref_len (e.g., ACGT -> A, ACGTCGT) XXX
// 3) extend phase blocks if possible, e.g., if ≥ 1 read supports the longer phase block
int make_variants(const call_var_opt_t *opt, bam_chunk_t *chunk, var_t **_var) {
    int n_cand_vars = chunk->n_cand_vars;
    if (n_cand_vars <= 0) return 0;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg;
    hts_pos_t active_reg_beg = chunk->reg_beg, active_reg_end = chunk->reg_end;
    cand_var_t *cand_vars = chunk->cand_vars;
    (*_var) = (var_t*)malloc(sizeof(var_t));
    var_t *var = *_var;
    var->n = 0; var->m = n_cand_vars;
    int hap1_idx=1, hap2_idx=2, hom_idx=0;
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    var->vars = (var1_t*)malloc(n_cand_vars * sizeof(var1_t));
    int i = 0, is_hom, hom_alt_is_set, hom_alle, hap1_alle, hap2_alle;
    int target_var_cate = LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR;
    if (opt->out_somatic == 1) target_var_cate |= LONGCALLD_CAND_SOMATIC_VAR;
    for (int cand_i = 0; cand_i < n_cand_vars; ++cand_i) {
        if ((chunk->var_i_to_cate[cand_i] & target_var_cate) == 0) continue;
        if (cand_vars[cand_i].var_type == BAM_CDEL || cand_vars[cand_i].var_type == BAM_CINS) {
            var->vars[i].pos = cand_vars[cand_i].pos-1;
            var->vars[i].ref_len = cand_vars[cand_i].ref_len+1;
        } else {
            var->vars[i].pos = cand_vars[cand_i].pos;
            var->vars[i].ref_len = cand_vars[cand_i].ref_len;
        }
        if (var->vars[i].pos < active_reg_beg || var->vars[i].pos > active_reg_end) continue;
        hom_alle = cand_vars[cand_i].hap_to_cons_alle[hom_idx];
        hap1_alle = cand_vars[cand_i].hap_to_cons_alle[hap1_idx];
        hap2_alle = cand_vars[cand_i].hap_to_cons_alle[hap2_idx];
        is_hom = 0; hom_alt_is_set = 0;
        if (hap1_alle == -1 && hap2_alle == -1) { // homozygous
            is_hom = 1;
            hap1_alle = hap2_alle = hom_alle;
            if (LONGCALLD_VERBOSE >= 2) _err_warning("Warning: both hap1 and hap2 are -1, set to hom_alle: %d var: %s:%" PRIi64 "\n", hom_alle, chunk->tname, cand_vars[cand_i].pos);
        } else if (hap1_alle == hap2_alle) is_hom = 1;
        if (hap1_alle == -1) hap1_alle = LONGCALLD_REF_ALLELE;
        if (hap2_alle == -1) hap2_alle = LONGCALLD_REF_ALLELE;

        var->vars[i].type = cand_vars[cand_i].var_type;
        if (chunk->var_i_to_cate[cand_i] == LONGCALLD_CAND_SOMATIC_VAR) var->vars[i].is_somatic = 1;
        else var->vars[i].is_somatic = 0;
        var->vars[i].tsd_len = cand_vars[cand_i].tsd_len;

        // retrotransposon info
        if (var->vars[i].tsd_len > 0) {
            var->vars[i].tsd_seq = (uint8_t*)malloc(var->vars[i].tsd_len * sizeof(uint8_t));
            for (int j = 0; j < var->vars[i].tsd_len; ++j) {
                var->vars[i].tsd_seq[j] = cand_vars[cand_i].tsd_seq[j];
            }
            var->vars[i].tsd_len = cand_vars[cand_i].tsd_len;
            var->vars[i].polya_len = cand_vars[cand_i].polya_len;
            var->vars[i].tsd_pos1 = cand_vars[cand_i].tsd_pos1;
            var->vars[i].tsd_pos2 = cand_vars[cand_i].tsd_pos2;
        }
        if (cand_vars[cand_i].te_seq_i >= 0) {
            var->vars[i].te_seq_i = cand_vars[cand_i].te_seq_i;
            var->vars[i].te_is_rev = cand_vars[cand_i].te_is_rev;
        } else var->vars[i].te_seq_i = -1;
        var->vars[i].PS = cand_vars[cand_i].phase_set;
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
        i++;
    }
    free(ovlp_b);
    var->n = i;
    return(var->m);
}

static void update_chunk_read_hap_phase_set1(bam_chunk_t *chunk) {
    // read haps
    if (chunk->flip_hap && chunk->flip_cur_PS != -1) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            if (chunk->haps[i] == 0) continue;
            if (chunk->phase_sets[i] == chunk->flip_cur_PS) {
                chunk->haps[i] = 3 - chunk->haps[i];
            }
        }
    }
    // read phase set
    hts_pos_t flip_pre_chunk_PS = chunk->flip_pre_PS, flip_cur_chunk_PS = chunk->flip_cur_PS;
    if (flip_pre_chunk_PS != -1 && flip_cur_chunk_PS != INT64_MAX) {
        for (int i = 0; i < chunk->n_reads; ++i) {
            if (chunk->phase_sets[i] == -1) continue;
            if (chunk->phase_sets[i] == flip_cur_chunk_PS) chunk->phase_sets[i] = flip_pre_chunk_PS;
        }
    }
}

static void update_chunk_var_hap_phase_set1(bam_chunk_t *chunk) {
    // var haps
    if (chunk->flip_hap && chunk->flip_cur_PS != -1) {
        for (int i = 0; i < chunk->n_cand_vars; ++i) {
            // if (chunk->cand_vars[i].hap_to_cons_alle[0] == 0) continue;
            if (chunk->cand_vars[i].phase_set == chunk->flip_cur_PS) {
                int tmp = chunk->cand_vars[i].hap_to_cons_alle[1];
                chunk->cand_vars[i].hap_to_cons_alle[1] = chunk->cand_vars[i].hap_to_cons_alle[2];
                chunk->cand_vars[i].hap_to_cons_alle[2] = tmp;
            }
        }
    }
    // var phase set
    hts_pos_t flip_pre_chunk_PS = chunk->flip_pre_PS, flip_cur_chunk_PS = chunk->flip_cur_PS;
    if (flip_pre_chunk_PS != -1 && flip_cur_chunk_PS != INT64_MAX) {
        for (int i = 0; i < chunk->n_cand_vars; ++i) {
            if (chunk->cand_vars[i].phase_set == -1) continue;
            if (chunk->cand_vars[i].phase_set == flip_cur_chunk_PS) chunk->cand_vars[i].phase_set = flip_pre_chunk_PS;
        }
    }
}

// XXX use overlapping reads to extend phase blocks
// 1. check if haplotype of variants in bam_chunk is inconsistent with the previous bam_chunk
// 2. extend phase blocks if possible (e.g., if ≥ 1 read supports the longer phase block)
void flip_variant_hap(call_var_opt_t *opt, bam_chunk_t *pre_chunk, bam_chunk_t *cur_chunk) {
    if (cur_chunk->tid != pre_chunk->tid) return;
    int n_cur_ovlp_reads = cur_chunk->n_up_ovlp_reads; 
    int n_pre_ovlp_reads = pre_chunk->n_down_ovlp_reads;
    // assert(n_cur_ovlp_reads == n_pre_ovlp_reads);
    if (n_cur_ovlp_reads != n_pre_ovlp_reads) {
        _err_error_exit("n_pre_ovlp_reads: %d (%s:%" PRIi64 "-%" PRIi64 "), n_cur_ovlp_reads: %d (%s:%" PRIi64 "-%" PRIi64 ")\n", n_pre_ovlp_reads, pre_chunk->tname, pre_chunk->reg_beg, pre_chunk->reg_end,
                        n_cur_ovlp_reads, cur_chunk->tname, cur_chunk->reg_beg, cur_chunk->reg_end);
    }
    if (n_cur_ovlp_reads <= 0) return;
    if (pre_chunk->n_cand_vars <= 0 || cur_chunk->n_cand_vars <= 0) return;
    // 1) find overlapping reads that ovlp with both prev and cur variants
    int *ovlp_read_i = cur_chunk->up_ovlp_read_i; // read_i in the previous bam_chunk
        int flip_hap_score = 0; hts_pos_t max_pre_read_PS = -1, min_cur_read_PS = INT64_MAX;
    for (int i = 0; i < n_cur_ovlp_reads; ++i) {
        int cur_read_i = cur_chunk->up_ovlp_read_i[i];
        int pre_read_i = pre_chunk->down_ovlp_read_i[i];
        // if (pre_read_i >= pre_chunk->n_reads) continue;
        if (LONGCALLD_VERBOSE >= 2) {
            if (strcmp(bam_get_qname(pre_chunk->reads[pre_read_i]), bam_get_qname(cur_chunk->reads[cur_read_i])) != 0) {
                _err_error_exit("OvlpRead not match: %d: cur_read: %s, pre_read: %s, %s:%" PRIi64 "-%" PRIi64 "\n", i, bam_get_qname(cur_chunk->reads[cur_read_i]), pre_chunk->reads[pre_read_i], cur_chunk->tname, cur_chunk->reg_beg, cur_chunk->reg_end);
            }
            fprintf(stderr, "read: %s pre_i: %d (%d)\n", bam_get_qname(pre_chunk->reads[pre_read_i]), pre_read_i, pre_chunk->n_reads);
            fprintf(stderr, "read: %s cur_i: %d (%d)\n", bam_get_qname(cur_chunk->reads[cur_read_i]), cur_read_i, cur_chunk->n_reads);
            fprintf(stderr, "pre_hap: %d\n", pre_chunk->haps[pre_read_i]);
            fprintf(stderr, "cur_hap: %d\n", cur_chunk->haps[cur_read_i]);
            fprintf(stderr, "pre_PS: %" PRIi64 "\n", pre_chunk->phase_sets[pre_read_i]);
            fprintf(stderr, "cur_PS: %" PRIi64 "\n", cur_chunk->phase_sets[cur_read_i]);
        }
        if (pre_chunk->is_skipped[pre_read_i] || pre_chunk->haps[pre_read_i] == 0 ||
            cur_chunk->is_skipped[cur_read_i] || cur_chunk->haps[cur_read_i] == 0) continue;
        int pre_read_hap = pre_chunk->haps[pre_read_i];
        hts_pos_t pre_read_PS = pre_chunk->phase_sets[pre_read_i];
        int cur_read_hap = cur_chunk->haps[cur_read_i];
        hts_pos_t cur_read_PS = cur_chunk->phase_sets[cur_read_i];
        // check if pre_read_hap is consistent with cur_chunk's vars/reads
        if (pre_read_hap == cur_read_hap) flip_hap_score -= 1;
        else flip_hap_score += 1;
        if (max_pre_read_PS < pre_read_PS) max_pre_read_PS = pre_read_PS;
        if (min_cur_read_PS > cur_read_PS) min_cur_read_PS = cur_read_PS;
    }
    if (flip_hap_score == 0) return; // no extension
    cur_chunk->flip_pre_PS = max_pre_read_PS;
    cur_chunk->flip_cur_PS = min_cur_read_PS;
    // fprintf(stderr, "pre_chunk_PS: %" PRId64 ", cur_chunk_PS: %" PRId64 "\n", max_pre_read_PS, min_cur_read_PS);
    int flip = flip_hap_score > 0 ? 1 : 0; // flip : keep
    cur_chunk->flip_hap = flip;
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Region: %s:%" PRIi64 "-%" PRIi64 ", flip_hap: %d (%d) pre_PS: %" PRIi64 ", cur_PS: %" PRIi64 "\n", cur_chunk->tname, cur_chunk->reg_beg, cur_chunk->reg_end, cur_chunk->flip_hap, flip_hap_score, max_pre_read_PS, min_cur_read_PS);
    update_chunk_var_hap_phase_set1(cur_chunk);
    if (opt->out_bam != NULL) update_chunk_read_hap_phase_set1(cur_chunk);
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


uint8_t collect_non_gap_base(uint8_t *ref_seq, int ref_pos) {
    while (ref_seq[ref_pos] == 5) ref_pos--;
    if (ref_pos < 0) return 4;
    return ref_seq[ref_pos];
}

int var_is_homopolymer_indel(bam_chunk_t *chunk, hts_pos_t ref_pos, int var_type, int ref_len, int alt_len, uint8_t *alt_seq) {
    if (var_type == BAM_CDIFF) return 0;
    if (var_type == BAM_CINS) { // INS
        // INS is HP base
        uint8_t ins_base0 = alt_seq[0];
        for (int i = 1; i < alt_len; ++i) {
            if (alt_seq[i] != ins_base0) return 0;
        }
        // region is HP
        for (int i = 0; i < 5; ++i) {
            if (chunk->ref_seq[ref_pos-chunk->ref_beg+i] != ins_base0) return 0;
        }
        return 1;
    } else { // DEL
        uint8_t ref_base0 = chunk->ref_seq[ref_pos-chunk->ref_beg];
        for (int i = 1; i < ref_len; ++i) {
            if (chunk->ref_seq[ref_pos-chunk->ref_beg+i] != ref_base0) return 0;
        }
        // region is HP
        for (int i = 0; i < 5; ++i) {
            if (chunk->ref_seq[ref_pos-chunk->ref_beg+i] != ref_base0) return 0;
        }
        return 1;
    }
}

void make_cand_vars0(cand_var_t *var, int tid, hts_pos_t ref_pos, int var_type, int n_alle, int ref_len, uint8_t ref_base, int alt_len, uint8_t *alt_seq,
                     int is_homopolymer_indel, int tsd_len, uint8_t *tsd_seq, hts_pos_t tsd_pos1, hts_pos_t tsd_pos2, int tsd_polya_len, int te_seq_i, int te_is_rev) {
    memset(var, 0, sizeof(cand_var_t));
    var->tid = tid; var->pos = ref_pos;
    var->var_type = var_type;
    var->n_uniq_alles = n_alle;
    var->alle_covs = (int*)calloc(n_alle, sizeof(int));
    var->ref_len = ref_len;
    var->is_homopolymer_indel = is_homopolymer_indel;
    if (var_type == BAM_CDIFF) var->ref_base = ref_base;
    if (alt_len > 0 && alt_seq != NULL) {
        var->alt_len = alt_len; 
        var->alt_seq = (uint8_t*)malloc(alt_len*sizeof(uint8_t));
        for (int i = 0; i < alt_len; ++i) var->alt_seq[i] = alt_seq[i];
    } else {
        var->alt_len = 0;
        var->alt_seq = NULL;
    }
    if (tsd_len > 0) {
        var->tsd_len = tsd_len; var->tsd_seq = tsd_seq; var->polya_len = tsd_polya_len;
        var->tsd_pos1 = tsd_pos1; var->tsd_pos2 = tsd_pos2;
    } else {
        var->tsd_len = 0; var->tsd_seq = NULL;
        var->tsd_pos1 = -1; var->tsd_pos2 = -1;
    }
    var->checked_tsd = 1;
    if (te_seq_i >= 0) {
        var->te_seq_i = te_seq_i; var->te_is_rev = te_is_rev;
    } else {
        var->te_seq_i = -1; var->te_is_rev = 0;
    }
}

// include vars: var.pos in within [reg_beg, reg_end]
// position with both gaps are filtered in advance
// for INS/DEL >= 50: look for TSD&polyA
// no_end_var: for somatic mode, not call variant for INDEL on both ends (partial clipping)
int make_cand_vars_from_baln0(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg,
                             uint8_t *ref_msa_seq, uint8_t *cons_msa_seq, int msa_len, cand_var_t **cand_vars, int **var_cate, int no_end_var) {
    *cand_vars = (cand_var_t *)malloc((msa_len + 1) * sizeof(cand_var_t));
    *var_cate = (int *)malloc((msa_len + 1) * sizeof(int));
    int n_vars = 0; int min_sv_len = opt->min_sv_len; // collect RetroTrans info
    hts_pos_t ref_pos = noisy_reg_beg; // int query_pos = 0
    int tid = chunk->tid;
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
            make_cand_vars0((*cand_vars)+n_vars, tid, ref_pos, BAM_CDIFF, 2, 1, ref_msa_seq[i], 1, cons_msa_seq+i, 0, 0, NULL, 0, 0, 0, -1, 0);
            n_vars++;
            i += 1; ref_pos += 1;
        } else if (ref_msa_seq[i] == 5) { // INS: 0->I
            int gap_len = 1;
            while (i+gap_len < msa_len && ref_msa_seq[i+gap_len] == 5 && cons_msa_seq[i+gap_len] != 5) gap_len++;
            // if no_end_var is set, both ends should be EQUAL bases
            if (no_end_var && (i-1 < 0 || i+gap_len >= msa_len || 
                               ref_msa_seq[i-1] == 5 || ref_msa_seq[i+gap_len] == 5 ||
                               cons_msa_seq[i-1] == 5 || cons_msa_seq[i+gap_len] == 5)) {
                i += gap_len;
                continue;
            }
            int is_homopolymer_indel = 0, tsd_len = 0, polya_len = 0; uint8_t *tsd_seq=NULL; hts_pos_t tsd_pos1=-1, tsd_pos2=-1; int te_seq_i=-1, te_is_rev=-1;
            if (gap_len >= min_sv_len) tsd_len = collect_te_info_from_cons(opt, chunk, ref_pos, i, BAM_CINS, gap_len, cons_msa_seq, &tsd_seq, &tsd_pos1, &tsd_pos2, &polya_len, &te_seq_i, &te_is_rev);
            else is_homopolymer_indel = var_is_homopolymer_indel(chunk, ref_pos, BAM_CINS, 0, gap_len, cons_msa_seq+i);
            make_cand_vars0((*cand_vars) + n_vars, tid, ref_pos, BAM_CINS, 2, 0, 0, gap_len, cons_msa_seq+i, is_homopolymer_indel, tsd_len, tsd_seq, tsd_pos1, tsd_pos2, polya_len, te_seq_i, te_is_rev);
            n_vars++;
            i += gap_len;
        } else if (cons_msa_seq[i] == 5) { // DEL: D->0
            int gap_len = 1;
            while (i+gap_len < msa_len && ref_msa_seq[i+gap_len] != 5 && cons_msa_seq[i+gap_len] == 5) gap_len++;
            // if no_end_var is set, both ends should be EQUAL bases
            if (no_end_var && (i-1 < 0 || i+gap_len >= msa_len || 
                               ref_msa_seq[i-1] == 5 || ref_msa_seq[i+gap_len] == 5 ||
                               cons_msa_seq[i-1] == 5 || cons_msa_seq[i+gap_len] == 5)) {
                i += gap_len; ref_pos += gap_len;
                continue;
            }
            int is_homopolymer_indel = 0, tsd_len = 0, polya_len = 0; uint8_t *tsd_seq=NULL; hts_pos_t tsd_pos1=-1, tsd_pos2=-1; int te_seq_i=-1, te_is_rev=-1;
            if (gap_len >= min_sv_len) tsd_len = collect_te_info_from_cons(opt, chunk, ref_pos, i, BAM_CDEL, gap_len, cons_msa_seq, &tsd_seq, &tsd_pos1, &tsd_pos2, &polya_len, &te_seq_i, &te_is_rev);
            else is_homopolymer_indel = var_is_homopolymer_indel(chunk, ref_pos, BAM_CDEL, gap_len, 0, NULL);
            make_cand_vars0((*cand_vars)+n_vars, tid, ref_pos, BAM_CDEL, 2, gap_len, 0, 0, NULL, is_homopolymer_indel, tsd_len, tsd_seq, tsd_pos1, tsd_pos2, polya_len, te_seq_i, te_is_rev);
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

int make_cand_vars_from_msa(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg,
                            uint8_t *ref_msa_seq, uint8_t *cons_msa_seq, int msa_len, cand_var_t **cand_vars, int **var_cate, int no_end_var) {
    uint8_t *_ref_msa_seq = (uint8_t *)malloc(msa_len * sizeof(uint8_t));
    uint8_t *_cons_msa_seq = (uint8_t *)malloc(msa_len * sizeof(uint8_t));
    int _msa_len = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (ref_msa_seq[i] != 5 || cons_msa_seq[i] != 5) {
            _ref_msa_seq[_msa_len] = ref_msa_seq[i];
            _cons_msa_seq[_msa_len++] = cons_msa_seq[i];
        }
    }
    int n_vars = make_cand_vars_from_baln0(opt, chunk, noisy_reg_beg, _ref_msa_seq, _cons_msa_seq, _msa_len, cand_vars, var_cate, no_end_var);
    free(_ref_msa_seq); free(_cons_msa_seq);
    return n_vars;
}

// var_site: {pos, var_type, ref_len, alt_len}
// order of variants with same pos: order by type, ref_len, alt_len
// only allow exact match
int exact_comp_var_site(const call_var_opt_t *opt, var_site_t *var1, var_site_t *var2) {
    if (var1->pos < var2->pos) return -1;
    if (var1->pos > var2->pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    if (var1->alt_len < var2->alt_len) return -1;
    if (var1->alt_len > var2->alt_len) return 1;

    if (var1->var_type == BAM_CDIFF || var1->var_type == BAM_CINS)
        return memcmp(var1->alt_seq, var2->alt_seq, var1->alt_len);
    return 0;
}

// allow fuzzy match for large INSs: consider same if length difference <20% of max.
// for other vars: only allow exact match
int exact_comp_var_site_ins(const call_var_opt_t *opt, var_site_t *var1, var_site_t *var2) {
    if (var1->pos < var2->pos) return -1;
    if (var1->pos > var2->pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    // for X: require same alt_seq
    // for D: require same ref_len
    // for I: small: require same alt_seq, large: require at most 10% difference in alt_len
    if (var1->var_type == BAM_CDIFF) {
        if (var1->alt_len < var2->alt_len) return -1;
        if (var1->alt_len > var2->alt_len) return 1;
        return memcmp(var1->alt_seq, var2->alt_seq, var1->alt_len); // for small INSs: only exact match is considered as a match, no fuzzy match
    } else if (var1->var_type == BAM_CINS) {
        if (var1->alt_len < opt->min_sv_len) {
            if (var1->alt_len < var2->alt_len) return -1;
            if (var1->alt_len > var2->alt_len) return 1;
            return memcmp(var1->alt_seq, var2->alt_seq, var1->alt_len); // for small INSs: only exact match is considered as a match, no fuzzy match
        } else { // for large INSs, similar length is considered as a match
            int min_len = var1->alt_len < var2->alt_len ? var1->alt_len : var2->alt_len;
            int max_len = var1->alt_len > var2->alt_len ? var1->alt_len : var2->alt_len;
            if (min_len >= max_len * 0.8) return 0; // similar length
            else return var1->alt_len - var2->alt_len; // order by length
        }
    }
    return 0;
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
    // int checked_base = 0;
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
            // checked_base++;
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
    int cur_pos = -1;
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
    int cur_pos = -1;
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

// for DELs
int get_full_cover_from_ref_cons_aln_str(aln_str_t *cons_aln_str, aln_str_t *ref_cons_aln_str, int beg_in_ref, int end_in_ref) {
    // colect target_pos and len using ref_cons_aln_str
    int cur_ref_pos=-1, cur_cons_pos=-1;
    int beg_in_cons = -1, end_in_cons = -1;
    int reach_end = 0;
    for (int i = 0; i < ref_cons_aln_str->aln_len; ++i) {
        if (ref_cons_aln_str->target_aln[i] != 5) cur_ref_pos++;
        if (ref_cons_aln_str->query_aln[i] != 5) cur_cons_pos++;
        if (i < ref_cons_aln_str->query_beg || i < ref_cons_aln_str->target_beg) continue;
        else if (i > ref_cons_aln_str->query_end || i > ref_cons_aln_str->target_end) break;

        if (cur_ref_pos == beg_in_ref && beg_in_cons == -1) {
            beg_in_cons = cur_cons_pos;
        }
        if (cur_ref_pos == end_in_ref) reach_end = 1;
        if (reach_end && ref_cons_aln_str->query_aln[i] != 5) {
            end_in_cons = cur_cons_pos;
            break;
        }
    }

    return is_cover_aln_str(cons_aln_str, beg_in_cons, end_in_cons-beg_in_cons+1);
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
        int full_cover = 0;
        allele_i = get_var_allele_i_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_ref_pos-delta_ref_alt, var_alt_len, cons_sim_thres, &full_cover);
        if (full_cover) {
            vars[i].total_cov++;
            if (allele_i != -1) vars[i].alle_covs[allele_i]++;
            update_read_var_profile_with_allele(i, allele_i, -1, p); // XXX read_base_pos: -1
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
void update_cand_var_profile_from_cons_aln_str21(int clu_idx, aln_str_t *cons_aln_str, aln_str_t *ref_cons_aln_str, hts_pos_t ref_pos_beg,
                                                 cand_var_t *vars, int n_vars, int *var_from_cons_idx, read_var_profile_t *p) {
    uint8_t *_cons_seq = cons_aln_str->target_aln; uint8_t *_cons_read_seq = cons_aln_str->query_aln; int _cons_aln_len = cons_aln_str->aln_len;
    float cons_sim_thres = 0.9; //, ref_sim_thres = 0.9; //, alt_sim_thres = 1.0;

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "cons:\t");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_seq[i]]); fprintf(stderr, "\nread:\t");
        for (int i = 0; i < _cons_aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[_cons_read_seq[i]]); fprintf(stderr, "\n");
    }
    int allele_i = -1, delta_ref_alt = 0;
    // collect var_ref_beg and var_ref_end for each var

    for (int i = 0; i < n_vars; ++i) {
        int var_beg_in_ref_str = vars[i].pos - ref_pos_beg, var_end_in_ref_str;
        int var_ref_len = vars[i].ref_len, var_alt_len = vars[i].alt_len;
        if (vars[i].var_type == BAM_CINS) var_end_in_ref_str = var_beg_in_ref_str;
        else var_end_in_ref_str = var_beg_in_ref_str + var_ref_len - 1;
        int full_cover = 0;
        if (var_from_cons_idx[i] & clu_idx) { // check if read contain the var (1)
            allele_i = get_var_allele_i_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_beg_in_ref_str-delta_ref_alt, var_alt_len, cons_sim_thres, &full_cover);
        } else { // var is from the other haplotype, only check if read covers the var
            if (vars[i].var_type != BAM_CDEL) full_cover = get_full_cover_from_cons_aln_str(cons_aln_str, vars[i].var_type, var_beg_in_ref_str-delta_ref_alt, var_ref_len);
            else full_cover = get_full_cover_from_ref_cons_aln_str(cons_aln_str, ref_cons_aln_str, var_beg_in_ref_str-1, var_end_in_ref_str+1);
            // fprintf(stderr, "FullCover: %d-%c %d\n", vars[i].ref_len, BAM_CIGAR_STR[vars[i].var_type], full_cover);
            allele_i = 0;
        }
        if (full_cover) {
            vars[i].total_cov++;
            if (allele_i != -1) vars[i].alle_covs[allele_i]++;
            update_read_var_profile_with_allele(i, allele_i, -1, p); // XXX read_base_pos: -1
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
int update_cand_var_profile_from_cons_aln_str2(const call_var_opt_t *opt, bam_chunk_t *chunk, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, 
                                          hts_pos_t noisy_reg_beg, cand_var_t *hap1_vars, int *hap1_var_cate, int n_hap1_vars, cand_var_t *hap2_vars, int *hap2_var_cate, int n_hap2_vars, 
                                          cand_var_t **noisy_vars, int **noisy_var_cate, read_var_profile_t **p) {
    if (n_hap1_vars + n_hap2_vars == 0) return 0;
    int n_vars = 0;
    (*noisy_vars) = (cand_var_t *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(cand_var_t));
    (*noisy_var_cate) = (int *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(int));
    int i1 = 0, i2 = 0;
    int *var_from_cons_idx = (int *)malloc((n_hap1_vars + n_hap2_vars) * sizeof(int));
    for (; i1 < n_hap1_vars && i2 < n_hap2_vars; ) {
        int ret = exact_comp_cand_var(opt, hap1_vars + i1, hap2_vars + i2);
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
            free_cand_vars1(hap2_vars + i2); i2++;
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
    for (int i = 0; i < 2; ++i) {
        aln_str_t *clu_aln_str = aln_strs[i];
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_id = clu_read_ids[i][j];
            read_var_profile_t *p1 = *p + read_id;
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%d: %s\n", read_id, bam_get_qname(chunk->reads[read_id]));
            aln_str_t *cons_aln_str = LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, j);
            aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
            update_cand_var_profile_from_cons_aln_str21(i+1, cons_aln_str, ref_cons_aln_str, noisy_reg_beg, *noisy_vars, n_vars, var_from_cons_idx, p1);
        }
    }
    free(var_from_cons_idx);
    return n_vars;
}
// XXX allow partial-aligned reads
// msa_seqs: ref + cons + clu_n_seqs
int make_vars_from_msa_cons_aln(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_reads, int *read_ids, hts_pos_t noisy_reg_beg,
                                int n_cons, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs,
                                cand_var_t **noisy_vars, int **noisy_var_cate, read_var_profile_t **p) {
    if (n_cons == 0) return 0;
    // WFA for ref vs cons1/2
    int n_vars = 0, n_hap1_vars = 0, n_hap2_vars = 0;
    cand_var_t *hap1_vars = NULL, *hap2_vars = NULL; int *hap1_var_cate=NULL, *hap2_var_cate=NULL;
    for (int i = 0; i < n_cons; ++i) {
        uint32_t *cigar_buf = NULL;
        if (n_cons == 1) {
            aln_str_t *clu_aln_str = aln_strs[0];
            aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
            uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
            n_vars = make_cand_vars_from_msa(opt, chunk, noisy_reg_beg,
                                             ref_seq_aln, cons_seq_aln, aln_len, noisy_vars, noisy_var_cate, 0);
            if (n_vars > 0) {
                *p = init_read_var_profile(chunk->n_reads, n_vars);
                update_cand_var_profile_from_cons_aln_str1(clu_n_seqs[0], clu_read_ids[0], clu_aln_str, noisy_reg_beg, *noisy_vars, n_vars, *p);
            } else {
                free(*noisy_vars); free(*noisy_var_cate);
            }
        } else {
            if (i == 0) {
                aln_str_t *clu_aln_str = aln_strs[0];
                aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
                uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
                n_hap1_vars = make_cand_vars_from_msa(opt, chunk, noisy_reg_beg,
                                                      ref_seq_aln, cons_seq_aln, aln_len, &hap1_vars, &hap1_var_cate, 0);
            } else if (i == 1) {
                aln_str_t *clu_aln_str = aln_strs[1];
                aln_str_t *ref_cons_aln_str = LONGCALLD_REF_CONS_ALN_STR(clu_aln_str);
                uint8_t *ref_seq_aln = ref_cons_aln_str->target_aln, *cons_seq_aln = ref_cons_aln_str->query_aln; int aln_len = ref_cons_aln_str->aln_len;
                n_hap2_vars = make_cand_vars_from_msa(opt, chunk, noisy_reg_beg,
                                                      ref_seq_aln, cons_seq_aln, aln_len, &hap2_vars, &hap2_var_cate, 0);
            }
        }
        if (cigar_buf != NULL) free(cigar_buf);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_vars; ++i)
            fprintf(stderr, "HOM: %" PRId64 ", %d-%c-%d\n", (*noisy_vars)[i].pos, (*noisy_vars)[i].ref_len, BAM_CIGAR_STR[(*noisy_vars)[i].var_type], (*noisy_vars)[i].alt_len);
    }
    if (n_cons == 2) {
        n_vars = update_cand_var_profile_from_cons_aln_str2(opt, chunk, clu_n_seqs, clu_read_ids, aln_strs, noisy_reg_beg,
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

// use <exact_comp_var_site_ins> for two somatic vars from noisy regions
int merge_somatic_vars2(const call_var_opt_t *opt, int n_vars1, cand_var_t *vars1, int *var_cate1,
                        int n_vars2, cand_var_t *vars2, int *var_cate2) {
    int i = 0, j = 0;
    for (; i < n_vars1 && j < n_vars2; ) {
        if (var_cate1[i] == LONGCALLD_NON_VAR || vars1[i].alle_covs[1] <= 0) {
            i++; continue;
        }
        if (var_cate2[j] == LONGCALLD_NON_VAR || vars2[j].alle_covs[1] <= 0) {
            j++; continue;
        }
        var_site_t var_site1 = make_var_site_from_cand_var(vars1+i);
        var_site_t var_site2 = make_var_site_from_cand_var(vars2+j);
        int ret = fuzzy_comp_var_site(opt, &var_site1, &var_site2);
        // int ret = exact_comp_var_site_ins(opt, &var_site1, &var_site2);
        // fprintf(stderr, "ExactInsCompCandVar ret: %d %" PRId64 " %d-%c-%d vs %" PRId64 " %d-%c-%d\n", ret,
                // vars1[i].pos, vars1[i].ref_len, BAM_CIGAR_STR[vars1[i].var_type], vars1[i].alt_len,
                // vars2[j].pos, vars2[j].ref_len, BAM_CIGAR_STR[vars2[j].var_type], vars2[j].alt_len);
        if (ret == 0) {
            vars1[i].alle_covs[1]++; vars2[j].alle_covs[1]--;
            i++; j++;
        } else if (ret < 0) i++;
        else j++;
    }
    return 0;
}

void sort_cand_vars(cand_var_t **vars, int **var_cates, int n_vars) {
    for (int i = 0; i < n_vars-1; i++) {
        for (int j = i+1; j < n_vars; j++) {
            if ((*vars)[i].pos > (*vars)[j].pos) {
                cand_var_t tmp = (*vars)[i]; (*vars)[i] = (*vars)[j]; (*vars)[j] = tmp;
                int tmp_cate = (*var_cates)[i]; (*var_cates)[i] = (*var_cates)[j]; (*var_cates)[j] = tmp_cate;
            }
        }
    }
}

int get_gap_len_from_aln_str(uint8_t *aln_str, int aln_len, int beg_pos, int ext_right) {
    if (beg_pos < 0 || beg_pos >= aln_len) return 0;
    if (ext_right) {
        for (int i = beg_pos; i < aln_len; ++i) {
            if (aln_str[i] != 5) return i - beg_pos;
        }
        return aln_len - beg_pos;
    } else {
        for (int i = beg_pos; i >= 0; --i) {
            if (aln_str[i] != 5) return beg_pos - i;
        }
        return beg_pos + 1;
    }
}

int confirm_somatic_var_based_on_aln_str(hts_pos_t noisy_reg_beg, cand_var_t *somatic_var, aln_str_t *ref_cons_aln_str, aln_str_t *cons_read_aln_str) {
    int ret = 0;
    int var_cons_pos = -1, var_ref_pos = -1, ref_cons_aln_var_pos = -1;
    for (int i = 0; i < ref_cons_aln_str->aln_len; ++i) {
        if (ref_cons_aln_str->target_aln[i] != 5) var_ref_pos++;
        if (ref_cons_aln_str->query_aln[i] != 5) var_cons_pos++;
        if (var_ref_pos + noisy_reg_beg == somatic_var->pos) {
            ref_cons_aln_var_pos = i;
            break;
        }
    }
    // assert(var_cons_pos > -1);
    int cons_pos = -1, read_pos = -1;
    for (int i = 0; i < cons_read_aln_str->aln_len; ++i) {
        if (cons_read_aln_str->target_aln[i] != 5) cons_pos++;
        if (cons_read_aln_str->query_aln[i] != 5) read_pos++;
        if (cons_pos == var_cons_pos) {
            if (somatic_var->var_type == BAM_CDIFF) {
                if (cons_pos >= 0 && read_pos >= 0 &&
                    ref_cons_aln_str->target_aln[ref_cons_aln_var_pos] == cons_read_aln_str->target_aln[i] &&
                    somatic_var->alt_seq[0] == cons_read_aln_str->query_aln[i]) {
                    ret = 1;
                }
            } else if (somatic_var->var_type == BAM_CINS) {
                int ins_len = get_gap_len_from_aln_str(cons_read_aln_str->target_aln, cons_read_aln_str->aln_len, i-1, 0);
                if (ins_len == somatic_var->alt_len) ret = 1;
            } else if (somatic_var->var_type == BAM_CDEL) {
                int del_len = get_gap_len_from_aln_str(cons_read_aln_str->query_aln, cons_read_aln_str->aln_len, i, 1);
                if (del_len == somatic_var->ref_len) ret = 1;
            }
            break;
        }
    }
    return ret;
}

// for somatic vs germline vars, we allow fuzzy matching
// as ref_read_aln_str is based on global alignment, so when merging candidate somatic vars: 
// we do not allow position-wise fuzzy matching, only allow fuzzy matching for large INS: <exact_comp_var_site_ins>
int make_somatic_vars_from_aln_str(const call_var_opt_t *opt, bam_chunk_t *chunk, int n_reads, int *read_ids, hts_pos_t noisy_reg_beg, int n_noisy_vars, cand_var_t *noisy_vars,
                                   int n_cons, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs,
                                   cand_var_t **noisy_somatic_vars, int **noisy_somatic_var_cate, read_var_profile_t **noisy_somatic_p) {
    int n_total_reads = 0, n_somatic_vars = 0, m_somatic_vars = 0;
    for (int i = 0; i < n_cons; ++i) {
        n_total_reads += clu_n_seqs[i];
        m_somatic_vars += aln_strs[i]->aln_len;
    }
    *noisy_somatic_vars = (cand_var_t *)malloc(m_somatic_vars * sizeof(cand_var_t));
    *noisy_somatic_var_cate = (int *)malloc(m_somatic_vars * sizeof(int));
    cand_var_t **read_wise_noisy_somatic_vars = (cand_var_t **)malloc(n_total_reads * sizeof(cand_var_t*));
    int *read_wise_n_noisy_somatic_vars = (int *)malloc(n_total_reads * sizeof(int));
    int **read_wise_noisy_somatic_var_cate = (int **)malloc(n_total_reads * sizeof(int*));
    int read_i = 0;
    for (int i = 0; i < n_cons; ++i) {
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            aln_str_t *ref_read_aln_str = LONGCALLD_REF_READ_ALN_STR(aln_strs[i], j);
            int read_id = clu_read_ids[i][j];
            int clu_n_seq = clu_n_seqs[i];
            // 1) collect all candidate somatic vars from msa_aln of read vs ref alignment (from POA, so could be wrong)
            int n_vars = make_cand_vars_from_msa(opt, chunk, noisy_reg_beg, ref_read_aln_str->target_aln, ref_read_aln_str->query_aln, ref_read_aln_str->aln_len,
                                                 &read_wise_noisy_somatic_vars[read_i], &read_wise_noisy_somatic_var_cate[read_i], 1);
            // fprintf(stderr, "read_id: %d, n_vars: %d\n", read_id, n_vars);
            read_wise_n_noisy_somatic_vars[read_i] = n_vars;
            for (int k = 0; k < n_vars; ++k) {
                if ((read_wise_noisy_somatic_vars[read_i][k].var_type == BAM_CINS && read_wise_noisy_somatic_vars[read_i][k].alt_len < opt->min_sv_len)
                 || (read_wise_noisy_somatic_vars[read_i][k].var_type == BAM_CDEL && read_wise_noisy_somatic_vars[read_i][k].ref_len < opt->min_sv_len)) {
                    read_wise_noisy_somatic_var_cate[read_i][k] = LONGCALLD_NON_VAR; // skip small indels
                    continue;
                }
                // check if the candidate somatic var can be confirmed by read vs cons alignment
                // int ret = confirm_somatic_var_based_on_aln_str(noisy_reg_beg, read_wise_noisy_somatic_vars[read_i] + k, ref_cons_aln_str, cons_read_aln_str);
                // fprintf(stderr, "Var: %" PRId64 ", %d-%c-%d %d\n", read_wise_noisy_somatic_vars[read_i][k].pos, read_wise_noisy_somatic_vars[read_i][k].ref_len,
                        // BAM_CIGAR_STR[read_wise_noisy_somatic_vars[read_i][k].var_type], read_wise_noisy_somatic_vars[read_i][k].alt_len, read_i);
                // if (ret == 0) { // inconsistent between readVSref and readVScons
                    // read_wise_noisy_somatic_var_cate[read_i][k] = LONGCALLD_NON_VAR; // skip
                // } else {
                read_wise_noisy_somatic_var_cate[read_i][k] = LONGCALLD_CAND_SOMATIC_VAR;
                read_wise_noisy_somatic_vars[read_i][k].alle_covs[1] = 1;
                // }
            }
            // check if somatic vars are fuzzy matching with germline vars
            for (int ii = 0; ii < n_vars; ++ii) {
                if (read_wise_noisy_somatic_var_cate[read_i][ii] == LONGCALLD_NON_VAR) continue;
                var_site_t var_site1 = make_var_site_from_cand_var(read_wise_noisy_somatic_vars[read_i]+ii);
                int var_len = var_site1.alt_len;
                if (var_site1.var_type == BAM_CDEL) var_len = var_site1.ref_len;
                int var_win = MAX_OF_TWO(500, var_len);
                for (int jj = 0; jj < n_noisy_vars; ++jj) {
                    var_site_t var_site2 = make_var_site_from_cand_var(noisy_vars+jj);
                    if (var_site2.var_type != var_site1.var_type) continue;
                    if (var_site2.pos < var_site1.pos - var_win) continue;
                    if (var_site2.pos > var_site1.pos + var_win) break; // no more overlap
                    int ret = fuzzy_comp_var_site(opt, &var_site1, &var_site2);
                    if (ret == 0)  {
                        read_wise_noisy_somatic_var_cate[read_i][ii] = LONGCALLD_NON_VAR;
                        break;
                    }
                }
            }
            read_i++;
        }
    }
    // 2) merge candidate somatic INS/DEL vars; call cons if >=3 candidate INS/DEL vars
    for (int i = 0; i < n_total_reads-1; ++i) {
        for (int j = i+1; j < n_total_reads; ++j) {
            // fprintf(stderr, "Merge: %d-%d, %d-%d\n", i, j, read_wise_n_noisy_somatic_vars[i], read_wise_n_noisy_somatic_vars[j]);
            merge_somatic_vars2(opt, read_wise_n_noisy_somatic_vars[i], read_wise_noisy_somatic_vars[i], read_wise_noisy_somatic_var_cate[i],
                                read_wise_n_noisy_somatic_vars[j], read_wise_noisy_somatic_vars[j], read_wise_noisy_somatic_var_cate[j]);
        }
    }
    // 3) filter by var_is_cand_somatic -> var_cate = LONGCALLD_CAND_SOMATIC_VAR
    for (int i = 0; i < n_total_reads; ++i) {
        for (int j = 0; j < read_wise_n_noisy_somatic_vars[i]; ++j) {
            if (read_wise_noisy_somatic_var_cate[i][j] == LONGCALLD_NON_VAR || read_wise_noisy_somatic_vars[i][j].alle_covs[1] <= 0) {
                free_cand_vars1(read_wise_noisy_somatic_vars[i] + j);
            } else {
                if (var_is_cand_somatic(chunk, opt, read_wise_noisy_somatic_vars[i]+j)) {
                    // if ((read_wise_noisy_somatic_vars[i][j].var_type == BAM_CINS || read_wise_noisy_somatic_vars[i][j].var_type == BAM_CDEL) && read_wise_noisy_somatic_vars[i][j].alle_covs[1] >= 3) {
                        // call cons XXX
                    // }
                    // fprintf(stderr, "n: %d, m: %d\n", n_somatic_vars, m_somatic_vars);
                    // fprintf(stderr, "%" PRId64 "\t%d-%c-%d\t%d\t%d\t%c\n", read_wise_noisy_somatic_vars[i][j].pos, read_wise_noisy_somatic_vars[i][j].ref_len, BAM_CIGAR_STR[read_wise_noisy_somatic_vars[i][j].var_type], read_wise_noisy_somatic_vars[i][j].alt_len,
                            // read_wise_noisy_somatic_vars[i][j].alle_covs[0], read_wise_noisy_somatic_vars[i][j].alle_covs[1], LONGCALLD_VAR_CATE_TYPE(read_wise_noisy_somatic_var_cate[i][j]));
                    (*noisy_somatic_vars)[n_somatic_vars] = read_wise_noisy_somatic_vars[i][j];
                    (*noisy_somatic_var_cate)[n_somatic_vars++] = LONGCALLD_CAND_SOMATIC_VAR;
                } else {
                    free_cand_vars1(read_wise_noisy_somatic_vars[i] + j);
                }
            }
        }
        free(read_wise_noisy_somatic_vars[i]); free(read_wise_noisy_somatic_var_cate[i]);
    }
    free(read_wise_noisy_somatic_var_cate); free(read_wise_noisy_somatic_vars); free(read_wise_n_noisy_somatic_vars);
    // 4) collect read_var_profile
    if (n_somatic_vars <= 0) {
        free(*noisy_somatic_vars); free(*noisy_somatic_var_cate);
        return 0;
    }
    sort_cand_vars(noisy_somatic_vars, noisy_somatic_var_cate, n_somatic_vars);
    for (int i = 0; i < n_somatic_vars; ++i) (*noisy_somatic_vars)[i].alle_covs[1] = 0;
    *noisy_somatic_p = init_read_var_profile(chunk->n_reads, n_somatic_vars);
    for (int i = 0; i < n_cons; ++i) {
        int start_var_i = 0;
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_id = clu_read_ids[i][j];
            start_var_i = update_read_vs_somatic_var_profile_from_digar(opt, chunk, chunk->digars+read_id, n_somatic_vars, *noisy_somatic_vars, start_var_i, (*noisy_somatic_p)+read_id);
            // get full_cover & allele_i: update total_cov and alle_covs
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "NoisyRegSomatic read_var_profile:\n");
                if ((*noisy_somatic_p)[read_id].start_var_idx >= 0) {
                    bam1_t *read = chunk->reads[read_id];
                    fprintf(stderr, "Read: %s, start_var_i: %d, end_var_i: %d\n", bam_get_qname(read), (*noisy_somatic_p)[read_id].start_var_idx, (*noisy_somatic_p)[read_id].end_var_idx);
                    for (int k = 0; k <= (*noisy_somatic_p)[read_id].end_var_idx-(*noisy_somatic_p)[read_id].start_var_idx; ++k) {
                        fprintf(stderr, "P\tVar: (%d) %" PRId64 " %c", k, (*noisy_somatic_vars)[k+(*noisy_somatic_p)[read_id].start_var_idx].pos, LONGCALLD_VAR_CATE_TYPE((*noisy_somatic_var_cate)[k+(*noisy_somatic_p)[read_id].start_var_idx]));
                        int alt_qi = (*noisy_somatic_p)[read_id].alt_qi[k];
                        int alt_qual = -1;
                        if (alt_qi != -1) alt_qual = chunk->digars[read_id].qual[alt_qi];
                        fprintf(stderr, " %d-%c-%d, allele: %d, alt_pos: %d, alt_qual: %d\n", (*noisy_somatic_vars)[k+(*noisy_somatic_p)[read_id].start_var_idx].ref_len, BAM_CIGAR_STR[(*noisy_somatic_vars)[k+(*noisy_somatic_p)[read_id].start_var_idx].var_type], (*noisy_somatic_vars)[k+(*noisy_somatic_p)[read_id].start_var_idx].alt_len, (*noisy_somatic_p)[read_id].alleles[k], alt_qi, alt_qual);
                    }
                }
            }
        }
    }
    return n_somatic_vars;
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
                    fprintf(stderr, "DIFF: old: %" PRIi64 " %d %c %d add: %" PRIi64 " %d %c %d\n", old_vars->vars[old_i].pos,
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
    int max_noisy_reg_len = opt->max_noisy_reg_len, max_noisy_reg_reads = opt->max_noisy_reg_reads;
    if (noisy_reg_end - noisy_reg_beg + 1 > max_noisy_reg_len) {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped long region: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " (>%d)\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, max_noisy_reg_len);
        free(ref_seq);
        return 0;
    }
    int *noisy_reads; int n_noisy_reads = collect_noisy_reg_reads1(chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, &noisy_reads);
    if (n_noisy_reads > max_noisy_reg_reads) {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped deep region: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        free(noisy_reads); free(ref_seq);
        return 0;
    }

    double realtime0 = realtime();

    if (LONGCALLD_VERBOSE >= 1) {
        fprintf(stderr, "NoisyReg: chunk_reg: %s:%" PRIi64 "-%" PRIi64 ", reg: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " (%d)\n", chunk->tname, chunk->reg_beg, chunk->reg_end, 
                        chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
    }
    // MSA and consensus calling
    int n_cons = 0; int *clu_n_seqs = (int*)calloc(2, sizeof(int)); int **clu_read_ids = (int**)malloc(2 * sizeof(int*));
    // for alignment of cons vs ref, read vs cons, read vs ref
    aln_str_t **aln_strs = (aln_str_t**)malloc(2 * sizeof(aln_str_t*));
    for (int i = 0; i < 2; ++i) {
        clu_read_ids[i] = NULL;
        aln_strs[i] = (aln_str_t*)malloc((1 + n_noisy_reads * 2) * sizeof(aln_str_t));
        for (int j = 0; j < 1+n_noisy_reads*2; ++j) {
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
            fprintf(stderr, "Skipped region: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        n_noisy_vars = -1;
        goto collect_noisy_vars1_end;
    }

    cand_var_t *noisy_vars = NULL; int *noisy_var_cate=NULL; read_var_profile_t *noisy_rvp = NULL;
    // 1. make variants and update cand_vars & read_var_profile
    //    if out_somatic==1: collect candidate somatic variants + read_var_profile for all reads
    // 2. collect RetroTrans info for cons: germline
    //    if out_somatic==1: collect RetroTans info for somatic variant reads
    n_noisy_vars = make_vars_from_msa_cons_aln(opt, chunk, n_noisy_reads, noisy_reads, noisy_reg_beg,
                                               n_cons, clu_n_seqs, clu_read_ids, aln_strs, &noisy_vars, &noisy_var_cate, &noisy_rvp);
    if (opt->out_somatic) {
        int n_noisy_somatic_var=0; cand_var_t *noisy_somatic_vars=NULL; int *noisy_somatic_var_cate=NULL; read_var_profile_t *noisy_somatic_rvp = NULL;
        n_noisy_somatic_var = make_somatic_vars_from_aln_str(opt, chunk, n_noisy_reads, noisy_reads, noisy_reg_beg, n_noisy_vars, noisy_vars,
                                                             n_cons, clu_n_seqs, clu_read_ids, aln_strs,
                                                             &noisy_somatic_vars, &noisy_somatic_var_cate, &noisy_somatic_rvp);
        merge_var_profile(opt, chunk, n_noisy_vars, noisy_vars, noisy_var_cate, noisy_rvp);
        merge_var_profile(opt, chunk, n_noisy_somatic_var, noisy_somatic_vars, noisy_somatic_var_cate, noisy_somatic_rvp);
    } else merge_var_profile(opt, chunk, n_noisy_vars, noisy_vars, noisy_var_cate, noisy_rvp);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "NoisyReg: chunk_reg: %s:%" PRIi64 "-%" PRIi64 ", reg: %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " (%d)\t", chunk->tname, chunk->reg_beg, chunk->reg_end, 
                        chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reads);
        fprintf(stderr, "Real time: %.3f sec.\n", realtime() - realtime0);
    }
collect_noisy_vars1_end:
    for (int i = 0; i < 2; ++i) {
        if (clu_read_ids[i] != NULL) free(clu_read_ids[i]);
        for (int j = 0; j < 1+n_noisy_reads*2; ++j) {
            if (aln_strs[i][j].target_aln != NULL) free(aln_strs[i][j].target_aln);
        }
        free(aln_strs[i]);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "DuringCollectNoisy\n");
        for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
            cand_var_t *var = chunk->cand_vars + var_i;
            hts_pos_t var_pos = var->pos;
            if (var->var_type != BAM_CDIFF) var_pos--;
            if (var_pos < chunk->reg_beg || var_pos > chunk->reg_end) continue;
            int var_cate = chunk->var_i_to_cate[var_i];
            if (var_cate == LONGCALLD_CAND_SOMATIC_VAR) {
                    fprintf(stderr, "SOMATIC: %" PRIi64 " %d-%c-%d %d\n", chunk->cand_vars[var_i].pos, 
                            chunk->cand_vars[var_i].ref_len, BAM_CIGAR_STR[chunk->cand_vars[var_i].var_type], chunk->cand_vars[var_i].alt_len, chunk->cand_vars[var_i].alle_covs[1]);
            }
        }
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

// read digars, check if there are dense seq errors, i.e. > max_n_seq_errors (mismatch, insertion, deletion) in a sliding window of size win
int read_has_dense_seq_error(int max_n_seq_errors, int win, digar_t *digars) {
    int n_digar = digars->n_digar;
    if (n_digar <= 0) return 0; // no digars, no seq errors
    // record the position of all seq errors (mismatch, insertion, deletion)
    hts_pos_t *seq_error_pos = (hts_pos_t*)malloc(n_digar * sizeof(hts_pos_t));
    int n_seq_errors = 0;
    for (int i = 0; i < n_digar; ++i) {
        if (digars->digars[i].type == BAM_CEQUAL || digars->digars[i].type == BAM_CREF_SKIP) continue;
        if (digars->digars[i].type == BAM_CSOFT_CLIP || digars->digars[i].type == BAM_CHARD_CLIP) continue; // skip clipping
        seq_error_pos[n_seq_errors++] = digars->digars[i].pos;
    }
    // check if there are dense seq errors in the sliding window
    int has_dense_error = 0;
    for (int i = max_n_seq_errors+1; i < n_seq_errors; ++i) {
        hts_pos_t cur_pos = seq_error_pos[i];
        hts_pos_t last_pos = seq_error_pos[i - max_n_seq_errors-1];
        if (cur_pos - last_pos <= win) {
            fprintf(stderr, "dense: %ld-%ld\n", last_pos, cur_pos);
            has_dense_error = 1; // found dense seq errors
            break;
        }
    }
    free(seq_error_pos);
    return has_dense_error;
}

// n_clean_non_low_comp_het_vars: (hom)
// 0. check if "nearby" n_clean_conflict_vars / n_clean_vars is too high, i.e., >= 0.3
// XXX should be moved to assign_hap.c, i.e., after digars are updated and candidate somatic vars are collected
// 1. check if cand. somatic vars are near long clippings
// 2. check if cand. somatic vars are near dense seq errors
void mark_invalid_somatic_reads(const call_var_opt_t *opt, bam_chunk_t *chunk) {
    int min_clip_len = 500;
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i]) continue;
        // 1. mark reads with long clipping ???
        // if only has long clipping, but no dense seq error, it should be good??
        // if (chunk->is_ont_palindrome[read_i] == 0) { // ignore clipping of ONT palindrome reads
        //     if ((chunk->digars[read_i].digars[0].type == BAM_CSOFT_CLIP || chunk->digars[read_i].digars[0].type == BAM_CHARD_CLIP) && chunk->digars[read_i].digars[0].len > min_clip_len) {
        //         if (LONGCALLD_VERBOSE >= 2) {
        //             fprintf(stderr, "read %s is skipped for somatic variant calling due to long clipping\n", bam_get_qname(chunk->reads[read_i]));
        //         }
        //         chunk->is_skipped_for_somatic[read_i] = 1;
        //         continue;
        //     }
        //     int digar_i = chunk->digars[read_i].n_digar - 1;
        //     if ((chunk->digars[read_i].digars[digar_i].type == BAM_CSOFT_CLIP || chunk->digars[read_i].digars[digar_i].type == BAM_CHARD_CLIP) && chunk->digars[read_i].digars[digar_i].len > min_clip_len) {
        //         chunk->is_skipped_for_somatic[read_i] = 1;
        //         if (LONGCALLD_VERBOSE >= 2) {
        //             fprintf(stderr, "read %s is skipped for somatic variant calling due to long clipping\n", bam_get_qname(chunk->reads[read_i]));
        //         }
        //         continue;
        //     }
        // }
        // 2. mark reads with dense seq errors, here, digars are already updated
        // int win, max_n_seq_errors;
        // if (opt->is_ont) { // ONT
        //     win=50; max_n_seq_errors=20;
        // } else {
        //     win=50; max_n_seq_errors=10;
        // }
        // if (read_has_dense_seq_error(max_n_seq_errors, win, chunk->digars+read_i)) {
        //     if (LONGCALLD_VERBOSE >= 2) {
        //         fprintf(stderr, "read %s is skipped for somatic variant calling due to dense seq errors\n", bam_get_qname(chunk->reads[read_i]));
        //     }
        //     chunk->is_skipped_for_somatic[read_i] = 1;
        //     continue;
        // }
        // 3. mark reads with conflicting phasing
        if (chunk->n_clean_conflict_snps[read_i] >= 2) {
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "read %s is skipped for somatic variant calling due to conflicting phasing (%d)\n", bam_get_qname(chunk->reads[read_i]), chunk->n_clean_conflict_snps[read_i]);
            }
            chunk->is_skipped_for_somatic[read_i] = 1;
            continue;
        }
    }
}
// collect somatic variants based on phased reads
// clean-reg somatic vars + noisy-reg somatic vars
// 1) require >= 2 (min_somatic_alt_dp) reads: <250 (min_somatic_te_len) bp INDEL or SNVs
// 2) require >= 1 (min_somatic_te_dp) reads: >=250 (min_somatic_te_len) bp INDEL, potential TEs
// LONGCALLD_CAND_GERMLINE_VAR:
// 3) previous candidate germline variants with refCall
// func: update var->phase_set, var->hap_to_cons_alle[1/2]
void collect_somatic_var(bam_chunk_t *chunk, const call_var_opt_t *opt) {
    // loop through all reads, mark non-somatic-var reads, i.e., dense seq error; long clipping; weak phasing;
    // for read with dense seq errors or long clipping, mark it as skipped for somatic variant calling, potentially sequencing artifact or alignment artifact
    mark_invalid_somatic_reads(opt, chunk);
    // loop through all candidate variants
    // 1) SOMATIC SNV; 2) SOMATIC SVs; 3) refCall germline variants
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        cand_var_t *var = chunk->cand_vars + var_i;
        hts_pos_t var_pos = var->pos;
        if (var->var_type != BAM_CDIFF) var_pos--;
        if (var_pos < chunk->reg_beg || var_pos > chunk->reg_end) continue;
        int var_cate = chunk->var_i_to_cate[var_i];

        if (LONGCALLD_VERBOSE >= 2) {
            if (var_cate == LONGCALLD_CAND_SOMATIC_VAR) {
                    fprintf(stderr, "SOMATIC: %" PRIi64 " %d-%c-%d %d\n", chunk->cand_vars[var_i].pos, 
                            chunk->cand_vars[var_i].ref_len, BAM_CIGAR_STR[chunk->cand_vars[var_i].var_type], chunk->cand_vars[var_i].alt_len, chunk->cand_vars[var_i].alle_covs[1]);
            } else if (var_cate & LONGCALLD_CAND_GERMLINE_VAR_CATE) {
                if (var->hap_to_cons_alle[1] == 0 && var->hap_to_cons_alle[2] == 0 && var_is_cand_somatic(chunk, opt, var)) { // refCall
                    fprintf(stderr, "RefCall: %" PRIi64 " %d-%c-%d %d\n", chunk->cand_vars[var_i].pos, 
                            chunk->cand_vars[var_i].ref_len, BAM_CIGAR_STR[chunk->cand_vars[var_i].var_type], chunk->cand_vars[var_i].alt_len, chunk->cand_vars[var_i].alle_covs[1]);
                }
            }
        }
        if (var_cate & LONGCALLD_CAND_GERMLINE_VAR_CATE) {
            if (var->hap_to_alle_profile != NULL) {
                if (var->hap_to_cons_alle[1] == 0 && var->hap_to_cons_alle[2] == 0 // refCall
                    && var_is_cand_somatic(chunk, opt, var))
                chunk->var_i_to_cate[var_i] = LONGCALLD_CAND_SOMATIC_VAR;
            } else {
                fprintf(stderr, "Warning: No hap_to_alle_profile for somatic var: %s:%" PRIi64 " %d-%c-%d %d\n", 
                        chunk->tname, chunk->cand_vars[var_i].pos, chunk->cand_vars[var_i].ref_len, BAM_CIGAR_STR[chunk->cand_vars[var_i].var_type], 
                        chunk->cand_vars[var_i].alt_len, chunk->cand_vars[var_i].alle_covs[1]);
            }
        }
    }
    assign_somatic_hap_based_on_phased_reads(opt, chunk, LONGCALLD_CAND_SOMATIC_VAR);
}

void collect_var_main(const call_var_pl_t *pl, bam_chunk_t *chunk) {
    call_var_opt_t *opt = pl->opt;
    // first round: easy-to-call SNPs (+indels)
    // 1.1 collect X/I/D sites from BAM
    collect_digars_from_bam(chunk, pl);

    // 1.2. merge all var sites from all reads, including low-depth ones, but not including nosiy-region ones
    var_site_t *var_sites = NULL; int n_var_sites;
    n_var_sites = collect_all_cand_var_sites(opt, chunk, &var_sites);
    // after collect_all_cand_var_sites: candidate somatic vars are not merged
    if (n_var_sites > 0) {
        // 1.3. collect all candidate variants, not including noisy-region ones
        // collect reference and alternative alleles for all var sites
        // all cand vars, including true/false germline/somatic variants
        // XXX for noisy/repeat regions, we need to carefully pick the candidate variants, based on supporting counts & re-alignments
        collect_cand_vars(opt, chunk, n_var_sites, var_sites);
    } free(var_sites);

    // 2.1. pre-process noisy regions (>5 XIDs in 100 win), not including low-complexity regions
    // current noisy regs are directly from BAM (large indels, VNTRs)
    pre_process_noisy_regs(chunk, opt);

    // 2.2. identify clean region vars
    // 2.3. additional noisy regions: low-complexity regions (hps)
    // (2.4) identify candidate somatic vars 
    // NOTE: after _classify_, n_cand_vars & cand_vars: only include clean region vars (may contain somatic vars)
    if (n_var_sites > 0)
        classify_cand_vars(chunk, n_var_sites, opt);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "AllNoisyRegions: %s:%" PRIi64 "-%" PRIi64 " (%ld)\n", chunk->tname, chunk->reg_beg, chunk->reg_end, chunk->chunk_noisy_regs ? chunk->chunk_noisy_regs->n_r : 0);
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "NoisyRegion: %s:%d-%d (%d)\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
    }
    if (chunk->n_cand_vars == 0 && (chunk->chunk_noisy_regs == NULL || chunk->chunk_noisy_regs->n_r == 0))
        return; // no variant to be called
    if (chunk->n_cand_vars > 0) { // XXX what if only have homozygous vars
        // 3.1. collect clean region read-wise var profiles
        chunk->read_var_profile = collect_read_var_profile(opt, chunk);

        // process candidate variants in the following order
        // within each round, use previously obtained phasing/haplotype information as intialization
        //   1. LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR
        //   2. LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR
        //   3. LONGCALLD_CAND_SOMATIC_VAR
        // 3.2. co-phasing and variant calling using clean region SNPs + indels
        assign_hap_based_on_germline_het_vars_kmeans(opt, chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR);
    }
    // 4. iteratively call variants in noisy regions and variant/read phasing
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
                        // assign_hap_based_on_het_vars(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR, pl->opt);
                        // assign_hap_based_on_germline_het_vars_kmeans(chunk, LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HET_VAR | LONGCALLD_NOISY_CAND_HOM_VAR, pl->opt);
                    } else {
                        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "No var: %s:%d-%d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, noisy_reg_i), cr_end(chunk->chunk_noisy_regs, noisy_reg_i));
                    }
                } // unable to resolve the region
            }
            if (new_var) {
                // re-assign haplotype information based on the all(old+newly) called variants
                assign_hap_based_on_germline_het_vars_kmeans(opt, chunk, LONGCALLD_CAND_GERMLINE_VAR_CATE);
            }
            if (new_region_is_done == 0) break;
        }
        free(sorted_noisy_regs); free(noisy_reg_is_done);
    }
    // 5. call somatic variants based on phased reads
    if (opt->out_somatic == 1) collect_somatic_var(chunk, opt);
}

// stitch ii-1 and ii
void stitch_var_main(call_var_step_t *step, bam_chunk_t *chunk) {
    call_var_pl_t *pl = step->pl; call_var_opt_t *opt = pl->opt;
    // extend phase set between two adjacent bam chunks
    for (int ii = 1; ii < step->n_chunks; ++ii) {
        flip_variant_hap(opt, chunk+ii-1, chunk+ii);
    }
}

void make_var_main(call_var_step_t *step, bam_chunk_t *chunk, var_t *var, long ii) {
    // generate variant-related information, e.g., GT, DP, etc.
    // merge variants based on ref_pos if needed
    call_var_pl_t *pl = step->pl; call_var_opt_t *opt = pl->opt;
    var_t *_vars;
    if (make_variants(opt, chunk, &_vars) > 0) {
        if (_vars->n > 0) merge_vars(var, _vars);
        free(_vars->vars); free(_vars);
    }
}