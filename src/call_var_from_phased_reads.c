#include "call_var_from_phased_reads.h"
#include "collect_var.h"
// #include "bam_utils.h"
// #include "cgranges.h"
//#include "htslib/kstring.h"
#include "htslib/kfunc.h"
// #include "htslib/sam.h"
// #include "call_var.h"
#include "align.h"

static inline double fisher_exact_test_p(int n11, int n12, int n21, int n22) {
    double left, right, two;
    return kt_fisher_exact(n11, n12, n21, n22, &left, &right, &two);
}

static inline double genotype_prob_from_phased_counts(int n10, int n11, int n20, int n21, double *p_ref_call, double *p_hom_call, double *p_h1_call, double *p_h2_call) {
    *p_ref_call = fisher_exact_test_p(n10, n11, n20, n21);
    *p_hom_call = fisher_exact_test_p(n11, n10, n21, n20);
    *p_h1_call = fisher_exact_test_p(n11, n10, n20, n21);
    *p_h2_call = fisher_exact_test_p(n21, n20, n10, n11);
}

// TODO: output p for each REP_HET_VAR/DENSE_REEEG_VAR
// see if true variants have lower p than false ones
// fisher-exact test for two haplotypes:
//   for each alt: (at most 2 alts)
//       calculate p using [H1.(n_alt, n_non_alt), H2.(n_alt, n_non_alt)]
// XXX no update of reads' HAP here, only call var
int call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, const call_var_opt_t *opt) {
    read_var_profile_t *p = bam_chunk->read_var_profile;
    cgranges_t *read_var_cr = bam_chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    int *var_idx = bam_chunk->var_cate_idx[var_cate_i];
    for (int i = 0; i < bam_chunk->var_cate_counts[var_cate_i]; ++i) {
        int var_i = var_idx[i];
        cand_var_t *var = bam_chunk->cand_vars + var_i;
        hts_pos_t pos = var->pos; int ref_len = var->ref_len;
        ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
        // for each alt allele
        for (int alt_i = 1; alt_i < var->n_uniq_alles; ++alt_i) {
            // count reads supporting each allele
            int n10 = 0, n11 = 0, n20 = 0, n21 = 0; // H1_non_alt, H1_alt, H2_non_alt, H2_alt
            double p_ref_call = 0, p_hom_call = 0, p_h1_call = 0, p_h2_call = 0;
            for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
                int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
                int read_var_idx = var_i - p[read_i].start_var_idx;
                if (bam_chunk->is_skipped[read_i] != 1 && p[read_i].var_is_used[read_var_idx] == 1) {
                    int hap = bam_chunk->haps[read_i];
                    int alle = p[read_i].alleles[read_var_idx]; // 0:ref, 1/2:alt
                    if (hap == 1) {
                        if (alle == alt_i) n11++;
                        else n10++;
                    } else if (hap == 2) {
                        if (alle == alt_i) n21++;
                        else n20++;
                    }
                }
            }
            genotype_prob_from_phased_counts(n10, n11, n20, n21, &p_ref_call, &p_hom_call, &p_h1_call, &p_h2_call);
            fprintf(stderr, "VarCate-%c: %s:%ld\t%c\t%d,%d,%d,%d\t", LONGCALLD_VAR_CATE_STR[var_cate_i], bam_chunk->tname, pos, BAM_CIGAR_STR[var->var_type], n10, n11, n20, n21);
            fprintf(stderr, "p_ref_call: %f, p_hom_call: %f, p_h1_call: %f, p_h2_call: %f\n", p_ref_call, p_hom_call, p_h1_call, p_h2_call);
            // double p = fisher_exact_test_p(n10, n11, n20, n21);
            // fprintf(stderr, "VarCate-%c: %s:%lld\t%c\t%f\n", LONGCALLD_VAR_CATE_STR[var_cate_i], bam_chunk->tname, pos, BAM_CIGAR_STR[var->var_type], p);
        }
    }
    free(ovlp_b); 
    return 0;
}

// for smaller number of alt. alleles within each phased cluster, we can use count-based method
// input: phased bam
// process: 
//  1) phased reads: count reads supporting ref/alt allele, pick the allele with higher count
//  2) unphased reads: collect consensus from phased ones, then compare with ref/alt cons. allele, poetential needs alignment
// output: genotype for each variant
// return: 0 if success, 1 if failed
int count_based_call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, const call_var_opt_t *opt) {
    read_var_profile_t *p = bam_chunk->read_var_profile;
    int *var_idx = bam_chunk->var_cate_idx[var_cate_i];
    for (int i = 0; i < bam_chunk->var_cate_counts[var_cate_i]; ++i) {
        int var_i = var_idx[i];
        cand_var_t *var = bam_chunk->cand_vars + var_i;
        hts_pos_t pos = var->pos; int ref_len = var->ref_len;
        // collect allele seq
        // count reads supporting each allele
    }
    return 0;
}

int digar_based_call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, const call_var_opt_t *opt) {
    return 0;
}

// only do this for hard-to-call variants, e.g., supporting count is not enough
// for easy-to-call ones, i.e., supporting count is enough, we can use count-based method
// for larger number of alt. alleles within each phased cluster, or top 2 haplotypes are too close, we need to use align-based method
// @input: phased reads, candidate variants
// @process: WFA/edlib-based, align phased reads to ref/alt allele sequences, pick the best allele with higher score
//           for each variant, collect all reads covering the variant, use phased reads first, then use unphased reads
// @output: genotype for each variant
//          update read-var-profile after alignment
int align_based_call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, kstring_t *ref_seq, const call_var_opt_t *opt) {
    read_var_profile_t *p = bam_chunk->read_var_profile;
    int *var_idx = bam_chunk->var_cate_idx[var_cate_i];
    for (int i = 0; i < bam_chunk->var_cate_counts[var_cate_i]; ++i) {
        int var_i = var_idx[i];
        cand_var_t *var = bam_chunk->cand_vars + var_i;
        hts_pos_t pos = var->pos, non_rep_beg = var->non_rep_beg, non_rep_end = var->non_rep_end;
        // collect ref/alt allele seq
        char *ref_allele_seq = ref_seq->s + non_rep_beg - 1; 
        int ref_allele_len = non_rep_end - non_rep_beg + 1;
        // for two haplotypes, collect reads seq
        // calculate S(read-seq|allele-seq)
    }
    return 0;
}

//
int msa_based_call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, kstring_t *ref_seq, const call_var_opt_t *opt) {
    read_var_profile_t *p = bam_chunk->read_var_profile;
    int *var_idx = bam_chunk->var_cate_idx[var_cate_i];
    for (int i = 0; i < bam_chunk->var_cate_counts[var_cate_i]; ++i) {
        int var_i = var_idx[i];
        cand_var_t *var = bam_chunk->cand_vars + var_i;
        hts_pos_t pos = var->pos, non_rep_beg = var->non_rep_beg, non_rep_end = var->non_rep_end;
        // collect ref/alt allele seq
        char *ref_allele_seq = ref_seq->s + non_rep_beg - 1; 
        int ref_allele_len = non_rep_end - non_rep_beg + 1;
        // for two haplotypes, collect reads seq
        // calculate S(read-seq|allele-seq)
    }
    return 0;
}