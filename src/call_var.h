#ifndef LONGCALLD_CALL_VAR_H
#define LONGCALLD_CALL_VAR_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"
#include "cgranges.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#define CALL_VAR_PL_THREAD_N 3
#define CALL_VAR_THREAD_N 8

#define LONGCALLD_MIN_CAND_MQ 5 // low-qual
#define LONGCALLD_MIN_CAND_BQ 0 // low-qual
#define LONGCALLD_MIN_CAND_DP 5 // total DP < 5: skipped
#define LONGCALLD_MIN_ALT_DP 2 // max alt depth < 2: skipped
#define LONGCALLD_MIN_SOMATIC_AF 0.05 // AF < 0.05: filtered out
#define LONGCALLD_MIN_CAND_AF 0.25 // AF < 0.25: not germline het.
#define LONGCALLD_MAX_CAND_AF 0.75 // AF > 0.75: not germline het.
#define LONGCALLD_MAX_LOW_QUAL_FRAC 0.25 // LOW_FRAC > 0.25: low-qual
#define LONGCALLD_DEF_PLOID 2 // diploid
#define LONGCALLD_REF_ALLELE 0
#define LONGCALLD_ALT_ALLELE1 1
#define LONGCALLD_ALT_ALLELE2 2
#define LONGCALLD_OTHER_ALT_ALLELE 3

#define LONGCALLD_DENSE_REG_MAX_XGAPS 10 // 5?  // dense X/gap region: more than n X/gap bases in a 100-bp window
#define LONGCALLD_DENSE_REG_SLIDE_WIN 100 //
#define LONGCALLD_DENSE_FLANK_WIN 25 // 100/(3+1)
#define LONGCALLD_NOISY_END_CLIP 100 // >= n bp clipping on both ends
#define LONGCALLD_NOISY_END_CLIP_WIN 100 // n bp flanking end-clipping region will be considered as low-quality region
#define LONGCALLD_INDEL_FLANK_WIN 0 // n bp around indel will be considered as low-quality region

#ifdef __cplusplus
extern "C" {
#endif

// for each vartiant:
// #site-level:
// * QUAL = -10 * log(1-p) (int), site-wise, p = P(data|ref)
// * FILTER: PASS, LowQual, RefCall, NoCall
// * INFO: DP/AD/AF/END
// #sample-level:
// * FORMAT: GT/GQ/DP/AD/PL/PS 
//     GP = P(G|D)
//     GL = P(D|G)
//   * PL (int) = -10 * log10(P(D|G)), Phred-scaled likelihoods of the possible genotypes.
//        "Normalized": the PL of the most likely genotype is set to 0, and the rest are scaled relative to this.
//   * GQ (int) = PL_sec - PL_lowest (int)
//   * GT = argmin(PL)
typedef struct {
    uint8_t type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    hts_pos_t pos, PS; // phase set
    uint8_t *ref_bases; int ref_len; // extra care needed for multi-allelic sites, e.g., GGTGT -> GGT,G
    int n_alt_allele; // alt allele: 1 or 2
    uint8_t **alt_bases; int *alt_len;
    int DP, AD[2]; uint8_t GT[2]; // DP/AD/GT
    int QUAL, GQ, PL[6]; // phred-scaled, QUAL/FILTER/GQ/PL
} var1_t;

typedef struct var_t {
    int n, m;
    var1_t *vars;
} var_t;

typedef struct call_var_opt_t {
    // input
    char *ref_fa_fn; char *rep_bed_fn; char *in_bam_fn; char *sample_name;
    char *region_list; uint8_t region_is_file; // for -R/--region/--region-file option
    // filters for variant calling
    int max_ploid, min_mq, min_bq, min_dp, min_alt_dp; 
    double min_somatic_af, min_af, max_af, max_low_qual_frac;
    int dens_reg_max_xgaps, dens_reg_slide_win, dens_reg_flank_win;
    int indel_flank_win;
    int end_clip_reg, end_clip_reg_flank_win;
    // general
    // int max_ploidy;
    int pl_threads, n_threads;
    // output
    htsFile *out_bam; // phased bam
    FILE *out_vcf; // phased vcf
} call_var_opt_t;

// shared data for all threads
typedef struct call_var_pl_t {
    // input files
    ref_seq_t *ref_seq; cgranges_t *rep_regs;
    samFile *bam; bam_hdr_t *header; hts_tpool *p;
    hts_idx_t *idx; hts_itr_t *iter; int use_iter;
    // parameters, output files
    const call_var_opt_t *opt;
    // m-threads
    int max_reads_per_chunk, max_reg_len_per_chunk;
    int n_threads;
} call_var_pl_t;

struct bam_chunk_t;
// separately/parallelly processed
typedef struct call_var_step_t {
    const call_var_pl_t *pl;
    int n_chunks, max_chunks; // 64
    struct bam_chunk_t *chunks;
    var_t *vars; // size: n_chunks; SNP/indel, SV, etc.
} call_var_step_t;

int call_var_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_call_var_H