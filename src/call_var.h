#ifndef LONGCALLD_CALL_VAR_H
#define LONGCALLD_CALL_VAR_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"

#define CALL_VAR_PL_THREAD_N 3
#define CALL_VAR_THREAD_N 4

#define LONGCALLD_MIN_CAND_SNP_MQ 0 // aln with MQ < 5 will be skipped
#define LONGCALLD_MIN_CAND_SNP_BQ 0 // base with BQ < 10 will be skipped, only when qual is available
#define LONGCALLD_MIN_CAND_SNP_DP 5 // site with high-qual DP < 5 will be skipped, high-qual: MQ >= 5, BQ >= 10
#define LONGCALLD_MIN_CAND_SNP_AF 0.25 // base with AF < 0.25 will be skipped
#define LONGCALLD_MAX_CAND_SNP_AF 0.75 // base with AF > 0.75 will be skipped
#define LONGCALLD_MAX_LOW_QUAL_FRAC 0.25 // base with LOW_FRAC > 0.25 will be skipped
#define LONGCALLD_DEF_PLOID 2 // diploid

// dense X/gap region: more than 5 X/gap bases in a 100-bp window
#define LONGCALLD_DENSE_REG_MAX_SITES 3   // 2
#define LONGCALLD_DENSE_REG_SLIDE_WIN 100 // 10
#define LONGCALLD_INDEL_FLANK_WIN_SIZE 3 // 3 bp around indel will be considered as low-quality region

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
    uint8_t type; // 0: SNP, 1: insertion, 2: deletion, etc.
    hts_pos_t pos, PS; int len; // PS
    uint8_t *ref_bases; int ref_len;
    uint8_t **alt_bases; int *alt_len;
    int n_allele; // including ref allele
    int DP, AD[2]; uint8_t GT[2]; // DP/AD/GT
    int QUAL, GQ, PL[6]; // phred-scaled, QUAL/FILTER/GQ/PL
} var1_t;

typedef struct var_t {
    int n, m;
    var1_t *vars;
} var_t;

typedef struct call_var_opt_t {
    // input files
    char *ref_fa_fn; char *in_bam_fn; char *sample_name;
    // int max_ploidy;
    char *region_list; uint8_t region_is_file; // for -R/--region/--region-file option

    // filters for cand. SNPs used in phasing
    int max_ploid, min_mq, min_bq, min_dp; 
    double min_af, max_af, max_low_qual_frac;
    int dens_reg_max_sites, dens_reg_slide_win, indel_flank_win_size;

    int pl_threads, n_threads;
    // output
    htsFile *out_bam; // phased bam
    FILE *out_vcf; // phased vcf

    // general
} call_var_opt_t;

// shared data for all threads
typedef struct call_var_pl_t {
    // input files
    ref_seq_t *ref_seq;
    samFile *bam; bam_hdr_t *header; hts_idx_t *idx; hts_itr_t *iter; int use_iter;
    // parameters, output files
    const call_var_opt_t *opt;
    // m-threads
    int max_reads_per_chunk, max_reg_len_per_chunk;
    int n_threads;
} call_var_pl_t;

struct bam_chunk_t;
// separately/parallelly processed
typedef struct {
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