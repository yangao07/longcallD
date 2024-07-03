#ifndef LONGCALLD_CALL_VAR_H
#define LONGCALLD_CALL_VAR_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"
#include "htslib/vcf.h"
#include "bam_utils.h"

#define CALL_VAR_PL_THREAD_N 3
#define CALL_VAR_THREAD_N 8


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

// // site wise summary data
// typedef struct {
//     // static information
//     hts_pos_t pos, phase_set;
//     int n_depth, n_uniq_bases; // ref+alt, used in variant calling & haplotype assignment
//     int n_low_depth; // including bases/regions with low quality, only count depth, not base
//     int8_t ref_base; // for VCF output
//     uint8_t *bases; int *base_covs, *base_to_i; // 'ACGTN.D', and corresponding coverage

//     // dynamic information, update during haplotype assignment
//     uint8_t *base_to_hap; // SNP-wise (base_i_to_hap): 0123456(ACGTN.D) -> 1:H1/2:H2/0:not set yet
//     int **hap_to_base_profile; // read-wise: 1:H1/2:H2 -> base_i -> read count
//     int *hap_to_cons_base; // HAP-wise (hap_to_cons_base_i): 1:H1/2:H2 -> 0123456(ACGTN.D)

//     // XXX
//     uint8_t is_low_qual; // read evidence is low quality to support this site as SNP
//                          // will be skipped or need extra operations to confirm
//     uint8_t is_skipped;  // skipped in VCF output
// } cand_snp_t;

typedef struct {
    uint8_t type; // 0: SNP, 1: insertion, 2: deletion, etc.
    hts_pos_t pos, phase_set; int len; // PS
    uint8_t *ref_bases; int ref_len;
    uint8_t **alt_bases; int *alt_len;
    int n_allele; // including ref allele
    int total_depth, depths[2]; uint8_t genotype[2]; // DP/AD/GT
    int qual; // phred-scaled, GQ
} var_t;

typedef struct call_var_opt_t {
    // input files
    char *ref_fa_fn;
    char *in_bam_fn; char *sample_name;
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
    // int verbose;
} call_var_opt_t;

// shared data for all threads
typedef struct {
    // input files
    ref_seq_t *ref_seq;
    samFile *bam; bam_hdr_t *header; hts_idx_t *idx; hts_itr_t *iter; int use_iter;
    // parameters, output files
    const call_var_opt_t *opt;
    // TODO: variant calling result
    // var_t *vars; int n_vars, m_vars; // SNPs, indels, complex variants (e.g. MNPs, SVs, duplications, inversions, translocations, etc.)
    // m-threads
    int max_reads_per_chunk, max_reg_len_per_chunk;
    int n_threads;
} call_var_pl_t;

// separately/parallelly processed
typedef struct {
    const call_var_pl_t *pl;
    int n_chunks, max_chunks; // 64
    bam_chunk_t *chunks;
} call_var_step_t;

int call_var_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_call_var_H