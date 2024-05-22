#ifndef LONGCALLD_PHASE_BAM_H
#define LONGCALLD_PHASE_BAM_H

#include <stdint.h>
#include <stdio.h>
#include "bam_utils.h"
#include "seq.h"
#include "htslib/vcf.h"

#define PHASE_BAM_PL_THREAD_N 3
#define PHASE_BAM_THREAD_N 8


#define LONGCALLD_MIN_CAND_SNP_MQ 0 // aln with MQ < 5 will be skipped
#define LONGCALLD_MIN_CAND_SNP_BQ 0 // base with BQ < 10 will be skipped, only when qual is available
#define LONGCALLD_MIN_CAND_SNP_DP 10 // base with high-qual DP < 5 will be skipped, high-qual: MQ >= 5, BQ >= 10
#define LONGCALLD_MIN_CAND_SNP_AF 0.25 // base with AF < 0.25 will be skipped
#define LONGCALLD_MAX_CAND_SNP_AF 0.75 // base with AF > 0.75 will be skipped
#define LONGCALLD_DEF_PLOID 2 // diploid

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    // input files
    char *ref_fa_fn;
    char *in_bam_fn;
    // int max_ploidy;
    char *region_list;
    uint8_t region_is_file;

    // filters for cand. SNPs used in phasing
    int max_ploid, min_mq, min_bq, min_dp; double min_af, max_af;

    int pl_threads, n_threads;
    // output
    htsFile *out_bam; // phased bam
    // vcfFile *out_vcf; // phased vcf
} phase_bam_opt_t;

// shared data for all threads
typedef struct {
    // papameters, output files
    const phase_bam_opt_t *pb_opt;
    // input files
    ref_seq_t *ref_seq;
    samFile *bam; bam_hdr_t *header; hts_idx_t *idx; // hts_itr_t *iter;
    // m-threads
    int max_reads_per_chunk, n_threads;
} phase_bam_pl_t;

// separately/parallelly processed
typedef struct {
    const phase_bam_pl_t *pl;
    int n_chunks, max_chunks; // 64
    bam_chunk_t *chunks;
} phase_bam_step_t;

int phase_bam_main(int argc, char *argv[]);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_PHASE_BAM_H