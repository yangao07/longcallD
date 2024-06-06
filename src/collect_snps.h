#ifndef LONGCALLD_COLLECT_SNPS_H
#define LONGCALLD_COLLECT_SNPS_H

#include "call_var.h"
#include "utils.h"

#ifdef __cplusplus
extern "C" {
#endif

// typedef struct {
//     // snps (mismatch bases)
//     int n_snps, m_snps;
//     int *snp_pos; // size: n_snps
//     int *n_snp_bases; uint8_t **snp_bases; // size: n_snps * n_snp_bases (A, C, G, T, N)
//     uint8_t *used_for_phasing_snps; // size: n_snps, 0: not used, 1: used
//     // reads
//     int n_reads; // init as 512 (100kb win, 10kb len. 50X cov -> 500 reads)
//     int *read_to_start_snp, *read_to_end_snp; // size: n_reads
//     // snps X reads
//     uint8_t **read_to_snp_map; // size: n_snps * n_reads
// } cand_snps_t;

int collect_snps_main(const call_var_pl_t *pl, bam_chunk_t *bam_chunk);

#ifdef __cplusplus
}
#endif

#endif // COLLECT_SNPS_H