#ifndef LONGCALLD_ASSIGN_HAP_H
#define LONGCALLD_ASSIGN_HAP_H

#include "bam_utils.h"

#ifdef __cplusplus
extern "C" {
#endif
// // input: tid, tname, start, end, reads, n_reads
// // output: HPs
// typedef struct {
//     int tid; char *tname; uint64_t start, end;
//     uint8_t bam_has_eqx_cigar, bam_has_md_tag;
//     int n_reads, m_reads;
//     bam1_t **reads; // size: m_reads
//     int *HPs; // size: m_reads
// } bam_chunk_t;

int assign_hap(read_snp_profile_t *p, int n_cand_snps, cand_snp_t *cand_snps, bam_chunk_t *bam_chunk);

#ifdef __cplusplus
}
#endif

#endif

