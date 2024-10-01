#ifndef LONGCALLD_PHASE_CALL_H
#define LONGCALLD_PHASE_CALL_H

#include "bam_utils.h"
#include "call_var.h"

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct call_var_opt_t;

int call_var_from_phased_reads(bam_chunk_t *bam_chunk, int var_cate_i, const call_var_opt_t *opt);
int count_based_call_var_from_phased_reads(struct bam_chunk_t *bam_chunk, int var_cate_i, const struct call_var_opt_t *opt);
int cigar_based_call_var_from_phased_reads(struct bam_chunk_t *bam_chunk, int var_cate_i, const struct call_var_opt_t *opt);
int align_based_call_var_from_phased_reads(struct bam_chunk_t *bam_chunk, int var_cate_i, kstring_t *ref_seq, const struct call_var_opt_t *opt);

#ifdef __cplusplus
}
#endif

#endif

