#ifndef LONGCALLD_ASSIGN_HAP_H
#define LONGCALLD_ASSIGN_HAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct call_var_opt_t;

int assign_hap_based_on_het_vars(struct bam_chunk_t *bam_chunk, int var_cate_i, const struct call_var_opt_t *opt);

#ifdef __cplusplus
}
#endif

#endif

