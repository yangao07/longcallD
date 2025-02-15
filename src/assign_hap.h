#ifndef LONGCALLD_ASSIGN_HAP_H
#define LONGCALLD_ASSIGN_HAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct call_var_opt_t;

int assign_hap_based_on_het_vars_kmeans(struct bam_chunk_t *chunk, int target_var_cate, struct call_var_opt_t *opt);

#ifdef __cplusplus
}
#endif

#endif

