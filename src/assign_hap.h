#ifndef LONGCALLD_ASSIGN_HAP_H
#define LONGCALLD_ASSIGN_HAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct call_var_opt_t;
struct cand_somatic_var_aux_info_t;

int assign_hap_based_on_germline_het_vars_kmeans(const struct call_var_opt_t *opt, struct bam_chunk_t *chunk, int target_var_cate);
int assign_somatic_hap_based_on_phased_reads(const struct call_var_opt_t *opt, struct bam_chunk_t *chunk, int target_var_cate);
void free_somatic_var_aux_info(struct cand_somatic_var_aux_info_t *aux_info);

#ifdef __cplusplus
}
#endif

#endif

