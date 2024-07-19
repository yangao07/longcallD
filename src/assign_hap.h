#ifndef LONGCALLD_ASSIGN_HAP_H
#define LONGCALLD_ASSIGN_HAP_H

#ifdef __cplusplus
extern "C" {
#endif

struct cand_snp_t;
struct read_snp_profile_t;
struct bam_chunk_t;

int assign_hap(struct read_snp_profile_t *p, int n_cand_snps, struct cand_snp_t *cand_snps, struct bam_chunk_t *bam_chunk);

#ifdef __cplusplus
}
#endif

#endif

