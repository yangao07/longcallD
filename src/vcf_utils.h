#ifndef LONGCALLD_VCF_UTILS_H
#define LONGCALLD_VCF_UTILS_H

#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

struct cand_var_t;
struct bam_chunk_t;
struct var_t;

int write_vcf_header(bam_hdr_t *hdr, struct call_var_opt_t *opt);
int write_var_to_vcf(struct var_t *vars, const struct call_var_opt_t *opt, char *chrom);

#ifdef __cplusplus
}
#endif

#endif