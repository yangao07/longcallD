#ifndef LONGCALLD_VCF_UTILS_H
#define LONGCALLD_VCF_UTILS_H

#include "bam_utils.h"
#include <stdio.h>

#ifdef __cplusplus
extern "C" {
#endif

int write_vcf_header(bam_hdr_t *hdr, FILE *out_vcf, char *sample_name);
int write_snp_to_vcf(cand_snp_t *cand_snps, int n_cand_snps, FILE *out_vcf, char *chrom);

#ifdef __cplusplus
}
#endif

#endif