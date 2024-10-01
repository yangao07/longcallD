#ifndef LONGCALLD_VCF_UTILS_H
#define LONGCALLD_VCF_UTILS_H

#ifdef __cplusplus
extern "C" {
#endif

struct cand_var_t;
struct bam_chunk_t;
struct var_t;

int write_vcf_header(bam_hdr_t *hdr, FILE *out_vcf, char *sample_name);
int write_var_to_vcf(struct var_t *vars, FILE *out_vcf, char *chrom);

#ifdef __cplusplus
}
#endif

#endif