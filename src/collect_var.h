#ifndef LONGCALLD_COLLECT_VAR_H
#define LONGCALLD_COLLECT_VAR_H

#include "htslib/sam.h"

#ifdef __cplusplus
extern "C" {
#endif

// site wise summary data
typedef struct cand_snp_t {
    // static information
    hts_pos_t pos, phase_set;
    int n_depth, n_uniq_bases; // ref+alt, used in variant calling & haplotype assignment
    int n_low_depth; // including bases/regions with low quality, only count depth, not base
    int8_t ref_base; // for VCF output
    uint8_t *bases; int *base_covs, *base_to_i; // 'ACGTN.D', and corresponding coverage

    // dynamic information, update during haplotype assignment
    uint8_t *base_to_hap; // SNP-wise (base_i_to_hap): 0123456(ACGTN.D) -> 1:H1/2:H2/0:not set yet
    int **hap_to_base_profile; // read-wise: 1:H1/2:H2 -> base_i -> read count
    int *hap_to_cons_base; // HAP-wise (hap_to_cons_base_i): 1:H1/2:H2 -> 0123456(ACGTN.D)

    // XXX
    uint8_t is_low_qual; // read evidence is low quality to support this site as SNP
                         // will be skipped or need extra operations to confirm
    uint8_t is_skipped; // skipped in VCF output
} cand_snp_t;

// read X snp
typedef struct read_snp_profile_t {
    int read_id; // 0 .. bam_chunk->n_read-1
    int start_snp_idx, end_snp_idx; // 0 .. n_total_cand_snps-1
    uint8_t *snp_is_used; // size: n_total_cand_snps
    uint8_t *snp_bases; // 0123456: ACGTN.D, set ref allele as '.' if ref seq is NA, D: deletion
    uint8_t *snp_qual;
} read_snp_profile_t;

typedef struct var_site_t {
    hts_pos_t pos; int var_type; // var_type: BAM_CINS/BAM_CDEL/BAM_CDIFF
    int ref_len; // SNP: 1, INS: 1, DEL: 1+del_len
} var_site_t;

// only include small variants, i.e., SNP, small indels
typedef struct cand_var_t {
    // static information
    hts_pos_t pos, phase_set;
    int var_type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    int n_depth, n_uniq_alles; // ref+alt, used in variant calling & haplotype assignment
    int n_low_depth; // including bases/regions with low quality, only count depth, not seq

    uint8_t *ref_seq; int ref_len; // REF in vcf file, could be NULL if not available
    uint8_t **alle_seqs; int *alle_seq_len, *alle_covs; // all candidate genotype sequences
    // XXX keep up to 2 alleles for each variant site

    // dynamic information, update during haplotype assignment
    uint8_t *alle_to_hap; // var-wise (alle_to_hap): alle_i -> 1:H1/2:H2/0:not set yet
    int **hap_to_alle_profile; // read-wise: 1:H1/2:H2 -> alle_i -> read count
    int *hap_to_cons_alle; // HAP-wise (hap_to_cons_alle_i): 1:H1/2:H2 -> alle_i

    uint8_t is_low_qual, is_skipped;
} cand_var_t;

// read X var
typedef struct read_var_profile_t {
    int read_id; // 0 .. bam_chunk->n_read-1
    int start_var_idx, end_var_idx; // 0 .. n_total_cand_vars-1
    uint8_t *var_is_used; // size: n_total_cand_vars
    // 0:ref, 1/2:alt
    uint8_t *var_bases; // 0123456: ACGTN.D, set ref allele as '.' if ref seq is NA, D: deletion
    uint8_t *var_qual;
} read_var_profile_t;

// complex variant: large SV, SNP/indel in complex region
typedef struct cand_complex_var_t {

} cand_complex_var_t;

struct bam_chunk_t;
struct var_t;
struct call_var_pl_t;

void collect_var_main(const struct call_var_pl_t *pl, struct bam_chunk_t *bam_chunk, struct var_t *var);
void stitch_var_main(const struct call_var_pl_t *pl, struct bam_chunk_t *bam_chunk, struct var_t *var, long ii);

#ifdef __cplusplus
}
#endif

#endif // LONGCALLD_COLLECT_VAR_H
