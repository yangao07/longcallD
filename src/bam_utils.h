#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"

#define LONGCALLD_BAM_BASE_STR "ACGTN.D"

#ifdef __cplusplus
extern "C" {
#endif

// site wise
typedef struct {
    hts_pos_t pos; int n_depth; // ref+alt
    int ref_base_cov; // ref_base could be '.' if ref_seq is not used
    int n_uniq_alt_bases; 
    uint8_t *alt_bases; int *alt_base_covs; // 'ACGTN.', and corresponding coverage
    // int *HPs;
} cand_snp_t;

// read base wise
typedef struct {
    hts_pos_t pos;
    uint8_t base, qual;
} mis_base_t;

// read X snp
typedef struct {
    int read_id; // 0 .. bam_chunk->n_read-1
    int start_snp_idx, end_snp_idx; // 0 .. n_total_cand_snps-1
    uint8_t *snp_bases; // 0123456: ACGTN.D, .: ref allele if no seq is used, D: deletion
    uint8_t *snp_qual; // size: n_span_snps
} read_snp_profile_t;

// input: tid, tname, start, end, reads, n_reads
// output: HPs
typedef struct {
    int tid; char *tname; uint64_t start, end;
    uint8_t bam_has_eqx_cigar, bam_has_md_tag;
    int n_reads, m_reads;
    bam1_t **reads; uint8_t *is_skipped; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    int *HPs; // size: m_reads
} bam_chunk_t;

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD);
int get_mis_bases_from_eqx_cigar(bam1_t *read, int min_bq, mis_base_t **mis_bases);
int update_mis_sites_from_eqx_cigar(bam1_t *read, int min_bq, int n_total_mis_pos, hts_pos_t *mis_pos, cand_snp_t *mis_sites);
int update_read_snp_profile_from_eqx_cigar(bam1_t *read, int n_total_mis_pos, cand_snp_t *mis_sites, int mis_start_i, read_snp_profile_t *read_snp_profile);
int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, mis_base_t **mis_bases);
int get_mis_bases_from_ref_seq(bam1_t *read, int min_bq, kstring_t *ref_seq, mis_base_t **mis_bases);
int bam_read_chunk(samFile *in_bam, bam_hdr_t *header, bam_chunk_t *chunk, int max_reads_per_chun);
read_snp_profile_t *init_read_snp_profile(int n_reads, int n_total_snps);
void free_read_snp_profile(read_snp_profile_t *p, int n_reads);

#ifdef __cplusplus
}
#endif

#endif
