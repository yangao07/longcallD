#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"

#define LONGCALLD_BAM_CHUNK_SIZE 4000 // 512
                           ///  012345678
#define LONGCALLD_BAM_BASE_STR "ACGTN.DIU"
#define LONGCALLD_BAM_BASE_N 9
#define LONGCALLD_BAM_REF_BASE '.'
#define LONGCALLD_BAM_REF_BASE_IDX 5
#define LONGCALLD_BAM_DEL_BASE 'D'
#define LONGCALLD_BAM_DEL_BASE_IDX 6
#define LONGCALLD_BAM_INS_BASE 'I'
#define LONGCALLD_BAM_INS_BASE_IDX 7
#define LONGCALLD_BAM_UNSET_BASE 'U'
#define LONGCALLD_BAM_UNSET_BASE_IDX 8

#ifdef __cplusplus
extern "C" {
#endif

// site wise summary data
typedef struct {
    // static information
    hts_pos_t pos, phase_set; int n_depth, n_uniq_bases; // ref+alt
    uint8_t *bases; int *base_covs, *base_to_i; // 'ACGTN.D', and corresponding coverage

    // dynamic information, update during haplotype assignment
    uint8_t *base_to_hap; // SNP-wise (base_i_to_hap): 0123456(ACGTN.D) -> 1:H1/2:H2/0:not set yet
    int **hap_to_base_profile; // read-wise: 1:H1/2:H2 -> base_i -> read count
    int *hap_to_cons_base; // HAP-wise (hap_to_cons_base_i): 1:H1/2:H2 -> 0123456(ACGTN.D)
} cand_snp_t;

// read/base wise X/I/D operations from CIGAR

typedef struct {
    hts_pos_t pos; int type, len, qi; // pos: 1-based ref position, qi: 0-based query position
    uint8_t base, qual;  // only for X
    // two rounds:
    // 1st: left_low_len, right_low_len 
    // 2nd: if len > left+right, then left+mid+right (split)
    //                           else all_low (no split)
    uint8_t is_low_qual; // low-qual sites will not be used for selecting candidate SNPs and haplotype assignment
} digar1_t;              // they may be: 1) in dense X/gap region
                         //              2) low base quality
                         //              3) X around indels
                         //              4) in repetitive region (tend to be dense X/gap region as well)
                         //              5) wrong mapping (tend to be dense X/indel region as well)
                         //              etc.
typedef struct {
    int n_digar, m_digar;
    digar1_t *digars;
} digar_t; // detailed CIGAR for each read

// typedef struct {
//     hts_pos_t pos;
//     uint8_t base, qual;          // only for X
//     uint8_t is_low_qual; // low-qual sites will not be used for selecting candidate SNPs and haplotype assignment
// } x_base_t;              // they may be: 1) in dense X/gap region
                         //              2) low base quality
                         //              3) X around indels
                         //              4) in repetitive region (tend to be dense X/gap region as well)
                         //              5) wrong mapping (tend to be dense X/indel region as well)
                         //              etc.


// read X snp
typedef struct {
    int read_id; // 0 .. bam_chunk->n_read-1
    int start_snp_idx, end_snp_idx; // 0 .. n_total_cand_snps-1
    uint8_t *snp_is_used; // size: n_total_cand_snps
    uint8_t *snp_bases; // 0123456: ACGTN.D, .: ref allele if no seq is used, D: deletion
    uint8_t *snp_qual;
} read_snp_profile_t;

struct call_var_opt_t;

typedef struct {
    // input
    int tid; char *tname; hts_pos_t start, end;
    uint8_t bam_has_eqx_cigar, bam_has_md_tag;
    int n_reads, m_reads;
    bam1_t **reads; 
    // intermidiate
    digar_t *digars;
    uint8_t *is_skipped; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    // output
    int *haps; // size: m_reads
} bam_chunk_t;

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD);
void collect_digar_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_MD_tag(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_ref_seq(bam1_t *read, const struct call_var_opt_t *opt, kstring_t *ref_seq, digar_t *digar);
int update_cand_snps_from_digar(digar_t *digar, int n_x_sites, hts_pos_t *x_sites, int start_i, cand_snp_t *cand_snps);
int update_read_snp_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_snps, cand_snp_t *cand_snps, int start_snp_i, read_snp_profile_t *read_snp_profile);
// int get_x_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, x_base_t **mis_bases);
// int update_x_sites_from_eqx_cigar(bam1_t *read, int min_bq, int n_total_mis_pos, hts_pos_t *mis_pos, cand_snp_t *mis_sites);
// int update_read_snp_profile_from_eqx_cigar(bam1_t *read, int n_total_mis_pos, cand_snp_t *mis_sites, int mis_start_i, read_snp_profile_t *read_snp_profile);


int bam_read_chunk(samFile *in_bam, bam_hdr_t *header, bam_chunk_t *chunk, int max_reads_per_chun);
void bam_read_chunk_free(bam_chunk_t *chunk);
read_snp_profile_t *init_read_snp_profile(int n_reads, int n_total_snps);
void free_read_snp_profile(read_snp_profile_t *p, int n_reads);

// int get_xid_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, xid_t **mis_bases);

// int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, x_base_t **mis_bases);

#ifdef __cplusplus
}
#endif

#endif