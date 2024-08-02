#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"

#define LONGCALLD_BAM_CHUNK_READ_COUNT 4000
#define LONGCALLD_BAM_CHUNK_REG_SIZE 1000000 // 1M
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

#define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

// read/base wise X/I/D operations from CIGAR
typedef struct {
    hts_pos_t pos; int type, len, qi; // pos: 1-based ref position, qi: 0-based query position
    // two rounds:
    // 1st: left_low_len, right_low_len 
    // 2nd: if len > left+right, then left+mid+right (split)
    //                           else all_low (no split)
    uint8_t is_low_qual; // low-qual sites will not be used for haplotype assignment,
} digar1_t;              // they may be: 1) in dense X/gap region
                         //                 a) in repetitive region
                         //                 b) wrong mapping
                         //              2) low base quality
                         //              3) X around indels
                         //              etc.
typedef struct {
    int n_digar, m_digar;
    digar1_t *digars;
    const uint8_t *bseq, *qual;
} digar_t; // detailed CIGAR for each read

typedef struct bam_chunk_t {
    // input
    int tid; char *tname; hts_pos_t beg, end;
    uint8_t bam_has_eqx_cigar, bam_has_md_tag;
    int n_reads, m_reads, n_ovlp_reads;
    bam1_t **reads; 
    // intermidiate
    digar_t *digars;
    uint8_t *is_skipped; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    // output
    uint8_t flip_hap; // XXX flip haplotype for this chunk, chromosome-wise global parameter, flip_hap ^= pre_flip
                      // so all the haplotypes need to be output sequentially to keep consistency
    int *haps; // size: m_reads
} bam_chunk_t;

struct call_var_opt_t;
struct cand_snp_t;
struct cand_var_t;
struct var_site_t;
struct read_snp_profile_t;

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD);
void collect_digar_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_MD_tag(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_ref_seq(bam1_t *read, const struct call_var_opt_t *opt, kstring_t *ref_seq, digar_t *digar);
int update_cand_snps_from_digar(digar_t *digar, bam1_t *read, int n_x_sites, hts_pos_t *x_sites, int start_i, struct cand_snp_t *cand_snps);
int update_cand_vars_from_digar(digar_t *digar, bam1_t *read, int n_var_sites, struct var_site_t *var_sites, int start_i, struct cand_var_t *cand_vars);
int update_read_snp_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_snps, struct cand_snp_t *cand_snps, int start_snp_i, struct read_snp_profile_t *read_snp_profile);

int collect_bam_chunk(samFile *in_bam, bam_hdr_t *header, hts_itr_t *iter, int use_iter, int max_reg_len_per_chunk, int **ovlp_read_i, int *n_ovlp_reads, bam_chunk_t *chunk);
void bam_chunk_free(bam_chunk_t *chunk);
struct read_snp_profile_t *init_read_snp_profile(int n_reads, int n_total_snps);
void free_read_snp_profile(struct read_snp_profile_t *p, int n_reads);
char *extract_sample_name_from_bam_header(bam_hdr_t *header);

// int get_xid_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, xid_t **mis_bases);

// int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, x_base_t **mis_bases);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_BAM_UTILS_H