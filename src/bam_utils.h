#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "cgranges.h"

#define LONGCALLD_BAM_CHUNK_READ_COUNT 4000
#define LONGCALLD_BAM_CHUNK_REG_SIZE 1000000 // 1M XXX smaller or bigger?
#define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

// include alignment-based simple variants, i.e., SNP, indels
// assume alignment is correct, no complex variants
// will be used for variant calling, haplotype assignment, phasing of variants and reads
typedef struct cand_var_t {
    // static information
    int tid; hts_pos_t pos, phase_set;
    // XXX check if is_rep_reg
    int is_rep_reg; hts_pos_t non_rep_beg, non_rep_end; // non-repetitive window, used for alignment-based variant calling 
    int var_type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    int n_depth; // ref+alt, used in variant calling & haplotype assignment, excluding low-qual bases
    int n_low_depth; // including bases/regions with low quality, only count depth, not seq
    int n_uniq_alles; // up ot 4: ref, alt1, alt2, minor_alt
                      // snp/ins: could be >2, del: ≤2
                      // minor_alt: not ref and not main alt alleles (mostly sequencing errors)
    int *alle_covs; // size: n_uniq_alles
    int ref_len, *alt_lens; uint8_t **alt_seqs; // all candidate alt alleles, ''/0 for minor_alt
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
    int *alleles; // 0:ref, 1/2:alt
} read_var_profile_t;

// read/base wise X/I/D operations from CIGAR
typedef struct digar1_t {
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
    // read-wise noisy region: active region for re-alignment
    cgranges_t *noisy_regs; // merge low_qual digar1_t if they are next to each other
    const uint8_t *bseq, *qual;
} digar_t; // detailed CIGAR for each read

// a chunk of BAM reads, and auxiliary information for variant calling/phasing
// the very basic unit where variants are called and haplotypes are assigned
typedef struct bam_chunk_t {
    // input
    int tid; char *tname; hts_pos_t beg, end; // for each chunk, only variants between (beg, end] will be considered
    uint8_t bam_has_eqx_cigar, bam_has_md_tag;
    int n_reads, m_reads;
    int n_up_ovlp_reads; // number of reads overlapping with upstream bam chunk
                         // for each chunk, only output [n_up_ovlp_reads, n_reads]'s reads
    int *up_ovlp_read_i; // read_i in upstream bam chunk; size: n_up_ovlp_reads
    bam1_t **reads; 
    // intermidiate
    uint8_t *is_ovlp; // size: m_reads, is_ovlp: overlap with other chunks
    uint8_t *is_skipped; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    digar_t *digars;
    // variant-related
    int n_cand_vars; cand_var_t *cand_vars;
    int **var_cate_idx, *var_cate_counts; // size: LONGCALLD_VAR_CATE_N
    int *var_i_to_cate; // size: n_cand_vars
    read_var_profile_t *read_var_profile; cgranges_t *read_var_cr;
    // output
    uint8_t flip_hap; // XXX flip haplotype for this chunk, chromosome-wise global parameter, flip_hap ^= pre_flip
                      // so all the haplotypes need to be output sequentially to keep consistency
    int *haps; // size: m_reads
} bam_chunk_t;

struct call_var_pl_t;
struct call_var_opt_t;
struct cand_var_t;
struct var_site_t;
struct read_var_profile_t;

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD);
void collect_digar_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_MD_tag(bam1_t *read, const struct call_var_opt_t *opt, digar_t *digar);
void collect_digar_from_ref_seq(bam1_t *read, const struct call_var_opt_t *opt, kstring_t *ref_seq, digar_t *digar);
int update_cand_vars_from_digar(digar_t *digar, bam1_t *read, int n_var_sites, struct var_site_t *var_sites, int start_i, struct cand_var_t *cand_vars);
int update_read_var_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_vars, struct cand_var_t *cand_vars, int start_var_i, struct read_var_profile_t *read_var_profile);

int collect_bam_chunk(struct call_var_pl_t *pl, int **last_chunk_read_i, int *n_last_chunk_reads, hts_pos_t *cur_reg_beg, bam_chunk_t *chunk);
void bam_chunk_free(bam_chunk_t *chunk);
void bam_chunks_free(bam_chunk_t *chunks, int n_chunks);
struct read_var_profile_t *init_read_var_profile(int n_reads, int n_total_vars);
void free_read_var_profile(struct read_var_profile_t *p, int n_reads);
char *extract_sample_name_from_bam_header(bam_hdr_t *header);

// int get_xid_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, xid_t **mis_bases);

// int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, x_base_t **mis_bases);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_BAM_UTILS_H