#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "cgranges.h"
#include "collect_var.h"

#define LONGCALLD_BAM_CHUNK_READ_COUNT 4000
#define LONGCALLD_BAM_CHUNK_REG_SIZE 500000 // 0.5M/1M


#define BAM_RECORD_LOW_QUAL 0x1
#define BAM_RECORD_WRONG_MAP 0x2
// #define BAM_RECORD_LARGE_CLIP 0x4

#define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

// read/base wise X/I/D operations from CIGAR
typedef struct digar1_t {
    hts_pos_t pos; int type, len, qi; // pos: 1-based ref position, qi: 0-based query position
    uint8_t *alt_seq; // only used for mismatch/insertion, deletion:NULL
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
    hts_pos_t beg, end; // [beg, end]
    int n_digar, m_digar;
    digar1_t *digars;
    // read-wise noisy region: active region for re-alignment
    uint8_t has_long_clip; // XXX long clip in the read, need to re-align see if there is a long ins/del, if not, probabely wrong mapping, skip this read
    cgranges_t *noisy_regs; // merge low_qual digar1_t if they are next to each other
    const uint8_t *bseq, *qual;
} digar_t; // detailed CIGAR for each read

// coordinate system: ACGT: 1234
//   samtools view/longCallD call: 1-based, [beg, end]: [1,4]
//   bam_chunk beg/end, ref_beg/ref_end: 1-based, [beg, end]:[1,4]
//   ref_reg_seq_t: 1-based, [beg, end]: [1,4]
//   cgranges: 0-based, (beg, end]: (0, 4]
// a chunk of BAM reads, and auxiliary information for variant calling/phasing
// the very basic working unit where variants are called and haplotypes are assigned
typedef struct bam_chunk_t {
    // input
    int tid; char *tname;
    // ref_seq: 
    //   reference sequence for this chunk, size: ref_end-ref_beg
    //   if reg_reg_cr contain >1 regions, ref_seq will be concatenated, including additonal gaps between regions
    char *ref_seq; uint8_t need_free_ref_seq; hts_pos_t ref_beg, ref_end; // [ref_beg, ref_end]
    // reg_cr: 
    //   reference region for this chunk. usually just 1 region, but could be multiple regions when regions are small, or reads are long
    //   only variants within regions will be considered, variants outside will be skipped
    //   for variants across boundaries/spaning multiple regions, they will be processed during stitching
    //   ref_beg <= reg_beg < reg_end <= ref_end, ref_seq may include additional flanking regions (1kb)
    cgranges_t *reg_cr, *low_comp_cr; hts_pos_t reg_beg, reg_end; // [reg_beg, reg_end]
    // hts_pos_t active_reg_beg, active_reg_end; // [active_reg_beg, active_reg_end]: only variants within this region will be considered
                                              // variants spaning this region will be re-processed during stitching
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
    int *var_i_to_cate; // size: n_cand_vars; HET/HOM/SOMATIC
    // XXX noisy region: 
    //  1) too many XIDs in non-repeat regions 
    //  2) higher depth than nearby normal regions
    //  3) normal reads with similar depth exist
    cgranges_t *chunk_noisy_regs; // merged noisy regions for all reads
    int **noisy_reg_to_reads, *noisy_reg_to_n_reads; // size: chunk->chunk_noisy_regs->n_r

    read_var_profile_t *read_var_profile; cgranges_t *read_var_cr;
    // output
    uint8_t flip_hap; // XXX flip haplotype for this chunk, chromosome-wise global parameter, flip_hap ^= pre_flip
                      // so all the haplotypes need to be output sequentially to keep consistency
    int *haps; hts_pos_t *PS; // size: m_reads
} bam_chunk_t;

struct call_var_pl_t;
struct call_var_opt_t;
struct cand_var_t;
struct var_site_t;
struct read_var_profile_t;

static inline int is_in_noisy_reg(hts_pos_t pos, cgranges_t *noisy_regs) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(noisy_regs, "cr", pos, pos+1, &ovlp_b, &max_b);
    free(ovlp_b);
    return (ovlp_n > 0);
}

static inline int is_overlap_cr(char *tname, hts_pos_t beg, hts_pos_t end, cgranges_t *reg_cr) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(reg_cr, tname, beg-1, end, &ovlp_b, &max_b);
    free(ovlp_b);
    return (ovlp_n > 0);
}

static inline int is_overlap_reg(hts_pos_t beg, hts_pos_t end, hts_pos_t reg_beg, hts_pos_t reg_end) {
    if (beg > reg_end || end < reg_beg) return 0;
    return 1;
}

int get_aux_int_from_bam(bam1_t *b, const char *tag);
char *get_aux_str_from_bam(bam1_t *b, const char *tag);
void print_digar1(digar1_t *digar, int n_digar, FILE *fp);
void print_digar(digar_t *digar, FILE *fp);
int collect_reg_digars_var_seqs(bam_chunk_t *chunk, int read_i, hts_pos_t reg_beg, hts_pos_t reg_end, digar1_t *reg_digars, uint8_t **reg_var_seqs, int *fully_cover);
void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD);
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_MD_tag(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_ref_seq(bam_chunk_t *chunk, bam1_t *read, const struct call_var_pl_t *pl, const struct call_var_opt_t *opt, digar_t *digar);
int update_cand_vars_from_digar(digar_t *digar, bam1_t *read, int n_var_sites, struct var_site_t *var_sites, int start_i, struct cand_var_t *cand_vars);
void update_read_var_profile_with_allele(int var_i, int allele_i, read_var_profile_t *read_var_profile);
int update_read_var_profile_from_digar(digar_t *digar, bam1_t *read, int n_cand_vars, struct cand_var_t *cand_vars, int start_var_i, struct read_var_profile_t *read_var_profile);

int collect_bam_chunk(struct call_var_pl_t *pl, bam_chunk_t *chunk, int chunk_i);
void copy_bam_chunk0(bam_chunk_t *from_chunk, bam_chunk_t *to_chunk);
void bam_chunk_free(bam_chunk_t *chunk);
void bam_ovlp_chunk_free(bam_chunk_t *chunk);
void bam_chunks_free(bam_chunk_t *chunks, int n_chunks);
struct read_var_profile_t *init_read_var_profile(int n_reads, int n_total_vars);
struct read_var_profile_t *init_read_var_profile_with_ids(int n_reads, int *read_ids, int n_total_vars);
void free_read_var_profile(struct read_var_profile_t *p, int n_reads);
char *extract_sample_name_from_bam_header(bam_hdr_t *header);

// int get_xid_from_eqx_cigar(bam1_t *read, const struct call_var_opt_t *opt, xid_t **mis_bases);

// int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, x_base_t **mis_bases);

// (beg, end]
static inline int is_in_reg(cgranges_t *reg_cr, char *tname, hts_pos_t beg, hts_pos_t end) {
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    ovlp_n = cr_overlap(reg_cr, tname, beg, end, &ovlp_b, &max_b);
    free(ovlp_b);
    // fprintf(stderr, "is_in_reg: %s:%d-%d : %d\n", tname, beg, end, ovlp_n);
    // for (int i = 0; i < reg_cr->n_r; ++i) {
        // fprintf(stderr, "reg_cr: %s:%d-%d\n", reg_cr->ctg[cr_label(reg_cr, i)].name, cr_start(reg_cr, i), cr_end(reg_cr, i));
    // }
    return (ovlp_n > 0);
}

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_BAM_UTILS_H