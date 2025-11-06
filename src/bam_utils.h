#ifndef LONGCALLD_BAM_UTILS_H
#define LONGCALLD_BAM_UTILS_H

#include "htslib/sam.h"
#include "htslib/kstring.h"
#include "cgranges.h"
#include "collect_var.h"

#define LONGCALLD_BAM_CHUNK_READ_COUNT 4000
#define LONGCALLD_BAM_CHUNK_REG_SIZE 500000 // 0.5M/1M
// region-based BAM_CHUNK: no split for a single chromosome
//                         multi regions may be merged into a single chunk if n_regs < 64
// #define LONGCALLD_BAM_MIN_REG_CHUNK_PER_RUN 64 // n_threads * 4


#define BAM_RECORD_LOW_QUAL 0x1
#define BAM_RECORD_WRONG_MAP 0x2
// #define BAM_RECORD_LARGE_CLIP 0x4

#define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

// read/base wise X/I/D operations from CIGAR
typedef struct digar1_t {
    // pos: 1-based ref position
    // qi: 0-based query position (for =XI); for DEL, qi is the first read base after the deletion
    hts_pos_t pos; int type, len, qi;
    uint8_t *alt_seq; // X/I: alt_seq; =/D: NULL
    uint8_t is_low_qual; // low-base-qual, effective for clean-region digars
} digar1_t;

typedef struct {
    hts_pos_t beg, end; // [beg, end]
    uint8_t is_rev;
    int n_digar, m_digar;
    digar1_t *digars;
    // read-wise noisy region: active region for re-alignment
    cgranges_t *noisy_regs; // merge low_qual digar1_t if they are next to each other
    int qlen; uint8_t *bseq, *qual; // copy from bam1_t
} digar_t; // detailed CIGAR for each read

typedef struct bam_chunk_t {
    // input
    // tid = pl->reg_chunks[reg_chunk_i].reg_tids[reg_i] 
    // reg_beg = pl->reg_chunks[reg_chunk_i].reg_begs[reg_i]
    // reg_end = pl->reg_chunks[reg_chunk_i].reg_ends[reg_i]
    int reg_chunk_i, reg_i, tid; hts_pos_t reg_beg, reg_end;
    int *qual_counts; int min_qual, first_quar_qual, median_qual, third_quar_qual, max_qual; // for each read, count the number of bases with qual
    char *tname;
    // ref_seq: 
    //   reference sequence for this chunk, size: ref_end-ref_beg+1
    //   if reg_reg_cr contain >1 regions, ref_seq will be concatenated, including additonal gaps between regions
    char *ref_seq; hts_pos_t ref_beg, ref_end, whole_ref_len; // [ref_beg, ref_end]
    // reg_cr: 
    //   reference region for this chunk. usually just 1 region, but could be multiple regions when regions are small, or reads are long
    //   only variants within regions will be considered, variants outside will be skipped
    //   for variants across boundaries/spaning multiple regions, they will be processed during stitching
    //   ref_beg <= reg_beg < reg_end <= ref_end, ref_seq may include additional flanking regions (50kb)
    // should always be 1 region XXX
    cgranges_t *low_comp_cr; // tandem_rep_cr;
    int n_reads, m_reads;
    int n_bam; // for each input bam, record the number of region-overlapping reads
    int *n_up_ovlp_reads, *n_down_ovlp_reads; // number of reads overlapping with up/downstream bam chunk
    int **up_ovlp_read_i, **down_ovlp_read_i;
    bam1_t **reads;
    // intermidiate
    int *n_clean_agree_snps, *n_clean_conflict_snps; // size: m_reads; XXX include both het and hom clean vars
    uint8_t *is_ont_palindrome; // size: m_reads, 1: palindromic read, 0: non-palindromic read
    digar_t *digars; uint8_t *is_skipped, *is_skipped_for_somatic; // size: m_reads, is_skipped: wrong mapping, low qual, etc.
    // variant-related
    // n_cand_vars: including candidate germline variants and somatic variants
    // for germline variants, including clean region and noisy region
    // for somatic, only non-noisy regions, we don't call somatic variants in noisy regions (for now)
    int n_cand_vars; cand_var_t *cand_vars; int *var_i_to_cate;
    read_var_profile_t *read_var_profile; cgranges_t *read_var_cr;
    // noisy regions
    // right now: not work with cooridinate > 2G (pow(2,31)), use noisy_reg_beg-reg_beg+1 instead for coordinates > 2G
    cgranges_t *chunk_noisy_regs; // merged noisy regions for all reads
    int **noisy_reg_to_reads, *noisy_reg_to_n_reads; // size: chunk->chunk_noisy_regs->n_r
    // output
    uint8_t flip_hap; // flip haplotype for this chunk, chromosome-wise global parameter
    hts_pos_t flip_pre_PS, flip_cur_PS; // so all the haplotypes need to be output sequentially to keep consistency
    // read-wise
    // phase_score: used to determine low-qual phased reads, which should not be used for somatic variant calling
    int *phase_scores, *haps; hts_pos_t *phase_sets; // size: m_reads 
} bam_chunk_t; // reg-based bam_chunk_t

struct call_var_pl_t;
struct call_var_opt_t;
struct cand_var_t;
struct var_site_t;
struct read_var_profile_t;
struct call_var_io_aux_t;

// static inline int double_check_digar(digar_t *digar) {
static inline int double_check_digar(digar1_t *digars, int n_digar) {
    if (n_digar == 0) return 0;
    for (int i = n_digar-2; i > 0; --i) {
        int qi, last_i = i-1; int last_qi = digars[last_i].qi;
        if (digars[last_i].type == BAM_CEQUAL ||
            digars[last_i].type == BAM_CMATCH ||
            digars[last_i].type == BAM_CDIFF ||
            digars[last_i].type == BAM_CINS ||
            digars[last_i].type == BAM_CSOFT_CLIP ||
            digars[last_i].type == BAM_CHARD_CLIP) {
                qi = last_qi + digars[last_i].len;
        } else qi = last_qi;
        if (qi != digars[i].qi) {
            // fprintf(stderr, "Error: digars[%d]=%d, digars[%d]=%d\n", i, digars[i].qi, last_i, qi);
            return 1;
        }
    }
    return 0;
}

static inline int digar2qlen(digar_t *digar) {
    if (digar->n_digar == 0) return 0;
    int qlen = digar->digars[digar->n_digar-1].qi;
    if (digar->digars[digar->n_digar-1].type == BAM_CEQUAL ||
        digar->digars[digar->n_digar-1].type == BAM_CMATCH ||
        digar->digars[digar->n_digar-1].type == BAM_CDIFF ||
        digar->digars[digar->n_digar-1].type == BAM_CINS ||
        digar->digars[digar->n_digar-1].type == BAM_CSOFT_CLIP ||
        digar->digars[digar->n_digar-1].type == BAM_CHARD_CLIP) {
        qlen += digar->digars[digar->n_digar-1].len;
    }
    return qlen;
}

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

digar1_t *push_digar_alt_seq(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d);
digar1_t *push_digar0(digar1_t *digar, int *n_digar, int *m_digar, digar1_t d);
void push_digar1(digar_t *digar, digar1_t d);
void free_digar1(digar1_t *digar1, int n_digar);
int get_aux_int_from_bam(bam1_t *b, const char *tag);
char *get_aux_str_from_bam(bam1_t *b, const char *tag);
void print_digar1(digar1_t *digar, int n_digar, FILE *fp);
void print_digar(digar_t *digar, FILE *fp);
int has_equal_X_in_bam_cigar(bam1_t *read);
int has_cs_in_bam(bam1_t *b);
int has_MD_in_bam(bam1_t *b);
int collect_reg_digars_var_seqs(bam_chunk_t *chunk, int read_i, hts_pos_t reg_beg, hts_pos_t reg_end, digar1_t *reg_digars, uint8_t **reg_var_seqs, int *fully_cover);
int collect_digar_from_eqx_cigar(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_cs_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_MD_tag(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int collect_digar_from_ref_seq(bam_chunk_t *chunk, int read_i, const struct call_var_opt_t *opt, digar_t *digar);
int update_cand_vars_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_var_sites, struct var_site_t *var_sites, struct cand_var_t *cand_vars);
void update_read_var_profile_with_allele(int var_i, int allele_i, int alt_qi, read_var_profile_t *read_var_profile);
int update_read_vs_all_var_profile_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, struct cand_var_t *cand_vars, int *var_i_to_cate, struct read_var_profile_t *read_var_profile);
int update_read_vs_somatic_var_profile_from_digar(const struct call_var_opt_t *opt, bam_chunk_t *chunk, digar_t *digar, int n_cand_vars, struct cand_var_t *cand_vars, struct read_var_profile_t *read_var_profile);

int collect_ref_seq_bam_main(const struct call_var_pl_t *pl, struct call_var_io_aux_t *io_aux, int reg_chunk_i, int reg_i, bam_chunk_t *chunks);
int write_read_to_bam(bam_chunk_t *chunk, const struct call_var_opt_t *opt, const struct call_var_io_aux_t *io_aux);
void bam_chunk_mid_free(bam_chunk_t *chunk);
void bam_chunks_mid_free(bam_chunk_t *chunks, int n_chunks);
void bam_chunk_post_free(bam_chunk_t *chunk);
void bam_chunks_post_free(bam_chunk_t *chunks, int n_chunks);
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