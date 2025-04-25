#ifndef LONGCALLD_ALIGN_H
#define LONGCALLD_ALIGN_H

#include "htslib/sam.h"

#define LONGCALLD_MATCH_SCORE 2
#define LONGCALLD_MISMATCH_SCORE 6 // 4
#define LONGCALLD_GAP_OPEN1_SCORE 6 // 4
#define LONGCALLD_GAP_EXT1_SCORE 2
#define LONGCALLD_GAP_OPEN2_SCORE 24
#define LONGCALLD_GAP_EXT2_SCORE 1

#define LONGCALLD_GAP_LEFT_ALN 1 // default: put gap at the left-most position
#define LONGCALLD_GAP_RIGHT_ALN 2

#define LONGCALLD_EXT_ALN_LEFT_TO_RIGHT 1
#define LONGCALLD_EXT_ALN_RIGHT_TO_LEFT  2

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct aln_str_t;
struct cand_var_t;

int collect_te_info_from_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *var);
int collect_te_info_from_cons(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t gap_ref_start, int msa_gap_start, int var_type, int gap_len, uint8_t *cons_msa_seq, 
                              uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev);
// int wfa_aln(int gap_pos, char *pattern, int plen, char *text, int tlen, uint32_t **cigar_buf);
int end2end_aln(const call_var_opt_t *opt, char *pattern, int plen, uint8_t *text, int tlen, uint32_t **cigar_buf);
int collect_reg_ref_cseq(bam_chunk_t *chunk, hts_pos_t *reg_beg, hts_pos_t *reg_end, char **ref_cseq);
int collect_reg_ref_bseq(bam_chunk_t *chunk, hts_pos_t *reg_beg, hts_pos_t *reg_end, uint8_t **ref_bseq);
int collect_noisy_reg_aln_strs(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, 
                               int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reads, uint8_t *ref_seq, int ref_seq_len,
                               int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs);

#ifdef __cplusplus
}
#endif

#endif
