#ifndef LONGCALLD_ALIGN_H
#define LONGCALLD_ALIGN_H

#include "htslib/sam.h"

#define LONGCALLD_NOISY_NO_COVER     0x0000
#define LONGCALLD_NOISY_RIGHT_GAP    0x0001
#define LONGCALLD_NOISY_LEFT_GAP     0x0002
#define LONGCALLD_NOISY_RIGHT_COVER  0x0004 // 4, 6
#define LONGCALLD_NOISY_LEFT_COVER   0x0008 // 8, 9
#define LONGCALLD_NOISY_BOTH_COVER   0x000C // 12

#define LONGCALLD_NOISY_IS_BOTH_COVER(c)    (((c) & LONGCALLD_NOISY_LEFT_COVER) && ((c) & LONGCALLD_NOISY_RIGHT_COVER))
#define LONGCALLD_NOISY_IS_LEFT_COVER(c)    (((c) & LONGCALLD_NOISY_LEFT_COVER) && ((c) & LONGCALLD_NOISY_RIGHT_COVER) == 0)
#define LONGCALLD_NOISY_IS_LEFT_GAP(c)      ((c) & LONGCALLD_NOISY_LEFT_GAP)
#define LONGCALLD_NOISY_IS_RIGHT_COVER(c)   (((c) & LONGCALLD_NOISY_LEFT_COVER) == 0 && ((c) & LONGCALLD_NOISY_RIGHT_COVER))
#define LONGCALLD_NOISY_IS_RIGHT_GAP(c)     ((c) & LONGCALLD_NOISY_RIGHT_GAP)
#define LONGCALLD_NOISY_IS_NOT_COVER(c)     (((c) & LONGCALLD_NOISY_LEFT_COVER) == 0 && ((c) & LONGCALLD_NOISY_RIGHT_COVER) == 0)


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


#define LONGCALLD_WFA_NO_HEURISTIC 0
#define LONGCALLD_WFA_ADAPTIVE 1
#define LONGCALLD_WFA_ZDROP 2

#define LONGCALLD_WFA_AFFINE_1P 0
#define LONGCALLD_WFA_AFFINE_2P 1

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;
struct aln_str_t;
struct cand_var_t;

int collect_te_info_from_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *var);
int collect_te_info_from_cons(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t gap_ref_start, int msa_gap_start, int var_type, int gap_len, uint8_t *cons_msa_seq, 
                              uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev);
int edlib_end2end_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid);
int edlib_infix_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid);
// int wfa_aln(int gap_pos, char *pattern, int plen, char *text, int tlen, uint32_t **cigar_buf);
int end2end_aln(const call_var_opt_t *opt, char *pattern, int plen, uint8_t *text, int tlen, uint32_t **cigar_buf);
int wfa_end2end_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen,
                    int gap_aln, int a, int b, int q, int e, int q2, int e2, int use_heuristic,
                    uint32_t **cigar_buf, int *cigar_length, uint8_t **pattern_alg, uint8_t **text_alg, int *alg_length);
int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen, int a, int b, int q, int e, int q2, int e2, int *n_eq, int *n_xid);

int collect_noisy_reg_aln_strs(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, 
                               int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reads, uint8_t *ref_seq, int ref_seq_len,
                               int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs);
int wfa_collect_diff_ins_seq(const call_var_opt_t *opt, uint8_t* large_seq, int large_len, uint8_t *small_seq, int small_len, uint8_t **diff_seq);

#ifdef __cplusplus
}
#endif

#endif
