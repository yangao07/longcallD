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

#define LONGCALLD_EXT_ALN_RIGHT 1
#define LONGCALLD_EXT_ALN_LEFT  2

#ifdef __cplusplus
extern "C" {
#endif

struct bam_chunk_t;

int test_wfa(char *pattern, char *text);
int test_abpoa(uint8_t **bseqs, int n_seqs, int *seq_lens);
int test_edlib(char *pattern, char *text);
int test_ksw2(char *pattern, char *text, const call_var_opt_t *opt);

int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen, int *n_eq, int *n_xid);

// int wfa_aln(int gap_pos, char *pattern, int plen, char *text, int tlen, uint32_t **cigar_buf);
int end2end_aln(const call_var_opt_t *opt, char *pattern, int plen, uint8_t *text, int tlen, uint32_t **cigar_buf);
int collect_ref_cons_msa_seqs(const call_var_opt_t *opt, uint8_t *ref_seq, int ref_seq_len, uint8_t **cons_seqs, int *cons_lens, int n_cons, uint8_t ***msa_seqs);
int collect_reg_ref_cseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, char **ref_cseq);
int collect_reg_ref_bseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, uint8_t **ref_bseq);
int collect_noisy_cons_seqs(bam_chunk_t *chunk, int noisy_reg_i, int **cons_lens, uint8_t ***cons_seqs);
int collect_noisy_reg_cons_seqs0(const call_var_opt_t *opt, bam_chunk_t *chunk, int noisy_reg_i, int *cons_lens, uint8_t **cons_seqs);
int collect_noisy_reg_msa_cons(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, 
                               int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reads, uint8_t *ref_seq, int ref_seq_len,
                               int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids, int *msa_len, uint8_t ***msa_seqs);

#ifdef __cplusplus
}
#endif

#endif
