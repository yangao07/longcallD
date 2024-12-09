#ifndef LONGCALLD_ALIGN_H
#define LONGCALLD_ALIGN_H

int test_wfa(char *pattern, char *text);
int test_abpoa(uint8_t **bseqs, int n_seqs, int *seq_lens);
int test_edlib(char *pattern, char *text);

int collect_wfa_aln(int gap_pos, char *ref_seq, int rlen, char *query, int qlen, char **ref_aln, char **query_aln);
int collect_msa_seqs(uint8_t *ref_seq, int ref_seq_len, uint8_t **cons_seqs, int *cons_lens, int n_cons, uint8_t ***msa_seqs);
int collect_reg_ref_cseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, char **ref_cseq);
int collect_reg_ref_bseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, uint8_t **ref_bseq);
int collect_noisy_cons_seqs(bam_chunk_t *chunk, int noisy_reg_i, int **cons_lens, uint8_t ***cons_seqs, int use_phase_info);
#endif
