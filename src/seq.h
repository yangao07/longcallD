#ifndef LONGCALLD_SEQ_H
#define LONGCALLD_SEQ_H
#include <stdint.h>
#include <zlib.h>
#include <assert.h>
#include "kseq.h"
#include "khash.h"
#include "cgranges.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
KSEQ_INIT(gzFile, gzread)

extern unsigned char nst_nt4_table[256];
extern unsigned char com_nst_nt4_table[256];
extern const uint8_t hash_nt4_table[6];
extern const uint8_t rc_nt4_table[6];
extern char n_char[6];

KHASH_MAP_INIT_STR(str, uint32_t)

#ifdef __cplusplus
extern "C" {
#endif

typedef struct {
    kstring_t name, seq, comment, qual;
} seq_t;

typedef struct {
    int n, m;
    kstring_t *seq, *name;
    khash_t(str) *h;
} ref_seq_t;

typedef struct {
    // int tid; 
    hts_pos_t beg, end; // 1-based, [beg, end]
    kstring_t seq;
} ref_reg_seq1_t;

typedef struct {
    int n, m;
    ref_reg_seq1_t *reg_seq; kstring_t *name;
    cgranges_t *reg_cr;
} ref_reg_seq_t;

uint32_t hash_key(uint8_t *bseq, int seq_len);
uint32_t hash_shift_key(uint32_t pre_key, uint8_t *bseq, int pre_i, int cur_i, int k);
uint8_t *get_bseq(char *seq, int seq_len);
char *get_rc_seq(char *seq, int seq_len);

ref_seq_t *ref_seq_init();
ref_seq_t *ref_seq_realloc(ref_seq_t *r);
ref_seq_t *read_ref_seq(const char *ref_fa_fn);
void ref_seq_free(ref_seq_t *r);

ref_reg_seq_t *ref_reg_seq_init(void);
ref_reg_seq_t *ref_reg_seq_realloc(ref_reg_seq_t *r);
void ref_reg_seq_free(ref_reg_seq_t *r);
void read_ref_reg_seq1(faidx_t *fai, ref_reg_seq_t *r, const char *rname, hts_pos_t beg, hts_pos_t end);
ref_reg_seq_t *read_ref_reg_seq(const char *ref_fa_fn);
// char get_ref_base_from_cr(ref_reg_seq_t *ref_seq, cgranges_t *chunk_cr, hts_pos_t pos);
// char *get_ref_bases_from_cr(ref_reg_seq_t *ref_seq, cgranges_t *chunk_cr, hts_pos_t beg, hts_pos_t end);

int ref_seq_name2id(ref_seq_t *r, const char *rname);
// int ref_seq_get_chr_i(ref_seq_t *r, const char *rname);

// char get_ref_base_from_cr(ref_reg_seq_t *ref_seq, cgranges_t *chunk_cr, hts_pos_t pos) {
//     int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
//     ovlp_n = cr_overlap(chunk_cr, "cr", pos-1, pos, &ovlp_b, &max_b);
//     assert(ovlp_n == 1);
//     int i = cr_label(chunk_cr, ovlp_b[0]);
//     assert(i < 0 || i >= ref_seq->n);
//     ref_reg_seq1_t *reg_seq = ref_seq->reg_seq+i;
//     assert(pos > reg_seq->beg && pos <= reg_seq->end);
//     char ref_base = reg_seq->seq.s[pos - reg_seq->beg - 1];
//     free(ovlp_b);
//     return ref_base;
// }

// // beg, end: 1-based [beg, end]
// static inline char *get_ref_bases_from_cr(ref_reg_seq_t *ref_seq, cgranges_t *chunk_cr, hts_pos_t beg, hts_pos_t end) {
//     int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
//     ovlp_n = cr_overlap(chunk_cr, "cr", beg, end, &ovlp_b, &max_b);
//     assert(ovlp_n == 1);
//     int i = cr_label(chunk_cr, ovlp_b[0]);
//     assert(i < 0 || i >= ref_seq->n);
//     ref_reg_seq1_t *reg_seq = ref_seq->reg_seq+i;
//     assert(beg > reg_seq->beg && end <= reg_seq->end);
//     char *ref_bases = reg_seq->seq.s + beg - 1 - reg_seq->beg;
//     free(ovlp_b);
//     return ref_bases;
// }

#ifdef __cplusplus
}
#endif

#endif