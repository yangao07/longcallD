#ifndef SEQ_H
#define SEQ_H
#include <stdint.h>
#include <zlib.h>
#include "kseq.h"
#include "khash.h"
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

uint32_t hash_key(uint8_t *bseq, int seq_len);
uint32_t hash_shift_key(uint32_t pre_key, uint8_t *bseq, int pre_i, int cur_i, int k);
uint8_t *get_bseq(char *seq, int seq_len);
char *get_rc_seq(char *seq, int seq_len);

ref_seq_t *ref_seq_init();
ref_seq_t *ref_seq_realloc(ref_seq_t *r);
ref_seq_t *read_ref_seq(const char *ref_fa_fn);
void ref_seq_free(ref_seq_t *r);
int ref_seq_name2id(ref_seq_t *r, const char *rname);
// int ref_seq_get_chr_i(ref_seq_t *r, const char *rname);

#ifdef __cplusplus
}
#endif

#endif