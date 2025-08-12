#ifndef _LONGCALLD_KMER_H
#define _LONGCALLD_KMER_H
#include <stdint.h>
#include <stdlib.h>
#include "khash.h"
#include "ksort.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef struct { uint32_t x, y; } kmer32_t;
typedef struct { size_t n, m; kmer32_t *a; } kmer32_v;


#define kmer32_hash(a) ((a))
#define kmer32_eq(a, b) ((a) == (b)) // key: up to 32bits (16-mer)
KHASH_INIT(kmer32, uint32_t, uint32_t, 1, kmer32_hash, kmer32_eq)
typedef khash_t(kmer32) kmer32_hash_t;



struct call_var_opt_t;
int make_te_kmer_idx(struct call_var_opt_t *opt);
int check_te_seq(const struct call_var_opt_t *opt, uint8_t *cand_te_seq, int cand_te_len, int *is_rev);
int test_te_kmer_query(struct call_var_opt_t *opt, char *query_fn);

#ifdef __cplusplus
}
#endif

#endif