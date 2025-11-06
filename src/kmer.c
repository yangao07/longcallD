#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include "kvec.h"
#include "kmer.h"
#include "seq.h"
#include "utils.h"
#include "call_var_main.h"
#include "sdust.h"

#define sort_key_64x(a) ((a).x)
KRADIX_SORT_INIT(64x, kmer32_t, sort_key_64x, 4) 

KSORT_INIT_GENERIC(uint32_t)

int not_simple_kmer(uint32_t kmer, int k) {
    // check if kmer is simple, i.e., AAAA..., CCCC..., GGGG..., TTTT...
    // 0: A, 1: C, 2: G, 3: T
    for (int i = 0; i < k; ++i) {
        if ((kmer & 3) != (kmer & 3) << 2 * i) return 1;
        kmer >>= 2;
    }
    return 0;
}

// collect consecutive kmers: 0, 1, 2, ...
int collect_rev_kmer(uint8_t *seq, int seq_len, int seq_id, int k, kmer32_v *a) {
    assert(k > 0 && k <= 16);
    void *km = 0; uint32_t shift1 = 2 * (k - 1);
    kv_resize(kmer32_t, km, *a, seq_len);
    // forawrd
    int i, l; uint32_t rc_hash_key = 0;
    kmer32_t m;
    for (i = l = 0; i < seq_len; ++i) {
        int c = nst_nt4_table[(int)seq[i]];
        if (c < 4) {
            rc_hash_key = rc_hash_key >> 2 | (c ^ 3) << shift1;
            l++;
            if (l >= k) {
                m.x = rc_hash_key;
                m.y = seq_id;
                if (not_simple_kmer(m.x, k)) {
                    kv_push(kmer32_t, km, *a, m);
                }
            }
        } else l = 0;
    }
    return a->n;
}

// collect consecutive kmers: 0, 1, 2, ...
int collect_kmer(uint8_t *seq, int seq_len, int seq_id, int k, kmer32_v *a) {
    assert(k > 0 && k <= 16);
    void *km = 0;
    kv_resize(kmer32_t, km, *a, seq_len);
    // forawrd
    int i, l; uint32_t hash_key = 0;
    kmer32_t m;
    for (i = l = 0; i < seq_len; ++i) {
        int c = nst_nt4_table[(int)seq[i]];
        if (c < 4) {
            hash_key = (hash_key << 2) | c;
            l++;
            if (l >= k) {
                m.x = hash_key & ((1 << 2 * k) - 1);
                m.y = seq_id;
                if (not_simple_kmer(m.x, k)) {
                    kv_push(kmer32_t, km, *a, m);
                }
            }
        } else l = 0;
    }
    return a->n;
}

// only keey kmers coming from one sequence
kmer32_hash_t *make_kmer_hash_tables(char *seq, int seq_len, int seq_id, int k, int rev) {
    kmer32_v a = {0, 0, 0};
    if (rev) {
        collect_rev_kmer((uint8_t*)seq, seq_len, seq_id, k, &a);
    } else {
        collect_kmer((uint8_t*)seq, seq_len, seq_id, k, &a);
    }
    if (a.n <= 0) return 0;
    // sort kmer
    radix_sort_64x(a.a, a.a + a.n);
    int j, n, n_keys;
    for (j = 1, n = 1, n_keys = 0; j <= a.n; ++j) {
		if (j == a.n || a.a[j].x != a.a[j-1].x) {
            // if (n == 1) 
            ++n_keys;
            n = 1;
		} else ++n;
	}
    kmer32_hash_t *h = kh_init(kmer32);
	kh_resize(kmer32, h, n_keys);

    // create hash table
    size_t start_a;
    for (j = 1, n = 1, start_a = 0; j <= a.n; ++j) {
        if (j == a.n || a.a[j].x != a.a[j-1].x) {
            // if (n == 1) {
                khint_t itr;
                int absent;
                kmer32_t *p = &a.a[j-1];
                itr = kh_put(kmer32, h, p->x, &absent);
                assert(absent && j == start_a + n);
                // kh_key(h, itr) |= 1;
                kh_val(h, itr) = p->y;
            // }
            start_a = j, n = 1;
        }
        else ++n;
    }
    free(a.a);
    // kh_destroy(kmer32, h);
    return h;
}

int make_te_kmer_idx(call_var_opt_t *opt) {
    gzFile f = gzopen(opt->te_seq_fn, "r");
    if (f == 0) {
        _err_error_exit("Cannot open TE sequence file: %s\n", opt->te_seq_fn);
    }
    int m_te_seqs = 3; // Alu, L1, SVA
    opt->te_seq_names = (char**)malloc(m_te_seqs*sizeof(char*));
    opt->te_for_h = (kmer32_hash_t**)malloc(m_te_seqs*sizeof(kmer32_hash_t*));
    opt->te_rev_h = (kmer32_hash_t**)malloc(m_te_seqs*sizeof(kmer32_hash_t*));
    kseq_t *ks = kseq_init(f);
    kseq_rewind(ks);
    int n_te_seqs = 0;
    while (kseq_read(ks) >= 0) {
        if (n_te_seqs >= m_te_seqs) {
            m_te_seqs = m_te_seqs * 2;
            opt->te_seq_names = (char**)realloc(opt->te_seq_names, m_te_seqs*sizeof(char*));
            opt->te_for_h = (kmer32_hash_t**)realloc(opt->te_for_h, m_te_seqs*sizeof(kmer32_hash_t*));
            opt->te_rev_h = (kmer32_hash_t**)realloc(opt->te_rev_h, m_te_seqs*sizeof(kmer32_hash_t*));
        }
        opt->te_seq_names[n_te_seqs] = strdup(ks->name.s);
        opt->te_for_h[n_te_seqs] = make_kmer_hash_tables(ks->seq.s, ks->seq.l, n_te_seqs, opt->te_kmer_len, 0);
        opt->te_rev_h[n_te_seqs] = make_kmer_hash_tables(ks->seq.s, ks->seq.l, n_te_seqs, opt->te_kmer_len, 1);
        n_te_seqs++;
    }
    opt->n_te_seqs = n_te_seqs;
    kseq_destroy(ks);
    gzclose(f);
    return 0;
}

// collect non-consecutive kmers: 0, k, 2k, ...
int collect_query_kmer(uint8_t *seq, int seq_len, int seq_id, int k, kmer32_v *a) {
    assert(k > 0 && k <= 16);
    void *km = 0;
    kv_resize(kmer32_t, km, *a, seq_len/k+1);
    // forawrd
    int i, l; uint32_t hash_key = 0;
    kmer32_t m;
    for (i = l = 0; i < seq_len; ++i) {
        int c = nst_nt4_table[(int)seq[i]];
        if (c < 4) {
            hash_key = (hash_key << 2) | c;
            l++;
            // only keep kmer at position 0, k, 2k, ...
            if (l == k) {
                m.x = hash_key & ((1 << 2 * k) - 1);
                m.y = seq_id;
                if (not_simple_kmer(m.x, k)) kv_push(kmer32_t, km, *a, m);
                l = 0;
            }
        } else l = 0;
    }
    return a->n;
}

int collect_kmer_hist(kmer32_v *read_kmers, kmer32_hash_t *h) {
    // collect kmer histogram, counts[]
    int count = 0;
    for (int i = 0; i < read_kmers->n; ++i) {
        kmer32_t *p = &read_kmers->a[i];
        khint_t itr = kh_get(kmer32, h, p->x);
        if (itr != kh_end(h)) {
            count++;
            // int ref_i = kh_val(h, itr);
            // assert(ref_i >= 0 && ref_i < n_ref_seqs);
            // counts[ref_i]++;
        }
    }
    return count;
}

// for somatic/mosaic SVs
// int check_tandem_dup(uint8_t *cand_query_seq, int cand_qlen, uint8_t *cand_target_seq, int cand_tlen) {
//     int q_count = 0, t_count = 0;
//     // collect all k-mers from cand_target_seq
//     // collect query_kmers from cand_query_seq
//     // check how fraction of query_kmers are in cand_target_seq
//     return 0;
// }

int check_low_complexity(uint8_t *cand_query_seq, int cand_qlen, float low_comp_frac) {
    if (cand_qlen < 50) return 0;
    // check if the query sequence is low complexity
    // return 1 if low complexity, 0 otherwise
    uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
    int total_low_comp_len = 0;
    r = sdust(0, cand_query_seq, cand_qlen, T, W, &n);
    for (int i = 0; i < n; ++i) {
        total_low_comp_len += ((int)r[i] - (int)(r[i]>>32));
    }
    free(r);
    if (total_low_comp_len >= cand_qlen * low_comp_frac) return 1;
    else return 0;
}

// database: all k-mers, 0,1,2,...
// query: non-consecutive k-mers, 0,k,2k,...
// XXX: use overlapped k-mers for query?
int check_te_seq(const call_var_opt_t *opt, uint8_t *cand_te_seq, int cand_te_len, int *is_rev) {
    int total_count, for_count, rev_count, max_for_count = 0, max_rev_count = 0, max_for_i = -1, max_rev_i = -1;
    int min_count = 3; // or 4
    kmer32_v a = {0, 0, 0};
    collect_query_kmer(cand_te_seq, cand_te_len, 0, opt->te_kmer_len, &a);
    total_count = a.n;
    if (total_count <= 0) {
        free(a.a);
        return -1;
    }
    for (int i = 0; i < opt->n_te_seqs; ++i) {
        for_count = collect_kmer_hist(&a, opt->te_for_h[i]);
        rev_count = collect_kmer_hist(&a, opt->te_rev_h[i]);
        if (for_count > max_for_count) max_for_count = for_count, max_for_i = i;
        if (rev_count > max_rev_count) max_rev_count = rev_count, max_rev_i = i;
    }
    free(a.a);
    if (max_for_count > max_rev_count) {
        *is_rev = 0;
        if (max_for_count >= min_count) // total_count / 2)
            return max_for_i;
        else
            return -1;
    } else {
        *is_rev = 1;
        if (max_rev_count >= min_count) // total_count / 2)
            return max_rev_i;
        else
            return -1;
    }
}

int test_te_kmer_query(call_var_opt_t *opt, char *query_fn) {
    gzFile f = gzopen(query_fn, "r");
    kseq_t *ks = kseq_init(f);
    kseq_rewind(ks);
    int for_count, rev_count;
    while (kseq_read(ks) >= 0) {
        kmer32_v a = {0, 0, 0};
        collect_query_kmer((uint8_t*)ks->seq.s, ks->seq.l, 0, opt->te_kmer_len, &a);
        fprintf(stderr, "Query: %s, %ld\n", ks->name.s, a.n);
        for (int j = 0; j < opt->n_te_seqs; ++j) {
            for_count = collect_kmer_hist(&a, opt->te_for_h[j]);
            rev_count = collect_kmer_hist(&a, opt->te_rev_h[j]);
            fprintf(stderr, "TE: %s, for_count: %d, rev_count: %d\n", opt->te_seq_names[j], for_count, rev_count);
        }
        if (a.n > 0) free(a.a);
    }
    kseq_destroy(ks);
    gzclose(f);
    return 0;
}
