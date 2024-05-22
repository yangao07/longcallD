#include "seq.h"
#include "utils.h"
#include "htslib/kstring.h"

#define nt_A 0
#define nt_C 1
#define nt_G 2
#define nt_T 3
#define nt_N 4

// ACGTN=>01234
unsigned char nst_nt4_table[256] = {
	0, 1, 2, 3,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 5 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

// ACGTN=>32104
unsigned char com_nst_nt4_table[256] = {
	3, 2, 1, 0,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4 /*'-'*/, 4, 4,
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 3, 4, 2,  4, 4, 4, 1,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  0, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
	4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4
};

//replace 'N' with 'G':           A  C  G  T  N->G
const uint8_t hash_nt4_table[6] = {0, 1, 2, 3, 2, 2};

//                               T  G  C  A  N->G
const uint8_t rc_nt4_table[6] = {3, 2, 1, 0, 2, 2};

char n_char[6] = {'A', 'C', 'G', 'T', 'N' };

uint32_t hash_key(uint8_t *bseq, int seq_len) {
    int i; uint32_t hash_key = 0;
    for (i = 0; i < seq_len; ++i) {
        // if (bseq[i] >= nt_N) err_printf("Error in bseq.\n");
        hash_key = (hash_key << 2) | bseq[i];
    }
    return hash_key;
}

uint32_t hash_shift_key(uint32_t pre_key, uint8_t *bseq, int pre_i, int cur_i, int k) {
    int i; uint32_t hash_key = pre_key;
    for (i = pre_i+1; i <= cur_i; ++i) {
        // if (bseq[i+k] >= nt_N) err_printf("Error in bseq.\n");
        hash_key = (hash_key << 2) | bseq[i];
    }
    return hash_key & ((1 << 2*k) - 1);
}

uint8_t *get_bseq(char *seq, int seq_len) {
    int i;
    uint8_t *bseq = (uint8_t*)_err_malloc(seq_len * sizeof(uint8_t));
    for (i = 0; i < seq_len; ++i) {
        // N(ambiguous base)
        // bseq[i] = hash_nt4_table[nst_nt4_table[(int)seq[i]]];
        bseq[i] = nst_nt4_table[(int)seq[i]];
    }
    return bseq;
}

char *get_rc_seq(char *seq, int seq_len) {
    int i;
    char *rc_seq = (char*)_err_malloc(sizeof(char) * (seq_len+1));
    for (i = 0; i < seq_len; ++i) {
        rc_seq[seq_len-i-1] = "ACGTN"[com_nst_nt4_table[(int)seq[i]]];
    } rc_seq[seq_len] = '\0';
    return rc_seq;
}

ref_seq_t *ref_seq_init(void) {
    ref_seq_t *r = (ref_seq_t*)_err_malloc(sizeof(ref_seq_t));
    r->n = r->m = 0; r->seq = 0, r->name = 0;
    r->h = kh_init(str);
    return r;
}

ref_seq_t *ref_seq_realloc(ref_seq_t *r) {
    if (r->n >= r->m) {
        r->m = r->m == 0 ? 1 : (r->m) << 1;
        r->seq = (kstring_t*)_err_realloc(r->seq, r->m * sizeof(kstring_t));
        r->name = (kstring_t*)_err_realloc(r->name, r->m * sizeof(kstring_t));
        int i;
        for (i = r->n; i < r->m; ++i) {
            r->seq[i] = (kstring_t){0,0,0};
            r->name[i] = (kstring_t){0,0,0};
        }
    }
    return r;
}

ref_seq_t *read_ref_seq(const char *ref_fa_fn) {
    if (!ref_fa_fn) _err_fatal("No reference FASTA file.");
    ref_seq_t *ref_seq = ref_seq_init(); gzFile ref_fp; kseq_t *ks; int absent;
    ref_fp = xzopen(ref_fa_fn, "r");
    ks = kseq_init(ref_fp); gzrewind(ref_fp); kseq_rewind(ks);
    while (kseq_read(ks) >= 0) {
        ref_seq_realloc(ref_seq);
        kputsn(ks->seq.s, ks->seq.l, ref_seq->seq+ref_seq->n);
        kputsn(ks->name.s, ks->name.l, ref_seq->name+ref_seq->n);
        khint_t pos = kh_put(str, ref_seq->h, ref_seq->name[ref_seq->n].s, &absent);
        if (absent) kh_val(ref_seq->h, pos) = ref_seq->n;
        else err_fatal(__func__, "Duplicated chromosome: \"%s\".", ks->name.s);
        ++ref_seq->n;
    }
    kseq_destroy(ks); err_gzclose(ref_fp);
    return ref_seq;
}

void ref_seq_free(ref_seq_t *r) {
    if (r->m > 0) {
        int i;
        for (i = 0; i < r->m; ++i) {
            if (r->seq[i].m) free(r->seq[i].s);
            if (r->name[i].m) free(r->name[i].s);
        }
        free(r->seq); free(r->name);
    }
    kh_destroy(str, r->h);
    free(r);
}

int ref_seq_name2id(ref_seq_t *r, const char *name) {
    if (!r || !r->h) _err_fatal("No reference sequence.");
    khint_t pos = kh_get(str, r->h, name);
    return pos == kh_end(r->h) ? -1 : kh_val(r->h, pos);
}