#include "seq.h"
#include "utils.h"
#include "htslib/kstring.h"
#include "htslib/faidx.h"

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

// 0,1,2,3,4,5=>A,C,G,T,N,-
// 65,97=>A, 67,99=>C, 71,103=>G, 84,85,116,117=>T
// else=>N
const char nt256_table[256] = {
       'A', 'C', 'G', 'T',  'N', '-', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', '-',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'A', 'N', 'C',  'N', 'N', 'N', 'G',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'T', 'T', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',
       'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N',  'N', 'N', 'N', 'N'
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

uint8_t get_bseq1(char *seq, hts_pos_t beg, hts_pos_t end, hts_pos_t pos) {
    if (pos < beg || pos > end) return 4;
    return nst_nt4_table[(int)seq[pos-beg]];
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

// // XXX load regional sequence if needed (similar to bam regions)
// ref_seq_t *read_ref_seq(const char *ref_fa_fn) {
//     if (!ref_fa_fn) _err_error_exit("No reference FASTA file.");
//     ref_seq_t *ref_seq = ref_seq_init(); gzFile ref_fp; kseq_t *ks; int absent;
//     ref_fp = xzopen(ref_fa_fn, "r");
//     ks = kseq_init(ref_fp); gzrewind(ref_fp); kseq_rewind(ks);
//     while (kseq_read(ks) >= 0) {
//         ref_seq_realloc(ref_seq);
//         kputsn(ks->seq.s, ks->seq.l, ref_seq->seq+ref_seq->n);
//         kputsn(ks->name.s, ks->name.l, ref_seq->name+ref_seq->n);
//         khint_t pos = kh_put(str, ref_seq->h, ref_seq->name[ref_seq->n].s, &absent);
//         if (absent) kh_val(ref_seq->h, pos) = ref_seq->n;
//         else _err_error_exit(__func__, "Duplicated chromosome: \"%s\".", ks->name.s);
//         ++ref_seq->n;
//     }
//     kseq_destroy(ks); err_gzclose(ref_fp);
//     return ref_seq;
// }

ref_seq_t *read_ref_seq(const char *ref_fa_fn) {
    ref_seq_t *ref_seq = ref_seq_init();
    htsFile *fp = hts_open(ref_fa_fn, "r");
    if (fp == NULL) _err_error_exit("Failed to open FASTA file from: %s\n", ref_fa_fn);
    kstring_t line = {0, 0, NULL}; int absent;
    while (hts_getline(fp, KS_SEP_LINE, &line) >= 0) {
        if (line.s[0] == '>') { // name
            ref_seq_realloc(ref_seq);
            ++ref_seq->n;
            int l = 0; for (l = 1; l < line.l; ++l) if (isspace(line.s[l])) break;
            kputsn(line.s+1, l-1, ref_seq->name+ref_seq->n-1);
            khint_t pos = kh_put(str, ref_seq->h, ref_seq->name[ref_seq->n-1].s, &absent);
            if (absent) kh_val(ref_seq->h, pos) = ref_seq->n-1;
            else _err_error_exit(__func__, "Duplicated chromosome: \"%s\".", ref_seq->name[ref_seq->n-1].s);
        } else { // seq
            kputsn(line.s, line.l, ref_seq->seq+ref_seq->n-1);
        }
    }
    free(line.s); hts_close(fp);
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
    if (!r || !r->h) _err_error_exit("No reference sequence.");
    khint_t pos = kh_get(str, r->h, name);
    return pos == kh_end(r->h) ? -1 : kh_val(r->h, pos);
}

ref_reg_seq_t *ref_reg_seq_init(void) {
    ref_reg_seq_t *r = (ref_reg_seq_t*)_err_malloc(sizeof(ref_reg_seq_t));
    r->n = r->m = 0;
    r->reg_seq = 0; r->name = 0;
    r->reg_cr = cr_init();
    return r;
}

ref_reg_seq_t *ref_reg_seq_realloc(ref_reg_seq_t *r) {
    if (r->n >= r->m) {
        r->m = r->m == 0 ? 1 : (r->m) << 1;
        r->reg_seq = (ref_reg_seq1_t*)_err_realloc(r->reg_seq, r->m * sizeof(ref_reg_seq1_t));
        r->name = (kstring_t*)_err_realloc(r->name, r->m * sizeof(kstring_t));
        int i;
        for (i = r->n; i < r->m; ++i) {
            r->reg_seq[i].seq = (kstring_t){0,0,0};
            r->name[i] = (kstring_t){0,0,0};
        }
    }
    return r;
}

void ref_reg_seq_free(ref_reg_seq_t *r) {
    if (r->m > 0) {
        int i;
        for (i = 0; i < r->m; ++i) {
            if (r->reg_seq[i].seq.m) free(r->reg_seq[i].seq.s);
            if (r->name[i].m) free(r->name[i].s);
        }
        free(r->reg_seq);
        free(r->name);
    }
    cr_destroy(r->reg_cr);
    free(r);
}

void read_ref_reg_seq1(faidx_t *fai, ref_reg_seq_t *r, const char *rname, hts_pos_t beg, hts_pos_t end) {
    if (beg >= end) _err_error_exit("Invalid region: %s:%d-%d\n", rname, beg, end);
    int len;
    int ref_seq_len = faidx_seq_len(fai, rname);
    int flank_len = 50000; // XXX use dynamimc flank length based on read length, needed for digar_from_ref_seq
    // hts_pos_t _beg = MAX_OF_TWO(0, beg-1000), _end = end+1000;
    int _beg = MAX_OF_TWO(flank_len, beg)-flank_len, _end = MIN_OF_TWO(ref_seq_len-flank_len, end)+flank_len;
    if (_beg >= _end || _beg >= ref_seq_len) _err_error_exit("Invalid region: %s:%d-%d\n", rname, beg, end);
    // fprintf(stderr, "rname: %s, beg: %" PRId64 ", end: %" PRId64 ", _beg: %d, _end: %d\n", rname, beg, end, _beg, _end);
    char *seq = faidx_fetch_seq(fai, rname, _beg, _end, &len);
    if (seq == NULL) _err_error_exit("Failed to fetch sequence (%s:%d-%d\n", rname, _beg, _end);
    ref_reg_seq_realloc(r);
    kputsn(seq, len, &r->reg_seq[r->n].seq); free(seq);
    kputsn(rname, strlen(rname), &r->name[r->n]);
    r->reg_seq[r->n].beg = _beg+1;
    r->reg_seq[r->n].end = _beg+len;
    // cr: (beg, end]
    // cr_add(r->reg_cr, rname, _beg, _end, r->n);
    cr_add(r->reg_cr, rname, beg, end > ref_seq_len?ref_seq_len : end, r->n);
    ++r->n;
}

ref_reg_seq_t *read_ref_reg_seq(const char *ref_fa_fn) {
    ref_reg_seq_t *r = ref_reg_seq_init();
    htsFile *fp = hts_open(ref_fa_fn, "r");
    if (fp == NULL) _err_error_exit("Failed to open FASTA file from: %s\n", ref_fa_fn);
    kstring_t line = {0, 0, NULL};
    while (hts_getline(fp, KS_SEP_LINE, &line) >= 0) {
        if (line.s[0] == '>') { // name
            ref_reg_seq_realloc(r);
            ++r->n;
            int l = 0; for (l = 1; l < line.l; ++l) if (isspace(line.s[l])) break;
            kputsn(line.s+1, l-1, r->name+r->n-1);
        } else { // seq
            kputsn(line.s, line.l, &(r->reg_seq[r->n-1].seq));
        }
    }
    for (int i = 0; i < r->n; ++i) {
        // cr: (beg, end]
        cr_add(r->reg_cr, r->name[i].s, 0, r->reg_seq[i].seq.l, i);
        r->reg_seq[i].beg = 1;
        r->reg_seq[i].end = r->reg_seq[i].seq.l;
    }
    cr_index(r->reg_cr);
    free(line.s); hts_close(fp);
    return r;
}

trans_ele_anno_t *trans_ele_anno_init(void) {
    trans_ele_anno_t *t = (trans_ele_anno_t*)_err_malloc(sizeof(trans_ele_anno_t));
    t->n = t->m = 0;
    t->seq = 0; t->dfam_ac = 0; t->rep_name = 0; t->rm_type = 0; t->rm_subtype = 0;
    t->seq_ver = 0; t->seq_lens = 0;
    t->h = kh_init(str);
    return t;
}

trans_ele_anno_t *trans_ele_anno_realloc(trans_ele_anno_t *t) {
    if (t->n >= t->m) {
        t->m = t->m == 0 ? 1 : (t->m) << 1;
        t->seq = (char**)_err_realloc(t->seq, t->m * sizeof(char*));
        t->seq_lens = (int*)_err_realloc(t->seq_lens, t->m * sizeof(int));
        t->dfam_ac = (char**)_err_realloc(t->dfam_ac, t->m * sizeof(char*));
        t->seq_ver = (char**)_err_realloc(t->seq_ver, t->m * sizeof(char*));
        t->rep_name = (char**)_err_realloc(t->rep_name, t->m * sizeof(char*));
        t->rm_type = (char**)_err_realloc(t->rm_type, t->m * sizeof(char*));
        t->rm_subtype = (char**)_err_realloc(t->rm_subtype, t->m * sizeof(char*));
        int i;
        for (i = t->n; i < t->m; ++i) {
            t->seq[i] = NULL;
            t->seq_lens[i] = 0;
            t->dfam_ac[i] = NULL;
            t->seq_ver[i] = NULL;
            t->rep_name[i] = NULL;
            t->rm_type[i] = NULL;
            t->rm_subtype[i] = NULL;
        }
    }
    return t;
}

trans_ele_anno_t *read_trans_ele_anno_from_embl(const char* fn) {
    gzFile fp = xzopen(fn, "r");
    if (fp == NULL) _err_error_exit("Failed to open Dfam file: %s\n", fn);
    trans_ele_anno_t *t = trans_ele_anno_init();
    char line[1024];
    int n_rep = 0, is_in_rep_mask_anno = 0;
    while (gzgets(fp, line, 1024) != NULL) {
        line[strcspn(line, "\r\n")] = '\0';
        if (strncmp("ID   ", line, 5) == 0) {
            t = trans_ele_anno_realloc(t);
            // extract the accession number
            char *semicolon = strchr(line+5, ';');
            if (semicolon == NULL) { _err_error_exit("Invalid Dfam file: %s %d\n", fn, __LINE__); }
            else *semicolon = '\0';
            t->dfam_ac[t->n] = strdup(line+5);
            // extract sequence version
            char *tmp = strchr(semicolon+5, ';');
            if (tmp == NULL) { _err_error_exit("Invalid Dfam file: %s %d\n", fn, __LINE__); }
            else *tmp = '\0';
            t->seq_ver[t->n] = strdup(semicolon+5);
        } else if (strncmp("NM   ", line, 5) == 0) {
            t->rep_name[t->n] = strdup(line+5);
        } else if (strncmp("CC   RepeatMasker Annotations", line, 29) == 0) {
            is_in_rep_mask_anno = 1;
        } else if (is_in_rep_mask_anno) {
            if (strncmp("CC        Type: ", line, 16) == 0) {
                t->rm_type[t->n] = strdup(line+16);
            } else if (strncmp("CC        SubType: ", line, 19) == 0) {
                t->rm_subtype[t->n] = strdup(line+19);
            } else if (strncmp("XX", line, 2) == 0) {
                is_in_rep_mask_anno = 0;
            }
        } else if (strncmp("SQ   Sequence ", line, 14) == 0) {
            int seq_len = atoi(line+14);
            t->seq_lens[t->n] = seq_len;
            t->seq[t->n] = (char*)_err_malloc(sizeof(char) * (seq_len+1));
            int i = 0;
            while (gzgets(fp, line, 1024) != NULL) {
                if (strncmp("//", line, 2) == 0) break;
                for (int j = 0; j < strlen(line); ++j) {
                    if (isalpha(line[j])) {
                        t->seq[t->n][i++] = line[j];
                    }
                }
            }
            t->seq[t->n][i] = '\0';
            ++t->n;
        }
    }
    // print all the sequences
    for (int i = 0; i < t->n; ++i) {
        fprintf(stderr, "%s\t%s\t%s\n", t->dfam_ac[i], t->rep_name[i], t->seq[i]);
    }
    err_gzclose(fp);
    return t;
}

void free_trans_ele_anno(trans_ele_anno_t *t) {
    if (t->m > 0) {
        int i;
        for (i = 0; i < t->m; ++i) {
            if (t->seq[i]) free(t->seq[i]);
            if (t->dfam_ac[i]) free(t->dfam_ac[i]);
            if (t->seq_ver[i]) free(t->seq_ver[i]);
            if (t->rep_name[i]) free(t->rep_name[i]);
            if (t->rm_type[i]) free(t->rm_type[i]);
            if (t->rm_subtype[i]) free(t->rm_subtype[i]);
        }
        free(t->seq);
        free(t->seq_lens);
        free(t->dfam_ac);
        free(t->seq_ver);
        free(t->rep_name);
        free(t->rm_type);
        free(t->rm_subtype);
    }
    kh_destroy(str, t->h);
    free(t);
}
