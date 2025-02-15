#include "wavefront/wavefront_align.h"
#include "alignment/cigar.h"
#include "seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"
#include "edlib.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "utils.h"
#include "align.h"

extern int LONGCALLD_VERBOSE;

int wfa_collect_pretty_alignment(cigar_t* const cigar,
                                 const uint8_t* const pattern,
                                 const int pattern_length,
                                 const uint8_t* const text,
                                 const int text_length,
                                 uint8_t **_pattern_alg, uint8_t **_text_alg) {
    // Parameters
    char* const operations = cigar->operations;
    const int begin_offset = cigar->begin_offset;
    const int end_offset = cigar->end_offset;
    // Allocate alignment buffers
    const int max_buffer_length = text_length + pattern_length + 1;
    uint8_t *mem = (uint8_t*)calloc(2*max_buffer_length,1);
    uint8_t *pattern_alg = mem;
    uint8_t *text_alg = mem + max_buffer_length;
    // Compute alignment buffers
    int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
    for (i=begin_offset;i<end_offset;++i) {
        switch (operations[i]) {
            case 'M':
                if (pattern[pattern_pos] != text[text_pos]) {
                    pattern_alg[alg_pos] = pattern[pattern_pos];
                    text_alg[alg_pos++] = text[text_pos];
                } else {
                    pattern_alg[alg_pos] = pattern[pattern_pos];
                    text_alg[alg_pos++] = text[text_pos];
                }
                pattern_pos++; text_pos++;
                break;
            case 'X':
                if (pattern[pattern_pos] != text[text_pos]) {
                    pattern_alg[alg_pos] = pattern[pattern_pos++];
                    text_alg[alg_pos++] = text[text_pos++];
                } else {
                    pattern_alg[alg_pos] = pattern[pattern_pos++];
                    text_alg[alg_pos++] = text[text_pos++];
                }
                break;
            case 'I':
                pattern_alg[alg_pos] = 5; //'-';
                text_alg[alg_pos++] = text[text_pos++];
                break;
            case 'D':
                pattern_alg[alg_pos] = pattern[pattern_pos++];
                text_alg[alg_pos++] = 5; //'-';
                break;
            default: break;
        }
    }
    assert(pattern_pos == pattern_length && text_pos == text_length);
    *_pattern_alg = pattern_alg; *_text_alg = text_alg;
    return alg_pos;
}

// return alignment score
int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen, int *n_eq, int *n_xid) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = 4;
    attributes.affine2p_penalties.gap_opening1 = 4;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 24;
    attributes.affine2p_penalties.gap_extension2 = 1;
    attributes.alignment_scope = compute_score;
    attributes.alignment_form.span = alignment_end2end;
    // xdrops
    attributes.heuristic.strategy = wf_heuristic_xdrop;
    attributes.heuristic.xdrop = MAX_OF_TWO(100, (int)((plen + tlen) * 0.1)); 
    attributes.heuristic.steps_between_cutoffs = 100;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    wavefront_align(wf_aligner, (const char*)pattern, plen, (const char*)text, tlen);
    // collect score
    int score = wf_aligner->cigar->score; //cigar_score_gap_affine2p(wf_aligner->cigar, &attributes.affine2p_penalties);
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        uint32_t *cigar; int cigar_length;
        // collect cigar
        cigar_get_CIGAR(wf_aligner->cigar, true, &cigar, &cigar_length);
        *n_eq = 0, *n_xid = 0;
        for (int i = 0; i < cigar_length; ++i) {
            int op = cigar[i] & 0xf; int len = cigar[i] >> 4;
            if (op == BAM_CEQUAL) *n_eq += len;
            else if (op == BAM_CDIFF || op == BAM_CINS || op == BAM_CDEL) *n_xid += len;
        }
    }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    return score;
}

// XXX filter out if the alignment score is too low
int wfa_ext_aln(int ext_direction, uint8_t *pattern, int plen, uint8_t *text, int tlen, uint32_t **cigar_buf) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = 4;
    attributes.affine2p_penalties.gap_opening1 = 4;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 24;
    attributes.affine2p_penalties.gap_extension2 = 1;
    attributes.alignment_scope = compute_alignment;
    attributes.heuristic.strategy = wf_heuristic_none;
    attributes.alignment_form.span = alignment_endsfree;
    if (ext_direction == LONGCALLD_EXT_ALN_RIGHT) { // left -> right
        attributes.alignment_form.pattern_begin_free = 0;
        attributes.alignment_form.text_begin_free = 0;
        if (plen > tlen) {
            attributes.alignment_form.pattern_end_free = plen-tlen;
            attributes.alignment_form.text_end_free = 0;
        } else {
            attributes.alignment_form.pattern_end_free = 0;
            attributes.alignment_form.text_end_free = tlen-plen;
        }
    } else { // right -> left
        attributes.alignment_form.pattern_end_free = 0;
        attributes.alignment_form.text_end_free = 0;
        if (plen > tlen) {
            attributes.alignment_form.pattern_begin_free = plen-tlen;
            attributes.alignment_form.text_begin_free = 0;
        } else {
            attributes.alignment_form.pattern_begin_free = 0;
            attributes.alignment_form.text_begin_free = tlen-plen;
        }
    }
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    wavefront_align(wf_aligner, (const char*)pattern, plen, (const char*)text, tlen);
    cigar_t *cigar = wf_aligner->cigar;
    if (LONGCALLD_VERBOSE >= 2) {
        char *p = (char*)malloc(plen+1); char *t = (char*)malloc(tlen+1);
        for (int i = 0; i < plen; ++i) p[i] = "ACGTN"[pattern[i]];
        for (int i = 0; i < tlen; ++i) t[i] = "ACGTN"[text[i]];
        cigar_print_pretty(stderr, cigar, p, plen, t, tlen);
        free(p); free(t);
    }
    // collect cigar
    uint32_t *tmp_cigar; int cigar_length;
    cigar_get_CIGAR(wf_aligner->cigar, true, &tmp_cigar, &cigar_length);
    *cigar_buf = (uint32_t*)malloc(cigar_length * sizeof(uint32_t));

    for (int i = 0; i < cigar_length; ++i) {
        (*cigar_buf)[i] = tmp_cigar[i];
    }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    return cigar_length;
}

int wfa_end2end_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen,
                    int gap_aln, int a, int b, int q, int e, int q2, int e2, uint32_t **cigar_buf, int *cigar_length,
                    uint8_t **pattern_alg, uint8_t **text_alg, int *alg_length) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0; // -a;
    attributes.affine2p_penalties.mismatch = b;
    attributes.affine2p_penalties.gap_opening1 = q;
    attributes.affine2p_penalties.gap_extension1 = e;
    attributes.affine2p_penalties.gap_opening2 = q2;
    attributes.affine2p_penalties.gap_extension2 = e2;
    attributes.alignment_scope = compute_alignment;
    attributes.alignment_form.span = alignment_end2end;
    attributes.heuristic.strategy = wf_heuristic_none;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    uint8_t *p = pattern, *t = text;
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse pattern and text
        p = (uint8_t*)malloc(plen); t = (uint8_t*)malloc(tlen);
        for (int i = 0; i < plen; ++i) p[i] = pattern[plen-i-1];
        for (int i = 0; i < tlen; ++i) t[i] = text[tlen-i-1];
    }

    // Align
    wavefront_align(wf_aligner, (const char*)p, plen, (const char*)t, tlen);
    cigar_t *cigar = wf_aligner->cigar;
    // if (LONGCALLD_VERBOSE >= 2) {
    //     char *p_seq = (char*)malloc(plen+1); char *t_seq = (char*)malloc(tlen+1);
    //     for (int i = 0; i < plen; ++i) p_seq[i] = "ACGTN"[p[i]];
    //     for (int i = 0; i < tlen; ++i) t_seq[i] = "ACGTN"[t[i]];
    //     cigar_print_pretty(stderr, cigar, p_seq, plen, t_seq, tlen);
    //     free(p_seq); free(t_seq);
    // }
    // collect cigar
    if (cigar_buf != NULL && cigar_length != NULL) {
        uint32_t *tmp_cigar;
        cigar_get_CIGAR(wf_aligner->cigar,true, &tmp_cigar, cigar_length);
        *cigar_buf = (uint32_t*)malloc(*cigar_length * sizeof(uint32_t));
        if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse cigar
            for (int i = 0; i < *cigar_length; ++i) {
                (*cigar_buf)[i] = tmp_cigar[*cigar_length-i-1];
            }
        } else {
            for (int i = 0; i < *cigar_length; ++i) {
                (*cigar_buf)[i] = tmp_cigar[i];
            }
        }
    }
    // collect alignment string
    if (pattern_alg != NULL && text_alg != NULL) {
        *alg_length = wfa_collect_pretty_alignment(cigar, p, plen, t, tlen, pattern_alg, text_alg);
        if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse alignment
            for (int i = 0; i < *alg_length/2; ++i) {
                uint8_t tmp = (*pattern_alg)[i];
                (*pattern_alg)[i] = (*pattern_alg)[*alg_length-i-1];
                (*pattern_alg)[*alg_length-i-1] = tmp;
                tmp = (*text_alg)[i];
                (*text_alg)[i] = (*text_alg)[*alg_length-i-1];
                (*text_alg)[*alg_length-i-1] = tmp;
            }
        }
    }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { free(p); free(t); }
    return 0;
}

int wfa_collect_aln_str(const call_var_opt_t *opt, uint8_t *target, int tlen, uint8_t *query, int qlen, aln_str_t *aln_str) {
    aln_str->target_aln = 0; aln_str->query_aln = 0; aln_str->aln_len = 0;
    wfa_end2end_aln(target, tlen, query, qlen,
                    opt->gap_aln, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2,
                    NULL, NULL, &aln_str->target_aln, &aln_str->query_aln, &aln_str->aln_len);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, ">target %d\n", tlen);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->target_aln[i]]);
        fprintf(stderr, "\n>query %d\n", qlen);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->query_aln[i]]);
        fprintf(stderr, "\n");
    }
    return 0;
}

int end2end_aln(const call_var_opt_t *opt, char *tseq, int tlen, uint8_t *qseq, int qlen, uint32_t **cigar_buf) {
    if (qlen <= 0 || tlen <= 0) return 0;
    int min_len = MIN_OF_TWO(tlen, qlen), max_len = MAX_OF_TWO(tlen, qlen);
    int delta_len = MAX_OF_TWO(1, max_len - min_len);

    uint8_t *tseq2 = (uint8_t*)malloc(max_len);
    for (int i = 0; i < tlen; ++i) tseq2[i] = nst_nt4_table[(uint8_t)tseq[i]];
    // use wfa if the length difference is small
    // if (max_len * delta_len + delta_len * delta_len < max_len * min_len) {
    int cigar_len = 0;
    wfa_end2end_aln(tseq2, tlen, qseq, qlen,
                    opt->gap_aln, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2,
                    cigar_buf, &cigar_len, NULL, NULL, NULL);
    // } else { // if (max_len - min_len > 1000) { // use ksw2 if the length difference is large
        // cigar_len = ksw2_aln(opt->gap_aln, tseq2, tlen, qseq, qlen, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, cigar_buf);
    // }
    free(tseq2);
    return cigar_len;
}

void collect_aln_beg_end(uint32_t *cigar_buf, int cigar_len, int ext_direction, int ref_len, int *ref_beg, int *ref_end, int read_len, int *read_beg, int *read_end) {
    *ref_beg = 1; *read_beg = 1; *ref_end = ref_len, *read_end = read_len;
    if (ext_direction == LONGCALLD_EXT_ALN_RIGHT) {
        int tmp_ref_end = 0, tmp_read_end = 0;
        for (int i = 0; i < cigar_len; ++i) { // collect ref_end & read_end with last M cigar
            int op = cigar_buf[i] & 0xf; int len = cigar_buf[i] >> 4;
            if (op == BAM_CEQUAL || op == BAM_CMATCH) {
                tmp_ref_end += len; tmp_read_end += len;
                *ref_end = tmp_ref_end; *read_end = tmp_read_end;
            } else if (op == BAM_CDIFF) {
                tmp_ref_end += len; tmp_read_end += len;
            } else if (op == BAM_CDEL) {
                tmp_ref_end += len;
            } else if (op == BAM_CINS) {
                tmp_read_end += len;
            } else continue;
        }
    } else {
        int tmp_ref_beg = 1, tmp_read_beg = 1;
        for (int i = 0; i < cigar_len; ++i) {
            int op = cigar_buf[i] & 0xf; int len = cigar_buf[i] >> 4;
            if (op == BAM_CEQUAL || op == BAM_CMATCH) {
                *ref_beg = tmp_ref_beg; *read_beg = tmp_read_beg;
                break;
                // tmp_ref_beg += len; tmp_read_beg += len;
            } else if (op == BAM_CDIFF) {
                tmp_ref_beg += len; tmp_read_beg += len;
            } else if (op == BAM_CDEL) {
                tmp_ref_beg += len;
            } else if (op == BAM_CINS) {
                tmp_read_beg += len;
            } else continue;
        }
    }
}

int cal_wfa_partial_aln_beg_end(int ext_direction, int ref_len, uint8_t *ref_seq, int read_len, uint8_t *read_seq, int *ref_beg, int *ref_end, int *read_beg, int *read_end) {
    uint32_t *cigar_buf = NULL;
    int ret = 1;
    int cigar_len = wfa_ext_aln(ext_direction, ref_seq, ref_len, read_seq, read_len, &cigar_buf);
    if (cigar_len == 0) ret = 0;
    else collect_aln_beg_end(cigar_buf, cigar_len, ext_direction, ref_len, ref_beg, ref_end, read_len, read_beg, read_end);
    if (cigar_buf != NULL) free(cigar_buf);
    return ret;
}

int collect_partial_aln_beg_end(int ref_len, uint8_t *ref_seq, int ref_full_cover, int read_len, uint8_t *read_seq, int read_full_cover, 
                                 int *ref_beg, int *ref_end, int *read_beg, int *read_end) {
    *ref_beg = 1, *ref_end = ref_len, *read_beg = 1, *read_end = read_len;
    int ret = 1;
    if (ref_full_cover == 3) { // ref is full cover
        if (read_full_cover == 3) {
            return 1;
        } else {
            if (read_full_cover == 1) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT, ref_len, ref_seq, read_len, read_seq, ref_beg, ref_end, read_beg, read_end);
            } else if (read_full_cover == 2) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT, ref_len, ref_seq, read_len, read_seq, ref_beg, ref_end, read_beg, read_end);
            }
        }
    } else if (ref_full_cover == 1) { // ref is left-cover
        if (read_full_cover != 1) _err_error_exit("Ref is left-cover but read is not left-cover\n");
        // assert(read_full_cover == 1);
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT, ref_len, ref_seq, read_len, read_seq, ref_beg, ref_end, read_beg, read_end);
    } else if (ref_full_cover == 2) { // ref is right-cover
        if (read_full_cover != 2) _err_error_exit("Ref is right-cover but read is not right-cover\n");
        // assert(read_full_cover == 2);
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT, ref_len, ref_seq, read_len, read_seq, ref_beg, ref_end, read_beg, read_end);
    } else _err_error_exit("Ref is not left or right-cover\n");
    return ret;
}

int abpoa_partial_aln_msa_cons(const call_var_opt_t *opt, abpoa_t *ab, int wb, int n_reads, int *read_ids, uint8_t **read_seqs, int *read_lens, int *read_full_cover, char **names,
                               int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids, int *msa_seq_lens, uint8_t **msa_seqs) {
    // abpoa_t *ab = abpoa_init();
    int needs_free_ab = 0;
    if (ab == NULL) {
        ab = abpoa_init();
        needs_free_ab = 1;
    }
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = wb;
    abpt->cons_algrm = ABPOA_MF;
    // abpt->sub_aln = 1;
    abpt->inc_path_score = 1;
    if (cons_lens != NULL && cons_seqs != NULL) abpt->out_cons = 1; else abpt->out_cons = 0;
    if (msa_seq_lens != NULL && msa_seqs != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    abpt->max_n_cons = max_n_cons;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;

    // msa
    if (msa_seq_lens != NULL && msa_seqs != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    if (cons_lens != NULL && cons_seqs != NULL) abpt->out_cons = 1; else abpt->out_cons = 0;
    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpoa_post_set_para(abpt);
    ab->abs->n_seq = n_reads;
    for (int i = 0; i < n_reads; ++i) {
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, ">%s %d %d\n", names[i], read_lens[i], read_full_cover[i]);
            for (int j = 0; j < read_lens[i]; ++j) {
                fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
            } fprintf(stderr, "\n");
        }
        abpoa_res_t res;
        res.graph_cigar = 0, res.n_cigar = 0;
        int exc_beg = 0, exc_end = 1, seq_beg_cut = 0, seq_end_cut = 0;
        if (i != 0) {
            int ref_beg, ref_end, read_beg, read_end, beg_id, end_id;
            if (collect_partial_aln_beg_end(read_lens[0], read_seqs[0], read_full_cover[0], read_lens[i], read_seqs[i], read_full_cover[i], &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Skipped in POA: %s\n", names[i]);
                continue;
            }
            beg_id = ref_beg+1, end_id = ref_end+1;
            seq_beg_cut = read_beg - 1, seq_end_cut = read_lens[i] - read_end;
            abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, &exc_beg, &exc_end);
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ExcBeg: %d, ExcEnd: %d, SeqBegCut: %d, SeqEndCut: %d\n", exc_beg, exc_end, seq_beg_cut, seq_end_cut);
        abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, read_lens[i]-seq_beg_cut-seq_end_cut, &res);
        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, NULL, read_lens[i]-seq_beg_cut-seq_end_cut, NULL, res, i, n_reads, 1);
        if (res.n_cigar) free(res.graph_cigar);
    }
    if (LONGCALLD_VERBOSE >= 2) abpoa_output(ab, abpt, stderr);
    else abpoa_output(ab, abpt, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    // cons bases
    if (cons_lens != NULL && cons_seqs != NULL) {
        if (abc->n_cons > 0) {
            for (int i = 0; i < abc->n_cons; ++i) {
                cons_lens[i] = abc->cons_len[i];
                cons_seqs[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
                for (int j = 0; j < abc->cons_len[i]; ++j) {
                    cons_seqs[i][j] = abc->cons_base[i][j];
                }
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
            }
            if (clu_n_seqs != NULL && clu_read_ids != NULL) {
                if (abc->n_cons > 1) {
                    for (int i = 0; i < abc->n_cons; ++i) {
                        clu_n_seqs[i] = abc->clu_n_seq[i];
                        clu_read_ids[i] = (int*)malloc(abc->clu_n_seq[i] * sizeof(int));
                        for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                            clu_read_ids[i][j] = read_ids[abc->clu_read_ids[i][j]];
                        }
                    }
                } else {
                    *clu_n_seqs = n_reads;
                    *clu_read_ids = (int*)malloc(n_reads * sizeof(int));
                    for (int i = 0; i < n_reads; ++i) (*clu_read_ids)[i] = read_ids[i];
                }
            }
            n_cons = abc->n_cons;
        } else {
            fprintf(stderr, "Unable to call consensus: %s\n", names[0]);
        }
    }

    // msa bases, include consensus sequences
    if (msa_seq_lens != NULL && msa_seqs != NULL) {
        *msa_seq_lens = abc->msa_len;
        for (int i = 0; i < abc->n_seq; ++i) {
            msa_seqs[i] = (uint8_t*)malloc(abc->msa_len * sizeof(uint8_t));
            for (int j = 0; j < abc->msa_len; ++j) {
                msa_seqs[i][j] = abc->msa_base[i][j];
            }
        }
    }
    if (needs_free_ab) abpoa_free(ab); 
    abpoa_free_para(abpt); 
    return n_cons;
 }

// reads are sorted by full-cover and length before calling abPOA
// use reads with full_cover == 3 as backbone for poa
// skip reads with full_cover == 0
int abpoa_partial_aln(const call_var_opt_t *opt, int n_reads, uint8_t **read_seqs, int *read_lens, int *read_full_cover, char **names, int max_n_cons, int *cons_lens, uint8_t **cons_seqs) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    // abpt->wb = -1;
    abpt->out_msa = 0;
    // abp->cons_algrm = ABPOA_MF;
    abpt->inc_path_score = 1;
    abpt->max_n_cons = max_n_cons;
    // msa
    abpt->out_cons = 1;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;
    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpoa_post_set_para(abpt);
    ab->abs->n_seq = n_reads;
    for (int i = 0; i < n_reads; ++i) {
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, ">%s %d %d\n", names[i], read_lens[i], read_full_cover[i]);
            for (int j = 0; j < read_lens[i]; ++j) {
                fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
            } fprintf(stderr, "\n");
        }
        abpoa_res_t res;
        res.graph_cigar = 0, res.n_cigar = 0;
        int exc_beg = 0, exc_end = 1, seq_beg_cut = 0, seq_end_cut = 0;
        if (i != 0) {
            int ref_beg, ref_end, read_beg, read_end, beg_id, end_id;
            if (collect_partial_aln_beg_end(read_lens[0], read_seqs[0], read_full_cover[0], read_lens[i], read_seqs[i], read_full_cover[i], &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Skipped in POA: %s\n", names[i]);
                continue;
            }
            beg_id = ref_beg+1, end_id = ref_end+1;
            seq_beg_cut = read_beg - 1, seq_end_cut = read_lens[i] - read_end;
            abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, &exc_beg, &exc_end);
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ExcBeg: %d, ExcEnd: %d, SeqBegCut: %d, SeqEndCut: %d\n", exc_beg, exc_end, seq_beg_cut, seq_end_cut);
        abpoa_align_sequence_to_subgraph(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, read_lens[i]-seq_beg_cut-seq_end_cut, &res);
        abpoa_add_subgraph_alignment(ab, abpt, exc_beg, exc_end, read_seqs[i]+seq_beg_cut, NULL, read_lens[i]-seq_beg_cut-seq_end_cut, NULL, res, i, n_reads, 0);
        if (res.n_cigar) free(res.graph_cigar);
    }
    if (LONGCALLD_VERBOSE >= 2) abpoa_output(ab, abpt, stderr);
    else abpoa_output(ab, abpt, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    if (abc->n_cons > 0) {
        for (int i = 0; i < abc->n_cons; ++i) {
            cons_lens[i] = abc->cons_len[i];
            cons_seqs[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
            for (int j = 0; j < abc->cons_len[i]; ++j) {
                cons_seqs[i][j] = abc->cons_base[i][j];
            }
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
        }
        n_cons = abc->n_cons;
    } else {
        fprintf(stderr, "Unable to call consensus: %s\n", names[0]);
    }
    abpoa_free(ab); abpoa_free_para(abpt);
    return n_cons;
 }

 // XXX limit abpoa memory usage, avoid memory allocation failure
 int abpoa_aln_msa_cons(const call_var_opt_t *opt, int wb, int n_reads, int *read_ids, uint8_t **read_seqs, int *read_lens, int max_n_cons,
                        int *cons_lens, uint8_t **cons_seqs,
                        int *clu_n_seqs, int **clu_read_ids, int *msa_seq_len, uint8_t ***msa_seq) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = wb;
    abpt->inc_path_score = 1;
    abpt->out_msa = 0; abpt->out_cons = 1;
    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->max_n_cons = max_n_cons;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;
    abpoa_post_set_para(abpt);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA (max %d cons): %d\n", max_n_cons, n_reads);
        abpoa_msa(ab, abpt, n_reads, NULL, read_lens, read_seqs, NULL, stderr);
    } else abpoa_msa(ab, abpt, n_reads, NULL, read_lens, read_seqs, NULL, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    // cons bases
    if (abc->n_cons > 0) {
        for (int i = 0; i < abc->n_cons; ++i) {
            cons_lens[i] = abc->cons_len[i];
            cons_seqs[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
            for (int j = 0; j < abc->cons_len[i]; ++j) {
                cons_seqs[i][j] = abc->cons_base[i][j];
            }
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
        }
        if (abc->n_cons == 2) { // collect read ids for each cluster
            for (int i = 0; i < 2; ++i) {
                clu_n_seqs[i] = abc->clu_n_seq[i];
                clu_read_ids[i] = (int*)malloc(abc->clu_n_seq[i] * sizeof(int));
                for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                    clu_read_ids[i][j] = read_ids[abc->clu_read_ids[i][j]];
                }
            }
        } else {
            *clu_n_seqs = n_reads;
            *clu_read_ids = (int*)malloc(n_reads * sizeof(int));
            for (int i = 0; i < n_reads; ++i) (*clu_read_ids)[i] = read_ids[i];
        }
        n_cons = abc->n_cons;
    } else {
        _err_error_exit("No Consensus: %ld\n");
    }
    // msa bases
    if (msa_seq != NULL && msa_seq_len != NULL) {
        // split msa and cons into 2 clusers
        for (int i = 0; i < n_cons; ++i) {
            msa_seq_len[i] = abc->msa_len;
            for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                msa_seq[i][j] = (uint8_t*)malloc(abc->msa_len * sizeof(uint8_t));
                for (int k = 0; k < abc->msa_len; ++k)
                    msa_seq[i][j][k] = abc->msa_base[abc->clu_read_ids[i][j]][k];
            }
            msa_seq[i][abc->clu_n_seq[i]] = (uint8_t*)malloc(abc->msa_len * sizeof(uint8_t));
            for (int j = 0; j < abc->msa_len; ++j)
                msa_seq[i][abc->clu_n_seq[i]][j] = abc->msa_base[abc->n_seq+i][j];
        }
    }
    abpoa_free_para(abpt); abpoa_free(ab); 
    return n_cons;
}

int collect_reg_ref_cseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, char **ref_cseq) {
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    *ref_cseq = (char*)malloc((reg_end - reg_beg + 1) * sizeof(char));
    for (int i = reg_beg; i <= reg_end; ++i) { // collect upper-case ref seq
        (*ref_cseq)[i-reg_beg] = toupper(ref_seq[i-ref_beg]);
    }
    return (reg_end - reg_beg + 1);
}

int collect_reg_ref_bseq(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, uint8_t **ref_bseq) {
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    *ref_bseq = (uint8_t*)malloc((reg_end - reg_beg + 1) * sizeof(uint8_t));
    for (int i = reg_beg; i <= reg_end; ++i) {
        (*ref_bseq)[i-reg_beg] = nst_nt4_table[(int)ref_seq[i-ref_beg]];
    }
    return (reg_end - reg_beg + 1);
}

void sort_by_full_cover_and_length(int n_reads, int *read_ids, int *read_lens, uint8_t **read_seqs, int *full_covers, char **names, int *haps, hts_pos_t *phase_sets) {
    for (int i = 0; i < n_reads-1; ++i) {
        for (int j = i+1; j < n_reads; ++j) {
            if (full_covers[i] < full_covers[j]) {
                int tmp_id = read_ids[i]; read_ids[i] = read_ids[j]; read_ids[j] = tmp_id;
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                int tmp_cover = full_covers[i]; full_covers[i] = full_covers[j]; full_covers[j] = tmp_cover;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
                int tmp_hap = haps[i]; haps[i] = haps[j]; haps[j] = tmp_hap;
                hts_pos_t tmp_ps = phase_sets[i]; phase_sets[i] = phase_sets[j]; phase_sets[j] = tmp_ps;
            } else if (full_covers[i] == full_covers[j] && read_lens[i] < read_lens[j]) {
                int tmp_id = read_ids[i]; read_ids[i] = read_ids[j]; read_ids[j] = tmp_id;
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
                int tmp_hap = haps[i]; haps[i] = haps[j]; haps[j] = tmp_hap;
                hts_pos_t tmp_ps = phase_sets[i]; phase_sets[i] = phase_sets[j]; phase_sets[j] = tmp_ps;
            }
        }
    }
}

int add_phase_set(hts_pos_t ps, hts_pos_t *uniq_phase_sets, int *n_uniq_phase_sets) {
    int i;
    for (i = 0; i < *n_uniq_phase_sets; ++i) {
        if (uniq_phase_sets[i] == ps) return i;
    }
    uniq_phase_sets[i] = ps;
    (*n_uniq_phase_sets)++;
    return i;
}

int determine_hap_based_on_noisy_cons(int *cons_lens, uint8_t **cons_seqs, int n_cons, int read_len, uint8_t *read_seq, int full_cover) {
    if (n_cons == 0) return 0;
    int scores[2] = {0, 0};
    for (int i = 0; i < n_cons; ++i) { // 0'th cons -> HAP1, 1'th cons -> HAP2
        int score = INT32_MIN;
        if (full_cover == 3) {
            score = wfa_heuristic_aln(cons_seqs[i], cons_lens[i], read_seq, read_len, NULL, NULL);
        } else if (full_cover == 1) {
            int len = MIN_OF_TWO(cons_lens[i], read_len);
            score = wfa_heuristic_aln(cons_seqs[i], len, read_seq, len, NULL, NULL);
        } else if (full_cover == 2) {
            int len = MIN_OF_TWO(cons_lens[i], read_len);
            score = wfa_heuristic_aln(cons_seqs[i]+cons_lens[i]-len, len, read_seq+read_len-len, len, NULL, NULL);
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "%d: %d\t", i+1, score);
            for (int j = 0; j < cons_lens[i]; ++j) {
                fprintf(stderr, "%c", "ACGTN"[cons_seqs[i][j]]);
            } fprintf(stderr, "\t");
            for (int j = 0; j < read_len; ++j) {
                fprintf(stderr, "%c", "ACGTN"[read_seq[j]]);
            } fprintf(stderr, "\n");
        }
        if (score != INT32_MIN) scores[i] = score;
    }
    if (scores[0] != scores[1]) {
        return (scores[0] > scores[1] ? 1 : 2);
    } else return 0;
}

int wfa_collect_noisy_aln_str_no_ps_hap(const call_var_opt_t *opt, int n_reads, int *read_ids, int *lens, uint8_t **seqs, char **qnames, int *fully_covers,
                                        uint8_t *ref_seq, int ref_seq_len, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs) {
    int *full_read_ids = (int*)malloc((n_reads+2) * sizeof(int));
    int *full_read_lens = (int*)malloc((n_reads+2) * sizeof(int));
    uint8_t **full_read_seqs = (uint8_t**)malloc((n_reads+2) * sizeof(uint8_t*));
    char **full_read_names = (char**)malloc((n_reads+2) * sizeof(char*));
    int *full_fully_covers = (int*)malloc((n_reads+2) * sizeof(int));
    int n_full_reads = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] <= 0 || fully_covers[i] != 3) continue;
        full_fully_covers[n_full_reads] = fully_covers[i];
        full_read_ids[n_full_reads] = i;
        full_read_lens[n_full_reads] = lens[i];
        full_read_seqs[n_full_reads] = seqs[i];
        full_read_names[n_full_reads] = qnames[i];
        n_full_reads++;
    }
    // abpoa
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    for (int i = 0; i < 2; ++i) cons_seqs[i] = NULL;
    int n_cons = 0;
    if (n_full_reads == 0) goto collect_noisy_msa_cons_no_ps_hap_end;
    n_cons = abpoa_aln_msa_cons(opt, -1, n_full_reads, full_read_ids, full_read_seqs, full_read_lens, 2,
                                cons_lens, cons_seqs, clu_n_seqs, clu_read_ids, NULL, NULL);

    // re-do POA with ref_seq and cons
    for (int i = 0; i < n_cons; ++i) {
        aln_str_t *clu_aln_str = aln_strs[i];
        wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[i], cons_lens[i], LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
        n_full_reads = 0;
        // full_read_ids[n_full_reads] = -1; full_fully_covers[n_full_reads] = 3; full_read_lens[n_full_reads] = ref_seq_len; full_read_seqs[n_full_reads] = ref_seq; full_read_names[n_full_reads++] = "REF";
        // full_read_ids[n_full_reads] = -2; full_fully_covers[n_full_reads] = 3; full_read_lens[n_full_reads] = cons_lens[i]; full_read_seqs[n_full_reads] = cons_seqs[i]; full_read_names[n_full_reads++] = "CONS";
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_i = clu_read_ids[i][j];
            int read_id = read_ids[read_i];
            clu_read_ids[i][j] = read_id;
            // cons vs read
            wfa_collect_aln_str(opt, cons_seqs[i], cons_lens[i], seqs[read_i], lens[read_i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads));
            // ref vs read
            wfa_collect_aln_str(opt, ref_seq, ref_seq_len, seqs[read_i], lens[read_i], LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_full_reads));
            // full_read_ids[n_full_reads] = read_id;
            // full_read_lens[n_full_reads] = lens[read_i];
            // full_read_seqs[n_full_reads] = seqs[read_i];
            // full_read_names[n_full_reads] = qnames[read_i];
            // full_fully_covers[n_full_reads] = fully_covers[read_i];
            n_full_reads++;
        }
    } 
    // XXX add clipped / partial aligned reads to the MSA
    // if (n_cons == 2) {
        // XXX update PS & hap if n_cons == 2
        // int n_old_haps[2] = {0, 0};
        // for (int hap = 1; hap <= 2; ++hap) {
        //     for (int i = 0; i < n_clu_reads[hap-1]; ++i) {
        //         int read_i = full_read_to_idx[clu_read_ids[hap-1][i]];
        //         phase_sets[read_i] = ps; haps[read_i] = hap;
        //         n_old_haps[hap-1] += 1;
        //         // fprintf(stderr, "UpdatePSHap: %s %d %d %d\n", qnames[read_i], hap, ps, clu_read_ids[hap-1][i]);
        //     }
        // }
        // // XXX for other reads, align to the consensus and then determine the haplotype
        // int assign_hap, n_new_haps[2]={0, 0};
        // for (int i = 0; i < n_reads; ++i) {
        //     if (lens[i] > 0 && (fully_covers[i] == 1 || fully_covers[i] == 2)) {
        //         assign_hap = determine_hap_based_on_noisy_cons(cons_lens, cons_seqs, n_cons, lens[i], seqs[i], fully_covers[i]);
        //         if (assign_hap > 0) {
        //             if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "AssignHap: %s %d\n", qnames[i], assign_hap);
        //             haps[i] = assign_hap; phase_sets[i] = ps;
        //             n_new_haps[assign_hap-1] += 1;
        //         }
        //     }
        // }
        // // re-do consensus calling using all reads
        // if (opt->disable_read_realign == 0) {
        //     for (int hap=1; hap<=2; ++hap) {
        //         if (n_new_haps[hap-1] >= n_old_haps[hap-1] / 4) {
        //             if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Re-Do Consensus Calling: %d %d\n", n_old_haps[hap-1], n_new_haps[hap-1]);
        //             free(cons_seqs[hap-1]);
        //             // XXX previous abpoa obj could be saved for 2nd round
        //             collect_noisy_cons_seqs_with_ps_hap0(opt, n_reads, lens, seqs, qnames, haps, phase_sets, fully_covers, ps, hap, cons_lens+hap-1, cons_seqs+hap-1);
        //         }
        //     }
        // }
        // for (int i = 0; i < 2; ++i) free(clu_read_ids[i]);
    // }
collect_noisy_msa_cons_no_ps_hap_end:
    free(full_read_lens); free(full_read_seqs); free(full_read_names); free(full_fully_covers); free(full_read_ids);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
    } free(cons_lens); free(cons_seqs);
    return n_cons;
}

// collect phase set with both haps having >= min_minor_hap_read_count reads
hts_pos_t collect_phase_set_with_both_haps(int n_reads, int *read_haps, hts_pos_t *phase_sets, int *fully_covers, int min_minor_hap_read_count) {
    int n_uniq_phase_sets = 0, phase_set_i = 0;
    hts_pos_t *uniq_phase_sets = (hts_pos_t*)calloc(n_reads, sizeof(hts_pos_t));
    int **phase_set_to_hap_read_count = (int**)malloc(n_reads * sizeof(int*));
    for (int i = 0; i < n_reads; ++i) {
        phase_set_to_hap_read_count[i] = (int*)calloc(2, sizeof(int));
    }
    // use one phase set if multiple phase sets exist, reads in other phase sets are considered as non-HAP
    for (int i = 0; i < n_reads; ++i) {
        if (read_haps[i] == 0) continue;
        if (fully_covers[i] != 3) continue;
        phase_set_i = add_phase_set(phase_sets[i], uniq_phase_sets, &n_uniq_phase_sets);
        phase_set_to_hap_read_count[phase_set_i][read_haps[i]-1]++;
    }
    hts_pos_t max_ps = -1; int max_ps_read_count1 = -1, max_ps_read_count2 = -1;
    for (int i = 0; i < n_uniq_phase_sets; ++i) {
        int phase_set_read_count1 = phase_set_to_hap_read_count[i][0] < phase_set_to_hap_read_count[i][1] ? phase_set_to_hap_read_count[i][0] : phase_set_to_hap_read_count[i][1];
        int phase_set_read_count2 = phase_set_to_hap_read_count[i][0] > phase_set_to_hap_read_count[i][1] ? phase_set_to_hap_read_count[i][0] : phase_set_to_hap_read_count[i][1];
        if (phase_set_read_count1 > max_ps_read_count1) {
            max_ps_read_count1 = phase_set_read_count1;
            max_ps_read_count2 = phase_set_read_count2;
            max_ps = uniq_phase_sets[i];
        } else if (phase_set_read_count1 == max_ps_read_count1 && phase_set_read_count2 > max_ps_read_count2) {
            max_ps_read_count1 = phase_set_read_count1;
            max_ps_read_count2 = phase_set_read_count2;
            max_ps = uniq_phase_sets[i];
        }
    }
    for (int i = 0; i < n_reads; ++i) free(phase_set_to_hap_read_count[i]);
    free(uniq_phase_sets); free(phase_set_to_hap_read_count);
    if (max_ps_read_count1 >= min_minor_hap_read_count) return max_ps;
    else return -1;
}

// pairwise alignment of cons and ref to append ref to the MSA
// total: 1+ n_reads*2
// 1. ref vs cons: 1
// 2. cons vs n_reads: n_reads
// 3. ref vs n_reads: n_reads
int wfa_collect_noisy_aln_str_with_ps_hap(const call_var_opt_t *opt, int n_reads, int *noisy_read_ids, int *lens, uint8_t **seqs, char **names,
                                          int *haps, hts_pos_t *phase_sets, int *fully_covers, hts_pos_t ps, uint8_t *ref_seq, int ref_seq_len,
                                          int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs) {
    // given specific PS, collect consensus sequences for each haplotype
    int n_cons = 0;
    int total_n_reads = n_reads+2, n_ps_hap_reads = 0;
    int *ps_hap_read_ids = (int*)malloc(total_n_reads * sizeof(int));
    int *ps_hap_read_lens = (int*)malloc(total_n_reads * sizeof(int));
    uint8_t **ps_hap_read_seqs = (uint8_t**)malloc(total_n_reads * sizeof(uint8_t*));
    char **ps_hap_read_names = (char**)malloc(total_n_reads * sizeof(char*));
    int *ps_hap_full_covers = (int*)calloc(total_n_reads, sizeof(int));
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    for (int hap=1; hap<=2; ++hap) {
        cons_seqs[hap-1] = NULL; cons_lens[hap-1] = 0;
        n_ps_hap_reads = 0;
        for (int i = 0; i < n_reads; ++i) {
            if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
            if (opt->disable_read_realign) {
                if (fully_covers[i] != 3) continue;
            }
            ps_hap_read_ids[n_ps_hap_reads] = noisy_read_ids[i];
            ps_hap_read_lens[n_ps_hap_reads] = lens[i];
            ps_hap_read_seqs[n_ps_hap_reads] = seqs[i];
            ps_hap_full_covers[n_ps_hap_reads] = fully_covers[i];
            ps_hap_read_names[n_ps_hap_reads] = names[i];
            n_ps_hap_reads++;
        }
        // only call consensus at first
        if (LONGCALLD_VERBOSE >=2 ) fprintf(stderr, "PS: %ld HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
        if (n_ps_hap_reads == 0) continue;
        // collect consensus sequences
        n_cons += abpoa_partial_aln_msa_cons(opt, NULL, 10, n_ps_hap_reads, ps_hap_read_ids, ps_hap_read_seqs, ps_hap_read_lens, ps_hap_full_covers, ps_hap_read_names,
                                             1, cons_lens+hap-1, cons_seqs+hap-1, clu_n_seqs+hap-1, clu_read_ids+hap-1, NULL, NULL);
    }

    if (n_cons != 2) n_cons = 0;
    // XXX add clipped / partial aligned reads to the MSA
    // if (opt->disable_read_realign == 0) { // for other reads with other PS or no PS, align to the consensus and then determine the haplotype
    //     int assign_hap, new_haps[2]={0, 0};
    //     for (int i = 0; i < n_reads; ++i) {
    //         if (lens[i] > 0 && fully_covers[i] != 0 && phase_sets[i] != ps) {
    //             assign_hap = determine_hap_based_on_noisy_cons(cons_lens, cons_seqs, n_cons, lens[i], seqs[i], fully_covers[i]);
    //             if (assign_hap > 0) {
    //                 if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "AssignHap: %s %d\n", names[i], assign_hap);
    //                 haps[i] = assign_hap; phase_sets[i] = ps;
    //                 new_haps[assign_hap-1] = 1;
    //             }
    //         }
    //     }
    //     // re-do consensus calling using all reads
    //     for (int hap=1; hap<=2; ++hap) {
    //         if (new_haps[hap-1]) {
    //             // free(cons_seqs[hap-1]);
    //             // XXX append new sequences to the MSA
    //             // collect_noisy_msa_cons_with_ps_hap0(opt, n_reads, noisy_read_ids, lens, seqs, names, haps, phase_sets, fully_covers, ps, hap,
    //                                                 // cons_lens+hap-1, cons_seqs+hap-1, clu_n_seq+hap-1, clu_read_ids+hap-1, msa_seq_lens+hap-1, msa_seqs[hap-1]);
    //         }
    //     }
    // }
    if (n_cons == 2) { // collect aln_strs
        for (int hap=1; hap<=2; ++hap) {
            // ref vs cons
            aln_str_t *clu_aln_str = aln_strs[hap-1];
            wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[hap-1], cons_lens[hap-1], LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
            n_ps_hap_reads = 0;
            for (int i = 0; i < n_reads; ++i) {
                if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
                if (opt->disable_read_realign) {
                    if (fully_covers[i] != 3) continue;
                }
                // cons vs read
                // use abpoa-MSA to derive aln_str of cons vs reads, ensure reads supporting the consensus also support all the variants derived from the consensus
                wfa_collect_aln_str(opt, cons_seqs[hap-1], cons_lens[hap-1], seqs[i], lens[i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                wfa_collect_aln_str(opt, ref_seq, ref_seq_len, seqs[i], lens[i], LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                // ps_hap_read_lens[n_ps_hap_reads] = lens[i];
                // ps_hap_read_seqs[n_ps_hap_reads] = seqs[i];
                // ps_hap_full_covers[n_ps_hap_reads] = fully_covers[i];
                n_ps_hap_reads++;
            }
            // only MSA
            if (LONGCALLD_VERBOSE >=2 ) fprintf(stderr, "With Ref+Cons PS: %ld HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
            if (n_ps_hap_reads == 0) continue;
        }
    }
    free(ps_hap_read_ids); free(ps_hap_read_lens); free(ps_hap_read_seqs); free(ps_hap_full_covers); free(ps_hap_read_names);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
    } free(cons_lens); free(cons_seqs);
    return n_cons;
}

int collect_noisy_read_info(bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reg_reads, int **read_lens, uint8_t ***read_seqs, char ***read_names, int **fully_covers, int **read_haps, hts_pos_t **phase_sets) {
    *read_lens = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_seqs = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    *read_names = (char**)malloc(n_noisy_reg_reads * sizeof(char*));
    *fully_covers = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_haps = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *phase_sets = (hts_pos_t*)calloc(n_noisy_reg_reads, sizeof(hts_pos_t));

    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        int read_i = noisy_reg_reads[i];
        digar_t *read_digars = chunk->digars+read_i; int n_digar = read_digars->n_digar; digar1_t *digars = read_digars->digars;
        hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
        int reg_read_beg = 0, reg_read_end = bam_cigar2qlen(chunk->reads[read_i]->core.n_cigar, bam_get_cigar(chunk->reads[read_i]))-1;
        (*read_names)[i] = bam_get_qname(chunk->reads[read_i]);
        for (int i = 0; i < n_digar; ++i) {
            hts_pos_t digar_beg = digars[i].pos, digar_end;
            int op = digars[i].type, len = digars[i].len, qi = digars[i].qi;
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) continue;
            if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_end = digar_beg + len - 1;
            else digar_end = digar_beg;
            if (digar_beg > reg_end) break;
            if (digar_end < reg_beg) continue;
            if (digar_beg <= reg_beg && digar_end >= reg_beg) {
                if (op == BAM_CDEL) {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi;
                } else {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi + (reg_beg - digar_beg);
                }
            }
            if (digar_beg <= reg_end && digar_end >= reg_end) {
                if (op == BAM_CDEL) {
                    reg_digar_end = reg_end;
                    reg_read_end = qi-1;
                } else {
                    reg_digar_end = reg_end;
                    reg_read_end = qi + (reg_end - digar_beg);
                }
            }
        }
        if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) (*fully_covers)[i] = 3;
        else if (reg_digar_beg == reg_beg) (*fully_covers)[i] = 1;
        else if (reg_digar_end == reg_end) (*fully_covers)[i] = 2;
        else (*fully_covers)[i] = 0;
        // if (2*(reg_read_end-reg_read_beg+1) < (reg_end-reg_beg+1)) return 0;
        (*read_seqs)[i] = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
        for (int j = reg_read_beg; j <= reg_read_end; ++j) {
            (*read_seqs)[i][j-reg_read_beg] = seq_nt16_int[bam_seqi(read_digars->bseq, j)];
        }
        (*read_lens)[i] = reg_read_end - reg_read_beg + 1;
        (*read_haps)[i] = chunk->haps[read_i];
        (*phase_sets)[i] = chunk->PS[read_i];
    }
    return 0;
}

// return n_cons
// 1. consensu calling; 2. WFA-based MSA
int collect_noisy_reg_aln_strs(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reads, uint8_t *ref_seq, int ref_seq_len,
                               int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs) {
    if (n_noisy_reg_reads <= 0) return 0;
    // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
    char **names = NULL; uint8_t **seqs = NULL; int *fully_covers = NULL, *lens = NULL, *haps = NULL; hts_pos_t *phase_sets = NULL;
    collect_noisy_read_info(chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, n_noisy_reg_reads, noisy_reads, &lens, &seqs, &names, &fully_covers, &haps, &phase_sets);

    sort_by_full_cover_and_length(n_noisy_reg_reads, noisy_reads, lens, seqs, fully_covers, names, haps, phase_sets);
    // at least 5 reads for each hap
    int min_hap_full_read_count = opt->min_hap_full_reads, min_no_hap_full_read_count = opt->min_no_hap_full_reads;
    hts_pos_t ps_with_both_haps = collect_phase_set_with_both_haps(n_noisy_reg_reads, haps, phase_sets, fully_covers, min_hap_full_read_count);
    int n_full_reads = 0;
    for (int i = 0; i < n_noisy_reg_reads; ++i) if (fully_covers[i] == 3) n_full_reads++;
    int n_cons = 0;
    // XXX only use fully-covered reads, including cliping reads (after re-align to backbone read)
    // two cases to call consensus sequences
    if (ps_with_both_haps > 0) { // call consensus sequences for each haplotype
        n_cons = wfa_collect_noisy_aln_str_with_ps_hap(opt, n_noisy_reg_reads, noisy_reads, lens, seqs, names, haps, phase_sets, fully_covers, ps_with_both_haps,
                                                       ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs);
    // } else if (n_full_cover_reads > n_noisy_reg_reads * 0.75 && reg_end - reg_beg + 1 <= 10000) { // de novo consensus calling using all reads, up to 2 consensus sequences
    } else if (ps_with_both_haps <= 0 && n_full_reads >= min_no_hap_full_read_count) {
        n_cons = wfa_collect_noisy_aln_str_no_ps_hap(opt, n_noisy_reg_reads, noisy_reads, lens, seqs, names, fully_covers,
                                                     ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs);
    } else { // if (n_full_cover_reads <= 0 || reg_end - reg_beg + 1 > 5000) {
        fprintf(stderr, "Skipped region: %s:%ld-%ld %ld %d reads (%d full)\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads);
    }
    for (int i = 0; i < n_noisy_reg_reads; ++i) free(seqs[i]);
    free(names); free(seqs); free(lens); free(fully_covers); free(haps); free(phase_sets);
    return n_cons;
}