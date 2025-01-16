#include "wavefront/wavefront_align.h"
#include "alignment/cigar.h"
#include "seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"
#include "edlib.h"
#include "ksw2.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "utils.h"
#include "align.h"

extern int LONGCALLD_VERBOSE;

// return alignment score
// ext_direction: 0 -> global, 1 -> left to right, 2 -> right to left
int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen) {
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

int wfa_aln(int gap_aln, uint8_t *pattern, int plen, uint8_t *text, int tlen, int a, int b, int q, int e, int q2, int e2, uint32_t **cigar_buf) {
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
        for (int i = 0; i < plen; ++i) {
            p[i] = pattern[plen-i-1];
        }
        for (int i = 0; i < tlen; ++i) {
            t[i] = text[tlen-i-1];
        }
    }

    // Align
    wavefront_align(wf_aligner, (const char*)p, plen, (const char*)t, tlen);
    cigar_t *cigar = wf_aligner->cigar;
    if (LONGCALLD_VERBOSE >= 2) {
        char *p_seq = (char*)malloc(plen+1); char *t_seq = (char*)malloc(tlen+1);
        for (int i = 0; i < plen; ++i) p_seq[i] = "ACGTN"[p[i]];
        for (int i = 0; i < tlen; ++i) t_seq[i] = "ACGTN"[t[i]];
        cigar_print_pretty(stderr, cigar, p_seq, plen, t_seq, tlen);
        free(p_seq); free(t_seq);
    }
    // collect cigar
    uint32_t *tmp_cigar; int cigar_length;
    cigar_get_CIGAR(wf_aligner->cigar,true, &tmp_cigar, &cigar_length);
    *cigar_buf = (uint32_t*)malloc(cigar_length * sizeof(uint32_t));
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { // reverse cigar
        for (int i = 0; i < cigar_length; ++i) {
            (*cigar_buf)[i] = tmp_cigar[cigar_length-i-1];
        }
    } else {
        for (int i = 0; i < cigar_length; ++i) {
            (*cigar_buf)[i] = tmp_cigar[i];
        }
    }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { free(p); free(t); }
    return cigar_length;
}

static void ksw_gen_simple_mat(int m, int8_t *mat, int8_t a, int8_t b) {
    int i, j;
    a = a < 0? -a : a;
    b = b > 0? -b : b;
    for (i = 0; i < m - 1; ++i) {
        for (j = 0; j < m - 1; ++j)
            mat[i * m + j] = i == j? a : b;
        mat[i * m + m - 1] = 0;
    }
    for (j = 0; j < m; ++j)
        mat[(m - 1) * m + j] = 0;
}

int ksw2_aln(int gap_aln, uint8_t *tseq, int tlen, uint8_t *qseq, int qlen, int a, int b, int q, int e, int q2, int e2, uint32_t **cigar_buf) {
    ksw_extz_t ez; uint8_t c[256];
    int w = -1, zdrop = -1;
    int8_t mat[25]; ksw_gen_simple_mat(5, mat, a, -b);
    memset(&ez, 0, sizeof(ksw_extz_t));
    int flag = KSW_EZ_EQX;
    if (gap_aln == LONGCALLD_GAP_RIGHT_ALN) flag |= KSW_EZ_RIGHT;
    ksw_extd2_sse(0, qlen, qseq, tlen, tseq, 5, mat, q, e, q2, e2, w, zdrop, 0, flag, &ez);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "ksw2-score: %d\n", ez.score);
        for (int i = 0; i < ez.n_cigar; ++i) {
            int op = ez.cigar[i] & 0xf; int len = ez.cigar[i] >> 4;
            fprintf(stderr, "%d%c", len, BAM_CIGAR_STR[op]);
        }
        fprintf(stderr, "\n");
        for (int i = 0; i < tlen; ++i) fprintf(stderr, "%c", "ACGTN"[tseq[i]]);
        fprintf(stderr, "\n");
        for (int i = 0; i < qlen; ++i) fprintf(stderr, "%c", "ACGTN"[qseq[i]]);
        fprintf(stderr, "\n");
    }
    *cigar_buf = ez.cigar;
    return ez.n_cigar;
}

int end2end_aln(const call_var_opt_t *opt, char *tseq, int tlen, uint8_t *qseq, int qlen, uint32_t **cigar_buf) {
    if (qlen <= 0 || tlen <= 0) return 0;
    int min_len = MIN_OF_TWO(tlen, qlen), max_len = MAX_OF_TWO(tlen, qlen);
    int delta_len = MAX_OF_TWO(1, max_len - min_len);

    uint8_t *tseq2 = (uint8_t*)malloc(max_len); int cigar_len = 0;
    for (int i = 0; i < tlen; ++i) tseq2[i] = nst_nt4_table[(uint8_t)tseq[i]];
    // use wfa if the length difference is small
    // if (max_len * delta_len + delta_len * delta_len < max_len * min_len) {
        cigar_len = wfa_aln(opt->gap_aln, tseq2, tlen, qseq, qlen, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, cigar_buf);
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

// reads are sorted by full-cover and length before calling abPOA
// use reads with full_cover == 3 as backbone for poa
// skip reads with full_cover == 0
int abpoa_partial_aln(int n_reads, uint8_t **read_seqs, int *read_lens, int *read_full_cover, char **names, int max_n_cons, int *cons_lens, uint8_t **cons_seqs) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    // abpt->wb = -1;
    abpt->out_msa = 0;
    abpt->cons_algrm = ABPOA_MF;
    abpt->sort_input_seq = 1;
    abpt->max_n_cons = max_n_cons;
    // msa
    abpt->out_cons = 1;
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
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Skip: %s\n", names[i]);
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

 int abpoa_aln(int n_reads, uint8_t **read_seqs, int *read_lens, int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int start_cons_i, int **n_clu_reads, int ***clu_read_ids) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abp = abpoa_init_para();
    // abpt->wb = -1;
    abp->out_msa = 0;
    if (LONGCALLD_VERBOSE >= 2) abp->out_msa = 1;
    abp->cons_algrm = ABPOA_MF;
    abp->max_n_cons = max_n_cons;
    // msa
    abp->out_cons = 1;
    abpoa_post_set_para(abp);

    // sort reads by length
    uint8_t **sorted_read_seqs = (uint8_t**)malloc(n_reads * sizeof(uint8_t*));
    int *sorted_read_lens = (int*)malloc(n_reads * sizeof(int));
    int *sorted_seq_to_ids = (int*)malloc(n_reads * sizeof(int));

    for (int i = 0; i < n_reads; ++i) {
        sorted_seq_to_ids[i] = i;
        sorted_read_lens[i] = read_lens[i];
        sorted_read_seqs[i] = read_seqs[i];
    }
    for (int i = 0; i < n_reads-1; ++i) {
        for (int j = i+1; j < n_reads; ++j) {
            if (sorted_read_lens[i] < sorted_read_lens[j]) {
                int tmp_len = sorted_read_lens[i]; sorted_read_lens[i] = sorted_read_lens[j]; sorted_read_lens[j] = tmp_len;
                uint8_t *tmp_seq = sorted_read_seqs[i]; sorted_read_seqs[i] = sorted_read_seqs[j]; sorted_read_seqs[j] = tmp_seq;
                int tmp_id = sorted_seq_to_ids[i]; sorted_seq_to_ids[i] = sorted_seq_to_ids[j]; sorted_seq_to_ids[j] = tmp_id;
            }
        }
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA (max %d cons): %d\n", max_n_cons, n_reads);
        abpoa_msa(ab, abp, n_reads, NULL, sorted_read_lens, sorted_read_seqs, NULL, stderr);
    } else abpoa_msa(ab, abp, n_reads, NULL, sorted_read_lens, sorted_read_seqs, NULL, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    int n_cons = 0;
    if (abc->n_cons > 0) {
        for (int i = 0; i < abc->n_cons; ++i) {
            cons_lens[start_cons_i+i] = abc->cons_len[i];
            cons_seqs[start_cons_i+i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
            for (int j = 0; j < abc->cons_len[i]; ++j) {
                cons_seqs[start_cons_i+i][j] = abc->cons_base[i][j];
            }
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[i]);
        }
        if (abc->n_cons == 2 && clu_read_ids != NULL) { // collect read ids for each cluster
            *clu_read_ids = (int**)malloc(2 * sizeof(int*));
            for (int i = 0; i < 2; ++i) {
                (*n_clu_reads)[i] = abc->clu_n_seq[i];
                (*clu_read_ids)[i] = (int*)malloc(abc->clu_n_seq[i] * sizeof(int));
                for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                    (*clu_read_ids)[i][j] = sorted_seq_to_ids[abc->clu_read_ids[i][j]];
                }
            }
        }
        n_cons = abc->n_cons;
    } else {
        _err_error_exit("No Consensus: %ld\n");
    }
    free(sorted_seq_to_ids); free(sorted_read_lens); free(sorted_read_seqs);
    abpoa_free(ab); abpoa_free_para(abp);
    return n_cons;
 }


int *collect_noisy_read_haps(bam_chunk_t *chunk, int noisy_reg_i) {
    int n_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    int *haps = (int*)calloc(n_reads, sizeof(int));
    int n_hap1 = 0, n_hap2 = 0;
    for (int j = 0; j < chunk->noisy_reg_to_n_reads[noisy_reg_i]; ++j) {
        int read_i = chunk->noisy_reg_to_reads[noisy_reg_i][j];
        haps[j] = chunk->haps[read_i];
        if (haps[j] == 1) n_hap1++;
        else if (haps[j] == 2) n_hap2++;
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "NOISY_REG: %d, n_reads: %d, n_hap1: %d, n_hap2: %d\n", noisy_reg_i, n_reads, n_hap1, n_hap2);
    }
    return haps;
}

hts_pos_t *collect_noisy_read_phase_sets(bam_chunk_t *chunk, int noisy_reg_i) {
    int n_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    hts_pos_t *phase_sets = (hts_pos_t*)calloc(n_reads, sizeof(hts_pos_t));
    for (int j = 0; j < chunk->noisy_reg_to_n_reads[noisy_reg_i]; ++j) {
        int read_i = chunk->noisy_reg_to_reads[noisy_reg_i][j];
        phase_sets[j] = chunk->PS[read_i];
    }
    return phase_sets;
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

int collect_reg_read_seq(bam_chunk_t *chunk, int read_i, hts_pos_t reg_beg, hts_pos_t reg_end, uint8_t **reg_seq, char **qname, int *fully_cover) {
    digar_t *read_digars = chunk->digars+read_i; int n_digar = read_digars->n_digar; digar1_t *digars = read_digars->digars;
    hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
    int reg_read_beg = 0, reg_read_end = bam_cigar2qlen(chunk->reads[read_i]->core.n_cigar, bam_get_cigar(chunk->reads[read_i]))-1;
    *qname = bam_get_qname(chunk->reads[read_i]);
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
    if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) *fully_cover = 3;
    else if (reg_digar_beg == reg_beg) *fully_cover = 1;
    else if (reg_digar_end == reg_end) *fully_cover = 2;
    else *fully_cover = 0;
    // if (2*(reg_read_end-reg_read_beg+1) < (reg_end-reg_beg+1)) return 0;
    *reg_seq = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
    for (int i = reg_read_beg; i <= reg_read_end; ++i) {
        (*reg_seq)[i-reg_read_beg] = seq_nt16_int[bam_seqi(read_digars->bseq, i)];
    }
    // fprintf(stderr, "%s\t%ld-%ld: %d-%d\t%d\tHAP%d\n", bam_get_qname(chunk->reads[read_i]), reg_beg, reg_end, reg_read_beg, reg_read_end, *fully_cover, chunk->haps[read_i]);
    // for (int i = 0; i < reg_read_end - reg_read_beg + 1; ++i) {
        // fprintf(stderr, "%c", "ACGTN"[(*reg_seq)[i]]);
    // } fprintf(stderr, "\n");

    return reg_read_end - reg_read_beg + 1;
}

int collect_msa_seqs(uint8_t *ref_seq, int ref_seq_len, uint8_t **cons_seqs, int *cons_lens, int n_cons, uint8_t ***msa_seqs) {
    int msa_seq_len = 0;
    int n_seqs = n_cons + 1;
    int *seq_lens = (int*)malloc(n_seqs * sizeof(int));
    uint8_t **seqs = (uint8_t**)malloc(n_seqs * sizeof(uint8_t*));
    for (int i = 0; i < n_cons+1; ++i) {
        if (i == 0) {
            seqs[i] = ref_seq;
            seq_lens[i] = ref_seq_len;
            // fprintf(stderr, "Ref: %d\n", ref_seq_len);
        } else {
            seqs[i] = cons_seqs[i-1];
            seq_lens[i] = cons_lens[i-1];
            // fprintf(stderr, "Cons%d: %d\n", i, cons_lens[i-1]);
        }
    }
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abp = abpoa_init_para();
    abp->wb = -1; // for variant calling, disable banded alignment
    abp->out_msa = 1;
    abp->out_cons = 0;
    abpoa_post_set_para(abp);
    if (LONGCALLD_VERBOSE >= 2) abpoa_msa(ab, abp, n_seqs, NULL, seq_lens, seqs, NULL, stderr);
    else abpoa_msa(ab, abp, n_seqs, NULL, seq_lens, seqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    *msa_seqs = (uint8_t**)malloc(n_seqs * sizeof(uint8_t**));
    msa_seq_len = abc->msa_len;
    for (int i = 0; i < n_seqs; ++i) {
        (*msa_seqs)[i] = (uint8_t*)malloc(msa_seq_len * sizeof(uint8_t));
        for (int j = 0; j < msa_seq_len; ++j) {
            (*msa_seqs)[i][j] = abc->msa_base[i][j];
        }
    }
    abpoa_free(ab); abpoa_free_para(abp); free(seqs); free(seq_lens);
    return msa_seq_len;
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

void sort_by_full_cover_and_length(int n_reads, int *read_lens, uint8_t **read_seqs, int *full_covers, char **names) {
    for (int i = 0; i < n_reads-1; ++i) {
        for (int j = i+1; j < n_reads; ++j) {
            if (full_covers[i] < full_covers[j]) {
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                int tmp_cover = full_covers[i]; full_covers[i] = full_covers[j]; full_covers[j] = tmp_cover;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
            } else if (full_covers[i] == full_covers[j] && read_lens[i] < read_lens[j]) {
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
            }
        }
    }
}

// given specific PS & hap, collect consensus sequence
int collect_noisy_cons_seqs_with_ps_hap0(const call_var_opt_t *opt, int n_reads, int *lens, uint8_t **seqs, char **names,
                                         int *read_haps, hts_pos_t *phase_sets, int *fully_covers,
                                         hts_pos_t ps, int hap, int *cons_lens, uint8_t **cons_seqs) {
    int n_ps_hap_reads = 0;
    int *ps_hap_read_lens = (int*)malloc(n_reads * sizeof(int));
    uint8_t **ps_hap_read_seqs = (uint8_t**)malloc(n_reads * sizeof(uint8_t*));
    char **ps_hap_read_names = (char**)malloc(n_reads * sizeof(char*));
    int *ps_hap_full_covers = (int*)calloc(n_reads, sizeof(int));
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] <= 0 || phase_sets[i] != ps || read_haps[i] != hap) continue;
        if (opt->disable_read_realign) {
            if (fully_covers[i] != 3) continue;
        }
        ps_hap_read_lens[n_ps_hap_reads] = lens[i];
        ps_hap_read_seqs[n_ps_hap_reads] = seqs[i];
        ps_hap_full_covers[n_ps_hap_reads] = fully_covers[i];
        ps_hap_read_names[n_ps_hap_reads] = names[i];
        n_ps_hap_reads++;
    }
    if (n_ps_hap_reads == 0) return 1;
    // sort by full cover and length
    sort_by_full_cover_and_length(n_ps_hap_reads, ps_hap_read_lens, ps_hap_read_seqs, ps_hap_full_covers, ps_hap_read_names);
    if (LONGCALLD_VERBOSE >=2 ) fprintf(stderr, "PS: %ld HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
    int n_cons = abpoa_partial_aln(n_ps_hap_reads, ps_hap_read_seqs, ps_hap_read_lens, ps_hap_full_covers, ps_hap_read_names, 1, cons_lens, cons_seqs);
    free(ps_hap_read_lens); free(ps_hap_read_seqs); free(ps_hap_full_covers); free(ps_hap_read_names);
    return n_cons;
}

int determine_hap_based_on_noisy_cons(int *cons_lens, uint8_t **cons_seqs, int n_cons, int read_len, uint8_t *read_seq, int full_cover) {
    if (n_cons == 0) return 0;
    int scores[2] = {0, 0};
    for (int i = 0; i < n_cons; ++i) { // 0'th cons -> HAP1, 1'th cons -> HAP2
        int score = INT32_MIN;
        if (full_cover == 3) {
            score = wfa_heuristic_aln(cons_seqs[i], cons_lens[i], read_seq, read_len);
        } else if (full_cover == 1) {
            int len = MIN_OF_TWO(cons_lens[i], read_len);
            score = wfa_heuristic_aln(cons_seqs[i], len, read_seq, len);
        } else if (full_cover == 2) {
            int len = MIN_OF_TWO(cons_lens[i], read_len);
            score = wfa_heuristic_aln(cons_seqs[i]+cons_lens[i]-len, len, read_seq+read_len-len, len);
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

// XXX de novo clustering based on read length if possible
// XXX only use fully-covered reads for now
int collect_noisy_cons_seq_no_ps_hap(const call_var_opt_t *opt, hts_pos_t ps, int n_reads, int *lens, uint8_t **seqs, char **qnames, int *haps, hts_pos_t *phase_sets, int *fully_covers, int *cons_lens, uint8_t **cons_seqs) {
    int *full_read_lens = (int*)malloc(n_reads * sizeof(int));
    uint8_t **full_read_seqs = (uint8_t**)malloc(n_reads * sizeof(uint8_t*));
    char **full_read_names = (char**)malloc(n_reads * sizeof(char*));
    int *full_read_to_idx = (int*)malloc(n_reads * sizeof(int));
    int n_full_reads = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] <= 0 || fully_covers[i] != 3) continue;
        full_read_lens[n_full_reads] = lens[i];
        full_read_seqs[n_full_reads] = seqs[i];
        full_read_to_idx[n_full_reads] = i;
        full_read_names[n_full_reads] = qnames[i];
        n_full_reads++;
    }
    // abpoa
    int **clu_read_ids = (int**)malloc(2 * sizeof(int*)), *n_clu_reads = (int*)calloc(2, sizeof(int));
    int n_cons = abpoa_aln(n_full_reads, full_read_seqs, full_read_lens, 2, cons_lens, cons_seqs, 0, &n_clu_reads, &clu_read_ids);
    if (n_cons == 2) {
        // XXX update PS & hap if n_cons == 2
        int n_old_haps[2] = {0, 0};
        for (int hap = 1; hap <= 2; ++hap) {
            for (int i = 0; i < n_clu_reads[hap-1]; ++i) {
                int read_i = full_read_to_idx[clu_read_ids[hap-1][i]];
                phase_sets[read_i] = ps; haps[read_i] = hap;
                n_old_haps[hap-1] += 1;
                // fprintf(stderr, "UpdatePSHap: %s %d %d %d\n", qnames[read_i], hap, ps, clu_read_ids[hap-1][i]);
            }
        }
        // XXX for other reads, align to the consensus and then determine the haplotype
        int assign_hap, n_new_haps[2]={0, 0};
        for (int i = 0; i < n_reads; ++i) {
            if (lens[i] > 0 && (fully_covers[i] == 1 || fully_covers[i] == 2)) {
                assign_hap = determine_hap_based_on_noisy_cons(cons_lens, cons_seqs, n_cons, lens[i], seqs[i], fully_covers[i]);
                if (assign_hap > 0) {
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "AssignHap: %s %d\n", qnames[i], assign_hap);
                    haps[i] = assign_hap; phase_sets[i] = ps;
                    n_new_haps[assign_hap-1] += 1;
                }
            }
        }
        // re-do consensus calling using all reads
        if (opt->disable_read_realign == 0) {
            for (int hap=1; hap<=2; ++hap) {
                if (n_new_haps[hap-1] >= n_old_haps[hap-1] / 4) {
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Re-Do Consensus Calling: %d %d\n", n_old_haps[hap-1], n_new_haps[hap-1]);
                    free(cons_seqs[hap-1]);
                    // XXX previous abpoa obj could be saved for 2nd round
                    collect_noisy_cons_seqs_with_ps_hap0(opt, n_reads, lens, seqs, qnames, haps, phase_sets, fully_covers, ps, hap, cons_lens+hap-1, cons_seqs+hap-1);
                }
            }
        }
        for (int i = 0; i < 2; ++i) free(clu_read_ids[i]);
    }
    free(clu_read_ids); free(n_clu_reads);
    free(full_read_lens); free(full_read_seqs);
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


int collect_noisy_read_uniq_ps(int n_noisy_reg_reads, hts_pos_t *phase_sets, hts_pos_t *uniq_ps) {
    int *uniq_ps_count = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    int n_uniq_ps = 0;
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        if (phase_sets[i] == 0) continue;
        int j, ps_i = -1;
        for (j = 0; j < n_uniq_ps; ++j) {
            if (uniq_ps[j] == phase_sets[i]) {
                ps_i = j;
                break;
            }
        }
        if (j == n_uniq_ps) {
            uniq_ps_count[n_uniq_ps] = 1; 
            uniq_ps[n_uniq_ps++] = phase_sets[i];
        } else {
            uniq_ps_count[ps_i]++;
        }
    }
    // sort by count
    for (int i = 0; i < n_uniq_ps-1; ++i) {
        for (int j = i+1; j < n_uniq_ps; ++j) {
            if (uniq_ps_count[i] < uniq_ps_count[j]) {
                int tmp_count = uniq_ps_count[i]; uniq_ps_count[i] = uniq_ps_count[j]; uniq_ps_count[j] = tmp_count;
                hts_pos_t tmp_ps = uniq_ps[i]; uniq_ps[i] = uniq_ps[j]; uniq_ps[j] = tmp_ps;
            }
        }
    }
    free(uniq_ps_count);
    return n_uniq_ps;
}

// given specific PS, collect consensus sequences for each haplotype
// for other reads with other PS, align to the consensus and then determine the haplotype
int collect_noisy_cons_seqs_with_ps_hap(const call_var_opt_t *opt, int n_reads, int *lens, uint8_t **seqs, char **names,
                                        int *haps, hts_pos_t *phase_sets, int *fully_covers,
                                        hts_pos_t ps, int *cons_lens, uint8_t **cons_seqs) {
    int n_cons = 0;
    // given specific PS, collect consensus sequences for each haplotype
    for (int hap=1; hap<=2; ++hap) {
        n_cons += collect_noisy_cons_seqs_with_ps_hap0(opt, n_reads, lens, seqs, names, haps, phase_sets, fully_covers, ps, hap, cons_lens+hap-1, cons_seqs+hap-1);
    }
    // for other reads with other PS or no PS, align to the consensus and then determine the haplotype
    int assign_hap, new_haps[2]={0, 0};
    assert(n_cons == 2);
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] > 0 && fully_covers[i] != 0 && phase_sets[i] != ps) {
            assign_hap = determine_hap_based_on_noisy_cons(cons_lens, cons_seqs, n_cons, lens[i], seqs[i], fully_covers[i]);
            if (assign_hap > 0) {
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "AssignHap: %s %d\n", names[i], assign_hap);
                haps[i] = assign_hap; phase_sets[i] = ps;
                new_haps[assign_hap-1] = 1;
            }
        }
    }
    // re-do consensus calling using all reads
    if (opt->disable_read_realign == 0) {
        for (int hap=1; hap<=2; ++hap) {
            if (new_haps[hap-1]) {
                free(cons_seqs[hap-1]);
                // XXX previous abpoa obj could be saved for 2nd round
                collect_noisy_cons_seqs_with_ps_hap0(opt, n_reads, lens, seqs, names, haps, phase_sets, fully_covers, ps, hap, cons_lens+hap-1, cons_seqs+hap-1);
            }
        }
    }
    return n_cons;
}

int collect_noisy_reads(bam_chunk_t *chunk, int noisy_reg_i, int **read_lens, uint8_t ***read_seqs, char ***read_names, int **fully_covers, int **read_haps, hts_pos_t **phase_sets) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int *noisy_reg_reads = chunk->noisy_reg_to_reads[noisy_reg_i], n_noisy_reg_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    hts_pos_t reg_beg = cr_start(noisy_regs, noisy_reg_i), reg_end = cr_end(noisy_regs, noisy_reg_i);
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

int collect_noisy_full_hap_reads(int n_noisy_reg_reads, int *haps, hts_pos_t *phase_sets, int *fully_covers, hts_pos_t *ps_with_full_hap_reads) {
    int n_full_cover_reads = 0; *ps_with_full_hap_reads = -1;
    int n_ps = 0; hts_pos_t *uniq_ps = (hts_pos_t*)calloc(n_noisy_reg_reads, sizeof(hts_pos_t));
    int *uniq_ps_hap = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        if (fully_covers[i] != 3) continue;
        n_full_cover_reads++;
        if (haps[i] == 0) continue;
        hts_pos_t ps = phase_sets[i]; int hap = haps[i];
        int j;
        for (j = 0; j < n_ps; ++j) {
            if (uniq_ps[j] == ps) {
                if (uniq_ps_hap[j] != hap) {
                    *ps_with_full_hap_reads = ps;
                    break;
                }
            }
        }
        if (j == n_ps) {
            uniq_ps[n_ps] = ps;
            uniq_ps_hap[n_ps++] = hap;
        }
        if (*ps_with_full_hap_reads != -1) break;
    }
    free(uniq_ps); free(uniq_ps_hap);
    return n_full_cover_reads;
}

void updated_ps_hap(bam_chunk_t *chunk, int noisy_reg_i, int n_noisy_reg_reads, int *haps, hts_pos_t *phase_sets) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int *noisy_reg_reads = chunk->noisy_reg_to_reads[noisy_reg_i];
    hts_pos_t reg_beg = cr_start(noisy_regs, noisy_reg_i), reg_end = cr_end(noisy_regs, noisy_reg_i);
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        int read_i = noisy_reg_reads[i];
        chunk->haps[read_i] = haps[i];
        chunk->PS[read_i] = phase_sets[i];
    }
}

// currently, we call consensus sequences only for two cases:
// 1. >=1 phase set with >=1 fully-covered read for each haplotype
// 2. no phase set, >= 1 fully-covered read, and region length < T, de novo consensus calling using all reads
// otherwise: we skip the region (return 0 consensus)
int collect_noisy_reg_cons_seqs0(const call_var_opt_t *opt, bam_chunk_t *chunk, int noisy_reg_i, int *cons_lens, uint8_t **cons_seqs) {
    int n_noisy_reg_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    if (n_noisy_reg_reads <= 0) return 0;
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t reg_beg = cr_start(noisy_regs, noisy_reg_i), reg_end = cr_end(noisy_regs, noisy_reg_i);
    // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
    char **names = NULL; uint8_t **seqs = NULL; int *fully_covers = NULL, *lens = NULL, *haps = NULL; hts_pos_t *phase_sets = NULL;
    // collect reads in the noisy region
    collect_noisy_reads(chunk, noisy_reg_i, &lens, &seqs, &names, &fully_covers, &haps, &phase_sets);
    // hts_pos_t ps_with_full_hap_reads;
    // int n_full_cover_reads = collect_noisy_full_hap_reads(n_noisy_reg_reads, haps, phase_sets, fully_covers, &ps_with_full_hap_reads);
    hts_pos_t ps_with_both_haps = collect_phase_set_with_both_haps(n_noisy_reg_reads, haps, phase_sets, fully_covers, 2);
    int n_full_reads = 0;
    for (int i = 0; i < n_noisy_reg_reads; ++i) if (fully_covers[i] == 3) n_full_reads++;
    int n_cons = 0;
    // XXX only use fully-covered reads, including cliping reads (after re-align to backbone read)
    // two cases to call consensus sequences
    if (ps_with_both_haps > 0) { // call consensus sequences for each haplotype
        // fprintf(stderr, "PerHap calling region: %s:%ld-%ld %ld %d (all)\n", chunk->tname, reg_beg, reg_end, reg_end-reg_beg+1, n_noisy_reg_reads);
        // n_cons = collect_noisy_cons_seqs_with_ps_hap(opt, n_noisy_reg_reads, lens, seqs, names, haps, phase_sets, fully_covers, ps_with_full_hap_reads, cons_lens, cons_seqs);
        n_cons = collect_noisy_cons_seqs_with_ps_hap(opt, n_noisy_reg_reads, lens, seqs, names, haps, phase_sets, fully_covers, ps_with_both_haps, cons_lens, cons_seqs);
    // } else if (n_full_cover_reads > n_noisy_reg_reads * 0.75 && reg_end - reg_beg + 1 <= 10000) { // de novo consensus calling using all reads, up to 2 consensus sequences
    } else if (n_full_reads >= 10) {
        // fprintf(stderr, "De nove calling region: %s:%ld-%ld %ld %d (all)\n", chunk->tname, reg_beg, reg_end, reg_end-reg_beg+1, n_noisy_reg_reads);
        n_cons = collect_noisy_cons_seq_no_ps_hap(opt, reg_beg, n_noisy_reg_reads, lens, seqs, names, haps, phase_sets, fully_covers, cons_lens, cons_seqs);
    } else { // if (n_full_cover_reads <= 0 || reg_end - reg_beg + 1 > 5000) {
        fprintf(stderr, "Skipped region: %s:%ld-%ld %ld %d reads (%d full)\n", chunk->tname, reg_beg, reg_end, reg_end-reg_beg+1, n_noisy_reg_reads, n_full_reads);
    }

    // update phase-set and HAP
    // updated_ps_hap(chunk, noisy_reg_i, n_noisy_reg_reads, haps, phase_sets);
    for (int i = 0; i < n_noisy_reg_reads; ++i) free(seqs[i]);
    free(names); free(seqs); free(lens); free(fully_covers); free(haps); free(phase_sets);
    return n_cons;
}

int test_edlib(char *pattern, char *text) {
    // 
    int k = -1;
    const EdlibAlignConfig config = edlibNewAlignConfig(k, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0);
    EdlibAlignResult result = edlibAlign(pattern, strlen(pattern), text, strlen(text), config);
    fprintf(stderr, "Edlib-Alignment returns score %d\n",result.editDistance);
    edlibFreeAlignResult(result);
    return 0;
}

int test_abpoa(uint8_t **bseqs, int n_seqs, int *seq_lens) {
    // Init abpoa graph
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->out_msa = 1;
    abpt->out_cons = 1;
    abpoa_post_set_para(abpt);
    abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, bseqs, NULL, stdout);

    abpoa_cons_t *abc = ab->abc;
    for (int i = 0; i < abc->n_cons; ++i) {
        fprintf(stdout, ">Consensus_sequence");
        if (abc->n_cons > 1) {
            fprintf(stdout, "_%d ", i+1);
            for (int j = 0; j < abc->clu_n_seq[i]; ++j) { // output read ids for each cluster/group
                fprintf(stdout, "%d", abc->clu_read_ids[i][j]);
                if (j != abc->clu_n_seq[i]-1) fprintf(stdout, ",");
            }
        }
        fprintf(stdout, "\n");
        for (int j = 0; j < abc->cons_len[i]; ++j)
            fprintf(stdout, "%c", "ACGTN"[abc->cons_base[i][j]]);
        fprintf(stdout, "\n");
    }
    // Free
    abpoa_free(ab); abpoa_free_para(abpt);
    return 0;
}

int test_ksw2(char *pattern, char *text, const call_var_opt_t *opt) {
    uint8_t *p = (uint8_t*)malloc(strlen(pattern)*sizeof(uint8_t));
    uint8_t *t = (uint8_t*)malloc(strlen(text)*sizeof(uint8_t));
    for (int i = 0; i < strlen(pattern); ++i) p[i] = nst_nt4_table[(int)pattern[i]];
    for (int i = 0; i < strlen(text); ++i) t[i] = nst_nt4_table[(int)text[i]];
    uint32_t *cigar=NULL;
    int cigar_len = ksw2_aln(opt->gap_aln, p, strlen(pattern), t, strlen(text), opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, &cigar);
    if (cigar_len > 0) {
        for (int i = 0; i < cigar_len; ++i) {
            int op = cigar[i]&0xf, len = cigar[i]>>4;
            fprintf(stderr, "%d%c", len, BAM_CIGAR_STR[op]);
        } fprintf(stderr, "\n");
        free(cigar);
    }
    return 0;
}