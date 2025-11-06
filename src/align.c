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
#include "math_utils.h"
#include "align.h"
#include "kmer.h"

extern int LONGCALLD_VERBOSE;


// INS/DEL : left-most
// short   :    XXXXX                      [TSD]XXXXX
// long    :    XXXXX[TSD]----------[polyA][TSD]XXXXX
// or:
// short   :    XXXXX                      [TSD]XXXXX
// long    :    XXXXX[TSD][polyT]----------[TSD]XXXXX
// INS/DEL : right-most (not considered in this function)
// short   :    XXXXX[TSD]                      XXXXX
// long    :    XXXXX[TSD]----------[polyA][TSD]XXXXX
// or:
// short   :    XXXXX[TSD]                      XXXXX
// long    :    XXXXX[TSD][polyT]----------[TSD]XXXXX
int collect_te_info(const call_var_opt_t *opt, int var_type, uint8_t *gap_seq, uint8_t *flank_ref_seq, int gap_len, hts_pos_t gap_pos,
                    uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev) {
    *tsd_seq = NULL; *tsd_pos1 = -1; *tsd_pos2 = -1; *tsd_polya_len = -1; *te_seq_i = -1; *te_is_rev = 0;
    int min_tsd_len = opt->min_tsd_len, max_tsd_len = opt->max_tsd_len;
    int min_polya_len = opt->min_polya_len; float min_polya_ratio = opt->min_polya_ratio;
    // compare long_seq[gap_start_i:] with ref_seq[gap_end_i+1:]
    int tsd_len = 0, n_mis = 0, max_allow_mis = 1;
    for (int i = 0; i < gap_len; i++) {
        uint8_t base1 = gap_seq[i], base2 = flank_ref_seq[i];
        if (base1 == base2) tsd_len = i+1;
        else {
            n_mis++;
            if (n_mis > max_allow_mis) break;
        }
        if (tsd_len > max_tsd_len) break;
    }
    int has_tsd = 0, has_polya = 0, max_search_polya_len = 20;
    if (tsd_len >= min_tsd_len && tsd_len <= max_tsd_len) {
        has_tsd = 1;
        // look for polyA
        for (int polya_len = 0, polya = 0, i = gap_len-1; i >= 0; i--) {
            polya_len++;
            if (gap_seq[i] == 0) { polya++;
                if (polya_len >= min_polya_len && polya >= min_polya_ratio*polya_len) {
                    has_polya = 1;
                    *tsd_polya_len = polya_len; // update tsd_polya_len
                }
            } else if (polya_len > max_search_polya_len) break; // stop searching if too long
        }
        if (has_polya == 0) { // look for polyT
            for (int polyt_len = 0, polyt = 0, i = tsd_len; i < gap_len; i++) {
                polyt_len++;
                if (gap_seq[i] == 3) {
                    polyt++;
                    if (polyt_len >= min_polya_len && polyt >= min_polya_ratio*polyt_len) {
                        has_polya = 1;
                        *tsd_polya_len = -polyt_len; // update tsd_polya_len
                    }
                } else if (polyt_len > max_search_polya_len) break;
            }
        }
    }
    // check if INS/DEL is TE
    if (has_tsd && has_polya) {
        if (opt->n_te_seqs > 0) *te_seq_i = check_te_seq(opt, gap_seq, gap_len, te_is_rev);
        *tsd_seq = (uint8_t*)malloc(tsd_len*sizeof(uint8_t));
        for (int i = 0; i < tsd_len; i++) (*tsd_seq)[i] = flank_ref_seq[i];
        *tsd_pos1 = gap_pos; 
        if (var_type == BAM_CDEL) *tsd_pos2 = gap_pos+gap_len; else *tsd_pos2 = -1;
        return tsd_len;
    } else return 0;
}


// for somatic TE variant:
int collect_te_info_from_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *var) {
    if (var->checked_tsd) return var->tsd_len; // already checked
    var->checked_tsd = 1; // mark as checked
    if (!((var->var_type == BAM_CINS && var->alt_len >= opt->min_sv_len) ||
        (var->var_type == BAM_CDEL && var->ref_len >= opt->min_sv_len)))
        return 0;
    // check if have TSD & polyA/T
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    uint8_t *gap_seq, *flank_ref_seq; int gap_len;
    if (var->var_type == BAM_CINS) {
        gap_len = var->alt_len;
        gap_seq = (uint8_t*)malloc(var->alt_len*sizeof(uint8_t));
        flank_ref_seq = (uint8_t*)malloc(var->alt_len*sizeof(uint8_t));
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = var->alt_seq[i];
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i);
        }
    } else {
        gap_len = var->ref_len;
        gap_seq = (uint8_t*)malloc(var->ref_len*sizeof(uint8_t));
        flank_ref_seq = (uint8_t*)malloc(var->ref_len*sizeof(uint8_t));
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i);
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, var->pos+i+gap_len);
        }
    }
    int tsd_len = 0, tsd_polya_len=0, te_seq_i=-1, te_is_rev=0;
    uint8_t *tsd_seq=NULL; hts_pos_t tsd_pos1=-1, tsd_pos2=-1; 
    tsd_len = collect_te_info(opt, var->var_type, gap_seq, flank_ref_seq, gap_len, var->pos,
                              &tsd_seq, &tsd_pos1, &tsd_pos2, &tsd_polya_len, &te_seq_i, &te_is_rev);
    free(gap_seq); free(flank_ref_seq);
    if (tsd_len > 0) { // copy tsd info to var
        var->tsd_seq = tsd_seq;
        var->tsd_len = tsd_len;
        var->polya_len = tsd_polya_len;
        var->tsd_pos1 = tsd_pos1;
        var->tsd_pos2 = tsd_pos2;
        var->te_seq_i = te_seq_i;
        var->te_is_rev = te_is_rev;
        return tsd_len;
    } else {
        var->tsd_seq = NULL; var->tsd_len = 0; var->tsd_pos1 = -1; var->tsd_pos2 = -1; var->polya_len = 0;
        var->te_seq_i = -1; var->te_is_rev = 0;
        return 0;
    }
}

// for germline TE variant:
// input  : cons' aln-str, or raw read' aln-str (mosaic variant)
// return : 0: not have TSD & polyA, 1: have TSD & polyA
//        tsd_len: length of TSD
//        tsd_target_pos1/2: position of 1st/2nd TSD on target
int collect_te_info_from_cons(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t gap_ref_start, int msa_gap_start, int var_type, int gap_len, uint8_t *cons_msa_seq, 
                              uint8_t **tsd_seq, hts_pos_t *tsd_pos1, hts_pos_t *tsd_pos2, int *tsd_polya_len, int *te_seq_i, int *te_is_rev) {
    xassert(var_type == BAM_CINS || var_type == BAM_CDEL, "Var_type should be INS/DEL\n");
    // in case flanking sequence in ref_msa_seq is not long enough
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    uint8_t *gap_seq = (uint8_t*)malloc(gap_len*sizeof(uint8_t));
    uint8_t *flank_ref_seq = (uint8_t*)malloc(gap_len*sizeof(uint8_t));
    if (var_type == BAM_CINS) {
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = cons_msa_seq[msa_gap_start+i];
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i);
        }
    } else {
        for (int i = 0; i < gap_len; ++i) {
            gap_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i);
            flank_ref_seq[i] = get_bseq1(ref_seq, ref_beg, ref_end, gap_ref_start+i+gap_len);
        }
    }
    int tsd_len = 0;
    tsd_len = collect_te_info(opt, var_type, gap_seq, flank_ref_seq, gap_len, gap_ref_start,
                    tsd_seq, tsd_pos1, tsd_pos2, tsd_polya_len, te_seq_i, te_is_rev);
    free(gap_seq); free(flank_ref_seq);
    return tsd_len;
}

int* edlibAlignmentToXID(const unsigned char* const alignment, const int alignmentLength) {
    // Maps move code from alignment to char in cigar.
    //                      0    1    2    3
    // moveCodeToChar[] : {'=', 'I', 'D', 'X'};
    int *xid = (int*) calloc(4, sizeof(int));

    int i;
    for (i = 0; i < alignmentLength; i++) {
        if (alignment[i] == EDLIB_EDOP_MATCH) {
            xid[0]++; // equal
        } else if (alignment[i] == EDLIB_EDOP_MISMATCH) {
            xid[1]++; // mismatch
        } else if (alignment[i] == EDLIB_EDOP_INSERT) {
            xid[2]++; // insertion
        } else if (alignment[i] == EDLIB_EDOP_DELETE) {
            xid[3]++; // deletion
        } else {
            free(xid);
            return NULL; // invalid alignment code
        }
    }
    return xid;
}

// count number of mismatches and gaps from edlib alignment
int edlibAlignmentToXGAPS(const unsigned char* const alignment, const int alignmentLength) {
    if (LONGCALLD_VERBOSE >= 3) {
        fprintf(stderr, "Edlib alignment: ");
        for (int i = 0; i < alignmentLength; i++) {
            fprintf(stderr, "%d", alignment[i]);
        }
        fprintf(stderr, "\n");
    }
    int n_gaps = 0, n_mismatch = 0;
    for (int i = 0; i < alignmentLength; i++) {
            if (alignment[i] == EDLIB_EDOP_MATCH) {
            continue;
        } else if (alignment[i] == EDLIB_EDOP_MISMATCH) {
            n_mismatch++;
        } else if (alignment[i] == EDLIB_EDOP_INSERT || alignment[i] == EDLIB_EDOP_DELETE) {
            if (i == 0 || alignment[i - 1] != alignment[i]) n_gaps++;
        }
    }
    return n_mismatch + n_gaps;
}

int edlib_edit_distance(uint8_t *target, int tlen, uint8_t *query, int qlen) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_DISTANCE, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

int edlib_xgaps(uint8_t *target, int tlen, uint8_t *query, int qlen) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    int n_x_gaps = edlibAlignmentToXGAPS(result.alignment, result.alignmentLength);
    edlibFreeAlignResult(result);
    return n_x_gaps;
}

int edlib_end2end_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
            *n_eq = xid[0]; *n_xid = xid[1] + xid[2] + xid[3];
            free(xid);
        } else {
            *n_eq = -1; *n_xid = -1;
        }
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

// infix/HW alignment using edlib
int edlib_infix_aln(uint8_t *target, int tlen, uint8_t *query, int qlen, int *n_eq, int *n_xid) {
    EdlibAlignResult result = edlibAlign((const char*)query, qlen, (const char*)target, tlen,
                                         edlibNewAlignConfig(-1, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0));
    if (result.status != EDLIB_STATUS_OK) {
        edlibFreeAlignResult(result);
        return -1;
    }
    if (n_eq != NULL && n_xid != NULL) { // collect equal and XID bases
        int *xid = edlibAlignmentToXID(result.alignment, result.alignmentLength);
        if (xid != 0) {
            *n_eq = xid[0]; *n_xid = xid[1] + xid[2] + xid[3];
            free(xid);
        } else {
            *n_eq = -1; *n_xid = -1;
        }
    }
    int score = result.editDistance;
    edlibFreeAlignResult(result);
    return score;
}

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
    // assert(pattern_pos == pattern_length && text_pos == text_length);
    *_pattern_alg = pattern_alg; *_text_alg = text_alg;
    return alg_pos;
}

// return alignment score
int wfa_heuristic_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen, 
                      int a, int b, int q, int e, int q2, int e2, int *n_eq, int *n_xid) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.memory_mode = wavefront_memory_ultralow;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = b;
    attributes.affine2p_penalties.gap_opening1 = q;
    attributes.affine2p_penalties.gap_extension1 = e;
    attributes.affine2p_penalties.gap_opening2 = q2;
    attributes.affine2p_penalties.gap_extension2 = e2;
    attributes.alignment_scope = compute_alignment;
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

// cons vs ref: ultra-low mem + dual-gap + no heuristic to ensure accurate alignment
// full read/cons vs full read: ultra-low mem + affine-gap + wf-adaptive to ensure accurate alignment and speed up
// full read/cons vs partial read: high mem (default) + affine-gap + heuristic (xdrop/zdrop) to speed up
int wfa_end2end_aln(uint8_t *pattern, int plen, uint8_t *text, int tlen,
                    int gap_aln, int b, int q, int e, int q2, int e2, int heuristic, int affine_gap, // heuristic: 0: no, 1: default, 2: zdrop
                    uint32_t **cigar_buf, int *cigar_length, uint8_t **pattern_alg, uint8_t **text_alg, int *alg_length) {
    // double realtime0 = realtime();
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    if (heuristic != LONGCALLD_WFA_ZDROP) attributes.memory_mode = wavefront_memory_ultralow;
    if (affine_gap == LONGCALLD_WFA_AFFINE_2P) {
        attributes.distance_metric = gap_affine_2p;
        attributes.affine2p_penalties.match = 0; // -a;
        attributes.affine2p_penalties.mismatch = b;
        attributes.affine2p_penalties.gap_opening1 = q;
        attributes.affine2p_penalties.gap_extension1 = e;
        attributes.affine2p_penalties.gap_opening2 = q2;
        attributes.affine2p_penalties.gap_extension2 = e2;
    } else {
        attributes.distance_metric = gap_affine;
        attributes.affine_penalties.match = 0; // -a;
        attributes.affine_penalties.mismatch = b;
        attributes.affine_penalties.gap_opening = q;
        attributes.affine_penalties.gap_extension = e;
    }
    attributes.alignment_scope = compute_alignment;
    attributes.alignment_form.span = alignment_end2end;
    if (heuristic == LONGCALLD_WFA_NO_HEURISTIC) 
        attributes.heuristic.strategy = wf_heuristic_none;
    if (heuristic == LONGCALLD_WFA_ADAPTIVE) // default heuristic
        attributes.heuristic.strategy = wf_heuristic_wfadaptive;
    else if (heuristic == LONGCALLD_WFA_ZDROP) { // Zdrop
        attributes.heuristic.strategy = wf_heuristic_zdrop;
        attributes.heuristic.zdrop = MIN_OF_TWO(500, (int)(MIN_OF_TWO(plen, tlen) * 0.1));
        attributes.heuristic.steps_between_cutoffs = 100;
    }
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
    if (LONGCALLD_VERBOSE >= 3 && cigar != NULL && cigar->cigar_length > 0) {
        char *p_seq = (char*)malloc(plen+1); char *t_seq = (char*)malloc(tlen+1);
        for (int i = 0; i < plen; ++i) p_seq[i] = "ACGTN"[p[i]];
        for (int i = 0; i < tlen; ++i) t_seq[i] = "ACGTN"[t[i]];
        cigar_print_pretty(stderr, cigar, p_seq, plen, t_seq, tlen);
        free(p_seq); free(t_seq);
    }
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
    // fprintf(stderr, "%d vs %d Real time: %.3f sec. Score: %d\n", plen, tlen, realtime() - realtime0, cigar_score_gap_affine2p(cigar, &attributes.affine2p_penalties));
    wavefront_aligner_delete(wf_aligner); 
    if (gap_aln == LONGCALLD_GAP_LEFT_ALN) { free(p); free(t); }
    return 0;
}


int wfa_collect_diff_ins_seq(const call_var_opt_t *opt, uint8_t* large_seq, int large_len, uint8_t *small_seq, int small_len, uint8_t **diff_seq) {
    int diff_ins_len = 0;
    uint8_t *large_aln_seq = NULL, *small_aln_seq = NULL; int aln_len;
    wfa_end2end_aln(large_seq, large_len, small_seq, small_len,
                    opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, // no heuristic, affine-2p
                    NULL, NULL, &large_aln_seq, &small_aln_seq, &aln_len);
    // collect the largest insertion sequence
    int largest_ins_len = 0, largest_ins_pos = -1;
    for (int i = 0; i < aln_len; ++i) {
        if (small_aln_seq[i] == 5 && large_aln_seq[i] != 5) { // insertion in large_seq
            int ins_len = 0, j = i;
            while (j < aln_len && small_aln_seq[j] == 5 && large_aln_seq[j] != 5) {
                ins_len++; j++;
            }
            if (ins_len > largest_ins_len) {
                largest_ins_len = ins_len; largest_ins_pos = i;
            }
            i = j-1; // skip the inserted bases
        }
    }
    if (largest_ins_len > 0) {
        diff_ins_len = largest_ins_len;
        *diff_seq = (uint8_t*)malloc(largest_ins_len * sizeof(uint8_t));
        for (int i = 0; i < largest_ins_len; ++i) {
            (*diff_seq)[i] = large_aln_seq[largest_ins_pos+i];
        }
    }
    free(large_aln_seq);
    return diff_ins_len;
}

// trim query bases if longer than target
// update query_beg/query_end if query is shorter than target
void wfa_trim_aln_str(int full_cover, aln_str_t *aln_str) {
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover) || LONGCALLD_NOISY_IS_BOTH_COVER(full_cover)) return;

    if ((LONGCALLD_NOISY_IS_LEFT_COVER(full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(full_cover)) ||
        (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(full_cover))) {
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = 0; aln_str->query_end = aln_str->aln_len-1;
        return;
    }
    if (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) { // left-cover: collect all aln_str until the very last = opertion
        int target_end = -1, query_end = -1;
        for (int i = aln_str->aln_len-1; i >= 0; --i) {
            if (query_end == -1) {
                if (aln_str->query_aln[i] != 5 && aln_str->target_aln[i] == aln_str->query_aln[i]) {
                    query_end = i;
                }
            }
            if (target_end == -1) {
                if (aln_str->target_aln[i] != 5) {
                    target_end = i;
                }
            }
            if (target_end != -1 && query_end != -1) break;
        }
        if (query_end == -1) query_end = target_end; // no = operation, query_end is the last position
        assert(query_end <= target_end);
        aln_str->aln_len = target_end+1;
        aln_str->target_beg = 0; aln_str->target_end = target_end;
        aln_str->query_beg = 0; aln_str->query_end = query_end;
        for (int i = query_end+1; i < aln_str->aln_len; ++i) {
            aln_str->query_aln[i] = 5;
        }
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover)) { // right-cover: collect all aln_str until the very first = opertion
        int query_start = -1, target_start = -1;
        for (int i = 0; i < aln_str->aln_len; ++i) {
            if (query_start == -1) {
                if (aln_str->query_aln[i] != 5 && aln_str->target_aln[i] == aln_str->query_aln[i]) {
                    query_start = i;
                }
            }
            if (target_start == -1) {
                if (aln_str->target_aln[i] != 5) {
                    target_start = i;
                }
            }
            if (target_start != -1 && query_start != -1) break;
        }
        if (query_start == -1) query_start = target_start; // no = operation, query_start is the first position
        assert(query_start >= target_start);
        aln_str->aln_len = aln_str->aln_len - target_start;
        if (target_start != 0) {
            uint8_t *_target_aln = (uint8_t*)malloc(aln_str->aln_len * 2 * sizeof(uint8_t));
            for (int i = 0; i < aln_str->aln_len; ++i)
                _target_aln[i] = aln_str->target_aln[i+target_start];
            for (int i = 0; i < aln_str->aln_len; ++i)
                _target_aln[i+aln_str->aln_len] = aln_str->query_aln[i+target_start];
            free(aln_str->target_aln);
            aln_str->target_aln = _target_aln;
            aln_str->query_aln = _target_aln + aln_str->aln_len;
        }
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = query_start-target_start; aln_str->query_end = aln_str->aln_len-1;
        for (int i = 0; i < aln_str->query_beg; ++i) {
            aln_str->query_aln[i] = 5;
        }
    }
}

// for read with full_cover as 1 or 2, collect beg/end positions of read mapped to target
int wfa_collect_aln_str(const call_var_opt_t *opt, uint8_t *target, int tlen, uint8_t *query, int qlen, int full_cover, int heuristic, int affine_gap, aln_str_t *aln_str) {
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover)) return 0;
    aln_str->target_aln = 0; aln_str->query_aln = 0; aln_str->aln_len = 0;
    int gap_aln = opt->gap_aln, b = opt->mismatch, q = opt->gap_open1, e = opt->gap_ext1, q2 = opt->gap_open2, e2 = opt->gap_ext2;
    if (LONGCALLD_NOISY_IS_BOTH_COVER(full_cover)) {
        wfa_end2end_aln(target, tlen, query, qlen, gap_aln, b, q, e, q2, e2, heuristic, affine_gap,
                        NULL, NULL, &aln_str->target_aln, &aln_str->query_aln, &aln_str->aln_len);
        aln_str->target_beg = 0; aln_str->target_end = aln_str->aln_len-1;
        aln_str->query_beg = 0; aln_str->query_end = aln_str->aln_len-1;
    } else { // trim aln_str if query is longer than target; 
             // update query_beg/query_end if query is shorter than target
        int ext_direction = (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) ? LONGCALLD_EXT_ALN_LEFT_TO_RIGHT : LONGCALLD_EXT_ALN_RIGHT_TO_LEFT;
        // trim target/query to opt->LONGCALLD_PARTIAL_ALN_RATIO
        double ratio = opt->partial_aln_ratio;
        int _tlen = tlen, _qlen = qlen;
        int _t_start = 0, _q_start = 0;
        if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
            if (tlen > qlen*ratio) _tlen = qlen*ratio;
            else if (qlen > tlen*ratio) _qlen = tlen*ratio;
        } else { // Right to Left
            if (tlen > qlen*ratio) { _tlen = qlen*ratio; _t_start = tlen - _tlen; }
            else if (qlen > tlen*ratio) { _qlen = tlen*ratio; _q_start = qlen - _qlen; }
        }
        if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
            gap_aln = (gap_aln == LONGCALLD_GAP_LEFT_ALN) ? LONGCALLD_GAP_RIGHT_ALN : LONGCALLD_GAP_LEFT_ALN;
        }
        // partial alignment
        wfa_end2end_aln(target+_t_start, _tlen, query+_q_start, _qlen, 
                        gap_aln, b, q, e, q2, e2, LONGCALLD_WFA_ZDROP, LONGCALLD_WFA_AFFINE_1P,
                        NULL, NULL, &aln_str->target_aln, &aln_str->query_aln, &aln_str->aln_len);
        // do not trim for end-gaps
        wfa_trim_aln_str(full_cover, aln_str);
        aln_str->target_beg += _t_start; aln_str->target_end += _t_start;
        aln_str->query_beg += _q_start; aln_str->query_end += _q_start;
    }
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, ">target %d:%d-%d\n", tlen, aln_str->target_beg, aln_str->target_end);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->target_aln[i]]);
        fprintf(stderr, "\n>query %d:%d-%d\n", qlen, aln_str->query_beg, aln_str->query_end);
        for (int i = 0; i < aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[aln_str->query_aln[i]]);
        fprintf(stderr, "\n");
    }
    return 0;
}

int end2end_aln(const call_var_opt_t *opt, char *tseq, int tlen, uint8_t *qseq, int qlen, uint32_t **cigar_buf) {
    if (qlen <= 0 || tlen <= 0) return 0;
    // int min_len = MIN_OF_TWO(tlen, qlen)
    int max_len = MAX_OF_TWO(tlen, qlen);
    // int delta_len = MAX_OF_TWO(1, max_len - min_len);

    uint8_t *tseq2 = (uint8_t*)malloc(max_len);
    for (int i = 0; i < tlen; ++i) tseq2[i] = nst_nt4_table[(uint8_t)tseq[i]];
    // use wfa if the length difference is small
    // if (max_len * delta_len + delta_len * delta_len < max_len * min_len) {
    int cigar_len = 0;
    wfa_end2end_aln(tseq2, tlen, qseq, qlen, opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P,
                    cigar_buf, &cigar_len, NULL, NULL, NULL);
    // } else { // if (max_len - min_len > 1000) { // use ksw2 if the length difference is large
        // cigar_len = ksw2_aln(opt->gap_aln, tseq2, tlen, qseq, qlen, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, cigar_buf);
    // }
    free(tseq2);
    return cigar_len;
}

void collect_aln_beg_end(uint32_t *cigar_buf, int cigar_len, int ext_direction, int ref_len, int *ref_beg, int *ref_end, int read_len, int *read_beg, int *read_end) {
    *ref_beg = 1; *read_beg = 1; *ref_end = ref_len, *read_end = read_len;
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
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
        int tmp_ref_beg = ref_len+1, tmp_read_beg = read_len+1;
        for (int i=cigar_len-1; i>=0; i--) {
            int op = cigar_buf[i] & 0xf; int len = cigar_buf[i] >> 4;
            if (op == BAM_CEQUAL || op == BAM_CMATCH) {
                tmp_ref_beg -= len; tmp_read_beg -= len;
                *ref_beg = tmp_ref_beg; *read_beg = tmp_read_beg;
            } else if (op == BAM_CDIFF) {
                tmp_ref_beg -= len; tmp_read_beg -= len;
            } else if (op == BAM_CDEL) {
                tmp_ref_beg -= len;
            } else if (op == BAM_CINS) {
                tmp_read_beg -= len;
            } else continue;
        }
    }
}

// calculate the beg/end positions of query mapped to target using WFA
// limit to (partial_aln_ratio*short_len vs short_len)
int cal_wfa_partial_aln_beg_end(int ext_direction, const call_var_opt_t *opt, uint8_t *_target, int _tlen, uint8_t *_query, int _qlen, int *target_beg, int *target_end, int *query_beg, int *query_end) {
    int gap_aln = opt->gap_aln, b = opt->mismatch, q = opt->gap_open1, e = opt->gap_ext1, q2 = opt->gap_open2, e2 = opt->gap_ext2;
    double ratio = opt->partial_aln_ratio;
    int tlen = _tlen, qlen = _qlen;
    uint8_t *target = _target, *query = _query;
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        if (_tlen > _qlen * ratio) tlen = (int)(_qlen * ratio);
        else if (_qlen > _tlen * ratio) qlen = (int)(_tlen * ratio);
    } else if (ext_direction == LONGCALLD_EXT_ALN_RIGHT_TO_LEFT) {
        if (_tlen > _qlen * ratio) {
            target = _target + _tlen - (int)(_qlen * ratio);
            tlen = (int)(_qlen * ratio);
        } else if (_qlen > _tlen * ratio) {
            query = _query + _qlen - (int)(_tlen * ratio);
            qlen = (int)(_tlen * ratio);
        }
    }
    uint32_t *cigar_buf = NULL;
    int ret = 1, cigar_len;
    // left to right
    if (ext_direction == LONGCALLD_EXT_ALN_LEFT_TO_RIGHT) {
        gap_aln = (gap_aln == LONGCALLD_GAP_RIGHT_ALN) ? LONGCALLD_GAP_LEFT_ALN : LONGCALLD_GAP_RIGHT_ALN;
    }
    // fprintf(stderr, "WFA partial aln: %d vs %d (Raw: %d vs %d)\n", tlen, qlen, _tlen, _qlen);
    wfa_end2end_aln(target, tlen, query, qlen, gap_aln, b, q, e, q2, e2, LONGCALLD_WFA_ZDROP, LONGCALLD_WFA_AFFINE_1P, &cigar_buf, &cigar_len, NULL, NULL, NULL);
    if (cigar_len == 0) ret = 0;
    else collect_aln_beg_end(cigar_buf, cigar_len, ext_direction, _tlen, target_beg, target_end, _qlen, query_beg, query_end);
    if (cigar_buf != NULL) free(cigar_buf);
    return ret;
}

int collect_partial_aln_beg_end(const call_var_opt_t *opt, int sampling_reads,
                                uint8_t *target, int tlen, int target_full_cover, uint8_t *query, int qlen, int query_full_cover, 
                                int *target_beg, int *target_end, int *query_beg, int *query_end) {
    *target_beg = 1, *target_end = tlen, *query_beg = 1, *query_end = qlen;
    int ret = 1;
    assert(LONGCALLD_NOISY_IS_BOTH_COVER(target_full_cover) != 0);
    if (LONGCALLD_NOISY_IS_BOTH_COVER(target_full_cover)) { // ref is full cover
        if (LONGCALLD_NOISY_IS_BOTH_COVER(query_full_cover) || 
           (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(query_full_cover)) ||
           (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(query_full_cover))) {
            if (sampling_reads) { // skip full-cover reads with >5% ED in sampling mode
                // int edit_dis = edlib_edit_distance(target, tlen, query, qlen);
                // if (edit_dis > MIN_OF_TWO(tlen, qlen) * 0.05) {
                int x_gaps = edlib_xgaps(target, tlen, query, qlen);
                if (x_gaps > MIN_OF_TWO(tlen, qlen) * 0.05) {
                    // if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "Skipped in POA due to high edit distance (%d vs %d): %d\n", tlen, qlen, edit_dis);
                    if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "Skipped in POA due to high # of X-gaps (%d vs %d): %d\n", tlen, qlen, x_gaps);
                    return 0;
                }
            } 
            return 1;
        } else {
            if (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT_TO_RIGHT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
            } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover)) {
                ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT_TO_LEFT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
            }
        }
    } else if (LONGCALLD_NOISY_IS_LEFT_COVER(target_full_cover)) { // ref is left-cover
        if (LONGCALLD_NOISY_IS_LEFT_COVER(query_full_cover) == 0) _err_error_exit("Target is left-cover but read is not left-cover\n");
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_LEFT_TO_RIGHT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(target_full_cover)) { // ref is right-cover
        if (LONGCALLD_NOISY_IS_RIGHT_COVER(query_full_cover) == 0) _err_error_exit("Ref is right-cover but read is not right-cover\n");
        ret = cal_wfa_partial_aln_beg_end(LONGCALLD_EXT_ALN_RIGHT_TO_LEFT, opt, target, tlen, query, qlen, target_beg, target_end, query_beg, query_end);
    } else _err_error_exit("Target is not left or right-cover\n");
    return ret;
}

int collect_exc_beg_end(const call_var_opt_t *opt, abpoa_t *ab, abpoa_para_t *abpt, int sampling_reads, int n_reads, char **names, uint8_t **read_seqs, 
                        int *read_lens, int *read_full_cover, int *exc_begs, int *exc_ends, int *seq_beg_cuts, int *seq_end_cuts) {
    for (int read_i = 1; read_i < n_reads; ++read_i) {
        int ref_beg, ref_end, read_beg, read_end, beg_id, end_id;
        if (collect_partial_aln_beg_end(opt, sampling_reads, read_seqs[0], read_lens[0], read_full_cover[0], read_seqs[read_i], read_lens[read_i], read_full_cover[read_i], &ref_beg, &ref_end, &read_beg, &read_end) == 0) {
            exc_begs[read_i] = -1, exc_ends[read_i] = -1;
            continue;
        }
        beg_id = ref_beg+1, end_id = ref_end+1;
        seq_beg_cuts[read_i] = read_beg - 1, seq_end_cuts[read_i] = read_lens[read_i] - read_end;
        abpoa_subgraph_nodes(ab, abpt, beg_id, end_id, exc_begs+read_i, exc_ends+read_i);
    }
    return 0;
}

int abpoa_partial_aln_msa_cons(const call_var_opt_t *opt, abpoa_t *ab, int sampling_reads, int n_reads, int *read_ids, uint8_t **read_seqs, uint8_t **read_quals, int *read_lens, int *read_full_cover, char **names,
                               int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids, int *msa_seq_lens, uint8_t **msa_seqs) {
    // abpoa_t *ab = abpoa_init();
    int needs_free_ab = 0;
    if (ab == NULL) {
        ab = abpoa_init(); needs_free_ab = 1;
    }
    abpoa_para_t *abpt = abpoa_init_para();
    if (opt->out_somatic) abpt->wf = 0.01; // more accurate for somatic variant calling
    else abpt->wf = 0.001; // limit memory usage for long sequences
    abpt->cons_algrm = ABPOA_MF;
    abpt->sub_aln = 1;
    abpt->inc_path_score = 1;
    if (cons_lens != NULL && cons_seqs != NULL) abpt->out_cons = 1; else abpt->out_cons = 0;
    if (msa_seq_lens != NULL && msa_seqs != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    abpt->max_n_cons = max_n_cons; abpt->min_freq = opt->min_af;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;

    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpoa_post_set_para(abpt);
    ab->abs->n_seq = n_reads;
    int *exc_begs = (int*)malloc(n_reads * sizeof(int)), *exc_ends = (int*)malloc(n_reads * sizeof(int));
    int *seq_beg_cuts = (int*)malloc(n_reads * sizeof(int)), *seq_end_cuts = (int*)malloc(n_reads * sizeof(int));
    for (int i = 0; i < n_reads; ++i) {
        exc_begs[i] = 0; exc_ends[i] = 1, seq_beg_cuts[i] = 0, seq_end_cuts[i] = 0;
    }
    for (int i = 0; i < n_reads; ++i) {
        if (LONGCALLD_VERBOSE >= 3) {
            fprintf(stderr, ">%s %d %d\n", names[i], read_lens[i], read_full_cover[i]);
            for (int j = 0; j < read_lens[i]; ++j) {
                fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
            } fprintf(stderr, "\n");
        }
        if (exc_begs[i] < 0 || exc_ends[i] < 0) continue;
        abpoa_res_t res; res.graph_cigar = 0, res.n_cigar = 0;
        if (LONGCALLD_VERBOSE >= 3) fprintf(stderr, "ExcBeg: %d, ExcEnd: %d, SeqBegCut: %d, SeqEndCut: %d, FullCover: %d\n", exc_begs[i], exc_ends[i], seq_beg_cuts[i], seq_end_cuts[i], read_full_cover[i]);
        abpoa_align_sequence_to_subgraph(ab, abpt, exc_begs[i], exc_ends[i], read_seqs[i]+seq_beg_cuts[i], read_lens[i]-seq_beg_cuts[i]-seq_end_cuts[i], &res);
        abpoa_add_subgraph_alignment(ab, abpt, exc_begs[i], exc_ends[i], read_seqs[i]+seq_beg_cuts[i], NULL, read_lens[i]-seq_beg_cuts[i]-seq_end_cuts[i], NULL, res, i, n_reads, 0);
        // collect exc_beg/exc_end for all reads with i>1
        if (i == 0) collect_exc_beg_end(opt, ab, abpt, sampling_reads, n_reads, names, read_seqs, read_lens, read_full_cover, exc_begs, exc_ends, seq_beg_cuts, seq_end_cuts);
        if (res.n_cigar) free(res.graph_cigar);
    }
    free(exc_begs); free(exc_ends); free(seq_beg_cuts); free(seq_end_cuts);
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
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Unable to call consensus: %s\n", names[0]);
        }
    }

    // msa bases, include consensus sequences
    if (msa_seq_lens != NULL && msa_seqs != NULL) {
        *msa_seq_lens = abc->msa_len;
        for (int i = 0; i < abc->n_seq+n_cons; ++i) {
            msa_seqs[i] = (uint8_t*)malloc(abc->msa_len * sizeof(uint8_t));
            for (int j = 0; j < abc->msa_len; ++j) {
                msa_seqs[i][j] = abc->msa_base[i][j];
            }
        }
    }
    abpoa_free_para(abpt); if (needs_free_ab) abpoa_free(ab);
    return n_cons;
 }

// TODO: for large deletions
// int two_sides_abpoa_partial_aln_msa_cons(const call_var_opt_t *opt, abpoa_t *ab, int n_reads, int *read_ids, uint8_t **read_seqs, uint8_t **read_quals, int *read_lens, int *read_full_cover, char **names,
//                                          int max_n_cons, int *cons_lens, uint8_t **cons_seqs, int *clu_n_seqs, int **clu_read_ids) {
//     int n_cons = 0;
//     // all left-cover: manually set full-cover of first read as 3, collect consensus
//     int n_left_cover = 0;
//     // all right-cover: manually set full-cover of first read as 3, collect consensus
//     int n_right_cover = 0;
//     // merge consensus
//     return n_cons;
//  }

  // XXX limit abpoa memory usage, avoid memory allocation failure
 int abpoa_aln_msa_cons(const call_var_opt_t *opt, int n_reads, int *read_ids, uint8_t **read_seqs, int *read_lens, int max_n_cons,
                        int *cons_lens, uint8_t **cons_seqs,
                        int *clu_n_seqs, int **clu_read_ids, int *msa_seq_len, uint8_t ***msa_seq) {
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    if (opt->out_somatic) abpt->wf = 0.01; // XXX for more accurate somatic variant calling
    else abpt->wf = 0.001; // limit memory usage for long sequences
    abpt->inc_path_score = 1;
    abpt->out_msa = 0; abpt->out_cons = 1;
    if (msa_seq != NULL && msa_seq_len != NULL) abpt->out_msa = 1; else abpt->out_msa = 0;
    if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->max_n_cons = max_n_cons; abpt->min_freq = opt->min_af;
    abpt->match = opt->match; abpt->mismatch = opt->mismatch;
    abpt->gap_open1 = opt->gap_open1; abpt->gap_ext1 = opt->gap_ext1;
    abpt->gap_open2 = opt->gap_open2; abpt->gap_ext2 = opt->gap_ext2;
    abpoa_post_set_para(abpt);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA (max %d cons, min_freq: %.2f): %d\n", max_n_cons, abpt->min_freq, n_reads);
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
        _err_error_exit("No Consensus.\n");
    }
    // msa bases
    if (msa_seq != NULL && msa_seq_len != NULL) {
        // split msa and cons into 2 clusers
        size_t msize = abc->msa_len * sizeof(uint8_t);
        for (int i = 0; i < n_cons; ++i) {
            msa_seq_len[i] = abc->msa_len;
            for (int j = 0; j < abc->clu_n_seq[i]; ++j) {
                msa_seq[i][j] = (uint8_t*)malloc(msize);
                for (int k = 0; k < abc->msa_len; ++k)
                    msa_seq[i][j][k] = abc->msa_base[abc->clu_read_ids[i][j]][k];
            }
            msa_seq[i][abc->clu_n_seq[i]] = (uint8_t*)malloc(msize);
            for (int j = 0; j < abc->msa_len; ++j)
                msa_seq[i][abc->clu_n_seq[i]][j] = abc->msa_base[abc->n_seq+i][j];
        }
    }
    abpoa_free_para(abpt); abpoa_free(ab); 
    return n_cons;
}

int full_cover_cmp(int cover1, int cover2) {
    if (cover1 == cover2) return 0;
    if (LONGCALLD_NOISY_IS_BOTH_COVER(cover1)) return 1;
    else if (LONGCALLD_NOISY_IS_BOTH_COVER(cover2)) return -1;
    if (LONGCALLD_NOISY_IS_LEFT_COVER(cover1) && LONGCALLD_NOISY_IS_LEFT_COVER(cover2)) return 0;
    if (LONGCALLD_NOISY_IS_RIGHT_COVER(cover1) && LONGCALLD_NOISY_IS_RIGHT_COVER(cover2)) return 0;
    return cover1-cover2;
}

// sort reads in noisy regions based on full-cover, [error-rate], read-length; error-rate is optional
void sort_noisy_region_reads(int n_reads, int *read_ids, int *read_lens, uint8_t **read_seqs, uint8_t **read_quals, uint8_t *strands, int *full_covers, char **names,
                             int *haps, hts_pos_t *phase_sets, int use_error_rate) {
    double *read_error_rates = NULL;
    if (use_error_rate) {
        read_error_rates = (double*)malloc(n_reads * sizeof(double));
        for (int i = 0; i < n_reads; ++i)
            read_error_rates[i] = calc_read_error_rate(read_lens[i], read_quals[i]);
    }
    for (int i = 0; i < n_reads-1; ++i) {
        for (int j = i+1; j < n_reads; ++j) {
            int cover_cmp = full_cover_cmp(full_covers[i], full_covers[j]);
            if (cover_cmp < 0 || 
               (cover_cmp == 0 && use_error_rate && read_error_rates[i] > read_error_rates[j]) ||
               (cover_cmp == 0 && ((use_error_rate && read_error_rates[i] == read_error_rates[j]) || !use_error_rate) && read_lens[i] < read_lens[j])) {
            // if (full_covers[i] < full_covers[j] || (full_covers[i] == full_covers[j] && read_lens[i] < read_lens[j])) {
                int tmp_full_cover = full_covers[i]; full_covers[i] = full_covers[j]; full_covers[j] = tmp_full_cover;
                int tmp_id = read_ids[i]; read_ids[i] = read_ids[j]; read_ids[j] = tmp_id;
                int tmp_len = read_lens[i]; read_lens[i] = read_lens[j]; read_lens[j] = tmp_len;
                uint8_t *tmp_seq = read_seqs[i]; read_seqs[i] = read_seqs[j]; read_seqs[j] = tmp_seq;
                uint8_t *tmp_qual = read_quals[i]; read_quals[i] = read_quals[j]; read_quals[j] = tmp_qual;
                uint8_t tmp_strand = strands[i]; strands[i] = strands[j]; strands[j] = tmp_strand;
                char *tmp_name = names[i]; names[i] = names[j]; names[j] = tmp_name;
                int tmp_hap = haps[i]; haps[i] = haps[j]; haps[j] = tmp_hap;
                hts_pos_t tmp_ps = phase_sets[i]; phase_sets[i] = phase_sets[j]; phase_sets[j] = tmp_ps;
                double tmp_error_rate = 0.0;
                if (use_error_rate) {
                    tmp_error_rate = read_error_rates[i]; read_error_rates[i] = read_error_rates[j]; read_error_rates[j] = tmp_error_rate;
                }
            }
        }
    }
    if (use_error_rate) free(read_error_rates);
}

static int add_phase_set(hts_pos_t ps, hts_pos_t *uniq_phase_sets, int *n_uniq_phase_sets) {
    int i;
    for (i = 0; i < *n_uniq_phase_sets; ++i) {
        if (uniq_phase_sets[i] == ps) return i;
    }
    uniq_phase_sets[i] = ps;
    (*n_uniq_phase_sets)++;
    return i;
}

// hp_start/end: flank_base
int is_homopolymer(uint8_t *seq, int seq_len, int flank_len, int *hp_start, int *hp_end, int *hp_len) {
    if (seq_len < 2*flank_len || seq_len > 2*flank_len + 50) return 0; // skip too short or too long cons
    int min_hp_len = 5; *hp_start = -1, *hp_end = -1, *hp_len = 0;
    for (int i = flank_len-1; i < seq_len-flank_len+1; ++i) {
        if (seq[i] == seq[i-1]) {
            if (*hp_start == -1) *hp_start = i-2;
            (*hp_len)++;
        } else {
            if (*hp_len >= min_hp_len) {
                *hp_end = i;
                return 1;
            } else {
                *hp_start = -1; *hp_len = 0;
            }
        }
    }
    if (*hp_len >= min_hp_len) {
        *hp_end = *hp_start + *hp_len + 1;
        return 1;
    }
    return 0;
}

// construct cons vs read alignment string based on abPOA's msa result
// read: A-CCGGT--A
// cons: ACCC-GT--A
// output:
// read: A-CCGGTA
// cons: ACCC-GTA
int make_cons_read_aln_str(const call_var_opt_t *opt, uint8_t *cons_str, uint8_t *read_str, int msa_len, int full_cover, aln_str_t *cons_read_aln_str) {
    int aln_len = 0;
    cons_read_aln_str->aln_len = 0;
    cons_read_aln_str->target_aln = (uint8_t*)malloc(msa_len * 2 * sizeof(uint8_t));
    cons_read_aln_str->query_aln = cons_read_aln_str->target_aln + msa_len;
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] != 5 || cons_str[i] != 5) { // not both gaps
            cons_read_aln_str->target_aln[aln_len] = cons_str[i];
            cons_read_aln_str->query_aln[aln_len] = read_str[i];
            aln_len++;
        }
    }
    cons_read_aln_str->aln_len = aln_len;
    cons_read_aln_str->target_beg = 0; cons_read_aln_str->target_end = aln_len - 1;
    cons_read_aln_str->query_beg = 0; cons_read_aln_str->query_end = aln_len - 1;
    wfa_trim_aln_str(full_cover, cons_read_aln_str);

    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, ">target %d-%d\n", cons_read_aln_str->target_beg, cons_read_aln_str->target_end);
        for (int i = 0; i < cons_read_aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->target_aln[i]]);
        fprintf(stderr, "\n>query %d-%d\n", cons_read_aln_str->query_beg, cons_read_aln_str->query_end);
        for (int i = 0; i < cons_read_aln_str->aln_len; ++i) fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[i]]);
        fprintf(stderr, "\n");
    }
    return aln_len;
}

int make_ref_read_aln_str(const call_var_opt_t *opt, aln_str_t *ref_cons_aln_str, aln_str_t *cons_read_aln_str, aln_str_t *ref_read_aln_str) {
    int aln_len = 0;
    int max_msa_len = ref_cons_aln_str->aln_len + cons_read_aln_str->aln_len;
    ref_read_aln_str->aln_len = 0;
    ref_read_aln_str->target_aln = (uint8_t*)malloc(max_msa_len * 2 * sizeof(uint8_t));
    ref_read_aln_str->query_aln = ref_read_aln_str->target_aln + max_msa_len;
    int i = 0, j = 0;
    while (i < ref_cons_aln_str->aln_len && j < cons_read_aln_str->aln_len) {
        if (ref_cons_aln_str->query_aln[i] == 5 && cons_read_aln_str->target_aln[j] == 5) { // both are gaps: extract the seqs for ref and read, until the non-gap base of cons, do WFA
            int ref_del_len = 1, read_del_len = 1;
            while (i+ref_del_len < ref_cons_aln_str->aln_len && ref_cons_aln_str->query_aln[i+ref_del_len] == 5) ref_del_len++;
            while (j+read_del_len < cons_read_aln_str->aln_len && cons_read_aln_str->target_aln[j+read_del_len] == 5) read_del_len++;
            uint8_t *ref_aln=0, *read_aln=0; int del_aln_len;
            // print ref_cons_aln_str->target_aln[i] and cons_read_aln_str->query_aln[j] for debugging
            if (LONGCALLD_VERBOSE >= 3) {
                fprintf(stderr, "Ref-Cons Gap: %d, Cons-Read Gap: %d\n", ref_del_len, read_del_len);
                fprintf(stderr, "Ref-Cons Target Aln: ");
                for (int k = i; k < i+ref_del_len; ++k) {
                    fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->target_aln[k]]);
                } fprintf(stderr, "\n");
                fprintf(stderr, "Cons-Read Query Aln: ");
                for (int k = j; k < j+read_del_len; ++k) {
                    fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[k]]);
                } fprintf(stderr, "\n");
            }
            wfa_end2end_aln(ref_cons_aln_str->target_aln+i, ref_del_len, cons_read_aln_str->query_aln+j, read_del_len, 
                            opt->gap_aln, opt->mismatch, opt->gap_open1, opt->gap_ext1, opt->gap_open2, opt->gap_ext2, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, // no heuristic, affine-2p
                            NULL, NULL, &ref_aln, &read_aln, &del_aln_len);
            for (int k = 0; k < del_aln_len; ++k) {
                ref_read_aln_str->target_aln[aln_len] = ref_aln[k];
                ref_read_aln_str->query_aln[aln_len] = read_aln[k];
                aln_len++;
            }
            i+=ref_del_len; j+=read_del_len;
            free(ref_aln);
        } else if (ref_cons_aln_str->query_aln[i] != 5 && cons_read_aln_str->target_aln[j] != 5) { // both are not gaps: match or mismatch
            ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
            ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
            aln_len++; i++; j++;
        } else if (ref_cons_aln_str->query_aln[i] == 5 && cons_read_aln_str->target_aln[j] != 5) { // ref-cons is gap, cons-read is not gap
            ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
            ref_read_aln_str->query_aln[aln_len] = 5;
            aln_len++; i++;
        } else if (ref_cons_aln_str->query_aln[i] != 5 && cons_read_aln_str->target_aln[j] == 5) { // ref-cons is not gap, cons-read is gap
            ref_read_aln_str->target_aln[aln_len] = 5;
            ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
            aln_len++; j++;
        } else {
            fprintf(stderr, "Error: ref_cons_aln_str->query_aln[%d]: %d cons_read_aln_str->target_aln[%d]: %d\n", i, ref_cons_aln_str->query_aln[i], j, cons_read_aln_str->target_aln[j]);
            exit(1);
        }
    }
    while (i < ref_cons_aln_str->aln_len) {
        ref_read_aln_str->target_aln[aln_len] = ref_cons_aln_str->target_aln[i];
        ref_read_aln_str->query_aln[aln_len] = 5;
        aln_len++; i++;
    }
    while (j < cons_read_aln_str->aln_len) {
        ref_read_aln_str->target_aln[aln_len] = 5;
        ref_read_aln_str->query_aln[aln_len] = cons_read_aln_str->query_aln[j];
        aln_len++; j++;
    }
    ref_read_aln_str->aln_len = aln_len;
    ref_read_aln_str->target_beg = ref_read_aln_str->target_end = ref_read_aln_str->query_beg = ref_read_aln_str->query_end = -1; // unset
    if (LONGCALLD_VERBOSE >= 3) { // print ref-read aln-str
        fprintf(stderr, "Ref-Cons AlnStr: len: %d\n", ref_cons_aln_str->aln_len);
        for (int k = 0; k < ref_cons_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < ref_cons_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_cons_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");

        fprintf(stderr, "Cons-Read AlnStr: len: %d\n", cons_read_aln_str->aln_len);
        for (int k = 0; k < cons_read_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < cons_read_aln_str->aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[cons_read_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");

        fprintf(stderr, "Ref-Read AlnStr: len: %d\n", aln_len);
        for (int k = 0; k < aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_read_aln_str->target_aln[k]]);
        } fprintf(stderr, "\n");
        for (int k = 0; k < aln_len; ++k) {
            fprintf(stderr, "%c", "ACGTN-"[ref_read_aln_str->query_aln[k]]);
        } fprintf(stderr, "\n");
    }
    return aln_len;
}

int wfa_collect_noisy_aln_str_no_ps_hap(const call_var_opt_t *opt, int n_reads, int *read_ids, int *lens, uint8_t **seqs, char **qnames, int *fully_covers,
                                        uint8_t *ref_seq, int ref_seq_len, int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, int collect_ref_read_aln_str) {
    int *full_read_ids = (int*)malloc((n_reads+2) * sizeof(int));
    int *full_read_lens = (int*)malloc((n_reads+2) * sizeof(int));
    uint8_t **full_read_seqs = (uint8_t**)malloc((n_reads+2) * sizeof(uint8_t*));
    char **full_read_names = (char**)malloc((n_reads+2) * sizeof(char*));
    int *full_fully_covers = (int*)malloc((n_reads+2) * sizeof(int));
    int n_full_reads = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (lens[i] <= 0 || LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
        full_fully_covers[n_full_reads] = fully_covers[i];
        full_read_ids[n_full_reads] = i;
        full_read_lens[n_full_reads] = lens[i];
        full_read_seqs[n_full_reads] = seqs[i];
        full_read_names[n_full_reads] = qnames[i];
        n_full_reads++;
    }
    // abpoa
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    int *msa_seq_lens = (int*)malloc(2 * sizeof(int)); uint8_t ***msa_seqs = (uint8_t***)malloc(2 * sizeof(uint8_t**));
    for (int i = 0; i < 2; ++i) {
        cons_seqs[i] = NULL; cons_lens[i] = 0;
        msa_seqs[i] = (uint8_t**)malloc((n_reads+1) * sizeof(uint8_t*));
        for (int j = 0; j < n_reads+1; ++j) msa_seqs[i][j] = NULL;
        msa_seq_lens[i] = 0;
    }
    int n_cons = 0;
    if (n_full_reads == 0) goto collect_noisy_msa_cons_no_ps_hap_end;
    else {
        if (full_read_lens[0] >= opt->max_noisy_reg_len) goto collect_noisy_msa_cons_no_ps_hap_end;
    }

    n_cons = abpoa_aln_msa_cons(opt, n_full_reads, full_read_ids, full_read_seqs, full_read_lens, 2,
                                cons_lens, cons_seqs, clu_n_seqs, clu_read_ids, msa_seq_lens, msa_seqs);

    // re-do POA with ref_seq and cons
    for (int i = 0; i < n_cons; ++i) {
        aln_str_t *clu_aln_str = aln_strs[i];
        wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[i], cons_lens[i], LONGCALLD_NOISY_BOTH_COVER, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
        n_full_reads = 0;
        for (int j = 0; j < clu_n_seqs[i]; ++j) {
            int read_i = clu_read_ids[i][j];
            int read_id = read_ids[read_i];
            clu_read_ids[i][j] = read_id;
            // cons vs read
            // wfa_collect_aln_str(opt, cons_seqs[i], cons_lens[i], seqs[read_i], lens[read_i], fully_covers[read_i], heuristic, affine_2p, LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads));
            make_cons_read_aln_str(opt, msa_seqs[i][clu_n_seqs[i]], msa_seqs[i][j], msa_seq_lens[i], fully_covers[i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads));
            if (collect_ref_read_aln_str)
                make_ref_read_aln_str(opt, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str), LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_full_reads), LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_full_reads));
            n_full_reads++;
        }
    } 
collect_noisy_msa_cons_no_ps_hap_end:
    free(full_read_lens); free(full_read_seqs); free(full_read_names); free(full_fully_covers); free(full_read_ids);
    for (int i = 0; i < 2; ++i) {
        for (int j = 0; j < n_reads+1; ++j) {
            if (msa_seqs[i][j] != NULL) free(msa_seqs[i][j]);
        } free(msa_seqs[i]);
    } free(msa_seqs); free(msa_seq_lens);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
    } free(cons_lens); free(cons_seqs);
    return n_cons;
}

int ps_has_enough_reads(int ps_i, int **phase_set_to_hap_full_read_count, int **phase_set_to_hap_left_read_count, int **phase_set_to_hap_right_read_count,
                        int min_hap_full_read_count, int min_hap_all_read_count, int min_hap_partial_read_count) {
    if (ps_i < 0) return 0; // invalid phase set
    for (int hap=0; hap < 2; ++hap) {
        if ((phase_set_to_hap_full_read_count[ps_i][hap] < min_hap_full_read_count || phase_set_to_hap_left_read_count[ps_i][hap] + phase_set_to_hap_right_read_count[ps_i][hap] < min_hap_all_read_count)
          && (phase_set_to_hap_left_read_count[ps_i][hap] < min_hap_partial_read_count || phase_set_to_hap_right_read_count[ps_i][hap] < min_hap_partial_read_count))
          return 0;
    }
    return 1;
}

// collect phase set with both haps having >= min_minor_hap_read_count reads
hts_pos_t collect_phase_set_with_both_haps(int n_reads, int *read_haps, int *read_lens, hts_pos_t *phase_sets, int *fully_covers, int min_minor_hap_full_read_count, int min_minor_hap_all_read_count) {
    int n_uniq_phase_sets = 0, phase_set_i = 0;
    hts_pos_t *uniq_phase_sets = (hts_pos_t*)calloc(n_reads, sizeof(hts_pos_t));
    int **phase_set_to_hap_full_read_count = (int**)malloc(n_reads * sizeof(int*));
    int **phase_set_to_hap_all_read_count = (int**)malloc(n_reads * sizeof(int*)); // including clipped reads
    int **phase_set_to_hap_min_full_read_len = (int**)malloc(n_reads * sizeof(int*));
    for (int i = 0; i < n_reads; ++i) {
        phase_set_to_hap_full_read_count[i] = (int*)calloc(2, sizeof(int));
        phase_set_to_hap_all_read_count[i] = (int*)calloc(2, sizeof(int));
        phase_set_to_hap_min_full_read_len[i] = (int*)malloc(2 * sizeof(int));
        for (int j = 0; j < 2; ++j) phase_set_to_hap_min_full_read_len[i][j] = INT32_MAX;
    }
    // use one phase set if multiple phase sets exist, reads in other phase sets are considered as non-HAP
    for (int i = 0; i < n_reads; ++i) {
        if (read_haps[i] == 0) continue;
        phase_set_i = add_phase_set(phase_sets[i], uniq_phase_sets, &n_uniq_phase_sets);
        if (LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i])) {
            phase_set_to_hap_full_read_count[phase_set_i][read_haps[i]-1]++;
            phase_set_to_hap_all_read_count[phase_set_i][read_haps[i]-1]++;
            if (phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1] > read_lens[i]) 
                phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1] = read_lens[i];
        } else if (LONGCALLD_NOISY_IS_LEFT_COVER(fully_covers[i]) || LONGCALLD_NOISY_IS_RIGHT_COVER(fully_covers[i])) {
            if (read_lens[i] >= phase_set_to_hap_min_full_read_len[phase_set_i][read_haps[i]-1]) {
                phase_set_to_hap_all_read_count[phase_set_i][read_haps[i]-1]++;
            }
        }
    }
    hts_pos_t max_ps = -1; int max_ps_i = -1; 
    int max_ps_full_read_count1 = -1, max_ps_full_read_count2 = -1;
    for (int i = 0; i < n_uniq_phase_sets; ++i) {
        int phase_set_full_read_count1 = phase_set_to_hap_full_read_count[i][0] < phase_set_to_hap_full_read_count[i][1] ? phase_set_to_hap_full_read_count[i][0] : phase_set_to_hap_full_read_count[i][1];
        int phase_set_full_read_count2 = phase_set_to_hap_full_read_count[i][0] > phase_set_to_hap_full_read_count[i][1] ? phase_set_to_hap_full_read_count[i][0] : phase_set_to_hap_full_read_count[i][1];
        if (phase_set_full_read_count1 > max_ps_full_read_count1) {
            max_ps_full_read_count1 = phase_set_full_read_count1;
            max_ps_full_read_count2 = phase_set_full_read_count2;
            max_ps = uniq_phase_sets[i];
            max_ps_i = i;
        } else if (phase_set_full_read_count1 == max_ps_full_read_count1 && phase_set_full_read_count2 > max_ps_full_read_count2) {
            max_ps_full_read_count2 = phase_set_full_read_count2;
            max_ps = uniq_phase_sets[i];
            max_ps_i = i;
        }
    }
    if (max_ps_full_read_count1 < min_minor_hap_full_read_count) max_ps = -1;
    if (max_ps != -1 && max_ps_i != -1) {
        if (phase_set_to_hap_all_read_count[max_ps_i][0] < min_minor_hap_all_read_count ||
            phase_set_to_hap_all_read_count[max_ps_i][1] < min_minor_hap_all_read_count) max_ps = -1;
    }
    for (int i = 0; i < n_reads; ++i) {
        free(phase_set_to_hap_full_read_count[i]);
        free(phase_set_to_hap_all_read_count[i]);
        free(phase_set_to_hap_min_full_read_len[i]);
    } free(uniq_phase_sets); free(phase_set_to_hap_full_read_count); free(phase_set_to_hap_all_read_count); free(phase_set_to_hap_min_full_read_len);
    return max_ps;
}

// pairwise alignment of cons and ref to append ref to the MSA
// total: 1+ n_reads*2
// 1. ref vs cons: 1
// 2. cons vs n_reads: n_reads
// 3. ref vs n_reads: n_reads
int wfa_collect_noisy_aln_str_with_ps_hap(const call_var_opt_t *opt, int sampling_reads, int n_reads, int *noisy_read_ids, int *lens, uint8_t **seqs, uint8_t *strands, uint8_t **quals, char **names,
                                          int *haps, hts_pos_t *phase_sets, int *fully_covers, hts_pos_t ps, int min_hap_full_reads, int min_hap_all_reads, uint8_t *ref_seq, int ref_seq_len,
                                          int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs, int collect_ref_read_aln_str) {
    // given specific phase_set, collect consensus sequences for each haplotype
    // optional: for long noisy regions, skip read with large edit distance to first read
    int n_cons = 0;
    int total_n_reads = n_reads+2, n_ps_hap_reads = 0;
    int *ps_hap_read_ids = (int*)malloc(total_n_reads * sizeof(int));
    int *ps_hap_read_lens = (int*)malloc(total_n_reads * sizeof(int));
    uint8_t **ps_hap_read_seqs = (uint8_t**)malloc(total_n_reads * sizeof(uint8_t*));
    uint8_t *ps_hap_read_strands = (uint8_t*)malloc(total_n_reads * sizeof(uint8_t));
    uint8_t **ps_hap_read_quals = (uint8_t**)malloc(total_n_reads * sizeof(uint8_t*));
    char **ps_hap_read_names = (char**)malloc(total_n_reads * sizeof(char*));
    int *ps_hap_full_covers = (int*)calloc(total_n_reads, sizeof(int));
    int *cons_lens = (int*)malloc(2 * sizeof(int)); uint8_t **cons_seqs = (uint8_t**)malloc(2 * sizeof(uint8_t*));
    int *msa_seq_lens = (int*)malloc(2 * sizeof(int)); uint8_t ***msa_seqs = (uint8_t***)malloc(2 * sizeof(uint8_t**));
    for (int i = 0; i < 2; ++i) {
        cons_seqs[i] = NULL; cons_lens[i] = 0;
        msa_seqs[i] = (uint8_t**)malloc((n_reads+1) * sizeof(uint8_t*));
        for (int j = 0; j < n_reads+1; ++j) msa_seqs[i][j] = NULL;
        msa_seq_lens[i] = 0;
    }

    int hp_flank_start, hp_flank_end, hp_len, use_non_full = 1;
    if (is_homopolymer(ref_seq, ref_seq_len, opt->noisy_reg_flank_len, &hp_flank_start, &hp_flank_end, &hp_len)) {
        use_non_full = 0;
    }
    for (int hap=1; hap<=2; ++hap) {
        // check if we have enough full-cover reads
        n_ps_hap_reads = 0;
        for (int i = 0; i < n_reads; ++i) {
            if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
            if (use_non_full == 0 && LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
            ps_hap_read_ids[n_ps_hap_reads] = noisy_read_ids[i];
            ps_hap_read_lens[n_ps_hap_reads] = lens[i];
            ps_hap_read_seqs[n_ps_hap_reads] = seqs[i];
            ps_hap_read_strands[n_ps_hap_reads] = strands[i];
            ps_hap_read_quals[n_ps_hap_reads] = quals[i];
            ps_hap_full_covers[n_ps_hap_reads] = fully_covers[i];
            ps_hap_read_names[n_ps_hap_reads] = names[i];
            n_ps_hap_reads++;
        }
        if (ps_hap_read_lens[0] >= opt->max_noisy_reg_len) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "SkipRegion: %" PRIi64 " %d %d %d\n", ps, hap, n_ps_hap_reads, ps_hap_read_lens[0]);
            break;
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "PS: %" PRIi64 " HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
        if (n_ps_hap_reads == 0) continue;
        // collect consensus sequences
        n_cons += abpoa_partial_aln_msa_cons(opt, NULL, sampling_reads, n_ps_hap_reads, ps_hap_read_ids, ps_hap_read_seqs, ps_hap_read_quals, ps_hap_read_lens, ps_hap_full_covers, ps_hap_read_names,
                                             1, cons_lens+hap-1, cons_seqs+hap-1, clu_n_seqs+hap-1, clu_read_ids+hap-1, msa_seq_lens+hap-1, msa_seqs[hap-1]);
    }

    if (n_cons != 2) n_cons = 0;
    else { // collect aln_strs
        for (int hap=1; hap<=2; ++hap) {
            // ref vs cons
            aln_str_t *clu_aln_str = aln_strs[hap-1];
            // fprintf(stderr, "WFA cons-ref align for HAP: %d %d vs %d\n", hap, ref_seq_len, cons_lens[hap-1]);
            wfa_collect_aln_str(opt, ref_seq, ref_seq_len, cons_seqs[hap-1], cons_lens[hap-1], LONGCALLD_NOISY_BOTH_COVER, LONGCALLD_WFA_NO_HEURISTIC, LONGCALLD_WFA_AFFINE_2P, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str));
            n_ps_hap_reads = 0;
            for (int i = 0; i < n_reads; ++i) {
                if (lens[i] <= 0 || phase_sets[i] != ps || haps[i] != hap) continue;
                if (use_non_full == 0 && LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i]) == 0) continue;
                // cons vs read XXX make?
                // heuristic = 1; affine_2p = 0; // 
                // wfa_collect_aln_str(opt, cons_seqs[hap-1], cons_lens[hap-1], seqs[i], lens[i], fully_covers[i], heuristic, affine_2p, LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                make_cons_read_aln_str(opt, msa_seqs[hap-1][clu_n_seqs[hap-1]], msa_seqs[hap-1][n_ps_hap_reads], msa_seq_lens[hap-1], fully_covers[i], LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                // ref vs read
                if (collect_ref_read_aln_str)
                    // fprintf(stderr, "Make ref-read aln str for HAP: %d read_id: %d, %d vs %d (%d)\n", hap, noisy_read_ids[i], lens[i], ref_seq_len, fully_covers[i]);
                    make_ref_read_aln_str(opt, LONGCALLD_REF_CONS_ALN_STR(clu_aln_str), LONGCALLD_CONS_READ_ALN_STR(clu_aln_str, n_ps_hap_reads), LONGCALLD_REF_READ_ALN_STR(clu_aln_str, n_ps_hap_reads));
                n_ps_hap_reads++;
            }
            if (LONGCALLD_VERBOSE >=2 ) fprintf(stderr, "With Ref+Cons PS: %" PRIi64 " HAP: %d n_reads: %d\n", ps, hap, n_ps_hap_reads);
            if (n_ps_hap_reads == 0) continue;
        }
    }
    free(ps_hap_read_ids); free(ps_hap_read_lens); free(ps_hap_read_seqs); free(ps_hap_read_strands); free(ps_hap_read_quals); free(ps_hap_full_covers); free(ps_hap_read_names);
    for (int i = 0; i < 2; ++i) {
        if (cons_seqs[i] != NULL) free(cons_seqs[i]);
        if (msa_seqs[i] != NULL) {
            for (int j = 0; j < n_reads+1; ++j) {
                if (msa_seqs[i][j] != NULL) free(msa_seqs[i][j]);
            }
            free(msa_seqs[i]);
        }
    } free(cons_lens); free(cons_seqs); free(msa_seq_lens); free(msa_seqs);
    return n_cons;
}

int collect_noisy_read_info(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t reg_beg, hts_pos_t reg_end, int noisy_reg_i, int n_noisy_reg_reads, int *noisy_reg_reads, int **read_lens,
                            uint8_t ***read_seqs, uint8_t **strands, uint8_t ***read_quals, char ***read_names, int **fully_covers, int **read_id_to_full_covers, 
                            int **read_reg_beg, int **read_reg_end, int **read_haps, hts_pos_t **phase_sets) {
    *read_lens = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_seqs = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    *strands = (uint8_t*)malloc(n_noisy_reg_reads * sizeof(uint8_t));
    *read_quals = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    *read_names = (char**)malloc(n_noisy_reg_reads * sizeof(char*));
    *fully_covers = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *read_id_to_full_covers = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_reg_beg = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_reg_end = (int*)calloc(chunk->n_reads, sizeof(int));
    *read_haps = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    *phase_sets = (hts_pos_t*)calloc(n_noisy_reg_reads, sizeof(hts_pos_t));

    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        int read_id = noisy_reg_reads[i];
        digar_t *read_digars = chunk->digars+read_id; int n_digar = read_digars->n_digar; digar1_t *digars = read_digars->digars;
        hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
        int reg_read_beg = 0, reg_read_end = digar2qlen(read_digars)-1;
        if (read_digars->digars[0].type == BAM_CHARD_CLIP) reg_read_beg = read_digars->digars[0].len;
        if (read_digars->digars[n_digar-1].type == BAM_CHARD_CLIP) reg_read_end = read_digars->digars[n_digar-1].qi - 1;
        if (LONGCALLD_VERBOSE >= 2) (*read_names)[i] = bam_get_qname(chunk->reads[read_id]);
        else (*read_names)[i] = NULL;
        (*strands)[i] = read_digars->is_rev;
        int beg_is_del = 0, end_is_del = 0, cover = 0;
        for (int digar_i = 0; digar_i < n_digar; ++digar_i) {
            hts_pos_t digar_beg = digars[digar_i].pos, digar_end;
            int op = digars[digar_i].type, len = digars[digar_i].len, qi = digars[digar_i].qi;
            if (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP) continue;
            if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_end = digar_beg + len - 1;
            else digar_end = digar_beg;
            if (digar_beg > reg_end) break;
            if (digar_end < reg_beg) continue;
            if (digar_beg <= reg_beg && digar_end >= reg_beg) {
                if (op == BAM_CDEL) {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi; // qi is on the right side of the DEL
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s noisyReg left boundary is DEL %s:%" PRId64 "-%" PRId64 "\n", (*read_names)[i], chunk->tname, reg_beg, reg_end);
                    if (len > opt->noisy_reg_flank_len) beg_is_del = 1;
                } else {
                    reg_digar_beg = reg_beg;
                    reg_read_beg = qi + (reg_beg - digar_beg);
                }
            }
            if (digar_beg <= reg_end && digar_end >= reg_end) {
                if (op == BAM_CDEL) {
                    reg_digar_end = reg_end;
                    reg_read_end = qi-1; // qi is the one the left side of the DEL
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s noisyReg right boundary is DEL %s:%" PRId64 "-%" PRId64 "\n", (*read_names)[i], chunk->tname, reg_beg, reg_end);
                    if (len > opt->noisy_reg_flank_len) end_is_del = 1;
                } else {
                    reg_digar_end = reg_end;
                    reg_read_end = qi + (reg_end - digar_beg);
                }
            }
        }
        if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) {
            if (beg_is_del == 0 && end_is_del == 0) cover = LONGCALLD_NOISY_LEFT_COVER | LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 3;
            else if (beg_is_del == 0 && end_is_del == 1) cover = LONGCALLD_NOISY_LEFT_COVER | LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 1;
            else if (beg_is_del == 1 && end_is_del == 0) cover = LONGCALLD_NOISY_LEFT_GAP | LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 2;
            else cover = LONGCALLD_NOISY_LEFT_GAP | LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 0;
        } else if (reg_digar_beg == reg_beg) {
            if (beg_is_del) cover = LONGCALLD_NOISY_LEFT_GAP; // (*fully_covers)[i] = 0;
            else cover = LONGCALLD_NOISY_LEFT_COVER; // (*fully_covers)[i] = 1;
        } else if (reg_digar_end == reg_end) {
            if (end_is_del) cover = LONGCALLD_NOISY_RIGHT_GAP; // (*fully_covers)[i] = 0;
            else cover = LONGCALLD_NOISY_RIGHT_COVER; // (*fully_covers)[i] = 2;
        } else cover = 0; // (*fully_covers)[i] = 0;
        // if (2*(reg_read_end-reg_read_beg+1) < (reg_end-reg_beg+1)) return 0;
        (*read_seqs)[i] = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
        (*read_quals)[i] = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
        for (int j = reg_read_beg; j <= reg_read_end; ++j) {
            (*read_seqs)[i][j-reg_read_beg] = seq_nt16_int[bam_seqi(read_digars->bseq, j)];
            (*read_quals)[i][j-reg_read_beg] = read_digars->qual[j];
        }
        (*read_lens)[i] = reg_read_end - reg_read_beg + 1;
        (*read_haps)[i] = chunk->haps[read_id];
        (*phase_sets)[i] = chunk->phase_sets[read_id];
        (*fully_covers)[i] = cover;
        (*read_id_to_full_covers)[read_id] = cover;
        (*read_reg_beg)[read_id] = reg_read_beg; (*read_reg_end)[read_id] = reg_read_end;
    }
    return 0;
}

digar1_t *collect_left_digars(digar_t *digar, int read_noisy_beg, hts_pos_t ref_noisy_beg, int *n_left_digars) {
    digar1_t *left_digars = NULL;
    *n_left_digars = 0; int m_left_digar = 0;
    for (int i = 0; i < digar->n_digar; ++i) {
        int op = digar->digars[i].type, qi = digar->digars[i].qi, digar_qi_end; 
        hts_pos_t digar_ref_beg = digar->digars[i].pos, digar_ref_end;
        // collect left clipping
        if (i == 0 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, digar->digars[i]);
            continue;
        }

        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) digar_qi_end = qi + digar->digars[i].len - 1;
        else digar_qi_end = qi;
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_ref_end = digar->digars[i].pos + digar->digars[i].len - 1;
        else digar_ref_end = digar->digars[i].pos;
        if (qi >= read_noisy_beg && digar_ref_beg >= ref_noisy_beg) break; // reach noisy region

        if (digar_qi_end < read_noisy_beg && digar_ref_end < ref_noisy_beg) // full digar
            left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, digar->digars[i]);
        else if (digar_qi_end >= read_noisy_beg || digar_ref_end >= ref_noisy_beg) { // partial digar, chop X=ID
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                digar1_t d = digar->digars[i]; d.len = read_noisy_beg - d.qi; // chop digar
                left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, d);
            } else if (op == BAM_CDEL) {
                digar1_t d = digar->digars[i]; d.len = ref_noisy_beg - d.pos; // chop digar
                left_digars = push_digar0(left_digars, n_left_digars, &m_left_digar, d);
            }
            break;
        } else {
            fprintf(stderr, "Error: digar_ref_end: %" PRId64 " ref_noisy_beg: %" PRId64 "\n", digar_ref_end, ref_noisy_beg);
            exit(1);
        }
    }
    return left_digars;
}

digar1_t *collect_right_digars(digar_t *digar, int read_noisy_end, hts_pos_t ref_noisy_end, int *n_right_digars) {
    digar1_t *right_digars = NULL; *n_right_digars = 0; int m_right_digar = 0;
    for (int i = 0; i < digar->n_digar; ++i) {
        int qi = digar->digars[i].qi, op = digar->digars[i].type, digar_qi_end;
        hts_pos_t digar_ref_beg = digar->digars[i].pos, digar_ref_end;
        // collect right clipping
        if (i == digar->n_digar-1 && (op == BAM_CSOFT_CLIP || op == BAM_CHARD_CLIP)) {
            right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, digar->digars[i]);
            continue;
        }
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CINS) digar_qi_end = qi + digar->digars[i].len - 1;
        else digar_qi_end = qi;
        if (op == BAM_CDIFF || op == BAM_CEQUAL || op == BAM_CDEL) digar_ref_end = digar->digars[i].pos + digar->digars[i].len - 1;
        else digar_ref_end = digar->digars[i].pos;
        if (digar_qi_end <= read_noisy_end && digar_ref_end <= ref_noisy_end) continue; // left side of noisy region

        if (qi > read_noisy_end && digar_ref_beg > ref_noisy_end) // full digar
            right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, digar->digars[i]);
        else if (qi <= read_noisy_end || digar_ref_beg <= ref_noisy_end) { // partial digar, chop IDX
            if (op == BAM_CINS || op == BAM_CEQUAL || op == BAM_CDIFF) {
                digar1_t d = digar->digars[i]; d.len = digar_qi_end - read_noisy_end; // chop digar
                d.qi = read_noisy_end + 1; // chop qi
                if (op == BAM_CINS) {
                    uint8_t *new_alt_seq = (uint8_t*)malloc(d.len * sizeof(uint8_t));
                    for (int j = 0; j < d.len; ++j) {
                        new_alt_seq[j] = digar->digars[i].alt_seq[j+read_noisy_end+1-qi];
                    }
                    free(d.alt_seq); d.alt_seq = new_alt_seq; // chop alt_seq
                } else d.pos = ref_noisy_end + 1; // chop pos
                right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, d);
            } else if (op == BAM_CDEL) {
                digar1_t d = digar->digars[i]; d.len = digar_ref_end - ref_noisy_end; // chop digar
                d.pos = ref_noisy_end + 1; // chop pos
                right_digars = push_digar0(right_digars, n_right_digars, &m_right_digar, d);
            }
        } else {
            fprintf(stderr, "Error: digar_ref_end: %" PRId64 " ref_noisy_end: %" PRId64 "\n", digar_ref_end, ref_noisy_end);
            exit(1);
        }
    }
    return right_digars;
}

digar1_t *collect_full_msa_digars(int read_beg, int read_end, hts_pos_t ref_beg, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg; hts_pos_t ref_pos = ref_beg;
    int left_read_start = 0, right_read_end = msa_len - 1;
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (read_str[i] == ref_str[i]) { // =
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                } else { // X
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                    d.alt_seq[0] = read_str[i];
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                }
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            ref_pos++;
        }
    }
    return msa_digars;
}

digar1_t *collect_left_msa_digars(int read_beg, int read_end, int qlen, hts_pos_t ref_beg, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg, read_end_pos = read_end; hts_pos_t ref_pos = ref_beg;
    int left_read_start = 0, right_read_end = msa_len - 1, right_skipped_read_base = 0;
    // set read_end_pos right_read_end
    read_end_pos = read_pos-1; int is_covered_by_ref = 0;
    for (int i = msa_len-1; i >= 0; --i) {
        if (ref_str[i] != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && read_str[i] != 5) {
            right_read_end = i; // find the last non-5 base
            break;
        } else if (!is_covered_by_ref && read_str[i] != 5) {
            right_skipped_read_base++;
        }
    }
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] != 5) read_end_pos++;
    }
    for (int i = 0; i < msa_len; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (i >= left_read_start && i <= right_read_end) {
                if (read_str[i] == ref_str[i]) { // =
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                } else { // X
                    digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                    d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                    d.alt_seq[0] = read_str[i];
                    msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
                }
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            if (i >= left_read_start && i <= right_read_end) {
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            ref_pos++;
        }
    }
    if (read_end_pos < qlen-1 || right_skipped_read_base > 0) {
        // push back S digar XXX do extend alignment
        digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_end_pos+1, .type = BAM_CSOFT_CLIP, .len = qlen-1-read_end_pos + right_skipped_read_base, .alt_seq = NULL, .is_low_qual = 0};
        msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
    }
    return msa_digars;
}

digar1_t *collect_right_msa_digars(int read_beg, int read_end, hts_pos_t ref_beg, hts_pos_t ref_end, int msa_len, uint8_t *read_str, uint8_t *ref_str, int *n_msa_digars) {
    digar1_t *msa_digars = NULL;
    *n_msa_digars = 0; int m_msa_digars = 0;
    if (msa_len <= 0) return NULL;
    int read_pos = read_beg; hts_pos_t _ref_pos, ref_pos;
    int left_read_start = 0, right_read_end = msa_len - 1, left_skipped_read_base = 0;
        // set read_pos && left_read_start
    read_pos = read_end+1; _ref_pos = ref_end+1; ref_pos = ref_beg;
    int is_covered_by_ref = 0;
    for (int i = 0; i < msa_len; ++i) {
        if (ref_str[i] != 5) is_covered_by_ref = 1;
        if (is_covered_by_ref && read_str[i] != 5) {
            left_read_start = i; // find the first non-5 base
            break;
        } else if (!is_covered_by_ref && read_str[i] != 5) {
            left_skipped_read_base++;
        }
    }
    for (int i = msa_len-1; i>=0; --i) {
        if (ref_str[i] != 5) _ref_pos--;
        if (read_str[i] != 5) {
            read_pos--;
            ref_pos = _ref_pos;
        }
    }
    if (read_pos > 0 || left_skipped_read_base > 0) { // push front S digar XXX do extend alignment
        digar1_t d = (digar1_t){.pos = ref_pos, .qi = 0, .type = BAM_CSOFT_CLIP, .len = read_pos+left_skipped_read_base, .alt_seq = NULL, .is_low_qual = 0};
        msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
    }
    read_pos += left_skipped_read_base;
    // for (int i = 0; i < msa_len; ++i) {
    for (int i = left_read_start; i <= right_read_end; ++i) {
        if (read_str[i] == 5 && ref_str[i] == 5) continue; // both are non-ACGTN
        if (read_str[i] != 5 && ref_str[i] != 5) {
            if (read_str[i] == ref_str[i]) { // =
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CEQUAL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            } else { // X
                digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDIFF, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
                d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
                d.alt_seq[0] = read_str[i];
                msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            }
            read_pos++; ref_pos++;
        } else if (read_str[i] != 5) { // INS: ref_str[i] == 5
            digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CINS, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
            d.alt_seq = (uint8_t*)malloc(1 * sizeof(uint8_t));
            d.alt_seq[0] = read_str[i];
            msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            read_pos++;
        } else { // DEL read_str[i] == 5, ref_str[i] != 5
            digar1_t d = (digar1_t){.pos = ref_pos, .qi = read_pos, .type = BAM_CDEL, .len = 1, .alt_seq = NULL, .is_low_qual = 0};
            msa_digars = push_digar0(msa_digars, n_msa_digars, &m_msa_digars, d);
            ref_pos++;
        }
    }
    return msa_digars;
}

void update_digars_from_msa1(char *tname, digar_t *digar, int msa_len, uint8_t *ref_str, uint8_t *read_str, int full_cover, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int read_beg, int read_end) {
    int new_n_digars = 0, new_m_digar = 0; digar1_t *new_digars = NULL;
    int old_left_n_digars, old_right_n_digars, msa_n_digars;
    digar1_t *old_left_digars=NULL, *old_right_digars=NULL, *msa_digars=NULL;
    if (LONGCALLD_NOISY_IS_NOT_COVER(full_cover)) return;
    else if (LONGCALLD_NOISY_IS_BOTH_COVER(full_cover) ||  // both left & right
             (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover) && LONGCALLD_NOISY_IS_RIGHT_GAP(full_cover)) ||
             (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover) && LONGCALLD_NOISY_IS_LEFT_GAP(full_cover))) {
        // chop digar
        old_left_digars = collect_left_digars(digar, read_beg, noisy_reg_beg, &old_left_n_digars);
        old_right_digars = collect_right_digars(digar, read_end, noisy_reg_end, &old_right_n_digars);
        msa_digars = collect_full_msa_digars(read_beg, read_end, noisy_reg_beg, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < old_left_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_left_digars[i]);
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
        for (int i = 0; i < old_right_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_right_digars[i]);
    } else if (LONGCALLD_NOISY_IS_LEFT_COVER(full_cover)) { // left-cover: update digar
        // chop digar
        old_left_digars = collect_left_digars(digar, read_beg, noisy_reg_beg, &old_left_n_digars);
        msa_digars = collect_left_msa_digars(read_beg, read_end, digar->qlen, noisy_reg_beg, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < old_left_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_left_digars[i]);
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
    } else if (LONGCALLD_NOISY_IS_RIGHT_COVER(full_cover)) { // right-cover: update digar and start pos
        // chop digar
        old_right_digars = collect_right_digars(digar, read_end, noisy_reg_end, &old_right_n_digars);
        msa_digars = collect_right_msa_digars(read_beg, read_end, noisy_reg_beg, noisy_reg_end, msa_len, read_str, ref_str, &msa_n_digars);
        // push to new digar
        for (int i = 0; i < msa_n_digars; ++i) new_digars = push_digar0(new_digars, &new_n_digars, &new_m_digar, msa_digars[i]);
        for (int i = 0; i < old_right_n_digars; ++i) new_digars = push_digar_alt_seq(new_digars, &new_n_digars, &new_m_digar, old_right_digars[i]);
    }
    if (double_check_digar(new_digars, new_n_digars)) {
        // fprintf(stderr, "%s: Error: after updating digar from MSA, invalid digar!\n", tname);
        // print_digar1(new_digars, new_n_digars, stderr);
        // fprintf(stderr, "old digar:\n");
        // print_digar1(digar->digars, digar->n_digar, stderr);
        free_digar1(new_digars, new_n_digars);
    } else {
        free_digar1(digar->digars, digar->n_digar);
        digar->n_digar = new_n_digars; digar->m_digar = new_m_digar; digar->digars = new_digars;
    }
    if (old_left_digars) free(old_left_digars); if (old_right_digars) free(old_right_digars); if (msa_digars) free(msa_digars);
}

void update_digars_from_aln_str(bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int *read_id_to_full_covers, int *read_reg_beg, int *read_reg_end,
                                aln_str_t **aln_str, int n_cons, int *clu_n_seqs, int **clu_read_ids) {
    for (int i = 0; i < n_cons; ++i) {
        int n_seqs = clu_n_seqs[i];
        for (int read_i = 0; read_i < n_seqs; ++read_i) {
            int read_id = clu_read_ids[i][read_i]; int read_beg = read_reg_beg[read_id], read_end = read_reg_end[read_id];
            aln_str_t *ref_read_aln_str = LONGCALLD_REF_READ_ALN_STR(aln_str[i], read_i);
            update_digars_from_msa1(chunk->tname, chunk->digars+read_id, ref_read_aln_str->aln_len, ref_read_aln_str->target_aln, ref_read_aln_str->query_aln, read_id_to_full_covers[read_id], noisy_reg_beg, noisy_reg_end, read_beg, read_end);
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%s digar->qlen: %d\n", bam_get_qname(chunk->reads[read_id]), digar2qlen(chunk->digars+read_id));
        }
    }
}

// return n_cons
// 1. consensu calling; 2. WFA-based MSA
int collect_noisy_reg_aln_strs(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t noisy_reg_beg, hts_pos_t noisy_reg_end, int noisy_reg_i,
                               int n_noisy_reg_reads, int *noisy_read_ids, uint8_t *ref_seq, int ref_seq_len,
                               int *clu_n_seqs, int **clu_read_ids, aln_str_t **aln_strs) {
    if (n_noisy_reg_reads <= 0) return 0;
    // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
    char **names = NULL; uint8_t **seqs = NULL; uint8_t *strands=NULL; int *fully_covers = NULL, *lens = NULL, *haps = NULL; hts_pos_t *phase_sets = NULL; uint8_t **base_quals = NULL;
    int *read_id_to_full_covers = NULL, *read_reg_beg = NULL, *read_reg_end = NULL;
    collect_noisy_read_info(opt, chunk, noisy_reg_beg, noisy_reg_end, noisy_reg_i, n_noisy_reg_reads, noisy_read_ids,
                            &lens, &seqs, &strands, &base_quals, &names, &fully_covers, &read_id_to_full_covers, 
                            &read_reg_beg, &read_reg_end, &haps, &phase_sets);

    // for large noisy regions, sort reads by fully_cover, error_base_rate, and length
    int sampling_reads = 0;
    if (noisy_reg_end - noisy_reg_beg + 1 >= opt->min_noisy_reg_size_to_sample_reads) sampling_reads = 1;
    sort_noisy_region_reads(n_noisy_reg_reads, noisy_read_ids, lens, seqs, base_quals, strands, fully_covers, names, haps, phase_sets, sampling_reads);
    // >= min_hap_full_read_count reads for each hap && >= min_hap_read_count reads (including not full-cover, but >= full-cover length) for each hap
    int min_hap_full_read_count = opt->min_hap_full_reads, min_hap_read_count = opt->min_hap_reads;
    int min_no_hap_full_read_count = opt->min_dp; // opt->min_no_hap_full_reads;
    // XXX find PS that fully cover the noisy region, allow 1 or 2 haps to have NO fully-covered reads, this can be resolved by POA-extension and consensus mering
    // if no valid PS was found, do haplotype-unaware local-assembly
    hts_pos_t ps_with_both_haps = collect_phase_set_with_both_haps(n_noisy_reg_reads, haps, lens, phase_sets, fully_covers, min_hap_full_read_count, min_hap_read_count);
    int n_full_reads = 0;
    for (int i = 0; i < n_noisy_reg_reads; ++i) if (LONGCALLD_NOISY_IS_BOTH_COVER(fully_covers[i])) n_full_reads++;
    int n_cons = 0;

    int collect_ref_read_aln_str = 0;
    if ((opt->refine_bam && opt->out_aln_fp != NULL) || opt->out_somatic) collect_ref_read_aln_str = 1;
    // XXX only use fully-covered reads, including cliping reads (after re-align to backbone read)
    // two cases to call consensus sequences
    if (ps_with_both_haps > 0) { // call consensus sequences for each haplotype
        n_cons = wfa_collect_noisy_aln_str_with_ps_hap(opt, sampling_reads, n_noisy_reg_reads, noisy_read_ids, lens, seqs, strands, base_quals, names, haps, phase_sets, fully_covers, ps_with_both_haps, min_hap_full_read_count, min_hap_read_count,
                                                       ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs, collect_ref_read_aln_str);
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Hap %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full) n_cons: %d\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads, n_cons);
    // >= min_no_hap_full_read_count full reads in total
    } else if (ps_with_both_haps <= 0 && n_full_reads >= min_no_hap_full_read_count) {
        // XXX do NOT de novo abPOA for homopolymer regions
        n_cons = wfa_collect_noisy_aln_str_no_ps_hap(opt, n_noisy_reg_reads, noisy_read_ids, lens, seqs, names, fully_covers,
                                                     ref_seq, ref_seq_len, clu_n_seqs, clu_read_ids, aln_strs, collect_ref_read_aln_str);
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "NoHap %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full) n_cons: %d\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads, n_cons);
    } else {
        if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "Skipped %s:%" PRIi64 "-%" PRIi64 " %" PRIi64 " %d reads (%d full)\n", chunk->tname, noisy_reg_beg, noisy_reg_end, noisy_reg_end-noisy_reg_beg+1, n_noisy_reg_reads, n_full_reads);
    }
    // update digar based on ref vs read in MSA
    if (n_cons > 0 && ((opt->refine_bam && opt->out_aln_fp != NULL) || opt->out_somatic)) {
        update_digars_from_aln_str(chunk, noisy_reg_beg, noisy_reg_end, read_id_to_full_covers, read_reg_beg, read_reg_end,
                                   aln_strs, n_cons, clu_n_seqs, clu_read_ids);
    }
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        free(seqs[i]); free(base_quals[i]);
    }
    free(names); free(strands); free(seqs); free(base_quals); free(lens); free(read_reg_beg); free(read_reg_end);
    free(fully_covers); free(read_id_to_full_covers); free(haps); free(phase_sets);
    return n_cons;
}