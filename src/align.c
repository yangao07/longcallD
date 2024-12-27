#include "wavefront/wavefront_align.h"
#include "seq.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include "abpoa.h"
#include "edlib.h"
#include "call_var.h"
#include "bam_utils.h"
#include "utils.h"

extern int LONGCALLD_VERBOSE;

int *collect_noisy_read_haps(bam_chunk_t *chunk, int noisy_reg_i) {
    int n_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    int *haps = (int*)calloc(n_reads, sizeof(int));
    for (int j = 0; j < chunk->noisy_reg_to_n_reads[noisy_reg_i]; ++j) {
        int read_i = chunk->noisy_reg_to_reads[noisy_reg_i][j];
        haps[j] = chunk->haps[read_i];
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
            if (op == BAM_CDEL) reg_read_beg = qi;
            else {
                reg_digar_beg = reg_beg;
                reg_read_beg = qi + (reg_beg - digar_beg);
            }
        }
        if (digar_beg <= reg_end && digar_end >= reg_end) {
            if (op == BAM_CDEL) reg_read_end = qi-1;
            else {
                reg_digar_end = reg_end;
                reg_read_end = qi + (reg_end - digar_beg);
            }
        }
    }
    if (reg_digar_beg == reg_beg && reg_digar_end == reg_end) *fully_cover = 3;
    else if (reg_digar_beg == reg_beg) *fully_cover = 1;
    else if (reg_digar_end == reg_end) *fully_cover = 2;
    else *fully_cover = 0;
    if (2*(reg_read_end-reg_read_beg+1) < (reg_end-reg_beg+1)) return 0;
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
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = -1;
    abpt->out_msa = 1;
    abpt->out_cons = 0;
    abpoa_post_set_para(abpt);
    if (LONGCALLD_VERBOSE >= 2) abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, seqs, NULL, stderr);
    else abpoa_msa(ab, abpt, n_seqs, NULL, seq_lens, seqs, NULL, NULL);

    abpoa_cons_t *abc = ab->abc;
    *msa_seqs = (uint8_t**)malloc(n_seqs * sizeof(uint8_t**));
    msa_seq_len = abc->msa_len;
    for (int i = 0; i < n_seqs; ++i) {
        (*msa_seqs)[i] = (uint8_t*)malloc(msa_seq_len * sizeof(uint8_t));
        for (int j = 0; j < msa_seq_len; ++j) {
            (*msa_seqs)[i][j] = abc->msa_base[i][j];
        }
    }
    abpoa_free(ab); abpoa_free_para(abpt); free(seqs); free(seq_lens);
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

int collect_noisy_cons_seq_no_hap(int n_reads, int *read_lens, uint8_t **read_seqs, char **qnames, int *fully_covers, int **cons_lens, uint8_t ***cons_seqs, char *reg_rname, hts_pos_t reg_beg, hts_pos_t reg_end) {
    int n_cons = 0;
    int *full_read_lens = (int*)malloc(n_reads * sizeof(int));
    uint8_t **full_read_seqs = (uint8_t**)malloc(n_reads * sizeof(uint8_t*));
    int n_full_reads = 0;
    for (int i = 0; i < n_reads; ++i) {
        if (read_lens[i] <= 0) continue;
        if (fully_covers[i] != 3) continue;
        full_read_lens[n_full_reads] = read_lens[i];
        full_read_seqs[n_full_reads] = read_seqs[i];
        n_full_reads++;
    }
    if (n_full_reads < 10) {
        fprintf(stderr, "NotEnoughFullReads: %s:%ld-%ld %d %d\n", reg_rname, reg_beg, reg_end, n_full_reads, n_reads);
        free(full_read_lens); free(full_read_seqs);
        return 0;
    }
    for (int i = 0; i < n_full_reads; ++i) {
        for (int j = i+1; j < n_full_reads; ++j) {
            if (full_read_lens[i] < full_read_lens[j]) {
                int tmp_len = full_read_lens[i]; full_read_lens[i] = full_read_lens[j]; full_read_lens[j] = tmp_len;
                uint8_t *tmp_seq = full_read_seqs[i]; full_read_seqs[i] = full_read_seqs[j]; full_read_seqs[j] = tmp_seq;
            }
        }
    }
    // abpoa
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = -1;
    abpt->out_msa = 0;
    abpt->cons_algrm = ABPOA_MF;
    abpt->sort_input_seq = 1;
    abpt->max_n_cons = 2;
    // if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpt->out_cons = 1;
    abpoa_post_set_para(abpt);
    // fprintf(stderr, "Hap%d: %d\n", hap, n_hap_reads);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA-no Hap: %d\n", n_full_reads);
        abpoa_msa(ab, abpt, n_full_reads, NULL, full_read_lens, full_read_seqs, NULL, stderr);
    } else abpoa_msa(ab, abpt, n_full_reads, NULL, full_read_lens, full_read_seqs, NULL, NULL);
    abpoa_cons_t *abc = ab->abc;
    
    if (abc->n_cons > 0) {
        for (int i = 0; i < abc->n_cons; ++i) {
            (*cons_lens)[i] = abc->cons_len[i];
            (*cons_seqs)[i] = (uint8_t*)malloc(abc->cons_len[i] * sizeof(uint8_t));
            for (int j = 0; j < abc->cons_len[i]; ++j) {
                (*cons_seqs)[i][j] = abc->cons_base[i][j];
            }
        }
        n_cons = abc->n_cons;
    } else {
        _err_error_exit("No Consensus: %ld\n");
    }
    abpoa_free(ab); abpoa_free_para(abpt); free(full_read_lens); free(full_read_seqs);
    return n_cons;
}

 int collect_noisy_cons_seq_with_hap(int hap, hts_pos_t ps, int n_reads, int *read_lens, uint8_t **read_seqs, char **qnames, int *read_haps, hts_pos_t *phase_sets, int *fully_covers, int **cons_lens, uint8_t ***cons_seqs, int cons_i, char *tname, hts_pos_t reg_beg, hts_pos_t reg_end) {
    int n_hap_reads = 0;
    int *hap_read_lens = (int*)malloc(n_reads * sizeof(int));
    int n_cons = 0;
    uint8_t **hap_read_seqs = (uint8_t**)malloc(n_reads * sizeof(uint8_t*));
    for (int i = 0; i < n_reads; ++i) {
        if (read_lens[i] <= 0 || read_haps[i] != hap || phase_sets[i] != ps) continue;
        if (fully_covers[i] != 3) {
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "NotFullyCovered: %s %d %d\t", qnames[i], fully_covers[i], read_lens[i]);
                for (int j = 0; j < read_lens[i]; ++j) {
                    fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
                } fprintf(stderr, "\n");
            }
            continue;
        }
        hap_read_lens[n_hap_reads] = read_lens[i];
        hap_read_seqs[n_hap_reads] = read_seqs[i];
        // fprintf(stderr, "%s\t%d\t", qnames[i], read_lens[i]);
        // for (int j = 0; j < read_lens[i]; ++j) {
        //     fprintf(stderr, "%c", "ACGTN"[read_seqs[i][j]]);
        // } fprintf(stderr, "\n");
        n_hap_reads++;
    }
    if (n_hap_reads < 3) {
        fprintf(stderr, "NotEnoughHapReads: %s:%ld-%ld %d\n", tname, reg_beg, reg_end, n_hap_reads);
        free(hap_read_lens); free(hap_read_seqs);
        return 0;
    }
    for (int i = 0; i < n_hap_reads-1; ++i) {
        for (int j = i+1; j < n_hap_reads; ++j) {
            if (hap_read_lens[i] < hap_read_lens[j]) {
                int tmp_len = hap_read_lens[i]; hap_read_lens[i] = hap_read_lens[j]; hap_read_lens[j] = tmp_len;
                uint8_t *tmp_seq = hap_read_seqs[i]; hap_read_seqs[i] = hap_read_seqs[j]; hap_read_seqs[j] = tmp_seq;
            }
        }
    }
    abpoa_t *ab = abpoa_init();
    abpoa_para_t *abpt = abpoa_init_para();
    abpt->wb = -1;
    // abpt->out_msa = 1;
    abpt->cons_algrm = ABPOA_MF;
    abpt->sort_input_seq = 1;
    // if (LONGCALLD_VERBOSE >= 2) abpt->out_msa = 1;
    abpt->out_cons = 1;
    abpoa_post_set_para(abpt);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "For abPOA: Hap%d: %d\n", hap, n_hap_reads);

        abpoa_msa(ab, abpt, n_hap_reads, NULL, hap_read_lens, hap_read_seqs, NULL, stderr);
    } else abpoa_msa(ab, abpt, n_hap_reads, NULL, hap_read_lens, hap_read_seqs, NULL, NULL);
    abpoa_cons_t *abc = ab->abc;
    if (abc->n_cons == 1) {
        (*cons_lens)[cons_i] = abc->cons_len[0];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "ConsLen: %d\n", abc->cons_len[0]);
        (*cons_seqs)[cons_i] = (uint8_t*)malloc(abc->cons_len[0] * sizeof(uint8_t));
        for (int i = 0; i < abc->cons_len[0]; ++i) {
            (*cons_seqs)[cons_i][i] = abc->cons_base[0][i];
        }
        n_cons = 1;
    } else {
        _err_error_exit("No Consensus: %ld\n", ps);
    }
    abpoa_free(ab); abpoa_free_para(abpt); free(hap_read_lens); free(hap_read_seqs);
    return n_cons;
 }

// return 2: HET, 1: HOM
int collect_noisy_cons_seqs(bam_chunk_t *chunk, int noisy_reg_i, int **cons_lens, uint8_t ***cons_seqs) {
    int n_noisy_reg_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    // fprintf(stderr, "Noisy region %d: %d reads\n", noisy_reg_i, n_noisy_reg_reads);
    if (n_noisy_reg_reads <= 0) return 0;
    int *read_haps = collect_noisy_read_haps(chunk, noisy_reg_i);
    hts_pos_t *phase_sets = collect_noisy_read_phase_sets(chunk, noisy_reg_i);
    // int n_reads_with_haps = 0, n_reads_hap1 = 0, n_reads_hap2 = 0;
    // for (int i = 0; i < n_noisy_reg_reads; ++i) {
    //     if (read_haps[i] == 1) {
    //         n_reads_hap1++; n_reads_with_haps++;
    //     } else if (read_haps[i] == 2) {
    //         n_reads_hap2++; n_reads_with_haps++;
    //     }
    // }
    // total haps, read count for each hap
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t reg_beg = cr_start(noisy_regs, noisy_reg_i), reg_end = cr_end(noisy_regs, noisy_reg_i);
    int *noisy_reg_reads = chunk->noisy_reg_to_reads[noisy_reg_i];
    uint8_t **read_seqs = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    char **read_names = (char**)malloc(n_noisy_reg_reads * sizeof(char*));
    int *read_lens = (int*)malloc(n_noisy_reg_reads * sizeof(int)), *fully_covers = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    // XXX clipping?
    // collect all read sequences in the noisy region
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        int read_i = noisy_reg_reads[i], fully_cover = 0;
        // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
        read_lens[i] = collect_reg_read_seq(chunk, read_i, reg_beg, reg_end, &read_seqs[i], &read_names[i], &fully_cover);
        // fprintf(stderr, "%s %d %d\n", read_names[i], read_lens[i], fully_cover);
        fully_covers[i] = fully_cover;
    }
    // haplotype wise consensus
    if (LONGCALLD_VERBOSE >= 2) {
        int hap_counts[3] = {0, 0, 0};
        for (int i = 0; i < n_noisy_reg_reads; ++i) hap_counts[read_haps[i]]++;
        fprintf(stderr, "Hap0: %d, Hap1: %d, Hap2: %d\n", hap_counts[0], hap_counts[1], hap_counts[2]);
    }
    int n_cons = 2;
    int sort_by_length = 1;
    *cons_lens = (int*)malloc(n_cons * sizeof(int));
    *cons_seqs = (uint8_t**)malloc(n_cons * sizeof(uint8_t*)); (*cons_seqs)[0] = NULL; (*cons_seqs)[1] = NULL;
    int n_uniq_phase_sets = 0, phase_set_i = 0;
    hts_pos_t *uniq_phase_sets = (hts_pos_t*)calloc(n_noisy_reg_reads, sizeof(hts_pos_t));
    int **phase_set_to_hap_read_count = (int**)malloc(n_noisy_reg_reads * sizeof(int*));
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        phase_set_to_hap_read_count[i] = (int*)calloc(2, sizeof(int));
    }
    // check if there are multiple phase sets
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        if (read_lens[i] <= 0) continue;
        if (read_haps[i] == 0) continue;
        if (fully_covers[i] != 3) continue;
        // fprintf(stderr, "HAP: %d PS: %ld\n", read_haps[i], phase_sets[i]);
        phase_set_i = add_phase_set(phase_sets[i], uniq_phase_sets, &n_uniq_phase_sets);
        phase_set_to_hap_read_count[phase_set_i][read_haps[i]-1]++;
    }
    hts_pos_t max_ps = -1; int max_ps_read_count1 = -1, max_ps_read_count2 = -1;
    int not_enough_hap_reads, min_hap_read_count = 3; // 5?
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
    not_enough_hap_reads = max_ps_read_count1 < min_hap_read_count ? 1 : 0;

    // XXX always collect two consensus sequences?
    int cons_i = 0;
    if (not_enough_hap_reads) { // call 2 consensus from all reads
        fprintf(stderr, "DeNovoHap: %s:%ld-%ld %d %d %d\n", chunk->tname, reg_beg, reg_end, n_noisy_reg_reads, max_ps_read_count1, max_ps_read_count2);
        cons_i += collect_noisy_cons_seq_no_hap(n_noisy_reg_reads, read_lens, read_seqs, read_names, fully_covers, cons_lens, cons_seqs, chunk->tname, reg_beg, reg_end);
        fprintf(stderr, "n_cons: %d\n", cons_i);
    } else { // call 1 consensus from each set of hap reads
        for (int hap = 1; hap <= 2; ++hap) {
            cons_i += collect_noisy_cons_seq_with_hap(hap, max_ps, n_noisy_reg_reads, read_lens, read_seqs, read_names, read_haps, phase_sets, fully_covers, cons_lens, cons_seqs, cons_i, chunk->tname, reg_end, reg_end);
        }
    }
    for (int i = 0; i < n_noisy_reg_reads; ++i) {
        if (read_lens[i] > 0) free(read_seqs[i]); 
        free(phase_set_to_hap_read_count[i]);
    }
    free(read_names); free(read_seqs); free(read_lens); free(fully_covers); free(read_haps); free(phase_sets);
    free(uniq_phase_sets); free(phase_set_to_hap_read_count);
    if (cons_i == 0) { // no consensus
        free(*cons_lens); free(*cons_seqs);
    }
    return cons_i;
}

hts_pos_t collect_one_phase_set(int n_reads, int *read_haps, hts_pos_t *phase_sets, int *fully_covers) {
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
    int min_hap_read_count = 3; // 5? XXX
    if (max_ps_read_count1 < min_hap_read_count) return -1;
    else return max_ps;
}

// two-rounds of abPOA when there are non-HAP and/or clipping reads
int collect_noisy_cons_seqs0(bam_chunk_t *chunk, int noisy_reg_i, int **cons_lens, uint8_t ***cons_seqs) {
    int n_noisy_reg_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    if (n_noisy_reg_reads <= 0) return 0;
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    int *noisy_reg_reads = chunk->noisy_reg_to_reads[noisy_reg_i];
    uint8_t **read_seqs = (uint8_t**)malloc(n_noisy_reg_reads * sizeof(uint8_t*));
    char **read_names = (char**)malloc(n_noisy_reg_reads * sizeof(char*));
    int *read_lens = (int*)malloc(n_noisy_reg_reads * sizeof(int)), *fully_covers = (int*)calloc(n_noisy_reg_reads, sizeof(int));
    // collect read sequences in the noisy region
    hts_pos_t reg_beg = cr_start(noisy_regs, noisy_reg_i), reg_end = cr_end(noisy_regs, noisy_reg_i);
    for (int i = 0; i < n_noisy_reg_reads; ++i)
        // fully_cover: 0 -> none, 1 -> left, 2 -> right, 3 -> both
        read_lens[i] = collect_reg_read_seq(chunk, noisy_reg_reads[i], reg_beg, reg_end, &read_seqs[i], &read_names[i], fully_covers+i);
    // collect HAP & phase_set information
    int *read_haps = collect_noisy_read_haps(chunk, noisy_reg_i);
    hts_pos_t *phase_sets = collect_noisy_read_phase_sets(chunk, noisy_reg_i);
    hts_pos_t max_ps = collect_one_phase_set(n_noisy_reg_reads, read_haps, phase_sets, fully_covers);

    // 1. abPOA with HAP, if there are enough reads, else (not enough) goto 4
    // 2. re-align NO-HAP or non-fully-covered reads to cons, if necessary
    // 3. re-abPOA with all reads
    // 4. abPOA with all fully-covered reads, without HAP info
    // 5. re-align non-fully-covered reads to cons, if necessary

    for (int i = 0; i < n_noisy_reg_reads; ++i) free(read_seqs[i]);
    free(read_names); free(read_seqs); free(read_lens); free(fully_covers); free(read_haps); free(phase_sets);
    return 0;
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

int collect_wfa_aln(int gap_pos, char *pattern, int plen, char *text, int tlen, char **pattern_alg, char **text_alg) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine_2p;
    attributes.affine2p_penalties.match = 0;
    attributes.affine2p_penalties.mismatch = 4;
    attributes.affine2p_penalties.gap_opening1 = 4;
    attributes.affine2p_penalties.gap_extension1 = 2;
    attributes.affine2p_penalties.gap_opening2 = 24;
    attributes.affine2p_penalties.gap_extension2 = 1;
    attributes.alignment_scope = compute_alignment;
    attributes.alignment_form.span = alignment_end2end;
    attributes.heuristic.strategy = wf_heuristic_none;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    int i, alg_pos = 0, pattern_pos = 0, text_pos = 0;
    int max_len = plen + tlen + 1;
    *pattern_alg = (char*)malloc(max_len * sizeof(char));
    *text_alg = (char*)malloc(max_len * sizeof(char));
    if (gap_pos == LONGCALLD_GAP_LEFT_ALN) { // reverse pattern and text
        for (i = 0; i < plen/2; ++i) {
            char tmp = pattern[i];
            pattern[i] = pattern[plen-i-1];
            pattern[plen-i-1] = tmp;
        }
        for (i = 0; i < tlen/2; ++i) {
            char tmp = text[i];
            text[i] = text[tlen-i-1];
            text[tlen-i-1] = tmp;
        }
    }

    wavefront_align(wf_aligner, pattern, plen, text, tlen);
    cigar_t *cigar = wf_aligner->cigar;
    if (LONGCALLD_VERBOSE >= 2) cigar_print_pretty(stderr, cigar, pattern, plen, text, tlen);
    // fprintf(stderr, "WFA-Alignment returns score %d\n",wf_aligner->cigar->score);
    for (i=cigar->begin_offset;i<cigar->end_offset;++i) {
        switch (cigar->operations[i]) {
        case 'M':
            if (pattern[pattern_pos] != text[text_pos]) {
                (*pattern_alg)[alg_pos] = pattern[pattern_pos];
                (*text_alg)[alg_pos++] = text[text_pos];
            } else {
                (*pattern_alg)[alg_pos] = pattern[pattern_pos];
                (*text_alg)[alg_pos++] = text[text_pos];
            }
            pattern_pos++;
            text_pos++;
            break;
        case 'X':
            if (pattern[pattern_pos] != text[text_pos]) {
                (*pattern_alg)[alg_pos] = pattern[pattern_pos++];
                (*text_alg)[alg_pos++] = text[text_pos++];
            } else {
                (*pattern_alg)[alg_pos] = pattern[pattern_pos++];
                (*text_alg)[alg_pos++] = text[text_pos++];
            }
            break;
        case 'I':
            (*pattern_alg)[alg_pos] = '-';
            (*text_alg)[alg_pos++] = text[text_pos++];
            break;
        case 'D':
            (*pattern_alg)[alg_pos] = pattern[pattern_pos++];
            (*text_alg)[alg_pos++] = '-';
            break;
        default:
            break;
        }
    }
    if (gap_pos == LONGCALLD_GAP_LEFT_ALN) { // reverse pattern and text
        for (i = 0; i < plen/2; ++i) {
            char tmp = pattern[i];
            pattern[i] = pattern[plen-i-1];
            pattern[plen-i-1] = tmp;
        }
        for (i = 0; i < tlen/2; ++i) {
            char tmp = text[i];
            text[i] = text[tlen-i-1];
            text[tlen-i-1] = tmp;
        }
        // reverse pattern and text
        for (i = 0; i < alg_pos/2; ++i) {
            char tmp = (*pattern_alg)[i];
            (*pattern_alg)[i] = (*pattern_alg)[alg_pos-i-1];
            (*pattern_alg)[alg_pos-i-1] = tmp;
            tmp = (*text_alg)[i];
            (*text_alg)[i] = (*text_alg)[alg_pos-i-1];
            (*text_alg)[alg_pos-i-1] = tmp;
        }
    }
    assert(pattern_pos == plen);
    assert(text_pos == tlen);
    // i = 0;
    // while (pattern_pos < plen) {
    //     (*pattern_alg)[alg_pos + i] = pattern[pattern_pos++];
    //     ++i;
    // }
    // i = 0;
    // while (text_pos < tlen) {
    //     (*text_alg)[alg_pos + i] = text[text_pos++];
    //     ++i;
    // }
    // Free
    wavefront_aligner_delete(wf_aligner); 
    return alg_pos;
}

int test_wfa(char *pattern, char *text) {
    wavefront_aligner_attr_t attributes = wavefront_aligner_attr_default;
    attributes.distance_metric = gap_affine;
    attributes.affine_penalties.match = 0;
    attributes.affine_penalties.mismatch = 4;
    attributes.affine_penalties.gap_opening = 6;
    attributes.affine_penalties.gap_extension = 2;
    // Initialize Wavefront Aligner
    wavefront_aligner_t* const wf_aligner = wavefront_aligner_new(&attributes);
    // Align
    wavefront_align(wf_aligner, pattern, strlen(pattern), text,strlen(text));
    fprintf(stderr, "WFA-Alignment returns score %d\n",wf_aligner->cigar->score);
    // Display alignment
    // Free
    wavefront_aligner_delete(wf_aligner);
    return 0;
}