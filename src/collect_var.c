#include <stdio.h>
#include <stdlib.h>
#include "call_var_from_phased_reads.h"
#include "collect_var.h"
#include "call_var.h"
#include "utils.h"
#include "bam_utils.h"
#include "seq.h"
#include "assign_hap.h"
#include "vcf_utils.h"
#include "kmedoids.h"

extern int LONGCALLD_VERBOSE;
char test_read_name[1024] = "m84039_230928_213653_s3/23396787/ccs";

// only 1 variant is stored in each cand_var_t, i.e., SNP/INS/DEL
// each read with be 
//   0: ref
//   1: alt
//   2: other alt. mostly sequencing errors, should be considered as non-informative during haplotype assignment, i.e. 2!=2
cand_var_t *init_cand_vars(int n_var_sites, var_site_t *var_sites) {
    cand_var_t *cand_vars = (cand_var_t*)malloc(n_var_sites * sizeof(cand_var_t));
    for (int i = 0; i < n_var_sites; ++i) {
        // static information
        // from var_sites
        cand_vars[i].tid = var_sites[i].tid;
        cand_vars[i].pos = var_sites[i].pos; 
        cand_vars[i].phase_set = 0; // unset
        cand_vars[i].var_type = var_sites[i].var_type;
        cand_vars[i].ref_len = var_sites[i].ref_len;

        cand_vars[i].n_depth = 0; cand_vars[i].n_low_depth = 0; 
        cand_vars[i].n_uniq_alles = 1; // ref/alt1/alt2/other_alt
        cand_vars[i].alle_covs = (int*)calloc(1, sizeof(int)); // ref/alt1/alt2/other_alt
        cand_vars[i].alt_lens = NULL; cand_vars[i].alt_seqs = NULL; // size: n_uniq_alles-1

        // dynamic information, allocate and update during haplotype assignment
        cand_vars[i].alle_to_hap = NULL; cand_vars[i].hap_to_alle_profile = NULL; cand_vars[i].hap_to_cons_alle = NULL;
        cand_vars[i].is_low_qual = 0; cand_vars[i].is_skipped = 0;
    }
    return cand_vars;
}

void free_cand_vars(cand_var_t *cand_vars, int m) {
    if (m > 0) {
        for (int i = 0; i < m; ++i) {
            if (cand_vars[i].alle_covs != NULL) free(cand_vars[i].alle_covs); 
            if (cand_vars[i].alt_lens != NULL) free(cand_vars[i].alt_lens);
            if (cand_vars[i].alt_seqs != NULL) {
                for (int j = 0; j < cand_vars[i].n_uniq_alles-1; ++j)
                    free(cand_vars[i].alt_seqs[j]);
                free(cand_vars[i].alt_seqs);
            }
            if (cand_vars[i].alle_to_hap != NULL) free(cand_vars[i].alle_to_hap);
            if (cand_vars[i].hap_to_alle_profile != NULL) {
                for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
                    free(cand_vars[i].hap_to_alle_profile[j]);
                } free(cand_vars[i].hap_to_alle_profile);
            }
            if (cand_vars[i].hap_to_cons_alle != NULL) free(cand_vars[i].hap_to_cons_alle);
        }
    }
    free(cand_vars);
}

int merge_x_sites(int n_total_x_sites, hts_pos_t **x_sites, int n_digar, digar1_t *digars) {
    int new_total_x_sites = 0;
    hts_pos_t *new_x_sites = (hts_pos_t*)malloc((n_total_x_sites + n_digar) * sizeof(hts_pos_t));
    int i, j;
    for (i = j = 0; i < n_total_x_sites && j < n_digar; ) {
        if (digars[j].type != BAM_CDIFF || digars[j].is_low_qual) { j++; continue; }
        if ((*x_sites)[i] < digars[j].pos) {
            new_x_sites[new_total_x_sites] = (*x_sites)[i];
            i++;
        } else if ((*x_sites)[i] > digars[j].pos) {
            new_x_sites[new_total_x_sites] = digars[j].pos;
            j++;
        } else {
            new_x_sites[new_total_x_sites] = (*x_sites)[i];
            i++; j++;
        }
        new_total_x_sites++;
    }
    for (; i < n_total_x_sites; ++i, ++new_total_x_sites) new_x_sites[new_total_x_sites] = (*x_sites)[i];
    for (; j < n_digar; ++j) {
        if (digars[j].type != BAM_CDIFF || digars[j].is_low_qual) continue;
        new_x_sites[new_total_x_sites++] = digars[j].pos;
    }
    free(*x_sites); *x_sites = new_x_sites;
    return new_total_x_sites;
}

int collect_cand_vars(bam_chunk_t *chunk, int n_var_sites, var_site_t *var_sites) {
    chunk->cand_vars = init_cand_vars(n_var_sites, var_sites);
    cand_var_t *cand_vars = chunk->cand_vars;
    // 2nd pass: update snp_sites, calculate the depth and allele frequency of each site
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        start_var_i = update_cand_vars_from_digar(chunk->digars+i, read, n_var_sites, var_sites, start_var_i, cand_vars);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_var_sites; ++i) {
            fprintf(stderr, "CandVar: %" PRId64 "\t", var_sites[i].pos);
            fprintf(stderr, "Type: %c\t", BAM_CIGAR_STR[var_sites[i].var_type]);
            fprintf(stderr, "RefLen: %d\t", var_sites[i].ref_len);
            fprintf(stderr, "Depth: %d\t", cand_vars[i].n_depth);
            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                fprintf(stderr, "%d ", j);
                if (j != 0) {
                    for (int k = 0; k < cand_vars[i].alt_lens[j-1]; ++k) {
                        fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seqs[j-1][k]]);
                    }
                } fprintf(stderr, ": ");
                fprintf(stderr, "%d\t", cand_vars[i].alle_covs[j]);
            }
            fprintf(stderr, "Low-Depth: %d\n", cand_vars[i].n_low_depth);
        }
    }
    return 0;
}

// XXX check the sequence around the variant site, not just left/right side
// filter: 
//   1) usable X/= over total non-low-qual depth >= threshold
//   2) total non-low-qual depth / total depth >= threshold

// 1-homopolymer: AAA 1*3
// 2-homopolymer: CGCGCG 2*3
// N-homopolymer: ACGACGACG N*3
int is_homopolymer(char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    // fprintf(stderr, "%d-%d: %s\n", ref_beg, ref_end, ref_seq);
    hts_pos_t start_pos, end_pos; // = var->pos, end_pos = var->pos;
    if (var->var_type == BAM_CDIFF) {
        start_pos = var->pos-1; end_pos = var->pos+1;
    } else if (var->var_type == BAM_CINS) {
        start_pos = var->pos-1; end_pos = var->pos;
    } else { // DEL
        start_pos = var->pos+var->ref_len-1; end_pos = var->pos;
    }
    int var_type = var->var_type;
    int is_homopolymer = 1, max_unit_len=6, n_check_copy_num = 3; // Short Tandem Repeat 1-6 bp
    uint8_t ref_bases[6];
    for (int i = 0; i < 6; ++i)
        ref_bases[i] = nst_nt4_table[(int) ref_seq[end_pos+i-ref_beg]];

    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (nst_nt4_table[(int) ref_seq[end_pos-ref_beg+i*rep_unit_len+j]] != ref_bases[j]) { is_homopolymer = 0; break; }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Homopolymer: %" PRId64 ", %c, unit: %d\n", var->pos, BAM_CIGAR_STR[var_type], rep_unit_len);
            break;
        }
    }
    if (is_homopolymer) return is_homopolymer;
    // reverse
    for (int i = 0; i < 6; ++i) {
        // fprintf(stderr, "%ld, %ld, %ld, %ld\n", start_pos, end_pos, ref_beg, ref_end);
        ref_bases[i] = nst_nt4_table[(int) ref_seq[start_pos-ref_beg-i]];
    }
    for (int rep_unit_len = 1; rep_unit_len <= max_unit_len; ++rep_unit_len) {
        is_homopolymer = 1;
        for (int i = 1; i < n_check_copy_num; ++i) {
            for (int j = 0; j < rep_unit_len; ++j) {
                if (nst_nt4_table[(int) ref_seq[start_pos-ref_beg-i*rep_unit_len-j]] != ref_bases[j]) { is_homopolymer = 0; break; }
            }
            if (is_homopolymer == 0) break;
        }
        if (is_homopolymer) {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Homopolymer: %" PRId64 ", %c, unit: %d (R)\n", var->pos, BAM_CIGAR_STR[var_type], rep_unit_len);
            break;
        }
    }
    return is_homopolymer;
}

// XXX skip detection for long ins/del
int is_repeat_region(char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end, cand_var_t *var) {
    hts_pos_t pos = var->pos; int ref_len = var->ref_len; int var_type = var->var_type;
    uint8_t *ref_bseq, *alt_bseq; int len = 0;
    int is_repeat = 1;
    if (var_type == BAM_CDEL) { // [pos, pos+del_len*N] vs [pos+del_len, pos+del_len*(N+1)]
        // nst_nt4_table['A'] -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4
        int del_len = ref_len;
        len = del_len * 3; // see if del seq is 3-fold repeat
        if (pos < ref_beg || pos+del_len+len >= ref_end) {
            fprintf(stderr, "DelLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", del_len, ref_len, pos, ref_beg, ref_end);
            return 0;
        }
        ref_bseq = get_bseq(ref_seq+pos-ref_beg, len);
        alt_bseq = get_bseq(ref_seq+pos-ref_beg+del_len, len);
        for (int i = 0; i < len; ++i) {
            if (ref_bseq[i] != alt_bseq[i]) { is_repeat = 0; break; }
        }
        free(ref_bseq); free(alt_bseq);
    } else { // if (var_type == BAM_CINS) { // [ins] * N vs [pos, pos+ ins_len * N]
        for (int i = 0; i < var->n_uniq_alles-1; ++i) {
            int ins_len = var->alt_lens[i];
            len = ins_len * 3; // see if ins seq is 3-fold repeat
            if (pos < ref_beg || pos+len >= ref_end) {
                fprintf(stderr, "InsLen: %d, RefLen: %d, Pos: %" PRId64 ", RefBeg: %" PRId64 ", RefEnd: %" PRId64 "\n", ins_len, ref_len, pos, ref_beg, ref_end);
                return 0;
            }
            ref_bseq = get_bseq(ref_seq+pos-ref_beg, len);
            alt_bseq = get_bseq(ref_seq+pos-ref_beg, len);
            for (int j = ins_len; j < len; ++j) alt_bseq[j] = alt_bseq[j-ins_len];
            for (int j = 0; j < ins_len; ++j) alt_bseq[j] = var->alt_seqs[i][j];
            is_repeat = 1;
            for (int k = 0; k < len; ++k) {
                if (ref_bseq[k] != alt_bseq[k]) { is_repeat = 0; break; }
            }
            free(ref_bseq); free(alt_bseq);
            if (is_repeat == 1) break;
        }
    }
    if (is_repeat) {
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "RepeatRegion: %" PRId64 ", type: %c, refLen: %d\n", pos, BAM_CIGAR_STR[var_type], ref_len);
    }
    return is_repeat;
}

// for each candidate variant, keep:
//   ref allele
//   1 alt alleles, potential HOM
//   up to 2 alt alleles, if both AF < max_af & > min_af, potential HET
void refactor_var(cand_var_t *var, int alt_dp1, double alt_af1, int alt_dp2, double alt_af2, 
                  double min_af_thres, double max_af_thres) {
    if (alt_af1 > max_af_thres || alt_af1 < min_af_thres) { // HOM_VAR/SOMATIC: only keep one with the highest dp
        int cov_i = 1, alt_alle_i = 0;
        for (int i = 1; i < var->n_uniq_alles; ++i) {
            if (var->alle_covs[i] == alt_dp1) {
                if (i != cov_i) {
                    var->alle_covs[cov_i] = alt_dp1;
                    var->alt_lens[alt_alle_i] = var->alt_lens[i-1];
                    uint8_t *tmp_seq = var->alt_seqs[alt_alle_i]; var->alt_seqs[alt_alle_i] = var->alt_seqs[i-1]; var->alt_seqs[i-1] = tmp_seq;
                }
                break;
            }
        }
        for (int i = alt_alle_i+1; i < var->n_uniq_alles-1; ++i) {
            free(var->alt_seqs[i]);
        }
        // if (var->n_uniq_alles == 2) var->alle_covs = (int*)realloc(var->alle_covs, 3 * sizeof(int)); // for minor_alt_allele
        var->n_uniq_alles = 2;
    } else { // HET, keep up to 2 alleles with AF < max_af & > min_af
        int cov_idx = 1, alt_alle_i = 0;
        int _alt_dp2 = -1; // alt_af1 < max_af & alt_af1 > min_af
        if (alt_af2 < max_af_thres && alt_af2 > min_af_thres) _alt_dp2 = alt_dp2;
        for (int i = 1; i < var->n_uniq_alles; ++i) {
            if (var->alle_covs[i] == alt_dp1 || var->alle_covs[i] == _alt_dp2) {
                if (i != cov_idx) {
                    var->alle_covs[cov_idx] = var->alle_covs[i];
                    var->alt_lens[alt_alle_i] = var->alt_lens[i-1];
                    uint8_t *tmp_seq = var->alt_seqs[alt_alle_i]; var->alt_seqs[alt_alle_i] = var->alt_seqs[i-1]; var->alt_seqs[i-1] = tmp_seq;
                }
                cov_idx++; alt_alle_i++;
            }
        }
        for (int i = alt_alle_i; i < var->n_uniq_alles-1; ++i) {
            free(var->alt_seqs[i]);
        }
        // if (alt_alle_i == var->n_uniq_alles) var->alle_covs = (int*)realloc(var->alle_covs, (alt_alle_i+1) * sizeof(int)); // for minor_alt_allele
        var->n_uniq_alles = alt_alle_i+1;
    }
}

// 6-tier classification: (5-round variant calling)
// classify all candidate variants into different categories, process each category sequentially
//   0) low total cov, skipped
//   1) easy-to-call germ+het(SNP or indel): high cov, min_af < AF < max_af - basic phasing info
//   2) repeat region germ+het: high cov, min_af < AF < max_af - additional phasing info
//   3) dense region (hom/het): high cov - (potential) additional phasing info
//   4) somatic: high cov, AF < min_af - require full phasing info (alt_af < min_af)
//   5) hom: high cov, AF > max_af - require full phasing info (alt_af > max_af)
// Note: category of variant can be changed during the process, i.e., from 0) to 1), or 0) to 4)
// XXX remove minor alleles in cand_var with low cov/AF
int classify_var_cate(char *ref_seq, hts_pos_t ref_beg, hts_pos_t ref_end,
                      cand_var_t *var, int min_dp_thres, int min_alt_dp_thres, double min_somatic_af_thres,
                      double min_af_thres, double max_af_thres, double max_low_qual_frac_thres) {
    // total depth < min_dp or total alt depth < min_alt_dp: (skip)
    if (var->n_depth + var->n_low_depth < min_dp_thres) return LONGCALLD_LOW_COV_VAR;
    // low-qual depth / total depth > max_low_qual_frac: too many low-qual bases, eg., dense X/gap region
    double af=0.0, alt_af1=0.0, alt_af2=0.0; // alt alleles
    int ad, alt_dp1=0, alt_dp2=0;
    for (int j = 1; j < var->n_uniq_alles; ++j) { // alt allele
        ad = var->alle_covs[j];
        af = (double) ad / var->n_depth;
        if (ad > alt_dp1) { alt_dp2 = alt_dp1; alt_dp1 = ad; }
        else if (ad > alt_dp2 && ad <= alt_dp1) alt_dp2 = ad;
        if (af > alt_af1) { alt_af2 = alt_af1; alt_af1 = af; }
        else if (af > alt_af2 && af <= alt_af1) alt_af2 = af;
    }
    if (alt_dp1 < min_alt_dp_thres || alt_af1 < min_somatic_af_thres) return LONGCALLD_LOW_COV_VAR;

    refactor_var(var, alt_dp1, alt_af1, alt_dp2, alt_af2, min_af_thres, max_af_thres);

    if ((double) var->n_low_depth / (var->n_low_depth + var->n_depth) > max_low_qual_frac_thres) return LONGCALLD_DENSE_REG_VAR; // require basic phasing info, MSA, provide additional phasing info
    // AF < min_af or AF > max_af
    if (alt_af1 > max_af_thres) return LONGCALLD_CAND_HOM_VAR;
    // if (af1 < min_af || af1 > max_af || af2 < min_af || af2 > max_af) return LONGCALLD_CAND_SOMA_VAR; // unlikely germline het., likely hom or somatic, require full phasing info
    // if (total_alt_af < min_af || total_alt_af > max_af || (1-total_alt_af) < min_af || (1-total_alt_af) > max_af) return LONGCALLD_CAND_SOMA_VAR; // unlikely germline het., likely hom or somatic, require full phasing info
    // indels in homo/repeat regions
    if ((var->var_type == BAM_CINS || var->var_type == BAM_CDEL) && (is_homopolymer(ref_seq, ref_beg, ref_end, var) || is_repeat_region(ref_seq, ref_beg, ref_end, var))) return LONGCALLD_REP_HET_VAR; // require basic phasing info, MSA, provide additional phasing info
    // not call somatic variant around homopolymer/repeat region
    if (alt_af1 < min_af_thres) return LONGCALLD_CAND_SOMA_VAR; // XXX could be het in homopolymer region
    if (var->var_type == BAM_CDIFF) return LONGCALLD_EASY_HET_SNP;
    else return LONGCALLD_EASY_HET_INDEL;
    // return LONGCALLD_EASY_HET_VAR;
}

void copy_var(cand_var_t *to_var, cand_var_t *from_var) {
    to_var->pos = from_var->pos; to_var->var_type = from_var->var_type;
    to_var->n_depth = from_var->n_depth; to_var->n_low_depth = from_var->n_low_depth;
    to_var->ref_len = from_var->ref_len;
    to_var->is_low_qual = from_var->is_low_qual; to_var->is_skipped = from_var->is_skipped;

    // allele_covs, alt_alle_seqs, alt_len
    if (to_var->alt_lens) free(to_var->alt_lens); 
    if (to_var->alt_seqs) {
        for (int i = 0; i < to_var->n_uniq_alles-1; ++i) free(to_var->alt_seqs[i]);
        free(to_var->alt_seqs);
    }
    to_var->alt_lens = (int*)malloc((from_var->n_uniq_alles-1) * sizeof(int));
    to_var->alt_seqs = (uint8_t**)malloc((from_var->n_uniq_alles-1) * sizeof(int8_t*));
    if (to_var->alle_covs) free(to_var->alle_covs);
    to_var->alle_covs = (int*)malloc(from_var->n_uniq_alles * sizeof(int));
    to_var->n_uniq_alles = from_var->n_uniq_alles;

    for (int j = 0; j < to_var->n_uniq_alles; ++j) {
        to_var->alle_covs[j] = from_var->alle_covs[j];
        if (j != 0) {
            to_var->alt_lens[j-1] = from_var->alt_lens[j-1];
            to_var->alt_seqs[j-1] = (uint8_t*)malloc(from_var->alt_lens[j-1] * sizeof(uint8_t));
            for (int k = 0; k < from_var->alt_lens[j-1]; ++k) {
                to_var->alt_seqs[j-1][k] = from_var->alt_seqs[j-1][k];
            }
        }
    }
}

int classify_cand_vars(bam_chunk_t *chunk, int n_var_sites, const call_var_opt_t *opt) {
    cand_var_t *cand_vars = chunk->cand_vars;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    // chunk->var_cate_counts = (int*)calloc(LONGCALLD_VAR_CATE_N, sizeof(int));
    // chunk->var_cate_idx = (int**)malloc(LONGCALLD_VAR_CATE_N * sizeof(int*));
    chunk->var_i_to_cate = (int*)malloc(n_var_sites * sizeof(int));
    // for (int i = 0; i < LONGCALLD_VAR_CATE_N; ++i) {
        // chunk->var_cate_idx[i] = (int*)malloc(n_var_sites * sizeof(int));
    // }
    // int **var_cate_idx = chunk->var_cate_idx; int *var_cate_counts = chunk->var_cate_counts; 
    int *var_i_to_cate = chunk->var_i_to_cate;
    int min_dp = opt->min_dp, min_alt_dp = opt->min_alt_dp;
    double min_somatic_af = opt->min_somatic_af, min_af = opt->min_af, max_af = opt->max_af, max_low_qual_frac = opt->max_low_qual_frac;
    int var_cate = -1, cand_var_i = 0;
    for (int i = 0; i < n_var_sites; ++i) {
        cand_var_t *var = cand_vars+i;
        var_cate = classify_var_cate(ref_seq, ref_beg, ref_end, var, min_dp, min_alt_dp, min_somatic_af, min_af, max_af, max_low_qual_frac);
        if (var_cate == LONGCALLD_LOW_COV_VAR) continue; // skipped
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "CandVarCate-%c: %s:%" PRId64 ", type: %c, refLen: %d, depth: %d\t", LONGCALLD_VAR_CATE_TYPE(var_cate), chunk->tname, cand_vars[i].pos, BAM_CIGAR_STR[cand_vars[i].var_type], cand_vars[i].ref_len, cand_vars[i].n_depth);
            fprintf(stderr, "Low-Depth: %d\t", cand_vars[i].n_low_depth);
            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                fprintf(stderr, "%d ", j);
                if (j != 0) {
                    for (int k = 0; k < cand_vars[i].alt_lens[j-1]; ++k) {
                        fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seqs[j-1][k]]);
                    }
                } fprintf(stderr, ": ");
                fprintf(stderr, "%d\t", cand_vars[i].alle_covs[j]);
            } fprintf(stderr, "\n");
        }
        // copy i'th to cand_var_i'th
        if (i != cand_var_i) copy_var(cand_vars+cand_var_i, cand_vars+i);
        // push var to var_cate's list
        var_i_to_cate[cand_var_i++] = var_cate;
        // var_cate_idx[var_cate][var_cate_counts[var_cate]++] = cand_var_i++;
    }
    for (int i = cand_var_i; i < n_var_sites; ++i) {
        free(cand_vars[i].alle_covs); free(cand_vars[i].alt_lens); 
        for (int j = 0; j < cand_vars[i].n_uniq_alles-1; ++j) free(cand_vars[i].alt_seqs[j]);
        free(cand_vars[i].alt_seqs);
    }
    return(chunk->n_cand_vars = cand_var_i);
}

/*int filter_cand_vars(bam_chunk_t *chunk, int n_var_sites, kstring_t *ref_seq, const call_var_opt_t *opt) {
    cand_var_t *cand_vars = chunk->cand_vars;
    int cand_var_i = 0;
    int min_dp = opt->min_dp; double min_af = opt->min_af, max_af = opt->max_af, max_low_qual_frac = opt->max_low_qual_frac;
    for (int i = 0; i < n_var_sites; ++i) {
        // 1. total non-low-qual depth < min_dp
        if (cand_vars[i].n_depth < min_dp) continue;
        // 2. low-qual depth / total depth > max_low_qual_frac
        if ((double) cand_vars[i].n_low_depth / (cand_vars[i].n_low_depth + cand_vars[i].n_depth) > max_low_qual_frac) continue;
        double _af1 = 0.0, _af2 = 0.0, _af = 0.0;
        for (int j = 0; j < 2; ++j) {
            _af = (double) cand_vars[i].alle_covs[j] / cand_vars[i].n_depth;
            if (_af > _af1) { _af2 = _af1; _af1 = _af; } 
            else if (_af > _af2 && _af <= _af1) _af2 = _af;
        }
        // 3. AF < min_af or AF > max_af
        if (_af1 < min_af || _af1 > max_af || _af2 < min_af || _af2 > max_af) continue;
        // 4. indels in repeat regions
        if ((cand_vars[i].var_type == BAM_CINS || cand_vars[i].var_type == BAM_CDEL) && is_repeat_region(ref_seq, cand_vars+i))
            continue;
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "FilteredCandVar %d: %lld, type: %c, refLen: %d, depth: %d, AF: %f|%f\t", cand_var_i, (long long) cand_vars[i].pos, BAM_CIGAR_STR[cand_vars[i].var_type], cand_vars[i].ref_len, cand_vars[i].n_depth, _af1, _af2);
            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                fprintf(stderr, "%d ", j);
                if (j == 1) {
                    for (int k = 0; k < cand_vars[i].alt_len; ++k) {
                        fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_seq[k]]);
                    }
                } fprintf(stderr, ": ");
                fprintf(stderr, "%d\t", cand_vars[i].alle_covs[j]);
            }
            fprintf(stderr, "Low-Depth: %d\n", cand_vars[i].n_low_depth);
        }
        // copy i'th to cand_var_i'th
        if (i != cand_var_i) {
            cand_vars[cand_var_i].pos = cand_vars[i].pos; cand_vars[cand_var_i].var_type = cand_vars[i].var_type;
            cand_vars[cand_var_i].n_depth = cand_vars[i].n_depth; cand_vars[cand_var_i].n_low_depth = cand_vars[i].n_low_depth;
            cand_vars[cand_var_i].ref_len = cand_vars[i].ref_len; 
            cand_vars[cand_var_i].is_low_qual = cand_vars[i].is_low_qual; cand_vars[cand_var_i].is_skipped = cand_vars[i].is_skipped;

            // allele_covs, alt_alle_seqs, alt_len
            if (cand_vars[i].n_uniq_alles > cand_vars[cand_var_i].n_uniq_alles) {
                cand_vars[cand_var_i].alle_covs = (int*)realloc(cand_vars[cand_var_i].alle_covs, cand_vars[i].n_uniq_alles * sizeof(int));
            }
            free(cand_vars[cand_var_i].alt_seq);
            cand_vars[cand_var_i].alt_seq = (uint8_t*)malloc(cand_vars[i].alt_len * sizeof(uint8_t*));

            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                cand_vars[cand_var_i].alle_covs[j] = cand_vars[i].alle_covs[j];
            }
            cand_vars[cand_var_i].alt_len = cand_vars[i].alt_len;
            for (int k = 0; k < cand_vars[i].alt_len; ++k) {
                cand_vars[cand_var_i].alt_seq[k] = cand_vars[i].alt_seq[k];
            }
            cand_vars[cand_var_i].n_uniq_alles = cand_vars[i].n_uniq_alles;
        }
        cand_var_i++;
    }
    for (int i = cand_var_i; i < n_var_sites; ++i) {
        free(cand_vars[i].alle_covs);
        free(cand_vars[i].alt_seq);
    }
    return (chunk->n_cand_vars = cand_var_i);
}*/

void collect_noisy_reg_reads(bam_chunk_t *chunk) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    chunk->noisy_reg_to_n_reads = (int*)calloc(noisy_regs->n_r, sizeof(int));
    chunk->noisy_reg_to_reads = (int**)malloc(noisy_regs->n_r * sizeof(int*));
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    hts_pos_t beg, end;
    
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            chunk->noisy_reg_to_n_reads[ovlp_b[ovlp_i]]++;
        }
    }
    for (int i = 0; i < noisy_regs->n_r; ++i) {
        chunk->noisy_reg_to_reads[i] = (int*)malloc(chunk->noisy_reg_to_n_reads[i] * sizeof(int));
        chunk->noisy_reg_to_n_reads[i] = 0;
    }
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        beg = chunk->digars[i].beg; end = chunk->digars[i].end;
        ovlp_n = cr_overlap(noisy_regs, "cr", beg-1, end, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            chunk->noisy_reg_to_reads[ovlp_b[ovlp_i]][chunk->noisy_reg_to_n_reads[ovlp_b[ovlp_i]]++] = i;
        }
    }
    // print read names for each noisy region
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < noisy_regs->n_r; ++i) {
            if (cr_label(noisy_regs, i) < 10) continue;
            fprintf(stderr, "NoisyRegReads: %s:%d-%d %d %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i) - cr_start(noisy_regs, i), cr_label(noisy_regs, i));
            for (int j = 0; j < chunk->noisy_reg_to_n_reads[i]; ++j) {
                fprintf(stderr, "Read: %s\n", bam_get_qname(chunk->reads[chunk->noisy_reg_to_reads[i][j]]));
            }
        }
    }
    free(ovlp_b);
}

void collect_clip_reads(bam_chunk_t *chunk) {

}

void collect_digars_from_bam(bam_chunk_t *chunk, const call_var_pl_t *pl) {
    chunk->chunk_noisy_regs = cr_init();
    chunk->large_gap_regs = cr_init();
    const call_var_opt_t *opt = pl->opt;
    // if (LONGCALLD_VERBOSE >= 2)
        // fprintf(stderr, "CHUNK: %s\tbeg: %" PRId64 ", end: %" PRId64 ", total_n: %d, ovlp_n: %d\n", chunk->tname, chunk->beg, chunk->end, chunk->n_reads, chunk->n_up_ovlp_reads);
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%d: qname: %s, flag: %d, pos: %" PRId64 ", end: %" PRId64 "\n", i, bam_get_qname(read), read->core.flag, read->core.pos+1, bam_endpos(read));
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        if (read->core.qual < opt->min_mq || read->core.flag & BAM_FSUPPLEMENTARY || read->core.flag & BAM_FSECONDARY) {
            chunk->is_skipped[i] = 1; 
            // fprintf(stderr, "Skip read %s with low mapping quality %d\n", bam_get_qname(read), read->core.qual);
            continue;
        }
        if (chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            collect_digar_from_eqx_cigar(chunk, read, opt, chunk->digars+i);
        } else if (chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
            collect_digar_from_MD_tag(chunk, read, opt, chunk->digars+i);
        // XXX TODO use cs tag
        } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
            collect_digar_from_ref_seq(chunk, read, opt, chunk->digars+i);
        }
    }
    cr_index(chunk->chunk_noisy_regs);
    cr_index(chunk->large_gap_regs);
    chunk->chunk_noisy_regs = cr_merge(chunk->chunk_noisy_regs);
    chunk->large_gap_regs = cr_merge(chunk->large_gap_regs);
    collect_noisy_reg_reads(chunk);
    // collect_clip_reads(chunk);
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "ChunkNoisyRegs-\tchr11:start-end 20\n");
        for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
            fprintf(stderr, "ChunkNoisyRegs: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->chunk_noisy_regs, i), cr_end(chunk->chunk_noisy_regs, i), cr_label(chunk->chunk_noisy_regs, i));
        }
        for (int i = 0; i < chunk->large_gap_regs->n_r; ++i) {
            fprintf(stderr, "LargeGapRegs: %s:%d-%d %d\n", chunk->tname, cr_start(chunk->large_gap_regs, i), cr_end(chunk->large_gap_regs, i), cr_label(chunk->large_gap_regs, i));
        }
    }
}

int update_digars_with_clip(bam_chunk_t *chunk, call_var_pl_t *pl) {
    return 0;
}

// XXX should be digar->pos-1 for INS/DEL
var_site_t make_var_site_from_digar(int tid, digar1_t *digar) {
    var_site_t var_site = {tid, digar->pos, digar->type, 1};
    if (digar->type == BAM_CINS) var_site.ref_len = 0;
    else if (digar->type == BAM_CDEL) var_site.ref_len = digar->len;
    return var_site;
}

var_site_t make_var_site_from_cand_var(cand_var_t *var) {
    var_site_t var_site = {var->tid, var->pos, var->var_type, var->ref_len};
    return var_site;
}

// check if two var_sites overlap based on their pos & ref_len
int ovlp_var_site(var_site_t *var1, var_site_t *var2) {
    int beg1 = var1->pos; int end1 = var1->pos + var1->ref_len;
    int beg2 = var2->pos; int end2 = var2->pos + var2->ref_len;
    if (var1->ref_len == 0 && var2->ref_len == 0) { // both are INS
        if (beg1 == beg2) return 1;
        else return 0;
    } else if (var1->ref_len == 0) { // var1 is INS, var2 is SNP/DEL
        if (beg1 > beg2 && end1 < end2) return 1;
        else return 0;
    } else if (var2->ref_len == 0) { // var2 is INS, var1 is SNP/DEL
        if (beg2 > beg1 && end2 < end1) return 1;
        else return 0;
    } else { // both are SNP/DEL
        if (beg1 >= end2 || beg2 >= end1) return 0; // no overlap
        return 1; // overlap
    }
}


// order of variants with same pos: order by type, ref_len, alt_len
int comp_var_site(var_site_t *var1, var_site_t *var2) {
    if (var1->pos < var2->pos) return -1;
    if (var1->pos > var2->pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    // if (var1->alt_len < var2->alt_len) return -1;
    // if (var1->alt_len > var2->alt_len) return 1;
    // alt_len == alt_len
    // return memcmp(var1->alt_seq, var2->alt_seq, var1->alt_len);
    return 0;
}


int comp_ovlp_var_site(var_site_t *var1, var_site_t *var2, int *is_ovlp) {
    *is_ovlp = ovlp_var_site(var1, var2);
    return comp_var_site(var1, var2);
}

// order of variants with same pos: order by type and ref_len
int merge_var_sites(int n_total_var_sites, var_site_t **var_sites, int tid, hts_pos_t reg_beg, hts_pos_t reg_end,
                    int n_digar, digar1_t *digars) {
    int new_total_var_sites = 0;
    var_site_t *new_var_sites = (var_site_t*)malloc((n_total_var_sites + n_digar) * sizeof(var_site_t));
    int i, j;
    for (i = j = 0; i < n_total_var_sites && j < n_digar; ) {
        hts_pos_t digar_pos = digars[j].pos;
        if (reg_beg != -1 && digar_pos < reg_beg) { j++; continue; }
        if (reg_end != -1 && digar_pos > reg_end) { break; }
        // if (!is_in_reg(reg_cr, tname, digar_pos-1, digar_pos)) { j++; continue; }
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
        // if (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL) { // keep noisy-regions variants
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            int ret = comp_var_site((*var_sites)+i, &digar_var_site);
            if (ret < 0) {
                new_var_sites[new_total_var_sites++] = (*var_sites)[i];
                i++;
            } else if (ret > 0) {
                new_var_sites[new_total_var_sites++] = digar_var_site;
                j++;
            } else { // merge
                new_var_sites[new_total_var_sites++] = (*var_sites)[i];
                i++; j++;
            }
        } else j++;
    }
    for (; i < n_total_var_sites; ++i) new_var_sites[new_total_var_sites++] = (*var_sites)[i];
    for (; j < n_digar; ++j) {
        hts_pos_t digar_pos = digars[j].pos;
        if (reg_beg != -1 && digar_pos < reg_beg) { j++; continue; }
        if (reg_end != -1 && digar_pos > reg_end) { break; }
        // if (!is_in_reg(reg_cr, tname, digar_pos-1, digar_pos)) continue;
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
            var_site_t digar_var_site = make_var_site_from_digar(tid, digars+j);
            new_var_sites[new_total_var_sites++] = digar_var_site;
        }
    }
    free(*var_sites); *var_sites = new_var_sites;
    return new_total_var_sites;
}

// collect all candidate variant sites from digars, excluding low-quality/noisy-region ones
int collect_all_cand_var_sites(bam_chunk_t *chunk, var_site_t **var_sites) {
    int n_total_var_sites = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
       n_total_var_sites = merge_var_sites(n_total_var_sites, var_sites, 
                                           chunk->tid, chunk->reg_beg, chunk->reg_end, //->beg, chunk->end, 
                                           chunk->digars[i].n_digar, chunk->digars[i].digars);
    }
    fprintf(stderr, "total_cand_var: %d\n", n_total_var_sites);
    // print cand_vars
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Total candidate variant sites: %d\n", n_total_var_sites);
        for (int i = 0; i < n_total_var_sites; ++i) {
            fprintf(stderr, "CandVarSite: %s:%" PRId64 "-%c\n", chunk->tname, (*var_sites)[i].pos, BAM_CIGAR_STR[(*var_sites)[i].var_type]);
        }
    }
    return n_total_var_sites;
}

read_var_profile_t *collect_read_var_profile(bam_chunk_t *chunk) {
    int n_cand_vars = chunk->n_cand_vars;
    cand_var_t *cand_vars = chunk->cand_vars;
    read_var_profile_t *p = init_read_var_profile(chunk->n_reads, n_cand_vars);
    cgranges_t *read_var_cr = cr_init();
    // 3rd pass: collect read-wise SNP profiles
    int start_var_i = 0;
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->is_skipped[i]) continue;
        bam1_t *read = chunk->reads[i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        start_var_i = update_read_var_profile_from_digar(chunk->digars+i, read, n_cand_vars, cand_vars, start_var_i, p+i);
        if (p[i].start_var_idx < 0 || p[i].end_var_idx < 0) continue;
        cr_add(read_var_cr, "cr", p[i].start_var_idx, p[i].end_var_idx+1, i);
        if (LONGCALLD_VERBOSE >= 2) {
            if (p[i].start_var_idx >= 0) {
                fprintf(stderr, "Read: %s, start_var_i: %d, end_var_i: %d\n", bam_get_qname(read), p[i].start_var_idx, p[i].end_var_idx);
                for (int j = 0; j <= p[i].end_var_idx-p[i].start_var_idx; ++j) {
                    if (p[i].var_is_used[j] == 0) continue;
                    fprintf(stderr, "P\tVar: (%d) %" PRId64 "", j, cand_vars[j+p[i].start_var_idx].pos);
                    fprintf(stderr, " %d-%c, allele: %d\n", cand_vars[j+p[i].start_var_idx].ref_len, BAM_CIGAR_STR[cand_vars[j+p[i].start_var_idx].var_type], p[i].alleles[j]);
                }
            }
        }
    }
    cr_index(read_var_cr); chunk->read_var_cr = read_var_cr;
    return p;
}

// collect variants based on hap_to_cons_alle
// 1) filter out low-quality variants, e.g., P(var|hap,phasing)
// 2) merge variants with the overlapping pos & ref_len (e.g., ACGT -> A, ACGTCGT) XXX
// 3) extend phase blocks if possible, e.g., if ≥ 1 read supports the longer phase block
int make_variants(bam_chunk_t *chunk, var_t *var) {
    int n_cand_vars = chunk->n_cand_vars;
    char *ref_seq = chunk->ref_seq; hts_pos_t ref_beg = chunk->ref_beg, ref_end = chunk->ref_end;
    if (n_cand_vars <= 0) return 0;
    cand_var_t *cand_vars = chunk->cand_vars;
    var->n = 0; var->m = n_cand_vars;
    int flip = chunk->flip_hap; int hap1_idx, hap2_idx;
    int64_t ovlp_n, *ovlp_b = 0, max_b = 0;
    if (flip) {
        hap1_idx = 2; hap2_idx = 1;
    } else {
        hap1_idx = 1; hap2_idx = 2;
    }
    var->vars = (var1_t*)malloc(n_cand_vars * sizeof(var1_t));
    int i = 0; int hap_alles[2];
    for (int cand_i = 0; cand_i < n_cand_vars; ++cand_i) {
        if (chunk->var_i_to_cate[cand_i] != LONGCALLD_EASY_HET_SNP &&
            chunk->var_i_to_cate[cand_i] != LONGCALLD_EASY_HET_INDEL) continue;
        hap_alles[0] = cand_vars[cand_i].hap_to_cons_alle[hap1_idx];
        // if (hap_alles[0] == -1) hap_alles[0] = LONGCALLD_REF_ALLELE;
        hap_alles[1] = cand_vars[cand_i].hap_to_cons_alle[hap2_idx];
        // if (hap_alles[1] == -1) hap_alles[1] = LONGCALLD_REF_ALLELE;
        // only keep het. vars
        if (hap_alles[0] == -1 && hap_alles[1] == -1) {
            fprintf(stderr, "skipped pos(-1): %" PRId64 ", %d-%c\n", cand_vars[cand_i].pos, cand_vars[cand_i].ref_len, BAM_CIGAR_STR[cand_vars[cand_i].var_type]);
            continue;
        }
        if (hap_alles[0] == hap_alles[1]) { // potential hom var
            fprintf(stderr, "skipped pos(==): %" PRId64 ", %d-%c\n", cand_vars[cand_i].pos, cand_vars[cand_i].ref_len, BAM_CIGAR_STR[cand_vars[cand_i].var_type]);
            continue;
        }

        if (hap_alles[0] == -1) hap_alles[0] = LONGCALLD_REF_ALLELE;
        if (hap_alles[1] == -1) hap_alles[1] = LONGCALLD_REF_ALLELE;

        var->vars[i].type = cand_vars[cand_i].var_type;
        var->vars[i].PS = cand_vars[cand_i].phase_set;
        var->vars[i].ref_len = cand_vars[cand_i].ref_len;
        if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) {
            var->vars[i].pos = cand_vars[cand_i].pos-1;
            var->vars[i].ref_len += 1;
        } else var->vars[i].pos = cand_vars[cand_i].pos;
        if (cr_overlap(chunk->reg_cr, chunk->tname, var->vars[i].pos-1, var->vars[i].pos+var->vars[i].ref_len-1, &ovlp_b, &max_b) <= 0)
            continue;
        // if (var->vars[i].pos == 49735305)
            // fprintf(stderr, "OK\n");
        var->vars[i].ref_bases = get_bseq(ref_seq+var->vars[i].pos-ref_beg, var->vars[i].ref_len);
        var->vars[i].alt_len = (int*)malloc(2 * sizeof(int));
        var->vars[i].alt_bases = (uint8_t**)malloc(2 * sizeof(uint8_t*));
        var->vars[i].n_alt_allele = 0;
        for (int hap=1; hap <= 2; ++hap) {
            int hap_alle = hap_alles[hap-1];
            if (hap_alle != 0) { // alt allele
                int alt_len = cand_vars[cand_i].alt_lens[hap_alle-1];
                var->vars[i].alt_bases[var->vars[i].n_alt_allele] = (uint8_t*)malloc((alt_len+1) * sizeof(uint8_t));
                if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) {
                    alt_len += 1;
                    var->vars[i].alt_bases[var->vars[i].n_alt_allele][0] = nst_nt4_table[(int)ref_seq[var->vars[i].pos-ref_beg]];
                    for (int j = 1; j < alt_len; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j] = cand_vars[cand_i].alt_seqs[hap_alle-1][j-1];
                    }
                } else {
                    for (int j = 0; j < alt_len; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j] = cand_vars[cand_i].alt_seqs[hap_alle-1][j];
                    }
                }
                var->vars[i].alt_len[var->vars[i].n_alt_allele] = alt_len;
                var->vars[i].GT[hap-1] = ++var->vars[i].n_alt_allele;
            } else var->vars[i].GT[hap-1] = 0;
        }
        var->vars[i].QUAL = 0;
        // var->vars[i].total_depth = cand_snps[i].n_depth;
        // var->vars[i].depths[0] = cand_snps[i].base_covs[cand_snps[i].base_to_i[cand_snps[i].ref_base]];
        // var->vars[i].depths[1] = cand_snps[i].n_depth - var->vars[i].depths[0];
        // var->vars[i].genotype[0] = 0; var->vars[i].genotype[1] = 0;
        i++;
    }
    free(ovlp_b);
    return(var->n = i);
}

// 1. check if haplotype of variants in bam_chunk is inconsistent with the previous bam_chunk
// 2. extend phase blocks if possible (e.g., if ≥ 1 read supports the longer phase block)
int flip_variant_hap(bam_chunk_t *prev_chunk, bam_chunk_t *cur_chunk) {
    int n_ovlp_reads = cur_chunk->n_up_ovlp_reads; if (n_ovlp_reads <= 0) return 0;
    if (prev_chunk->n_cand_vars <= 0 || cur_chunk->n_cand_vars <= 0) return 0;
    // 1) find overlapping reads that ovlp with both prev and cur variants
    int *ovlp_read_i = cur_chunk->up_ovlp_read_i; // read_i in the previous bam_chunk
    int *used_ovlp_read_i = (int*)calloc(n_ovlp_reads, sizeof(int)); int n_used_ovlp_reads = 0;
    // hts_pos_t ovlp_start_pos = -1;
    hts_pos_t ovlp_end_pos = 0;
    // int pre_var_start_i=0, pre_var_end_i=prev_chunk->n_cand_vars-1;
    int cur_var_start_i=0, cur_var_end_i=cur_chunk->n_cand_vars-1; // variants covered by overlapping reads
    // hts_pos_t pre_var_end_pos = prev_chunk->cand_vars[pre_var_end_i].pos;
    hts_pos_t cur_var_start_pos = cur_chunk->cand_vars[cur_var_start_i].pos;
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        if (prev_chunk->is_skipped[read_i]) continue;
        // hts_pos_t start_pos = prev_chunk->reads[read_i]->core.pos+1;
        hts_pos_t end_pos = bam_endpos(prev_chunk->reads[read_i]);
        // if (start_pos > pre_var_end_pos) continue;
        if (end_pos < cur_var_start_pos) continue;
        // if (ovlp_start_pos == -1) ovlp_start_pos = start_pos;
        if (end_pos > ovlp_end_pos) ovlp_end_pos = end_pos;
        used_ovlp_read_i[i] = 1; n_used_ovlp_reads++;
    }
    // collect ovlp variants from the previous bam_chunk
    // for (int i = pre_var_end_i; i >= 0; --i) {
    //     cand_var_t *var = prev_chunk->cand_vars + i;
    //    if (var->pos >= ovlp_start_pos) pre_var_start_i = i; else break;
    // }
    // collect ovlp variants from the current bam_chunk
    for (int i = cur_var_start_i; i < cur_chunk->n_cand_vars; ++i) {
        cand_var_t *var = cur_chunk->cand_vars + i;
        if (var->pos <= ovlp_end_pos)
            cur_var_end_i = i;
        else
            break;
    }
    // 2) collect read_var_profile for the overlapping reads & overlapping variants
    // 3) check if haplotype needs to be flipped
    //    var_hap == read_hap: +1
    //    var_hap != read_hap: -1
    //    var_hap == 0: 0
    int keep_hap_score = 0;
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        if (prev_chunk->is_skipped[read_i] || used_ovlp_read_i[i] == 0)
            continue;
        int read_hap = prev_chunk->haps[read_i];
        if (read_hap == 0)
            continue;
        // bam1_t *read = prev_chunk->reads[read_i];
        // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        // int pre_var_start_i0 = prev_chunk->read_var_profile[read_i].start_var_idx;
        // int pre_var_end_i0 = prev_chunk->read_var_profile[read_i].end_var_idx;
        // for (int j = pre_var_start_i; j <= pre_var_end_i; ++j) {
        // if (j < pre_var_start_i0 || j > pre_var_end_i0 || prev_chunk->read_var_profile[read_i].var_is_used[j-pre_var_start_i0] == 0) continue;
        // int allele_i = prev_chunk->read_var_profile[read_i].alleles[j-pre_var_start_i0];
        // fprintf(stderr, "\tPre-Var: %lld, %c, allele: %d\t", (long long) prev_chunk->cand_vars[j].pos, BAM_CIGAR_STR[prev_chunk->cand_vars[j].var_type], allele_i);
        // cand_var_t *var = prev_chunk->cand_vars + j;
        // if (var->hap_to_cons_alle[1] == allele_i) fprintf(stderr, "H1\n");
        // else if (var->hap_to_cons_alle[2] == allele_i) fprintf(stderr, "H2\n");
        // else fprintf(stderr, "none\n");
        // }
        int cur_var_start_i0 = cur_chunk->read_var_profile[i].start_var_idx;
        int cur_var_end_i0 = cur_chunk->read_var_profile[i].end_var_idx;
        for (int j = cur_var_start_i; j <= cur_var_end_i; ++j)
        {
            if (j < cur_var_start_i0 || j > cur_var_end_i0 || cur_chunk->read_var_profile[i].var_is_used[j - cur_var_start_i0] == 0)
                continue;
            int allele_i = cur_chunk->read_var_profile[i].alleles[j - cur_var_start_i0];
            // fprintf(stderr, "\tCur-Var: %lld, %c, allele: %d\n", (long long) cur_chunk->cand_vars[j].pos, BAM_CIGAR_STR[cur_chunk->cand_vars[j].var_type], allele_i);
            cand_var_t *var = cur_chunk->cand_vars + j;
            if (var->hap_to_cons_alle[read_hap] == allele_i)
                keep_hap_score++;
            else if (var->hap_to_cons_alle[3 - read_hap] == allele_i)
                keep_hap_score--;
        }
    }
    free(used_ovlp_read_i);
    // fprintf(stderr, "pos: %" PRId64 ", keep_hap_score: %d\n", prev_chunk->ref_end, keep_hap_score);
    int flip = (keep_hap_score > 0 ? 0 : 1);
    // 4) update HP tag for the overlapping reads (in prev_bam_chunk, not cur_bam_chunk, to maintain the order in output bam file)
    // since only part of the variants (either prev or cur variants) are used to determine the HP in previous step
    // XXX currently only use HP & variant from the cur_bam_chunk
    // char test_read_name[1024] = "m84039_231005_222902_s1/234751166/ccs";
    for (int i = 0; i < n_ovlp_reads; ++i) {
        int read_i = ovlp_read_i[i];
        // bam1_t *read = prev_chunk->reads[read_i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
        // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        if (prev_chunk->haps[read_i] != 0 || cur_chunk->haps[i] == 0)
            continue;
        if (flip)
            prev_chunk->haps[read_i] = cur_chunk->haps[i] ^ 3;
        else
            prev_chunk->haps[read_i] = cur_chunk->haps[i];
    }
    // 5) update phase block for the current bam_chunk
    return flip;
}

// update phase set within bam_chunk
int extend_phase_set(bam_chunk_t *chunk) {
    if (chunk->n_cand_vars <= 1)
        return 0;
    hts_pos_t phase_set = -1;
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
    }
return 0;
}

void print_var_seqs(digar1_t *digars, uint8_t **reg_var_seqs, int reg_n_digar, FILE *fp) {
    fprintf(fp, "pos\ttype\tlen\tqi\tbase\tis_low_qual\n");
    for (int i = 0; i < reg_n_digar; ++i) {
        fprintf(fp, "%" PRId64 "\t%c\t%d\t%d\t%d", digars[i].pos, BAM_CIGAR_STR[digars[i].type], digars[i].len, digars[i].qi, digars[i].is_low_qual);
        if (reg_var_seqs[i]) {
            int seq_len = digars[i].len;
            if (digars[i].type == BAM_CDIFF) seq_len *= 2;
            fprintf(fp, "\t");
            for (int j = 0; j < seq_len; ++j) {
                fprintf(fp, "%c", "ACGT"[reg_var_seqs[i][j]]);
            }
        }
        fprintf(fp, "\n");
    }
}

int update_read_var_profile_with_phased_info(bam_chunk_t *chunk, int target_var_cate) {
    return 0;
}

int *collect_noisy_read_haps(bam_chunk_t *chunk, int noisy_reg_i) {
    int n_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    int *haps = (int*)calloc(n_reads, sizeof(int));
    for (int j = 0; j < chunk->noisy_reg_to_n_reads[noisy_reg_i]; ++j) {
        int read_i = chunk->noisy_reg_to_reads[noisy_reg_i][j];
        haps[j] = chunk->haps[read_i];
    }
    return haps;
}

int **collect_read_xid_profile(bam_chunk_t *chunk, int noisy_reg_i) {
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    hts_pos_t reg_beg = cr_start(chunk->chunk_noisy_regs, noisy_reg_i), reg_end = cr_end(chunk->chunk_noisy_regs, noisy_reg_i);
    int n_reads = chunk->noisy_reg_to_n_reads[noisy_reg_i];
    int **xid_profile = (int**)malloc(n_reads * sizeof(int*));
    for (int j = 0; j < chunk->noisy_reg_to_n_reads[noisy_reg_i]; ++j) {
        xid_profile[j] = (int*)calloc(20, sizeof(int));
        int read_i = chunk->noisy_reg_to_reads[noisy_reg_i][j];
        digar1_t *reg_digars = (digar1_t*)malloc(chunk->digars[read_i].m_digar * sizeof(digar1_t)); int fully_cover = 0;
        uint8_t **reg_var_seqs = (uint8_t**)malloc(chunk->digars[read_i].m_digar * sizeof(uint8_t*));
        int reg_n_digar = collect_reg_digars_var_seqs(chunk, read_i, reg_beg, reg_end, reg_digars, reg_var_seqs, &fully_cover);
        if (!fully_cover) {
            xid_profile[j][0] = -1;
            goto end_loop; // only collect reads that fully cover the noisy region
        }
        for (int k = 0; k < reg_n_digar; ++k) {
            if (reg_digars[k].type == BAM_CEQUAL) continue;
            else if (reg_digars[k].type == BAM_CDIFF) {
                for (int l = 0; l < reg_digars[k].len; ++l) {
                    xid_profile[j][reg_var_seqs[k][l*2]] += 1;
                    xid_profile[j][5+reg_var_seqs[k][l*2+1]] += 1;
                }
            } else if (reg_digars[k].type == BAM_CINS) {
                for (int l = 0; l < reg_digars[k].len; ++l)
                    xid_profile[j][10+reg_var_seqs[k][l]]++;
            } else if (reg_digars[k].type == BAM_CDEL) {
                for (int l = 0; l < reg_digars[k].len; ++l)
                    xid_profile[j][15+reg_var_seqs[k][l]]++;
            }
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "Read: %s\n", bam_get_qname(chunk->reads[read_i]));
            // print XID profile
            for (int k = 0; k < 20; ++k) {
                fprintf(stderr, "%d\t", xid_profile[j][k]);
            } fprintf(stderr, "\n");
            print_var_seqs(reg_digars, reg_var_seqs, reg_n_digar, stderr);
        }
    end_loop:
        free(reg_digars);
        for (int k = 0; k < reg_n_digar; ++k) {
            if (reg_var_seqs[k]) free(reg_var_seqs[k]);
        } free(reg_var_seqs);
    }
    return xid_profile;
}


int make_vars_from_mediod_reads(bam_chunk_t *chunk, int noisy_reg_i, int *medoids, int n_mediods, cand_var_t **cand_vars) {
    return 0;
}

read_var_profile_t *collect_noisy_read_var_profile(bam_chunk_t *chunk, int noisy_reg_i, int n_cand_vars, cand_var_t *cand_vars) {
    return NULL;
}

// excluding non-fully covering reads: not long enough, or clipped
// keep up to 2 medoid reads for each haplotype, in case the phasing was wrong
int *collect_medoid_noisy_reads(bam_chunk_t *chunk, int noisy_reg_i, int *n_medoids, int use_phase_info) {
    // check if alignment-based medoid selection method is needed
    // if not, use XID-profile-based medoid selection method
    return NULL;
}

// phase info is used as initial cluster assignment for noisy region variants
int update_read_var_profile(bam_chunk_t *chunk, int target_var_cate, int use_phase_info) {
    // 1) update read_var_profile for target_var_cate using re-alignment
    // merge_replace

    // 2) collect candidate variants in noisy region, insert into cand_vars, then collect read_var_profile
    // merge_insert
    int n_vars = 0;
    cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
    for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
        if (cr_label(noisy_regs, i) < 10) continue;
        // collect candicate variants in noisy region based on XID-profile
        hts_pos_t reg_beg = cr_start(noisy_regs, i), reg_end = cr_end(noisy_regs, i);
        // collect medoid reads and candidate genotype sequences for each noisy region
        int n_medoids = 0;
        int *medoid_read_idxs = collect_medoid_noisy_reads(chunk, i, &n_medoids, use_phase_info);
        cand_var_t *cand_vars = NULL;
        int n_cand_vars = make_vars_from_mediod_reads(chunk, i, medoids, n_medoids, &cand_vars);
        // free XID-profile & medoids
        if (medoids!= NULL) free(medoids);
        
        read_var_profile_t *p = collect_noisy_read_var_profile(chunk, i, n_cand_vars, cand_vars);
        // insert cand_vars into chunk->cand_vars
        // insert p into chunk->read_var_profile
    }
    return 0;
}

int old_update_read_var_profile(bam_chunk_t *chunk, int target_var_cate, int use_phase_info) {
    // 1) update read_var_profile for target_var_cate using re-alignment
    // merge_replace

    // 2) collect candidate variants in noisy region, insert into cand_vars, then collect read_var_profile
    // merge_insert
    int n_vars = 0;
    for (int i = 0; i < chunk->chunk_noisy_regs->n_r; ++i) {
        cgranges_t *noisy_regs = chunk->chunk_noisy_regs;
        if (cr_label(noisy_regs, i) < 10) continue;
        // collect candicate variants in noisy region based on XID-profile
        hts_pos_t reg_beg = cr_start(noisy_regs, i), reg_end = cr_end(noisy_regs, i);
        // collect candidate genotype sequences for each noisy region
        // xid_profile[read_i][0]==-1: read_i does not cover the noisy region, skip in k-medoids
        int **xid_profile = collect_read_xid_profile(chunk, i);
        int n_reads = chunk->noisy_reg_to_n_reads[i];
        int *haps = collect_noisy_read_haps(chunk, i);
        int n_medoids = 3;
        int *medoids = xid_profile_2medoids(xid_profile, 20, n_medoids, haps, n_reads, use_phase_info); // 20 = RefX, ReadX, Ins, Del * (ACGTN)
        int tmp_hp_tag = -1;
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "NoisyRegs: %s:%d-%d %d %d\n", chunk->tname, cr_start(noisy_regs, i), cr_end(noisy_regs, i), cr_end(noisy_regs, i) - cr_start(noisy_regs, i), cr_label(noisy_regs, i));
            fprintf(stderr, "Medoid-reads\t");
            for (int j = 0; j < n_medoids; ++j) {
                int read_i = chunk->noisy_reg_to_reads[i][medoids[j]];
                int hp_tag = get_aux_int_from_bam(chunk->reads[read_i], "HP");
                fprintf(stderr, "%s HP:%d\t", bam_get_qname(chunk->reads[read_i]), hp_tag);
            } fprintf(stderr, "\n");
        }
        cand_var_t *cand_vars = NULL;
        int n_cand_vars = make_vars_from_mediod_reads(chunk, i, medoids, n_medoids, &cand_vars);
        // free XID-profile & medoids
        if (medoids!= NULL) free(medoids);
        for (int j = 0; j < chunk->noisy_reg_to_n_reads[i]; ++j) free(xid_profile[j]); free(xid_profile);
        
        // read_var_profile_t *p = init_read_var_profile(chunk->n_reads, n_cand_vars);

        read_var_profile_t *p = collect_noisy_read_var_profile(chunk, i, n_cand_vars, cand_vars);
        // insert cand_vars into chunk->cand_vars
        // insert p into chunk->read_var_profile
        free(haps);
    }
    return 0;
}

void collect_var_main(const call_var_pl_t *pl, bam_chunk_t *chunk, var_t *var) {
    // first round: easy-to-call SNPs (+indels)
    // 1. collect X/I/D sites from BAM
    collect_digars_from_bam(chunk, pl);

    // 2. merge all var sites from all reads, including low-depth ones, but not including nosiy-region ones
    var_site_t *var_sites = NULL; int n_var_sites;
    if ((n_var_sites = collect_all_cand_var_sites(chunk, &var_sites)) <= 0) return;

    // 3. collect all candidate variants, not including noisy-region ones
    // collect reference and alternative alleles for all var sites
    // all cand vars, including true/false germline/somatic variants
    // XXX for noisy/repeat regions, we need to carefully pick the candidate variants, based on supporting counts & re-alignments
    // so collect_cand_vars and classify_cand_vars should be run simultaneously
    collect_cand_vars(chunk, n_var_sites, var_sites); free(var_sites);

    // 4. filter out low-depth ones;
    //    idenitfy repeat region;
    if (classify_cand_vars(chunk, n_var_sites, pl->opt) <= 0) return;
    // collect all read_var_profile, including repeat/noisy regions

    // 5. collect read-wise var profiles
    chunk->read_var_profile = collect_read_var_profile(chunk);

    // process candidate variants in the following order
    // within each round, use previously obtained phasing/haplotype information as intialization
    //   1. LONGCALLD_EASY_HET_SNP | LONGCALLD_EASY_HET_INDEL
    //   2. LONGCALLD_REP_HET_VAR
    //   3. LONGCALLD_DENSE_REG_VAR + ALL VARS IN NOISY_REG
    //   4. LONGCALLD_CAND_HOM_VAR
    //   5. LONGCALLD_CAND_SOMATIC_VAR
    // 1st round: easy-to-call het.
    //   0) easy-to-call germline het. variants here to determine the haplotype
    // Both SNP & INDEL are considered, we can use only SNP, need to test the performance difference
    // 6. co-phasing and variant calling using easy-to-call SNPs (+ indels)
    assign_hap_based_on_het_vars(chunk, LONGCALLD_EASY_HET_SNP | LONGCALLD_EASY_HET_INDEL, pl->opt);
    // co-update cand_vars & read_var_profile based on the assigned haplotype for hard-to-call regions
    // for LONGCALLD_REP_HET_VAR, update support read count (read_var_profile) based on re-alignment
    // for noisy regions, collect candidate variants and insert into cand_vars, then collect read_var_profile
    // second-round: update variant & profile based on re-alignment
    // 1. clipping around large indels & noisy region
    // 2. repeat-region
    // 3. noisy-region
    update_read_var_profile(chunk, LONGCALLD_REP_HET_VAR | LONGCALLD_DENSE_REG_VAR, 1);
    // 2nd round of co-phasing & haplotype assignment using all variants
    // assign_hap_based_on_all_vars()
    // variant calling using haplotype & supporting read counts
    // 5th round: extend phase block within bam_chunk
}

// stitch ii and ii+1
void stitch_var_main(const call_var_step_t *step, bam_chunk_t *chunk, var_t *var, long ii) {
    const call_var_pl_t *pl = step->pl;
    // ref_seq_t *ref_seq = pl->ref_seq;
    // if (ii == 0) { // extend phase set between two adjacent bam chunks
        // bam_chunk_t *prev_bam_chunk = step->chunks+ii-1;
        // chunk->flip_hap = prev_chunk->flip_hap ^ flip_variant_hap(prev_bam_chunk, bam_chunk);
    // }

    // generate variant-related information, e.g., GT, DP, etc.
    // merge variants based on ref_pos if needed
    make_variants(chunk, var);
}