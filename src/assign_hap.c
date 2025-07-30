#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include "align.h"
#include "sdust.h"
#include "math_utils.h"
#include <inttypes.h>
#include <stdlib.h>

extern int LONGCALLD_VERBOSE;

// init reads' hap & ps
void read_init_hap_phase_set(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->n_reads; ++i) {
        chunk->haps[i] = 0; chunk->phase_sets[i] = -1; chunk->phase_scores[i] = 0;
    }
}

static int get_var_init_max_cov_allele(int is_ont, cand_var_t *var) {
    // XXX do this for hifi as well ???
    if (is_ont == 1 && var->is_homopolymer_indel) return -1; // not confident
    int max_cov = 0, max_cov_alle_i = -1, total_cov = 0;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        total_cov += var->alle_covs[i];
        if (var->alle_covs[i] > max_cov) { // prefer ref allele
            max_cov = var->alle_covs[i];
            max_cov_alle_i = i;
        } 
    }
    return max_cov_alle_i;
}

// all candidate vars, including heterozygous and homozygous SNPs:
// hap_to_alle_profile: HAP 0/1/2: -> allele: 0/1 -> count: 0
// hap_to_cons_alle: HAP 0/1/2: HOM -> -1/1/1, HET -> -1/-1/-1
void var_init_hap_profile_cons_allele(const call_var_opt_t *opt, cand_var_t *cand_vars, int *var_idx, int n_cand_vars, int *var_i_to_cate) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
            var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
            var->hap_to_cons_alle[0] = get_var_init_max_cov_allele(opt->is_ont, var); // hom_idx
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        } else {
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j)
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            var->hap_to_cons_alle[0] = get_var_init_max_cov_allele(opt->is_ont, var); // hom_idx
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

void var_init_hap_profile_cons_allele1(const call_var_opt_t *opt, cand_var_t *var) {
    if (var->hap_to_alle_profile == NULL) {
        var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
        for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
        var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
        var->hap_to_cons_alle[0] = get_var_init_max_cov_allele(opt->is_ont, var); // hom_idx
        for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
    } else {
        for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
        var->hap_to_cons_alle[0] = get_var_init_max_cov_allele(opt->is_ont, var); // hom_idx
        for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
    }
}


void var_init_hap_to_alle_profile(cand_var_t *cand_vars, int *var_idx, int n_cand_vars) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
        } else {
            for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            }
        }
    }
}

int select_init_var(cand_var_t *cand_vars, int *var_idx, int n_cand_vars, int *var_i_to_cate) {
    int clean_het_snp_i = -1, clean_het_indel_i = -1, noisy_het_snp_i = -1, noisy_het_indel_i = -1;
    int clean_het_snp_depth = 0, clean_het_indel_depth = 0, noisy_het_snp_depth = 0, noisy_het_indel_depth = 0;
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars + var_i;
        if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_SNP) {
            if (clean_het_snp_i == -1 || clean_het_snp_depth < var->total_cov) { 
                clean_het_snp_i = _var_i; clean_het_snp_depth = var->total_cov;
            }
        } else if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_INDEL) {
            if (clean_het_indel_i == -1 || clean_het_indel_depth < var->total_cov) { 
                clean_het_indel_i = _var_i; clean_het_indel_depth = var->total_cov;
            }
        } else if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR) {
            if (cand_vars[var_i].var_type == BAM_CDIFF) {
                if (noisy_het_snp_i == -1 || noisy_het_snp_depth < var->total_cov) {
                    noisy_het_snp_i = _var_i; noisy_het_snp_depth = var->total_cov;
                }
            } else if (var->is_homopolymer_indel == 0) { // not homopolymer indel
                if (noisy_het_indel_i == -1 || noisy_het_indel_depth < var->total_cov) {
                    noisy_het_indel_i = _var_i; noisy_het_indel_depth = var->total_cov;
                }
            }
        }
    }
    if (clean_het_snp_i != -1) return clean_het_snp_i;
    else if (clean_het_indel_i != -1) return clean_het_indel_i;
    else if (noisy_het_snp_i != -1) return noisy_het_snp_i;
    else if (noisy_het_indel_i != -1) return noisy_het_indel_i;
    else return -1; // no het variants
}

int read_to_cons_allele_score(int read_i, int hap, cand_var_t *var, int var_cate, int allele_i) {
    assert(hap == 1 || hap == 2);
    int var_score = 1;
    if (var_cate == LONGCALLD_CLEAN_HET_SNP) var_score = 2;
    else if (var_cate == LONGCALLD_CLEAN_HET_INDEL) var_score = 2;
        // if (var->var_type == BAM_CDIFF) var_score = 4;
        // else var_score = 3;
    // } else if (var_cate == LONGCALLD_NOISY_CAND_HET_VAR) {
        // var_score = 1;
        // if (var->var_type == BAM_CDIFF) var_score = 2;
        // else var_score = 1;
    // }
    if (var->hap_to_cons_alle[hap] == -1 && var->hap_to_cons_alle[3-hap] == -1) return 0;
    else {
        if (var->hap_to_cons_alle[hap] == -1) var->hap_to_cons_alle[hap] = 1-var->hap_to_cons_alle[3-hap];
        if (var->hap_to_cons_alle[3-hap] == -1) var->hap_to_cons_alle[3-hap] = 1-var->hap_to_cons_alle[hap];
    }
    if (var->hap_to_cons_alle[hap] == allele_i) return var_score;
    else if (var->hap_to_cons_alle[hap] == -1) return 0;
    else return -var_score;
}

// int assign_read_hap_based_on_pre_reads(int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
// no check for reject counts, as this is the initialization step
int init_assign_read_hap_based_on_cons_alle(bam_chunk_t *chunk, int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int *hap_scores = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int n_vars_used[3] = {0, 0, 0};

    chunk->n_clean_agree_vars[read_i] = chunk->n_clean_conflict_vars[read_i] = 0;

    int n_clean_agree_vars[3] = {0, 0, 0}; // for clean vars
    int n_clean_conflict_vars[3] = {0, 0, 0}; // for clean vars
    for (int var_i = p[read_i].start_var_idx; var_i <= p[read_i].end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = cand_vars+var_i;
        if (var->is_homopolymer_indel == 1 || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) continue; 
        int read_var_idx = var_i - p[read_i].start_var_idx;
        if (p[read_i].alleles[read_var_idx] < 0) continue;
        int score0;
        for (int hap = 1; hap <= 2; ++hap) {
            score0 = read_to_cons_allele_score(read_i, hap, var, var_i_to_cate[var_i], p[read_i].alleles[read_var_idx]);
            if (score0 != 0) {
                n_vars_used[hap]++;
                if (score0 > 0) {
                    if (!!(var_i_to_cate[var_i] & LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE)) {
                        n_clean_agree_vars[hap]++;
                    }
                } else {
                    if (!!(var_i_to_cate[var_i] & LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE)) {
                        n_clean_conflict_vars[hap]++;
                    }
                }
            }
            hap_scores[hap] += score0;
        }
    }
    int max_hap = 0, max_score = 0, min_hap = 0, min_score = 0;
    for (int hap = 1; hap <= 2; ++hap) {
        if (hap_scores[hap] > max_score) {
            max_hap = hap; max_score = hap_scores[hap];
        } else if (hap_scores[hap] < min_score) {
            min_hap = hap; min_score = hap_scores[hap];
        }
    }
    free(hap_scores);
    if (n_vars_used[1] == 0 && n_vars_used[2] == 0) return -1; // no used vars
    else if (max_score == 0 && min_score == 0) return 0; // tied
    else if (max_score > 0) {
        chunk->n_clean_agree_vars[read_i] = n_clean_agree_vars[max_hap]; chunk->n_clean_conflict_vars[read_i] = n_clean_conflict_vars[max_hap];
        return max_hap;
    } else return 3-min_hap;
    // return max_hap;
}

// Assign read haplotype based on consensus allele of candidate variants
// check for reject counts, reject the haplotype has too many clean rejects (>=2)
int assign_read_hap_based_on_cons_alle(int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int *hap_scores = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int n_vars_used[3] = {0, 0, 0};

    int n_clean_reject[2] = {0, 0};
    for (int var_i = p[read_i].start_var_idx; var_i <= p[read_i].end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = cand_vars+var_i;
        if (var->is_homopolymer_indel == 1 || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) continue; 
        int read_var_idx = var_i - p[read_i].start_var_idx;
        if (p[read_i].alleles[read_var_idx] < 0) continue;
        int score0;
        for (int hap = 1; hap <= 2; ++hap) {
            score0 = read_to_cons_allele_score(read_i, hap, var, var_i_to_cate[var_i], p[read_i].alleles[read_var_idx]);
            if (score0 != 0) n_vars_used[hap]++;
            if (score0 < 0 && !!(var_i_to_cate[var_i] & LONGCALLD_CAND_GERMLINE_CLEAN_VAR_CATE)) n_clean_reject[hap-1] += 1;
            hap_scores[hap] += score0;
        }
    }
    int max_hap = 0, max_score = 0, min_hap = 0, min_score = 0;
    for (int hap = 1; hap <= 2; ++hap) {
        if (hap_scores[hap] > max_score) {
            max_hap = hap; max_score = hap_scores[hap];
        } else if (hap_scores[hap] < min_score) {
            min_hap = hap; min_score = hap_scores[hap];
        }
    }
    free(hap_scores);
    if (n_vars_used[1] == 0 && n_vars_used[2] == 0) return -1; // no used vars
    else if (max_score == 0 && min_score == 0) return 0; // tied
    else if (max_score > 0) {
        // fprintf(stderr, "Error, should not go to this line. n_clean_reject: %d,%d, hap_scores: %d,%d, max_hap: %d, %d\n", n_clean_reject[0], n_clean_reject[1], hap_scores[1], hap_scores[2], max_hap, max_score);
        // XXX for ONT this is not good!!!
        if (n_clean_reject[max_hap-1] >= 2) return -1; // reject if both haplotype have too many clean rejects
        else return max_hap;
    } else {
        // fprintf(stderr, "Error, should not go to this line. n_clean_reject: %d,%d, hap_scores: %d,%d, min_hap: %d, %d\n", n_clean_reject[0], n_clean_reject[1], hap_scores[1], hap_scores[2], min_hap, min_score);
        if (n_clean_reject[3-min_hap-1] >= 2) return -1; // reject if both haplotype have too many clean rejects
        else return 3-min_hap;
    }
}

int update_var_hap_to_cons_alle(const call_var_opt_t *opt, cand_var_t *var, int var_cate, int hap) {
    if (hap == 0) return 0;
    int max_cov = 0, max_cov_alle_i = -1;
    // int n_max_alt = 0;
    // if (var_cate & (LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_NOISY_CAND_HET_VAR)) n_max_alt = 1;
    // else if (var_cate & (LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HOM_VAR)) n_max_alt = 2;

    int total_cov = 0;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        total_cov += var->hap_to_alle_profile[hap][i];
        if (var->hap_to_alle_profile[hap][i] > max_cov) { // prefer ref allele
            max_cov = var->hap_to_alle_profile[hap][i];
            max_cov_alle_i = i;
        } 
        // else if (var->hap_to_alle_profile[hap][i] == max_cov && i != 0) { // alt allele
        //     if (n_max_alt > 0) {
        //         max_cov_alle_i = i; n_max_alt--;
        //     }
        // }
    }
    if (opt->is_ont && var->is_homopolymer_indel == 1 && max_cov < total_cov * 0.67) max_cov_alle_i = -1; // not confident
    // fprintf(stderr, "%ld: hap: %d, max_allele_i: %d, total_cov: %d, max_cov: %d, var_cate: %d\n", var->pos, hap, max_cov_alle_i, total_cov, max_cov, var_cate);
    var->hap_to_cons_alle[hap] = max_cov_alle_i;
    return 0;
}

int update_var_hap_profile_cons_alle_based_on_read_hap(const call_var_opt_t *opt, int read_i, int hap, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "var_i: %d, %" PRIi64 " allele_i: %d\n", var_i, cand_vars[var_i].pos, allele_i);
        if (allele_i < 0) continue;
        cand_var_t *var = cand_vars+var_i; int var_cate = var_i_to_cate[var_i];
        if (hap == 0) {
            for (int i = 1; i <= 2; ++i) {
                var->hap_to_alle_profile[i][allele_i] += 1;
                update_var_hap_to_cons_alle(opt, var, var_cate, i);
            }
        } else {
            var->hap_to_alle_profile[hap][allele_i] += 1;
            update_var_hap_to_cons_alle(opt, var, var_cate, hap);
        }
    }
    return 0;
}

int update_var_hap_profile_based_on_read_hap(int read_i, int hap, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (allele_i < 0) continue;
        if (hap == 0) {
            cand_vars[var_i].hap_to_alle_profile[1][allele_i] += 1;
            cand_vars[var_i].hap_to_alle_profile[2][allele_i] += 1;
        } else cand_vars[var_i].hap_to_alle_profile[hap][allele_i] += 1;
    }
    return 0;
}

int check_agree_haps(int read_i, int hap, cand_var_t *vars, int var1, int var2, read_var_profile_t *p) {
    if (var1 < p[read_i].start_var_idx || var2 > p[read_i].end_var_idx) return -1;
    if (hap == 0) return -1; // always agree for hap==0
    int allele_i1 = p[read_i].alleles[var1 - p[read_i].start_var_idx];
    int allele_i2 = p[read_i].alleles[var2 - p[read_i].start_var_idx];
    if (allele_i1 < 0 || allele_i2 < 0) return -1;

    int agree = 0, conflict = 0, agree_hap = hap, conflict_hap = 3-hap;
    if (vars[var1].hap_to_cons_alle[agree_hap] == allele_i1 && vars[var2].hap_to_cons_alle[agree_hap] == allele_i2) agree = 1;
    if (vars[var1].hap_to_cons_alle[hap] == allele_i1 && vars[var2].hap_to_cons_alle[conflict_hap] == allele_i2) conflict = 1;

    if (agree) return 1; else if (conflict) return 0;
    else return -1;
}

void update_read_phase_set(bam_chunk_t *chunk, int *var_is_valid, read_var_profile_t *p, cand_var_t *cand_vars) {
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i]) continue;
        if (p[read_i].start_var_idx == -1) continue;
        hts_pos_t phase_set = -1;
        for (int var_i = p[read_i].start_var_idx; var_i <= p[read_i].end_var_idx; ++var_i) {
            if (var_is_valid[var_i] == 0) continue;
            cand_var_t *var = cand_vars+var_i;
            // check if var is het
            if (var->hap_to_cons_alle[1] != -1 && var->hap_to_cons_alle[2] != -1 && var->hap_to_cons_alle[1] != var->hap_to_cons_alle[2]) {
                phase_set = var->phase_set;
            }
            if (phase_set != -1) break;
        }
        chunk->phase_sets[read_i] = phase_set;
    }
}

// record n_agree/n_conflit for each read as well
// XXX assign each read a HaplotypeScore based on the consensus allele of candidate variants
// could be used for somatic variant calling
// reads->haps & haps_to_cons_alle are UpToDate
int iter_update_var_hap_cons_phase_set(bam_chunk_t *chunk, int *var_idx, read_var_profile_t *p, cand_var_t *cand_vars, int n_cand_vars, int *var_i_to_cate) {
    int *het_var_idx = (int*)malloc(n_cand_vars * sizeof(int));
    int n_het_vars = 0;
    int *is_het = (int*)calloc(n_cand_vars, sizeof(int));
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars+var_i;
        // XXX only use clean het vars, do not use homopolymer indels
        // if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR) continue;
        if (var->hap_to_cons_alle[1] != -1 && var->hap_to_cons_alle[2] != -1 
            && var->hap_to_cons_alle[1] != var->hap_to_cons_alle[2] && var->is_homopolymer_indel == 0) {
            is_het[_var_i] = 1;
            het_var_idx[n_het_vars++] = _var_i;
        }
    }

    cgranges_t *read_var_cr = chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    int *n_agree = (int*)calloc(n_cand_vars, sizeof(int)); // number of reads supporting current haplotype for each var
    int *n_conflict = (int*)calloc(n_cand_vars, sizeof(int)); // number of reads conflicting current haplotype for each var
    // for each two adjacent het vars, check if they are in the same phase set
    // for (int _var_i = 1; _var_i < n_cand_vars; ++_var_i) {
        // int var_i = var_idx[_var_i];
    for (int __var_i = 1; __var_i < n_het_vars; ++__var_i) {
        int _var_i = het_var_idx[__var_i];
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars+var_i;
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "var_i %d %" PRIi64 " %d-%c-%d %d reads\n", var_i, var->pos, var->ref_len, BAM_CIGAR_STR[cand_vars[var_i].var_type], cand_vars[var_i].alt_len, cand_vars[var_i].total_cov);
        ovlp_n = cr_overlap(read_var_cr, "cr", var_idx[het_var_idx[__var_i-1]], var_i+1, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) { // assign each read a haplotype if it is not assigned yet
            int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
            if (chunk->is_skipped[read_i]) continue;
            // XXX only use var with valid target_var_cate
            int agree = check_agree_haps(read_i, chunk->haps[read_i], cand_vars, var_idx[het_var_idx[__var_i-1]], var_i, p);
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read: %s agree: %d\n", bam_get_qname(chunk->reads[read_i]), agree);
            if (agree > 0) n_agree[_var_i]++;
            else if (agree == 0) n_conflict[_var_i]++;
        }
    }
    // new phase set if n_agree < THRESHOLD & n_conflict < THRESHOLD
    // flip hap_to_cons_alle if n_agree < n_conflict
    int changed = 0;
    int flip = 0; hts_pos_t phase_set = -1;
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars+var_i;
        if (_var_i == 0) {
            if (var->var_type == BAM_CDIFF) phase_set = var->pos;
            else phase_set = var->pos - 1;
            var->phase_set = phase_set;
            continue;
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%" PRIi64 " %d %d\n", var->pos, n_agree[_var_i], n_conflict[_var_i]);
        if (is_het[_var_i] == 1) {
            if (n_agree[_var_i] < 2 && n_conflict[_var_i] < 2) { // new phase set
                if (var->var_type == BAM_CDIFF) phase_set = var->pos;
                else phase_set = var->pos - 1;
                // flip ^= 0;
            } else if (n_conflict[_var_i] > n_agree[_var_i]) {
                flip ^= 1;
            } 
            // else if (n_agree[_var_i] > n_conflict[_var_i]) {
                // flip ^= 0;
            // }
            if (flip == 1) { // flip hap_to_cons_alle
                changed = 1;
                for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
                    int tmp = var->hap_to_cons_alle[hap];
                    var->hap_to_cons_alle[hap] = var->hap_to_cons_alle[3-hap];
                    var->hap_to_cons_alle[3-hap] = tmp;
                }
            }
        }
        var->phase_set = phase_set;
    }
    free(ovlp_b); free(n_agree); free(n_conflict); free(het_var_idx); free(is_het);
    return changed;
}

// update read's haps & var's hap_to_cons_alle
int iter_update_var_hap_to_cons_alle(const call_var_opt_t *opt, bam_chunk_t *chunk, int *var_idx, int n_cand_vars, read_var_profile_t *p, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate) {
    int **cur_hap_to_cons_alle = (int**)malloc(n_cand_vars * sizeof(int*));
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        cur_hap_to_cons_alle[_var_i] = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
        for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
            cur_hap_to_cons_alle[_var_i][hap] = var->hap_to_cons_alle[hap];
        }
    }
    // var-wise loop, update phase set XXX
    // update hap_to_alle_profile + hap_to_cons_alle
    var_init_hap_to_alle_profile(cand_vars, var_idx, n_cand_vars);
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i]) continue;
        int hap;
        hap = init_assign_read_hap_based_on_cons_alle(chunk, read_i, cand_vars, p, var_i_to_cate, target_var_cate);
        // if (opt->is_ont) hap = init_assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
        // else hap = assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
        if (hap == -1) hap = 0;
        chunk->haps[read_i] = hap;
        update_var_hap_profile_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
    }
    // update hap_to_cons_alle
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
            update_var_hap_to_cons_alle(opt, var, var_i_to_cate[var_i], hap);
        }
    }
    // check if any hap_to_cons_alle changed
    int changed = 0;
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
            if (var->hap_to_cons_alle[hap] != cur_hap_to_cons_alle[_var_i][hap]) {
                changed = 1; break;
            }
        }
    }
    for (int i = 0; i < n_cand_vars; ++i) free(cur_hap_to_cons_alle[i]); free(cur_hap_to_cons_alle);
    return changed;
}

// this is for germline variant calling only
// read's PS was assigned using a fake HET (HOM) var
// goal: 1) assign haplotype + phase set to all reads
//       2) variant calling based on clustered reads
int assign_hap_based_on_germline_het_vars_kmeans(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    read_var_profile_t *p = chunk->read_var_profile;
    int n_valid_vars = 0, *valid_var_idx = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    int *var_is_valid = (int*)calloc(chunk->n_cand_vars, sizeof(int));
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        if ((chunk->var_i_to_cate[i] & target_var_cate) == 0) continue;
        valid_var_idx[n_valid_vars++] = i;
        var_is_valid[i] = 1;
    }
    if (n_valid_vars == 0) {
        free(valid_var_idx); free(var_is_valid);
        return 0;
    }
    int *var_i_to_cate = chunk->var_i_to_cate;
    cand_var_t *cand_vars = chunk->cand_vars;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;

    read_init_hap_phase_set(chunk);
    var_init_hap_profile_cons_allele(opt, cand_vars, valid_var_idx, n_valid_vars, var_i_to_cate);

    // 1st loop: var-wise loop
    int init_var_i = select_init_var(cand_vars, valid_var_idx, n_valid_vars, var_i_to_cate);
    if (init_var_i != -1) {
        int *var_ii = (int*)malloc(n_valid_vars * sizeof(int));
        var_ii[0] = init_var_i;
        for (int var_i = init_var_i-1; var_i >= 0; --var_i) var_ii[init_var_i-var_i] = var_i;
        for (int var_i = init_var_i+1; var_i < n_valid_vars; ++var_i) var_ii[var_i] = var_i;

        // for each var: assign reads covering the var to HAP1 or HAP2
        for (int _var_i = 0; _var_i < n_valid_vars; ++_var_i) {
            int var_i = valid_var_idx[var_ii[_var_i]];
            // HOM vars not used in the 1st round
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) // || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR)
                continue;
            cand_var_t *var = cand_vars+var_i;
            if (LONGCALLD_VERBOSE >= 2)
                fprintf(stderr, "var_i %d %" PRIi64 " %d-%c-%d %d reads\n", var_i, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len, var->total_cov);
            ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
            for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) { // assign each read a haplotype if it is not assigned yet
                int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
                if (chunk->is_skipped[read_i] || chunk->haps[read_i] != 0) continue;
                // this round of assignment may not be correct for all reads, will be updated in the following rounds
                // 1/2: hap, 0: tied, no hap, -1: no var can be used
                int hap = init_assign_read_hap_based_on_cons_alle(chunk, read_i, cand_vars, p, var_i_to_cate, target_var_cate);
                if (hap == -1) { // no used vars, new phase set
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "new PS: %" PRIi64 " %s\n", cand_vars[var_i].pos, bam_get_qname(chunk->reads[read_i]));
                    hap = 1;
                }
                chunk->haps[read_i] = hap;
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read: %s hap: %d\n", bam_get_qname(chunk->reads[read_i]), hap);
                // update_var_hap_profile_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
                update_var_hap_profile_cons_alle_based_on_read_hap(opt, read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
            }
        } free(var_ii); free(ovlp_b);
    }
    // after 1st round, read_haps/hap_to_alle_profile/hap_to_cons_alle are upToDate and will be used in the following rounds
    int max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "iter: %d\n", i_iter);
        // XXX var-wise loop update PS
        // update var->hap_to_cons_alle & var->phase_set
        int changed_hap1 = iter_update_var_hap_cons_phase_set(chunk, valid_var_idx, p, cand_vars, n_valid_vars, var_i_to_cate);
        // update
        int changed_hap2 = iter_update_var_hap_to_cons_alle(opt, chunk, valid_var_idx, n_valid_vars, p, cand_vars, var_i_to_cate, target_var_cate);
        if (changed_hap1 == 0 && changed_hap2 == 0) break;
        else {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "changed_hap1: %d changed_hap2: %d\n", changed_hap1, changed_hap2);
        }
    }
    // update PS for read & var after all iterations
    update_read_phase_set(chunk, var_is_valid, p, cand_vars);
    free(valid_var_idx); free(var_is_valid);
    return 0;
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

void sort_phase_sets(hts_pos_t *uniq_phase_set, int *phase_set_to_total_depth, int n_uniq_phase_set, int ***phase_set_to_hap_alle_profile) {
    for (int i = 0; i < n_uniq_phase_set-1; ++i) {
        for (int j = i+1; j < n_uniq_phase_set; ++j) {
            if (phase_set_to_total_depth[i] < phase_set_to_total_depth[j]) {
                int tmp = phase_set_to_total_depth[i]; phase_set_to_total_depth[i] = phase_set_to_total_depth[j]; phase_set_to_total_depth[j] = tmp;
                hts_pos_t tmp_ps = uniq_phase_set[i]; uniq_phase_set[i] = uniq_phase_set[j]; uniq_phase_set[j] = tmp_ps;
                int **tmp_profile = phase_set_to_hap_alle_profile[i]; phase_set_to_hap_alle_profile[i] = phase_set_to_hap_alle_profile[j]; phase_set_to_hap_alle_profile[j] = tmp_profile;
            }
        }
    }
}

// return:
// >= 0: phase set index
// -1: no valid phase set
// -2: not phased
static int select_somatic_phase_set0(hts_pos_t *uniq_phase_set, int n_uniq_phase_set, int ***phase_set_to_hap_alle_profile, int min_somatic_hap_depth) {
    int ps_i = -1;
    int *phase_set_to_total_depth = (int*)calloc(n_uniq_phase_set, sizeof(int));
    for (int i = 0; i < n_uniq_phase_set; ++i) {
        for (int j = 1; j < 3; ++j) {
            phase_set_to_total_depth[i] += (phase_set_to_hap_alle_profile[i][j][0] + phase_set_to_hap_alle_profile[i][j][1]);
        }
    }
    sort_phase_sets(uniq_phase_set, phase_set_to_total_depth, n_uniq_phase_set, phase_set_to_hap_alle_profile);
    free(phase_set_to_total_depth);
    int has_valid_ps = 0;
    for (int i = 0; i < n_uniq_phase_set; ++i) {
        int skip = 0, n_hap_with_alt = 0, n_valid_hap_with_alt = 0; //, n_total_alt = 0, n_total_ref = 0;
        for (int hap = 1; hap <= 2; ++hap) {
            // check if this hap has alt reads
            int n_alt = phase_set_to_hap_alle_profile[i][hap][1], n_ref = phase_set_to_hap_alle_profile[i][hap][0];
            // n_total_alt += n_alt; n_total_ref += n_ref;
            if (n_alt > 0) {
                n_hap_with_alt++;
                if (n_alt <= n_ref && n_alt+n_ref >= min_somatic_hap_depth)
                    n_valid_hap_with_alt++;
            }
        }
        if (n_hap_with_alt == 1 && n_valid_hap_with_alt == 1) { // XXX alt only show up in one haplotype
            has_valid_ps = 1;
            ps_i = i; break;
        }
    }
    if (has_valid_ps == 0) return -2;
    else return ps_i;
}

cand_somatic_var_aux_info_t *init_somatic_aux_info(int max_alt_dp) {
    cand_somatic_var_aux_info_t *aux_info = (cand_somatic_var_aux_info_t*)malloc(sizeof(cand_somatic_var_aux_info_t));
    memset(aux_info, 0, sizeof(cand_somatic_var_aux_info_t));
    aux_info->alt_read_ids = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->alt_quals = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->win_low_qual = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->dis_to_indel_error = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->no_dense_error = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->is_not_homopolymer_error = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->low_comp_reg_has_no_error = (int*)malloc(max_alt_dp * sizeof(int));
    return aux_info;
}

void free_somatic_var_aux_info(cand_somatic_var_aux_info_t *aux_info) {
    if (aux_info == NULL) return;
    free(aux_info->alt_read_ids);
    free(aux_info->alt_quals);
    free(aux_info->win_low_qual);
    free(aux_info->dis_to_indel_error);
    free(aux_info->no_dense_error);
    free(aux_info->is_not_homopolymer_error);
    free(aux_info->low_comp_reg_has_no_error);
    free(aux_info);
}

int get_min_dis_to_var(bam_chunk_t *chunk, int var_i) {
    int min_dis = INT32_MAX;
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    for (int i = var_i-1; i >= 0; --i) {
        if ((var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE) == 0) continue;
        min_dis = MIN(min_dis, cand_vars[var_i].pos - cand_vars[i].pos - cand_vars[i].ref_len);
        break;
    }
    for (int i = var_i+1; i < chunk->n_cand_vars; ++i) {
        if ((var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE) == 0) continue;
        min_dis = MIN(min_dis, cand_vars[i].pos - cand_vars[var_i].pos - cand_vars[var_i].ref_len);
        break;
    }
    return min_dis;
}

int get_read_win_min_qual(const call_var_opt_t *opt, digar_t *digar, int alt_qi) {
    uint8_t *qual = digar->qual; int qlen = digar->qlen;
    int flank_win_size = 3; // XXX opt->read_win_size; 7? 3+1+3
    int start = MAX(0, alt_qi - flank_win_size);
    int end = MIN(alt_qi + flank_win_size, qlen - 1);
    int min_qual = INT32_MAX; // default min qual
    for (int i = start; i <= end; ++i) {
        if (qual[i] < min_qual) min_qual = qual[i];
    }
    return min_qual;
}

int get_read_win_low_qual(const call_var_opt_t *opt, digar_t *digar, int alt_qi) {
    uint8_t *qual = digar->qual; int qlen = digar->qlen;
    int flank_win_size = 3; // XXX opt->read_win_size; 7? 3+1+3
    int start = MAX(0, alt_qi - flank_win_size);
    int end = MIN(alt_qi + flank_win_size, qlen - 1);
    int low_qual = INT32_MAX; // default low qual
    for (int i = start; i <= end; ++i) {
        if (qual[i] < low_qual) low_qual = qual[i];
    }
    return low_qual;
}

int digar_is_germline_var(bam_chunk_t *chunk, int var_i, digar1_t *digar1) {
    // check if the digar1 is a candidate variant
    if (digar1->type == BAM_CEQUAL) return 0; // skip equal bases
    for (int i = var_i; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        if (var->var_type == digar1->type && digar1->pos == var->pos && 
            ((var->var_type == BAM_CINS && var->alt_len == digar1->len) || (var->var_type == BAM_CDEL && var->ref_len == digar1->len)))
            return 1;
        else if (var->pos > digar1->pos) break; // no need to check further
    }
    for (int i = var_i-1; i >= 0; --i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        if (var->var_type == digar1->type && digar1->pos == var->pos &&
            ((var->var_type == BAM_CINS && var->alt_len == digar1->len) || (var->var_type == BAM_CDEL && var->ref_len == digar1->len)))
            return 1;
        else if (var->pos < digar1->pos) break; // no need to check further
    }
    return 0;
}

int digar_is_var(bam_chunk_t *chunk, int var_i, digar1_t *digar1) {
    // check if the digar1 is a candidate variant
    if (digar1->type == BAM_CEQUAL) return 0; // skip equal bases
    for (int i = var_i; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE) &&
            !(chunk->var_i_to_cate[i] & LONGCALLD_CAND_SOMATIC_VAR)) continue;
        if (var->var_type == digar1->type && digar1->pos == var->pos && 
            (var->var_type == BAM_CDIFF || (var->var_type == BAM_CINS && var->alt_len == digar1->len) || (var->var_type == BAM_CDEL && var->ref_len == digar1->len)))
            return 1;
        else if (var->pos > digar1->pos) break; // no need to check further
    }
    for (int i = var_i-1; i >= 0; --i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE) &&
            !(chunk->var_i_to_cate[i] & LONGCALLD_CAND_SOMATIC_VAR)) continue;
        if (var->var_type == digar1->type && digar1->pos == var->pos &&
            (var->var_type == BAM_CDIFF || (var->var_type == BAM_CINS && var->alt_len == digar1->len) || (var->var_type == BAM_CDEL && var->ref_len == digar1->len)))
            return 1;
        else if (var->pos < digar1->pos) break; // no need to check further
    }
    return 0;
}

// check if the read has indel error which is not a candidate variant
// XXX add digar_i to read_var_allele_profile???
int get_dis_to_seq_error(bam_chunk_t *chunk, int var_i, digar_t *digar, int var_type, int var_len, int alt_qi, int only_indel) {
    int dis_to_seq_error = 10; // maximum distance to indel error
    for (int i = 0; i < digar->n_digar; ++i) {
        digar1_t *digar1 = digar->digars + i;
        if (only_indel && digar1->type != BAM_CINS && digar1->type != BAM_CDEL) continue; // only check indels
        if (!only_indel && digar1->type != BAM_CINS && digar1->type != BAM_CDEL && digar1->type != BAM_CDIFF) continue; // skip equal bases
        if (alt_qi - digar1->qi > dis_to_seq_error) continue;
        // if (alt_qi == digar1->qi && digar1->type == var_type && digar1->len == var_len) continue; // skip itself
        if (digar1->qi - alt_qi > dis_to_seq_error) break; // no need to check further
        if (digar_is_var(chunk, var_i, digar1)) continue;
        dis_to_seq_error = MIN(dis_to_seq_error, ABS(digar1->qi - alt_qi));
        if (digar1->type == BAM_CINS) dis_to_seq_error = MIN(dis_to_seq_error, ABS(digar1->qi+digar1->len - alt_qi));
    }
    return dis_to_seq_error;
}

int digar1_is_homopolymer_error(bam_chunk_t *chunk, digar_t *digar, digar1_t *digar1, int *qi_hp_start, int *qi_hp_end) {
    if (digar1->type == BAM_CINS) {
        // INS is HP base
        uint8_t ins_base0 = digar1->alt_seq[0];
        for (int i = 1; i < digar1->len; ++i) {
            if (digar1->alt_seq[i] != ins_base0) return 0; // not a homopolymer
        }
        // INS is in homopolymer region
        *qi_hp_start = digar1->qi; *qi_hp_end = digar1->qi + digar1->len - 1;
        int i = digar1->qi-1;
        while (i >= 0) {
            if (seq_nt16_int[bam_seqi(digar->bseq, i)] != ins_base0) break; // not a homopolymer
            *qi_hp_start = i;
            i--;
        }
        i = digar1->qi+ digar1->len;
        while (i < digar->qlen) {
            if (seq_nt16_int[bam_seqi(digar->bseq, i)] != ins_base0) break; // not a homopolymer
            *qi_hp_end = i;
            i++;
        }
        if (*qi_hp_end - *qi_hp_start + 1 - digar1->len >= 3) return 1;
        else return 0;
    } else if (digar1->type == BAM_CDEL) {
        // DEL is HP
        char del_char0 = chunk->ref_seq[digar1->pos-chunk->ref_beg];
        uint8_t del_base0 = nst_nt4_table[(int)del_char0];
        for (int i = 1; i < digar1->len; ++i) {
            if (chunk->ref_seq[digar1->pos-chunk->ref_beg+i] != del_char0) return 0; // not a homopolymer
        }
        // DEL is in homopolymer region
        *qi_hp_start = digar1->qi; *qi_hp_end = digar1->qi;
        int i = digar1->qi-1;
        while (i >= 0) {
            if (seq_nt16_int[bam_seqi(digar->bseq, i)] != del_base0) break; // not a homopolymer
            *qi_hp_start = i;
            i--;
        }
        i = digar1->qi;
        while (i < digar->qlen) {
            if (seq_nt16_int[bam_seqi(digar->bseq, i)] != del_base0) break; // not a homopolymer
            *qi_hp_end = i;
            i++;
        }
        if (*qi_hp_end - *qi_hp_start + 1 >= 3) return 1;
        else return 0;
    } else { // MISMATCH

    }
    return 0;
}

// 1. get dis to homopolymer
// 2. check if the read has indel error, which is not a candidate variant
// 3. if yes, return dis, else return 10
// int get_dis_to_homopolymer_error(bam_chunk_t *chunk, int var_i, digar_t *digar, int alt_qi) {
//     int dis_to_homopolymer = 5; // maximum distance to homopolymer
//     for (int i = 0; i < digar->n_digar; ++i) {
//         digar1_t *digar1 = digar->digars + i;
//         // if (digar1->type != BAM_CINS && digar1->type != BAM_CDEL) continue; // only check indels
//         int hp_start = -1, hp_end = -1;
//         // print digar
//         if (!digar1_is_homopolymer_error(chunk, digar, digar1, &hp_start, &hp_end)) continue; // skip non-homopolymer errors
//         // fprintf(stderr, "alt_qi: %d %c i: %d hp_start: %d hp_end: %d %c\n", alt_qi, "ACGTN"[seq_nt16_int[bam_seqi(digar->bseq, alt_qi)]], i, hp_start, hp_end, "ACGTN"[seq_nt16_int[bam_seqi(digar->bseq, hp_start)]]);
//         // fprintf(stderr, "digar1: %ld %d %c %d\n", digar1->pos, digar1->qi, BAM_CIGAR_STR[digar1->type], digar1->len);
//         if (alt_qi - hp_end > dis_to_homopolymer) continue;
//         if (hp_start - alt_qi > dis_to_homopolymer) break; // no need to check further
//         if (digar_is_var(chunk, var_i, digar1)) continue;
//         dis_to_homopolymer = MIN(ABS(alt_qi - hp_end), ABS(alt_qi - hp_start));
//     }
//     return dis_to_homopolymer;
// }


int collect_noisy_read_info1(const call_var_opt_t *opt, bam_chunk_t *chunk, int read_id, hts_pos_t reg_beg, hts_pos_t reg_end, int *read_len,
                             uint8_t **read_seqs, int *fully_covers, int *read_reg_beg, int *read_reg_end) {
    digar_t *read_digars = chunk->digars+read_id; int n_digar = read_digars->n_digar; digar1_t *digars = read_digars->digars;
    hts_pos_t reg_digar_beg = -1, reg_digar_end = -1;
    int reg_read_beg = 0, reg_read_end = digar2qlen(read_digars)-1;
    if (read_digars->digars[0].type == BAM_CHARD_CLIP) reg_read_beg = read_digars->digars[0].len;
    if (read_digars->digars[n_digar-1].type == BAM_CHARD_CLIP) reg_read_end = read_digars->digars[n_digar-1].qi - 1;
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
    (*read_seqs) = (uint8_t*)malloc((reg_read_end - reg_read_beg + 1) * sizeof(uint8_t));
    for (int j = reg_read_beg; j <= reg_read_end; ++j) {
        (*read_seqs)[j-reg_read_beg] = seq_nt16_int[bam_seqi(read_digars->bseq, j)];
    }
    (*read_len) = reg_read_end - reg_read_beg + 1;
    (*fully_covers) = cover;
    (*read_reg_beg) = reg_read_beg; (*read_reg_end) = reg_read_end;
    return 0;
}

// >= 3 same base
int is_1mer_hp(uint8_t *seq, int i, int len) {
    if (i < 0 || i >= len) return 0; // out of bounds
    uint8_t base0 = seq[i];
    int hp_len = 1;
    for (int j = i+1; j < len; ++j) {
        if (seq[j] != base0) break; // not a homopolymer
        hp_len++;
    }
    if (hp_len >= 3) return hp_len; // >= 2 same base
    else return 0; // not a homopolymer
}

// >= 3 2-mers, e.g., ACACAC or GTGTGT
int is_2mer_hp(uint8_t *seq, int i, int len) {
    if (i < 0 || i >= len-1) return 0; // out of bounds
    uint8_t base0 = seq[i], base1 = seq[i+1];
    if (base0 == base1) return 0; // not a 2-mer homopolymer
    int hp_len = 2;
    for (int j = i+2; j < len; j += 2) {
        if (j >= len) break; // out of bounds
        if (seq[j] != base0) break;
        if (j+1 >= len) break;
        if (seq[j+1] != base1) break; // not a 2-mer homopolymer
        hp_len += 2;
    }
    if (hp_len >= 6) return hp_len; // >= 2 2-mers
    else return 0;
}

// compressed 1-mer and 2-mer homopolymer sequence
int get_hp_compressed_seq(uint8_t *seq, int len, uint8_t *hp_seq, int *hp_seq_len) {
    if (len == 0) {
        hp_seq[0] = '\0'; return 0;
    }
    int hp_len = 0, hp_len0, i = 0;
    while (i < len) {
        if ((hp_len0 = is_2mer_hp(seq, i, len)) > 0) {
            hp_seq[hp_len++] = seq[i]; hp_seq[hp_len++] = seq[i+1];
            hp_seq_len[hp_len-2] = -1; hp_seq_len[hp_len-1] = hp_len0; // -1 for 2-mer
            i+= hp_len0; 
        } else {
            if ((hp_len0 = is_1mer_hp(seq, i, len)) > 0) {
                hp_seq[hp_len++] = seq[i]; // not a homopolymer, add the base
                hp_seq_len[hp_len-1] = hp_len0; // 1-mer
                i += hp_len0;
            } else { // move by 1 base
                hp_seq[hp_len++] = seq[i]; // not a homopolymer, add the base
                hp_seq_len[hp_len-1] = 1; // 1-mer
                i++;
            }
        }
    }
    return hp_len;
}

int is_hp_compressed_match(uint8_t *seq1, int len1, uint8_t *seq2, int len2) {
    uint8_t *hp_seq1 = (uint8_t*)malloc(len1 + 1); int *hp_seq_len1 = (int*)malloc((len1 + 1) * sizeof(int));
    uint8_t *hp_seq2 = (uint8_t*)malloc(len2 + 1); int *hp_seq_len2 = (int*)malloc((len2 + 1) * sizeof(int));
    int hp_len1 = get_hp_compressed_seq(seq1, len1, hp_seq1, hp_seq_len1);
    int hp_len2 = get_hp_compressed_seq(seq2, len2, hp_seq2, hp_seq_len2);
    // fprintf(stderr, "hp_len1: %d hp_len2: %d\n", hp_len1, hp_len2);
    // fprintf(stderr, "hp_seq1: ");
    // for (int i = 0; i < hp_len1; ++i) {
    //     fprintf(stderr, "%d%c", hp_seq_len1[i], "ACGTN-"[hp_seq1[i]]);
    // } fprintf(stderr, "\n");
    // fprintf(stderr, "hp_seq2: ");
    // for (int i = 0; i < hp_len2; ++i) {
    //     fprintf(stderr, "%d%c", hp_seq_len2[i], "ACGTN-"[hp_seq2[i]]);
    // } fprintf(stderr, "\n");

    int ret = 1;
    if (hp_len1 != hp_len2 || hp_len1 == 0) ret = 0;
    else {
        for (int i = 0; i < hp_len1; ++i) {
            if (hp_seq1[i] != hp_seq2[i] || (hp_seq_len1[i] > 0 && hp_seq_len2[i] < 0) || (hp_seq_len1[i] < 0 && hp_seq_len2[i] > 0)) {
                ret = 0; break;
            }
        }
    }
    free(hp_seq1); free(hp_seq2); free(hp_seq_len1); free(hp_seq_len2);
    return ret;
}

int var_low_comp_reg_has_error(bam_chunk_t *chunk, int var_i, digar_t *digar, hts_pos_t low_comp_beg, hts_pos_t low_comp_end) {
    // fprintf(stderr, "var %ld low_comp_beg: %ld low_comp_end: %ld\n", chunk->cand_vars[var_i].pos, low_comp_beg, low_comp_end);
    for (int i = 0; i < digar->n_digar; ++i) {
        digar1_t *digar1 = digar->digars + i;
        if (digar1->type != BAM_CINS && digar1->type != BAM_CDEL && digar1->type != BAM_CDIFF) continue; // skip equal bases
        hts_pos_t digar1_end = digar1->pos;
        if (digar1->type == BAM_CDEL) digar1_end += digar1->len - 1;
        if (digar1_end < low_comp_beg) continue;
        if (digar1->pos > low_comp_end) break; // no need to check further
        if (digar_is_var(chunk, var_i, digar1)) continue;
        // fprintf(stderr, "digar1: %ld %d %c %d\n", digar1->pos, digar1->qi, BAM_CIGAR_STR[digar1->type], digar1->len);
        return 1;
    }
    return 0;
}

int is_diff_between_ref_hap_aln(const call_var_opt_t *opt, uint8_t *read_reg_seq, int read_reg_len, uint8_t *hap_reg_seq, int hap_reg_len, uint8_t *ref_reg_seq, int ref_reg_len, int alt_ref_pos) {
    int is_diff = 0;
    // wfa alignment between read_reg_seq and hap_reg_seq
    uint8_t *read_hap_aln_str, *read_ref_aln_str, *hap_aln_str, *ref_aln_str;
    int hap_aln_len = 0, ref_aln_len = 0;
    wfa_end2end_aln(hap_reg_seq, hap_reg_len, read_reg_seq, read_reg_len, opt->gap_aln, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1,
                    opt->gap_open2, opt->gap_ext2, NULL, NULL, &hap_aln_str, &read_hap_aln_str, &hap_aln_len);
    // wfa alignment between read_reg_seq and ref_reg_seq
    wfa_end2end_aln(ref_reg_seq, ref_reg_len, read_reg_seq, read_reg_len, opt->gap_aln, opt->match, opt->mismatch, opt->gap_open1, opt->gap_ext1,
                    opt->gap_open2, opt->gap_ext2, NULL, NULL, &ref_aln_str, &read_ref_aln_str, &ref_aln_len);
    // check if the alignment at position alt_qi is the same or not
    int hap_aln_i = -1, ref_aln_i = -1;
    int read_aln_i = -1, alt_read_pos = -1;
    for (int i = 0; i < ref_aln_len; ++i) {
        if (ref_aln_str[i] != 5) ref_aln_i++;
        if (ref_aln_i == alt_ref_pos) { alt_read_pos = i; break; }
    }
    read_aln_i = -1;
    for (int i = 0; i < hap_aln_len; ++i) {
        if (read_hap_aln_str[i] != 5) read_aln_i++;
        if (read_aln_i == alt_read_pos) { hap_aln_i = i; break; }
    }
    if (hap_aln_i < 0 || ref_aln_i < 0) {
        // fprintf(stderr, "alt_qi: %d hap_aln_i: %d ref_aln_i: %d\n", alt_qi, hap_aln_i, ref_aln_i);
        is_diff = 1; // cannot find the position in the alignment, consider as different
    } else {
        // fprintf(stderr, "alt_qi: %d hap_aln_i: %d ref_aln_i: %d\n", alt_qi, hap_aln_i, ref_aln_i);
        if (hap_aln_str[hap_aln_i] != ref_aln_str[ref_aln_i]) is_diff = 1;
    }
    free(hap_aln_str); free(ref_aln_str);
    return is_diff;
}

// check if 
// 1) HP-compressed seq of read is the same as the haplotype seq (ref+HapGermlineVars)
// or 
// 2) alignment between read-haplotype seq and read-reference seq is the same at the candidate variant position
//   return 1 if match (seq error), 0 if not match (potential var)
int var_is_homopolymer_error(const call_var_opt_t *opt, bam_chunk_t *chunk, int hap, int read_i, hts_pos_t low_comp_beg, hts_pos_t low_comp_end, hts_pos_t var_pos) {
    int read_reg_len = 0, read_reg_beg, read_reg_end, full_cover; uint8_t *read_reg_seq = NULL;
    hts_pos_t ref_reg_beg = low_comp_beg - opt->noisy_reg_flank_len, ref_reg_end = low_comp_end + opt->noisy_reg_flank_len;
    int alt_ref_pos = var_pos - ref_reg_beg; // alt_qi in read_reg_seq
    collect_noisy_read_info1(opt, chunk, read_i, ref_reg_beg, ref_reg_end, &read_reg_len, &read_reg_seq, &full_cover, &read_reg_beg, &read_reg_end);
    // add var->alt_seq to ref_reg_seq
    int all_ins_len = 0;
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
        int var_cate = chunk->var_i_to_cate[i];
        if (!(var_cate & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        hts_pos_t var_pos = var->pos; hts_pos_t var_end = var_pos + var->ref_len - 1;
        if (var_pos < ref_reg_beg) continue;
        if (var_end > ref_reg_end) break;
        int allele_i = var->hap_to_cons_alle[hap];
        if (allele_i != 1) continue;
        if (var->var_type == BAM_CINS) {
            all_ins_len += var->alt_len;
        } else if (var->var_type == BAM_CDEL) {
            all_ins_len -= var->ref_len;
        }
    }
    uint8_t *hap_reg_seq = (uint8_t*)malloc(ref_reg_end - ref_reg_beg + all_ins_len + 1);
    uint8_t *ref_reg_seq = (uint8_t*)malloc(ref_reg_end - ref_reg_beg + 1);
    hts_pos_t last_pos = ref_reg_beg; int hap_reg_len = 0;
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
        int var_cate = chunk->var_i_to_cate[i];
        if (!(var_cate & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        hts_pos_t var_pos = var->pos; hts_pos_t var_end = var_pos + var->ref_len - 1;
        if (var_pos < ref_reg_beg) continue;
        if (var_end > ref_reg_end) break;
        int allele_i = var->hap_to_cons_alle[hap];
        if (allele_i != 1) continue;
        // last_pos ... var_pos-1
        for (hts_pos_t j = last_pos; j < var_pos; ++j) hap_reg_seq[hap_reg_len++] = nst_nt4_table[(int)chunk->ref_seq[j - chunk->ref_beg]];
        // alt_seq
        for (int j = 0; j < var->alt_len; ++j) hap_reg_seq[hap_reg_len++] = var->alt_seq[j];
        // update last_pos as var_end + 1
        last_pos = var_end + 1;
    }
    // last_pos ... ref_reg_end
    for (hts_pos_t j = last_pos; j <= ref_reg_end; ++j) {
        hap_reg_seq[hap_reg_len++] = nst_nt4_table[(int)chunk->ref_seq[j - chunk->ref_beg]];
    }
    int ref_reg_len = 0;
    for (hts_pos_t j = ref_reg_beg; j <= ref_reg_end; ++j) {
        ref_reg_seq[ref_reg_len++] = nst_nt4_table[(int)chunk->ref_seq[j - chunk->ref_beg]];
    }
    int is_hp_compresed = is_hp_compressed_match(read_reg_seq, read_reg_len, hap_reg_seq, hap_reg_len);
    int is_diff_between_ref_cons_aln = 0;
    if (is_hp_compresed == 0)
       is_diff_between_ref_cons_aln = is_diff_between_ref_hap_aln(opt, read_reg_seq, read_reg_len, hap_reg_seq, hap_reg_len, ref_reg_seq, ref_reg_len, alt_ref_pos);
    // fprintf(stderr, "ref_seq: ");
    // for (int i = 0; i < ref_reg_len; ++i) {
    //     fprintf(stderr, "%c", "ACGTN-"[ref_reg_seq[i]]);
    // } fprintf(stderr, "\n");
    // fprintf(stderr, "read_seq: ");
    // for (int i = 0; i < read_reg_len; ++i) {
    //     fprintf(stderr, "%c", "ACGTN-"[read_reg_seq[i]]);
    // } fprintf(stderr, "\n");
    // fprintf(stderr, "HP-match-ret: %d\n", ret);

    free(read_reg_seq); free(ref_reg_seq); free(hap_reg_seq);
    if (is_hp_compresed || is_diff_between_ref_cons_aln) return 1; // seq error
    else return 0; // potential var
}

int get_read_win_median_qual(digar_t *d, int alt_qi, int len) {
    int *quals = (int*)malloc(len * sizeof(int));
    for (int i = 0; i < len; ++i) {
        int pos = alt_qi + i;
        if (pos < 0 || pos >= d->qlen) {
            quals[i] = 0; // default min qual
        } else {
            quals[i] = d->qual[pos];
        }
    }
    int median_qual = median_int(quals, len);
    free(quals);
    return median_qual;
}

int get_alt_qual(digar_t *d, int var_type, int var_len, int alt_qi) {
    if (var_type == BAM_CDIFF) return d->qual[alt_qi];
    else if (var_type == BAM_CINS) return get_read_win_median_qual(d, alt_qi, var_len);
    else if (var_type == BAM_CDEL) return get_read_win_median_qual(d, alt_qi-1, 2);
    else return d->qual[alt_qi];
}

// return:
//   1: in low-complexity region
//   2: near low-complexity region
//   0: none
int var_is_low_comp_reg(bam_chunk_t *chunk, cand_var_t *var, hts_pos_t *low_comp_beg, hts_pos_t *low_comp_end) {
    // check if the variant is in low-complexity regions
    int64_t low_comp_n = 0, max_low_comp_n = 0;
    if (chunk->low_comp_cr == NULL || chunk->low_comp_cr->n_r == 0) return 0;
    int64_t *low_comp_b = 0; 
    low_comp_n = cr_overlap(chunk->low_comp_cr, "cr", var->pos, var->pos+var->ref_len-1, &low_comp_b, &max_low_comp_n);
    if (low_comp_n >= 1) { // var is in low-complexity region
        *low_comp_beg = cr_start(chunk->low_comp_cr, low_comp_b[0])+1;
        *low_comp_end = cr_end(chunk->low_comp_cr, low_comp_b[low_comp_n-1]);
        free(low_comp_b);
        return 1;
    } else { // check if var is near low-complexity region
        int flank_len = 5;
        low_comp_n = cr_overlap(chunk->low_comp_cr, "cr", var->pos-flank_len, var->pos+var->ref_len+flank_len-1, &low_comp_b, &max_low_comp_n);
        if (low_comp_n >= 1) {
            *low_comp_beg = cr_start(chunk->low_comp_cr, low_comp_b[0])+1-flank_len;
            *low_comp_end = cr_end(chunk->low_comp_cr, low_comp_b[low_comp_n-1])+flank_len;
            free(low_comp_b);
            return 2;
        }
    }
    free(low_comp_b);
    return 0;
}

// check if the variant is in dense-error region, > max_err errors in a window of win bases
int has_dense_error(bam_chunk_t *chunk, cand_var_t *var, int var_i, digar_t *digar, int alt_qi) {
    int win = 50, max_err = 3;
    int left_has_dense_error = 0, right_has_dense_error = 0;
    int n_err = 0, i = 0;
    for (i = 0; i < digar->n_digar; ++i) {
        digar1_t *digar1 = digar->digars + i;
        if (digar1->type != BAM_CINS && digar1->type != BAM_CDEL && digar1->type != BAM_CDIFF) continue;
        if (alt_qi - digar1->qi > win) continue;
        if (digar1->qi == alt_qi && digar1->type == var->var_type) break;
        if (digar_is_var(chunk, var_i, digar1)) continue;
        n_err++;
        if (n_err >= max_err) {
            left_has_dense_error = 1;
            break;
        }
    }
    n_err = 0;
    for (; i < digar->n_digar; ++i) {
        digar1_t *digar1 = digar->digars + i;
        if (digar1->type != BAM_CINS && digar1->type != BAM_CDEL && digar1->type != BAM_CDIFF) continue;
        if (digar1->qi - alt_qi > win) break;
        if (digar1->qi == alt_qi && digar1->type == var->var_type) continue;
        if (digar_is_var(chunk, var_i, digar1)) continue;
        n_err++;
        if (n_err >= max_err) {
            right_has_dense_error = 1;
            break;
        }
    }
    if (left_has_dense_error == 1 || right_has_dense_error == 1) return 1;
    return 0;
}

// consider both SNV and INDEL
cand_somatic_var_aux_info_t *collect_somatic_var_aux_info(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t ps, int alt_hap, int var_i, int64_t ovlp_n, int64_t *ovlp_b) {
    cand_somatic_var_aux_info_t *aux_info = init_somatic_aux_info(ovlp_n);
    cand_var_t *var = chunk->cand_vars + var_i; cgranges_t *read_var_cr = chunk->read_var_cr; read_var_profile_t *p = chunk->read_var_profile;
    int64_t ovlp_i; hts_pos_t low_comp_beg, low_comp_end;

    aux_info->is_low_comp = var_is_low_comp_reg(chunk, var, &low_comp_beg, &low_comp_end);
    aux_info->min_dis_to_var = get_min_dis_to_var(chunk, var_i);
    int alt_i = 0;
    for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
        int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
        if (chunk->is_skipped[read_i]) continue;
        aux_info->total_dp += 1;
        int hap = chunk->haps[read_i]; hts_pos_t phase_set = chunk->phase_sets[read_i];
        if (ps != -1 && alt_hap != 0) {
            if (hap != alt_hap || phase_set != ps) continue; // only collect reads with the same haplotype and phase set
        }
        read_var_profile_t *p1 = p + read_i; int var_idx = var_i - p1->start_var_idx; int alle_i = p1->alleles[var_idx];
        int alt_qi = p1->alt_qi[var_idx];
        digar_t *d = chunk->digars + read_i;
        // if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read %s %" PRIi64 " hap: %d phase_set: %" PRId64 " alle_i: %d\n", bam_get_qname(chunk->reads[read_i]), var->pos, hap, phase_set, alle_i);
        aux_info->hap_total_dp += 1;
        if (alle_i == 1 && alt_qi != -1) {
            if (d->is_rev) aux_info->hap_alt_rev_cov += 1; else aux_info->hap_alt_for_cov += 1;
            aux_info->alt_read_ids[alt_i] = read_i;
            aux_info->hap_alt_dp += 1;
            aux_info->alt_quals[alt_i] = get_alt_qual(d, var->var_type, var->alt_len, alt_qi);
            aux_info->win_low_qual[alt_i] = get_read_win_low_qual(opt, d, alt_qi);
            aux_info->dis_to_indel_error[alt_i] = get_dis_to_seq_error(chunk, var_i, d, var->var_type, var->alt_len, alt_qi, 1);
            aux_info->no_dense_error[alt_i] = 1-has_dense_error(chunk, var, var_i, d, alt_qi);
            // fprintf(stderr, "%" PRIi64 "\n", var->pos);
            if (var->var_type == BAM_CDIFF) {
                if (aux_info->is_low_comp > 0) aux_info->is_not_homopolymer_error[alt_i] = 1-var_is_homopolymer_error(opt, chunk, hap, read_i, low_comp_beg, low_comp_end, var->pos);
                else aux_info->is_not_homopolymer_error[alt_i] = 1-var_is_homopolymer_error(opt, chunk, hap, read_i, var->pos, var->pos+var->ref_len-1, var->pos);
            } else aux_info->is_not_homopolymer_error[alt_i] = 1; // default value, no need to check
            // var is in low-complexity region, check if the read has error in low-complexity region
            if (aux_info->is_low_comp == 1) aux_info->low_comp_reg_has_no_error[alt_i] = 1-var_low_comp_reg_has_error(chunk, var_i, d, low_comp_beg, low_comp_end);
            else aux_info->low_comp_reg_has_no_error[alt_i] = 1; // default value, no need to check
            alt_i++;
        } else { // ref
            if (d->is_rev) aux_info->hap_ref_rev_cov += 1; else aux_info->hap_ref_for_cov += 1;
        }
    }
    return aux_info;
}

// comp two INS/DEL variants that fuzzy_ovlp with each other
// return 0: fuzzy match; non-0: not fuzzy match
int vntr_fuzzy_comp_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *var1, cand_var_t *var2) {
    if (var1->var_type == BAM_CDEL && var2->var_type == BAM_CDEL) {
        int min_len = MIN_OF_TWO(var1->ref_len, var2->ref_len);
        int max_len = MAX_OF_TWO(var1->ref_len, var2->ref_len);
        if (min_len < max_len * 0.8) return 1;
        if ((vntr_fuzzy_comp_seq(opt, (uint8_t*)chunk->ref_seq+var1->pos-chunk->ref_beg, var1->ref_len, (uint8_t*)chunk->ref_seq+var2->pos-chunk->ref_beg, var2->ref_len) == 0))
            return 0;

    } else if (var1->var_type == BAM_CINS && var2->var_type == BAM_CINS) {
        int min_len = MIN_OF_TWO(var1->alt_len, var2->alt_len);
        int max_len = MAX_OF_TWO(var1->alt_len, var2->alt_len);
        if (min_len < max_len * 0.8) return 1;
        if ((vntr_fuzzy_comp_seq(opt, var1->alt_seq, var1->alt_len, var2->alt_seq, var2->alt_len) == 0)) return 0;
    }
    return 1;
}

// comp two INS variants that have different length (difference is large)
// return 0: difference is low-complexity
// return 1: difference is not low-complexity
int low_comp_ins_comp_var(const call_var_opt_t *opt, bam_chunk_t *chunk, cand_var_t *large_var, cand_var_t *small_var) {
    if (large_var->var_type == BAM_CINS && small_var->var_type == BAM_CINS) {
        int large_ins_len = large_var->alt_len, small_ins_len = small_var->alt_len;
        if (large_ins_len < small_ins_len) return 0; // not low_comp_ins
        // if ((vntr_fuzzy_comp_seq(opt, large_var->alt_seq, large_var->alt_len, small_var->alt_seq, small_var->alt_len) == 0)) return 0;
        int diff_len = 0; uint8_t *diff_seq = 0;
        diff_len = wfa_collect_diff_ins_seq(opt, large_var->alt_seq, large_var->alt_len, small_var->alt_seq, small_var->alt_len, &diff_seq);
        uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
        r = sdust(0, diff_seq, diff_len, T, W, &n);
        int low_comp_len = 0;
        for (int i = 0; i < n; ++i) {
            // fprintf(stderr, "low_comp_cr: %s %d-%d\n", chunk->tname, chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1);
            low_comp_len += ((int)r[i]- (int)(r[i]>>32));
        }
        free(r); if (diff_len > 0) free(diff_seq);
        if (low_comp_len > diff_len * 0.8) return 1;
    }
    return 0;
}

// check if var 
// 1. fuzzily match germline variant, if so: potentially artifact
// 2. different from germline variant but the difference is low-complexity, if so: potentially artifact
int var_is_germline(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    int *var_i_to_cate = chunk->var_i_to_cate;
    assert(var->var_type == BAM_CINS || var->var_type == BAM_CDEL);
    hts_pos_t var_beg, var_end; int var_len;
    var_beg = var->pos;
    if (var->var_type == BAM_CDEL) {
        var_len = var->ref_len; var_end = var->pos + var->ref_len - 1;
    } else {
        var_len = var->alt_len; var_end = var->pos;
    }
    int var_win = MAX(500, var_len);
    // 1. check if fuzzy match
    for (int i = var_i+1; i < chunk->n_cand_vars; ++i) {
        if (!(var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        cand_var_t *var1 = chunk->cand_vars + i;
        if (var1->pos - var_end > var_win) break; // no need to check further

        if (vntr_fuzzy_comp_var(opt, chunk, var, var1) == 0) return 1;
    }
    for (int i = var_i-1; i >= 0; --i) {
        if (!(var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        cand_var_t *var1 = chunk->cand_vars + i;
        hts_pos_t var1_end = var1->pos;
        if (var->var_type == BAM_CDEL) var1_end += var1->ref_len - 1;
        if (var_beg - var1_end > var_win) break; // no need to check further
        if (vntr_fuzzy_comp_var(opt, chunk, var, var1) == 0) return 1;
    }
    // 2. check if the variant has artifact low-complexity insertion sequences
    var_win = 50; // only check 50 bases around the variant
    int not_germline = 0;
    for (int i = var_i+1; i < chunk->n_cand_vars; ++i) {
        if (!(var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        cand_var_t *var1 = chunk->cand_vars + i;
        if (var1->pos - var_end > var_win) break; // no need to check further
        if (low_comp_ins_comp_var(opt, chunk, var, var1) == 1) return 1;
    }
    for (int i = var_i-1; i >= 0; --i) {
        if (!(var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip non-germline vars
        cand_var_t *var1 = chunk->cand_vars + i;
        hts_pos_t var1_end = var1->pos;
        if (var->var_type == BAM_CDEL) var1_end += var1->ref_len - 1;
        if (var_beg - var1_end > var_win) break; // no need to check further
        if (low_comp_ins_comp_var(opt, chunk, var, var1) == 1) return 1;
    }
    return 0;
}

int sv_is_te(cand_var_t *var) {
    int tsd_len = var->tsd_len, polya_len = var->polya_len, te_seq_i = var->te_seq_i;
    int min_tsd_len = 5, min_polya_len = 20;
    int n_criteria = 0;
    if (tsd_len >= min_tsd_len) n_criteria++;
    if (abs(polya_len) >= min_polya_len) n_criteria++;
    if (te_seq_i >= 0) n_criteria++;
    if (n_criteria >= 2) return 1; // at least two criteria are met
    else return 0; // not a TE
}

int somatic_var_seq_is_low_comp(bam_chunk_t *chunk, int var_i) {
    cand_var_t *var = chunk->cand_vars + var_i;
    if (var->var_type == BAM_CDIFF) return 0; // SNVs are not low complexity
    uint8_t *var_seq; int var_len;
    if (var->var_type == BAM_CINS) {
        var_len = var->alt_len;
        var_seq = var->alt_seq;
    } else {
        var_len = var->ref_len;
        var_seq = (uint8_t*)chunk->ref_seq + var->pos - chunk->ref_beg;
    }
    if (var->tsd_len > 0 && abs(var->polya_len) > 0) {
        if (var->tsd_len + abs(var->polya_len) > var_len * 0.8) return 1;
    }
    uint64_t *r; int n=0, T=LONGCALLD_SDUST_T, W=LONGCALLD_SDUST_W;
    r = sdust(0, var_seq, var_len, T, W, &n);
    int low_comp_len = 0;
    for (int i = 0; i < n; ++i) {
        // fprintf(stderr, "low_comp_cr: %s %d-%d\n", chunk->tname, chunk->reg_beg+(int)(r[i]>>32)-1, chunk->reg_beg+(int)r[i]-1);
        low_comp_len += ((int)r[i]- (int)(r[i]>>32));
    }
    free(r);
    if (low_comp_len > var_len * 0.8) return 1;
    else return 0;
}

int phased_sv_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    // check TE infomation
    cand_var_t *var = chunk->cand_vars + var_i;
    if (var->alle_covs[1]  < opt->min_somatic_alt_dp) {
        if (var->alle_covs[1] < opt->min_somatic_te_dp || sv_is_te(chunk->cand_vars + var_i) == 0)
            return 0;
    }
    if (var_is_germline(opt, chunk, var_i, aux_info)) return 0;
    if (var->alle_covs[1] == 1 && somatic_var_seq_is_low_comp(chunk, var_i)) return 0;
    if (median_int(aux_info->no_dense_error, aux_info->hap_alt_dp) == 0) return 0;
    // XXX For SVs, low_comp_reg should be VNTR regions, not homopolymer regions
    // if (min_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp) == 0) return 0; // low complexity region has error, not somatic
    return 1;
}

int no_phase_sv_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    // check TE infomation
    cand_var_t *var = chunk->cand_vars + var_i;
    if (var->alle_covs[1]  < opt->min_somatic_alt_dp) {
        if (var->alle_covs[1] < opt->min_somatic_te_dp || sv_is_te(chunk->cand_vars + var_i) == 0)
            return 0;
    }
    if (var_is_germline(opt, chunk, var_i, aux_info)) return 0;
    if (var->alle_covs[1] == 1 && somatic_var_seq_is_low_comp(chunk, var_i)) return 0;
    if (median_int(aux_info->no_dense_error, aux_info->hap_alt_dp) == 0) return 0;
    // XXX For SVs, low_comp_reg should be VNTR regions, not homopolymer regions
    // if (min_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp) == 0) return 0; // low complexity region has error, not somatic
    return 1;
}

// check all filters
int phased_snv_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    int is_somatic = 1;
    fprintf(stderr, "%s\t%" PRIi64 "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", chunk->tname, var->pos,
            aux_info->total_dp, aux_info->hap_total_dp, aux_info->hap_alt_dp, 
            aux_info->hap_ref_for_cov, aux_info->hap_ref_rev_cov, aux_info->hap_alt_for_cov, aux_info->hap_alt_rev_cov,
            aux_info->hap_alt_dp * 1.0 / aux_info->hap_total_dp,
            aux_info->min_dis_to_var, median_int(aux_info->alt_quals, aux_info->hap_alt_dp), median_int(aux_info->win_low_qual, aux_info->hap_alt_dp),
            median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp), median_int(aux_info->no_dense_error, aux_info->hap_alt_dp),
            median_int(aux_info->is_not_homopolymer_error, aux_info->hap_alt_dp), median_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp));

    if (aux_info->hap_alt_dp < opt->min_somatic_alt_dp) return 0;
    if (aux_info->min_dis_to_var < opt->min_somatic_dis_to_var) return 0;
    // if (median_int(aux_info->alt_quals, aux_info->hap_alt_dp) < 27) return 0;
    if (median_int(aux_info->alt_quals, aux_info->hap_alt_dp) < chunk->median_qual) return 0;
    // if (median_int(aux_info->win_low_qual, aux_info->hap_alt_dp) < 17) return 0;
    if (median_int(aux_info->win_low_qual, aux_info->hap_alt_dp) < chunk->first_quar_qual) return 0;
    if (median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp) < opt->min_somatic_dis_to_seq_error) return 0;
    if (median_int(aux_info->no_dense_error, aux_info->hap_alt_dp) == 0) return 0;
    if (min_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp) == 0) return 0; // low complexity region has error, not somatic
    if (min_int(aux_info->is_not_homopolymer_error, aux_info->hap_alt_dp) == 0) return 0;
    return is_somatic;
}

int no_phase_snv_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    int is_somatic = 1;
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "%s\t%" PRIi64 "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n", chunk->tname, var->pos,
                aux_info->total_dp, aux_info->hap_total_dp, aux_info->hap_alt_dp, 
                aux_info->hap_ref_for_cov, aux_info->hap_ref_rev_cov, aux_info->hap_alt_for_cov, aux_info->hap_alt_rev_cov,
                aux_info->hap_alt_dp * 1.0 / aux_info->hap_total_dp,
                aux_info->min_dis_to_var, median_int(aux_info->alt_quals, aux_info->hap_alt_dp), median_int(aux_info->win_low_qual, aux_info->hap_alt_dp),
                median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp), 
                median_int(aux_info->no_dense_error, aux_info->hap_alt_dp), median_int(aux_info->is_not_homopolymer_error, aux_info->hap_alt_dp), median_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp));
    }

    if (aux_info->hap_alt_dp < opt->min_somatic_alt_dp) return 0;
    if (aux_info->min_dis_to_var < opt->min_somatic_dis_to_var) return 0;
    if (median_int(aux_info->alt_quals, aux_info->hap_alt_dp) < chunk->third_quar_qual) return 0;
    if (median_int(aux_info->win_low_qual, aux_info->hap_alt_dp) < chunk->median_qual) return 0;
    if (median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp) < opt->min_somatic_dis_to_seq_error) return 0;
    if (median_int(aux_info->no_dense_error, aux_info->hap_alt_dp) == 0) return 0;
    if (min_int(aux_info->low_comp_reg_has_no_error, aux_info->hap_alt_dp) == 0) return 0; // low complexity region has error, not somatic
    if (min_int(aux_info->is_not_homopolymer_error, aux_info->hap_alt_dp) == 0) return 0;
    return is_somatic;
}

// check all filters
int phased_var_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    // XXX no small indels for now
    if (var->var_type == BAM_CDIFF) return phased_snv_is_somatic(opt, chunk, var_i, aux_info);
    else return phased_sv_is_somatic(opt, chunk, var_i, aux_info);
}

int no_phase_var_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    // XXX no small indels for now
    if (var->var_type == BAM_CDIFF) return no_phase_snv_is_somatic(opt, chunk, var_i, aux_info);
    else return no_phase_sv_is_somatic(opt, chunk, var_i, aux_info);
}

void mark_somatic_var(cand_var_t *var, hts_pos_t phase_set, int alt_hap) {
    var->phase_set = phase_set; // set phase set
    var->hap_to_cons_alle[alt_hap] = 1; // alt allele
    var->hap_to_cons_alle[3-alt_hap] = 0; // ref allele
}

void mark_non_somatic_var(cand_var_t *var) {
    var->phase_set = -1; // no phase set
    var->hap_to_cons_alle[1] = var->hap_to_cons_alle[2] = 0;
}

void mark_skip_somatic_reads_based_on_invalid_var(bam_chunk_t *chunk, int var_i) {
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
    for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
        int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
        // fprintf(stderr, "mark read %d as skipped for somatic variant calling\n", read_i);
        chunk->is_skipped_for_somatic[read_i] = 1; // tag read as skipped for somatic variant calling
    }
}

void mark_invalid_somatic_vars(int somatic_win, int somatic_max_var, hts_pos_t *somatic_var_beg, hts_pos_t *somatic_var_end, int n_somatic_vars, int *is_invalid_somatic) {
    if (n_somatic_vars <= somatic_max_var) return; // no need to mark somatic variants
    // for any sliding window of size somatic_dense_win, if there are more than somatic_dense_max_var somatic variants, mark them as invalid
    for (int i = 0; i < n_somatic_vars - somatic_max_var; ++i) {
        int j = i + somatic_max_var;
        if (somatic_var_beg[j] - somatic_var_end[i] < somatic_win) {
            for (int k = i; k <= j; ++k) {
                is_invalid_somatic[k] = 1; // mark as invalid
            }
        }
    }
}

void mark_somatic_reads(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    int somatic_win = opt->somatic_win; // default 1000-bp window
    int somatic_max_var = opt->somatic_win_max_vars; // maximum number of somatic variants in a dense window
    hts_pos_t *somatic_var_beg = (hts_pos_t*)malloc(chunk->n_cand_vars * sizeof(hts_pos_t));
    hts_pos_t *somatic_var_end = (hts_pos_t*)malloc(chunk->n_cand_vars * sizeof(hts_pos_t));
    int *somatic_var_i = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    int *is_invalid_somatic_var = (int*)calloc(chunk->n_cand_vars, sizeof(int)); // mark somatic variants
    int n_somatic_vars = 0;
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue; // skip non-target variants
        cand_var_t *var = cand_vars + var_i;
        if (var->hap_to_cons_alle[1] == 0 && var->hap_to_cons_alle[2] == 0) continue; // skip variants without alt alleles
        somatic_var_i[n_somatic_vars] = var_i; // collect somatic variant indices
        somatic_var_beg[n_somatic_vars++] = var->pos; // collect somatic variant positions
        somatic_var_end[n_somatic_vars-1] = var->pos + var->ref_len-1; // collect somatic variant end positions
    }
    // var-wise: mark invalid var & reads
    mark_invalid_somatic_vars(somatic_win, somatic_max_var, somatic_var_beg, somatic_var_end, n_somatic_vars, is_invalid_somatic_var);
    for (int i = 0; i < n_somatic_vars; ++i) {
        int var_i = somatic_var_i[i];
        if (is_invalid_somatic_var[i] == 1) mark_skip_somatic_reads_based_on_invalid_var(chunk, var_i);
    }
    free(somatic_var_beg); free(somatic_var_end);
    free(somatic_var_i); free(is_invalid_somatic_var);
}

void mark_somatic_vars(bam_chunk_t *chunk, int target_var_cate) {
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue; // skip non-target variants
        cand_var_t *var = cand_vars + var_i;
        if (var->hap_to_cons_alle[1] == 0 && var->hap_to_cons_alle[2] == 0) continue; // skip variants without alt alleles
        cand_somatic_var_aux_info_t *aux_info = var->somatic_aux_info;
        for (int i = 0; i < aux_info->hap_alt_dp; ++i) {
            int read_i = aux_info->alt_read_ids[i];
            // XXX here any of the somatic var's covering reads is marked as skipped for somatic var calling, mark the somatic var as non-somatic var
            if (chunk->is_skipped_for_somatic[read_i]) { 
                mark_non_somatic_var(var); // mark as non-somatic variant
                break;
            }
        }
    }
}

// post-process somatic variants:
// 1) tag reads with dense somatic variants, i.e., >= 2 somatic variants in 1000-bp window
//    re-tag somatic variants based on tagged reads
// 2) collect consensus for overlapping somatic variants, i.e., INSs
int post_process_somatic_vars(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    mark_somatic_reads(opt, chunk, target_var_cate);
    mark_somatic_vars(chunk, target_var_cate);
    return 0;
}

// somatic SNP:
// 1) alt reads all come from one haplotype, no alt reads come from ref haplotype or non-haplotype
// 2) alt freqency in alt haplotype < 0.5
// 3) total depth in alt haplotype >= 10, total depth in ref haplotype >= 10
// somatic SVs, look for TE
// pick phase_set & alt_hap that have no conflict haplotype for somatic variant
hts_pos_t select_somatic_phase_set_alt_hap(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, read_var_profile_t *p, int *alt_hap) { 
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    hts_pos_t ps = -1;
    cand_var_t *var = chunk->cand_vars + var_i;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);

    int min_somatic_hap_depth = opt->min_somatic_hap_dp;
    hts_pos_t *uniq_phase_set = (hts_pos_t*)calloc(ovlp_n, sizeof(hts_pos_t)); int n_uniq_phase_set = 0;
    int ***phase_set_to_hap_alle_profile = (int***)calloc(ovlp_n, sizeof(int**));
    for (int i = 0; i < ovlp_n; ++i) {
        phase_set_to_hap_alle_profile[i] = (int**)calloc(3, sizeof(int*));
        for (int j = 0; j < 3; ++j) phase_set_to_hap_alle_profile[i][j] = (int*)calloc(var->n_uniq_alles, sizeof(int));
    }
    for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
        int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
        if (chunk->is_skipped[read_i]) continue;
        int hap = chunk->haps[read_i]; hts_pos_t phase_set = chunk->phase_sets[read_i];
        if (hap == 0 || phase_set == -1) continue; // skip unphased reads
        read_var_profile_t *p1 = p + read_i; int var_idx = var_i - p1->start_var_idx; int alle_i = p1->alleles[var_idx];
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "read %s %" PRIi64 " hap: %d phase_set: %" PRId64 " alle_i: %d\n", bam_get_qname(chunk->reads[read_i]), var->pos, hap, phase_set, alle_i);
        }
        if (alle_i != 1) alle_i = 0; // change all non-alt alleles to 0
        int phase_set_i = add_phase_set(phase_set, uniq_phase_set, &n_uniq_phase_set);
        phase_set_to_hap_alle_profile[phase_set_i][hap][alle_i] += 1;
    }
    // pick phase_set
    int phase_set_i = select_somatic_phase_set0(uniq_phase_set, n_uniq_phase_set, phase_set_to_hap_alle_profile, min_somatic_hap_depth);
    if (phase_set_i >= 0) {
        ps = uniq_phase_set[phase_set_i];
        int _alt_hap = 0; // alt haplotype
        for (int hap=1; hap <= 2; ++hap) {
            if (phase_set_to_hap_alle_profile[phase_set_i][hap][1]) {
                if (_alt_hap != 0) _alt_hap = 0; // more than one alt haplotype
                else _alt_hap = hap; // alt hap
            }
        }
        // check if all filters are passed
        if (_alt_hap != 0) { // no conflict haplotype
            cand_somatic_var_aux_info_t *aux_info = collect_somatic_var_aux_info(opt, chunk, ps, _alt_hap, var_i, ovlp_n, ovlp_b);
            if (phased_var_is_somatic(opt, chunk, var_i, aux_info)) *alt_hap = _alt_hap;
            var->somatic_aux_info = aux_info; // store aux info for somatic variant
        }
    } else if (phase_set_i == -2){ // no phase set to use: use simple count-based method, ignore phasing
        *alt_hap = 0;
        cand_somatic_var_aux_info_t *aux_info = collect_somatic_var_aux_info(opt, chunk, -1, -1, var_i, ovlp_n, ovlp_b);
        if (no_phase_var_is_somatic(opt, chunk, var_i, aux_info)) ps = 0;
        var->somatic_aux_info = aux_info; // store aux info for somatic variant
    }
    for (int i = 0; i < ovlp_n; ++i) {
        for (int j = 0; j < 3; ++j) free(phase_set_to_hap_alle_profile[i][j]);
        free(phase_set_to_hap_alle_profile[i]);
    } free(phase_set_to_hap_alle_profile); free(uniq_phase_set);
    if (ovlp_b != NULL) free(ovlp_b);
    return ps;
}

int assign_somatic_hap_based_on_phased_reads(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = cand_vars+var_i;
        var_init_hap_profile_cons_allele1(opt, var);
        read_var_profile_t *p = chunk->read_var_profile;
        // collect phase_set & alt_hap for candidate somatic variant
        hts_pos_t phase_set = -1; int alt_hap = 0; // haplotype with alt somatic var
        phase_set = select_somatic_phase_set_alt_hap(opt, chunk, var_i, p, &alt_hap); // phase_set & alt_hap: no confilct for alt_hap in phase_set
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "CandSomaticVar: %s\t%" PRIi64 "\t%d\t%c\t%d\t%d\t%d\t%ld\t%d\n", chunk->tname, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len, var->total_cov, var->alle_covs[1], phase_set, alt_hap);
        // update hap_to_cons_alle based on phase_set_to_use
        if (phase_set > 0 && alt_hap != 0) { // phased somatic var
            assert(alt_hap == 1 || alt_hap == 2);
            mark_somatic_var(var, phase_set, alt_hap);
        } else if (phase_set == 0) { // non-phased somatic var
            mark_somatic_var(var, phase_set, 2);
        } else {
            mark_non_somatic_var(var);
        }
    }
    // post-process somatic variants
    post_process_somatic_vars(opt, chunk, target_var_cate);
    return 0;
}