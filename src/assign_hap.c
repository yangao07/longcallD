#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include "math_utils.h"
#include <inttypes.h>
#include <stdlib.h>

extern int LONGCALLD_VERBOSE;

// init reads' hap & ps
void read_init_hap_phase_set(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->n_reads; ++i) {
        chunk->haps[i] = 0; chunk->phase_sets[i] = -1;
    }
}

static int get_var_max_cov_allele(cand_var_t *var) {
    int max_cov = 0, max_cov_alle_i = -1;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
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
void var_init_hap_profile_cons_allele(cand_var_t *cand_vars, int *var_idx, int n_cand_vars, int *var_i_to_cate) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
            var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
            var->hap_to_cons_alle[0] = get_var_max_cov_allele(var); // hom_idx
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        } else {
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j)
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            var->hap_to_cons_alle[0] = get_var_max_cov_allele(var); // hom_idx
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

void var_init_hap_profile_cons_allele1(cand_var_t *var) {
    if (var->hap_to_alle_profile == NULL) {
        var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
        for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
        var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
        var->hap_to_cons_alle[0] = get_var_max_cov_allele(var); // hom_idx
        for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
    } else {
        for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
        var->hap_to_cons_alle[0] = get_var_max_cov_allele(var); // hom_idx
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
        if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_SNP) {
            if (clean_het_snp_i == -1 || clean_het_snp_depth < cand_vars[var_i].total_cov) { 
                clean_het_snp_i = _var_i; clean_het_snp_depth = cand_vars[var_i].total_cov;
            }
        } else if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_INDEL) {
            if (clean_het_indel_i == -1 || clean_het_indel_depth < cand_vars[var_i].total_cov) { 
                clean_het_indel_i = _var_i; clean_het_indel_depth = cand_vars[var_i].total_cov;
            }
        } else if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR) {
            if (cand_vars[var_i].var_type == BAM_CDIFF) {
                if (noisy_het_snp_i == -1 || noisy_het_snp_depth < cand_vars[var_i].total_cov) {
                    noisy_het_snp_i = _var_i; noisy_het_snp_depth = cand_vars[var_i].total_cov;
                }
            } else {
                if (noisy_het_indel_i == -1 || noisy_het_indel_depth < cand_vars[var_i].total_cov) {
                    noisy_het_indel_i = _var_i; noisy_het_indel_depth = cand_vars[var_i].total_cov;
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
int init_assign_read_hap_based_on_cons_alle(int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int *hap_scores = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int n_vars_used[3] = {0, 0, 0};

    for (int var_i = p[read_i].start_var_idx; var_i <= p[read_i].end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) continue; 
        int read_var_idx = var_i - p[read_i].start_var_idx;
        if (p[read_i].alleles[read_var_idx] < 0) continue;
        int score0;
        cand_var_t *var = cand_vars+var_i;
        for (int hap = 1; hap <= 2; ++hap) {
            score0 = read_to_cons_allele_score(read_i, hap, var, var_i_to_cate[var_i], p[read_i].alleles[read_var_idx]);
            if (score0 != 0) n_vars_used[hap]++;
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
    else if (max_score > 0) return max_hap;
    else return 3-min_hap;
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
        if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CLEAN_HOM_VAR) continue; 
        int read_var_idx = var_i - p[read_i].start_var_idx;
        if (p[read_i].alleles[read_var_idx] < 0) continue;
        int score0;
        cand_var_t *var = cand_vars+var_i;
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

int update_var_hap_to_cons_alle(cand_var_t *var, int var_cate, int hap) {
    if (hap == 0) return 0;
    int max_cov = 0, max_cov_alle_i = -1;
    // int n_max_alt = 0;
    // if (var_cate & (LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_NOISY_CAND_HET_VAR)) n_max_alt = 1;
    // else if (var_cate & (LONGCALLD_CLEAN_HOM_VAR | LONGCALLD_NOISY_CAND_HOM_VAR)) n_max_alt = 2;

    for (int i = 0; i < var->n_uniq_alles; ++i) {
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
    var->hap_to_cons_alle[hap] = max_cov_alle_i;
    return 0;
}

int update_var_hap_profile_cons_alle_based_on_read_hap(int read_i, int hap, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
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
                update_var_hap_to_cons_alle(var, var_cate, i);
            }
        } else {
            var->hap_to_alle_profile[hap][allele_i] += 1;
            update_var_hap_to_cons_alle(var, var_cate, hap);
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

// reads->haps & haps_to_cons_alle are UpToDate
int iter_update_var_hap_cons_phase_set(bam_chunk_t *chunk, int *var_idx, read_var_profile_t *p, cand_var_t *cand_vars, int n_cand_vars, int *var_i_to_cate) {
    int *het_var_idx = (int*)malloc(n_cand_vars * sizeof(int));
    int n_het_vars = 0;
    int *is_het = (int*)calloc(n_cand_vars, sizeof(int));
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars+var_i;
        // XXX only use clean het vars
        // if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR) continue;
        if (var->hap_to_cons_alle[1] != -1 && var->hap_to_cons_alle[2] != -1 
            && var->hap_to_cons_alle[1] != var->hap_to_cons_alle[2]) {
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
        if (opt->is_ont) hap = init_assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
        else hap = assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
        if (hap == -1) hap = 0;
        chunk->haps[read_i] = hap;
        update_var_hap_profile_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
    }
    // update hap_to_cons_alle
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
            update_var_hap_to_cons_alle(var, var_i_to_cate[var_i], hap);
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

// read's PS was assigned using a fake HET (HOM) var
// goal: 1) assign haplotype + phase set to all reads
//       2) variant calling based on clustered reads
int assign_hap_based_on_het_vars_kmeans(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
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
    var_init_hap_profile_cons_allele(cand_vars, valid_var_idx, n_valid_vars, var_i_to_cate);

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
                int hap = init_assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
                if (hap == -1) { // no used vars, new phase set
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "new PS: %" PRIi64 " %s\n", cand_vars[var_i].pos, bam_get_qname(chunk->reads[read_i]));
                    hap = 1;
                }
                chunk->haps[read_i] = hap;
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read: %s hap: %d\n", bam_get_qname(chunk->reads[read_i]), hap);
                // update_var_hap_profile_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
                update_var_hap_profile_cons_alle_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
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
    for (int i = 0; i < n_uniq_phase_set; ++i) {
        int skip = 0, alt_hap = 0;
        for (int hap = 1; hap <= 2; ++hap) {
            // check if this hap has alt reads
            if (phase_set_to_hap_alle_profile[i][hap][1] > 0) {
                alt_hap++;
                // n_alt > 2*n_ref
                // if (phase_set_to_hap_alle_profile[i][hap][0]*2 < phase_set_to_hap_alle_profile[i][hap][1]) {
                    // skip = 1; break;
                // }
            }
            // hap depth
            if (phase_set_to_hap_alle_profile[i][hap][0] + phase_set_to_hap_alle_profile[i][hap][1] < min_somatic_hap_depth) {
                skip = 1; break;
            }
            //
        }
        if (skip == 0 && alt_hap == 1) { // XXX alt only show up in one haplotype
            ps_i = i; break;
        }
    }
    return ps_i;
}

cand_somatic_var_aux_info_t *init_somatic_aux_info(int max_alt_dp) {
    cand_somatic_var_aux_info_t *aux_info = (cand_somatic_var_aux_info_t*)malloc(sizeof(cand_somatic_var_aux_info_t));
    memset(aux_info, 0, sizeof(cand_somatic_var_aux_info_t));
    aux_info->alt_read_ids = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->alt_quals = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->min_win_quals = (int*)malloc(max_alt_dp * sizeof(int));
    aux_info->dis_to_indel_error = (int*)malloc(max_alt_dp * sizeof(int));
    return aux_info;
}

void free_somatic_var_aux_info(cand_somatic_var_aux_info_t *aux_info) {
    if (aux_info == NULL) return;
    free(aux_info->alt_read_ids);
    free(aux_info->alt_quals);
    free(aux_info->min_win_quals);
    free(aux_info->dis_to_indel_error);
    free(aux_info);
}

int get_min_dis_to_het_var(bam_chunk_t *chunk, int var_i) {
    int min_dis = INT32_MAX;
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    for (int i = var_i-1; i >= 0; --i) {
        if ((var_i_to_cate[i] & LONGCALLD_CAND_HET_VAR_CATE) == 0) continue;
        min_dis = MIN(min_dis, cand_vars[var_i].pos - cand_vars[i].pos - cand_vars[i].ref_len);
        break;
    }
    for (int i = var_i+1; i < chunk->n_cand_vars; ++i) {
        if ((var_i_to_cate[i] & LONGCALLD_CAND_HET_VAR_CATE) == 0) continue;
        min_dis = MIN(min_dis, cand_vars[i].pos - cand_vars[var_i].pos - cand_vars[var_i].ref_len);
        break;
    }
    return min_dis;
}

int get_read_win_min_qual(const call_var_opt_t *opt, digar_t *digar, int alt_qi) {
    uint8_t *qual = digar->qual; int qlen = digar->qlen;
    int flank_win_size = 3; // XXX opt->read_win_size; 5? 3+1+3
    int start = MAX(0, alt_qi - flank_win_size);
    int end = MIN(alt_qi + flank_win_size, qlen - 1);
    int min_qual = INT32_MAX; // default min qual
    for (int i = start; i <= end; ++i) {
        if (qual[i] < min_qual) min_qual = qual[i];
    }
    return min_qual;
}

int digar_is_var(bam_chunk_t *chunk, int var_i, digar1_t *digar1) {
    // check if the digar1 is a candidate variant
    if (digar1->type == BAM_CEQUAL) return 0; // skip equal bases
    for (int i = var_i; i < chunk->n_cand_vars; ++i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip germline vars
        if (var->var_type == digar1->type && digar1->pos == var->pos) return 1;
        else if (var->pos > digar1->pos) break; // no need to check further
    }
    for (int i = var_i-1; i >= 0; --i) {
        cand_var_t *var = chunk->cand_vars + i;
        if (!(chunk->var_i_to_cate[i] & LONGCALLD_CAND_GERMLINE_VAR_CATE)) continue; // skip germline vars
        if (var->var_type == digar1->type && digar1->pos == var->pos) return 1;
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
        if (alt_qi == digar1->qi && digar1->type == var_type && digar1->len == var_len) continue; // skip itself
        if (digar1->qi - alt_qi > dis_to_seq_error) break; // no need to check further
        if (digar_is_var(chunk, var_i, digar1)) continue;
        dis_to_seq_error = MIN(dis_to_seq_error, ABS(digar1->qi - alt_qi));
    }
    return dis_to_seq_error;
}
// int_cmp

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

// consider both SNV and INDEL
cand_somatic_var_aux_info_t *collect_somatic_var_aux_info(const call_var_opt_t *opt, bam_chunk_t *chunk, hts_pos_t ps, int alt_hap, int var_i, int64_t ovlp_n, int64_t *ovlp_b) {
    cand_somatic_var_aux_info_t *aux_info = init_somatic_aux_info(ovlp_n);
    cand_var_t *var = chunk->cand_vars + var_i; cgranges_t *read_var_cr = chunk->read_var_cr; read_var_profile_t *p = chunk->read_var_profile;
    int64_t ovlp_i;

    aux_info->min_dis_to_het_var = get_min_dis_to_het_var(chunk, var_i);
    int alt_i = 0;
    for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
        int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
        if (chunk->is_skipped[read_i]) continue;
        aux_info->total_dp += 1;
        int hap = chunk->haps[read_i]; hts_pos_t phase_set = chunk->phase_sets[read_i];
        if (hap != alt_hap || phase_set != ps) continue; // only collect reads with the same haplotype and phase set
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
            aux_info->min_win_quals[alt_i] = get_read_win_min_qual(opt, d, alt_qi);
            aux_info->dis_to_indel_error[alt_i] = get_dis_to_seq_error(chunk, var_i, d, var->var_type, var->alt_len, alt_qi, 1);
            // dis_to_homopolymer
            alt_i++;
        } else { // ref
            if (d->is_rev) aux_info->hap_ref_rev_cov += 1; else aux_info->hap_ref_for_cov += 1;
        }
    }
    return aux_info;
}

// check all filters
int var_is_somatic(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, cand_somatic_var_aux_info_t *aux_info) {
    cand_var_t *var = chunk->cand_vars + var_i;
    // if (var->var_type != BAM_CDIFF) return 1; // XXX only check SNVs
    int is_somatic = 1;
    // print aux_info
    // fprintf(stderr, "chr\tpos\ttotal_dp\thap_total_dp\thap_alt_dp\thap_ref_for_cov\thap_ref_rev_cov\thap_alt_for_cov\thap_alt_rev_cov\tstrand_fisher\t");
    // fprintf(stderr, "dis_to_het_var\tmedian_alt_qual\tmedian_min_win_qual\tmedian_dis_to_indel_error\n");

    fprintf(stderr, "%s\t%" PRIi64 "\t%d\t%d\t%d\t%d\t%d\t%d\t%d\t%.4f\t%d\t%d\t%d\t%d\n", chunk->tname, var->pos,
            aux_info->total_dp, aux_info->hap_total_dp, aux_info->hap_alt_dp, 
            aux_info->hap_ref_for_cov, aux_info->hap_ref_rev_cov, aux_info->hap_alt_for_cov, aux_info->hap_alt_rev_cov,
            log_betabinom_pmf(aux_info->hap_alt_dp, aux_info->hap_total_dp, opt->somatic_beta_alpha, opt->somatic_beta_beta, opt),
            // fisher_exact_test(aux_info->hap_ref_for_cov, aux_info->hap_ref_rev_cov, aux_info->hap_alt_for_cov, aux_info->hap_alt_rev_cov, opt),
            aux_info->min_dis_to_het_var, median_int(aux_info->alt_quals, aux_info->hap_alt_dp), median_int(aux_info->min_win_quals, aux_info->hap_alt_dp), median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp));

    // if (aux_info->hap_total_dp * opt->max_somatic_alt_af < aux_info->hap_alt_dp) return 0;
    if (aux_info->min_dis_to_het_var < opt->min_somatic_dis_to_het_var) return 0;
    if (median_int(aux_info->alt_quals, aux_info->hap_alt_dp) < opt->min_somatic_alt_qual) return 0;
    if (median_int(aux_info->min_win_quals, aux_info->hap_alt_dp) < opt->min_somatic_bq) return 0;
    if (median_int(aux_info->dis_to_indel_error, aux_info->hap_alt_dp) < opt->min_somatic_dis_to_seq_error) return 0;
    // for INDELs, not check beta-binomial
    if (var->var_type == BAM_CDIFF && log_betabinom_pmf(aux_info->hap_alt_dp, aux_info->hap_total_dp, opt->somatic_beta_alpha, opt->somatic_beta_beta, opt) < opt->min_somatic_log_beta_binom) return 0;
    // if (fisher_exact_test(aux_info->hap_ref_for_cov, aux_info->hap_ref_rev_cov, aux_info->hap_alt_for_cov, aux_info->hap_alt_rev_cov, opt) < opt->min_somatic_fisher_pval) return 0;
    return is_somatic;
}

hts_pos_t select_somatic_phase_set(const call_var_opt_t *opt, bam_chunk_t *chunk, int var_i, read_var_profile_t *p, int *alt_hap) { 
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
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read %d %" PRIi64 " hap: %d phase_set: %" PRId64 " alle_i: %d\n", read_i, var->pos, hap, phase_set, alle_i);
        if (alle_i != 1) alle_i = 0; // change all non-alt alleles to 0
        int phase_set_i = add_phase_set(phase_set, uniq_phase_set, &n_uniq_phase_set);
        phase_set_to_hap_alle_profile[phase_set_i][hap][alle_i] += 1;
    }
    int phase_set_i = select_somatic_phase_set0(uniq_phase_set, n_uniq_phase_set, phase_set_to_hap_alle_profile, min_somatic_hap_depth);
    if (phase_set_i != -1) {
        ps = uniq_phase_set[phase_set_i];
        int _alt_hap = 0; // alt haplotype
        for (int hap=1; hap <= 2; ++hap) {
            if (phase_set_to_hap_alle_profile[phase_set_i][hap][1]) {
                if (_alt_hap != 0) _alt_hap = 0; // more than one alt haplotype
                else _alt_hap = hap; // alt hap
            }
        }
        // check if all filters are passed
        if (_alt_hap != 0) {
            // if (var->var_type == BAM_CDIFF) { // only for somatic SNVs
            cand_somatic_var_aux_info_t *aux_info = collect_somatic_var_aux_info(opt, chunk, ps, _alt_hap, var_i, ovlp_n, ovlp_b);
            if (var_is_somatic(opt, chunk, var_i, aux_info)) *alt_hap = _alt_hap;
            var->somatic_aux_info = aux_info; // store aux info for somatic variant
                // free_somatic_var_aux_info(aux_info);
            // } else *alt_hap = _alt_hap; // for other types of variants, i.e., INDELs, we don't have somatic filters
        }
    }
    for (int i = 0; i < ovlp_n; ++i) {
        for (int j = 0; j < 3; ++j) free(phase_set_to_hap_alle_profile[i][j]);
        free(phase_set_to_hap_alle_profile[i]);
    } free(phase_set_to_hap_alle_profile); free(uniq_phase_set);
    if (ovlp_b != NULL) free(ovlp_b);
    return ps;
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

void mark_skip_somatic_reads(bam_chunk_t *chunk, int var_i) {
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
    for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
        int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
        // fprintf(stderr, "mark read %d as skipped for somatic variant calling\n", read_i);
        chunk->is_skipped_for_somatic[read_i] = 1; // tag read as skipped for somatic variant calling
    }
}

void mark_invalid_somatic(int somatic_dense_win, int somatic_dense_max_var, hts_pos_t *somatic_var_pos, int n_somatic_vars, int *is_invalid_somatic) {
    if (n_somatic_vars < somatic_dense_max_var) return; // no need to mark somatic variants
    // for any sliding window of size somatic_dense_win, if there are more than somatic_dense_max_var somatic variants, mark them as invalid
    for (int i = 0; i <= n_somatic_vars - somatic_dense_max_var; ++i) {
        int j = i + somatic_dense_max_var - 1;
        if (somatic_var_pos[j] - somatic_var_pos[i] < somatic_dense_win) {
            for (int k = i; k <= j; ++k) {
                is_invalid_somatic[k] = 1; // mark as invalid
            }
        }
    }
}

void mark_somatic_reads(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars;
    int somatic_dense_win = opt->somatic_dense_win; // default 1000-bp window
    int somatic_dense_max_var = opt->somatic_dense_win_max_vars; // maximum number of somatic variants in a dense window
    hts_pos_t *somatic_var_pos = (hts_pos_t*)malloc(chunk->n_cand_vars * sizeof(hts_pos_t));
    int *somatic_var_i = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    int *is_invalid_somatic = (int*)calloc(chunk->n_cand_vars, sizeof(int)); // mark somatic variants
    int n_somatic_vars = 0;
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue; // skip non-target variants
        cand_var_t *var = cand_vars + var_i;
        if (var->hap_to_cons_alle[1] == 0 && var->hap_to_cons_alle[2] == 0) continue; // skip variants without alt alleles
        somatic_var_i[n_somatic_vars] = var_i; // collect somatic variant indices
        somatic_var_pos[n_somatic_vars++] = var->pos; // collect somatic variant positions
    }
    mark_invalid_somatic(somatic_dense_win, somatic_dense_max_var, somatic_var_pos, n_somatic_vars, is_invalid_somatic);
    for (int i = 0; i < n_somatic_vars; ++i) {
        int var_i = somatic_var_i[i];
        if (is_invalid_somatic[i] == 1) mark_skip_somatic_reads(chunk, var_i);
    }
    free(somatic_var_pos); free(somatic_var_i); free(is_invalid_somatic);
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
int assign_somatic_hap_based_on_phased_reads(const call_var_opt_t *opt, bam_chunk_t *chunk, int target_var_cate) {
    int *var_i_to_cate = chunk->var_i_to_cate; cand_var_t *cand_vars = chunk->cand_vars; cgranges_t *read_var_cr = chunk->read_var_cr;
    for (int var_i = 0; var_i < chunk->n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = cand_vars+var_i;
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "CandSomaticVar: %s\t%" PRIi64 "\t%d\t%d\t%d\t%d\n", chunk->tname, var->pos, var->ref_len, var->var_type, var->alt_len, var->total_cov);
        var_init_hap_profile_cons_allele1(var);
        read_var_profile_t *p = chunk->read_var_profile;
        // collect phase_set for somatic variant calling
        hts_pos_t phase_set = -1; int alt_hap = 0; // haplotype with alt somatic var
        phase_set = select_somatic_phase_set(opt, chunk, var_i, p, &alt_hap);
        // update hap_to_cons_alle based on phase_set_to_use
        if (phase_set != -1 && alt_hap != 0) {
            assert(alt_hap == 1 || alt_hap == 2);
            mark_somatic_var(var, phase_set, alt_hap);
        } else {
            mark_non_somatic_var(var);
        }
    }
    // post-process somatic variants
    post_process_somatic_vars(opt, chunk, target_var_cate);
    return 0;
}