#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include <inttypes.h>

extern int LONGCALLD_VERBOSE;


void read_init_hap_ps(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->n_reads; ++i) {
        chunk->haps[i] = 0; chunk->PS[i] = 0;
    }
}

// all candidate vars, including heterozygous and homozygous SNPs:
// alle_to_hap: allele: 0/1 -> HAP: 0
// hap_to_alle_profile: HAP 0/1/2: -> allele: 0/1 -> count: 0
// hap_to_cons_alle: HAP 0/1/2: -1
void var_init_hap_profile_cons_allele(cand_var_t *cand_vars, int *var_idx, int n_cand_vars, int *var_i_to_cate) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->alle_to_hap = (uint8_t*)calloc(var->n_uniq_alles, sizeof(uint8_t)); // +1: minor_alt_allele
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
            var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        } else {
            memset(var->alle_to_hap, 0, var->n_uniq_alles * sizeof(uint8_t));
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j)
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

void var_init_hap_profile(cand_var_t *cand_vars, int *var_idx, int n_cand_vars) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->alle_to_hap = (uint8_t*)calloc(var->n_uniq_alles, sizeof(uint8_t)); // +1: minor_alt_allele
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
        } else {
            memset(var->alle_to_hap, 0, var->n_uniq_alles * sizeof(uint8_t));
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            }
        }
    }
}


void var_init_hap_cons_alle0(cand_var_t *var, int min_alt_dp) {
    // select the most common allele as the consensus allele, based on hap_to_alle_profile
    for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
        int max_cov = 0, max_cov_alle_i = -1;
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (var->hap_to_alle_profile[hap][i] > max_cov) {
                max_cov = var->hap_to_alle_profile[hap][i];
                if (max_cov >= min_alt_dp) max_cov_alle_i = i;
            }
        }
        var->hap_to_cons_alle[hap] = max_cov_alle_i;
        var->is_skipped = 0;
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

int update_var_hap_profile_based_on_changed_hap(int new_hap, int old_hap, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate, int min_alt_dp, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        if (var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR) continue;
        int read_var_idx = var_i - start_var_idx;
        // if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
        cand_var_t *var = cand_vars+var_i;
        // only update Het-Var
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "pos: %" PRId64 ", old_hap: %d, new_hap: %d, var: %d\n", var->pos, old_hap, new_hap, allele_i);
        if (old_hap == 0) {
            var->hap_to_alle_profile[1][allele_i] -= 1;
            var->hap_to_alle_profile[2][allele_i] -= 1;
        } else var->hap_to_alle_profile[old_hap][allele_i] -= 1;
        if (new_hap == 0) {
            var->hap_to_alle_profile[1][allele_i] += 1;
            var->hap_to_alle_profile[2][allele_i] += 1;
        } else var->hap_to_alle_profile[new_hap][allele_i] += 1;
        var_init_hap_cons_alle0(var, min_alt_dp);
    }
    return 0;
}

// int assign_read_hap_based_on_pre_reads(int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
int assign_read_hap_phase_set_based_on_cons_alle(int read_i, hts_pos_t *phase_set, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int *hap_scores = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    *phase_set = -1;
    int n_vars_used[3] = {0, 0, 0};

    for (int var_i = p[read_i].start_var_idx; var_i <= p[read_i].end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) continue; 
        int read_var_idx = var_i - p[read_i].start_var_idx;
        if (p[read_i].alleles[read_var_idx] == -1) continue;
        int score0;
        cand_var_t *var = cand_vars+var_i;
        for (int hap = 1; hap <= 2; ++hap) {
            score0 = read_to_cons_allele_score(read_i, hap, var, var_i_to_cate[var_i], p[read_i].alleles[read_var_idx]);
            if (score0 != 0) n_vars_used[hap]++;
            if (*phase_set == -1) *phase_set = var->phase_set;
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

int update_var_hap_to_cons_alle(cand_var_t *var, int var_cate, int hap) {
    if (hap == 0) return 0;
    int max_cov = 0, max_cov_alle_i = -1;
    int n_max_alt = 0;
    if (var_cate & (LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_NOISY_CAND_HET_VAR)) n_max_alt = 1;
    else if (var_cate & (LONGCALLD_CAND_HOM_VAR | LONGCALLD_NOISY_CAND_HOM_VAR)) n_max_alt = 2;

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
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "var_i: %d, %ld allele_i: %d\n", var_i, cand_vars[var_i].pos, allele_i);
        if (allele_i == -1) continue;
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
        if (allele_i == -1) continue;
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
    if (allele_i1 == -1 || allele_i2 == -1) return -1;

    int agree = 0, conflict = 0, agree_hap = hap, conflict_hap = 3-hap;
    if (vars[var1].hap_to_cons_alle[agree_hap] == allele_i1 && vars[var2].hap_to_cons_alle[agree_hap] == allele_i2) agree = 1;
    if (vars[var1].hap_to_cons_alle[hap] == allele_i1 && vars[var2].hap_to_cons_alle[conflict_hap] == allele_i2) conflict = 1;

    if (agree) return 1; else if (conflict) return 0;
    else return -1;
}

// reads->haps/haps_to_cons_alle are UpToDate
int update_var_hap_cons_phase_set(bam_chunk_t *chunk, int *var_idx, read_var_profile_t *p, cand_var_t *cand_vars, int n_cand_vars) {
    int *het_var_idx = (int*)malloc(n_cand_vars * sizeof(int));
    int n_het_vars = 0;
    int *is_het = (int*)calloc(n_cand_vars, sizeof(int));
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i];
        cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_cons_alle[1] != -1 && var->hap_to_cons_alle[2] != -1 
            && var->hap_to_cons_alle[1] != var->hap_to_cons_alle[2]) {
            // fprintf(stderr, "is_het: %d %ld %d-%c-%d\n", var_i, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len);
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
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "var_i %d %ld %d-%c-%d %d reads\n", var_i, var->pos, var->ref_len, BAM_CIGAR_STR[cand_vars[var_i].var_type], cand_vars[var_i].alt_len, cand_vars[var_i].total_cov);
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
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%ld %d %d\n", var->pos, n_agree[_var_i], n_conflict[_var_i]);
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


int collect_tmp_hap_cons_allele_by_deduct_read(cand_var_t *var, int hap, int allele_i, int *tmp_hap_to_cons_alle) {
    for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
        tmp_hap_to_cons_alle[i] = -1;
        if (i != hap || var->hap_to_cons_alle[i] != allele_i) { // no change
            tmp_hap_to_cons_alle[i] = var->hap_to_cons_alle[i];
        } else { // i==hap && var->hap_to_cons_alle[i] == allele_i
            int cov, max_cov = 0, max_cov_alle_i = -1;
            for (int j = 0; j < var->n_uniq_alles; ++j) {
                cov = var->hap_to_alle_profile[i][j];
                if (j == allele_i) cov -= 1;
                if (cov > max_cov) {
                    max_cov = cov; max_cov_alle_i = j;
                }
            }
            if (max_cov_alle_i == -1 && LONGCALLD_VERBOSE >= 2) _err_func_printf("No HAP %d allele in SNP: %" PRId64 ", %d-%c\n", i, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
            tmp_hap_to_cons_alle[i] = max_cov_alle_i;
        }
    }
    return 0;
}

// update hap_cons_profile based on hap_to_base_profile
// update haplotype for a read based on SNP profiles of all other overlapping reads
// input: target read, SNP profiles of all reads
// output: updated haplotype of the target read
int update_var_aln_hap1(int target_read_i, hts_pos_t *phase_set, int cur_hap,  bam_chunk_t *chunk, read_var_profile_t *p, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate) {
    int start_var_idx = p[target_read_i].start_var_idx, end_var_idx = p[target_read_i].end_var_idx;
    // deduct target read from hap_to_alle_profile, then compare target read's var profile with hap_cons_alle
    int *hap_match_cnt = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int *tmp_hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
    *phase_set = -1;

    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        cand_var_t *var = cand_vars+var_i;
        int var_weight = 1; // XXX
        if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_SNP)  var_weight = 4;
        else if (var_i_to_cate[var_i] == LONGCALLD_CLEAN_HET_INDEL)  var_weight = 3;
        else if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR) {
            if (var->var_type == BAM_CDIFF) var_weight = 2;
            else var_weight = 1;
        } 

        // if (p[target_read_i].var_is_used[read_var_idx] == 0) continue;
        if (var->is_skipped) continue;
        int allele_i = p[target_read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
        collect_tmp_hap_cons_allele_by_deduct_read(var, cur_hap, allele_i, tmp_hap_to_cons_alle);
        if (*phase_set == -1) *phase_set = var->phase_set;
        for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
            if (tmp_hap_to_cons_alle[i] == allele_i) {
                hap_match_cnt[i] += var_weight;
            } else if (tmp_hap_to_cons_alle[i] != -1) {
                hap_match_cnt[i] -= var_weight;
            }
        }
    }
    int max_cnt=0, max_hap=0, sec_cnt=0;
    for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
        if (hap_match_cnt[i] > max_cnt) {
            sec_cnt = max_cnt;
            max_cnt = hap_match_cnt[i]; max_hap = i;
        } else if (hap_match_cnt[i] > sec_cnt) {
            sec_cnt = hap_match_cnt[i];
        }
    }
    free(hap_match_cnt); free(tmp_hap_to_cons_alle);
    if (max_cnt == 0) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == 0 (pos: %" PRId64 ")\n", bam_get_qname(chunk->reads[target_read_i]), chunk->reads[target_read_i]->core.pos);
        return 0; // unknown
    } else if (max_cnt == sec_cnt) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == sec_cnt (%" PRId64 ")\n", bam_get_qname(chunk->reads[target_read_i]), chunk->reads[target_read_i]->core.pos);
        max_hap = cur_hap;
    }
    return max_hap;
}

int update_read_hap_phase_set(const call_var_opt_t *opt, bam_chunk_t *chunk, read_var_profile_t *p, cand_var_t *cand_vars, int n_cand_vars, int *var_i_to_cate, int target_var_cate) {
    // read-wise loop
    // re-calculate read-wise haplotype (Var-wise 1. Hap, 2. Allele, 3. Read, 4. AlleleCons are all up-to-date)
    int changed = 0;
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i]) continue;
        int cur_hap = chunk->haps[read_i];
        // XXX TODO: potential local optima
        hts_pos_t phase_set = -1;
        int new_hap = update_var_aln_hap1(read_i, &phase_set, cur_hap, chunk, p, cand_vars, var_i_to_cate, target_var_cate);
        if (new_hap != cur_hap) { // update bam_haps, hap_to_base_profile, hap_to_cons_base
            if (LONGCALLD_VERBOSE >= 2) {
                fprintf(stderr, "read (%d): %s, pos: %" PRId64 "\t", read_i, bam_get_qname(chunk->reads[read_i]), chunk->reads[read_i]->core.pos);
                fprintf(stderr, "\t\t cur_hap: %d, new_hap: %d\n", cur_hap, new_hap);
            }
            changed = 1;
            chunk->haps[read_i] = new_hap; // update intermediately
            // bam_aux_append(chunk->reads[read_i], "XT", 'i', 4, (uint8_t*)&(chunk->haps[read_i]));
            update_var_hap_profile_based_on_changed_hap(new_hap, cur_hap, cand_vars, var_i_to_cate, target_var_cate, opt->min_alt_dp, p, read_i);
        }
        if (phase_set != -1) chunk->PS[read_i] = phase_set;
    }
    return changed;
}

int update_hap_to_cons_alle(bam_chunk_t *chunk, int *var_idx, int n_cand_vars, read_var_profile_t *p, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate) {
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
    var_init_hap_profile(cand_vars, var_idx, n_cand_vars);
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i]) continue;
        // hap == -1 ?? XXX
        hts_pos_t phase_set = -1;
        int hap = assign_read_hap_phase_set_based_on_cons_alle(read_i, &phase_set, cand_vars, p, var_i_to_cate, target_var_cate);
        if (hap == -1) hap = 0;
        chunk->haps[read_i] = hap; chunk->PS[read_i] = phase_set;
        // fprintf(stderr, "update read : %s hap: %d\n", bam_get_qname(chunk->reads[read_i]), hap);
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



// goal: 1) assign haplotype + phase set to all reads
//       2) variant calling based on clustered reads
int assign_hap_based_on_het_vars_kmeans(bam_chunk_t *chunk, int target_var_cate, call_var_opt_t *opt) {
    read_var_profile_t *p = chunk->read_var_profile;
    int n_valid_vars = 0, *valid_var_idx = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        if ((chunk->var_i_to_cate[i] & target_var_cate) == 0) continue;
        valid_var_idx[n_valid_vars++] = i;
    }
    int *var_i_to_cate = chunk->var_i_to_cate;
    cand_var_t *cand_vars = chunk->cand_vars;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;

    read_init_hap_ps(chunk);
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
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR)
                continue;
            cand_var_t *var = cand_vars+var_i;
            if (LONGCALLD_VERBOSE >= 2)
                fprintf(stderr, "var_i %d %ld %d-%c-%d %d reads\n", var_i, cand_vars[var_i].pos, cand_vars[var_i].ref_len, 
                                BAM_CIGAR_STR[cand_vars[var_i].var_type], cand_vars[var_i].alt_len, cand_vars[var_i].total_cov);
            ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
            for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) { // assign each read a haplotype if it is not assigned yet
                int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
                if (chunk->is_skipped[read_i] || chunk->haps[read_i] != 0) continue;
                // this round of assignment may not be correct for all reads, will be updated in the following rounds
                hts_pos_t phase_set = -1;
                int hap = assign_read_hap_phase_set_based_on_cons_alle(read_i, &phase_set, cand_vars, p, var_i_to_cate, target_var_cate);
                if (hap == -1) { // no used vars, new phase set
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "new PS: %ld %s\n", cand_vars[var_i].pos, bam_get_qname(chunk->reads[read_i]));
                    hap = 1;
                }
                chunk->haps[read_i] = hap; chunk->PS[read_i] = phase_set;
                if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "read: %s hap: %d\n", bam_get_qname(chunk->reads[read_i]), hap);
                // update_var_hap_profile_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
                update_var_hap_profile_cons_alle_based_on_read_hap(read_i, hap, cand_vars, p, var_i_to_cate, target_var_cate);
            }
        } free(var_ii); free(ovlp_b);
    }
    // update hap_to_cons_alle for both het and hom vars
    // for (int _var_i = 0; _var_i < n_valid_vars; ++_var_i) {
    //     int var_i = valid_var_idx[_var_i];
    //     for (int hap = 0; hap <= LONGCALLD_DEF_PLOID; ++hap) {
    //         update_var_hap_to_cons_alle(cand_vars+var_i, hap);
    //     }
    // }
    // after 1st round, read_haps/hap_to_alle_profile/hap_to_cons_alle are upToDate and will be used in the following rounds
    int max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "iter: %d\n", i_iter);
        // XXX var-wise loop update PS
        // read-wise loop: update read_hap
        int changed_hap1 = update_var_hap_cons_phase_set(chunk, valid_var_idx, p, cand_vars, n_valid_vars);
        int changed_hap2 = update_hap_to_cons_alle(chunk, valid_var_idx, n_valid_vars, p, cand_vars, var_i_to_cate, target_var_cate);
        // int changed_hap2 = update_read_hap_phase_set(chunk, valid_var_idx, n_valid_vars, p, cand_vars, var_i_to_cate, target_var_cate);
        if (changed_hap1 == 0 && changed_hap2 == 0) break;
        else {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "changed_hap1: %d changed_hap2: %d\n", changed_hap1, changed_hap2);
        }
    }
    free(valid_var_idx);
    return 0;
}

// alle_to_hap: is only used for HET vars, not used for HOM vars
int collect_max_cov_allele(cand_var_t *var) {
    int max_cov = 0, max_cov_alle_i = -1;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        if (var->alle_covs[i] > max_cov) {
            max_cov = var->alle_covs[i]; max_cov_alle_i = i;
        }
    }
    return max_cov_alle_i;
}

hts_pos_t assign_var_init_hom_hap(cand_var_t *var) {
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Init Hom-Var hap: %" PRId64 ", %d-%c\n", var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        var->alle_to_hap[i] = 0;
    }
    return var->phase_set;
}

// 1st round operations: update base_to_hap -> {1/2/0}
// assign haplotype to a SNP, when no other information avaliable
// most common base -> 1, second common base -> 2, others -> 0
// potential start of a PhaseSet, could be merged with others (intra- or inter-blocks)
hts_pos_t assign_var_init_hap(cand_var_t *var) {
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Init Het-Var hap: %" PRId64 ", %d-%c\n", var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
    if (var->var_type == BAM_CDIFF) var->phase_set = var->pos; // potential start of a PhaseSet
    else var->phase_set = var->pos - 1;
    int hap1_alle_i = -1, hap2_alle_i = -1, hap1_cov = 0, hap2_cov = 0;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        var->alle_to_hap[i] = 0;
        if (var->alle_covs[i] > hap1_cov) {
            hap2_alle_i = hap1_alle_i; hap2_cov = hap1_cov;
            hap1_alle_i = i; hap1_cov = var->alle_covs[i];
        } else if (var->alle_covs[i] > hap2_cov) {
            hap2_alle_i = i; hap2_cov = var->alle_covs[i];
        }
    }
    if (hap1_alle_i == -1) {
        fprintf(stderr, "No candidate allele in Var: %d-%" PRId64 ", %d-%c\n", var->tid, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
        return -1;
    }
    if (hap2_alle_i == -1) { // only one allele: low coverage or homozygous
        var->alle_to_hap[hap1_alle_i] = 1; // var->alle_to_hap[hap1_alle_i] = 2;
    } else {
        if (hap2_alle_i == 0) { // ref allele
            var->alle_to_hap[hap1_alle_i] = 2; var->alle_to_hap[hap2_alle_i] = 1;
        } else {
            var->alle_to_hap[hap1_alle_i] = 1; var->alle_to_hap[hap2_alle_i] = 2;
        }
    }
    return var->phase_set;
}


// assign haplotype to a SNP based on read SNP profiles
// 1. pick the most common haplotype and corresponding most common base 
// 2. assign most common base to the most common haplotype, and other bases to the other haplotype
hts_pos_t assign_var_hap_based_on_pre_reads1(cand_var_t *var) {
    int first_hap=0, first_hap_cnt=0, first_hap_alle_i=-1;
    int sec_hap=0, sec_hap_cnt=0, sec_hap_alle_i=-1;

    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < var->n_uniq_alles; ++i) { // ref + alt_allele; minor_alt_allele is not considered
            if (var->hap_to_alle_profile[j][i] > first_hap_cnt) {
                sec_hap_cnt = first_hap_cnt; sec_hap = first_hap; sec_hap_alle_i = first_hap_alle_i;
                first_hap_cnt = var->hap_to_alle_profile[j][i]; first_hap = j; first_hap_alle_i = i;
            } else if (var->hap_to_alle_profile[j][i] > sec_hap_cnt) {
                sec_hap_cnt = var->hap_to_alle_profile[j][i]; sec_hap = j; sec_hap_alle_i = i;
            }
        }
    }
    if (first_hap == 0) _err_error_exit("major haplotype is not set yet\n"); 
    if (sec_hap == 0) {
        if (first_hap == 1) sec_hap = 2; else sec_hap = 1;
        // set the most common bases other than first_hap_base to sec_hap
        int sec_hap_alle_cov = 0;
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (var->alle_covs[i] > sec_hap_alle_cov && i != first_hap_alle_i) {
                sec_hap_alle_i = i; sec_hap_alle_cov = var->alle_covs[i];
            }
        }
    }
    if (first_hap == sec_hap) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Var: %" PRId64 ", %d-%c, first_hap: %d (%d: %d), sec_hap: %d (%d: %d)\n", var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], first_hap, first_hap_alle_i, first_hap_cnt, sec_hap, sec_hap_alle_i, sec_hap_cnt);
        return assign_var_init_hap(var);
    // } else if (first_hap_alle_i == sec_hap_alle_i) { // homozygous
    } else {
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (i == first_hap_alle_i) var->alle_to_hap[i] = first_hap;
            else if (i == sec_hap_alle_i) var->alle_to_hap[i] = sec_hap;
            else var->alle_to_hap[i] = 0;
        }
    }
    return var->phase_set;
} 

// after a read is assigned with hap, update hap of all other SNPs covered by this read
// including homozygous and heterozygous SNPs
// void update_var_hap_profile_based_on_aln_hap(int hap, hts_pos_t phase_set, cand_var_t *var, int *var_i_to_cate, int target_var_cate, read_var_profile_t *p, int read_i) {
// haps of vars may be conflicted, will be updated in the following rounds
void update_var_hap_profile_based_on_aln_hap(int hap, cand_var_t *var, int *var_i_to_cate, int target_var_cate, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        // if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
        var[var_i].hap_to_alle_profile[hap][allele_i] += 1;
        // update HOM var's phase_set as well ??
        // if (var[var_i].phase_set == 0 || phase_set <= var[var_i].pos)
            // var[var_i].phase_set = phase_set;
    }
}

void init_hap_cons_alle(const call_var_opt_t *opt, bam_chunk_t *chunk, read_var_profile_t *p, cand_var_t *cand_vars, int n_cand_vars, int *var_i_to_cate, int target_var_cate) {
    // update hap_to_cons_profile for reads with HAP==0
    for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
        if (chunk->is_skipped[read_i] || chunk->haps[read_i] != 0) continue;
        for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
            update_var_hap_profile_based_on_aln_hap(hap, cand_vars, var_i_to_cate, target_var_cate, p, read_i);
        }
    }
    // init hap_to_cons_alle for all vars, including HOM and HET vars
    for (int var_i = 0; var_i < n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = cand_vars+var_i;
        var_init_hap_cons_alle0(var, opt->min_alt_dp);
    }
}

int var_hap_profile_cov(cand_var_t *var) {
    int cov = 0;
    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            cov += var->hap_to_alle_profile[j][i];
        }
    }
    return cov;
}

// hts_pos_t assign_var_hap_based_on_pre_reads(cand_var_t *var, int min_dp) {
int assign_var_hap_based_on_pre_reads(cand_var_t *var, int min_dp) {
    if (var_hap_profile_cov(var) < min_dp) 
        return assign_var_init_hap(var);
    else {
        return assign_var_hap_based_on_pre_reads1(var);
    }
}


int *sort_cand_vars(cand_var_t *cand_vars, int n_cand_vars, int *var_i_to_cate, int target_var_cate) {
    int *sorted_cand_var_i = (int*)malloc(n_cand_vars * sizeof(int));
    for (int i = 0; i < n_cand_vars; ++i) {
        sorted_cand_var_i[i] = i;
    }
    // sort var by type: clean SNP -> clean indel -> noisy SNP -> noisy INDEL
    for (int i = 0; i < n_cand_vars-1; ++i) {
        for (int j = i+1; j < n_cand_vars; ++j) {
            if (var_i_to_cate[sorted_cand_var_i[i]] > var_i_to_cate[sorted_cand_var_i[j]]) {
                int tmp = sorted_cand_var_i[i]; sorted_cand_var_i[i] = sorted_cand_var_i[j]; sorted_cand_var_i[j] = tmp;
            } else if (var_i_to_cate[sorted_cand_var_i[i]] == var_i_to_cate[sorted_cand_var_i[j]]) {
                if (cand_vars[sorted_cand_var_i[i]].var_type != cand_vars[sorted_cand_var_i[j]].var_type) {
                    if (cand_vars[sorted_cand_var_i[j]].var_type == BAM_CDIFF && cand_vars[sorted_cand_var_i[i]].var_type != BAM_CDIFF) {
                        int tmp = sorted_cand_var_i[i]; sorted_cand_var_i[i] = sorted_cand_var_i[j]; sorted_cand_var_i[j] = tmp;
                    }
                } else {
                    if (cand_vars[sorted_cand_var_i[i]].pos > cand_vars[sorted_cand_var_i[j]].pos) {
                        int tmp = sorted_cand_var_i[i]; sorted_cand_var_i[i] = sorted_cand_var_i[j]; sorted_cand_var_i[j] = tmp;
                    } else if (cand_vars[sorted_cand_var_i[i]].pos == cand_vars[sorted_cand_var_i[j]].pos) {
                        if (cand_vars[sorted_cand_var_i[i]].ref_len > cand_vars[sorted_cand_var_i[j]].ref_len) {
                            int tmp = sorted_cand_var_i[i]; sorted_cand_var_i[i] = sorted_cand_var_i[j]; sorted_cand_var_i[j] = tmp;
                        }
                    }
                }
            }
        }
    }
    return sorted_cand_var_i;
}

// XXX hap_cons_base not upated yet
// start with the very first SNP
// 1st round: assign hap to SNP's base, then assign hap to reads covering this SNP,
// 2~N rounds: re-assign haplotypes to reads based on haplotype clusters in previous rounds, 
//             until no changes to any reads
// output chunk->haps[i] to 1 or 2, 0: unknown
int assign_hap_based_on_het_vars(bam_chunk_t *chunk, int target_var_cate, call_var_opt_t *opt) {
    read_var_profile_t *p = chunk->read_var_profile;
    // int n_cand_vars = chunk->var_cate_counts[target_var_cate];
    int n_cand_vars = chunk->n_cand_vars;
    // int *var_cate_idx = chunk->var_cate_idx[target_var_cate];
    int *var_i_to_cate = chunk->var_i_to_cate;
    cand_var_t *cand_vars = chunk->cand_vars;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;

    int n_valid_vars = 0, *valid_var_idx = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        if ((chunk->var_i_to_cate[i] & target_var_cate) == 0) continue;
        valid_var_idx[n_valid_vars++] = i;
    }

    // 1st loop: var-wise loop
    read_init_hap_ps(chunk);
    var_init_hap_profile_cons_allele(cand_vars, valid_var_idx, n_valid_vars, var_i_to_cate);
    // sort var by type: clean SNP -> clean indel -> noisy SNP -> noisy INDEL
    int *sorted_cand_var_i = sort_cand_vars(cand_vars, n_cand_vars, var_i_to_cate, target_var_cate);
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = sorted_cand_var_i[_var_i];
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        if (var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR) {
            assign_var_init_hom_hap(cand_vars+var_i);
            continue;
        }
        cand_var_t *var = cand_vars+var_i;
        ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
        // here we highly rely on reads that were assigned with HAPs previously, >= 2 reads is enough
        hts_pos_t phase_set = assign_var_hap_based_on_pre_reads(var, 2); //opt->min_dp); // update alle_to_hap
        if (phase_set == -1) continue;
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
            int read_var_idx = var_i - p[read_i].start_var_idx;
            if (chunk->haps[read_i] == 0 && chunk->is_skipped[read_i] != 1) { // && p[read_i].var_is_used[read_var_idx] == 1) {
                int var_alle_i = p[read_i].alleles[read_var_idx];
                if (var_alle_i == -1) continue;
                int hap = var->alle_to_hap[var_alle_i];
                // XXX for hap == 0, due to alle_to_hap was not updated yet
                if (hap != 0) {
                    // first time assign hap to the read (update bam_haps)
                    chunk->haps[read_i] = hap;
                    if (LONGCALLD_VERBOSE >= 2)
                        fprintf(stderr, "read: %s, cur_var: %" PRId64 ", %d-%c, alle: %d, hap: %d\n", bam_get_qname(chunk->reads[read_i]), var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var_alle_i, hap);
                    // update hap_to_alle_profile for all Vars covered by this read, based on its assigned haplotype
                    // udpated profile will then be used for following Vars (assign_var_hap_based_on_pre_reads)
                    // update_var_hap_profile_based_on_aln_hap(hap, phase_set, cand_vars, var_i_to_cate, target_var_cate, p, read_i);
                    // XXX only reads with HAP==1 or HAP==2 will be updated
                    update_var_hap_profile_based_on_aln_hap(hap, cand_vars, var_i_to_cate, target_var_cate, p, read_i);
                    // update PS for the read
                    if (chunk->PS[read_i] == 0 || chunk->PS[read_i] > phase_set) chunk->PS[read_i] = phase_set;
                }
            }
        }
        // var_init_hap_cons_alle0(var, opt->min_alt_dp); // update hap_to_cons_alle
    } 
    init_hap_cons_alle(opt, chunk, p, cand_vars, n_cand_vars, var_i_to_cate, target_var_cate);
    
    // after first round, they are upToDate and will be used in the following rounds
        // bam_haps: 0/1/2
        // hap_to_alle_profile: HAP 1/2 -> allele 0/1 -> count, no HAP 0
        // hap_to_cons_alle: HAP 1/2 -> allele 0/1, no HAP 0
        // alle_to_hap will not be used (may be NOT upToDate)
    // 2nd loop: read-wise iterative loop
    int changed_hap1, changed_hap2, max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "iter: %d\n", i_iter);
        // for each var: flip hap_to_cons_alle if n_agree < n_conflict, update phase_set
        changed_hap1 = update_var_hap_cons_phase_set(chunk, valid_var_idx, p, cand_vars, n_valid_vars);
        // for each read: update hap & phase set
        // re-assign read to hap, check if cons_alle changed
        // changed_hap2 = update_hap_to_cons_alle(chunk, valid_var_idx, n_valid_vars, p, cand_vars, var_i_to_cate, target_var_cate);
        changed_hap2 = update_read_hap_phase_set(opt, chunk, p, cand_vars, n_cand_vars, var_i_to_cate, target_var_cate);
        if (changed_hap1 == 0 && changed_hap2 == 0) break;
        else {
            if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "changed_hap1: %d, changed_hap2: %d\n", changed_hap1, changed_hap2);
        }
    }
    if (LONGCALLD_VERBOSE >= 2) _err_info("Iteration: %d\n", i_iter);
    free(ovlp_b); free(sorted_cand_var_i); free(valid_var_idx);
    return 0;
}