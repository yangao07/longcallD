#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include <inttypes.h>

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
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        } else {
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j)
                memset(var->hap_to_alle_profile[j], 0, var->n_uniq_alles * sizeof(int));
            var->hap_to_cons_alle[0] = get_var_max_cov_allele(var); // hom_idx
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = 1;
            } else {
                for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) var->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

void var_init_hap_to_alle_profile(cand_var_t *cand_vars, int *var_idx, int n_cand_vars) {
    for (int _var_i = 0; _var_i < n_cand_vars; ++_var_i) {
        int var_i = var_idx[_var_i]; cand_var_t *var = cand_vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
        } else {
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
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
int assign_read_hap_based_on_cons_alle(int read_i, cand_var_t *cand_vars, read_var_profile_t *p, int *var_i_to_cate, int target_var_cate) {
    int *hap_scores = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
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
    // int n_max_alt = 0;
    // if (var_cate & (LONGCALLD_CLEAN_HET_SNP | LONGCALLD_CLEAN_HET_INDEL | LONGCALLD_NOISY_CAND_HET_VAR)) n_max_alt = 1;
    // else if (var_cate & (LONGCALLD_CAND_HOM_VAR | LONGCALLD_NOISY_CAND_HOM_VAR)) n_max_alt = 2;

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

int iter_update_var_hap_to_cons_alle(bam_chunk_t *chunk, int *var_idx, int n_cand_vars, read_var_profile_t *p, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate) {
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
        // hap == -1 XXX
        int hap = assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
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
int assign_hap_based_on_het_vars_kmeans(bam_chunk_t *chunk, int target_var_cate, call_var_opt_t *opt) {
    read_var_profile_t *p = chunk->read_var_profile;
    int n_valid_vars = 0, *valid_var_idx = (int*)malloc(chunk->n_cand_vars * sizeof(int));
    int *var_is_valid = (int*)calloc(chunk->n_cand_vars, sizeof(int));
    for (int i = 0; i < chunk->n_cand_vars; ++i) {
        if ((chunk->var_i_to_cate[i] & target_var_cate) == 0) continue;
        valid_var_idx[n_valid_vars++] = i;
        var_is_valid[i] = 1;
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
            if (var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HOM_VAR || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) // || var_i_to_cate[var_i] == LONGCALLD_NOISY_CAND_HET_VAR)
                continue;
            cand_var_t *var = cand_vars+var_i;
            if (LONGCALLD_VERBOSE >= 2)
                fprintf(stderr, "var_i %d %ld %d-%c-%d %d reads\n", var_i, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], var->alt_len, var->total_cov);
            ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
            for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) { // assign each read a haplotype if it is not assigned yet
                int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
                if (chunk->is_skipped[read_i] || chunk->haps[read_i] != 0) continue;
                // this round of assignment may not be correct for all reads, will be updated in the following rounds
                // 1/2: hap, 0: tied, no hap, -1: no var can be used
                int hap = assign_read_hap_based_on_cons_alle(read_i, cand_vars, p, var_i_to_cate, target_var_cate);
                if (hap == -1) { // no used vars, new phase set
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "new PS: %ld %s\n", cand_vars[var_i].pos, bam_get_qname(chunk->reads[read_i]));
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
        int changed_hap2 = iter_update_var_hap_to_cons_alle(chunk, valid_var_idx, n_valid_vars, p, cand_vars, var_i_to_cate, target_var_cate);
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