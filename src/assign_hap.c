#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var.h"
#include "bam_utils.h"
#include "collect_var.h"

extern int LONGCALLD_VERBOSE;
// 1st round operations: update base_to_hap -> {1/2/0}
// assign haplotype to a SNP, when no other information avaliable
// most common base -> 1, second common base -> 2, others -> 0
// potential start of a PhaseSet, could be merged with others (intra- or inter-blocks)
hts_pos_t assign_snp_init_hap(cand_snp_t *snp) {
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Init SNP hap: %ld\n", snp->pos);
    snp->phase_set = snp->pos; // potential start of a PhaseSet
    // if (snp->pos == 10737188)
        // printf("ok");
    int hap1_base_i = -1, hap2_base_i = -1, hap1_cov = 0, hap2_cov = 0;
    for (int i = 0; i < snp->n_uniq_bases; ++i) {
        snp->base_to_hap[i] = 0;
        if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
        if (snp->base_covs[i] > hap1_cov) {
            hap2_base_i = hap1_base_i; hap2_cov = hap1_cov;
            hap1_base_i = i; hap1_cov = snp->base_covs[i];
        } else if (snp->base_covs[i] > hap2_cov) {
            hap2_base_i = i; hap2_cov = snp->base_covs[i];
        }
    }
    if (hap2_base_i == -1) _err_error_exit("Only one base in SNP: %ld\n", snp->pos);
    if (snp->bases[hap2_base_i] == LONGCALLD_BAM_REF_BASE_IDX) {
        snp->base_to_hap[hap1_base_i] = 2; snp->base_to_hap[hap2_base_i] = 1;
    } else {
        snp->base_to_hap[hap1_base_i] = 1; snp->base_to_hap[hap2_base_i] = 2;
    }
    return snp->phase_set;
}

hts_pos_t assign_var_init_hap(cand_var_t *var) {
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Init Var hap: %ld\n", var->pos);
    var->phase_set = var->pos; // potential start of a PhaseSet
    // if (snp->pos == 10737188)
        // printf("ok");
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
    if (hap2_alle_i == -1) _err_error_exit("Only one allele in Var: %ld\n", var->pos);
    if (hap2_alle_i == 0) { // ref allele
        var->alle_to_hap[hap1_alle_i] = 2; var->alle_to_hap[hap2_alle_i] = 1;
    } else {
        var->alle_to_hap[hap1_alle_i] = 1; var->alle_to_hap[hap2_alle_i] = 2;
    }
    return var->phase_set;
}

// assign haplotype to a SNP based on read SNP profiles
// 1. pick the most common haplotype and corresponding most common base 
// 2. assign most common base to the most common haplotype, and other bases to the other haplotype
hts_pos_t assign_snp_hap_based_on_pre_reads1(cand_snp_t *snp) {
    int first_hap=0, first_hap_cnt=0, first_hap_base_i=-1;
    int sec_hap=0, sec_hap_cnt=0, sec_hap_base_i=-1;

    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
            if (snp->hap_to_base_profile[j][i] > first_hap_cnt) {
                sec_hap_cnt = first_hap_cnt; sec_hap = first_hap; sec_hap_base_i = first_hap_base_i;
                first_hap_cnt = snp->hap_to_base_profile[j][i]; first_hap = j; first_hap_base_i = i;
            } else if (snp->hap_to_base_profile[j][i] > sec_hap_cnt) {
                sec_hap_cnt = snp->hap_to_base_profile[j][i]; sec_hap = j; sec_hap_base_i = i;
            }
        }
    }
    if (first_hap == 0) _err_error_exit("major haplotype is not set yet\n"); 
    if (sec_hap == 0) {
        if (first_hap == 1) sec_hap = 2; else sec_hap = 1;
        // set the most common bases other than first_hap_base to sec_hap
        int sec_hap_base_cov = 0;
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
            if (snp->base_covs[i] > sec_hap_base_cov && i != first_hap_base_i) {
                sec_hap_base_i = i; sec_hap_base_cov = snp->base_covs[i];
            }
        }
    }
    if (first_hap == sec_hap || first_hap_base_i == sec_hap_base_i) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("SNP: %ld, first_hap: %d (%c: %d), sec_hap: %d (%c: %d)\n", snp->pos, first_hap, LONGCALLD_BAM_BASE_STR[snp->bases[first_hap_base_i]], first_hap_cnt, sec_hap, LONGCALLD_BAM_BASE_STR[snp->bases[sec_hap_base_i]], sec_hap_cnt);
        assign_snp_init_hap(snp);
    } else {
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            if (i == first_hap_base_i) snp->base_to_hap[i] = first_hap;
            else if (i == sec_hap_base_i) snp->base_to_hap[i] = sec_hap;
            else snp->base_to_hap[i] = 0;
        }
    }
    return snp->phase_set;
} 

hts_pos_t assign_var_hap_based_on_pre_reads1(cand_var_t *var) {
    int first_hap=0, first_hap_cnt=0, first_hap_alle_i=-1;
    int sec_hap=0, sec_hap_cnt=0, sec_hap_alle_i=-1;

    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < var->n_uniq_alles; ++i) {
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
    if (first_hap == sec_hap || first_hap_alle_i == sec_hap_alle_i) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Var: %ld, first_hap: %d (%d: %d), sec_hap: %d (%d: %d)\n", var->pos, first_hap, first_hap_alle_i, first_hap_cnt, sec_hap, sec_hap_alle_i, sec_hap_cnt);
        assign_var_init_hap(var);
    } else {
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (i == first_hap_alle_i) var->alle_to_hap[i] = first_hap;
            else if (i == sec_hap_alle_i) var->alle_to_hap[i] = sec_hap;
            else var->alle_to_hap[i] = 0;
        }
    }
    return var->phase_set;
} 

void snp_init_hap_profile(cand_snp_t *snps, int n_cand_snps) {
    for (int i = 0; i < n_cand_snps; ++i) {
        cand_snp_t *snp = snps+i;
        if (snp->hap_to_base_profile == NULL) {
            snp->base_to_hap = (uint8_t*)calloc(snp->n_uniq_bases, sizeof(uint8_t));
            snp->hap_to_base_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) snp->hap_to_base_profile[i] = (int*)calloc(snp->n_uniq_bases, sizeof(int));
            snp->hap_to_cons_base = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
        }
    }
}

void var_init_hap_profile(cand_var_t *vars, int n_cand_vars) {
    for (int i = 0; i < n_cand_vars; ++i) {
        cand_var_t *var = vars+i;
        if (var->hap_to_alle_profile == NULL) {
            var->alle_to_hap = (uint8_t*)calloc(var->n_uniq_alles, sizeof(uint8_t));
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
            var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
        }
    }
}

void snp_init_hap_cons_base(cand_snp_t *snps, int n_cand_snps) {
    for (int snp_i = 0; snp_i < n_cand_snps; ++snp_i) {
        cand_snp_t *snp = snps+snp_i;
        for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
            // select the most common base as the consensus base, based on hap_to_base_profile
            int max_cov = 0, max_cov_base_i = -1;
            for (int i = 0; i < snp->n_uniq_bases; ++i) {
                if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
                if (snp->hap_to_base_profile[j][i] > max_cov) {
                    max_cov = snp->hap_to_base_profile[j][i]; max_cov_base_i = i;
                }
            }
            if (max_cov_base_i == -1) {
                _err_func_printf("No HAP %d base in SNP: %ld\n", j, snp->pos);
                snp->is_skipped = 1;
            } else {
                snp->is_skipped = 0;
                snp->hap_to_cons_base[j] = max_cov_base_i;
            }
        }
    }
}

void snp_init_hap_cons_base0(cand_snp_t *snp) {
    // select the most common base as the consensus base, based on hap_to_base_profile
    for (int hap = 1; hap <= LONGCALLD_DEF_PLOID ; ++hap) {
        int max_cov = 0, max_cov_base_i = -1;
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
            if (snp->hap_to_base_profile[hap][i] > max_cov) {
                max_cov = snp->hap_to_base_profile[hap][i];
                max_cov_base_i = i;
            }
        }
        if (max_cov_base_i == -1) {
            _err_func_printf("No HAP %d base in SNP: %ld\n", hap, snp->pos);
            snp->is_skipped = 1;
        } else {
            snp->hap_to_cons_base[hap] = max_cov_base_i;
            snp->is_skipped = 0;
        }
    }
}

void var_init_hap_cons_alle0(cand_var_t *var) {
    // select the most common allele as the consensus allele, based on hap_to_alle_profile
    for (int hap = 1; hap <= LONGCALLD_DEF_PLOID ; ++hap) {
        int max_cov = 0, max_cov_alle_i = -1;
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (var->hap_to_alle_profile[hap][i] > max_cov) {
                max_cov = var->hap_to_alle_profile[hap][i];
                max_cov_alle_i = i;
            }
        }
        if (max_cov_alle_i == -1) {
            if (LONGCALLD_VERBOSE >= 2) _err_func_printf("No HAP allele %d in Var: %ld\n", hap, var->pos);
            var->is_skipped = 1;
        } else {
            var->hap_to_cons_alle[hap] = max_cov_alle_i;
            var->is_skipped = 0;
        }
    }
}

void snp_init_hap_cons_base1(cand_snp_t *snp, int hap) {
    // select the most common base as the consensus base, based on hap_to_base_profile
    int max_cov = 0, max_cov_base_i = -1;
    for (int i = 0; i < snp->n_uniq_bases; ++i) {
        if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX)
            continue;
        if (snp->hap_to_base_profile[hap][i] > max_cov) {
            max_cov = snp->hap_to_base_profile[hap][i];
            max_cov_base_i = i;
        }
    }
    if (max_cov_base_i == -1) _err_func_printf("No HAP %d base in SNP: %ld\n", hap, snp->pos);
    snp->hap_to_cons_base[hap] = max_cov_base_i;
}

int snp_hap_profile_cov(cand_snp_t *snp) {
    int cov = 0;
    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            if (snp->bases[i] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
            cov += snp->hap_to_base_profile[j][i];
        }
    }
    return cov;
}

hts_pos_t assign_snp_hap_based_on_pre_reads(cand_snp_t *snp) {
    if (snp_hap_profile_cov(snp) < LONGCALLD_MIN_CAND_SNP_DP) 
        return assign_snp_init_hap(snp);
    else {
        return assign_snp_hap_based_on_pre_reads1(snp);
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

hts_pos_t assign_var_hap_based_on_pre_reads(cand_var_t *var) {
    if (var_hap_profile_cov(var) < LONGCALLD_MIN_CAND_SNP_DP) 
        return assign_var_init_hap(var);
    else {
        return assign_var_hap_based_on_pre_reads1(var);
    }
}

void reset_snp_hap_profile(cand_snp_t *snp) {
    for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
        for (int i = 0; i < snp->n_uniq_bases; ++i) {
            snp->hap_to_base_profile[j][i] = 0;
        }
    }
}

// after a read is assigned with hap, update hap of all other SNPs covered by this read
void update_snp_hap_profile_based_on_aln_hap(int hap, hts_pos_t ps, cand_snp_t *snp, read_snp_profile_t *p, int read_i) {
    int start_snp_idx = p[read_i].start_snp_idx, end_snp_idx = p[read_i].end_snp_idx;
    for (int snp_i = start_snp_idx; snp_i <= end_snp_idx; ++snp_i) {
        int read_snp_idx = snp_i - start_snp_idx;
        if (p[read_i].snp_is_used[read_snp_idx] == 0) continue;
        uint8_t snp_base = p[read_i].snp_bases[read_snp_idx];
        int snp_base_i = snp[snp_i].base_to_i[snp_base];
        snp[snp_i].hap_to_base_profile[hap][snp_base_i] += 1;
        if (snp[snp_i].phase_set == 0 || ps <= snp[snp_i].pos)
            snp[snp_i].phase_set = ps;
    }
}

void update_var_hap_profile_based_on_aln_hap(int hap, hts_pos_t phase_set, cand_var_t *var, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        int read_var_idx = var_i - start_var_idx;
        if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        var[var_i].hap_to_alle_profile[hap][allele_i] += 1;
        if (var[var_i].phase_set == 0 || phase_set <= var[var_i].pos)
            var[var_i].phase_set = phase_set;
    }
}

// 2-n round operations
// assign haplotype to all other SNPs covered by a read
//  void assign_hap_to_read_snps(read_snp_profile_t *p, int read_i, cand_snp_t *snps, int hap) {
//     for (int snp_i = p[read_i].start_snp_idx; snp_i <= p[read_i].end_snp_idx; ++snp_i) {
//         if (snps[snp_i]->base_to_hap[0] == -1) 
//     }
//  }

// if the first read' hap is 2, then flip all haps
void flip_aln_hap(bam_chunk_t *bam_chunk) {
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->haps[i] == 1 || bam_chunk->haps[i] == 2) {
            if (bam_chunk->haps[i] == 1) return;
            else break;
        } else continue;
    }
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->haps[i] == 1) bam_chunk->haps[i] = 2;
        else if (bam_chunk->haps[i] == 2) bam_chunk->haps[i] = 1;
    }
}

int flip_snp_hap(cand_snp_t *snp) {
    return 0;
}

int collect_tmp_hap_cons_by_deduct_read(cand_snp_t *snp, int hap, int snp_base_i, int *tmp_hap_to_cons_base) {
    for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
        if (i != hap || snp->hap_to_cons_base[i] != snp_base_i) { // no change
            tmp_hap_to_cons_base[i] = snp->hap_to_cons_base[i];
        } else { // i==hap && snp->hap_to_cons_base[i] == snp_base_i
            int cov, max_cov = 0, max_cov_base_i = -1;
            for (int j = 0; j < snp->n_uniq_bases; ++j) {
                if (snp->bases[j] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
                cov = snp->hap_to_base_profile[i][j];
                if (j == snp_base_i) cov -= 1;
                if (cov > max_cov) {
                    max_cov = cov; max_cov_base_i = j;
                }
            }
            if (max_cov_base_i == -1) _err_func_printf("No HAP %d base in SNP: %ld\n", i, snp->pos);
            tmp_hap_to_cons_base[i] = max_cov_base_i;
        }
    }
    return 0;
}

int collect_tmp_hap_cons_allele_by_deduct_read(cand_var_t *var, int hap, int allele_i, int *tmp_hap_to_cons_alle) {
    for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
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
            if (max_cov_alle_i == -1 && LONGCALLD_VERBOSE >= 2) _err_func_printf("No HAP %d allele in SNP: %ld\n", i, var->pos);
            tmp_hap_to_cons_alle[i] = max_cov_alle_i;
        }
    }
    return 0;
}

// update hap_cons_profile based on hap_to_base_profile
// update haplotype for a read based on SNP profiles of all other overlapping reads
// input: target read, SNP profiles of all reads
// output: updated haplotype of the target read
int update_aln_hap1(int target_read_i, int cur_hap,  bam_chunk_t *bam_chunk, read_snp_profile_t *p, cand_snp_t *cand_snps) {
    int start_snp_idx = p[target_read_i].start_snp_idx, end_snp_idx = p[target_read_i].end_snp_idx;
    // deduct target read's SNP profile from hap_cons_profile, then compare target read's SNP profile with hap_cons_profile
    int *hap_match_cnt = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int *tmp_hap_to_cons_base = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));

    for (int snp_i = start_snp_idx; snp_i <= end_snp_idx; ++snp_i) {
        int read_snp_idx = snp_i - start_snp_idx;
        if (p[target_read_i].snp_is_used[read_snp_idx] == 0) continue;

        cand_snp_t *snp = cand_snps+snp_i;
        int snp_base_i = snp->base_to_i[p[target_read_i].snp_bases[read_snp_idx]];
        collect_tmp_hap_cons_by_deduct_read(snp, cur_hap, snp_base_i, tmp_hap_to_cons_base);
        for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
            if (tmp_hap_to_cons_base[i] == snp_base_i) {
                hap_match_cnt[i] += 1;
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
    free(hap_match_cnt); free(tmp_hap_to_cons_base);
    if (max_cnt == 0) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == 0 (pos: %ld)\n", bam_get_qname(bam_chunk->reads[target_read_i]), bam_chunk->reads[target_read_i]->core.pos);
        return 0; // unknown
    } else if (max_cnt == sec_cnt) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == sec_cnt (%ld)\n", bam_get_qname(bam_chunk->reads[target_read_i]), bam_chunk->reads[target_read_i]->core.pos);
        max_hap = cur_hap;
    }
    return max_hap;
}

int update_var_aln_hap1(int target_read_i, int cur_hap,  bam_chunk_t *bam_chunk, read_var_profile_t *p, cand_var_t *cand_vars) {
    int start_var_idx = p[target_read_i].start_var_idx, end_var_idx = p[target_read_i].end_var_idx;
    // deduct target read from hap_to_alle_profile, then compare target read's var profile with hap_cons_alle
    int *hap_match_cnt = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int *tmp_hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));

    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        int read_var_idx = var_i - start_var_idx;
        if (p[target_read_i].var_is_used[read_var_idx] == 0) continue;

        cand_var_t *var = cand_vars+var_i;
        int allele_i =p[target_read_i].alleles[read_var_idx];
        collect_tmp_hap_cons_allele_by_deduct_read(var, cur_hap, allele_i, tmp_hap_to_cons_alle);
        for (int i = 1; i <= LONGCALLD_DEF_PLOID; ++i) {
            if (tmp_hap_to_cons_alle[i] == allele_i) {
                hap_match_cnt[i] += 1;
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
            _err_func_printf("Read %s max_cnt == 0 (pos: %ld)\n", bam_get_qname(bam_chunk->reads[target_read_i]), bam_chunk->reads[target_read_i]->core.pos);
        return 0; // unknown
    } else if (max_cnt == sec_cnt) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == sec_cnt (%ld)\n", bam_get_qname(bam_chunk->reads[target_read_i]), bam_chunk->reads[target_read_i]->core.pos);
        max_hap = cur_hap;
    }
    return max_hap;
}

int update_snp_hap_profile_based_on_changed_hap(int new_hap, int old_hap, cand_snp_t *cand_snps, read_snp_profile_t *p, int read_i) {
    int start_snp_idx = p[read_i].start_snp_idx, end_snp_idx = p[read_i].end_snp_idx;
    for (int snp_i = start_snp_idx; snp_i <= end_snp_idx; ++snp_i) {
        int read_snp_idx = snp_i - start_snp_idx;
        if (p[read_i].snp_is_used[read_snp_idx] == 0) continue;
        uint8_t snp_base = p[read_i].snp_bases[read_snp_idx];
        cand_snp_t *snp = cand_snps+snp_i;
        int snp_base_i = snp->base_to_i[snp_base];
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "pos: %ld, old_hap: %d, new_hap: %d, snp: %d\n", snp->pos, old_hap, new_hap, snp_base_i);
        snp->hap_to_base_profile[old_hap][snp_base_i] -= 1;
        snp->hap_to_base_profile[new_hap][snp_base_i] += 1;
        snp_init_hap_cons_base0(snp);
    }
    return 0;
}

int update_var_hap_profile_based_on_changed_hap(int new_hap, int old_hap, cand_var_t *cand_vars, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        int read_var_idx = var_i - start_var_idx;
        if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        cand_var_t *var = cand_vars+var_i;
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "pos: %ld, old_hap: %d, new_hap: %d, var: %d\n", var->pos, old_hap, new_hap, allele_i);
        var->hap_to_alle_profile[old_hap][allele_i] -= 1;
        var->hap_to_alle_profile[new_hap][allele_i] += 1;
        var_init_hap_cons_alle0(var);
    }
    return 0;
}

// XXX hap_cons_base not upated yet
void update_snp_hap_cons_base1(int read_i, int new_hap, bam_chunk_t *bam_chunk, read_snp_profile_t *p, cand_snp_t *cand_snps) {
    int start_snp_idx = p[read_i].start_snp_idx, end_snp_idx = p[read_i].end_snp_idx;
    for (int snp_i = start_snp_idx; snp_i <= end_snp_idx; ++snp_i) {
        int snp_base_i = cand_snps[snp_i].base_to_i[p[read_i].snp_bases[snp_i-start_snp_idx]];
        cand_snp_t *snp = cand_snps+snp_i;
        snp->hap_to_base_profile[new_hap][snp_base_i] += 1;
        snp_init_hap_cons_base1(snp, new_hap);
    }
}

int *sort_snps_by_cov(int n_cand_snps, cand_snp_t *cand_snps) {
    int *ordered_snps = (int *)malloc(n_cand_snps * sizeof(int));
    for (int i = 0; i < n_cand_snps; ++i) ordered_snps[i] = i;
    for (int i = 0; i < n_cand_snps; ++i) {
        for (int j = i + 1; j < n_cand_snps; ++j) {
            if (cand_snps[ordered_snps[i]].n_depth < cand_snps[ordered_snps[j]].n_depth) {
                int tmp = ordered_snps[i];
                ordered_snps[i] = ordered_snps[j];
                ordered_snps[j] = tmp;
            }
        }
    }
    return ordered_snps;
}

// start with the very first SNP
// 1st round: assign hap to SNP's base, then assign hap to reads covering this SNP,
// 2~N rounds: re-assign haplotypes to reads based on haplotype clusters in previous rounds, 
//             until no changes to any reads
// output bam_chunk->haps[i] to 1 or 2, 0: unknown
char read_name[1024] = "m84039_231005_222902_s1/80479720/ccs";

int assign_hap_based_on_cand_vars(bam_chunk_t *bam_chunk) {
    read_var_profile_t *p = bam_chunk->read_var_profile;
    int n_cand_vars = bam_chunk->n_cand_vars;
    cand_var_t *cand_vars = bam_chunk->cand_vars;
    // init read_snp_cr for read overlapping query
    cgranges_t *read_var_cr = cr_init(); int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        if (p[i].start_var_idx < 0 || p[i].end_var_idx < 0) continue;
        cr_add(read_var_cr, "cr", p[i].start_var_idx, p[i].end_var_idx+1, i); // [start, end): 0-based
    } cr_index(read_var_cr);

    // 1st loop: var-wise loop
    var_init_hap_profile(cand_vars, n_cand_vars);
    for (int var_i = 0; var_i < n_cand_vars; ++var_i) {
        cand_var_t *var = cand_vars+var_i;
        ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
        hts_pos_t phase_set = assign_var_hap_based_on_pre_reads(var); // update alle_to_hap
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
            int read_var_idx = var_i - p[read_i].start_var_idx;
            if (bam_chunk->haps[read_i] == 0 && bam_chunk->is_skipped[read_i] != 1 && p[read_i].var_is_used[read_var_idx] == 1) {
                int var_alle_i = p[read_i].alleles[read_var_idx];
                int hap = cand_vars[var_i].alle_to_hap[var_alle_i];
                // XXX for hap == 0, due to alle_to_hap was not updated yet
                if (hap != 0) {
                    // first time assign hap to the read (update bam_haps)
                    bam_chunk->haps[read_i] = hap;
                    if (LONGCALLD_VERBOSE >= 2)
                        fprintf(stderr, "read: %s, cur_var: %ld, alle: %d, hap: %d\n", bam_get_qname(bam_chunk->reads[read_i]), var->pos, var_alle_i, hap);
                    // update hap_to_alle_profile for all Vars covered by this read, based on its assigned haplotype
                    // udpated profile will then be used for following Vars (assign_var_hap_based_on_pre_reads)
                    update_var_hap_profile_based_on_aln_hap(hap, phase_set, cand_vars, p, read_i);
                }
            }
        }
        var_init_hap_cons_alle0(var); // update hap_to_cons_alle
    } // after first round, 
      // bam_haps/hap_to_alle_profile/hap_to_cons_alle are upToDate and will be used in the following rounds
      // alle_to_hap will not be used (may be NOT upToDate)
    // 2nd loop: read-wise iterative loop
    int changed_hap, max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "iter: %d\n", i_iter);
        changed_hap = 0;
        // read-wise loop
        // re-calculate read-wise haplotype (Var-wise 1. Hap, 2. Allele, 3. Read, 4. AlleleCons are all up-to-date)
        for (int read_i = 0; read_i < bam_chunk->n_reads; ++read_i) {
            if (bam_chunk->is_skipped[read_i] || bam_chunk->haps[read_i] == 0) continue;
            // if (strcmp(read_name, bam_get_qname(bam_chunk->reads[read_i])) == 0) 
                // fprintf(stderr, "ok\n");
            int cur_hap = bam_chunk->haps[read_i];
            // XXX TODO: potential local optima
            int new_hap = update_var_aln_hap1(read_i, cur_hap, bam_chunk, p, cand_vars);
            if (new_hap != cur_hap) { // update bam_haps, hap_to_base_profile, hap_to_cons_base
                if (LONGCALLD_VERBOSE >= 2) {
                    fprintf(stderr, "read (%d): %s, pos: %ld\t", read_i, bam_get_qname(bam_chunk->reads[read_i]), bam_chunk->reads[read_i]->core.pos);
                    fprintf(stderr, "\t\t cur_hap: %d, new_hap: %d\n", cur_hap, new_hap);
                }
                changed_hap = 1;
                bam_chunk->haps[read_i] = new_hap; // update intermediately
                // bam_aux_append(bam_chunk->reads[read_i], "XT", 'i', 4, (uint8_t*)&(bam_chunk->haps[read_i]));
                update_var_hap_profile_based_on_changed_hap(new_hap, cur_hap, cand_vars, p, read_i);
            }
        } if (changed_hap == 0) break;
    }
    if (LONGCALLD_VERBOSE >= 2) _err_info("Iteration: %d\n", i_iter);
    free(ovlp_b); cr_destroy(read_var_cr);
    return 0;
}

int assign_hap(read_snp_profile_t *p, int n_cand_snps, cand_snp_t *cand_snps, bam_chunk_t *bam_chunk) {
    // init read_snp_cr for read overlapping query
    cgranges_t *read_snp_cr = cr_init(); int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        if (p[i].start_snp_idx < 0 || p[i].end_snp_idx < 0) continue;
        cr_add(read_snp_cr, "cr", p[i].start_snp_idx, p[i].end_snp_idx+1, i); // [start, end): 0-based
    } cr_index(read_snp_cr);

    // 1st loop: SNP wise loop
    snp_init_hap_profile(cand_snps, n_cand_snps);
    for (int snp_i = 0; snp_i < n_cand_snps; ++snp_i) {
        cand_snp_t *snp = cand_snps+snp_i;
        ovlp_n = cr_overlap(read_snp_cr, "cr", snp_i, snp_i+1, &ovlp_b, &max_b);
        hts_pos_t ps = assign_snp_hap_based_on_pre_reads(snp); // update base_to_hap
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int read_i = cr_label(read_snp_cr, ovlp_b[ovlp_i]);
            int read_snp_idx = snp_i - p[read_i].start_snp_idx;
            if (bam_chunk->haps[read_i] == 0 && bam_chunk->is_skipped[read_i] != 1 && p[read_i].snp_is_used[read_snp_idx] == 1) {
                int snp_base_i = cand_snps[snp_i].base_to_i[p[read_i].snp_bases[read_snp_idx]];
                int hap = cand_snps[snp_i].base_to_hap[snp_base_i];
                // XXX for hap == 0, due to base_to_hap was not updated yet
                if (hap != 0) {
                    // first time assign hap to the read (update bam_haps)
                    bam_chunk->haps[read_i] = hap;
                    if (LONGCALLD_VERBOSE >= 2)
                        fprintf(stderr, "read: %s, cur_snp: %ld, base: %d(%c), hap: %d\n", bam_get_qname(bam_chunk->reads[read_i]), snp->pos, p[read_i].snp_bases[snp_i-p[read_i].start_snp_idx], LONGCALLD_BAM_BASE_STR[p[read_i].snp_bases[snp_i-p[read_i].start_snp_idx]], hap);
                    // update hap_to_base_profile for all SNPs covered by this read, based on its assigned haplotype
                    // udpated profile will then be used for following SNPs (assign_snp_hap_based_on_pre_reads, L273)
                    update_snp_hap_profile_based_on_aln_hap(hap, ps, cand_snps, p, read_i);
                }
            }
        }
        snp_init_hap_cons_base0(snp); // update hap_to_cons_base
    } // after first round, 
      // bam_haps/hap_to_base_profile/hap_to_cons_base are upToDate and will be used in the following rounds
      // base_to_hap will not be used (may be NOT upToDate)
    // 2nd loop: read-wise iterative loop
    int changed_hap, max_iter = 10, i_iter=0;
    while (i_iter++ < max_iter) {
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "iter: %d\n", i_iter);
        changed_hap = 0;
        // read-wise loop, or SNP-wise loop again?
        // re-calculate read-wise haplotype (SNP-wise 1. Hap, 2. Base, 3. Read, 4. Cons are all up-to-date)
        for (int read_i = 0; read_i < bam_chunk->n_reads; ++read_i) {
            if (bam_chunk->is_skipped[read_i] || bam_chunk->haps[read_i] == 0) continue;
            // if (strcmp(read_name, bam_get_qname(bam_chunk->reads[read_i])) == 0) 
                // fprintf(stderr, "ok\n");
            int cur_hap = bam_chunk->haps[read_i];
            // XXX TODO: potential local optima
            int new_hap = update_aln_hap1(read_i, cur_hap, bam_chunk, p, cand_snps);
            if (new_hap != cur_hap) { // update bam_haps, hap_to_base_profile, hap_to_cons_base
                if (LONGCALLD_VERBOSE >= 2) {
                    fprintf(stderr, "read (%d): %s, pos: %ld\t", read_i, bam_get_qname(bam_chunk->reads[read_i]), bam_chunk->reads[read_i]->core.pos);
                    fprintf(stderr, "\t\t cur_hap: %d, new_hap: %d\n", cur_hap, new_hap);
                }
                changed_hap = 1;
                bam_chunk->haps[read_i] = new_hap; // update intermediately
                // bam_aux_append(bam_chunk->reads[read_i], "XT", 'i', 4, (uint8_t*)&(bam_chunk->haps[read_i]));
                update_snp_hap_profile_based_on_changed_hap(new_hap, cur_hap, cand_snps, p, read_i);
            }
        } if (changed_hap == 0) break;
    }
    if (LONGCALLD_VERBOSE >= 2) _err_info("Iteration: %d\n", i_iter);
    free(ovlp_b); cr_destroy(read_snp_cr);
    return 0;
}