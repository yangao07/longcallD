#include "assign_hap.h"
#include "cgranges.h"
#include "utils.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include <inttypes.h>

extern int LONGCALLD_VERBOSE;

int collect_max_cov_allele(cand_var_t *var) {
    int max_cov = 0, max_cov_alle_i = -1;
    for (int i = 0; i < var->n_uniq_alles; ++i) {
        if (var->alle_covs[i] > max_cov) {
            max_cov = var->alle_covs[i]; max_cov_alle_i = i;
        }
    }
    return max_cov_alle_i;
}

// 1st round operations: update base_to_hap -> {1/2/0}
// assign haplotype to a SNP, when no other information avaliable
// most common base -> 1, second common base -> 2, others -> 0
// potential start of a PhaseSet, could be merged with others (intra- or inter-blocks)
hts_pos_t assign_var_init_hap(cand_var_t *var) {
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "Init Var hap: %" PRId64 ", %d-%c\n", var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
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
    if (hap1_alle_i == -1) _err_error_exit("No candidate allele in Var: %d-%" PRId64 ", %d-%c\n", var->tid, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
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
    if (first_hap == sec_hap || first_hap_alle_i == sec_hap_alle_i) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Var: %" PRId64 ", %d-%c, first_hap: %d (%d: %d), sec_hap: %d (%d: %d)\n", var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type], first_hap, first_hap_alle_i, first_hap_cnt, sec_hap, sec_hap_alle_i, sec_hap_cnt);
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

// all candidate vars, including heterozygous and homozygous SNPs:
// init haplotype profile: HAP: 0 -> max_cons_allele_i, 1/2: -1
void var_init_hap_profile(cand_var_t *vars, int n_cand_vars, int *var_i_to_cate, int target_var_cate) {
// void var_init_hap_profile(cand_var_t *vars, int n_cand_vars, int *var_cate_list) {
    // for (int i = 0; i < n_cand_vars; ++i) {
        // int var_i = var_cate_list[i];
    for (int var_i = 0; var_i < n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        cand_var_t *var = vars+var_i;
        if (var->hap_to_alle_profile == NULL) {
            var->alle_to_hap = (uint8_t*)calloc(var->n_uniq_alles, sizeof(uint8_t)); // +1: minor_alt_allele
            var->hap_to_alle_profile = (int**)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int*));
            for (int i = 0; i <= LONGCALLD_DEF_PLOID; ++i) var->hap_to_alle_profile[i] = (int*)calloc(var->n_uniq_alles, sizeof(int));
            var->hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));
            var->hap_to_cons_alle[0] = collect_max_cov_allele(var);
            for (int j = 1; j <= LONGCALLD_DEF_PLOID; ++j) {
                var->hap_to_cons_alle[j] = -1;
            }
        }
    }
}

void var_init_hap_cons_alle0(cand_var_t *var, int min_alt_dp) {
    // select the most common allele as the consensus allele, based on hap_to_alle_profile
    // int n_depth = var->n_depth;
    for (int hap = 1; hap <= LONGCALLD_DEF_PLOID; ++hap) {
        int max_cov = 0, max_cov_alle_i = -1;
        for (int i = 0; i < var->n_uniq_alles; ++i) {
            if (var->hap_to_alle_profile[hap][i] > max_cov) {
                max_cov = var->hap_to_alle_profile[hap][i];
                if (max_cov >= min_alt_dp) max_cov_alle_i = i;
            }
        }
        // if (max_cov_alle_i == -1) {
            // if (LONGCALLD_VERBOSE >= 2) _err_func_printf("No HAP allele %d in Var: %ld-%d-%c\n", hap, var->pos, var->ref_len, BAM_CIGAR_STR[var->var_type]);
            // var->is_skipped = 1;
        // } else {
        var->hap_to_cons_alle[hap] = max_cov_alle_i;
        var->is_skipped = 0;
        // }
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

hts_pos_t assign_var_hap_based_on_pre_reads(cand_var_t *var, int min_dp) {
    if (var_hap_profile_cov(var) < min_dp) 
        return assign_var_init_hap(var);
    else {
        return assign_var_hap_based_on_pre_reads1(var);
    }
}

// after a read is assigned with hap, update hap of all other SNPs covered by this read
// including homozygous and heterozygous SNPs
void update_var_hap_profile_based_on_aln_hap(int hap, hts_pos_t phase_set, cand_var_t *var, int *var_i_to_cate, int target_var_cate, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
        var[var_i].hap_to_alle_profile[hap][allele_i] += 1;
        if (var[var_i].phase_set == 0 || phase_set <= var[var_i].pos)
            var[var_i].phase_set = phase_set;
    }
}

// if the first read' hap is 2, then flip all haps
void flip_aln_hap(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->haps[i] == 1 || chunk->haps[i] == 2) {
            if (chunk->haps[i] == 1) return;
            else break;
        } else continue;
    }
    for (int i = 0; i < chunk->n_reads; ++i) {
        if (chunk->haps[i] == 1) chunk->haps[i] = 2;
        else if (chunk->haps[i] == 2) chunk->haps[i] = 1;
    }
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
int update_var_aln_hap1(int target_read_i, int cur_hap,  bam_chunk_t *chunk, read_var_profile_t *p, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate) {
    int start_var_idx = p[target_read_i].start_var_idx, end_var_idx = p[target_read_i].end_var_idx;
    // deduct target read from hap_to_alle_profile, then compare target read's var profile with hap_cons_alle
    int *hap_match_cnt = (int*)calloc((LONGCALLD_DEF_PLOID+1), sizeof(int));
    int *tmp_hap_to_cons_alle = (int*)malloc((LONGCALLD_DEF_PLOID+1) * sizeof(int));

    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        if (p[target_read_i].var_is_used[read_var_idx] == 0) continue;

        cand_var_t *var = cand_vars+var_i;
        if (var->is_skipped) continue;
        int allele_i = p[target_read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
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
            _err_func_printf("Read %s max_cnt == 0 (pos: %" PRId64 ")\n", bam_get_qname(chunk->reads[target_read_i]), chunk->reads[target_read_i]->core.pos);
        return 0; // unknown
    } else if (max_cnt == sec_cnt) {
        if (LONGCALLD_VERBOSE >= 2)
            _err_func_printf("Read %s max_cnt == sec_cnt (%" PRId64 ")\n", bam_get_qname(chunk->reads[target_read_i]), chunk->reads[target_read_i]->core.pos);
        max_hap = cur_hap;
    }
    return max_hap;
}

int update_var_hap_profile_based_on_changed_hap(int new_hap, int old_hap, cand_var_t *cand_vars, int *var_i_to_cate, int target_var_cate, int min_alt_dp, read_var_profile_t *p, int read_i) {
    int start_var_idx = p[read_i].start_var_idx, end_var_idx = p[read_i].end_var_idx;
    for (int var_i = start_var_idx; var_i <= end_var_idx; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0) continue;
        int read_var_idx = var_i - start_var_idx;
        if (p[read_i].var_is_used[read_var_idx] == 0) continue;
        int allele_i = p[read_i].alleles[read_var_idx];
        if (allele_i == -1) continue;
        cand_var_t *var = cand_vars+var_i;
        if (LONGCALLD_VERBOSE >= 2)
            fprintf(stderr, "pos: %" PRId64 ", old_hap: %d, new_hap: %d, var: %d\n", var->pos, old_hap, new_hap, allele_i);
        var->hap_to_alle_profile[old_hap][allele_i] -= 1;
        var->hap_to_alle_profile[new_hap][allele_i] += 1;
        var_init_hap_cons_alle0(var, min_alt_dp);
    }
    return 0;
}

// XXX hap_cons_base not upated yet
// start with the very first SNP
// 1st round: assign hap to SNP's base, then assign hap to reads covering this SNP,
// 2~N rounds: re-assign haplotypes to reads based on haplotype clusters in previous rounds, 
//             until no changes to any reads
// output chunk->haps[i] to 1 or 2, 0: unknown
char read_name[1024] = "m84039_231005_222902_s1/80479720/ccs";
int assign_hap_based_on_het_vars(bam_chunk_t *chunk, int target_var_cate, call_var_opt_t *opt) {
    read_var_profile_t *p = chunk->read_var_profile;
    // int n_cand_vars = chunk->var_cate_counts[target_var_cate];
    int n_cand_vars = chunk->n_cand_vars;
    // int *var_cate_idx = chunk->var_cate_idx[target_var_cate];
    int *var_i_to_cate = chunk->var_i_to_cate;
    cand_var_t *cand_vars = chunk->cand_vars;
    cgranges_t *read_var_cr = chunk->read_var_cr;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;

    // 1st loop: var-wise loop
    var_init_hap_profile(cand_vars, n_cand_vars, var_i_to_cate, target_var_cate); // var_cate_idx);
    // for (int i = 0; i < n_cand_vars; ++i) {
        // int var_i = var_cate_idx[i];
    for (int var_i = 0; var_i < n_cand_vars; ++var_i) {
        if ((var_i_to_cate[var_i] & target_var_cate) == 0 || var_i_to_cate[var_i] == LONGCALLD_CAND_HOM_VAR) // skip cand homozygous SNPs
            continue;
        cand_var_t *var = cand_vars+var_i;
        ovlp_n = cr_overlap(read_var_cr, "cr", var_i, var_i+1, &ovlp_b, &max_b);
        hts_pos_t phase_set = assign_var_hap_based_on_pre_reads(var, opt->min_dp); // update alle_to_hap
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int read_i = cr_label(read_var_cr, ovlp_b[ovlp_i]);
            int read_var_idx = var_i - p[read_i].start_var_idx;
            if (chunk->haps[read_i] == 0 && chunk->is_skipped[read_i] != 1 && p[read_i].var_is_used[read_var_idx] == 1) {
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
                    update_var_hap_profile_based_on_aln_hap(hap, phase_set, cand_vars, var_i_to_cate, target_var_cate, p, read_i);
                    // update PS for the read
                    if (chunk->PS[read_i] == 0 || chunk->PS[read_i] > phase_set) chunk->PS[read_i] = phase_set;
                }
            }
        }
        var_init_hap_cons_alle0(var, opt->min_alt_dp); // update hap_to_cons_alle
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
        for (int read_i = 0; read_i < chunk->n_reads; ++read_i) {
            if (chunk->is_skipped[read_i] || chunk->haps[read_i] == 0) continue;
            // if (strcmp(read_name, bam_get_qname(chunk->reads[read_i])) == 0) 
                // fprintf(stderr, "ok\n");
            int cur_hap = chunk->haps[read_i];
            // XXX TODO: potential local optima
            int new_hap = update_var_aln_hap1(read_i, cur_hap, chunk, p, cand_vars, var_i_to_cate, target_var_cate);
            if (new_hap != cur_hap) { // update bam_haps, hap_to_base_profile, hap_to_cons_base
                if (LONGCALLD_VERBOSE >= 2) {
                    fprintf(stderr, "read (%d): %s, pos: %" PRId64 "\t", read_i, bam_get_qname(chunk->reads[read_i]), chunk->reads[read_i]->core.pos);
                    fprintf(stderr, "\t\t cur_hap: %d, new_hap: %d\n", cur_hap, new_hap);
                }
                changed_hap = 1;
                chunk->haps[read_i] = new_hap; // update intermediately
                // bam_aux_append(chunk->reads[read_i], "XT", 'i', 4, (uint8_t*)&(chunk->haps[read_i]));
                update_var_hap_profile_based_on_changed_hap(new_hap, cur_hap, cand_vars, var_i_to_cate, target_var_cate, opt->min_alt_dp, p, read_i);
            }
        } if (changed_hap == 0) break;
    }
    if (LONGCALLD_VERBOSE >= 2) _err_info("Iteration: %d\n", i_iter);
    free(ovlp_b); 
    // cr_destroy(read_var_cr);
    return 0;
}