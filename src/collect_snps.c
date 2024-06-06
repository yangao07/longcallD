#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "collect_snps.h"
#include "utils.h"
#include "bam_utils.h"
#include "seq.h"
#include "assign_hap.h"
#include "vcf_utils.h"

extern int LONGCALLD_VERBOSE;
char test_read_name[1024] = "m84039_230928_213653_s3/45487002/ccs";
// init cand_snps_t
// cand_snps_t *init_cand_snps(void) {
//     cand_snps_t *snps_t = malloc(sizeof(cand_snps_t));
//     int m_snps = 64, n_reads = 512;
//     // snps
//     snps_t->n_snps = 0; snps_t->m_snps = m_snps;
//     snps_t->snp_pos = malloc(snps_t->m_snps * sizeof(int));
//     snps_t->n_snp_bases = malloc(snps_t->m_snps * sizeof(int));
//     snps_t->snp_bases = malloc(snps_t->m_snps * sizeof(uint8_t*));
//     snps_t->used_for_phasing_snps = malloc(snps_t->m_snps * sizeof(uint8_t));
//     // reads 
//     snps_t->n_reads = n_reads;
//     snps_t->read_to_start_snp = malloc(n_reads * sizeof(int));
//     snps_t->read_to_end_snp = malloc(n_reads * sizeof(int));
//     // reads X snps
//     snps_t->read_to_snp_map = malloc(n_reads * sizeof(uint8_t*));
//     for (int i = 0; i < n_reads; i++) {
//         snps_t->read_to_snp_map[i] = malloc(m_snps * sizeof(uint8_t));
//     }
//     return snps_t;
// }

// // free cand_snps_t
// void free_cand_snps(cand_snps_t *snps_t) {
//     for (int i = 0; i < snps_t->n_reads; i++) {
//         free(snps_t->read_to_snp_map[i]);
//     }
//     free(snps_t->read_to_snp_map);
//     free(snps_t->read_to_end_snp);
//     free(snps_t->read_to_start_snp);
//     free(snps_t->used_for_phasing_snps);
//     for (int i = 0; i < snps_t->n_snps; i++) {
//         free(snps_t->snp_bases[i]);
//     }
//     free(snps_t->snp_bases);
//     free(snps_t->n_snp_bases);
//     free(snps_t);
// }

// init cand_snp_t
cand_snp_t *init_cand_snps(int n_mis_pos, hts_pos_t *mis_pos) {
    cand_snp_t *mis_sites = (cand_snp_t*)malloc(n_mis_pos * sizeof(cand_snp_t));
    for (int i = 0; i < n_mis_pos; ++i) {
        // static information
        mis_sites[i].pos = mis_pos[i];
        mis_sites[i].n_depth = 0; mis_sites[i].n_uniq_bases = 0;
        mis_sites[i].bases = (uint8_t*)malloc(LONGCALLD_BAM_BASE_N * sizeof(uint8_t)); 
        mis_sites[i].base_covs = (int*)malloc(LONGCALLD_BAM_BASE_N * sizeof(int));
        mis_sites[i].base_to_i = (int*)malloc(LONGCALLD_BAM_BASE_N * sizeof(int));
        // dynamic information, allocate and update during haplotype assignment
        mis_sites[i].base_to_hap = NULL; mis_sites[i].hap_to_base_profile = NULL; mis_sites[i].hap_to_cons_base = NULL;
    }
    return mis_sites;
}

// void copy_mis_t(cand_snp_t *mis_sites, int mis_i, cand_snp_t mis) {
//     mis_sites[mis_i].pos = mis.pos;
//     mis_sites[mis_i].n_depth = mis.n_depth; mis_sites[mis_i].n_uniq_bases = mis.n_uniq_bases;
//     for (int i = 0; i < mis.n_uniq_bases; ++i) {
//         mis_sites[mis_i].bases[i] = mis.bases[i];
//         mis_sites[mis_i].base_covs[i] = mis.base_covs[i];
//     }
// }

// void copy_mis(cand_snp_t *mis_sites, int mis_i, x_t mis_base) {
//     mis_sites[mis_i].pos = mis_base.pos;
//     mis_sites[mis_i].n_depth = 1; mis_sites[mis_i].n_uniq_bases = 1;
//     mis_sites[mis_i].bases[0] = mis_base.base; mis_sites[mis_i].base_covs[0] = 1;
// }

// void merge_mist1(cand_snp_t *mis_sites, int mis_i, cand_snp_t snp1, x_t mis_base2) {
//     mis_sites[mis_i].pos = mis_base2.pos;
//     mis_sites[mis_i].n_depth = snp1.n_depth + 1;

//     int base_i = snp1.n_uniq_bases, exist=0;
//     mis_sites[mis_i].n_uniq_bases = snp1.n_uniq_bases;
//     for (int i = 0; i < snp1.n_uniq_bases; ++i) {
//         if (snp1.bases[i] == mis_base2.base) {
//             exist = 1; base_i = i;
//         }
//         mis_sites[mis_i].bases[i] = snp1.bases[i];
//         mis_sites[mis_i].base_covs[i] = snp1.base_covs[i];
//     }
//     if (exist == 0) {
//         mis_sites[mis_i].bases[base_i] = mis_base2.base;
//         mis_sites[mis_i].base_covs[base_i] = 1;
//         mis_sites[mis_i].n_uniq_bases += 1;
//     } else {
//         mis_sites[mis_i].base_covs[base_i] += 1;
//     }
// }

void free_cand_snps(cand_snp_t *mis_sites, int m) {
    if (m > 0) {
        for (int i = 0; i < m; ++i) {
            free(mis_sites[i].bases); free(mis_sites[i].base_covs); free(mis_sites[i].base_to_i);
            if (mis_sites[i].base_to_hap != NULL) free(mis_sites[i].base_to_hap);
            if (mis_sites[i].hap_to_base_profile != NULL) {
                for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
                    free(mis_sites[i].hap_to_base_profile[j]);
                } free(mis_sites[i].hap_to_base_profile);
            }
            if (mis_sites[i].hap_to_cons_base != NULL) free(mis_sites[i].hap_to_cons_base);
        }
        free(mis_sites);
    }
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

// // read_id: read index in bam_chunk
// int merge_mis_sites(int read_id, int *n_total_mis_sites, int *m_total_cand_snps, cand_snp_t **mis_sites, 
//                      int n_mis_bases, x_t *mis_bases) {
//     int new_total_mis_sites = 0;
//     cand_snp_t *new_mis_sites = malloc((*n_total_mis_sites + n_mis_bases) * sizeof(cand_snp_t));
//     for (int i = 0; i < *n_total_mis_sites + n_mis_bases; ++i) {
//         new_mis_sites[i].bases = malloc(5 * sizeof(uint8_t));
//         new_mis_sites[i].base_covs = malloc(5 * sizeof(int));
//     }
//     // merge sort cand_snp_pos and mis_pos, store new positions in new_cand_snp_pos
//     int i, j;
//     for (i = j = 0; i < *n_total_mis_sites && j < n_mis_bases; ) {
//         if ((*mis_sites)[i].pos < mis_bases[j].pos) {
//             copy_mis_t(new_mis_sites, new_total_mis_sites, (*mis_sites)[i]);
//             i++;
//         } else if ((*mis_sites)[i].pos > mis_bases[j].pos) {
//             copy_mis(new_mis_sites, new_total_mis_sites, mis_bases[j]);
//             j++;
//         } else { // merge
//             merge_mist1(new_mis_sites, new_total_mis_sites, (*mis_sites)[i], mis_bases[j]);
//             i++; j++;
//         }
//         new_total_mis_sites++;
//     }
//     for (; i < *n_total_mis_sites; ++i, ++new_total_mis_sites) copy_mis_t(new_mis_sites, new_total_mis_sites, (*mis_sites)[i]);
//     for (; j < n_mis_bases; ++j, ++new_total_mis_sites) copy_mis(new_mis_sites, new_total_mis_sites, mis_bases[j]);
//     // copy new_cand_snp_pos to cand_snp_pos
//     mis_sites_free(*mis_sites, *m_total_cand_snps); *mis_sites = new_mis_sites; 
//     *m_total_cand_snps = *n_total_mis_sites + n_mis_bases; *n_total_mis_sites = new_total_mis_sites;
//     return 0;
// }

// cand_snp_t *collect_cand_snps(int n_x_pos, hts_pos_t *x_pos, bam_chunk_t *bam_chunk) {
//     cand_snp_t *cand_snps = init_cand_snps(n_x_pos, x_pos);
//     int start_x_i = 0;
//     for (int read_i = 0; read_i < bam_chunk->n_reads; ++read_i) {
//         if (bam_chunk->is_skipped[read_i]) continue;
//         bam1_t *read = bam_chunk->reads[read_i];
//         start_x_i += update_cand_snps_from_readx(bam_chunk->read_xs+read_i, n_x_pos-start_x_i, x_pos+start_x_i, cand_snps+start_x_i);
//     }
//     return cand_snps;
// }

cand_snp_t *collect_cand_snps(bam_chunk_t *bam_chunk, int n_x_sites, hts_pos_t *x_sites) {
    cand_snp_t *cand_snps = init_cand_snps(n_x_sites, x_sites);
    // 2nd pass: update mis_sites, calculate the depth and allele frequency of each site
    int start_x_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        start_x_i = update_cand_snps_from_digar(bam_chunk->digars+i, n_x_sites, x_sites, start_x_i, cand_snps);
        // start_x_i += update_x_sites_from_eqx_cigar(read, opt->min_bq, n_x_sites-start_x_i, x_sites+start_x_i, cand_snps+start_x_i);
        // else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
        //     start_mis_i += update_mis_sites_from_MD_tag(read, opt->min_bq, n_total_mis_pos, mis_pos, mis_sites);
        // } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
        //     if (pl->ref_seq == NULL) _err_fatal("Reference file required.");
        //     char *chrom = pl->header->target_name[read->core.tid];
        //     start_mis_i += update_mis_sites_from_ref_seq(read, opt->min_bq, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, chrom), n_total_mis_pos, mis_pos, mis_sites);
        // }
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_x_sites; ++i) {
            fprintf(stderr, "Mis: %lld\t", (long long) x_sites[i]);
            fprintf(stderr, "Depth: %d\t", cand_snps[i].n_depth);
            for (int j = 0; j < cand_snps[i].n_uniq_bases; ++j) {
                fprintf(stderr, "%c: %d\t", LONGCALLD_BAM_BASE_STR[cand_snps[i].bases[j]], cand_snps[i].base_covs[j]);
            } fprintf(stderr, "\n");
        }
    }
    return cand_snps;
}

// given a CHUNK of reads, collect all the valid x sites
// int collect_x_sites(hts_pos_t **x_sites, const call_var_pl_t *pl, bam_chunk_t *bam_chunk) {
//     // 1st pass: collect candidate SNP sites, only check mismatch bases
//     const call_var_opt_t *opt = pl->opt;
//     int n_total_x_sites = 0;
//     // cand_snp_t *cand_sites = malloc(m_total_cand_sites * sizeof(cand_snp_t));
//     for (int i = 0; i < bam_chunk->n_reads; ++i) {
//         if (bam_chunk->is_skipped[i]) continue;
//         int read_n_x_sites = 0; x_base_t *read_x_sites = NULL;
//         // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
//             // printf("ok\n");
//         if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
//             read_n_x_sites = get_x_from_eqx_cigar(read, opt, &read_x_sites);
//         } 
//         // else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
//         //     read_n_x_sites = get_mis_bases_from_MD_tag(read, opt->min_bq, &read_x_sites);
//         // } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
//         //     if (pl->ref_seq == NULL) _err_fatal("Reference file required.");
//         //     char *chrom = pl->header->target_name[read->core.tid];
//         //     read_n_x_sites = get_mis_bases_from_ref_seq(read, opt->min_bq, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, chrom), &read_x_sites);
//         // }
//         // printf("Read: %s, n_mis_bases: %d\n", bam_get_qname(read), n_mis_bases);
//         // for (int j = 0; j < n_mis_bases; ++j) {
//             // printf("Mis: %lld, base: %c, is_dense: %c\n", (long long) mis_bases[j].pos, LONGCALLD_BAM_BASE_STR[mis_bases[j].base], "NY"[mis_bases[j].is_low_qual]);
//         // }
//         // bam_chunk->haps[i] = n_mis_bases;
//         n_total_x_sites = merge_x_sites(n_total_x_sites, x_sites, read_n_x_sites, read_x_sites);
//         if (read_n_x_sites > 0) free(read_x_sites);
//     }
//     // print cand_snps
//     fprintf(stderr, "Total candidate SNP sites: %d\n", n_total_x_sites);
//     for (int i = 0; i < n_total_x_sites; ++i) {
//         fprintf(stderr, "Mis: %lld\n", (long long) (*x_sites)[i]);
//     }
//     return n_total_x_sites;
// }

// read_snp_profile_t *collect_read_snp_profile(int n_cand_snps, cand_snp_t *cand_snp_sites, const call_var_pl_t *pl, const call_var_opt_t *opt, bam_chunk_t *bam_chunk) {
//     read_snp_profile_t *p = init_read_snp_profile(bam_chunk->n_reads, n_cand_snps);
//     // 3rd pass: collect read-wise SNP profiles
//     int start_snp_i = 0;
//     for (int i = 0; i < bam_chunk->n_reads; ++i) {
//         if (bam_chunk->is_skipped[i]) continue;
//         bam1_t *read = bam_chunk->reads[i];
//         // check if MQ >= min_mq
//         if (read->core.qual < opt->min_mq) continue;
//         // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
//             // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
//         // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
//         if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
//             start_snp_i = update_read_snp_profile_from_eqx_cigar(read, n_cand_snps, cand_snp_sites, start_snp_i, p+i);
//             // fprintf(stderr, "start_snp_i: %d\n", start_snp_i);
//         }
//         // if (p[i].start_snp_idx >= 0) {
//         //     // fprintf(stderr, "Read: %s, start_snp_i: %d, end_snp_i: %d\n", bam_get_qname(read), p[i].start_snp_idx, p[i].end_snp_idx);
//         //     for (int j = p[i].start_snp_idx; j <= p[i].end_snp_idx; ++j) {
//         //         fprintf(stderr, "\tSNP: (%d) %lld, base: %c, qual: %d\n", j, (long long) cand_snp_sites[j].pos, 
//         //         LONGCALLD_BAM_BASE_STR[p[i].snp_bases[j-p[i].start_snp_idx]], p[i].snp_qual[j-p[i].start_snp_idx]);
//         //     }
//         // }
//     }
//     return p;
// }

// XXX filter: usable X/= over total non-low-qual depth >= threshold
int filter_cand_snps(cand_snp_t *cand_snps, int n_x_sites, const call_var_opt_t *opt) {
    // no need to malloc for cand_snp_pos
    int cand_snp_i = 0;
    int min_dp = opt->min_dp; double min_af = opt->min_af, max_af = opt->max_af;
    for (int i = 0; i < n_x_sites; ++i) {
        if (cand_snps[i].n_depth < min_dp) continue;
        double _af1 = 0.0, _af2 = 0.0, _af = 0.0;
        for (int j = 0; j < cand_snps[i].n_uniq_bases; ++j) {
            if (cand_snps[i].bases[j] == LONGCALLD_BAM_DEL_BASE_IDX) continue;
            _af = (double) cand_snps[i].base_covs[j] / cand_snps[i].n_depth;
            if (_af > _af1) { _af2 = _af1; _af1 = _af; } 
            else if (_af > _af2 && _af <= _af1) _af2 = _af;
        }
        if (_af1 < min_af || _af1 > max_af || _af2 < min_af || _af2 > max_af) continue; // XXX
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "Cand SNP %d: %lld, depth: %d, AF: %f|%f\t", cand_snp_i, (long long) cand_snps[i].pos, cand_snps[i].n_depth, _af1, _af2);
            for (int j = 0; j < cand_snps[i].n_uniq_bases; ++j) {
                fprintf(stderr, "%c: %d\t", LONGCALLD_BAM_BASE_STR[cand_snps[i].bases[j]], cand_snps[i].base_covs[j]);
            } fprintf(stderr, "\n");
        }
        // copy i'th to cand_snp_i'th
        if (i != cand_snp_i) {
            cand_snps[cand_snp_i].pos = cand_snps[i].pos;
            cand_snps[cand_snp_i].n_depth = cand_snps[i].n_depth;
            cand_snps[cand_snp_i].n_uniq_bases = cand_snps[i].n_uniq_bases;
            for (int j = 0; j < cand_snps[i].n_uniq_bases; ++j) {
                cand_snps[cand_snp_i].bases[j] = cand_snps[i].bases[j];
                cand_snps[cand_snp_i].base_covs[j] = cand_snps[i].base_covs[j];
            }
            for (int j = 0; j < LONGCALLD_BAM_BASE_N; ++j) cand_snps[cand_snp_i].base_to_i[j] = cand_snps[i].base_to_i[j];
        }
        cand_snp_i++;
        // cand_snp_pos[cand_snp_i++] = cand_snps[i].pos;
    }
    for (int i = cand_snp_i; i < n_x_sites; ++i) {
        free(cand_snps[i].bases); free(cand_snps[i].base_covs); free(cand_snps[i].base_to_i);
    }
    return cand_snp_i;
}

void collect_digars_from_bam(bam_chunk_t *bam_chunk, const call_var_pl_t *pl) {
    const call_var_opt_t *opt = pl->opt;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        if (read->core.qual < opt->min_mq) {
            bam_chunk->is_skipped[i] = 1; continue;
        }
        if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            collect_digar_from_eqx_cigar(read, opt, bam_chunk->digars+i);
        } else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
            collect_digar_from_MD_tag(read, opt, bam_chunk->digars+i);
        } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
            if (pl->ref_seq == NULL) _err_fatal("Reference genome FASTA file is required for BAM without =/X cigars and MD tags.");
            char *chrom = pl->header->target_name[read->core.tid];
            collect_digar_from_ref_seq(read, opt, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, chrom), bam_chunk->digars+i);
        }
    }
}

int collect_x_sites(bam_chunk_t *bam_chunk, hts_pos_t **x_sites) {
    int n_total_x_sites = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
       n_total_x_sites = merge_x_sites(n_total_x_sites, x_sites, bam_chunk->digars[i].n_digar, bam_chunk->digars[i].digars);
    }
    // print cand_snps
    // fprintf(stderr, "Total candidate SNP sites: %d\n", n_total_x_sites);
    // for (int i = 0; i < n_total_x_sites; ++i) {
    //     fprintf(stderr, "Mis: %lld\n", (long long) (*x_sites)[i]);
    // }
    return n_total_x_sites;
}

read_snp_profile_t *collect_read_snp_profile(bam_chunk_t *bam_chunk, int n_cand_snps, cand_snp_t *cand_snps) {
    read_snp_profile_t *p = init_read_snp_profile(bam_chunk->n_reads, n_cand_snps);
    // 3rd pass: collect read-wise SNP profiles
    int start_snp_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        start_snp_i = update_read_snp_profile_from_digar(bam_chunk->digars+i, read, n_cand_snps, cand_snps, start_snp_i, p+i);
            // fprintf(stderr, "start_snp_i: %d\n", start_snp_i);
        if (LONGCALLD_VERBOSE >= 2) {
            if (p[i].start_snp_idx >= 0) {
                fprintf(stderr, "Read: %s, start_snp_i: %d, end_snp_i: %d\n", bam_get_qname(read), p[i].start_snp_idx, p[i].end_snp_idx);
                for (int j = 0; j <= p[i].end_snp_idx-p[i].start_snp_idx; ++j) {
                    // fprintf(stderr, "\tSNP: (%d) %lld, base: %c, qual: %d\n", j, (long long) cand_snps[j].pos, 
                    // LONGCALLD_BAM_BASE_STR[p[i].snp_bases[j-p[i].start_snp_idx]], p[i].snp_qual[j-p[i].start_snp_idx]);
                    if (p[i].snp_is_used[j] == 0) continue;
                    fprintf(stderr, "\tSNP: (%d) %lld", j, (long long) cand_snps[j+p[i].start_snp_idx].pos);
                    fprintf(stderr, ", base: %c", LONGCALLD_BAM_BASE_STR[p[i].snp_bases[j]]);
                    fprintf(stderr, ", qual: %d\n", p[i].snp_qual[j]);
                }
            }
        }
    }
    return p;
}

int collect_snps_main(const call_var_pl_t *pl, bam_chunk_t *bam_chunk) {
    // collect X/I/D sites from BAM
    collect_digars_from_bam(bam_chunk, pl);
    
    // merge X sites from all reads
    hts_pos_t *x_sites = NULL; int n_x_sites;
    if ((n_x_sites = collect_x_sites(bam_chunk, &x_sites)) <= 0) return 0;
    
    // collect reference and alternative alleles for all X sites 
    cand_snp_t *cand_snps = collect_cand_snps(bam_chunk, n_x_sites, x_sites); free(x_sites);
    
    // filter candidate SNPs based depth, allele frequency, etc.
    int n_cand_snps;
    if ((n_cand_snps = filter_cand_snps(cand_snps, n_x_sites, pl->opt)) <= 0) return 0;
    
    // collect read-wise snp profiles
    read_snp_profile_t *p = collect_read_snp_profile(bam_chunk, n_cand_snps, cand_snps);

    // assign hap to each read and SNP
    assign_hap(p, n_cand_snps, cand_snps, bam_chunk); free_read_snp_profile(p, bam_chunk->n_reads); 

    // write snps to VCF
    write_snp_to_vcf(cand_snps, n_cand_snps, pl->opt->out_vcf, pl->header->target_name[bam_chunk->tid]); free_cand_snps(cand_snps, n_cand_snps); 
    return 0;
}