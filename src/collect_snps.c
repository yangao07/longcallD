#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "collect_snps.h"
#include "utils.h"
#include "bam_utils.h"
#include "seq.h"
#include "assign_aln_hap.h"

// init cand_snps_t
// cand_snps_t *init_cand_snps(void) {
//     cand_snps_t *snps_t = _err_malloc(sizeof(cand_snps_t));
//     int m_snps = 64, n_reads = 512;
//     // snps
//     snps_t->n_snps = 0; snps_t->m_snps = m_snps;
//     snps_t->snp_pos = _err_malloc(snps_t->m_snps * sizeof(int));
//     snps_t->n_snp_bases = _err_malloc(snps_t->m_snps * sizeof(int));
//     snps_t->snp_bases = _err_malloc(snps_t->m_snps * sizeof(uint8_t*));
//     snps_t->used_for_phasing_snps = _err_malloc(snps_t->m_snps * sizeof(uint8_t));
//     // reads 
//     snps_t->n_reads = n_reads;
//     snps_t->read_to_start_snp = _err_malloc(n_reads * sizeof(int));
//     snps_t->read_to_end_snp = _err_malloc(n_reads * sizeof(int));
//     // reads X snps
//     snps_t->read_to_snp_map = _err_malloc(n_reads * sizeof(uint8_t*));
//     for (int i = 0; i < n_reads; i++) {
//         snps_t->read_to_snp_map[i] = _err_malloc(m_snps * sizeof(uint8_t));
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
    cand_snp_t *mis_sites = _err_malloc(n_mis_pos * sizeof(cand_snp_t));
    for (int i = 0; i < n_mis_pos; ++i) {
        mis_sites[i].pos = mis_pos[i];
        mis_sites[i].n_depth = 0; mis_sites[i].ref_base_cov = 0; mis_sites[i].n_uniq_alt_bases = 0;
        mis_sites[i].alt_bases = _err_malloc(5 * sizeof(uint8_t)); mis_sites[i].alt_base_covs = _err_malloc(5 * sizeof(int));
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

// void copy_mis(cand_snp_t *mis_sites, int mis_i, mis_base_t mis_base) {
//     mis_sites[mis_i].pos = mis_base.pos;
//     mis_sites[mis_i].n_depth = 1; mis_sites[mis_i].n_uniq_bases = 1;
//     mis_sites[mis_i].bases[0] = mis_base.base; mis_sites[mis_i].base_covs[0] = 1;
// }

// void merge_mist1(cand_snp_t *mis_sites, int mis_i, cand_snp_t snp1, mis_base_t mis_base2) {
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

void free_mis_sites(cand_snp_t *mis_sites, int m) {
    if (m > 0) {
        for (int i = 0; i < m; ++i) {
            free(mis_sites[i].alt_bases); free(mis_sites[i].alt_base_covs);
        }
        free(mis_sites);
    }
}

int merge_mis_pos(int n_total_mis_sites, hts_pos_t **mis_pos, int n_mis_bases, mis_base_t *mis_bases) {
    int new_total_mis_sites = 0;
    hts_pos_t *new_mis_pos = _err_malloc((n_total_mis_sites + n_mis_bases) * sizeof(hts_pos_t));
    int i, j;
    for (i = j =0; i < n_total_mis_sites && j < n_mis_bases; ) {
        if ((*mis_pos)[i] < mis_bases[j].pos) {
            new_mis_pos[new_total_mis_sites] = (*mis_pos)[i];
            i++;
        } else if ((*mis_pos)[i] > mis_bases[j].pos) {
            new_mis_pos[new_total_mis_sites] = mis_bases[j].pos;
            j++;
        } else {
            new_mis_pos[new_total_mis_sites] = (*mis_pos)[i];
            i++; j++;
        }
        new_total_mis_sites++;
    }
    for (; i < n_total_mis_sites; ++i, ++new_total_mis_sites) new_mis_pos[new_total_mis_sites] = (*mis_pos)[i];
    for (; j < n_mis_bases; ++j, ++new_total_mis_sites) new_mis_pos[new_total_mis_sites] = mis_bases[j].pos;
    free(*mis_pos); *mis_pos = new_mis_pos;
    return new_total_mis_sites;
}

// // read_id: read index in bam_chunk
// int merge_mis_sites(int read_id, int *n_total_mis_sites, int *m_total_cand_snps, cand_snp_t **mis_sites, 
//                      int n_mis_bases, mis_base_t *mis_bases) {
//     int new_total_mis_sites = 0;
//     cand_snp_t *new_mis_sites = _err_malloc((*n_total_mis_sites + n_mis_bases) * sizeof(cand_snp_t));
//     for (int i = 0; i < *n_total_mis_sites + n_mis_bases; ++i) {
//         new_mis_sites[i].bases = _err_malloc(5 * sizeof(uint8_t));
//         new_mis_sites[i].base_covs = _err_malloc(5 * sizeof(int));
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

int *collect_cand_snp_sites(int n_total_mis_sites, cand_snp_t *mis_sites) {
    return 0;
}

// input: sorted sites of mismatch bases
// working: update mis_sites, calculate the depth and allele frequency of each site
//          considering both reference allele and alternative alleles
// output: updated mis_sites
int update_cand_snps(cand_snp_t *and_snps, int n_total_mis_pos, hts_pos_t *mis_pos, const phase_bam_pl_t *pl, const phase_bam_opt_t *opt, bam_chunk_t *bam_chunk) {
    // 2nd pass: update mis_sites, calculate the depth and allele frequency of each site
    int start_mis_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        bam1_t *read = bam_chunk->reads[i];
        // check if MQ >= min_mq
        if (read->core.qual < opt->min_mq) continue;
        // printf("Read name: %s, start_mis_i: %d\n", bam_get_qname(read), start_mis_i);
        if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            start_mis_i += update_mis_sites_from_eqx_cigar(read, opt->min_bq, n_total_mis_pos-start_mis_i, mis_pos+start_mis_i, and_snps+start_mis_i);
        } 
        // else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
        //     start_mis_i += update_mis_sites_from_MD_tag(read, opt->min_bq, n_total_mis_pos, mis_pos, mis_sites);
        // } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
        //     if (pl->ref_seq == NULL) _err_fatal("Reference file required.");
        //     char *chrom = pl->header->target_name[read->core.tid];
        //     start_mis_i += update_mis_sites_from_ref_seq(read, opt->min_bq, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, chrom), n_total_mis_pos, mis_pos, mis_sites);
        // }
    }
    // for (int i = 0; i < n_total_mis_pos; ++i) {
    //     printf("Mis: %lld\t", (long long) mis_pos[i]);
    //     printf("Depth: %d\t", mis_sites[i].n_depth);
    //     printf("Ref: %d\tAlt: ", mis_sites[i].ref_base_cov);
    //     for (int j = 0; j < mis_sites[i].n_uniq_alt_bases; ++j) {
    //         printf("%c: %d\t", "ACGTN"[mis_sites[i].alt_bases[j]], mis_sites[i].alt_base_covs[j]);
    //     } printf("\n");
    // }
    return 0;
}

// given a CHUNK of reads, collect the candidate SNP sites
int collect_mis_sites(hts_pos_t **mis_pos, const phase_bam_pl_t *pl, const phase_bam_opt_t *opt, bam_chunk_t *bam_chunk) {
    // 1st pass: collect candidate SNP sites, only check mismatch bases
    int n_total_mis_sites = 0;
    // cand_snp_t *cand_sites = _err_malloc(m_total_cand_sites * sizeof(cand_snp_t));
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        bam1_t *read = bam_chunk->reads[i];
        // check if MQ >= min_mq
        if (read->core.qual < opt->min_mq) continue;
        int n_mis_bases = 0; mis_base_t *mis_bases = NULL;
        if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            n_mis_bases = get_mis_bases_from_eqx_cigar(read, opt->min_bq, &mis_bases);
        } else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
            n_mis_bases = get_mis_bases_from_MD_tag(read, opt->min_bq, &mis_bases);
        } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
            if (pl->ref_seq == NULL) _err_fatal("Reference file required.");
            char *chrom = pl->header->target_name[read->core.tid];
            n_mis_bases = get_mis_bases_from_ref_seq(read, opt->min_bq, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, chrom), &mis_bases);
        }
        n_total_mis_sites = merge_mis_pos(n_total_mis_sites, mis_pos, n_mis_bases, mis_bases);
        if (n_mis_bases > 0) free(mis_bases);
    }
    // print cand_snps
    // printf("Total candidate SNP sites: %d\n", n_total_mis_sites);
    // for (int i = 0; i < n_total_mis_sites; ++i) {
    //     printf("Mis: %lld\n", (long long) (*mis_pos)[i]);
    // }
    return n_total_mis_sites;
}

char test_read_name[1024] = "m84039_231005_222902_s1/80349975/ccs";
read_snp_profile_t *collect_read_snp_profile(int n_cand_snps, cand_snp_t *cand_snp_sites, const phase_bam_pl_t *pl, const phase_bam_opt_t *opt, bam_chunk_t *bam_chunk) {
    read_snp_profile_t *p = init_read_snp_profile(bam_chunk->n_reads, n_cand_snps);
    // 3rd pass: collect read-wise SNP profiles
    int start_snp_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        bam1_t *read = bam_chunk->reads[i];
        // check if MQ >= min_mq
        if (read->core.qual < opt->min_mq) continue;
        if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            printf("Read: %s\n", bam_get_qname(read));
        printf("Read: %s\n", bam_get_qname(read));
        if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            start_snp_i = update_read_snp_profile_from_eqx_cigar(read, n_cand_snps, cand_snp_sites, start_snp_i, p+i);
            // printf("start_snp_i: %d\n", start_snp_i);
        }
        if (p[i].start_snp_idx >= 0) {
            // printf("Read: %s, start_snp_i: %d, end_snp_i: %d\n", bam_get_qname(read), p[i].start_snp_idx, p[i].end_snp_idx);
            for (int j = p[i].start_snp_idx; j <= p[i].end_snp_idx; ++j) {
                printf("\tSNP: (%d) %lld, base: %c, qual: %d\n", j, (long long) cand_snp_sites[j].pos, 
                LONGCALLD_BAM_BASE_STR[p[i].snp_bases[j-p[i].start_snp_idx]], p[i].snp_qual[j-p[i].start_snp_idx]);
            }
        }
    }
    return p;
}

int filter_cand_snps(cand_snp_t *mis_sites, int n_mis_pos, const phase_bam_opt_t *opt) {
    // no need to malloc for cand_snp_pos
    int cand_snp_i = 0;
    int min_dp = opt->min_dp; double min_af = opt->min_af, max_af = opt->max_af;
    for (int i = 0; i < n_mis_pos; ++i) {
        if (mis_sites[i].n_depth < min_dp) continue;
        int ref_base_cov = mis_sites[i].ref_base_cov;
        double _af = (double) ref_base_cov / mis_sites[i].n_depth;
        for (int j = 0; j < mis_sites[i].n_uniq_alt_bases; ++j) {
            _af = MAX_OF_TWO(_af, (double) mis_sites[i].alt_base_covs[j] / mis_sites[i].n_depth);
        }
        if (_af < min_af || _af > max_af) continue; // XXX
        printf("Cand SNP: %lld, n_depth: %d, af: %f, Ref: %d\t", (long long) mis_sites[i].pos, mis_sites[i].n_depth, _af, ref_base_cov);
        for (int j = 0; j < mis_sites[i].n_uniq_alt_bases; ++j) {
            printf("%c: %d\t", LONGCALLD_BAM_BASE_STR[mis_sites[i].alt_bases[j]], mis_sites[i].alt_base_covs[j]);
        } printf("\n");
        // copy i'th to cand_snp_i'th
        if (i != cand_snp_i) {
            mis_sites[cand_snp_i].pos = mis_sites[i].pos;
            mis_sites[cand_snp_i].n_depth = mis_sites[i].n_depth;
            mis_sites[cand_snp_i].ref_base_cov = mis_sites[i].ref_base_cov;
            mis_sites[cand_snp_i].n_uniq_alt_bases = mis_sites[i].n_uniq_alt_bases;
            for (int j = 0; j < mis_sites[i].n_uniq_alt_bases; ++j) {
                mis_sites[cand_snp_i].alt_bases[j] = mis_sites[i].alt_bases[j];
                mis_sites[cand_snp_i].alt_base_covs[j] = mis_sites[i].alt_base_covs[j];
            }
        }
        cand_snp_i++;
        // cand_snp_pos[cand_snp_i++] = mis_sites[i].pos;
    }
    for (int i = cand_snp_i; i < n_mis_pos; ++i) {
        free(mis_sites[i].alt_bases); free(mis_sites[i].alt_base_covs);
    }
    return cand_snp_i;
}

int collect_snps_main(const phase_bam_pl_t *pl, bam_chunk_t *bam_chunk) {
    hts_pos_t *mis_pos = NULL;
    // collect mismatch positions
    int n_mis_pos = collect_mis_sites(&mis_pos, pl, pl->pb_opt, bam_chunk);
    if (n_mis_pos <= 0) return 0;

    cand_snp_t *cand_snps = init_cand_snps(n_mis_pos, mis_pos);
    // collect reference and alternative alleles for all mismatch positions
    update_cand_snps(cand_snps, n_mis_pos, mis_pos, pl, pl->pb_opt, bam_chunk);

    // filter candidate SNPs based depth, allele frequency, etc.
    int n_cand_snps = filter_cand_snps(cand_snps, n_mis_pos, pl->pb_opt);

    // collect read-wise snp profiles
    read_snp_profile_t *p = collect_read_snp_profile(n_cand_snps, cand_snps, pl, pl->pb_opt, bam_chunk);

    // assign_aln_hap(p, n_cand_snps, cand_snps, bam_chunk);
    if (n_mis_pos > 0) free(mis_pos);
    free_read_snp_profile(p, bam_chunk->n_reads); free_mis_sites(cand_snps, n_cand_snps); 
    return 0;
}
