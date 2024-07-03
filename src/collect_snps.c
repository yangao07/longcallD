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

// init cand_snp_t
cand_snp_t *init_cand_snps(int n_mis_pos, hts_pos_t *mis_pos) {
    cand_snp_t *cand_snps = (cand_snp_t*)malloc(n_mis_pos * sizeof(cand_snp_t));
    for (int i = 0; i < n_mis_pos; ++i) {
        // static information
        cand_snps[i].ref_base = -1; // unset
        cand_snps[i].pos = mis_pos[i]; cand_snps[i].phase_set = 0;
        cand_snps[i].n_depth = 0; cand_snps[i].n_low_depth = 0; cand_snps[i].n_uniq_bases = 0;
        cand_snps[i].bases = (uint8_t*)malloc(LONGCALLD_BAM_BASE_N * sizeof(uint8_t)); 
        cand_snps[i].base_covs = (int*)malloc(LONGCALLD_BAM_BASE_N * sizeof(int));
        cand_snps[i].base_to_i = (int*)malloc(LONGCALLD_BAM_BASE_N * sizeof(int));

        // dynamic information, allocate and update during haplotype assignment
        cand_snps[i].base_to_hap = NULL; cand_snps[i].hap_to_base_profile = NULL; cand_snps[i].hap_to_cons_base = NULL;
        cand_snps[i].is_low_qual = 0; cand_snps[i].is_skipped = 0;
    }
    return cand_snps;
}

void free_cand_snps(cand_snp_t *snp_sites, int m) {
    if (m > 0) {
        for (int i = 0; i < m; ++i) {
            free(snp_sites[i].bases); free(snp_sites[i].base_covs); free(snp_sites[i].base_to_i);
            if (snp_sites[i].base_to_hap != NULL) free(snp_sites[i].base_to_hap);
            if (snp_sites[i].hap_to_base_profile != NULL) {
                for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
                    free(snp_sites[i].hap_to_base_profile[j]);
                } free(snp_sites[i].hap_to_base_profile);
            }
            if (snp_sites[i].hap_to_cons_base != NULL) free(snp_sites[i].hap_to_cons_base);
        }
        free(snp_sites);
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

cand_snp_t *collect_cand_snps(bam_chunk_t *bam_chunk, int n_x_sites, hts_pos_t *x_sites) {
    cand_snp_t *cand_snps = init_cand_snps(n_x_sites, x_sites);
    // 2nd pass: update snp_sites, calculate the depth and allele frequency of each site
    int start_x_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        start_x_i = update_cand_snps_from_digar(bam_chunk->digars+i, read, n_x_sites, x_sites, start_x_i, cand_snps);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_x_sites; ++i) {
            fprintf(stderr, "Mis: %lld\t", (long long) x_sites[i]);
            fprintf(stderr, "Depth: %d\t", cand_snps[i].n_depth);
            for (int j = 0; j < cand_snps[i].n_uniq_bases; ++j) {
                fprintf(stderr, "%c: %d\t", LONGCALLD_BAM_BASE_STR[cand_snps[i].bases[j]], cand_snps[i].base_covs[j]);
            }
            fprintf(stderr, "Low-Depth: %d\n", cand_snps[i].n_low_depth);
        }
    }
    return cand_snps;
}

// filter: 
//   1) usable X/= over total non-low-qual depth >= threshold
//   2) total non-low-qual depth / total depth >= threshold
int filter_cand_snps(cand_snp_t *cand_snps, int n_x_sites, const call_var_opt_t *opt) {
    // no need to malloc for cand_snp_pos
    int cand_snp_i = 0;
    int min_dp = opt->min_dp; double min_af = opt->min_af, max_af = opt->max_af, max_low_qual_frac = opt->max_low_qual_frac;
    for (int i = 0; i < n_x_sites; ++i) {
        if (cand_snps[i].n_depth < min_dp) continue;
        if ((double) cand_snps[i].n_low_depth / (cand_snps[i].n_low_depth + cand_snps[i].n_depth) > max_low_qual_frac) continue;
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
            } 
            fprintf(stderr, "Low-Depth: %d\n", cand_snps[i].n_low_depth);
        }
        // copy i'th to cand_snp_i'th
        if (i != cand_snp_i) {
            cand_snps[cand_snp_i].pos = cand_snps[i].pos;
            cand_snps[cand_snp_i].n_depth = cand_snps[i].n_depth;
            cand_snps[cand_snp_i].n_uniq_bases = cand_snps[i].n_uniq_bases;
            cand_snps[cand_snp_i].ref_base = cand_snps[i].ref_base;
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
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "CHUNK: %s\tbeg: %ld, end: %ld, ovlp_n: %d\n", bam_chunk->tname, (long) bam_chunk->beg, (long) bam_chunk->end, bam_chunk->n_ovlp_reads);
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "qname: %s, pos: %ld, end: %ld\n", bam_get_qname(read), (long) read->core.pos+1, (long) bam_endpos(read));
        if (read->core.qual < opt->min_mq) {
            bam_chunk->is_skipped[i] = 1; continue;
        }
        if (bam_chunk->bam_has_eqx_cigar) { // 1) look for Xs in cigar if =/X in cigar
            collect_digar_from_eqx_cigar(read, opt, bam_chunk->digars+i);
        } else if (bam_chunk->bam_has_md_tag) { // 2) look for mismatches in MD tag
            collect_digar_from_MD_tag(read, opt, bam_chunk->digars+i);
        } else { // 3) no =/X in cigar and no MD tag, compare bases with ref_seq
            if (pl->ref_seq == NULL) _err_error_exit("Reference genome FASTA file is required for BAM without =/X cigars and MD tags.");
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
    assign_hap(p, n_cand_snps, cand_snps, bam_chunk); 

    // write snps to VCF
    write_snp_to_vcf(cand_snps, n_cand_snps, pl->opt->out_vcf, pl->header->target_name[bam_chunk->tid]);
    free_read_snp_profile(p, bam_chunk->n_reads); free_cand_snps(cand_snps, n_cand_snps); 
    return 0;
}