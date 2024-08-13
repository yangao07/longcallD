#include <stdio.h>
#include <stdlib.h>
#include "htslib/sam.h"
#include "collect_var.h"
#include "call_var.h"
#include "utils.h"
#include "bam_utils.h"
#include "seq.h"
#include "assign_hap.h"
#include "vcf_utils.h"

extern int LONGCALLD_VERBOSE;
char test_read_name[1024] = "m84039_230928_213653_s3/23396787/ccs";

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

cand_var_t *init_cand_vars(int n_var_sites, var_site_t *var_sites) {
    cand_var_t *cand_vars = (cand_var_t*)malloc(n_var_sites * sizeof(cand_var_t));
    for (int i = 0; i < n_var_sites; ++i) {
        // static information
        // from var_sites
        cand_vars[i].tid = var_sites[i].tid;
        cand_vars[i].pos = var_sites[i].pos; 
        cand_vars[i].var_type = var_sites[i].var_type;
        cand_vars[i].ref_len = var_sites[i].ref_len;

        cand_vars[i].phase_set = 0; // unset
        cand_vars[i].n_depth = 0; cand_vars[i].n_low_depth = 0; 
        cand_vars[i].n_uniq_alles = 1; // ref allele
        cand_vars[i].alle_covs = (int*)calloc(1, sizeof(int)); // ref+alt
        cand_vars[i].alt_alle_seqs = NULL; cand_vars[i].alt_alle_seq_len = NULL; // only alt

        // dynamic information, allocate and update during haplotype assignment
        cand_vars[i].alle_to_hap = NULL; cand_vars[i].hap_to_alle_profile = NULL; cand_vars[i].hap_to_cons_alle = NULL;
        cand_vars[i].is_low_qual = 0; cand_vars[i].is_skipped = 0;
    }
    return cand_vars;
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

void free_cand_vars(cand_var_t *cand_vars, int m) {
    if (m > 0) {
        for (int i = 0; i < m; ++i) {
            free(cand_vars[i].alle_covs); free(cand_vars[i].alt_alle_seq_len);
            for (int j = 0; j < cand_vars[i].n_uniq_alles-1; ++j) {
                if (cand_vars[i].alt_alle_seqs[j] != NULL) free(cand_vars[i].alt_alle_seqs[j]);
            } free(cand_vars[i].alt_alle_seqs);
            if (cand_vars[i].alle_to_hap != NULL) free(cand_vars[i].alle_to_hap);
            if (cand_vars[i].hap_to_alle_profile != NULL) {
                for (int j = 0; j <= LONGCALLD_DEF_PLOID; ++j) {
                    free(cand_vars[i].hap_to_alle_profile[j]);
                } free(cand_vars[i].hap_to_alle_profile);
            }
            if (cand_vars[i].hap_to_cons_alle != NULL) free(cand_vars[i].hap_to_cons_alle);
        }
        free(cand_vars);
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

cand_var_t *collect_cand_vars(bam_chunk_t *bam_chunk, int n_var_sites, var_site_t *var_sites) {
    cand_var_t *cand_vars = init_cand_vars(n_var_sites, var_sites);
    // 2nd pass: update snp_sites, calculate the depth and allele frequency of each site
    int start_var_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        start_var_i = update_cand_vars_from_digar(bam_chunk->digars+i, read, n_var_sites, var_sites, start_var_i, cand_vars);
    }
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < n_var_sites; ++i) {
            fprintf(stderr, "CandVar: %lld\t", (long long) var_sites[i].pos);
            fprintf(stderr, "Type: %c\t", BAM_CIGAR_STR[var_sites[i].var_type]);
            fprintf(stderr, "RefLen: %d\t", var_sites[i].ref_len);
            fprintf(stderr, "Depth: %d\t", cand_vars[i].n_depth);
            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                fprintf(stderr, "%d ", j);
                if (j != 0) {
                    for (int k = 0; k < cand_vars[i].alt_alle_seq_len[j-1]; ++k) {
                        fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_alle_seqs[j-1][k]]);
                    }
                } fprintf(stderr, ": ");
                fprintf(stderr, "%d\t", cand_vars[i].alle_covs[j]);
            }
            fprintf(stderr, "Low-Depth: %d\n", cand_vars[i].n_low_depth);
        }
    }
    return cand_vars;
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

// XXX skip detection for long ins/del
int is_repeat_region(kstring_t *ref_seq, cand_var_t *var) {
    hts_pos_t pos = var->pos; int ref_len = var->ref_len; int var_type = var->var_type;
    uint8_t *ref_bseq, *alt_bseq; int len = 0;
    int is_repeat = 1;
    if (var_type == BAM_CDEL) { // [pos, pos+del_len*N] vs [pos+del_len, pos+del_len*(N+1)]
        // nst_nt4_table['A'] -> 0, 'C' -> 1, 'G' -> 2, 'T' -> 3, 'N' -> 4
        int del_len = ref_len-1;
        len = del_len * 3; // see if del seq is 3-fold repeat
        ref_bseq = get_bseq(ref_seq->s+pos-1, len);
        alt_bseq = get_bseq(ref_seq->s+pos-1+del_len, len);
        for (int i = 0; i < len; ++i) {
            if (ref_bseq[i] != alt_bseq[i]) { is_repeat = 0; break; }
        }
        free(ref_bseq); free(alt_bseq);
    } else { // if (var_type == BAM_CINS) { // [ins] * N vs [pos, pos+ ins_len * N]
        for (int i = 1; i < var->n_uniq_alles; ++i) {
            int ins_len = var->alt_alle_seq_len[i-1];
            len = ins_len * 3; // see if ins seq is 3-fold repeat
            ref_bseq = get_bseq(ref_seq->s+pos-1, len);
            alt_bseq = get_bseq(ref_seq->s+pos-1-ins_len, len);
            for (int j = 0; j < ins_len; ++j) alt_bseq[j] = var->alt_alle_seqs[i-1][j];
            is_repeat = 1;
            for (int k = 0; k < len; ++k) {
                if (ref_bseq[k] != alt_bseq[k]) { is_repeat = 0; break; }
            }
            free(ref_bseq); free(alt_bseq);
            if (is_repeat == 1) break;
        }
    }
    if (is_repeat) {
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "RepeatRegion: %lld, type: %c, refLen: %d\n", (long long) pos, BAM_CIGAR_STR[var_type], ref_len);
    }
    return is_repeat;
}

// XXX filter out var sites with too many alt alleles? (repeats)
// filter indels in repeat regions
// only keep het. vars
int filter_cand_vars(cand_var_t *cand_vars, int n_var_sites, kstring_t *ref_seq, const call_var_opt_t *opt) {
    // no need to malloc for cand_snp_pos
    int cand_var_i = 0;
    int min_dp = opt->min_dp; double min_af = opt->min_af, max_af = opt->max_af, max_low_qual_frac = opt->max_low_qual_frac;
    for (int i = 0; i < n_var_sites; ++i) {
        if (cand_vars[i].n_depth < min_dp) continue;
        if ((double) cand_vars[i].n_low_depth / (cand_vars[i].n_low_depth + cand_vars[i].n_depth) > max_low_qual_frac) continue;
        double _af1 = 0.0, _af2 = 0.0, _af = 0.0;
        for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
            _af = (double) cand_vars[i].alle_covs[j] / cand_vars[i].n_depth;
            if (_af > _af1) { _af2 = _af1; _af1 = _af; } 
            else if (_af > _af2 && _af <= _af1) _af2 = _af;
        }
        if (_af1 < min_af || _af1 > max_af || _af2 < min_af || _af2 > max_af) continue; // XXX
        if (cand_vars[i].var_type == BAM_CINS || cand_vars[i].var_type == BAM_CDEL) {
            if (is_repeat_region(ref_seq, cand_vars+i)) continue;
        }
        if (LONGCALLD_VERBOSE >= 2) {
            fprintf(stderr, "FilteredCandVar %d: %lld, type: %c, refLen: %d, depth: %d, AF: %f|%f\t", cand_var_i, (long long) cand_vars[i].pos, BAM_CIGAR_STR[cand_vars[i].var_type], cand_vars[i].ref_len, cand_vars[i].n_depth, _af1, _af2);
            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                fprintf(stderr, "%d ", j);
                if (j != 0) {
                    for (int k = 0; k < cand_vars[i].alt_alle_seq_len[j-1]; ++k) {
                        fprintf(stderr, "%c", "ACGTN"[cand_vars[i].alt_alle_seqs[j-1][k]]);
                    }
                } fprintf(stderr, ": ");
                fprintf(stderr, "%d\t", cand_vars[i].alle_covs[j]);
            }
            fprintf(stderr, "Low-Depth: %d\n", cand_vars[i].n_low_depth);
        }
        // copy i'th to cand_snp_i'th
        if (i != cand_var_i) {
            cand_vars[cand_var_i].pos = cand_vars[i].pos; cand_vars[cand_var_i].var_type = cand_vars[i].var_type;
            cand_vars[cand_var_i].n_depth = cand_vars[i].n_depth; cand_vars[cand_var_i].n_low_depth = cand_vars[i].n_low_depth;
            cand_vars[cand_var_i].ref_len = cand_vars[i].ref_len; 
            cand_vars[cand_var_i].is_low_qual = cand_vars[i].is_low_qual; cand_vars[cand_var_i].is_skipped = cand_vars[i].is_skipped;

            // allele_covs, alt_alle_seqs, alt_alle_seq_len
            if (cand_vars[i].n_uniq_alles > cand_vars[cand_var_i].n_uniq_alles) {
                cand_vars[cand_var_i].alle_covs = (int*)realloc(cand_vars[cand_var_i].alle_covs, cand_vars[i].n_uniq_alles * sizeof(int));
                cand_vars[cand_var_i].alt_alle_seq_len = (int*)realloc(cand_vars[cand_var_i].alt_alle_seq_len, (cand_vars[i].n_uniq_alles-1) * sizeof(int));
            }
            for (int j = 0; j < cand_vars[cand_var_i].n_uniq_alles-1; ++j) {
                free(cand_vars[cand_var_i].alt_alle_seqs[j]); 
            } free(cand_vars[cand_var_i].alt_alle_seqs);
            cand_vars[cand_var_i].alt_alle_seqs = (uint8_t**)malloc((cand_vars[i].n_uniq_alles-1) * sizeof(uint8_t*));

            for (int j = 0; j < cand_vars[i].n_uniq_alles; ++j) {
                cand_vars[cand_var_i].alle_covs[j] = cand_vars[i].alle_covs[j];
            }
            for (int j = 0; j < cand_vars[i].n_uniq_alles-1; j++) {
                cand_vars[cand_var_i].alt_alle_seq_len[j] = cand_vars[i].alt_alle_seq_len[j];
                cand_vars[cand_var_i].alt_alle_seqs[j] = (uint8_t*)malloc(cand_vars[i].alt_alle_seq_len[j] * sizeof(uint8_t));
                for (int k = 0; k < cand_vars[i].alt_alle_seq_len[j]; ++k) {
                    cand_vars[cand_var_i].alt_alle_seqs[j][k] = cand_vars[i].alt_alle_seqs[j][k];
                }
            } 
            cand_vars[cand_var_i].n_uniq_alles = cand_vars[i].n_uniq_alles;
        }
        cand_var_i++;
        // cand_snp_pos[cand_snp_i++] = cand_snps[i].pos;
    }
    for (int i = cand_var_i; i < n_var_sites; ++i) {
        free(cand_vars[i].alle_covs); free(cand_vars[i].alt_alle_seq_len); 
        for (int j = 0; j < cand_vars[i].n_uniq_alles-1; ++j) free(cand_vars[i].alt_alle_seqs[j]);
        free(cand_vars[i].alt_alle_seqs);
    }
    return cand_var_i;
}
void collect_digars_from_bam(bam_chunk_t *bam_chunk, const call_var_pl_t *pl) {
    const call_var_opt_t *opt = pl->opt;
    if (LONGCALLD_VERBOSE >= 2)
        fprintf(stderr, "CHUNK: %s\tbeg: %ld, end: %ld, total_n: %d, ovlp_n: %d\n", bam_chunk->tname, (long) bam_chunk->beg, (long) bam_chunk->end, bam_chunk->n_reads, bam_chunk->n_up_ovlp_reads);
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "%d: qname: %s, flag: %d, pos: %ld, end: %ld\n", i, bam_get_qname(read), read->core.flag, (long) read->core.pos+1, (long) bam_endpos(read));
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        if (read->core.qual < opt->min_mq) {
            bam_chunk->is_skipped[i] = 1; 
            fprintf(stderr, "Skip read %s with low mapping quality %d\n", bam_get_qname(read), read->core.qual);
            continue;
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

// XXX should be digar->pos-1 for INS/DEL
var_site_t make_var_site(int tid, digar1_t *digar) {
    var_site_t var_site = {tid, digar->pos, digar->type, 0};
    if (digar->type == BAM_CDIFF) {
        var_site.ref_len = 1;
    } else if (digar->type == BAM_CINS) {
        var_site.ref_len = 1;
    } else if (digar->type == BAM_CDEL) {
        var_site.ref_len = 1 + digar->len;
    }
    return var_site;
}

int comp_var_site(var_site_t *var1, var_site_t *var2) {
    if (var1->pos < var2->pos) return -1;
    if (var1->pos > var2->pos) return 1;
    if (var1->var_type < var2->var_type) return -1;
    if (var1->var_type > var2->var_type) return 1;
    if (var1->ref_len < var2->ref_len) return -1;
    if (var1->ref_len > var2->ref_len) return 1;
    return 0;
}

// order of variants with same pos: order by type and ref_len
int merge_var_sites(int n_total_var_sites, var_site_t **var_sites, int tid, hts_pos_t beg, hts_pos_t end, int n_digar, digar1_t *digars) {
    int new_total_var_sites = 0;
    var_site_t *new_var_sites = (var_site_t*)malloc((n_total_var_sites + n_digar) * sizeof(var_site_t));
    int i, j;
    for (i = j = 0; i < n_total_var_sites && j < n_digar; ) {
        // if (digars[j].pos == 10829597 || (*var_sites)[i].pos == 10829597)
            // fprintf(stderr, "i: %d, j: %d\n", i, j);
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
            var_site_t digar_var_site = make_var_site(tid, digars+j);
            if (beg != -1 && digar_var_site.pos < beg) { j++; continue; }
            else if (end != -1 && digar_var_site.pos > end) break;
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
        if (!digars[j].is_low_qual && (digars[j].type == BAM_CDIFF || digars[j].type == BAM_CINS || digars[j].type == BAM_CDEL)) {
            var_site_t digar_var_site = make_var_site(tid, digars+j);
            new_var_sites[new_total_var_sites++] = digar_var_site;
        }
    }
    free(*var_sites); *var_sites = new_var_sites;
    return new_total_var_sites;
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

int collect_var_sites(bam_chunk_t *bam_chunk, var_site_t **var_sites) {
    int n_total_var_sites = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
       n_total_var_sites = merge_var_sites(n_total_var_sites, var_sites, bam_chunk->tid, bam_chunk->beg, bam_chunk->end, bam_chunk->digars[i].n_digar, bam_chunk->digars[i].digars);
    }
    // print cand_snps
    if (LONGCALLD_VERBOSE >= 2) {
        fprintf(stderr, "Total candidate variant sites: %d\n", n_total_var_sites);
        for (int i = 0; i < n_total_var_sites; ++i) {
            fprintf(stderr, "CandVarSite: %s:%lld\t%d\t%c\n", bam_chunk->tname, (long long) (*var_sites)[i].pos, (*var_sites)[i].ref_len, BAM_CIGAR_STR[(*var_sites)[i].var_type]);
        }
    }
    return n_total_var_sites;
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

read_var_profile_t *collect_read_var_profile(bam_chunk_t *bam_chunk) {
    int n_cand_vars = bam_chunk->n_cand_vars;
    cand_var_t *cand_vars = bam_chunk->cand_vars;
    read_var_profile_t *p = init_read_var_profile(bam_chunk->n_reads, n_cand_vars);
    // 3rd pass: collect read-wise SNP profiles
    int start_var_i = 0;
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        if (bam_chunk->is_skipped[i]) continue;
        bam1_t *read = bam_chunk->reads[i];
        // if (strcmp(test_read_name, bam_get_qname(read)) == 0)
            // fprintf(stderr, "Read: %s\n", bam_get_qname(read));
        start_var_i = update_read_var_profile_from_digar(bam_chunk->digars+i, read, n_cand_vars, cand_vars, start_var_i, p+i);
            // fprintf(stderr, "start_snp_i: %d\n", start_snp_i);
        if (LONGCALLD_VERBOSE >= 2) {
            if (p[i].start_var_idx >= 0) {
                fprintf(stderr, "Read: %s, start_var_i: %d, end_var_i: %d\n", bam_get_qname(read), p[i].start_var_idx, p[i].end_var_idx);
                for (int j = 0; j <= p[i].end_var_idx-p[i].start_var_idx; ++j) {
                    // fprintf(stderr, "\tSNP: (%d) %lld, base: %c, qual: %d\n", j, (long long) cand_snps[j].pos, 
                    // LONGCALLD_BAM_BASE_STR[p[i].snp_bases[j-p[i].start_snp_idx]], p[i].snp_qual[j-p[i].start_snp_idx]);
                    if (p[i].var_is_used[j] == 0) continue;
                    fprintf(stderr, "\tVar: (%d) %lld", j, (long long) cand_vars[j+p[i].start_var_idx].pos);
                    fprintf(stderr, ",%c, allele: %d\n", BAM_CIGAR_STR[cand_vars[j+p[i].start_var_idx].var_type], p[i].alleles[j]);
                    // fprintf(stderr, ", qual: %d\n", p[i].snp_qual[j]);
                }
            }
        }
    }
    return p;
}

// collect variants based on hap_to_cons_alle
// filter out low-quality variants, e.g., P(var|hap,phasing)
int make_variants(cand_var_t *cand_vars, int n_cand_vars, kstring_t *ref_seq, var_t *var) {
    if (n_cand_vars <= 0) return 0;
    var->n = 0; var->m = n_cand_vars;
    int i = 0; int hap_alles[2];
    var->vars = (var1_t*)malloc(n_cand_vars * sizeof(var1_t));
    for (int cand_i = 0; cand_i < n_cand_vars; ++cand_i) {
        hap_alles[0] = cand_vars[cand_i].hap_to_cons_alle[1];
        hap_alles[1] = cand_vars[cand_i].hap_to_cons_alle[2];
        // only keep het. vars
        if (hap_alles[0] == -1 || hap_alles[1] == -1) continue;
        if (hap_alles[0] == hap_alles[1]) continue;

        var->vars[i].type = cand_vars[cand_i].var_type;
        var->vars[i].PS = cand_vars[cand_i].phase_set;
        var->vars[i].ref_len = cand_vars[cand_i].ref_len;
        if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) var->vars[i].pos = cand_vars[cand_i].pos-1;
        else var->vars[i].pos = cand_vars[cand_i].pos;
        var->vars[i].ref_bases = get_bseq(ref_seq->s+var->vars[i].pos-1, var->vars[i].ref_len);
        var->vars[i].alt_len = (int*)malloc(2 * sizeof(int));
        var->vars[i].alt_bases = (uint8_t**)malloc(2 * sizeof(uint8_t*));
        var->vars[i].n_alt_allele = 0;
        for (int hap=1; hap <= 2; ++hap) {
            int hap_alle = hap_alles[hap-1];
            if (hap_alle != 0) {
                int alt_len = cand_vars[cand_i].alt_alle_seq_len[hap_alle-1];
                if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) alt_len += 1;
                var->vars[i].alt_len[var->vars[i].n_alt_allele] = alt_len;
                var->vars[i].alt_bases[var->vars[i].n_alt_allele] = (uint8_t*)malloc(alt_len * sizeof(uint8_t));
                if (var->vars[i].type == BAM_CDEL || var->vars[i].type == BAM_CINS) {
                    var->vars[i].alt_bases[var->vars[i].n_alt_allele][0] = nst_nt4_table[(int)ref_seq->s[var->vars[i].pos-1]];
                    for (int j = 0; j < alt_len-1; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j+1] = cand_vars[cand_i].alt_alle_seqs[hap_alle-1][j];
                    }
                } else {
                    for (int j = 0; j < alt_len; ++j) {
                        var->vars[i].alt_bases[var->vars[i].n_alt_allele][j] = cand_vars[cand_i].alt_alle_seqs[hap_alle-1][j];
                    }
                }
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
    var->n = i;
    return i;
}

void collect_var_main(const call_var_pl_t *pl, bam_chunk_t *bam_chunk, var_t *var) {
    // collect X/I/D sites from BAM
    collect_digars_from_bam(bam_chunk, pl);
    
    // merge var sites from all reads
    var_site_t *var_sites = NULL; int n_var_sites;
    if ((n_var_sites = collect_var_sites(bam_chunk, &var_sites)) <= 0) return;
    
    // collect reference and alternative alleles for all var sites
    bam_chunk->cand_vars = collect_cand_vars(bam_chunk, n_var_sites, var_sites); free(var_sites);
    
    // filter candidate vars based depth, allele frequency, etc.
    if ((bam_chunk->n_cand_vars = filter_cand_vars(bam_chunk->cand_vars, n_var_sites, pl->ref_seq->seq + ref_seq_name2id(pl->ref_seq, bam_chunk->tname), pl->opt)) <= 0) return;
    
    // collect read-wise var profiles
    bam_chunk->read_var_profile = collect_read_var_profile(bam_chunk);

    // assign hap to each read and SNP
    assign_hap_based_on_cand_vars(bam_chunk);
}

// stitch ii and ii+1
void stitch_var_main(const call_var_step_t *step, bam_chunk_t *bam_chunk, var_t *var, long ii) {
    const call_var_pl_t *pl = step->pl;
    ref_seq_t *ref_seq = pl->ref_seq;
    if (ii == 0) {
        // no stitching needed, copy variants directly
        // collect variants: only keep high-quality het. variants, filter out low-quality variants (depth, allele frequency, phasing)
        // within each bam chunk (chr, start, end), all variants called are based on all supporting reads, no other bam chunks are needed
        // between two adjacent bam chunks, we need to stitch the variants (flip haplotypes) later using the overlapping reads
        make_variants(bam_chunk->cand_vars, bam_chunk->n_cand_vars, ref_seq->seq + ref_seq_name2id(ref_seq, bam_chunk->tname), var);
        return;
    } else {
        // stitch variants of ii-1 and ii's bam chunk
        // 1) find overlapping reads between two adjacent bam chunks
        bam_chunk_t *prev_bam_chunk = step->chunks+ii-1;
        make_variants(bam_chunk->cand_vars, bam_chunk->n_cand_vars, ref_seq->seq + ref_seq_name2id(ref_seq, bam_chunk->tname), var);
        // stick_variants(prev_bam_chunk, bam_chunk, var);
        // 2) flip haplotypes for variants in the overlapping region
        // 3) merge variants from two adjacent bam chunks
    }
    // 4) write variants to VCF
    // write_var_to_vcf(cand_vars, n_cand_vars, pl->opt->out_vcf, pl->header->target_name[bam_chunk->tid]);
}