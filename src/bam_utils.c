#include <ctype.h>
#include "bam_utils.h"
#include "utils.h"


read_snp_profile_t *init_read_snp_profile(int n_reads, int n_total_snps) {
    read_snp_profile_t *p = _err_malloc(n_reads * sizeof(read_snp_profile_t));
    for (int i = 0; i < n_reads; ++i) {
        p[i].read_id = i;
        p[i].start_snp_idx = -1; p[i].end_snp_idx = -2;
        p[i].snp_bases = (uint8_t*) _err_malloc(n_total_snps * sizeof(uint8_t));
        p[i].snp_qual = (uint8_t*) _err_malloc(n_total_snps * sizeof(uint8_t));
    }
    return p;
}

void free_read_snp_profile(read_snp_profile_t *p, int n_reads) {
    for (int i = 0; i < n_reads; ++i) {
        if (p[i].snp_bases) free(p[i].snp_bases);
        if (p[i].snp_qual) free(p[i].snp_qual);
    }
    free(p);
}

int has_equal_X_in_bam_cigar(bam1_t *read) {
    // Get the CIGAR string for the read
    const uint32_t *cigar = bam_get_cigar(read);
    int n_cigar = read->core.n_cigar;
    // Check if '=' or 'X' is in the CIGAR operations
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        if (op == BAM_CEQUAL || op == BAM_CDIFF) {
            return 1;
        } else if (op == BAM_CMATCH) {
            return 0;
        }
    }
    return 2;
}

int has_MD_in_bam(bam1_t *b) {
    uint8_t *s = bam_aux_get(b, "MD");
    return s != NULL;
}

void check_eqx_cigar_MD_tag(samFile *in_bam, bam_hdr_t *header, uint8_t *has_eqx, uint8_t *has_MD) {
    int n_to_check_reads = 10;
    // if any of n_to_check_reads reads do NOT have eqx or MD tag
    // set tag has NO
    bam1_t *b = bam_init1();
    int r, n_checked_reads = 0;
    *has_eqx = 1; *has_MD = 1;
    while ((r = sam_read1(in_bam, header, b)) >= 0) {
        if (has_equal_X_in_bam_cigar(b) == 0) *has_eqx = 0;
        if (has_MD_in_bam(b) == 0) *has_MD = 0;
        n_checked_reads++;
        if ((*has_eqx == 0 && *has_MD == 0) || n_checked_reads >= n_to_check_reads)
            break;
    }
    bam_destroy1(b);
    // rewind the file
}

int get_mis_site_index(hts_pos_t *mis_pos, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos) {
    int mis_i = -1;
    for (int i = cur_mis_i; i < n_total_mis_pos; i++) {
        if (mis_pos[i] == pos) {
            mis_i = i;
            break;
        }
    }
    if (mis_i == -1) _err_fatal("Mismatch position not found: %ld", pos);
    return mis_i;
}

int check_mis_site_index(hts_pos_t *mis_pos, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos, int *hit) {
    int mis_i = -1; *hit = 0;
    for (int i = cur_mis_i; i < n_total_mis_pos; i++) {
        if (mis_pos[i] >= pos) {
            mis_i = i;
            if (mis_pos[i] == pos) *hit = 1;
            break;
        }
    }
    if (mis_i == -1) mis_i = n_total_mis_pos;
        // _err_fatal("Mismatch position not found: %ld", pos);
    return mis_i;
}

int get_mis_pos_start(hts_pos_t *mis_pos, int n_total_mis_pos, hts_pos_t start) {
    int i;
    for (i = 0; i < n_total_mis_pos; ++i) {
        if (mis_pos[i] >= start) return i;
    }
    return i;
}

int get_cur_mis_pos_start(cand_snp_t *mis_sites, int cur_mis_i, int n_total_mis_pos, hts_pos_t start) {
    int i;
    for (i = cur_mis_i; i < n_total_mis_pos; ++i) {
        if (mis_sites[i].pos >= start) return i;
    }
    return i;
}

void merge_mis_site1(cand_snp_t *mis_site, uint8_t base) {
    mis_site->n_depth++;
    int base_i = mis_site->n_uniq_alt_bases, exist=0;
    for (int i = 0; i < mis_site->n_uniq_alt_bases; ++i) {
        if (mis_site->alt_bases[i] == base) {
            exist = 1; base_i = i;
        }
    }
    if (exist == 0) {
        mis_site->alt_bases[base_i] = base;
        mis_site->alt_base_covs[base_i] = 1;
        mis_site->n_uniq_alt_bases += 1;
    } else {
        mis_site->alt_base_covs[base_i] += 1;
    }
}

void merge_ref_site1(cand_snp_t *mis_site) {
    mis_site->n_depth++; mis_site->ref_base_cov += 1;
}

int merge_ref_sites(cand_snp_t *mis_sites, hts_pos_t *mis_pos, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos_start, hts_pos_t pos_end) {
    int new_mis_i = cur_mis_i;
    for (int i = cur_mis_i; i < n_total_mis_pos; i++) {
        if (mis_pos[i] >= pos_start && mis_pos[i] < pos_end) {
            merge_ref_site1(mis_sites+i);
            // printf("Merge ref site: %ld (%d): %d\n", mis_sites[i].pos, i, mis_sites[i].n_depth);
            new_mis_i = i+1;
        } else if (mis_pos[i] >= pos_end) {
            break;
        }
    }
    return new_mis_i;
}

int update_read_snp_profile_from_D_sites(cand_snp_t *mis_sites, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos_start, hts_pos_t pos_end, uint8_t qual, read_snp_profile_t *read_snp_profile, int *snp_idx) {
    int new_mis_i = cur_mis_i;
    for (int i = cur_mis_i; i < n_total_mis_pos; ++i) {
        if (mis_sites[i].pos >= pos_start && mis_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
            /* 6(D): ref allele */
            read_snp_profile->snp_bases[*snp_idx] = 6; 
            read_snp_profile->snp_qual[*snp_idx] = qual;
            (*snp_idx)++;
            new_mis_i = i+1;
        } else if (mis_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_mis_i;
}

int udpate_read_snp_profile_from_ref_sites(cand_snp_t *mis_sites, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos_start, hts_pos_t pos_end, const uint8_t *qual, read_snp_profile_t *read_snp_profile, int *snp_idx) {
    int new_mis_i = cur_mis_i;
    for (int i = cur_mis_i; i < n_total_mis_pos; ++i) {
        if (mis_sites[i].pos >= pos_start && mis_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
             /* 5(.): ref allele */
            read_snp_profile->snp_bases[*snp_idx] = 5; 
            read_snp_profile->snp_qual[*snp_idx] = qual[mis_sites[i].pos-pos_start];
            (*snp_idx)++;
            new_mis_i = i+1;
        } else if (mis_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_mis_i;
}

int update_read_snp_profile_from_mis_sites(cand_snp_t *mis_sites, int cur_mis_i, int n_total_mis_pos, hts_pos_t pos_start, hts_pos_t pos_end, int qi, const uint8_t *bam_seq, const uint8_t *qual, read_snp_profile_t *read_snp_profile, int *snp_idx) {
    int new_mis_i = cur_mis_i;
    for (int i = cur_mis_i; i < n_total_mis_pos; ++i) {
        if (mis_sites[i].pos >= pos_start && mis_sites[i].pos < pos_end) {
            if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = i;
            read_snp_profile->end_snp_idx = i;
            read_snp_profile->snp_bases[*snp_idx] = seq_nt16_int[bam_seqi(bam_seq, qi+mis_sites[i].pos-pos_start)];
            read_snp_profile->snp_qual[*snp_idx] = qual[qi+mis_sites[i].pos-pos_start];
            (*snp_idx)++;
            new_mis_i = i+1;
        } else if (mis_sites[i].pos >= pos_end) {
            break;
        }
    }
    return new_mis_i;
}

// input: mis_pos, n_total_mis_pos: start from 0
//        mis_start_i: previous start index
// return: index in mis_pos that is covered by current read (should be the start for next read)
int update_read_snp_profile_from_eqx_cigar(bam1_t *read, int n_total_mis_pos, cand_snp_t *mis_sites, int mis_start_i, read_snp_profile_t *read_snp_profile) {
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    const uint8_t *qual = bam_get_qual(read), *bam_seq = bam_get_seq(read);
    hts_pos_t pos = read->core.pos+1, qi = 0;

    int cur_start_i = get_cur_mis_pos_start(mis_sites, mis_start_i, n_total_mis_pos, pos);
    int snp_idx = 0;

    int mis_i = cur_start_i;
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        // printf("%d%c\t", len, BAM_CIGAR_STR[bam_cigar_op(cigar[i])]);
        if (op == BAM_CDIFF) {
            // for (int j = 0; j < len; ++j) {
            //     mis_i = check_mis_site_index(mis_pos, mis_i, n_total_mis_pos, pos, &hit);
            //     if (hit) {
            //         if (read_snp_profile->start_snp_idx == -1) read_snp_profile->start_snp_idx = mis_i;
            //         read_snp_profile->end_snp_idx = mis_i;
            //         read_snp_profile->snp_bases[snp_idx] = seq_nt16_int[bam_seqi(bam_get_seq(read), qi)];
            //         read_snp_profile->snp_qual[snp_idx] = qual[qi];
            //         // printf("SNP: %ld : %d\n", mis_sites[mis_i].pos, read_snp_profile->snp_bases[snp_idx]);
            //         snp_idx++;
            //         if (snp_idx != read_snp_profile->end_snp_idx - read_snp_profile->start_snp_idx + 1)
            //             printf("SNP index error: %d, %d, %d\n", snp_idx, read_snp_profile->end_snp_idx, read_snp_profile->start_snp_idx);
            //             // _err_fatal("SNP index error: %d, %d, %d", snp_idx, read_snp_profile->end_snp_idx, read_snp_profile->start_snp_idx);
            //     }
            //     pos++; qi++;
            // }
            mis_i = update_read_snp_profile_from_mis_sites(mis_sites, mis_i, n_total_mis_pos, pos, pos+len, qi, bam_seq, qual, read_snp_profile, &snp_idx);
            pos+=len; qi+=len;
        }
        if (op == BAM_CEQUAL) {
            mis_i = udpate_read_snp_profile_from_ref_sites(mis_sites, mis_i, n_total_mis_pos, pos, pos+len, qual+qi, read_snp_profile, &snp_idx);
            pos += len; qi += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            uint8_t _qual;
            if (qi > 0) _qual = (qual[qi] + qual[qi-1]) / 2; else _qual = qual[qi];
            mis_i = update_read_snp_profile_from_D_sites(mis_sites, mis_i, n_total_mis_pos, pos, pos+len, _qual, read_snp_profile, &snp_idx);
            pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) { // XXX
            qi += len;
        } else if (op == BAM_CMATCH) {
            _err_fatal("CIGAR operation 'M' is not expected in EQX CIGAR: %s", bam_get_qname(read));
        }
        // printf("op: %c, mis_i: %d\n", BAM_CIGAR_STR[op], mis_i);
        if (pos > mis_sites[n_total_mis_pos-1].pos || mis_i >= n_total_mis_pos) break;
    }
    return cur_start_i;
}

int update_mis_sites_from_eqx_cigar(bam1_t *read, int min_bq, int n_total_mis_pos, hts_pos_t *mis_pos, cand_snp_t *mis_sites) {
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    const uint8_t *qual = bam_get_qual(read);
    hts_pos_t pos = read->core.pos+1, qi = 0;

    int n_skipped_mis = get_mis_pos_start(mis_pos, n_total_mis_pos, pos);

    int mis_i = n_skipped_mis;
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]), len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                if (qual[qi] >= min_bq) {
                    // update
                    mis_i = get_mis_site_index(mis_pos, mis_i, n_total_mis_pos, pos);
                    // merge the current mismatch base to the mis_site
                    merge_mis_site1(mis_sites+mis_i, seq_nt16_int[bam_seqi(bam_get_seq(read), qi)]);
                    // printf("Merge mis site: %ld (%d): %d\n", mis_sites[mis_i].pos, mis_i, mis_sites[mis_i].n_depth);
                }
                pos++; qi++;
            }
        }
        if (op == BAM_CEQUAL) {
            // if mis_pos is within [pos,pos+len]
            // update
            mis_i = merge_ref_sites(mis_sites, mis_pos, mis_i, n_total_mis_pos, pos, pos+len);
            pos += len; qi += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            // XXX no action if D cover any mismatch sites
            pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qi += len;
        } else if (op == BAM_CMATCH) {
            _err_fatal("CIGAR operation 'M' is not expected in EQX CIGAR: %s", bam_get_qname(read));
        }
        // printf("op: %c, mis_i: %d\n", BAM_CIGAR_STR[op], mis_i);
        if (pos > mis_pos[n_total_mis_pos-1] || mis_i >= n_total_mis_pos-1) break;
    }
    return n_skipped_mis;
}
// collect postions for all X bases from EQX CIGAR
// store positions into mis_pos, return the number of X bases
int get_mis_bases_from_eqx_cigar(bam1_t *read, int min_bq, mis_base_t **mis_bases) {
    const uint32_t *cigar = bam_get_cigar(read); int n_cigar = read->core.n_cigar;
    const uint8_t *qual = bam_get_qual(read);
    int n_mis_bases = 0, m_mis_bases = 0;
    hts_pos_t pos = read->core.pos+1, qi = 0;

    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        
        if (op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                if (qual[qi] >= min_bq) {
                    if (mis_bases != NULL) {
                        // Reallocate memory for mis_pos
                        _uni_realloc(*mis_bases, n_mis_bases, m_mis_bases, mis_base_t);
                        // Store the current position
                        (*mis_bases)[n_mis_bases].pos = pos;
                        (*mis_bases)[n_mis_bases].base = seq_nt16_int[bam_seqi(bam_get_seq(read), qi)];
                        (*mis_bases)[n_mis_bases].qual = qual[qi];
                    }
                    n_mis_bases++;
                }
                pos++; qi++;
            }
        }
        if (op == BAM_CMATCH || op == BAM_CEQUAL) {
            pos += len; qi += len;
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) {
            pos += len;
        } else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qi += len;
        }
    }
    return n_mis_bases;
}

int get_mis_bases_from_MD_tag(bam1_t *read, int min_bq, mis_base_t **mis_bases) {
    int n_mis_bases = 0, m_mis_bases = 0;
    hts_pos_t pos = read->core.pos+1, qi = 0;

    uint8_t *s = bam_aux_get(read, "MD");
    if (s == NULL) {
        err_fatal_core(__func__, "MD tag not found in the BAM file: %s", bam_get_qname(read));
        return 0;
    }

    const uint8_t *qual = bam_get_qual(read);
    uint8_t *qual_without_insertions = NULL;
    int n_cigar = read->core.n_cigar;
    const uint32_t *cigar = bam_get_cigar(read);
    qual_without_insertions = (uint8_t*)malloc(read->core.l_qseq);
    int qi_without_insertions = 0;
    for (int i = 0; i < n_cigar; i++) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                qual_without_insertions[qi_without_insertions] = qual[qi];
                qi_without_insertions++;
                qi++;
            }
        }
        if (op == BAM_CINS || op == BAM_CSOFT_CLIP) {
            qi += len;
        }
    }

    char *md = bam_aux2Z(s);
    // printf("MD: %s\n", md);
    int md_i = 0; qi = 0;
    while (md[md_i]) {
        if (isdigit(md[md_i])) { // This is a run of matches
            int eq_len = strtol(&md[md_i], &md, 10);
            pos += eq_len;
            qi += eq_len;
            md_i = 0; // start of the next run
        } else if (isalpha(md[md_i])) { // This is a mismatch
            if (qual_without_insertions[qi] >= min_bq) {
                if (mis_bases != NULL) {
                    // Reallocate memory for mis_pos
                    _uni_realloc(*mis_bases, n_mis_bases, m_mis_bases, mis_base_t);
                    // Store the current position
                    (*mis_bases)[n_mis_bases].pos = pos;
                    (*mis_bases)[n_mis_bases].base = seq_nt16_int[bam_seqi(bam_get_seq(read), qi)];
                    (*mis_bases)[n_mis_bases].qual = qual[qi];
                }
                n_mis_bases++;
            }
            pos++; qi++; md_i++;
        } else { // ^: start of deletions
            md_i++;
            while (md[md_i] && isalpha(md[md_i])) {
                pos++; md_i++;
            }
        }
    }
    free(qual_without_insertions);
    return n_mis_bases;
}

int get_mis_bases_from_ref_seq(bam1_t *read, int min_bq, kstring_t *ref_seq, mis_base_t **mis_bases) {
    int n_mis_bases = 0, m_mis_bases = 0;
    int n_cigar = read->core.n_cigar;
    const uint32_t *cigar = bam_get_cigar(read);
    const uint8_t *qual = bam_get_qual(read);
    hts_pos_t pos = read->core.pos+1, qi = 0;

    for (int i = 0; i < n_cigar; ++i) {
        int op = bam_cigar_op(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);
        if (op == BAM_CMATCH || op == BAM_CEQUAL || op == BAM_CDIFF) {
            for (int j = 0; j < len; ++j) {
                if (qual[qi] >= min_bq) {
                    // Get the reference base
                    if (pos > ref_seq->l) err_fatal_core(__func__, "Read exceed reference sequence length: %s", bam_get_qname(read));
                    char ref_base = ref_seq->s[pos-1];
                    // Get the read base
                    char read_base = seq_nt16_str[bam_seqi(bam_get_seq(read), qi)];
                    if (ref_base != read_base) {
                        if (mis_bases != NULL) {
                            // Reallocate memory for mis_pos
                            _uni_realloc(*mis_bases, n_mis_bases, m_mis_bases, mis_base_t);
                            // Store the current position
                            (*mis_bases)[n_mis_bases].pos = pos;
                            (*mis_bases)[n_mis_bases].base = seq_nt16_int[bam_seqi(bam_get_seq(read), qi)];
                            (*mis_bases)[n_mis_bases].qual = qual[qi];
                        }
                        n_mis_bases++;
                    }
                }
                pos++; qi++;
            }
        } else if (op == BAM_CDEL || op == BAM_CREF_SKIP) pos += len;
        else if (op == BAM_CINS || op == BAM_CSOFT_CLIP) qi += len;
    }

    return n_mis_bases;
}

int bam_read_chunk_init(bam_chunk_t *chunk, int n_reads) {
    chunk->reads = (bam1_t**)malloc(n_reads * sizeof(bam1_t*));
    chunk->is_skipped = (uint8_t*)calloc(n_reads, sizeof(uint8_t));
    chunk->HPs = (int*)malloc(n_reads * sizeof(int));
    for (int i = 0; i < n_reads; i++) {
        chunk->reads[i] = bam_init1();
    }
    chunk->n_reads = 0; chunk->m_reads = n_reads;
    return 0;
}

int bam_read_chunk_free(bam_chunk_t *chunk) {
    for (int i = 0; i < chunk->m_reads; i++) {
        bam_destroy1(chunk->reads[i]);
    }
    free(chunk->reads); free(chunk->is_skipped); free(chunk->HPs);
    return 0;
}

int bam_read_chunk(samFile *in_bam, bam_hdr_t *header, bam_chunk_t *chunk, int max_reads_per_chunk) {
    int r, has_eqx_cigar=0, has_MD=0;
    bam_read_chunk_init(chunk, max_reads_per_chunk);
    int tid=-1; uint64_t beg=0, _end=0, end; char *tname=NULL;
    while (chunk->n_reads < max_reads_per_chunk) {
        if ((r = sam_read1(in_bam, header, chunk->reads[chunk->n_reads])) < 0) break;
        // skip unmapped/secondary alignments
        if (chunk->reads[chunk->n_reads]->core.flag & (BAM_FUNMAP | BAM_FSECONDARY)) continue; // BAM_FSUPPLEMENTARY
        if (tid == -1) {
            tid = chunk->reads[chunk->n_reads]->core.tid;
            tname = header->target_name[tid];
            beg = chunk->reads[chunk->n_reads]->core.pos;
        } else if (tid != chunk->reads[chunk->n_reads]->core.tid) break; // different chromosome
        else { // same chromosome
            end = chunk->reads[chunk->n_reads]->core.pos + bam_cigar2rlen(chunk->reads[chunk->n_reads]->core.n_cigar, bam_get_cigar(chunk->reads[chunk->n_reads]));
            if (end > _end) _end = end;
        }
        // check eqx cigar or MD tag
        if (has_eqx_cigar == 0) has_eqx_cigar = has_equal_X_in_bam_cigar(chunk->reads[chunk->n_reads]);
        if (has_MD == 0) has_MD = has_MD_in_bam(chunk->reads[chunk->n_reads]);
        chunk->n_reads++;
    }
    if (chunk->n_reads > 0) {
        chunk->tid = tid; chunk->tname = tname;
        chunk->start = beg; chunk->end = _end;
        chunk->bam_has_eqx_cigar = has_eqx_cigar;
        chunk->bam_has_md_tag = has_MD;
    } else {
        bam_read_chunk_free(chunk);
    }
    return chunk->n_reads;
}