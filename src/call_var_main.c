#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#ifdef __linux__
#include <sys/sysinfo.h>  // Linux-specific
#elif __APPLE__
#include <sys/types.h>
#include <sys/sysctl.h>    // macOS-specific
#include <mach/mach.h>     // For memory information
#else
#error "Unsupported platform"
#endif
#include "main.h"
#include "call_var_main.h"
#include "bam_utils.h"
#include "vcf_utils.h"
#include "utils.h"
#include "seq.h"
#include "collect_var.h"
#include "kthread.h"
#include "align.h"
#include "math_utils.h"
#include "kmer.h"

extern int LONGCALLD_VERBOSE;

const struct option call_var_opt [] = {
    // long options
    { "amb-base", 0, NULL, 0},
    { "hifi", 0, NULL, 0},
    { "ont", 0, NULL, 0},
    { "region-file", 1, NULL, 0},
    { "regions-file", 1, NULL, 0},
    { "all-ctg", 0, NULL, 0},
    { "autosome-XY", 0, NULL, 0},
    { "autosome", 0, NULL, 0},
    { "mosaic", 0, NULL, 0},
    { "refine-aln", 0, NULL, 0},
    { "alpha", 1, NULL, 0},
    { "beta", 1, NULL, 0},
    { "max-somvar", 1, NULL, 0},
    { "som-alt", 1, NULL, 0},
    { "som-mei-alt", 1, NULL, 0},

    { "exclude-ctg", 1, NULL, 'E'},
    // { "ont-hp-prof", 1, NULL, 0},

    // short options
    { "ref-idx", 1, NULL, 'r' },
    { "out-vcf", 1, NULL, 'o'},
    { "out-type", 1, NULL, 'O'},
    { "min-sv-len", 1, NULL, 'l'},
    { "out-sam", 1, NULL, 'S' },
    { "out-bam", 1, NULL, 'b' },
    { "out-cram", 1, NULL, 'C' },
    { "no-vcf-header", 0, NULL, 'H'},
    { "somatic", 0, NULL, 's'},
    { "methylation", 0, NULL, 'm'},
    { "gap-aln", 1, NULL, 'g'},
    // { "out-bam", 1, NULL, 'b' },
    { "sample-name", 1, NULL, 'n'},
    { "trans-elem", 1, NULL, 'T'},

    { "min-cov", 1, NULL, 'c' },
    { "alt-cov", 1, NULL, 'd' },
    { "alt-ratio", 1, NULL, 'a' },
    { "min-mapq", 1, NULL, 'M'}, // map quality
    { "min-bq", 1, NULL, 'B'}, // base quality
    // { "max-ploidy", 1, NULL, 'p' },

    { "max-xgap", 1, NULL, 'x' },
    { "win-size", 1, NULL, 'w'},
    { "noisy-rat", 1, NULL, 'j' },
    { "noisy-flank", 1, NULL, 'f' },
    { "merge-dis", 1, NULL, 'D'},
    { "end-clip", 1, NULL, 'c' },
    { "clip-flank", 1, NULL, 'F' },
    { "hap-read", 1, NULL, 'p' },
    { "full-read", 1, NULL, 'f' },
    { "no-re-aln", 1, NULL, 'N'},

    { "threads", 1, NULL, 't' },
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    { "verbose", 1, NULL, 'V' },
    { 0, 0, 0, 0}
};

int get_num_processors() {
#ifdef __linux__
    return get_nprocs();
#elif __APPLE__
    int nm[2] = {CTL_HW, HW_AVAILCPU};
    int ncpu;
    size_t len = sizeof(ncpu);

    if (sysctl(nm, 2, &ncpu, &len, NULL, 0) == -1 || ncpu < 1) {
        nm[1] = HW_NCPU;
        sysctl(nm, 2, &ncpu, &len, NULL, 0);
        if (ncpu < 1) ncpu = 1;
    }
    return ncpu;
#else
    #error "Unsupported platform"
#endif
}

void set_hifi_opt(call_var_opt_t *opt) {
    opt->is_pb_hifi = 1; opt->is_ont = 0;
    // HiFi: 10/200
    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS; //10
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN; // 200bp slide window
}

void set_ont_opt(call_var_opt_t *opt) {
    opt->is_ont = 1; opt->is_pb_hifi = 0;
    opt->strand_bias_pval = LONGCALLD_STRAND_BIAS_PVAL_ONT;
    // ONT: 10/50
    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS; // 10
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_ONT_SLIDE_WIN; // 50bp slide window
}

// extreme low allele frequency < 5% or even 1%
// void set_mosaic_opt(call_var_opt_t *opt) {
//     opt->somatic_beta_alpha = 1; opt->somatic_beta_beta = 10;
//     opt->min_somatic_log_beta_binom = -3.0; // log10(0.001)
// }

// // moderate low allele frequency, > 5%, (tumor purity: 40~100% -> HET var: 20%~50% allele frequency)
// void set_tumor_somatic_opt(call_var_opt_t *opt) {
//     opt->somatic_beta_alpha = 2; opt->somatic_beta_beta = 5;
//     opt->min_somatic_log_beta_binom = -4.0; // log10(0.0001) ???
// }

call_var_opt_t *call_var_init_para(void) {
    call_var_opt_t *opt = (call_var_opt_t*)_err_malloc(sizeof(call_var_opt_t));

    opt->sample_name = NULL;
    opt->ref_fa_fn = NULL; opt->ref_fa_fai_fn = NULL; opt->reg_bed_fn = NULL;

    opt->is_pb_hifi = 1; opt->is_ont = 0;
    opt->only_autosome = 0; opt->only_autosome_XY = 1; 
    opt->n_exc_tnames = 0; opt->exc_tnames = NULL;

    // opt->region_list = NULL; opt->region_is_file = 0;
    opt->max_ploid = LONGCALLD_DEF_PLOID;
    opt->min_mq = LONGCALLD_MIN_CAND_MQ;
    opt->min_bq = LONGCALLD_MIN_CAND_BQ;
    opt->min_dp = LONGCALLD_MIN_CAND_DP;
    opt->min_alt_dp = LONGCALLD_MIN_ALT_DP;
    opt->min_af = LONGCALLD_MIN_CAND_AF;
    opt->max_af = LONGCALLD_MAX_CAND_AF;
    // somatic/mosaic var
    opt->min_somatic_dis_to_var = LONGCALLD_MIN_SOMATIC_DIS_TO_VAR;
    opt->min_somatic_dis_to_homopolymer_indel_error = LONGCALLD_MIN_SOMATIC_DIS_TO_HP_INDEL_ERROR;
    opt->min_somatic_dis_to_seq_error = LONGCALLD_MIN_SOMATIC_DIS_TO_SEQ_ERROR;
    opt->min_somatic_fisher_pval = LONGCALLD_MIN_SOMATIC_FISHER_PVAL;
    opt->min_somatic_alt_dp = LONGCALLD_MIN_SOMATIC_ALT_DP;
    opt->min_somatic_hap_dp = LONGCALLD_MIN_SOMATIC_HAP_READS;
    opt->min_somatic_te_dp = LONGCALLD_MIN_SOMATIC_TE_ALT_DP;
    // opt->min_somatic_af = LONGCALLD_MIN_SOMATIC_AF;
    // opt->max_somatic_alt_af = LONGCALLD_MAX_SOMATIC_ALT_AF;
    // opt->min_somatic_log_beta_binom = LONGCALLD_MIN_SOMATIC_LOG_BETA_BINOM;
    // opt->somatic_beta_alpha = LONGCALLD_SOMATIC_BETA_ALPHA;
    // opt->somatic_beta_beta = LONGCALLD_SOMATIC_BETA_BETA;
    opt->somatic_win = LONGCALLD_SOMATIC_WIN;
    opt->somatic_win_max_vars = LONGCALLD_SOMATIC_WIN_MAX_VARS;

    opt->noisy_reg_merge_dis = LONGCALLD_NOISY_REG_MERGE_DIS;
    opt->noisy_reg_flank_len = LONGCALLD_NOISY_REG_FLANK_LEN;

    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS;
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN;

    opt->end_clip_reg = LONGCALLD_NOISY_END_CLIP;
    opt->end_clip_reg_flank_win = LONGCALLD_NOISY_END_CLIP_WIN;

    opt->max_var_ratio_per_read = LONGCALLD_MAX_VAR_RATIO_PER_READ;
    opt->max_noisy_reg_reads = LONGCALLD_MAX_NOISY_REG_READS;
    opt->max_noisy_reg_len  = LONGCALLD_MAX_NOISY_REG_LEN;
    opt->min_noisy_reg_reads = LONGCALLD_NOISY_REG_READS;
    // opt->min_noisy_reg_ratio = LONGCALLD_NOISY_REG_RATIO;
    opt->max_noisy_frac_per_read = LONGCALLD_MAX_NOISY_FRAC_PER_READ;
    opt->min_hap_full_reads = LONGCALLD_MIN_HAP_FULL_READS;
    opt->min_hap_reads = LONGCALLD_MIN_HAP_READS;
    // opt->min_no_hap_full_reads = LONGCALLD_MIN_NO_HAP_FULL_READS;

    opt->match = LONGCALLD_MATCH_SCORE;
    opt->mismatch = LONGCALLD_MISMATCH_SCORE;
    opt->gap_open1 = LONGCALLD_GAP_OPEN1_SCORE;
    opt->gap_ext1 = LONGCALLD_GAP_EXT1_SCORE;
    opt->gap_open2 = LONGCALLD_GAP_OPEN2_SCORE;
    opt->gap_ext2 = LONGCALLD_GAP_EXT2_SCORE;
    opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
    opt->disable_read_realign = 1; // XXX right now, always disable read realignment

    opt->pl_threads = MIN_OF_TWO(CALL_VAR_PL_THREAD_N, get_num_processors());
    opt->n_threads = MIN_OF_TWO(CALL_VAR_THREAD_N, get_num_processors());

    opt->min_sv_len = LONGCALLD_MIN_SV_LEN;
    opt->min_tsd_len = LONGCALLD_MIN_TSD_LEN; opt->max_tsd_len = LONGCALLD_MAX_TSD_LEN;
    opt->min_polya_len = LONGCALLD_MIN_POLYA_LEN; opt->min_polya_ratio = LONGCALLD_MIN_POLYA_RATIO;

    initialize_lgamma_cache(opt);
    opt->te_seq_fn = NULL; opt->n_te_seqs = 0; opt->te_kmer_len = 15; opt->te_for_h = NULL; opt->te_rev_h = NULL;

    opt->p_error = 0.001; opt->log_p = -3.0; opt->log_1p = log10(1-opt->p_error); opt->log_2 = 0.301023;
    opt->max_gq = 60; opt->max_qual = 60;
    opt->out_vcf = NULL; opt->vcf_hdr = NULL; opt->out_vcf_fn = NULL; opt->out_vcf_type = 'v'; opt->no_vcf_header = 0; opt->out_amb_base = 0;
    opt->out_bam = NULL; opt->out_is_cram = 0; opt->refine_bam = 0;
    opt->out_somatic = 0; opt->out_methylation = 0;
    // opt->verbose = 0;
    return opt;
}

void call_var_free_para(call_var_opt_t *opt) {
    if (opt->ref_fa_fn) free(opt->ref_fa_fn);
    if (opt->in_bam_fn) free(opt->in_bam_fn);
    if (opt->sample_name) free(opt->sample_name);
    // if (opt->region_list) free(opt->region_list);
    if (opt->n_exc_tnames > 0) {
        for (int i = 0; i < opt->n_exc_tnames; ++i) free(opt->exc_tnames[i]);
        free(opt->exc_tnames);
    }
    if (opt->out_vcf_fn != NULL) free(opt->out_vcf_fn);
    if (opt->te_seq_fn != NULL) {
        free(opt->te_seq_fn);
        for (int i = 0; i < opt->n_te_seqs; ++i) {
            free(opt->te_seq_names[i]);
            kh_destroy(kmer32, opt->te_for_h[i]);
            kh_destroy(kmer32, opt->te_rev_h[i]);
        }
        free(opt->te_seq_names); free(opt->te_for_h); free(opt->te_rev_h);
    }
    free(opt);
}

void reg_chunks_free(reg_chunks_t *reg_chunks, int m_reg_chunks) {
    for (int i = 0; i < m_reg_chunks; ++i) {
        if (reg_chunks[i].reg_tids) free(reg_chunks[i].reg_tids);
        if (reg_chunks[i].reg_begs) free(reg_chunks[i].reg_begs);
        if (reg_chunks[i].reg_ends) free(reg_chunks[i].reg_ends);
    }
    free(reg_chunks);
}

void call_var_io_aux_free(call_var_io_aux_t *aux, int n) {
    for (int i = 0; i < n; ++i) {
        if (aux[i].fai) fai_destroy(aux[i].fai);
        if (aux[i].bam) {
            if (aux[i].idx) hts_idx_destroy(aux[i].idx);
            bam_hdr_destroy(aux[i].header);
            sam_close(aux[i].bam); 
        }
    }
    free(aux);
}
void call_var_free_pl(call_var_pl_t pl) {
    call_var_io_aux_free(pl.io_aux, pl.n_threads);
    reg_chunks_free(pl.reg_chunks, pl.m_reg_chunks);
}

void var_free(var_t *v) {
    if (v->vars) {
        for (int i = 0; i < v->n; ++i) {
            if (v->vars[i].ref_bases) free(v->vars[i].ref_bases);
            if (v->vars[i].alt_len) free(v->vars[i].alt_len);
            if (v->vars[i].alt_bases) {
                for (int j = 0; j < v->vars[i].n_alt_allele; ++j) free(v->vars[i].alt_bases[j]);
                free(v->vars[i].alt_bases);
            }
            if (v->vars[i].tsd_len > 0) {
                free(v->vars[i].tsd_seq);
            }
        }
        free(v->vars);
    }
}

void var1_free(var1_t *v) {
    if (v->ref_bases) free(v->ref_bases);
    if (v->alt_len) free(v->alt_len);
    if (v->alt_bases) {
        for (int j = 0; j < v->n_alt_allele; ++j) free(v->alt_bases[j]);
        free(v->alt_bases);
    }
    if (v->tsd_len > 0) {
        free(v->tsd_seq);
    }
}

static void collect_bam_call_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks + ii;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d) ... \n", __func__, tid, ii, step->n_chunks);
    // if (collect_ref_seq_bam_main(step->pl, step->pl->io_aux+tid, step->pl->reg_chunk_i, ii, c) > 0) {
    collect_ref_seq_bam_main(step->pl, step->pl->io_aux+tid, step->pl->reg_chunk_i, ii, c);
    collect_var_main(step->pl, c);
    // }
    bam_chunk_mid_free(c);
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d) ... done\n", __func__, tid, ii, step->n_chunks);
}

// merge variants from two chunks
static void make_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks+ii; var_t *var = step->vars+ii;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] thread-id: %d, chunk: %ld (%d), n_reads: %d\n", __func__, tid, ii, step->n_chunks, c->n_reads);
    make_var_main(step, c, var, ii);
}

static void reg_chunks_realloc(call_var_pl_t *pl) {
    if (pl->n_reg_chunks >= pl->m_reg_chunks) {
        pl->m_reg_chunks *= 2;
        pl->reg_chunks = (reg_chunks_t*)realloc(pl->reg_chunks, pl->m_reg_chunks * sizeof(reg_chunks_t));
        for (int k = pl->n_reg_chunks; k < pl->m_reg_chunks; k++) {
            pl->reg_chunks[k].n_regions = 0;
            pl->reg_chunks[k].m_regions = 512;
            pl->reg_chunks[k].reg_tids = (int*)malloc(pl->reg_chunks[k].m_regions * sizeof(int));
            pl->reg_chunks[k].reg_begs = (hts_pos_t*)malloc(pl->reg_chunks[k].m_regions * sizeof(hts_pos_t));
            pl->reg_chunks[k].reg_ends = (hts_pos_t*)malloc(pl->reg_chunks[k].m_regions * sizeof(hts_pos_t));
        }
    }
}

void set_exclude_ctg(call_var_opt_t *opt, char *ctg) {
    opt->exc_tnames = (char**)realloc(opt->exc_tnames, (opt->n_exc_tnames+1) * sizeof(char*));
    opt->exc_tnames[opt->n_exc_tnames] = strdup(ctg);
    opt->n_exc_tnames++;
}

// Enum for chromosome types
typedef enum {
    AUTOSOME,        // Chromosomes 1-22
    SEX_CHROMOSOME,  // X or Y
    OTHER            // Mitochondrial (MT/M) or unplaced/unusual chromosomes
} ChromosomeType;

ChromosomeType classify_chromosome(const char *chr) {
    char chr_copy[1000];  // Buffer to hold the extracted chromosome name
    strncpy(chr_copy, chr, sizeof(chr_copy) - 1);
    chr_copy[sizeof(chr_copy) - 1] = '\0';

    // Find first occurrence of ':', indicating a genomic region
    char *colon_pos = strchr(chr_copy, ':');
    if (colon_pos != NULL) {
        *colon_pos = '\0';  // Truncate at ':'
    }

    // Remove "chr" prefix if present
    if (strncmp(chr_copy, "chr", 3) == 0) {
        memmove(chr_copy, chr_copy + 3, strlen(chr_copy) - 2);  // Shift left
    }

    // Check for sex chromosomes
    if (strcmp(chr_copy, "X") == 0 || strcmp(chr_copy, "Y") == 0) {
        return SEX_CHROMOSOME;
    }

    // Check for other chromosomes (e.g., Mitochondrial DNA)
    if (strcmp(chr_copy, "MT") == 0 || strcmp(chr_copy, "M") == 0) {
        return OTHER;
    }

    // Convert string to integer and check for autosome range
    char *endptr;
    long num = strtol(chr_copy, &endptr, 10);

    // if (*endptr == '\0' && num >= 1 && num <= 22) {
    if (*endptr == '\0' && num >= 1) {
        return AUTOSOME;
    }
    return OTHER; // Anything else falls under "OTHER"
}

static void reg_chunks_region_realloc(reg_chunks_t *reg_chunks) {
    if (reg_chunks->n_regions >= reg_chunks->m_regions) {
        reg_chunks->m_regions *= 2;
        reg_chunks->reg_tids = (int*)realloc(reg_chunks->reg_tids, reg_chunks->m_regions * sizeof(int));
        reg_chunks->reg_begs = (hts_pos_t*)realloc(reg_chunks->reg_begs, reg_chunks->m_regions * sizeof(hts_pos_t));
        reg_chunks->reg_ends = (hts_pos_t*)realloc(reg_chunks->reg_ends, reg_chunks->m_regions * sizeof(hts_pos_t));
    }
}

static int exclude_target_ctg(call_var_opt_t *opt, char *tname) {
    for (int i = 0; i < opt->n_exc_tnames; ++i) {
        if (strcmp(tname, opt->exc_tnames[i]) == 0) return 1;
    }
    return 0;
}

static int skip_target_region(call_var_opt_t *opt, char *tname) {
    ChromosomeType t = classify_chromosome(tname);
    if (opt->only_autosome && t != AUTOSOME) return 1;
    if (opt->only_autosome_XY && t != AUTOSOME && t != SEX_CHROMOSOME) return 1;
    if (exclude_target_ctg(opt, tname)) return 1;
    return 0;
}

static void collect_regions_from_region_list(call_var_opt_t *opt, call_var_pl_t *pl, hts_itr_t *iter, int n_regions, char **regions) {
    int last_tid = -1; bam_hdr_t *hdr = pl->io_aux[0].header;
    int min_reg_chunks_per_run = pl->min_reg_chunks_per_run, max_reg_len_per_chunk = pl->max_reg_len_per_chunk;
    for (int i = 0; i < iter->n_reg; ++i) {
        for (int j = 0; j < iter->reg_list[i].count; ++j) { // ACGT:1234, (beg, end]:(0,4]
            int tid = iter->reg_list[i].tid; hts_pos_t chr_len = hdr->target_len[tid]; char *tname = hdr->target_name[tid];
            if (skip_target_region(opt, tname)) continue;
            // fprintf(stderr, "Region: %s:%" PRId64 "-%" PRId64 "\n", hdr->target_name[tid], iter->reg_list[i].intervals[j].beg, iter->reg_list[i].intervals[j].end);
            if (last_tid != -1 && tid != last_tid && pl->reg_chunks[pl->n_reg_chunks].n_regions >= min_reg_chunks_per_run)
                pl->n_reg_chunks++;
            reg_chunks_realloc(pl);
            hts_pos_t reg_beg = MAX_OF_TWO(1, iter->reg_list[i].intervals[j].beg + 1); // 0-base
            hts_pos_t reg_end = MIN_OF_TWO(iter->reg_list[i].intervals[j].end, chr_len);
            hts_pos_t n_regions = (reg_end - reg_beg + max_reg_len_per_chunk) / max_reg_len_per_chunk;
            reg_chunks_t *reg_chunks = &pl->reg_chunks[pl->n_reg_chunks];
            for (int reg_i = 0; reg_i < n_regions; ++ reg_i) {
                hts_pos_t beg = (hts_pos_t) reg_i * max_reg_len_per_chunk + reg_beg;
                hts_pos_t end = MIN_OF_TWO((hts_pos_t) (reg_i + 1) * max_reg_len_per_chunk + reg_beg-1, reg_end);
                // add region
                reg_chunks_region_realloc(reg_chunks);
                reg_chunks->reg_tids[reg_chunks->n_regions] = tid;
                reg_chunks->reg_begs[reg_chunks->n_regions] = beg; // [beg, end]
                reg_chunks->reg_ends[reg_chunks->n_regions] = end;
                reg_chunks->n_regions++;
            }
            last_tid = tid;
        }
    }
    if (pl->reg_chunks[pl->n_reg_chunks].n_regions > 0) pl->n_reg_chunks++;
    hts_itr_destroy(iter);
    // Print the region_chunks
    if (LONGCALLD_VERBOSE >= 1) {
        for (int i = 0; i < pl->n_reg_chunks; i++) {
            fprintf(stderr, "Chunk %d:\n", i);
            for (int j = 0; j < pl->reg_chunks[i].n_regions; j++) {
                fprintf(stderr, "\tRegion %d %s:%" PRIi64 "-%" PRIi64 "\n", j, pl->io_aux[0].header->target_name[pl->reg_chunks[i].reg_tids[j]], pl->reg_chunks[i].reg_begs[j], pl->reg_chunks[i].reg_ends[j]);
            }
        }
    }
}

static void collect_regions_from_bed_file(call_var_opt_t *opt, call_var_pl_t *pl) {
    FILE *fp = fopen(opt->reg_bed_fn, "r");
    if (fp == NULL) {
        _err_error("Failed to open region bed file: %s\n", opt->reg_bed_fn);
        return;
    }
    char line[1024]; int needs_sort = 0;
    // collect regions
    bam_hdr_t *hdr = pl->io_aux[0].header;
    int n_regions = 0, m_regions = 1024;
    char **regions = (char**)malloc(m_regions * sizeof(char*));
    for (int i = 0; i < m_regions; ++i) regions[i] = (char*)malloc(1024);
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue; // skip header/comment
        // exclude \n in line if necessary
        if (line[strlen(line)-1] == '\n') line[strlen(line)-1] = '\0';
        char *token;
        token = strtok(line, "\t");
        if (token == NULL) continue;
        char *tname = token; int tid = bam_name2id(hdr, tname);
        if (tid == -1) {
            _err_warning("Skip unknonw region: %s\n", line);
            continue;
        }
        if (skip_target_region(opt, tname)) continue;
        hts_pos_t beg = 1, end = hdr->target_len[tid];
        token = strtok(NULL, "\t");
        if (token != NULL) {
            beg = atoi(token)+1;
            token = strtok(NULL, "\t");
            if (token != NULL) end = atoi(token);
        }
        if (beg > end || beg <= 0 || end <= 0) continue;
        if (n_regions == m_regions) {
            n_regions *= 2;
            regions = (char**)realloc(regions, m_regions * sizeof(char*));
            for (int i = n_regions; i < m_regions; ++i) regions[i] = (char*)malloc(1024);
        }
        sprintf(regions[n_regions], "%s:%" PRIi64 "-%" PRIi64 "", tname, beg, end);
        n_regions++;
    }
    hts_idx_t *idx = pl->io_aux[0].idx;
    hts_itr_t *iter = sam_itr_regarray(idx, hdr, regions, n_regions);
    if (iter != NULL) collect_regions_from_region_list(opt, pl, iter, n_regions, regions);
}

static void collect_regions(call_var_pl_t *pl, call_var_opt_t *opt, int n_regions, char **regions) {
    int min_reg_chunks_per_run = pl->min_reg_chunks_per_run;
    int max_reg_len_per_chunk = pl->max_reg_len_per_chunk;
    pl->n_reg_chunks = 0; pl->m_reg_chunks = 100;
    pl->reg_chunks = (reg_chunks_t*)malloc(pl->m_reg_chunks * sizeof(reg_chunks_t));
    for (int i = 0; i < pl->m_reg_chunks; i++) {
        pl->reg_chunks[i].n_regions = 0;
        pl->reg_chunks[i].m_regions = 512;
        pl->reg_chunks[i].reg_tids = (int*)malloc(pl->reg_chunks[i].m_regions * sizeof(int));
        pl->reg_chunks[i].reg_begs = (hts_pos_t*)malloc(pl->reg_chunks[i].m_regions * sizeof(hts_pos_t));
        pl->reg_chunks[i].reg_ends = (hts_pos_t*)malloc(pl->reg_chunks[i].m_regions * sizeof(hts_pos_t));
    }
    // check if region(s) is provided
    if (n_regions > 0) {
        opt->only_autosome = 0, opt->only_autosome_XY=0;
        if (opt->reg_bed_fn != NULL) {
            _err_error("Both region(s) and region bed file are provided. Only region(s) will be used.\n");
        }
        hts_idx_t *idx = pl->io_aux[0].idx; bam_hdr_t *hdr = pl->io_aux[0].header;
        hts_itr_t *iter = sam_itr_regarray(idx, hdr, regions, n_regions);
        if (iter == NULL) {
            _err_error("Failed to parse provided regions: ");
            for (int i = 0; i < n_regions; ++i) fprintf(stderr, "%s ", regions[i]);
            fprintf(stderr, "\n");
            _err_error("Will process the whole alignment file.\n");
        } else {
            return collect_regions_from_region_list(opt, pl, iter, n_regions, regions);
        }
    } else if (opt->reg_bed_fn != NULL) { // check if region bed file is provided
        opt->only_autosome = 0, opt->only_autosome_XY=0;
        collect_regions_from_bed_file(opt, pl);
        if (pl->n_reg_chunks == 0) {
            _err_error("Failed to parse provided region bed file: %s\n", opt->reg_bed_fn);
            _err_error("Will process the whole alignment file.\n");
        } else return;
    }
    // whole genome, plit chromosomes into region_chunks
    // collect_regions_from_whole_genome(opt, pl, min_reg_chunks_per_run, max_reg_len_per_chunk);
    for (int i = 0; i < pl->io_aux[0].header->n_targets; i++) {
        reg_chunks_realloc(pl);
        // all regions in a chromosome will be in one reg_chunk
        char *tname = pl->io_aux[0].header->target_name[i];
        if (skip_target_region(opt, tname)) continue;
        int tid = i; hts_pos_t chr_len = pl->io_aux[0].header->target_len[i];
        int n_regions_chr = (chr_len + max_reg_len_per_chunk - 1) / max_reg_len_per_chunk;
        reg_chunks_t *reg_chunks = &pl->reg_chunks[pl->n_reg_chunks];
        for (int reg_i = 0; reg_i < n_regions_chr; reg_i++) {
            hts_pos_t beg = (hts_pos_t) reg_i * max_reg_len_per_chunk + 1;
            hts_pos_t end = MIN_OF_TWO((hts_pos_t) (reg_i + 1) * max_reg_len_per_chunk, chr_len);
            // add region
            reg_chunks_region_realloc(reg_chunks);
            reg_chunks->reg_tids[reg_chunks->n_regions] = tid;
            reg_chunks->reg_begs[reg_chunks->n_regions] = beg;
            reg_chunks->reg_ends[reg_chunks->n_regions] = end;
            reg_chunks->n_regions++;
        }
        if (reg_chunks->n_regions >= min_reg_chunks_per_run) {
            pl->n_reg_chunks++;
        }
    }
    if (pl->reg_chunks[pl->n_reg_chunks].n_regions > 0) pl->n_reg_chunks++;

    // Print the region_chunks
    if (LONGCALLD_VERBOSE >= 2) {
        for (int i = 0; i < pl->n_reg_chunks; i++) {
            fprintf(stderr, "Chunk %d:\n", i);
            for (int j = 0; j < pl->reg_chunks[i].n_regions; j++) {
                fprintf(stderr, "\tRegion %d %s:%" PRIi64 "-%" PRIi64 "\n", j, pl->io_aux[0].header->target_name[pl->reg_chunks[i].reg_tids[j]], pl->reg_chunks[i].reg_begs[j], pl->reg_chunks[i].reg_ends[j]);
            }
        }
    }
}

// only works with sorted index BAM/CRAM
static void call_var_pl_open_fa_bam(call_var_opt_t *opt, call_var_pl_t *pl, char **regions, int n_regions) {
    // input BAM file
    if (strcmp(opt->in_bam_fn, "-") == 0) { _err_error_exit("Input from pipe/stdin is not supported\n"); }
    else { _err_info("Opening alignment file: %s\n", opt->in_bam_fn); }
    pl->io_aux = (call_var_io_aux_t*)malloc(pl->n_threads * sizeof(call_var_io_aux_t));
    // multi-threading
    for (int i = 0; i < pl->n_threads; ++i){
        pl->io_aux[i].bam = sam_open(opt->in_bam_fn, "r"); if (pl->io_aux[i].bam == NULL) _err_error_exit("Failed to open alignment file \'%s\'\n", opt->in_bam_fn);
        const htsFormat *fmt = hts_get_format(pl->io_aux[i].bam); if (!fmt) _err_error_exit("Failed to get format of alignment file \'%s\'\n", opt->in_bam_fn);
        if (fmt->format != bam && fmt->format != cram) {
            sam_close(pl->io_aux[i].bam); _err_error_exit("Input alignment file must be BAM or CRAM format.\n");
        }
        pl->io_aux[i].hts_fmt = fmt->format;
        pl->io_aux[i].header = sam_hdr_read(pl->io_aux[i].bam);
        if (pl->io_aux[i].header == NULL) {
            sam_close(pl->io_aux[i].bam);
            _err_error_exit("Failed to read alignment file header \'%s\'\n", opt->in_bam_fn);
        }
        pl->io_aux[i].fai = fai_load3(opt->ref_fa_fn, opt->ref_fa_fai_fn, NULL, FAI_CREATE);
        if (pl->io_aux[i].fai == NULL) {
            sam_hdr_destroy(pl->io_aux[i].header); sam_close(pl->io_aux[i].bam);
            _err_error_exit("Failed to load/build reference fasta index: %s\n", opt->ref_fa_fn);
        }
        if (fmt->format == cram) {
            if (hts_set_fai_filename(pl->io_aux[i].bam, opt->ref_fa_fn) != 0) {
                fai_destroy(pl->io_aux[i].fai); sam_hdr_destroy(pl->io_aux[i].header); sam_close(pl->io_aux[i].bam);
                _err_error_exit("Failed to set reference file for CRAM decoding: %s %s\n", opt->in_bam_fn, opt->ref_fa_fn);
            }
        }
        pl->io_aux[i].idx = sam_index_load(pl->io_aux[i].bam, opt->in_bam_fn);
        if (pl->io_aux[i].idx == NULL) { // attempt to create index if not exist
            _err_warning("Index not found for \'%s\', creating now ...\n", opt->in_bam_fn);
            if (sam_index_build(opt->in_bam_fn, 0) < 0) {
                fai_destroy(pl->io_aux[i].fai); sam_hdr_destroy(pl->io_aux[i].header); sam_close(pl->io_aux[i].bam);
                _err_error_exit("Failed to build index for \'%s\'\n", opt->in_bam_fn);
            } else {
                pl->io_aux[i].idx = sam_index_load(pl->io_aux[i].bam, opt->in_bam_fn);
                if (pl->io_aux[i].idx == NULL) {
                    fai_destroy(pl->io_aux[i].fai); sam_hdr_destroy(pl->io_aux[i].header); sam_close(pl->io_aux[i].bam);
                    _err_error_exit("Failed to load index for \'%s\'\n", opt->in_bam_fn);
                }
            }
        }
    }
    // collect sample name (SM) from BAM header
    if (opt->sample_name == NULL) opt->sample_name = extract_sample_name_from_bam_header(pl->io_aux[0].header); // extract sample names from BAM header
    if (opt->sample_name == NULL) opt->sample_name = strdup(opt->in_bam_fn);
    // collect regions
    collect_regions(pl, opt, n_regions, regions);
}

static void call_var_pl_write_bam_header(call_var_opt_t *opt, bam_hdr_t *header) {
    if (opt->out_is_cram) {
        if (hts_set_fai_filename(opt->out_bam, opt->ref_fa_fn) != 0) _err_error_exit("Failed to set reference file for output CRAM encoding: %s\n", opt->ref_fa_fn);
    }
    if (sam_hdr_add_pg(header, PROG, "VN", LONGCALLD_VERSION, "CL", CMD, NULL) < 0) _err_error_exit("Fail to add PG line to bam header.\n");
    if (sam_hdr_write(opt->out_bam, header) < 0) _err_error_exit("Failed to write BAM header.\n");
}

// work with sorted SAM/BAM/CRAM
static void *call_var_worker_pipeline(void *shared, int step, void *in) { // kt_pipeline() callback
    call_var_pl_t *pl = (call_var_pl_t*)shared;
    if (step == 0) { // step 0: read bam records into BAM chunks, call variants
        if (pl->reg_chunk_i >= pl->n_reg_chunks) return 0;
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 0: load BAM & call variants.\n");
        // _err_info("Processing contig(s): %d/%d ...\n", pl->reg_chunk_i+1, pl->n_reg_chunks);
        call_var_step_t *s = calloc(1, sizeof(call_var_step_t));
        s->pl = pl;
        s->n_chunks = pl->reg_chunks[pl->reg_chunk_i].n_regions;
        s->chunks = calloc(s->n_chunks, sizeof(bam_chunk_t));
        // _err_info("Processing ctg-id: %d, n_chunks: %d (%d/%d)\n", pl->reg_chunks[pl->reg_chunk_i].reg_tids[0], s->n_chunks, pl->reg_chunk_i+1, pl->n_reg_chunks);
        kt_for(pl->n_threads, collect_bam_call_var_worker_for, s, s->n_chunks);
        pl->reg_chunk_i++;
        return s;
    } else if (step == 1) { // step 1: stitch the results
                            //         merge variants from different chunks
                            //         1) merge overlapping variants
                            //         2) update phase set and haplotype (flip if needed)
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 2: stitch & make variants\n");
        if (((call_var_step_t*)in)->n_chunks > 0) {
            ((call_var_step_t*)in)->vars = calloc(((call_var_step_t*)in)->n_chunks, sizeof(var_t));
            stitch_var_main((call_var_step_t*)in, ((call_var_step_t*)in)->chunks);
            kt_for(pl->n_threads, make_var_worker_for, in, ((call_var_step_t*)in)->n_chunks);
        }
        int64_t n_processed_reads = 0;
        for (int i = 0; i < ((call_var_step_t*)in)->n_chunks; ++i) {
            n_processed_reads += (((call_var_step_t*)in)->chunks[i].n_reads - ((call_var_step_t*)in)->chunks[i].n_up_ovlp_reads);
        }
        if (n_processed_reads > 0) _err_info("Processed %" PRIi64 " reads, %d/%d chunks\n", n_processed_reads, pl->reg_chunk_i, pl->n_reg_chunks);
        return in;
    } else if (step == 2) { // step 3: write the buffer to output
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 3: output variants (& phased bam)\n");
        call_var_step_t *s = (call_var_step_t*)in;
        int n_out_vars = 0, n_out_reads = 0;
        for (int i = 0; i < s->n_chunks; ++i) {
            bam_chunk_t *c = s->chunks + i;
            n_out_vars += write_var_to_vcf(s->vars+i, pl->opt, c->tname);
            if (pl->opt->out_bam != NULL) n_out_reads += write_read_to_bam(c, pl->opt);
            var_free(s->vars + i);  // free output
        }
        if (n_out_vars > 0) _err_info("Output %d variants to VCF\n", n_out_vars);
        if (n_out_reads > 0) {
            if (pl->opt->out_is_cram) _err_info("Output %d reads to CRAM\n", n_out_reads);
            else _err_info("Output %d reads to BAM\n", n_out_reads);
        }
        bam_chunks_post_free(s->chunks, s->n_chunks); // free input
        free(s->vars); free(s);
    }
    return 0;
}

static void call_var_usage(void) {//main usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s (%s)\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", LONGCALLD_VERSION, CONTACT); 

    fprintf(stderr, "Usage: %s call [options] <ref.fa> <input.bam/cram> [region ...] > phased.vcf\n", PROG);
    fprintf(stderr, "       Note: \'ref.fa\' should contain the same contig/chromosome names as \'input.bam/cram\'\n");
    fprintf(stderr, "             \'input.bam/cram' should be sorted and indexed\n");
    // fprintf(stderr, "                  \'ref.fa\' is not needed if 1) BAM contains CIGARS with =/X, or 2) BAM contain MD tags\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Input:\n");
    fprintf(stderr, "    --hifi                HiFi reads [Default]\n");
    fprintf(stderr, "    --ont                 ONT reads [False]\n");
    fprintf(stderr, "    --region-file   FILE  region bed file [NULL]\n");
    fprintf(stderr, "                          each line is a region, e.g., chr1         (whole chromosome)\n");
    fprintf(stderr, "                                                       chr1 100     (chr1:100-END)\n");
    fprintf(stderr, "                                                       chr1 100 200 (chr1:100-200)\n");
    fprintf(stderr, "    --autosome-XY         only call variants on autosomes and sex chromosomes [Default]\n");
    fprintf(stderr, "                          e.g., human: chr{1-22,XY} or 1-22,XY\n");
    fprintf(stderr, "    --all-ctg             call variants on all chromosomes/contigs [False]\n");
    fprintf(stderr, "    --autosome            only call variants on autosomes [False]\n");
    fprintf(stderr, "                          e.g., human: chr{1-22} or 1-22\n");
    fprintf(stderr, "    -E --exclude-ctg STR  exclude contig/chromosome [NULL]\n");
    fprintf(stderr, "                          can be used multiple times, e.g., -E hs37d5 -E chrM\n");
    fprintf(stderr, "    -n --sample-name STR  sample name. 'RG/SM' in BAM header will be used if not provided [NULL]\n");
    fprintf(stderr, "                          BAM file name will be used if not provided and no 'RG/SM' in BAM header\n");
    fprintf(stderr, "    -r --ref-idx    FILE  .fai index file for reference genome FASTA file, detect automaticaly if not provided [NULL]\n");
    fprintf(stderr, "    -T --trans-elem FILE  transposable element sequence [NULL]\n");
    fprintf(stderr, "  Output:\n");
    fprintf(stderr, "    -o --out-vcf    FILE  output phased VCF file [stdout]\n");
    fprintf(stderr, "    -O --out-type    STR  v/z: un/compressed VCF [v]\n");
    fprintf(stderr, "    -l --min-sv-len  INT  min. length to be considered as SV [%d]\n", LONGCALLD_MIN_SV_LEN);
    fprintf(stderr, "                          SV-related information will be added to INFO field, e.g., SVLEN/SVTYPE/TSD\n");
    fprintf(stderr, "    -H --no-vcf-header    do NOT output VCF header [False]\n");
    fprintf(stderr, "    -s --mosaic/somatic   output low allele-frequency mosaic/somatic variants [False]\n");
    // fprintf(stderr, "    -m --methylation    output methylation site information [False]\n");
    // fprintf(stderr, "                          if present, MM/ML tags will be used to calculate methylation level\n");
    fprintf(stderr, "       --amb-base         output variant with ambiguous base (N) [False]\n");
    fprintf(stderr, "    -b --out-bam    FILE  output phased BAM file [NULL]\n");
    fprintf(stderr, "    -S --out-sam    FILE  output phased SAM file [NULL]\n");
    fprintf(stderr, "    -C --out-cram   FILE  output phased CRAM file [NULL]\n");
    fprintf(stderr, "    --refine-aln          refine alignment in output SAM/BAM/CRAM [False]\n");
    fprintf(stderr, "                          Note: output SAM/BAM/CRAM may be unsorted when --refine-aln is set\n");
    // fprintf(stderr, "    -g --gap-aln     STR  put gap on the \'left\' or \'right\' side in alignment [left/l]\n");
    // fprintf(stderr, "                          \'left\':  ATTTG\n");
    // fprintf(stderr, "                                   | |||\n");
    // fprintf(stderr, "                                   A-TTG\n");
    // fprintf(stderr, "                          \'right\': ATTTG\n");
    // fprintf(stderr, "                                   ||| |\n");
    // fprintf(stderr, "                                   ATT-G\n");
    // fprintf(stderr, "\n");
    fprintf(stderr, "  Variant calling:\n");
    fprintf(stderr, "    -c --min-cov     INT  min. total read coverage for candidate variant [%d]\n", LONGCALLD_MIN_CAND_DP);
    fprintf(stderr, "    -d --alt-cov     INT  min. alt. read coverage for candidate variant [%d]\n", LONGCALLD_MIN_ALT_DP);
    fprintf(stderr, "    -a --alt-ratio FLOAT  min. alt. read ratio for candidate variant [%.2f]\n", LONGCALLD_MIN_CAND_AF);
    fprintf(stderr, "    -M --min-mapq    INT  min. mapping quality for long-read alignment to be used [%d]\n", LONGCALLD_MIN_CAND_MQ);
    // fprintf(stderr, "    -B --min-bq      INT  filter out base with base quality < -B/--min-bq [%d]\n", LONGCALLD_MIN_CAND_BQ);
    // fprintf(stderr, "    -p --max-ploidy  INT  max. ploidy [%d]\n", LONGCALLD_DEF_PLOID);
    fprintf(stderr, "  Low allele-frequency mosaic/somatic variant calling: (effective when -s/--mosaic/--somatic is used)\n");
    // fprintf(stderr, "    --alpha          INT  alpha value of beta-binomial distribution [%d]\n", LONGCALLD_SOMATIC_BETA_ALPHA);
    // fprintf(stderr, "    --beta           INT  beta value of beta-binomial distribution [%d]\n", LONGCALLD_SOMATIC_BETA_BETA);
    fprintf(stderr, "    --som-alt        INT  min. alt. read coverage for somatic variant [%d]\n", LONGCALLD_MIN_SOMATIC_ALT_DP);
    fprintf(stderr, "    --som-mei-alt    INT  min. alt. read coverage for somatic MEI variant [%d]\n", LONGCALLD_MIN_SOMATIC_TE_ALT_DP);
    fprintf(stderr, "    --max-somvar INT,INT  max. number of mosaic/somatic variants allowed in a window (m,w) [%d,%d]\n", LONGCALLD_SOMATIC_WIN_MAX_VARS, LONGCALLD_SOMATIC_WIN);
    // fprintf(stderr, "  Variant calling in noisy regions:\n");
    // fprintf(stderr, "    -x --max-xgap    INT  max. number of allowed substitutions/gap-bases in a sliding window(-w/--win-size) [%d]\n", LONGCALLD_NOISY_REG_MAX_XGAPS);
    // fprintf(stderr, "                          window with more than -x subs/gap-bases will be considered as noisy region\n");
    // fprintf(stderr, "    -w --win-size    INT  window size for searching noisy region [%d]\n", LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN);
    // fprintf(stderr, "    -D --merge-dis   INT  default max. distance to merge two potential noisy SV regions [%d]\n", LONGCALLD_NOISY_REG_MERGE_DIS);
    // fprintf(stderr, "    -j --noisy-rat FLOAT  min. ratio of noisy reads in a window to call a noisy region [%.2f]\n", LONGCALLD_NOISY_REG_RATIO);
    // fprintf(stderr, "    -f --noisy-flank INT  flanking mask window size for noisy region [%d]\n", LONGCALLD_DENSE_FLANK_WIN);
    // fprintf(stderr, "    -c --end-clip    INT  max. number of clipping bases on both ends [%d]\n", LONGCALLD_NOISY_END_CLIP);
    // fprintf(stderr, "                          end-clipping region with more than -c bases will be considered as noisy clipping region\n");
    // fprintf(stderr, "    -F --clip-flank  INT  flanking mask window size for noisy clipping region [%d]\n", LONGCALLD_NOISY_END_CLIP_WIN);
    // fprintf(stderr, "\n");
    // fprintf(stderr, "    -p --hap-read    INT  when haplotype is available, min. number of full-spanning reads for each haplotype in noisy region to call a variant [%d]\n", LONGCALLD_MIN_HAP_FULL_READS);
    // fprintf(stderr, "    -f --full-read   INT  when haplotype is not available, min. number of full-spanning reads in noisy region to call a variant [%d]\n", LONGCALLD_MIN_NO_HAP_FULL_READS);

    // fprintf(stderr, "    -N --no-re-aln        disable read re-alignment\n");
    // fprintf(stderr, "\n");
    fprintf(stderr, "  General:\n");
    fprintf(stderr, "    -t --threads     INT  number of threads to use [%d]\n", MIN_OF_TWO(CALL_VAR_THREAD_N, get_num_processors()));
    fprintf(stderr, "    -h --help             print this help usage\n");
    fprintf(stderr, "    -v --version          print version number\n");
    fprintf(stderr, "    -V --verbose     INT  verbose level (0-2). 0: none, 1: information, 2: debug [0]\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    exit(1);
}

int call_var_main(int argc, char *argv[]) {
    // _err_cmd("%s\n", CMD);
    const char *opt_str = "r:T:E:o:O:ml:Hb:C:S:sc:d:M:B:a:n:x:w:D:j:L:f:p:g:Nt:hvV:";
    int c, op_idx; call_var_opt_t *opt = call_var_init_para(); char *s;
    double realtime0 = realtime();
    // first round of getopt_long() to parse preset options
    while ((c = getopt_long(argc, argv, opt_str, call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 0: 
                if (strcmp(call_var_opt[op_idx].name, "hifi") == 0) set_hifi_opt(opt);
                else if (strcmp(call_var_opt[op_idx].name, "ont") == 0) set_ont_opt(opt);
                break;
        }
    }
    optind = 1; // reset optind for next getopt_long() call
    // second round of getopt_long() to parse other options
    while ((c = getopt_long(argc, argv, opt_str, call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 'r': opt->ref_fa_fai_fn = strdup(optarg); break;
            // case 'b': cgp->var_block_size = atoi(optarg); break;
            case 'o': opt->out_vcf_fn = strdup(optarg); break; // = fopen(optarg, "w"); break;
            case 'O': if (strcmp(optarg, "v") == 0 || strcmp(optarg, "z") == 0) opt->out_vcf_type = optarg[0];
                      else _err_error_exit("\'-O/--out-type\' can only be \'v\' or \'z\'\n");
                      break;
            case 'H': opt->no_vcf_header = 1; break;
            case 'l': opt->min_sv_len = atoi(optarg); break;
            case 'T': opt->te_seq_fn = strdup(optarg); make_te_kmer_idx(opt); break;
            case 0: if (strcmp(call_var_opt[op_idx].name, "amb-base") == 0) opt->out_amb_base = 1; 
                    // else if (strcmp(call_var_opt[op_idx].name, "hifi") == 0) set_hifi_opt(opt);
                    // else if (strcmp(call_var_opt[op_idx].name, "ont") == 0) set_ont_opt(opt);
                    else if (strcmp(call_var_opt[op_idx].name, "region-file") == 0 || 
                             strcmp(call_var_opt[op_idx].name, "regions-file") == 0) opt->reg_bed_fn = strdup(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "all-ctg") == 0) opt->only_autosome = 0, opt->only_autosome_XY=0;
                    else if (strcmp(call_var_opt[op_idx].name, "autosome") == 0) opt->only_autosome = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "autosome-XY") == 0) opt->only_autosome_XY = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "mosaic") == 0) opt->out_somatic = 1;
                    else if (strcmp(call_var_opt[op_idx].name, "alpha") == 0) opt->somatic_beta_alpha = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "beta") == 0) opt->somatic_beta_beta = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "som-alt") == 0) opt->min_somatic_alt_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "som-mei-alt") == 0) opt->min_somatic_te_dp = atoi(optarg);
                    else if (strcmp(call_var_opt[op_idx].name, "max-somvar") == 0) {
                        opt->somatic_win_max_vars = strtol(optarg, &s, 10); 
                        if (*s == ',') opt->somatic_win = strtol(s+1, &s, 10); break;
                    } else if (strcmp(call_var_opt[op_idx].name, "refine-aln") == 0) opt->refine_bam = 1;
                    break;
            case 's': opt->out_somatic = 1; break;
            case 'm': opt->out_methylation = 1; break;
            case 'E': set_exclude_ctg(opt, optarg); break;
            case 'b': opt->out_bam = hts_open(optarg, "wb"); opt->out_is_cram = 0; break;
            case 'C': opt->out_bam = hts_open(optarg, "wc"); opt->out_is_cram = 1; break;
            case 'S': opt->out_bam = hts_open(optarg, "w"); opt->out_is_cram = 0; break;
            case 'c': opt->min_dp = atoi(optarg); break;
            case 'd': opt->min_alt_dp = atoi(optarg); break;
            case 'a': opt->min_af = atof(optarg); break;
            case 'M': opt->min_mq = atoi(optarg); break;
            case 'B': opt->min_bq = atoi(optarg); break;
            case 'n': opt->sample_name = strdup(optarg); break;
            case 'x': opt->noisy_reg_max_xgaps = atoi(optarg); break;
            case 'w': opt->noisy_reg_slide_win = atoi(optarg); break;
            case 'p': opt->min_hap_full_reads = atoi(optarg); break;
            // case 'f': opt->min_no_hap_full_reads = atoi(optarg); break;
            // case 'D': opt->noisy_reg_merge_dis = atoi(optarg); break;
            // case 'j': opt->min_noisy_reg_ratio = atof(optarg); break;
            case 'L': opt->max_noisy_reg_len = atoi(optarg); break;
            // case 'c': opt->end_clip_reg = atoi(optarg); break;
            // case 'F': opt->end_clip_reg_flank_win = atoi(optarg); break;
            case 'g': if (strcmp(optarg, "right") == 0 || strcmp(optarg, "r") == 0) opt->gap_aln = LONGCALLD_GAP_RIGHT_ALN;
                      else if (strcmp(optarg, "left") == 0 || strcmp(optarg, "l") == 0) opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
                      else _err_error_exit("\'-a/--gap-aln\' can only be \'left\'/\'l\' or \'right\'/\'r\'\n"); // call_var_usage();
            case 'N': opt->disable_read_realign = 1; break;
            case 't': opt->n_threads = atoi(optarg); break;
            case 'h': call_var_free_para(opt); call_var_usage();
            case 'v': fprintf(stdout, "%s\n", LONGCALLD_VERSION); call_var_free_para(opt); return 0;
            case 'V': LONGCALLD_VERBOSE = atoi(optarg); break;
            default: call_var_free_para(opt); return 0;
        }
    }

    if (argc - optind < 2) {
        call_var_free_para(opt); call_var_usage();
    // } else if (argc - optind < 2) {
        // _err_error_exit("Reference genome FASTA and alignment BAM are required.\n");
    } else {
        opt->ref_fa_fn = argv[optind++];
        opt->in_bam_fn = argv[optind++];
    }
    // check if input is (short) URL
    opt->ref_fa_fn = retrieve_full_url(opt->ref_fa_fn); opt->in_bam_fn = retrieve_full_url(opt->in_bam_fn);
    // check if hifi and ont are both set
    if (opt->is_pb_hifi && opt->is_ont) _err_error_exit("Cannot set both --hifi and --ont\n");
    if (opt->only_autosome && opt->only_autosome_XY) _err_error_exit("Cannot set both --autosome and --autosome-XY\n");
    // threads
    if (opt->n_threads <= 0) {
        opt->n_threads = 1; opt->pl_threads = 1;
    } else {
        opt->n_threads = MIN_OF_TWO(opt->n_threads, get_num_processors());
        opt->pl_threads = MIN_OF_TWO(2, opt->n_threads); // MIN_OF_TWO(opt->n_threads, CALL_VAR_PL_THREAD_N);
    }
    if (opt->out_vcf_fn == NULL) opt->out_vcf_fn = strdup("-");
    if (opt->out_vcf_type == 'z') opt->out_vcf = hts_open(opt->out_vcf_fn, "wz");
    else opt->out_vcf = hts_open(opt->out_vcf_fn, "w");
    // set up pipeline for multi-threading
    call_var_pl_t pl;
    memset(&pl, 0, sizeof(call_var_pl_t));
    pl.max_reg_len_per_chunk = LONGCALLD_BAM_CHUNK_REG_SIZE; // pl.max_reads_per_chunk = LONGCALLD_BAM_CHUNK_READ_COUNT; 
    pl.n_threads = opt->n_threads; pl.min_reg_chunks_per_run = opt->n_threads * 4;
    // open BAM file & reference genome
    call_var_pl_open_fa_bam(opt, &pl, argv+optind, argc-optind);

    // write VCF/BAM header
    if (opt->out_bam != NULL) call_var_pl_write_bam_header(opt, pl.io_aux[0].header);
    if (opt->out_vcf != NULL && (opt->out_vcf_type == 'z' || opt->no_vcf_header == 0)) write_vcf_header(pl.io_aux[0].header, opt);
    // start to work !!!
    pl.opt = opt;
    kt_pipeline(opt->pl_threads, call_var_worker_pipeline, &pl, 3);
    if (opt->out_bam != NULL) hts_close(opt->out_bam);
    if (opt->out_vcf != NULL) {
        if (opt->vcf_hdr != NULL) bcf_hdr_destroy(opt->vcf_hdr);
        hts_close(opt->out_vcf);
    }
    call_var_free_pl(pl); call_var_free_para(opt); 
    // finish
    _err_info("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.\n", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    _err_success("%s\n", CMD);
    return 0;
}