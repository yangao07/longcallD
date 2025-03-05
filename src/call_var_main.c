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

extern int LONGCALLD_VERBOSE;

const struct option call_var_opt [] = {
    // long options
    { "amb-base", 0, NULL, 0},
    { "hifi", 0, NULL, 0},
    { "ont", 0, NULL, 0},
    { "ont-hp-prof", 1, NULL, 0},

    // short options
    { "ref-idx", 1, NULL, 'r' },
    { "out-vcf", 1, NULL, 'o'},
    { "out-bam", 1, NULL, 'b' },
    { "no-vcf-header", 0, NULL, 'H'},
    { "gap-aln", 1, NULL, 'g'},
    // { "out-bam", 1, NULL, 'b' },
    { "sample-name", 1, NULL, 'n'},

    { "min-depth", 1, NULL, 'c' },
    { "alt-depth", 1, NULL, 'd' },
    { "alt-ratio", 1, NULL, 'a' },
    // { "max-ploidy", 1, NULL, 'p' },

    { "max-xgap", 1, NULL, 'x' },
    { "win-size", 1, NULL, 'w'},
    { "noisy-rat", 1, NULL, 'j' },
    { "noisy-flank", 1, NULL, 'f' },
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
    opt->is_pb_hifi = 1;
    opt->noisy_reg_max_xgaps = 5;
}

void set_ont_opt(call_var_opt_t *opt) {
    opt->is_ont = 1;
    opt->noisy_reg_max_xgaps = 20;
}

// profile format
// CTA     5       +       5       10661016.0      0.8805382103122837
// beg/hp/end   ref_len strand alt_len count prob
void set_ont_hp_prof(call_var_opt_t *opt, char *prof_fn) {
    opt->ont_hp_profile = (ont_hp_profile_t***)malloc(4*4*4*sizeof(ont_hp_profile_t**));
    for (int i = 0; i < 4*4*4; ++i) {
        opt->ont_hp_profile[i] = (ont_hp_profile_t**)malloc(50*sizeof(ont_hp_profile_t*));
        for (int j = 0; j < 50; ++j) {
            opt->ont_hp_profile[i][j] = (ont_hp_profile_t*)malloc(2*sizeof(ont_hp_profile_t));
            for (int k = 0; k < 2; ++k) {
                opt->ont_hp_profile[i][j][k].hp_len_to_prob = (double*)malloc(20*sizeof(double));
                opt->ont_hp_profile[i][j][k].alt_hp_lens = (int*)malloc(20*sizeof(int));
                opt->ont_hp_profile[i][j][k].beg_flank_base = i/16;
                opt->ont_hp_profile[i][j][k].hp_base = (i%16)/4;
                opt->ont_hp_profile[i][j][k].end_flank_base = i%4;
                opt->ont_hp_profile[i][j][k].ref_hp_len = j+1;
                opt->ont_hp_profile[i][j][k].strand = k;
                opt->ont_hp_profile[i][j][k].n_alt_hp_lens = 0;
            }
        }
    }
    // load profile
    FILE *fp = fopen(prof_fn, "r");
    if (fp == NULL) _err_error_exit("Failed to open HP profile file: %s", prof_fn);
    char line[1024];
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '#') continue; // skip header/comment
        char *token;
        token = strtok(line, "\t");
        char beg_base = token[0], hp_base = token[1], end_base = token[2];
        token = strtok(NULL, "\t");
        int ref_len = atoi(token);
        if (ref_len < 5 || ref_len > 50) continue; // skip too short or too long ref length
        token = strtok(NULL, "\t");
        char strand = token[0];
        // int hp_idx = nst_nt4_table[(int)beg_base]*16 + nst_nt4_table[(int)hp_base]*4 + nst_nt4_table[(int)end_base];
        int hp_base_idx = nst_nt4_table[(int)beg_base]*16 + nst_nt4_table[(int)hp_base]*4 + nst_nt4_table[(int)end_base];
        int hp_ref_idx = ref_len-1;
        int hp_strand_idx = (strand == '+') ? 0 : 1;
        // fprintf(stderr, "HP idx: %d, %c-%c-%c\n", hp_idx, beg_base, hp_base, end_base);
        ont_hp_profile_t *hp_profile = &(opt->ont_hp_profile[hp_base_idx][hp_ref_idx][hp_strand_idx]);
        if (hp_profile->n_alt_hp_lens >= 20) continue;

        token = strtok(NULL, "\t");
        int alt_len = atoi(token);
        token = strtok(NULL, "\t");
        double count = atof(token);
        token = strtok(NULL, "\t");
        double prob = atof(token);

        hp_profile->ref_hp_len = ref_len;
        hp_profile->alt_hp_lens[hp_profile->n_alt_hp_lens] = alt_len;
        hp_profile->hp_len_to_prob[hp_profile->n_alt_hp_lens] = prob;
        hp_profile->n_alt_hp_lens++;
    }
    fclose(fp);
}

call_var_opt_t *call_var_init_para(void) {
    call_var_opt_t *opt = (call_var_opt_t*)_err_malloc(sizeof(call_var_opt_t));

    opt->sample_name = NULL;
    opt->ref_fa_fn = NULL; opt->ref_fa_fai_fn = NULL;

    opt->is_pb_hifi = 0; opt->is_ont = 0;

    opt->region_list = NULL; opt->region_is_file = 0;

    opt->max_ploid = LONGCALLD_DEF_PLOID;
    opt->min_mq = LONGCALLD_MIN_CAND_MQ;
    opt->min_bq = LONGCALLD_MIN_CAND_BQ;
    opt->min_dp = LONGCALLD_MIN_CAND_DP;
    opt->min_alt_dp = LONGCALLD_MIN_ALT_DP;

    opt->noisy_reg_flank_len = LONGCALLD_NOISY_REG_FLANK_LEN;

    opt->noisy_reg_max_xgaps = LONGCALLD_NOISY_REG_MAX_XGAPS;
    opt->noisy_reg_slide_win = LONGCALLD_NOISY_REG_SLIDE_WIN;

    opt->end_clip_reg = LONGCALLD_NOISY_END_CLIP;
    opt->end_clip_reg_flank_win = LONGCALLD_NOISY_END_CLIP_WIN;

    opt->max_var_ratio_per_read = LONGCALLD_MAX_VAR_RATIO_PER_READ;
    opt->max_noisy_reg_reads = LONGCALLD_MAX_NOISY_REG_READS;
    opt->max_noisy_reg_len  = LONGCALLD_MAX_NOISY_REG_LEN;
    opt->min_noisy_reg_reads = LONGCALLD_NOISY_REG_READS;
    opt->min_noisy_reg_ratio = LONGCALLD_NOISY_REG_RATIO;
    opt->max_noisy_frac_per_read = LONGCALLD_MAX_NOISY_FRAC_PER_READ;
    opt->min_hap_full_reads = LONGCALLD_MIN_HAP_FULL_READS;
    opt->min_hap_reads = LONGCALLD_MIN_HAP_READS;
    opt->min_no_hap_full_reads = LONGCALLD_MIN_NO_HAP_FULL_READS;

    opt->min_somatic_af = LONGCALLD_MIN_SOMATIC_AF;
    opt->min_af = LONGCALLD_MIN_CAND_AF;
    opt->max_af = LONGCALLD_MAX_CAND_AF;

    opt->match = LONGCALLD_MATCH_SCORE;
    opt->mismatch = LONGCALLD_MISMATCH_SCORE;
    opt->gap_open1 = LONGCALLD_GAP_OPEN1_SCORE;
    opt->gap_ext1 = LONGCALLD_GAP_EXT1_SCORE;
    opt->gap_open2 = LONGCALLD_GAP_OPEN2_SCORE;
    opt->gap_ext2 = LONGCALLD_GAP_EXT2_SCORE;
    opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
    opt->disable_read_realign = 1; // XXX right now, always disable read realignment
    opt->ont_hp_profile = NULL;

    opt->pl_threads = MIN_OF_TWO(CALL_VAR_PL_THREAD_N, get_num_processors());
    opt->n_threads = MIN_OF_TWO(CALL_VAR_THREAD_N, get_num_processors());

    opt->p_error = 0.001; opt->log_p = -3.0; opt->log_1p = log10(1-opt->p_error); opt->log_2 = 0.301023;
    opt->max_gq = 60; opt->max_qual = 60;
    opt->out_vcf = NULL; opt->no_vcf_header = 0; opt->out_amb_base = 0;
    opt->out_bam = NULL; opt->no_bam_header = 0;
    // opt->verbose = 0;
    return opt;
}

void call_var_free_para(call_var_opt_t *opt) {
    if (opt->ref_fa_fn) free(opt->ref_fa_fn);
    if (opt->in_bam_fn) free(opt->in_bam_fn);
    if (opt->sample_name) free(opt->sample_name);
    if (opt->region_list) free(opt->region_list);
    if (opt->ont_hp_profile != NULL) {
        for (int i = 0; i < 4*4*4; ++i) {
            for (int j = 0; j < 50; ++j) {
                for (int k = 0; k < 2; ++k) {
                    free(opt->ont_hp_profile[i][j][k].alt_hp_lens);
                    free(opt->ont_hp_profile[i][j][k].hp_len_to_prob);
                }
                free(opt->ont_hp_profile[i][j]);
            }
            free(opt->ont_hp_profile[i]);
        }
        free(opt->ont_hp_profile);
    }
    free(opt);
}

void call_var_free_pl(call_var_pl_t pl) {
    // if (pl.ref_seq) ref_seq_free(pl.ref_seq);
    if (pl.ref_reg_seq) ref_reg_seq_free(pl.ref_reg_seq);
    if (pl.fai) fai_destroy(pl.fai);
    if (pl.bam) {
        if (pl.iter) hts_itr_destroy(pl.iter);
        if (pl.idx) hts_idx_destroy(pl.idx);
        bam_hdr_destroy(pl.header);
        sam_close(pl.bam); 
        hts_tpool_destroy(pl.p);
    }
    bam_ovlp_chunk_free(pl.last_chunk); free(pl.last_chunk);
    if (pl.last_chunk_read_i != NULL) free(pl.last_chunk_read_i);
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
}

void call_var_pl_open_ref_reg_fa0(const char *ref_fa_fn, call_var_pl_t *pl) {
    if (ref_fa_fn) {
        _err_info("Loading reference genome: %s\n", ref_fa_fn);
        pl->ref_reg_seq = read_ref_reg_seq(ref_fa_fn);
        if (pl->ref_reg_seq == NULL) _err_error_exit("Failed to read reference genome: %s\n", ref_fa_fn);
        _err_info("Loading done: %s\n", ref_fa_fn);
    } else {
        _err_error_exit("No reference genome provided\n");
    }
}

ref_reg_seq_t *call_var_pl_open_ref_reg_fa(const char *ref_fa_fn, faidx_t *fai, hts_itr_t *iter, bam_hdr_t *header) {
    ref_reg_seq_t *ref_reg_seq = ref_reg_seq_init();
    _err_info("Loading region(s) from reference genome: %s\n", ref_fa_fn);
    for (int i = 0; i < iter->n_reg; ++i) {
        for (int j = 0; j < iter->reg_list[i].count; ++j) { // ACGT:1234, (beg, end]:(0,4]
            read_ref_reg_seq1(fai, ref_reg_seq, header->target_name[iter->reg_list[i].tid], iter->reg_list[i].intervals[j].beg, iter->reg_list[i].intervals[j].end);
        }
    }
    cr_index(ref_reg_seq->reg_cr);
    // _err_info("Loading regions done!\n");
    return ref_reg_seq;
} 

// call variants for each chunk
static void call_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks + ii; // var_t *v = step->vars + ii;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] tid: %d, chunk: %ld n_reads: %d\n", __func__, tid, ii, c->n_reads);
    collect_var_main(step->pl, c);
}

// merge variants from two chunks
static void stitch_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks+ii; var_t *var = step->vars+ii;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] tid: %d, chunk: %ld (%d), n_reads: %d\n", __func__, tid, ii, step->n_chunks, c->n_reads);
    stitch_var_main(step, c, var, ii);
}

static void call_var_pl_open_fa_bam(call_var_opt_t *opt, call_var_pl_t *pl, char **regions, int n_regions) {
    // input BAM file
    if (strcmp(opt->in_bam_fn, "-") == 0) _err_info("Opening BAM file from pipe/stdin\n");
    else _err_info("Opening BAM file: %s\n", opt->in_bam_fn);
    pl->bam = sam_open(opt->in_bam_fn, "r"); if (pl->bam == NULL) _err_error_exit("Failed to open alignment file \'%s\'\n", opt->in_bam_fn);
    const htsFormat *fmt = hts_get_format(pl->bam); if (!fmt) _err_error_exit("Failed to get format of alignment file \'%s\'\n", opt->in_bam_fn);
    // multi-threading
    pl->p = hts_tpool_init(MAX_OF_TWO(1, pl->n_threads));
    // pl->p = hts_tpool_init(3);
    htsThreadPool thread_pool = {pl->p, 0};
    if (hts_set_thread_pool(pl->bam, &thread_pool) != 0) _err_error_exit("Failed to set thread pool.\n");

    pl->header = sam_hdr_read(pl->bam); if (pl->header == NULL) {
        hts_tpool_destroy(pl->p); sam_close(pl->bam);
        _err_error_exit("Failed to read BAM header \'%s\'\n", opt->in_bam_fn);
    }
    // init chunk
    pl->last_chunk = calloc(1, sizeof(bam_chunk_t)); pl->last_chunk_read_i = NULL; pl->n_last_chunk_reads = 0; pl->cur_active_reg_beg = -1; pl->cur_active_reg_beg = -1;
    // input reference fasta file
    // pl->fai = fai_load(opt->ref_fa_fn);
    pl->fai = fai_load3(opt->ref_fa_fn, opt->ref_fa_fai_fn, NULL, FAI_CREATE);
    if (pl->fai == NULL) {
        hts_tpool_destroy(pl->p); sam_hdr_destroy(pl->header); sam_close(pl->bam);
        _err_error_exit("Failed to load/build reference fasta index: %s\n", opt->ref_fa_fn);
    }
    if (fmt->format == cram) {
        if (hts_set_fai_filename(pl->bam, opt->ref_fa_fn) != 0) {
            fai_destroy(pl->fai); hts_tpool_destroy(pl->p); sam_hdr_destroy(pl->header); sam_close(pl->bam);
            _err_error_exit("Failed to set reference file for CRAM decoding: %s %s\n", opt->in_bam_fn, opt->ref_fa_fn);
        }
    }
    // collect sample name (SM) from BAM header
    if (opt->sample_name == NULL) opt->sample_name = extract_sample_name_from_bam_header(pl->header); // extract sample names from BAM header
    if (opt->sample_name == NULL) opt->sample_name = strdup(opt->in_bam_fn);
    pl->use_iter = 0;
    if (n_regions > 0) {
        if (strcmp(opt->in_bam_fn, "-") == 0) _err_error_exit("When provided with region(s), input BAM file cannot be pipe/stdin(-).\n");
        pl->idx = sam_index_load(pl->bam, opt->in_bam_fn);
        if (pl->idx == NULL) { // attempt to create index if not exist
            _err_warning("BAM index not found for \'%s\', creating now ...\n", opt->in_bam_fn);
            if (sam_index_build(opt->in_bam_fn, 0) < 0) {
                bam_hdr_destroy(pl->header);
                sam_close(pl->bam);
                _err_error_exit("Failed to build BAM index \'%s\'\n", opt->in_bam_fn);
            } else {
                pl->idx = sam_index_load(pl->bam, opt->in_bam_fn);
                if (pl->idx == NULL) {
                    bam_hdr_destroy(pl->header);
                    sam_close(pl->bam);
                    _err_error_exit("Failed to load BAM index \'%s\'\n", opt->in_bam_fn);
                }
            }
        }
        pl->iter = sam_itr_regarray(pl->idx, pl->header, regions, n_regions);
        if (pl->iter == NULL) {
            _err_error("Failed to parse provided regions: ");
            for (int i = 0; i < n_regions; ++i) fprintf(stderr, "%s ", regions[i]);
            fprintf(stderr, "\n");
            _err_error("Will process the whole BAM file.\n");
            call_var_pl_open_ref_reg_fa0(opt->ref_fa_fn, pl);
        } else {
            pl->use_iter = 1;
            // read reference fasta based on intervals
            pl->ref_reg_seq = call_var_pl_open_ref_reg_fa(opt->ref_fa_fn, pl->fai, pl->iter, pl->header);
        }
    } else {
        pl->iter = NULL;
        // read whole reference fasta
        call_var_pl_open_ref_reg_fa0(opt->ref_fa_fn, pl);
    }
}

static void call_var_pl_write_bam_header(samFile *out_bam, bam_hdr_t *header) {
    if (sam_hdr_add_pg(header, PROG, "VN", VERSION, "CL", CMD, NULL) < 0) _err_error_exit("Fail to add PG line to bam header.\n");
    if (sam_hdr_write(out_bam, header) < 0) _err_error_exit("Failed to write BAM header.\n");
}

static void *call_var_worker_pipeline(void *shared, int step, void *in) { // kt_pipeline() callback
    call_var_pl_t *p = (call_var_pl_t*)shared;
    if (step == 0) { // step 0: read bam records into BAM chunks
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 0: read BAM into chunks.\n");
        call_var_step_t *s;
        s = calloc(1, sizeof(call_var_step_t));
        s->pl = p;
        s->max_chunks = 4*p->opt->n_threads; s->chunks = calloc(s->max_chunks, sizeof(bam_chunk_t));
        int r;
        // for each round, collect s->max_chunks of reads
        // TODO:
        // XXX use the variant calling result of current round to guide the reg_beg of next round 
        // adjacent rounds: try to overlap at least one high-quality variant
        // XXX (maybe not necessary) if overlapping part is too large, enlarge the total chunk size for the next round
        // if there is a variant gap, simply start with the first read of the next round
        // this enables the phasing across multiple rounds
        // XXX add a tmp_bam_chunk to store the last chunk of last round, stitch it with the first chunk of current round
        while (!p->reach_bam_end) {
            r = collect_bam_chunk(p, s->chunks, s->n_chunks);
            if (s->chunks[s->n_chunks].n_reads == 0) break;
            if (++s->n_chunks >= s->max_chunks) break;
            if (r < 0) {
                p->reach_bam_end = 1;
                break;
            }
        }
        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "n_chunks: %d\n", s->n_chunks);
        // XXX merge chunks if too few reads
        if (s->n_chunks > 0) {
            // copy last chunk to the first chunk
            if (p->n_last_chunk_reads > 0) {
                bam_ovlp_chunk_free(p->last_chunk);
                copy_bam_chunk0(s->chunks+s->n_chunks-1, p->last_chunk);
            }
            return s;
        } else {
            free(s->chunks); free(s);
            return 0;
        }
    } else if (step == 1) { // step 1: work on the BAM chunks
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 1: call variants.\n");
        if (((call_var_step_t*)in)->n_chunks > 0) {
            ((call_var_step_t*)in)->vars = calloc(((call_var_step_t*)in)->n_chunks, sizeof(var_t));
            // if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "n_chunks: %d\n", ((call_var_step_t*)in)->n_chunks);
            kt_for(p->n_threads, call_var_worker_for, in, ((call_var_step_t*)in)->n_chunks);
        }
        return in;
    } else if (step == 2) { // step 2: stitch the results
                            //         merge variants from different chunks
                            //         1) merge overlapping variants
                            //         2) update phase set and haplotype (flip if needed)
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 2: stitch variants.\n");
        // XXX read BAM in kt_for, create iterater for each chunk
        kt_for(p->n_threads, stitch_var_worker_for, in, ((call_var_step_t*)in)->n_chunks);
        return in;
    } else if (step == 3) { // step 3: write the buffer to output
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 3: output variants.\n");
        call_var_step_t *s = (call_var_step_t*)in;
        for (int i = 0; i < s->n_chunks; ++i) {
            // if (i != 0) continue;
            bam_chunk_t *c = s->chunks + i;
            if (p->opt->out_bam != NULL) {
                bam_chunk_t *pre_c = NULL;
                if (i > 0) pre_c = s->chunks + i - 1;
                else pre_c = p->last_chunk;
                for (int j = 0; j < c->n_reads; ++j) {
                    bam1_t *b = c->reads[j];
                    int hap = c->haps[j]; hts_pos_t ps = c->PS[j]; // Add HP+PS tag
                    if (j < c->n_up_ovlp_reads) { // ovlp reads
                        int pre_read_i = c->up_ovlp_read_i[j];
                        int pre_hap = pre_c->haps[pre_read_i]; hts_pos_t pre_ps = pre_c->PS[pre_read_i];
                        if (hap == 0) hap = pre_hap;
                        if (ps == -1) ps = pre_ps;
                    }
                    if (hap != 0) {
                        // if (c->flip_hap) hap ^= 3; // check if HP tag exists
                        uint8_t *hp = bam_aux_get(b, "HP");
                        if (hp != NULL) {
                            if (hap != bam_aux2i(hp)) {
                                bam_aux_del(b, hp);
                                bam_aux_append(b, "HP", 'i', 4, (uint8_t*)&(hap));
                            }
                        } else bam_aux_append(b, "HP", 'i', 4, (uint8_t*)&(hap));
                    }
                    if (ps != -1) { // check if PS tag exists
                        uint8_t *ps_tag = bam_aux_get(b, "PS");
                        if (ps_tag != NULL) {
                            if (ps != bam_aux2i(ps_tag)) {
                                bam_aux_del(b, ps_tag);
                                bam_aux_append(b, "PS", 'i', 4, (uint8_t *)&(ps));
                            }
                        } else bam_aux_append(b, "PS", 'i', 4, (uint8_t *)&(ps));
                    }
                    if (sam_write1(p->opt->out_bam, p->header, b) < 0) _err_error_exit("Failed to write BAM record.");
                }
            }
            write_var_to_vcf(s->vars+i, p->opt, c->tname);
            var_free(s->vars + i);  // free output
        }
        int64_t n_processed_reads = 0;
        for (int i = 0; i < s->n_chunks; ++i) {
            n_processed_reads += s->chunks[i].n_reads;
        }
        _err_info("Processed %ld reads.\n", n_processed_reads);
        bam_chunks_free(s->chunks, s->n_chunks); // free input
        free(s->vars); free(s);
    }
    return 0;
}

static void call_var_usage(void) {//main usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s (%s)\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", VERSION, CONTACT); 

    fprintf(stderr, "Usage: %s call [options] <ref.fa> <input.sam/bam/cram> [region ...] > phased.vcf\n", PROG);
    fprintf(stderr, "       Note: \'ref.fa\' should contain the same contig/chromosome names as \'input.sam/bam/cram\'\n");
    fprintf(stderr, "             \'input.sam/bam/cram' should be sorted\n");
    // fprintf(stderr, "                  \'ref.fa\' is not needed if 1) BAM contains CIGARS with =/X, or 2) BAM contain MD tags\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Intput and output:\n");
    fprintf(stderr, "    -r  --ref-idx    STR  .fai index file for reference genome FASTA file, detect automaticaly if not provided [NULL]\n");
    // fprintf(stderr, "    -b --out-bam     STR  output phased BAM file [NULL]\n");
    fprintf(stderr, "    -n --sample-name STR  sample name. 'RG/SM' in BAM header will be used if not provided [NULL]\n");
    fprintf(stderr, "                          BAM file name will be used if not provided and not 'RG/SM' in BAM header [NULL]\n");
    fprintf(stderr, "    -o --out-vcf     STR  output phased VCF file [stdout]\n");
    fprintf(stderr, "    -H --no-vcf-header    do NOT output VCF header\n");
    fprintf(stderr, "       --amb-base         output variant with ambiguous base [False]\n");
    // fprintf(stderr, "    -b --out-bam     STR  output phased BAM file [NULL]\n");
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
    // fprintf(stderr, "    -p --max-ploidy  INT  max. ploidy [%d]\n", LONGCALLD_DEF_PLOID);
    fprintf(stderr, "  Variant calling in noisy regions:\n");
    fprintf(stderr, "    -x --max-xgap    INT  max. number of allowed substitutions/gap-bases in a sliding window(-w/--win-size) [%d]\n", LONGCALLD_NOISY_REG_MAX_XGAPS);
    fprintf(stderr, "                          window with more than -x subs/gap-bases will be considered as noisy region\n");
    fprintf(stderr, "    -w --win-size    INT  window size for searching noisy region [%d]\n", LONGCALLD_NOISY_REG_SLIDE_WIN);
    fprintf(stderr, "    -j --noisy-rat FLOAT  min. ratio of noisy reads in a window to call a noisy region [%.2f]\n", LONGCALLD_NOISY_REG_RATIO);
    // fprintf(stderr, "    -f --noisy-flank INT  flanking mask window size for noisy region [%d]\n", LONGCALLD_DENSE_FLANK_WIN);
    // fprintf(stderr, "    -c --end-clip    INT  max. number of clipping bases on both ends [%d]\n", LONGCALLD_NOISY_END_CLIP);
    // fprintf(stderr, "                          end-clipping region with more than -c bases will be considered as noisy clipping region\n");
    // fprintf(stderr, "    -F --clip-flank  INT  flanking mask window size for noisy clipping region [%d]\n", LONGCALLD_NOISY_END_CLIP_WIN);
    // fprintf(stderr, "\n");
    fprintf(stderr, "    -p --hap-read    INT  when haplotype is available, min. number of full-spanning reads for each haplotype in noisy region to call a variant [%d]\n", LONGCALLD_MIN_HAP_FULL_READS);
    fprintf(stderr, "    -f --full-read   INT  when haplotype is not available, min. number of full-spanning reads in noisy region to call a variant [%d]\n", LONGCALLD_MIN_NO_HAP_FULL_READS);

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
    int c, op_idx; call_var_opt_t *opt = call_var_init_para();
    double realtime0 = realtime();
    while ((c = getopt_long(argc, argv, "r:o:Hb:c:d:a:n:x:w:j:L:f:p:g:Nt:hvV:", call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 'r': opt->ref_fa_fai_fn = strdup(optarg); break;
            // case 'b': cgp->var_block_size = atoi(optarg); break;
            case 'o': opt->out_vcf = fopen(optarg, "w"); break;
            case 'H': opt->no_vcf_header = 1; break;
            case 0: if (strcmp(call_var_opt[op_idx].name, "amb-base") == 0) opt->out_amb_base = 1; 
                    else if (strcmp(call_var_opt[op_idx].name, "hifi") == 0) set_hifi_opt(opt);
                    else if (strcmp(call_var_opt[op_idx].name, "ont") == 0) set_ont_opt(opt);
                    else if (strcmp(call_var_opt[op_idx].name, "ont-hp-prof") == 0) set_ont_hp_prof(opt, optarg);
                    break;
            case 'b': opt->out_bam = hts_open(optarg, "wb"); break;
            case 'c': opt->min_dp = atoi(optarg); break;
            case 'd': opt->min_alt_dp = atoi(optarg); break;
            case 'a': opt->min_af = atof(optarg); break;
            case 'n': opt->sample_name = strdup(optarg); break;
            case 'x': opt->noisy_reg_max_xgaps = atoi(optarg); break;
            case 'w': opt->noisy_reg_slide_win = atoi(optarg); break;
            case 'p': opt->min_hap_full_reads = atoi(optarg); break;
            case 'f': opt->min_no_hap_full_reads = atoi(optarg); break;
            case 'j': opt->min_noisy_reg_ratio = atof(optarg); break;
            case 'L': opt->max_noisy_reg_len = atoi(optarg); break;
            // case 'f': opt->dens_reg_flank_win = atoi(optarg); break;
            // case 'c': opt->end_clip_reg = atoi(optarg); break;
            // case 'F': opt->end_clip_reg_flank_win = atoi(optarg); break;
            case 'g': if (strcmp(optarg, "right") == 0 || strcmp(optarg, "r") == 0) opt->gap_aln = LONGCALLD_GAP_RIGHT_ALN;
                      else if (strcmp(optarg, "left") == 0 || strcmp(optarg, "l") == 0) opt->gap_aln = LONGCALLD_GAP_LEFT_ALN;
                      else _err_error_exit("\'-a/--gap-aln\' can only be \'left\'/\'l\' or \'right\'/\'r\'\n"); // call_var_usage();
            case 'N': opt->disable_read_realign = 1; break;
            case 't': opt->n_threads = atoi(optarg); break;
            case 'h': call_var_free_para(opt); call_var_usage();
            case 'v': fprintf(stderr, "%s\n", VERSION); call_var_free_para(opt); return 0;
            case 'V': LONGCALLD_VERBOSE = atoi(optarg); break;
            default: call_var_free_para(opt); call_var_usage();
        }
    }

    if (!isatty(STDIN_FILENO) && argc - optind == 1) { // bam file is piped, only ref.fa is provided, not working with region when bam is from pipe
        opt->in_bam_fn = "-";
        if (argc - optind < 1) {
            _err_error("Reference FASTA is required.\n");
            call_var_free_para(opt); call_var_usage();
        } else {
            opt->ref_fa_fn = argv[optind++];
        }
    } else {
        if (argc - optind < 2) {
            _err_error("Reference genome FASTA and alignment BAM are required.\n");
            call_var_free_para(opt); call_var_usage();
        } else {
            opt->ref_fa_fn = argv[optind++];
            opt->in_bam_fn = argv[optind++];
        }
    }
    // check if input is (short) URL
    opt->ref_fa_fn = retrieve_full_url(opt->ref_fa_fn); opt->in_bam_fn = retrieve_full_url(opt->in_bam_fn);
    // check if hifi and ont are both set
    if (opt->is_pb_hifi && opt->is_ont) _err_error_exit("Cannot set both --hifi and --ont\n");
    // threads
    if (opt->n_threads <= 0) {
        opt->n_threads = 1; opt->pl_threads = 1;
    } else {
        opt->n_threads = MIN_OF_TWO(opt->n_threads, get_num_processors());
        opt->pl_threads = MIN_OF_TWO(2, opt->n_threads); // MIN_OF_TWO(opt->n_threads, CALL_VAR_PL_THREAD_N);
    }
    if (opt->out_vcf == NULL) opt->out_vcf = stdout;
    // set up pipeline for multi-threading
    call_var_pl_t pl;
    memset(&pl, 0, sizeof(call_var_pl_t));
    pl.max_reg_len_per_chunk = LONGCALLD_BAM_CHUNK_REG_SIZE; // pl.max_reads_per_chunk = LONGCALLD_BAM_CHUNK_READ_COUNT; 
    pl.ovlp_region_len = LONGCALLD_BAM_CHUNK_OVLP_REG_SIZE;
    pl.n_threads = opt->n_threads;
    // open BAM file & reference genome
    call_var_pl_open_fa_bam(opt, &pl, argv+optind, argc-optind);

    // write VCF/BAM header
    if (opt->out_bam != NULL) call_var_pl_write_bam_header(opt->out_bam, pl.header);
    if (opt->out_vcf != NULL && opt->no_vcf_header == 0) write_vcf_header(pl.header, opt->out_vcf, opt->sample_name);
    // start to work !!!
    pl.opt = opt;
    kt_pipeline(opt->pl_threads, call_var_worker_pipeline, &pl, 4);
    if (opt->out_bam != NULL) hts_close(opt->out_bam);
    if (opt->out_vcf != NULL) fclose(opt->out_vcf);
    call_var_free_pl(pl); call_var_free_para(opt); 
    // finish
    _err_info("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.\n", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    _err_success("%s\n", CMD);
    return 0;
}