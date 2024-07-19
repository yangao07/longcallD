#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/sysinfo.h>
#include "main.h"
#include "call_var.h"
#include "bam_utils.h"
#include "vcf_utils.h"
#include "utils.h"
#include "seq.h"
#include "collect_var.h"
#include "kthread.h"

extern int LONGCALLD_VERBOSE;

const struct option call_var_opt [] = {
    { "ref-fa", 1, NULL, 'r' },
    { "out-vcf", 1, NULL, 'o'},
    { "out-bam", 1, NULL, 'b' },
    { "sample-name", 1, NULL, 'n'},
    { "min-depth", 1, NULL, 'd' },
    { "max-ploidy", 1, NULL, 'p' },
    { "max-sites", 1, NULL, 's' },
    { "win-size", 1, NULL, 'w'},
    { "threads", 1, NULL, 't' },
    { "help", 0, NULL, 'h' },
    { "version", 0, NULL, 'v' },
    // { "verbose", 1, NULL, 'V' },
    { 0, 0, 0, 0}
};

call_var_opt_t *call_var_init_para(void) {
    call_var_opt_t *opt = (call_var_opt_t*)_err_malloc(sizeof(call_var_opt_t));

    opt->sample_name = NULL;
    opt->ref_fa_fn = NULL;

    opt->region_list = NULL; opt->region_is_file = 0;

    opt->max_ploid = LONGCALLD_DEF_PLOID;
    opt->min_mq = LONGCALLD_MIN_CAND_SNP_MQ;
    opt->min_bq = LONGCALLD_MIN_CAND_SNP_BQ;
    opt->min_dp = LONGCALLD_MIN_CAND_SNP_DP;

    opt->dens_reg_max_sites = LONGCALLD_DENSE_REG_MAX_SITES;
    opt->dens_reg_slide_win = LONGCALLD_DENSE_REG_SLIDE_WIN;
    opt->indel_flank_win_size = LONGCALLD_INDEL_FLANK_WIN_SIZE;

    opt->min_af = LONGCALLD_MIN_CAND_SNP_AF;
    opt->max_af = LONGCALLD_MAX_CAND_SNP_AF;
    opt->max_low_qual_frac = LONGCALLD_MAX_LOW_QUAL_FRAC;

    opt->pl_threads = MIN_OF_TWO(CALL_VAR_PL_THREAD_N, get_nprocs());
    opt->n_threads = MIN_OF_TWO(CALL_VAR_THREAD_N, get_nprocs());

    opt->out_vcf = NULL; opt->out_bam = NULL;
    // opt->verbose = 0;
    return opt;
}

void call_var_free_para(call_var_opt_t *opt) {
    if (opt->ref_fa_fn) free(opt->ref_fa_fn);
    if (opt->sample_name) free(opt->sample_name);
    if (opt->region_list) free(opt->region_list);
    free(opt);
}

void call_var_free_pl(call_var_pl_t pl) {
    if (pl.ref_seq) ref_seq_free(pl.ref_seq);
    if (pl.bam) {
        if (pl.iter) hts_itr_destroy(pl.iter);
        if (pl.idx) hts_idx_destroy(pl.idx);
        bam_hdr_destroy(pl.header);
        sam_close(pl.bam); 
    }
}

void var_free(var_t *v) {
    if (v->vars) {
        for (int i = 0; i < v->m; ++i) {
            if (v->vars[i].ref_bases) free(v->vars[i].ref_bases);
            if (v->vars[i].alt_bases) {
                for (int j = 0; j < v->vars[i].n_allele-1; ++j)
                    free(v->vars[i].alt_bases[j]);
                free(v->vars[i].alt_bases);
                free(v->vars[i].alt_len);
            }
        }
        free(v->vars);
    }
}

// call variants for each chunk
static void call_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks + ii; var_t *v = step->vars + ii;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] tid: %d, chunk: %ld n_reads: %d\n", __func__, tid, ii, c->n_reads);
    collect_var_main(step->pl, c, v);
    // Add HP tag
    if (step->pl->opt->out_bam != NULL) {
        for (int i = 0; i < c->n_reads; ++i) {
            bam1_t *b = c->reads[i];
            if (c->is_skipped[i] || c->haps[i] <= 0) continue;
            bam_aux_append(b, "XP", 'i', 4, (uint8_t*)&(c->haps[i]));
        }
    }
}

// merge variants from two chunks
static void stitch_var_worker_for(void *_data, long ii, int tid) {
    call_var_step_t *step = (call_var_step_t*)_data;
    bam_chunk_t *c = step->chunks; var_t *var = step->vars;
    if (LONGCALLD_VERBOSE >= 1) fprintf(stderr, "[%s] tid: %d, chunk: %ld (%d), n_reads: %d\n", __func__, tid, ii, step->n_chunks, (c+ii)->n_reads);
    stitch_var_main(step->pl, c, var, ii);
}

static void call_var_pl_open_bam(call_var_opt_t *opt, call_var_pl_t *pl, char **regions, int n_regions) {
    // input BAM file
    pl->bam = sam_open(opt->in_bam_fn, "r"); if (pl->bam == NULL) _err_error_exit("Failed to open BAM file \'%s\'\n", opt->in_bam_fn);
    pl->header = sam_hdr_read(pl->bam); if (pl->header == NULL) {
        sam_close(pl->bam); _err_error_exit("Failed to read BAM header \'%s\'\n", opt->in_bam_fn);
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
                if (pl->idx == NULL)
                {
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
        } else pl->use_iter = 1;
    } else pl->iter = NULL;
}

static void call_var_pl_write_bam_header(samFile *out_bam, bam_hdr_t *header) {
    if (sam_hdr_add_pg(header, PROG, "VN", VERSION, "CL", CMD, NULL) < 0) _err_error_exit("Fail to add PG line to bam header.\n");
    if (sam_hdr_write(out_bam, header) < 0) _err_error_exit("Failed to write BAM header.\n");
}

static void call_var_pl_open_ref_fa(const char *ref_fa_fn, call_var_pl_t *pl) {
    // if (pl->bam_has_eqx_cigar || pl->bam_has_md_tag) return;
    if (ref_fa_fn) {
        pl->ref_seq = read_ref_seq(ref_fa_fn);
        if (pl->ref_seq == NULL) _err_error_exit("Failed to read reference genome: %s.\n", ref_fa_fn);
    } else {
        _err_error_exit("No reference genome provided.\n");
    }
}

static void *call_var_worker_pipeline(void *shared, int step, void *in) { // kt_pipeline() callback
    call_var_pl_t *p = (call_var_pl_t*)shared;
    if (step == 0) { // step 0: read bam records into BAM chunks
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 0: read BAM into chunks.\n");
        call_var_step_t *s;
        s = calloc(1, sizeof(call_var_step_t));
        s->pl = p;
        s->max_chunks = 64; s->chunks = calloc(s->max_chunks, sizeof(bam_chunk_t));
        int r, n_ovlp_reads=0, *ovlp_read_i=NULL;
        while ((r = collect_bam_chunk(p->bam, p->header, p->iter, p->use_iter, p->max_reg_len_per_chunk, &ovlp_read_i, &n_ovlp_reads, s->chunks+s->n_chunks)) > 0) {
            if (++s->n_chunks >= s->max_chunks) break;
        }
        if (ovlp_read_i != NULL) free(ovlp_read_i);
        // XXX merge chunks if too few reads
        if (s->n_chunks > 0) return s;
        else { free(s->chunks); free(s); return 0; }
    } else if (step == 1) { // step 1: work on the BAM chunks
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 1: call variants.\n");
        if (((call_var_step_t*)in)->n_chunks > 0) {
            ((call_var_step_t*)in)->vars = calloc(((call_var_step_t*)in)->n_chunks, sizeof(var_t));
            kt_for(p->n_threads, call_var_worker_for, in, ((call_var_step_t*)in)->n_chunks);
        }
        return in;
    } else if (step == 2) { // step 2: stitch the results
                            //         merge variants from different chunks
                            //         1) merge overlapping variants
                            //         2) update phase set and haplotype (flip if needed)
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 2: stitch variants.\n");
        // XXX stitch all results using 1 thread?
        if (((call_var_step_t*)in)->n_chunks > 1) {
            kt_for(p->n_threads, stitch_var_worker_for, in, ((call_var_step_t*)in)->n_chunks-1);
        }
        // stitch all the results, or not necessary ?
        return in;
    } else if (step == 3) { // step 3: write the buffer to output
        if (LONGCALLD_VERBOSE >= 1) _err_info("Step 3: output variants.\n");
        call_var_step_t *s = (call_var_step_t*)in;
        for (int i = 0; i < s->n_chunks; ++i) {
            bam_chunk_t *c = s->chunks + i;
            for (int j = 0; j < c->n_reads; ++j) {
                // chrome, pos, qname, hap, rlen
                if (LONGCALLD_VERBOSE >= 2)
                    fprintf(stderr, "%d\t%s\t%ld\t%s\t%d\t%d\t%ld\n", j, c->tname, c->reads[j]->core.pos+1, bam_get_qname(c->reads[j]), c->reads[j]->core.flag, c->haps[j], bam_cigar2rlen(c->reads[j]->core.n_cigar, bam_get_cigar(c->reads[j])));
                if (p->opt->out_bam != NULL)
                    if (sam_write1(p->opt->out_bam, p->header, c->reads[j]) < 0) _err_error_exit("Failed to write BAM record.");
            }
            bam_chunk_free(c); // free input
            var_free(s->vars + i);  // free output
        }
        free(s->chunks); free(s->vars); free(s);
    }
    return 0;
}

static int call_var_usage(void) {//main usage
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s (%s)\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", VERSION, CONTACT); 

    fprintf(stderr, "Usage:   %s call [options] input.bam [region ...] > phased.vcf\n", PROG);
    // fprintf(stderr, "            Note: \'ref.fa\' needs to contain the same contig/chromosome names as \'input.bam\'\n");
    // fprintf(stderr, "                  \'ref.fa\' is not needed if 1) BAM contains CIGARS with =/X, or 2) BAM contain MD tags\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Intput and output:\n");
    fprintf(stderr, "    -r --ref-fa      STR    reference genome, FASTA(.gz) file\n");
    fprintf(stderr, "                            Note: \'ref.fa\' must to contain the same contig/chromosome names as \'input.bam\'\n");
    fprintf(stderr, "                                  \'ref.fa\' is not needed if input BAM contains 1) =/X CIGARS, or 2) MD tags\n");
    fprintf(stderr, "    -n --sample-name STR  sample name. if not provided, input BAM file name will be used as sample name in output VCF [NULL]\n");
    fprintf(stderr, "    -o --out-vcf     STR  output phased VCF file [stdout]\n");
    fprintf(stderr, "    -b --out-bam     STR  output phased BAM file [NULL]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Variant:\n");
    fprintf(stderr, "    -d --min-depth   INT  minimum depth to call a SNP [%d]\n", LONGCALLD_MIN_CAND_SNP_DP);
    fprintf(stderr, "    -p --max-ploidy  INT  maximum ploidy [%d]\n", LONGCALLD_DEF_PLOID);
    fprintf(stderr, "    -s --max-sites   INT  maximum number of substitutions/gaps in a window(-w/--win-size) [%d]\n", LONGCALLD_DENSE_REG_MAX_SITES);
    fprintf(stderr, "    -w --win-size    INT  window size for noisy region [%d]\n", LONGCALLD_DENSE_REG_SLIDE_WIN);
    fprintf(stderr, "                          noisy region with more than -s subs/gaps in a window of -w bases will be skipped for initial haplotype assignment\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  General:\n");
    fprintf(stderr, "    -t --thread      INT  number of threads to use [%d]\n", MIN_OF_TWO(CALL_VAR_THREAD_N, get_nprocs()));
    fprintf(stderr, "    -h --help             print this help usage\n");
    fprintf(stderr, "    -v --version          print version number\n");
    fprintf(stderr, "    -V --verbose     INT  verbose level (0-2). 0: none, 1: information, 2: debug [0]\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    return 1;
}

int call_var_main(int argc, char *argv[]) {
    _err_cmd("%s\n", CMD);
    int c, op_idx; call_var_opt_t *opt = call_var_init_para();
    double realtime0 = realtime();
    while ((c = getopt_long(argc, argv, "r:o:b:d:n:s:w:t:hvV:", call_var_opt, &op_idx)) >= 0) {
        switch(c) {
            case 'r': opt->ref_fa_fn = strdup(optarg); break;
            // case 'b': cgp->var_block_size = atoi(optarg); break;
            case 'o': opt->out_vcf = fopen(optarg, "w"); break;
            case 'b': opt->out_bam = hts_open(optarg, "wb"); break;
            case 'd': opt->min_dp = atoi(optarg); break;
            case 'n': opt->sample_name = strdup(optarg); break;
            case 's': opt->dens_reg_max_sites = atoi(optarg); break;
            case 'w': opt->dens_reg_slide_win = atoi(optarg); break;
            case 't': opt->n_threads = atoi(optarg); break;
            case 'h': call_var_usage(); call_var_free_para(opt); return 0;
            case 'v': fprintf(stderr, "%s\n", VERSION); call_var_free_para(opt); return 0;
            case 'V': LONGCALLD_VERBOSE = atoi(optarg); break;
            default: call_var_free_para(opt); return 0; // call_var_usage();
        }
    }

    if (!isatty(STDIN_FILENO)) {
        opt->in_bam_fn = "-";
    } else {
        if (argc - optind < 1) {
            _err_error("Input BAM is required.\n");
            call_var_free_para(opt); return call_var_usage();
        } else {
            opt->in_bam_fn = argv[optind++];
            // opt->ref_fa_fn = argv[optind+1];
        }
    }
    // threads
    if (opt->n_threads <= 0) {
        opt->n_threads = 1; opt->pl_threads = 1;
    } else {
        opt->n_threads = MIN_OF_TWO(opt->n_threads, get_nprocs());
        opt->pl_threads = MIN_OF_TWO(opt->n_threads, CALL_VAR_PL_THREAD_N);
    }
    // if (opt->out_bam == NULL) // output to stdout
        // opt->out_bam = hts_open("-", "wb");
    if (opt->out_vcf == NULL) opt->out_vcf = stdout;
    // set up m-threading
    call_var_pl_t pl;
    memset(&pl, 0, sizeof(call_var_pl_t));
    // read reference genome
    if (opt->ref_fa_fn) {
        pl.ref_seq = read_ref_seq(opt->ref_fa_fn);
        if (pl.ref_seq == NULL) {
            _err_error("Failed to read reference genome: %s.\n", opt->ref_fa_fn);
            call_var_free_para(opt); hts_close(opt->out_bam); call_var_free_pl(pl);
            return 1;
        }
    }
    // open BAM file
    call_var_pl_open_bam(opt, &pl, argv+optind, argc-optind);
    // call_var_pl_open_ref_fa(opt->ref_fa_fn, &pl); // XXX


    // write VCF/BAM header
    if (opt->out_bam != NULL) call_var_pl_write_bam_header(opt->out_bam, pl.header);
    if (opt->out_vcf != NULL) write_vcf_header(pl.header, opt->out_vcf, opt->sample_name);
    // start to work !!!
    pl.opt = opt;
    pl.max_reg_len_per_chunk = LONGCALLD_BAM_CHUNK_REG_SIZE; // pl.max_reads_per_chunk = LONGCALLD_BAM_CHUNK_READ_COUNT; 
    pl.n_threads = opt->n_threads;
    kt_pipeline(opt->pl_threads, call_var_worker_pipeline, &pl, 4);
    if (opt->out_bam != NULL) hts_close(opt->out_bam);
    if (opt->out_vcf != NULL) fclose(opt->out_vcf);
    call_var_free_pl(pl); call_var_free_para(opt); 
    // finish
    _err_cmd("%s\n", CMD);
    _err_info("Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.\n", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    _err_success("Done.\n");
    return 0;
}