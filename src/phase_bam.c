#include <stdio.h>
#include <stdlib.h>
#include <getopt.h>
#include <string.h>
#include <sys/sysinfo.h>
#include "main.h"
#include "phase_bam.h"
#include "bam_utils.h"
#include "utils.h"
#include "seq.h"
#include "collect_snps.h"
#include "kthread.h"

const struct option phase_bam_opt [] = {
    { "ref-fa", 1, NULL, 'r' },
    { "out-bam", 1, NULL, 'o' },
    { "min-depth", 1, NULL, 'd' },
    { "max-ploidy", 1, NULL, 'p' },
    { "threads", 1, NULL, 't' },
    { "help", 0, NULL, 'h' },
    { 0, 0, 0, 0}
};

phase_bam_opt_t *phase_bam_init_para(void) {
    phase_bam_opt_t *pb_opt = (phase_bam_opt_t*)_err_malloc(sizeof(phase_bam_opt_t));

    pb_opt->ref_fa_fn = NULL;

    pb_opt->region_list = NULL; pb_opt->region_is_file = 0;

    pb_opt->max_ploid = LONGCALLD_DEF_PLOID;
    pb_opt->min_mq = LONGCALLD_MIN_CAND_SNP_MQ;
    pb_opt->min_bq = LONGCALLD_MIN_CAND_SNP_BQ;
    pb_opt->min_dp = LONGCALLD_MIN_CAND_SNP_DP;
    pb_opt->min_af = LONGCALLD_MIN_CAND_SNP_AF;
    pb_opt->max_af = LONGCALLD_MAX_CAND_SNP_AF;

    pb_opt->pl_threads = MIN_OF_TWO(PHASE_BAM_PL_THREAD_N, get_nprocs());
    pb_opt->n_threads = MIN_OF_TWO(PHASE_BAM_THREAD_N, get_nprocs());

    pb_opt->out_bam = NULL;
    return pb_opt;
}

void phase_bam_free_para(phase_bam_opt_t *pb_opt) {
    if (pb_opt->region_list) free(pb_opt->region_list);
    free(pb_opt);
}

void phase_bam_free_pl(phase_bam_pl_t pl) {
    if (pl.ref_seq) ref_seq_free(pl.ref_seq);
    if (pl.bam) {
        sam_close(pl.bam);
        bam_hdr_destroy(pl.header);
        hts_idx_destroy(pl.idx);
    }
}

static void phase_bam_worker_for(void *_data, long ii, int tid) // kt_for() callback
{
    phase_bam_step_t *step = (phase_bam_step_t*)_data;
    bam_chunk_t *c = step->chunks + ii;
    // fprintf(stderr, "[M::%s] tid: %d, n_reads: %d\n", __func__, tid, c->n_reads);
    collect_snps_main(step->pl, c);
    // Add HP tag
    for (int i = 0; i < c->n_reads; ++i) {
        bam1_t *b = c->reads[i];
        bam_aux_append(b, "HP", 'i', 4, (uint8_t*)&(c->HPs[i]));
    }
}

static void phase_bam_pl_open_bam(const char *in_bam_fn, phase_bam_pl_t *pl) {
    // output BAM file
    pl->bam = sam_open(in_bam_fn, "r"); if (pl->bam == NULL) _err_fatal("Error: failed to open BAM file \'%s\'", in_bam_fn);
    pl->header = sam_hdr_read(pl->bam); if (pl->header == NULL) {
        sam_close(pl->bam); _err_fatal("Error: failed to read BAM header \'%s\'", in_bam_fn);
    }
    pl->idx = sam_index_load(pl->bam, in_bam_fn);
    // if (pl->idx == NULL) {
        // bam_hdr_destroy(pl->header); sam_close(pl->bam);
        // _err_fatal("Error: failed to load BAM index \'%s\'", in_bam_fn);
    // }
    // pl->iter = sam_itr_querys(pl->idx, pl->header, region);
    // if (pl->iter == NULL) {
        // hts_idx_destroy(pl->idx); bam_hdr_destroy(pl->header); sam_close(pl->bam);
        // _err_fatal("Error: failed to get iterator for region %s\n", region);
    // }
    // check =/X in CIGAR or MD tag
    // check_eqx_cigar_MD_tag(pl->bam, pl->header, &(pl->bam_has_eqx_cigar), &(pl->bam_has_md_tag));
}

static void phase_bam_pl_write_bam_header(samFile *out_bam, bam_hdr_t *header) {
    if (sam_hdr_add_pg(header, PROG, "VN", VERSION, "CL", CMD, NULL) < 0) _err_fatal("Fail to add PG line to bam header.");
    if (sam_hdr_write(out_bam, header) < 0) _err_fatal("Failed to write BAM header.");
}

static void phase_bam_pl_open_ref_fa(const char *ref_fa_fn, phase_bam_pl_t *pl) {
    // if (pl->bam_has_eqx_cigar || pl->bam_has_md_tag) return;
    if (ref_fa_fn) {
        pl->ref_seq = read_ref_seq(ref_fa_fn);
        if (pl->ref_seq == NULL) _err_fatal("Error: failed to read reference genome: %s.", ref_fa_fn);
    } else {
        _err_fatal("Error: reference genome is required.");
    }
}

static void *phase_bam_worker_pipeline(void *shared, int step, void *in) // kt_pipeline() callback
{
    phase_bam_pl_t *p = (phase_bam_pl_t*)shared;
    if (step == 0) { // step 0: read bam records into BAM chunks
        phase_bam_step_t *s;
        // fprintf(stderr, "[M::%s] step 0\n", __func__);
        s = calloc(1, sizeof(phase_bam_step_t));
        s->pl = p;
        s->max_chunks = 64; s->chunks = calloc(s->max_chunks, sizeof(bam_chunk_t));
        int r;
        while ((r = bam_read_chunk(p->bam, p->header, s->chunks+s->n_chunks, p->max_reads_per_chunk)) > 0) {
            if (++s->n_chunks >= s->max_chunks) break;
        }
        if (s->n_chunks > 0) return s;
        else { free(s->chunks); free(s); return 0; }
    } else if (step == 1) { // step 1: work on the BAM chunks
        kt_for(p->n_threads, phase_bam_worker_for, in, ((phase_bam_step_t*)in)->n_chunks);
        return in;
    } else if (step == 2) { // step 2: write the buffer to output
        // fprintf(stderr, "[M::%s] step 2\n", __func__);
        phase_bam_step_t *s = (phase_bam_step_t*)in;
        for (int i = 0; i < s->n_chunks; ++i) {
            bam_chunk_t *c = s->chunks + i;
            for (int j = 0; j < c->n_reads; ++j) {
                // printf("%s, HP: %d\n", bam_get_qname(c->reads[j]), c->HPs[j]);
                // if (sam_write1(p->pb_opt->out_bam, p->header, c->reads[j]) < 0)
                    // _err_fatal("Failed to write BAM record.");
                bam_destroy1(c->reads[j]);
            }
            for (int j = c->n_reads; j < p->max_reads_per_chunk; ++j) bam_destroy1(c->reads[j]);
            free(c->HPs); free(c->reads);
        }
        free(s->chunks); free(s);
    }
    return 0;
}

static int phase_bam_usage(void)	//main usage
{
    fprintf(stderr, "\n");
    fprintf(stderr, "Program: %s (%s)\n", PROG, DESCRIP);
    fprintf(stderr, "Version: %s\tContact: %s\n\n", VERSION, CONTACT); 

    fprintf(stderr, "Usage:   %s phase-bam [options] input.bam [region ...] > phased.bam\n", PROG);
    // fprintf(stderr, "            Note: \'ref.fa\' needs to contain the same contig/chromosome names as \'input.bam\'\n");
    // fprintf(stderr, "                  \'ref.fa\' is not needed if 1) BAM contains CIGARS with =/X, or 2) BAM contain MD tags\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "Options:\n");
    fprintf(stderr, "  Intput and output:\n");
    fprintf(stderr, "    -r --ref-fa    STR    reference genome, FASTA(.gz) file (not needed if 1) BAM contains CIGARS with =/X, or 2) BAM contain MD tags) \n");
    fprintf(stderr, "                   Note: \'ref.fa\' needs to contain the same contig/chromosome names as \'input.bam\'\n");
    fprintf(stderr, "                         \'ref.fa\' is not needed if BAM contains 1) =/X CIGARS, or 2) MD tags\n");
    fprintf(stderr, "    -o --out-bam   STR    output phased BAM file [stdout]\n");
    fprintf(stderr, "\n");
    fprintf(stderr, "  Variant:\n");
    fprintf(stderr, "    -d --min-depth INT    minimum depth to call a SNP [%d]\n", LONGCALLD_MIN_CAND_SNP_DP);
    fprintf(stderr, "    -p --max-ploidy INT   maximum ploidy [%d]\n", LONGCALLD_DEF_PLOID);
    fprintf(stderr, "\n");
    fprintf(stderr, "  General:\n");
    fprintf(stderr, "    -t --thread    INT    number of threads to use [%d]\n", MIN_OF_TWO(PHASE_BAM_THREAD_N, get_nprocs()));
    fprintf(stderr, "    -h --help             print this help usage\n");
    fprintf(stderr, "\n");

    fprintf(stderr, "\n");
    return 1;
}

int phase_bam_main(int argc, char *argv[]) {
    int c, op_idx; phase_bam_opt_t *pb_opt = phase_bam_init_para();
    // for (i = 0; i < argc; ++i) fprintf(stderr, "%s ", argv[i]); fprintf(stderr, "\n");
    double realtime0 = realtime();
    while ((c = getopt_long(argc, argv, "r:o:d:t:h", phase_bam_opt, &op_idx)) >= 0) {
        switch(c) {
            case 'r': pb_opt->ref_fa_fn = optarg; break;
            // case 'b': cgp->var_block_size = atoi(optarg); break;
            case 'o': pb_opt->out_bam = hts_open(optarg, "wb"); break;
            case 'd': pb_opt->min_dp = atoi(optarg); break;
            case 't': pb_opt->n_threads = atoi(optarg); break;
            case 'h': phase_bam_usage(); phase_bam_free_para(pb_opt); return 0;
            default: err_printf("\n[main] Error: unknown option: -%c(%d) %s.\n", c, c, optarg);
                     phase_bam_free_para(pb_opt); return phase_bam_usage();
        }
    }
    if (argc - optind < 1) {
        err_func_printf(__func__, "Error: reference FASTA and input BAM are required.");
        phase_bam_free_para(pb_opt); return phase_bam_usage();
    } else {
        pb_opt->in_bam_fn = argv[optind];
        // pb_opt->ref_fa_fn = argv[optind+1];
    }
    // BAM regions
    // if (argc - optind > 2) {
    //     int i;
    //     kstring_t tmp = {0,0,0};
    //     kputs(argv[optind+1], &tmp);
    //     for (i=optind+2; i<argc; i++) { kputc(',',&tmp); kputs(argv[i],&tmp); }
    //     if ( bcf_sr_set_regions(args->files, tmp.s, 0)<0 )
    //         error("Failed to read the regions: %s\n", tmp.s);
    //     free(tmp.s);
    // }

    // threads
    if (pb_opt->n_threads <= 0) {
        pb_opt->n_threads = 1; pb_opt->pl_threads = 1;
    } else {
        pb_opt->n_threads = MIN_OF_TWO(pb_opt->n_threads, get_nprocs());
        pb_opt->pl_threads = MIN_OF_TWO(pb_opt->n_threads, PHASE_BAM_PL_THREAD_N);
    }
    if (pb_opt->out_bam == NULL) // output to stdout
        pb_opt->out_bam = hts_open("-", "wb");
    // set up m-threading
    phase_bam_pl_t pl;
    memset(&pl, 0, sizeof(phase_bam_pl_t));
    // read reference genome
    if (pb_opt->ref_fa_fn) {
        pl.ref_seq = read_ref_seq(pb_opt->ref_fa_fn);
        if (pl.ref_seq == NULL) {
            err_func_printf(__func__, "Error: failed to read reference genome: %s.", pb_opt->ref_fa_fn);
            phase_bam_free_para(pb_opt); hts_close(pb_opt->out_bam); phase_bam_free_pl(pl);
            return 1;
        }
    }
    // open BAM file
    phase_bam_pl_open_bam(pb_opt->in_bam_fn, &pl);
    // phase_bam_pl_open_ref_fa(pb_opt->ref_fa_fn, &pl); // XXX

    // write BAM header
    phase_bam_pl_write_bam_header(pb_opt->out_bam, pl.header);
    // start to work !!!
    pl.pb_opt = pb_opt;
    pl.max_reads_per_chunk = 512; pl.n_threads = pb_opt->n_threads;
    kt_pipeline(pb_opt->pl_threads, phase_bam_worker_pipeline, &pl, 3);
    err_func_printf(__func__, "Real time: %.3f sec; CPU: %.3f sec; Peak RSS: %.3f GB.", realtime() - realtime0, cputime(), peakrss() / 1024.0 / 1024.0 / 1024.0);
    hts_close(pb_opt->out_bam); phase_bam_free_pl(pl); phase_bam_free_para(pb_opt); 
    return 0;
}