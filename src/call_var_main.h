#ifndef LONGCALLD_CALL_VAR_MAIN_H
#define LONGCALLD_CALL_VAR_MAIN_H

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "seq.h"
#include "kmer.h"
#include "cgranges.h"
#include "wavefront/wavefront_align.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"
#include "htslib/thread_pool.h"

#define CALL_VAR_PL_THREAD_N 2 // num of threads for pipeline
#define CALL_VAR_THREAD_N 8    // num of threads for variant calling

#define LONGCALLD_MIN_CAND_MQ 30 // ignore reads with MAPQ < 30
#define LONGCALLD_MIN_CAND_BQ 10 // ignore bases with BQ < 10
#define LONGCALLD_MIN_CAND_DP 5 // total depth < 5: skipped
#define LONGCALLD_MIN_ALT_DP 2 // alt depth < 2: skipped
#define LONGCALLD_MIN_CAND_AF 0.20 // AF < 0.20: not germline het.
#define LONGCALLD_MAX_CAND_AF 0.80 // AF > 0.80: not germline het.
#define LONGCALLD_DEF_PLOID 2 // diploid
#define LONGCALLD_REF_ALLELE 0
#define LONGCALLD_ALT_ALLELE1 1
#define LONGCALLD_ALT_ALLELE2 2

// define noisy region parameters:
#define LONGCALLD_NOISY_REG_MAX_XGAPS 5 //5 // or 10; dense X/gap region: more than n X/gap bases in a 100-bp window
#define LONGCALLD_NOISY_REG_HIFI_SLIDE_WIN 100 //200
#define LONGCALLD_NOISY_REG_ONT_SLIDE_WIN 25
#define LONGCALLD_MAX_NOISY_FRAC_PER_READ 0.5 // skip reads with more than 50% bases in noisy region
#define LONGCALLD_MAX_VAR_RATIO_PER_READ 1.00 // skip reads with n_var / ref_span > 5% 
#define LONGCALLD_MAX_READ_DEPTH 500 // vars with >500 reads will be skipped
#define LONGCALLD_MAX_NOISY_REG_READS 1000 // regions with >1000 reads will be skipped
#define LONGCALLD_NOISY_END_CLIP 30 // 100 // >= 100 bp clipping on both ends will be considered as long clipping
#define LONGCALLD_NOISY_END_CLIP_WIN 100 // 100 bp next to the long end-clipping will be considered as low-quality region
#define LONGCALLD_NOISY_REG_MERGE_DIS 500 // 500 bp, merge noisy regions within 500 bp
#define LONGCALLD_NOISY_REG_FLANK_LEN 10 // during re-alignment, include 10-bp flanking region for both ends of noisy region

#define LONGCALLD_MAX_NOISY_REG_LEN 50000 // >50kb noisy region will be skipped
#define LONGCALLD_NOISY_REG_READS 2 // >= 5 reads supporting noisy region
// #define LONGCALLD_NOISY_REG_RATIO 0.20 // >= 25% reads supporting noisy region

#define LONGCALLD_MIN_HAP_FULL_READS 1 // >= 1full read supporting each haplotype
#define LONGCALLD_MIN_HAP_READS 2 // >= 3 reads supporting each haplotype, including partial/clipped reads, call consensus from >= 3 reads
// #define LONGCALLD_MIN_NO_HAP_FULL_READS 10 // >10 total full reads in noisy region
#define LONGCALLD_MIN_READ_TO_HAP_CONS_SIM 0.9 // for reads with >= 90% equal bases, assign haplotype
// SV related parameters
#define LONGCALLD_MIN_SV_LEN 30 // size >= 30 bp -> SV
#define LONGCALLD_MIN_TSD_LEN 2 // TSD >= 2 bp
#define LONGCALLD_MAX_TSD_LEN 100 // TSD <= 100 bp
#define LONGCALLD_MIN_POLYA_LEN 10 // polyA >= 10 bp
#define LONGCALLD_MIN_POLYA_RATIO 0.8 // polyA >= 80% of the total length
// SOMATIC/MOSAIC Var
#define LONGCALLD_MIN_SOMATIC_DIS_TO_VAR 5
#define LONGCALLD_MIN_SOMATIC_DIS_TO_HP_INDEL_ERROR 3
#define LONGCALLD_MIN_SOMATIC_DIS_TO_INDEL_VAR 5
#define LONGCALLD_MIN_SOMATIC_DIS_TO_SEQ_ERROR 5
// XXX MOSAIC (low) vs SOMATIC (high)
// #define LONGCALLD_MIN_SOMATIC_AF 0.01
// #define LONGCALLD_MAX_SOMATIC_ALT_AF 0.1
#define LONGCALLD_SOMATIC_WIN 1000
#define LONGCALLD_SOMATIC_WIN_MAX_VARS 5 // >= 2 somatic vars in a 1000-bp window, consider as artifact, tag both var & read as artifact
#define LONGCALLD_MIN_SOMATIC_HAP_READS 5 // >= 10 reads supporting each haplotype
#define LONGCALLD_MIN_SOMATIC_ALT_DP 2 // >= 2 reads supporting somatic variant
#define LONGCALLD_MIN_SOMATIC_TE_ALT_DP 1 // >= 1 read supporting somatic TE variant
// fisher is NOT used for PacBio-HiFi
#define LONGCALLD_MIN_SOMATIC_FISHER_PVAL 0.05
#define LONGCALLD_STRAND_BIAS_PVAL_ONT 0.01
// BLT50: 1/10/3 works better
// #define LONGCALLD_SOMATIC_BETA_ALPHA 2 // beta prior for somatic variant calling
// #define LONGCALLD_SOMATIC_BETA_BETA 10 // beta prior for somatic variant calling
// #define LONGCALLD_MIN_SOMATIC_LOG_BETA_BINOM -4 // min log10(p-value) for beta-binomial test


// for sdust
#define LONGCALLD_SDUST_T 5  // 10
#define LONGCALLD_SDUST_W 20 // 50

// for math_utils
#define LONGCALLD_LGAMMA_MAX_I 500




#ifdef __cplusplus
extern "C" {
#endif

// for each vartiant:
// #site-level:
// * QUAL = -10 * log(1-p) (int), site-wise, p = P(data|ref)
// * FILTER: PASS, LowQual, RefCall, NoCall
// * INFO: DP/AD/AF/END
// #sample-level:
// * FORMAT: GT/GQ/DP/AD/PL/PS 
//     GP = P(G|D)
//     GL = P(D|G)
//   * PL (int) = -10 * log10(P(D|G)), Phred-scaled likelihoods of the possible genotypes.
//        "Normalized": the PL of the most likely genotype is set to 0, and the rest are scaled relative to this.
//   * GQ (int) = PL_sec - PL_lowest (int)
//   * GT = argmin(PL)
typedef struct {
    uint8_t type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    hts_pos_t pos, PS; // phase set
    uint8_t *ref_bases; int ref_len;
    uint8_t is_somatic;
    uint8_t *tsd_seq; int tsd_len, polya_len; hts_pos_t tsd_pos1, tsd_pos2; // target site duplication, 2 TSDs for DEL
    int te_seq_i, te_is_rev; // char *rep_name, *rep_family, *rep_class; use te_seq_i to retrieve TE sequence info

    int n_alt_allele; // currently, only one alt allele
    uint8_t **alt_bases; int *alt_len;
    int DP, AD[2]; uint8_t GT[2]; // DP/AD/GT
    int QUAL, GQ, PL[6]; // phred-scaled, QUAL/FILTER/GQ/PL
} var1_t;

typedef struct var_t {
    int n, m;
    var1_t *vars;
} var_t;

typedef struct call_var_opt_t {
    // input
    char *ref_fa_fn, *ref_fa_fai_fn;
    char *reg_bed_fn;
    char *in_bam_fn; char *sample_name;
    uint8_t is_pb_hifi, is_ont; float strand_bias_pval; // for ONT reads
    // variant calling regions
    uint8_t only_autosome, only_autosome_XY;
    char **exc_tnames; int n_exc_tnames;
    // char *region_list; uint8_t region_is_file; // for -R/--region/--region-file option
    // filters for variant calling
    int max_ploid, min_mq, min_bq, min_dp, min_alt_dp;
    double min_af, max_af;
    // somatic/mosaic variant
    int min_somatic_dis_to_var, min_somatic_dis_to_homopolymer_indel_error, min_somatic_dis_to_seq_error, min_somatic_alt_dp, min_somatic_te_dp, min_somatic_hap_dp;
    double min_somatic_log_beta_binom, min_somatic_fisher_pval; //, max_somatic_alt_af; // 0.01~0.1
    int somatic_beta_alpha, somatic_beta_beta; // beta prior for somatic variant calling
    int somatic_win, somatic_win_max_vars; // somatic variant calling window size, max vars in the window

    int noisy_reg_max_xgaps, noisy_reg_slide_win;
    int end_clip_reg, end_clip_reg_flank_win;
    int noisy_reg_merge_dis, noisy_reg_flank_len; // noisy_reg_merge_win; // for re-alignment
    // filters for noisy region, i.e., coverage/ratio
    int max_noisy_reg_reads, max_noisy_reg_len, min_noisy_reg_reads; 
    double max_var_ratio_per_read, max_noisy_frac_per_read; //, min_noisy_reg_ratio;
    int min_hap_full_reads, min_hap_reads; //, min_no_hap_full_reads;
    // alignment
    int match, mismatch, gap_open1, gap_ext1, gap_open2, gap_ext2;
    int gap_aln; // default: 1: left (minimap2, abpoa), 2: right (wfa2)
    double min_read_to_hap_cons_sim;
    int disable_read_realign; // disable re-alignment/MSA, only use variant calling from consensus sequence
                              // phasing is based on read cluster of the current haplotype, so not error-robust
    // TSD & polyA/T
    int min_tsd_len, max_tsd_len, min_polya_len; float min_polya_ratio;
    // Alu/L1/SVA sequences
    char *te_seq_fn; char **te_seq_names; int n_te_seqs;
    int te_kmer_len; kmer32_hash_t **te_for_h, **te_rev_h;
    // general
    // int max_ploidy;
    int pl_threads, n_threads;
    // math utils
    double lgamma_cache[LONGCALLD_LGAMMA_MAX_I+1]; int min_lgamma_i, max_lgamma_i; // 0, 999

    // output
    int min_sv_len; // classify as SV if length >= min_sv_len (50)
    htsFile *out_bam; uint8_t out_is_cram; uint8_t refine_bam; // phased bam
    htsFile *out_vcf; bcf_hdr_t *vcf_hdr; char *out_vcf_fn; char out_vcf_type; // u/b/v/z
    double p_error, log_p, log_1p, log_2; int max_gq; int max_qual;
    int8_t no_vcf_header, out_amb_base, out_somatic, out_methylation;
} call_var_opt_t;

typedef struct {
    int n_regions, m_regions;
    int *reg_tids; hts_pos_t *reg_begs, *reg_ends;
} reg_chunks_t;

typedef struct call_var_io_aux_t {
    faidx_t *fai;
    int hts_fmt; samFile *bam; bam_hdr_t *header; hts_idx_t *idx;
} call_var_io_aux_t; // per thread

// shared data for all threads
typedef struct call_var_pl_t {
    // input files
    call_var_io_aux_t *io_aux;
    // TODO mv aligner_obj (WFA/abPOA) here
    // parameters, output files
    struct call_var_opt_t *opt;
    // m-threads
    // struct bam_chunk_t *last_chunk; int n_last_chunk_reads, *last_chunk_read_i; hts_pos_t cur_active_reg_beg;
    // int max_reads_per_chunk, 
    int min_reg_chunks_per_run, max_reg_len_per_chunk;
    int reg_chunk_i, n_reg_chunks, m_reg_chunks; reg_chunks_t *reg_chunks;
    int n_threads;
} call_var_pl_t;

struct bam_chunk_t;
// separately/parallelly processed
typedef struct call_var_step_t {
    struct call_var_pl_t *pl;
    int n_chunks, max_chunks;
    struct bam_chunk_t *chunks; // input, size: n_chunks
    var_t *vars; // output, size: n_chunks
} call_var_step_t;

int call_var_main(int argc, char *argv[]);
void var_free(var_t *v);
void var1_free(var1_t *v);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_call_var_H