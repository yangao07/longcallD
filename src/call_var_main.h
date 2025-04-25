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
#define LONGCALLD_MIN_CAND_SOMATIC_BQ 30 // ignor bases with BQ < 30
#define LONGCALLD_MIN_CAND_DP 5 // total depth < 5: skipped
#define LONGCALLD_MIN_ALT_DP 2 // max alt depth < 2: skipped
#define LONGCALLD_MIN_SOMATIC_AF 0.05 // AF < 0.05: filtered out, 0.05~0.25: candidate somatic
#define LONGCALLD_MIN_CAND_AF 0.20 // AF < 0.20: not germline het.
#define LONGCALLD_MAX_CAND_AF 0.80 // AF > 0.80: not germline het.
#define LONGCALLD_DEF_PLOID 2 // diploid
#define LONGCALLD_REF_ALLELE 0
#define LONGCALLD_ALT_ALLELE1 1
#define LONGCALLD_ALT_ALLELE2 2

#define LONGCALLD_NOISY_REG_MAX_XGAPS 5 // or 10; dense X/gap region: more than n X/gap bases in a 100-bp window
#define LONGCALLD_NOISY_REG_SLIDE_WIN 100
#define LONGCALLD_MAX_NOISY_FRAC_PER_READ 0.5 // skip reads with more than 50% bases in noisy region
#define LONGCALLD_MAX_VAR_RATIO_PER_READ 0.05 // skip reads with n_var / ref_span > 5% 
#define LONGCALLD_MAX_READ_DEPTH 500 // vars with >500 reads will be skipped
#define LONGCALLD_MAX_NOISY_REG_READS 1000 // regions with >1000 reads will be skipped
#define LONGCALLD_NOISY_END_CLIP 100 // >= 100 bp clipping on both ends will be considered as long clipping
#define LONGCALLD_NOISY_END_CLIP_WIN 100 // 100 bp next to the long end-clipping will be considered as low-quality region
#define LONGCALLD_NOISY_REG_FLANK_LEN 10 // during re-alignment, include 10-bp flanking region for both ends of noisy region

#define LONGCALLD_MAX_NOISY_REG_LEN 50000 // >50kb noisy region will be skipped
#define LONGCALLD_NOISY_REG_READS 2 // >= 5 reads supporting noisy region
#define LONGCALLD_NOISY_REG_RATIO 0.20 // >= 25% reads supporting noisy region

#define LONGCALLD_MIN_HAP_FULL_READS 1 // >= 1full read supporting each haplotype
#define LONGCALLD_MIN_HAP_READS 2 // >= 3 reads supporting each haplotype, including partial/clipped reads, call consensus from >= 3 reads
// #define LONGCALLD_MIN_NO_HAP_FULL_READS 10 // >10 total full reads in noisy region
#define LONGCALLD_MIN_READ_TO_HAP_CONS_SIM 0.9 // for reads with >= 90% equal bases, assign haplotype
// MOSAIC SNPs
#define LONGCALLD_MIN_SOMATIC_HAP_READS 10 // >= 10 reads supporting each haplotype
#define LONGCALLD_MIN_SOMATIC_ALT_DP 2 // >= 2 reads supporting somatic variant
#define LONGCALLD_MIN_SOMATIC_TE_ALT_DP 1 // >= 1 read supporting somatic TE variant
#define LONGCALLD_MIN_SOMATIC_TE_LEN 250 // >= 250 bp INDEL to be considered as somatic TE

#define LONGCALLD_MIN_SV_LEN 30 // diff >= 30 bp: classify as SV
#define LONGCALLD_MIN_TSD_LEN 2 // TSD >= 2 bp
#define LONGCALLD_MAX_TSD_LEN 25 // TSD <= 25 bp
#define LONGCALLD_MIN_POLYA_LEN 10 // polyA >= 10 bp
#define LONGCALLD_MIN_POLYA_RATIO 0.8 // polyA >= 80% of the total length

// for sdust
#define LONGCALLD_SDUST_T 5
#define LONGCALLD_SDUST_W 20





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

typedef struct ont_hp_profile_t {
    uint8_t beg_flank_base, hp_base, end_flank_base; // 0: A, 1: C, 2: G, 3: T
    int ref_hp_len; uint8_t strand; int *alt_hp_lens; int n_alt_hp_lens;
    double *hp_len_to_prob; // -> len -> prob
} ont_hp_profile_t; // 4x4x4 * 50 * 2 = 6400

typedef struct call_var_opt_t {
    // input
    char *ref_fa_fn, *ref_fa_fai_fn;
    char *reg_bed_fn;
    char *in_bam_fn; char *sample_name;
    uint8_t is_pb_hifi, is_ont;
    // variant calling regions
    uint8_t only_autosome, only_autosome_XY;
    char **exc_tnames; int n_exc_tnames;
    // char *region_list; uint8_t region_is_file; // for -R/--region/--region-file option
    // filters for variant calling
    int max_ploid, min_mq, min_bq, min_dp, min_alt_dp;
    double min_af, max_af;
    int min_somatic_bq, min_somatic_alt_dp, min_somatic_te_dp, min_somatic_hap_dp;
    int noisy_reg_max_xgaps, noisy_reg_slide_win;
    int end_clip_reg, end_clip_reg_flank_win;
    int noisy_reg_flank_len; // noisy_reg_merge_win; // for re-alignment
    // filters for noisy region, i.e., coverage/ratio
    int max_noisy_reg_reads, max_noisy_reg_len, min_noisy_reg_reads; 
    double max_var_ratio_per_read, max_noisy_frac_per_read, min_noisy_reg_ratio;
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
    ont_hp_profile_t ***ont_hp_profile; // for ONT data, use HP profile to call variants, size: 4x4x4
    // general
    // int max_ploidy;
    int pl_threads, n_threads;
    // math utils
    double lgamma_cache[100]; int min_lgamma_i, max_lgamma_i; // 0, 999

    // output
    int min_sv_len; // classify as SV if length >= min_sv_len (50)
    htsFile *out_bam; uint8_t out_is_cram; // phased bam
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