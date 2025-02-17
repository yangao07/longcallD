#ifndef LONGCALLD_COLLECT_VAR_H
#define LONGCALLD_COLLECT_VAR_H

#include <math.h>
#include "htslib/sam.h"

// category of candidate variants
// #define LONGCALLD_VAR_CATE_N 8
#define LONGCALLD_VAR_CATE_STR "LENIRXSHehl"

#define LONGCALLD_LOW_COV_VAR        0x001 // "L"
// #define LONGCALLD_CLEAN_HET_VAR      0x002 // "E" // not used for now
#define LONGCALLD_CLEAN_HET_SNP      0x004 // "N"
#define LONGCALLD_CLEAN_HET_INDEL    0x008 // "I"
#define LONGCALLD_REP_HET_VAR        0x010 // "R"
#define LONGCALLD_NOISY_REG_VAR      0x020 // "X"
#define LONGCALLD_CAND_SOMA_VAR      0x040 // "S"
#define LONGCALLD_CAND_HOM_VAR       0x080 // "H"
#define LONGCALLD_NOISY_CAND_HET_VAR 0x100 // "e"
#define LONGCALLD_NOISY_CAND_HOM_VAR 0x200 // "h"
#define LONGCALLD_LOW_AF_VAR         0x400 // "l"

#define LONGCALLD_VAR_CATE_TYPE(var_cate) LONGCALLD_VAR_CATE_STR[(int)(log2(var_cate))]

#define LONGCALLD_REF_CONS_ALN_STR(clu_aln_strs) clu_aln_strs
#define LONGCALLD_CONS_READ_ALN_STR(clu_aln_strs, read_i) clu_aln_strs+read_i+1
#define LONGCALLD_REF_READ_ALN_STR(clu_aln_strs, read_i) clu_aln_strs+read_i+1
// #define LONGCALLD_CONS_READ_ALN_STR(clu_aln_strs, read_i) clu_aln_strs+(read_i+1)*2-1
// #define LONGCALLD_REF_READ_ALN_STR(clu_aln_strs, read_i) clu_aln_strs+(read_i+1)*2

#ifdef __cplusplus
extern "C" {
#endif

typedef struct var_site_t {
    int tid; hts_pos_t pos; // var_type: BAM_CINS/BAM_CDEL/BAM_CDIFF
    int var_type, ref_len, alt_len; // SNP: 1:1, INS: 0:I, DEL: D:1
    uint8_t *alt_seq; // only used for mismatch/insertion, deletion:NULL
} var_site_t;

// XXX each cand_var only have one alt_allele/var/seq, previously we keep multiple insertions in one var
// XXX for insertion/deletion, not include the first reference base, will add it when output to VCF
// alignment-based simple variant: SNP, small indel in non-repetitive region
// assume alignment is correct, no complex variants
// will be used for variant calling, haplotype assignment, phasing of variants and reads
typedef struct cand_var_t {
    // static information
    int tid; hts_pos_t pos, phase_set;
    int var_type; // BAM_CINS/BAM_CDEL/BAM_CDIFF
    int total_cov; // ref+alt, used in variant calling & haplotype assignment, excluding low-qual bases
    int low_qual_cov; // including bases/regions with low quality, only count depth, not seq
    int n_uniq_alles; // up ot 4: ref, alt1, alt2, minor_alt
                      // snp/ins: could be >2, del: ≤2
                      // minor_alt: not ref and not main alt alleles (mostly sequencing errors)
    int *alle_covs; // size: n_uniq_alles
    int ref_len; uint8_t ref_base; // 1-base ref_base, only used for X
    int alt_len; uint8_t *alt_seq; // only used for mismatch/insertion, deletion:NULL

    // dynamic information, update during haplotype assignment
    int **hap_to_alle_profile; // read-wise: 1:H1/2:H2 -> alle_i -> read count
    int *hap_to_cons_alle; // HAP-wise (hap_to_cons_alle_i): 1:H1/2:H2 -> alle_i

    uint8_t is_low_qual, is_skipped;
} cand_var_t;

// read X var
typedef struct read_var_profile_t {
    int read_id; // 0 .. bam_chunk->n_read-1 XXX for noisy region
    int start_var_idx, end_var_idx; // 0 .. n_total_cand_vars-1
    // uint8_t *var_is_used; // size: n_total_cand_vars
    int *alleles; // 0:ref, 1:alt
} read_var_profile_t;

typedef struct aln_str_t {
    uint8_t *target_aln;
    uint8_t *query_aln;
    int aln_len;
    // for partially aligned reads
    int target_beg, target_end, query_beg, query_end; // 0..aln_len-1, including '-'/5
} aln_str_t;

struct bam_chunk_t;
struct var_t;
struct call_var_pl_t;
struct call_var_step_t;
struct digar1_t;
struct cand_var_t;

int comp_var_site(var_site_t *var1, var_site_t *var2);
int comp_ovlp_var_site(var_site_t *var1, var_site_t *var2, int *is_ovlp);
var_site_t make_var_site_from_digar(int tid, struct digar1_t *digar);
var_site_t make_var_site_from_cand_var(struct cand_var_t *cand_var);
void collect_var_main(const struct call_var_pl_t *pl, struct bam_chunk_t *bam_chunk);
void stitch_var_main(struct call_var_step_t *step, struct bam_chunk_t *bam_chunk, struct var_t *var, long ii);
void free_cand_vars(struct cand_var_t *cand_vars, int m);

#ifdef __cplusplus
}
#endif

#endif // LONGCALLD_COLLECT_VAR_H
