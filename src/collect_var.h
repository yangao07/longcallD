#ifndef LONGCALLD_COLLECT_VAR_H
#define LONGCALLD_COLLECT_VAR_H

#include <math.h>
#include "htslib/sam.h"

// category of candidate variants
// #define LONGCALLD_VAR_CATE_N 8
#define LONGCALLD_VAR_CATE_STR "LENIRXSHehl"

#define LONGCALLD_LOW_COV_VAR        0x001 // "L"
#define LONGCALLD_CLEAN_HET_VAR      0x002 // "E" // not used for now
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


#ifdef __cplusplus
extern "C" {
#endif

typedef struct var_site_t {
    int tid; hts_pos_t pos; // var_type: BAM_CINS/BAM_CDEL/BAM_CDIFF
    int var_type, ref_len, alt_len; // SNP: 1:1, INS: 0:I, DEL: D:1
    uint8_t *alt_seq; // only used for mismatch/insertion, deletion:NULL
} var_site_t;

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
