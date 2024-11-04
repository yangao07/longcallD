#ifndef LONGCALLD_KMEDOIDS_H
#define LONGCALLD_KMEDOIDS_H

// #include "cgranges.h"

// #define bam_bseq2base(bseq, qi) seq_nt16_str[bam_seqi(bseq, qi)]

#ifdef __cplusplus
extern "C" {
#endif

int *xid_profile_2medoids(int **xid_profile, int n_dim, int n_medoids, int *haps, int n_reads, int use_hap_info);

#ifdef __cplusplus
}
#endif

#endif // end of LONGCALLD_BAM_UTILS_H