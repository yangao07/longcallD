#include "assign_aln_hap.h"
#include "cgranges.h"
#include "utils.h"

int *sort_snps_by_cov(int n_cand_snps, cand_snp_t *cand_snps) {
    int *ordered_snps = (int *)malloc(n_cand_snps * sizeof(int));
    for (int i = 0; i < n_cand_snps; ++i) ordered_snps[i] = i;
    for (int i = 0; i < n_cand_snps; ++i) {
        for (int j = i + 1; j < n_cand_snps; ++j) {
            if (cand_snps[ordered_snps[i]].n_depth < cand_snps[ordered_snps[j]].n_depth) {
                int tmp = ordered_snps[i];
                ordered_snps[i] = ordered_snps[j];
                ordered_snps[j] = tmp;
            }
        }
    }
    return ordered_snps;
}

// assign haplotype to each read based on the read_snp_profile_t
// 1st round: select SNP with highest coverage, assign haplotypes to reads covering this SNP,
//            for all remaining un-assigned reads, select another highest-coverage SNP and assign haplotype, 
//            consistent with previous ones if necessary, until no more reads can be assigned
// 2~n rounds: follow the order in 1st round, re-assign haplotypes to reads based on haplotype clusters in previous rounds, 
//             until no changes to any reads
// set bam_chunk->HPs[i] to 1 or 2, 0: unknown
int assign_aln_hap(read_snp_profile_t *p, int n_cand_snps, cand_snp_t *cand_snps, bam_chunk_t *bam_chunk) {
    int *ordered_reads = (int *)malloc(bam_chunk->n_reads * sizeof(int));
    // initialize haps for all reads
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        ordered_reads[i] = i; bam_chunk->HPs[i] = -1; // -1: unassigned, 1: hap1, 2: hap2, 0: unknown
    }
    // sort SNPs by coverage
    int *ordered_snps = sort_snps_by_cov(n_cand_snps, cand_snps);
    printf("ordered_snps: ");
    for (int i = 0; i < n_cand_snps; ++i) printf("%d ", ordered_snps[i]);
    printf("\n");
    // init cr obj for read overlapping query
    cgranges_t *read_snp_cr = cr_init();
    for (int i = 0; i < bam_chunk->n_reads; ++i) {
        // [start, end): 0-based
        if (p[i].start_snp_idx < 0 || p[i].end_snp_idx < 0) continue;
        cr_add(read_snp_cr, "cr", p[i].start_snp_idx, p[i].end_snp_idx+1, i);
        printf("read_snp_cr: %d %d %d\n", p[i].start_snp_idx, p[i].end_snp_idx+1, i);
    }
    cr_index(read_snp_cr);

    // 1st round
    int n_remaining_snps = n_cand_snps, snp_i = 0;
    int64_t ovlp_i, ovlp_n, *ovlp_b = 0, max_b = 0;
    // while (n_remaining_snps > 0) {
        int cur_snp_i = ordered_snps[snp_i];
        // assign hap to cur_snp_i
        // assign_snp_hap(cand_snps+cur_snp_i);
        printf(_YELLOW("SNP: %ld (%d)") "\n", cand_snps[cur_snp_i].pos, cur_snp_i);
        // retrieve reads overlapping with current SNP
        ovlp_n = cr_overlap(read_snp_cr, "cr", cur_snp_i, cur_snp_i+1, &ovlp_b, &max_b);
        for (ovlp_i = 0; ovlp_i < ovlp_n; ++ovlp_i) {
            int read_i = cr_label(read_snp_cr, ovlp_b[ovlp_i]);
            if (bam_chunk->HPs[read_i] == -1) {
                printf("read_i: %d\n", read_i);
                // bam_chunk->HPs[read_i] = 1;
            }
        }
        snp_i++;
    // } 
    free(ovlp_b);

    // iterative updating haplotypes
    // int changed_hap = 0, read_i = 0;
    // while (1) {
    //     int read_i; // Add declaration statement
    //     for (int i = 0; i < bam_chunk->n_reads; ++i) {
    //         read_i = ordered_reads[i];
    //         int cur_hap = bam_chunk->HPs[read_i];
    //         int new_hap = update_aln_hap(p, n_cand_snps, cand_snps, read_i, bam_chunk);
    //         if (new_hap != cur_hap) {
    //             bam_chunk->HPs[read_i] = new_hap; // update intermediately
    //             changed_hap = 1;
    //         }
    //     }
    //     if (changed_hap == 0) break;
    // }

    free(ordered_snps); cr_destroy(read_snp_cr);
    return 0;
}