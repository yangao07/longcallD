#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // Add the missing import statement for the 'time' library
#include "main.h"
#include "bam_utils.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"

// Write VCF header to FILE
void write_vcf_header(bam_hdr_t *hdr, FILE *out_vcf, char *sample_name) {
    fprintf(out_vcf, "##fileformat=VCFv4.2\n");
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    // Get current date
    time_t t = time(NULL); struct tm *tm = localtime(&t);
    char date[11]; strftime(date, sizeof(date), "%Y%m%d", tm);
    fprintf(out_vcf, "##fileDate=%s\n", date);
    fprintf(out_vcf, "##source=%s version=%s\n", PROG, VERSION);
    fprintf(out_vcf, "##CL=%s\n", CMD);
    // write reference sequence information
    for (int i = 0; i < hdr->n_targets; i++)
        fprintf(out_vcf, "##contig=<ID=%s,length=%d>\n", hdr->target_name[i], hdr->target_len[i]);
    
    // INFO field
    // fprintf(out_vcf, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");

    // FORMAT field
    fprintf(out_vcf, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total Depth\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase Set\">\n");

    // header
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sample_name);
}

int write_snp_to_vcf(cand_snp_t *cand_snps, int n_cand_snps, FILE *out_vcf, char *chrom) {
    int i;
    int ret = 0;
    // Write each SNP to VCF
    // XXX correct phase set if needed: set it as leftmost SNP position within the same phase block
    for (i = 0; i < n_cand_snps; i++) {
        cand_snp_t snp = cand_snps[i];
        if (snp.is_skipped) continue;
        // Determine the genotype based on hap_to_cons_base
        char gt[4] = "0|0"; uint8_t allele_bases[2];
        char ref_base='.', alt_bases[4]; int n_alt_bases=0;
        
        if (snp.ref_base != -1) ref_base = LONGCALLD_BAM_BASE_STR[snp.ref_base];

        for (int i = 0; i < 2; ++i) {
            int base_i = snp.hap_to_cons_base[i+1];
            if (base_i < 0 || base_i >= snp.n_uniq_bases) {
                fprintf(stderr, "Error: hap_to_cons_base[%d] = %d, out of range [0, %d)\n", i, base_i, snp.n_uniq_bases);
                ret = 1; break;
            }
            allele_bases[i] = snp.bases[base_i];
            if (allele_bases[i] == LONGCALLD_BAM_REF_BASE_IDX) {
                gt[i*2] = '0';
            } else {
                alt_bases[n_alt_bases*2] = LONGCALLD_BAM_BASE_STR[allele_bases[i]];
                if (n_alt_bases > 0) alt_bases[n_alt_bases*2-1] = ',';
                // check if alt_bases is already in the list
                int j, add_alt_base = 1;
                for (j = 0; j < n_alt_bases; j++) {
                    if (allele_bases[i] == allele_bases[j]) {
                        add_alt_base = 0;
                        break;
                    }
                }
                n_alt_bases += add_alt_base;
                gt[i*2] = n_alt_bases + '0';
            }
        }
        if (n_alt_bases == 0) continue;
        alt_bases[n_alt_bases*2-1] = '\0'; gt[3] = '\0';

        // Write SNP information to VCF
        fprintf(out_vcf, "%s\t%ld\t.\t%c\t%s\t.\t.\t.\tGT:PS\t%s:%ld\n", 
                          chrom, snp.pos, ref_base, alt_bases, gt, snp.phase_set);
    }
    return ret;
}