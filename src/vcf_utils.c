#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h> // Add the missing import statement for the 'time' library
#include "main.h"
#include "bam_utils.h"
#include "collect_var.h"
#include "call_var.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"

// Write VCF header to FILE
void write_vcf_header(bam_hdr_t *hdr, FILE *out_vcf, char *sample_name) {
    fprintf(out_vcf, "##fileformat=VCFv4.3\n");
    // Get current date
    time_t t = time(NULL); struct tm *tm = localtime(&t);
    char date[11]; strftime(date, sizeof(date), "%Y%m%d", tm);
    fprintf(out_vcf, "##fileDate=%s\n", date);
    fprintf(out_vcf, "##source=%s version=%s\n", PROG, VERSION);
    fprintf(out_vcf, "##CL=%s\n", CMD);
    // write reference sequence information
    for (int i = 0; i < hdr->n_targets; i++)
        fprintf(out_vcf, "##contig=<ID=%s,length=%d>\n", hdr->target_name[i], hdr->target_len[i]);
    
    // FILTER field
    fprintf(out_vcf, "##FILTER=<ID=PASS,Description=\"All filters passed\">\n");
    fprintf(out_vcf, "##FILTER=<ID=LowQual,Description=\"Low quality variant\">\n");
    fprintf(out_vcf, "##FILTER=<ID=RefCall,Description=\"Reference call\">\n");
    fprintf(out_vcf, "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">\n");

    // INFO field
    // fprintf(out_vcf, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Combined depth across samples\">\n");
    // fprintf(out_vcf, "##INFO=<ID=AD,Number=R,Type=Integer,Description=\"Total read depth for each allele\">\n");
    // fprintf(out_vcf, "##INFO=<ID=AF,Number=1,Type=Float,Description=\"Allele Frequency\">\n");
    fprintf(out_vcf, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position on CHROM\">\n");

    // FORMAT field
    fprintf(out_vcf, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Conditional genotype quality\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">\n");
    fprintf(out_vcf, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">\n");

    // field header
    fprintf(out_vcf, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t%s\n", sample_name);
}

int write_var_to_vcf(var_t *vars, FILE *out_vcf, char *chrom) {
    int n_vars = vars->n;
    // fprintf(stdout, "n: %d\n", n_vars);
    int ret = 0;
    // XXX correct phase set if needed: set it as leftmost position within the same phase block
    for (int i = 0; i < n_vars; i++) {
        var1_t var = vars->vars[i];
        // if (var.pos == 49735305)
            // fprintf(stderr, "OK\n");
        if (var.n_alt_allele == 0) continue;
        fprintf(out_vcf, "%s\t%" PRId64 "\t.\t", chrom, var.pos);
        // ref bases
        for (int j = 0; j < var.ref_len; j++) {
            fprintf(out_vcf, "%c", "ACGTN"[var.ref_bases[j]]);
        }
        // alt bases
        fprintf(out_vcf, "\t");
        for (int j = 0; j < var.n_alt_allele; j++) {
            for (int k = 0; k < var.alt_len[j]; k++) {
                fprintf(out_vcf, "%c", "ACGTN"[var.alt_bases[j][k]]);
            }
            if (j < var.n_alt_allele - 1) fprintf(out_vcf, ",");
        }
        // QUAL, FILTER, INFO
        fprintf(out_vcf, "\t%d\tPASS\tEND=%" PRId64 "\tGT:PS\t%d|%d:%" PRId64 "\n", var.QUAL, var.pos + var.ref_len - 1, var.GT[0], var.GT[1], var.PS);
    }
    return ret;
}