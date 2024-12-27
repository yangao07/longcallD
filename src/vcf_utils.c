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

extern int LONGCALLD_VERBOSE;

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
    fprintf(out_vcf, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">\n");
    fprintf(out_vcf, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">\n");

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

int write_var_to_vcf(var_t *vars, const struct call_var_opt_t *opt, char *chrom) {
    FILE *out_vcf = opt->out_vcf;
    int n_vars = vars->n;
    // fprintf(stdout, "n: %d\n", n_vars);
    int ret = 0;
    // XXX correct phase set if needed: set it as leftmost position within the same phase block
    for (int i = 0; i < n_vars; i++) {
        var1_t var = vars->vars[i];
        if (var.n_alt_allele == 0) continue;
        if (opt->out_amb_base == 0) {
            uint8_t skip = 0;
            for (int j = 0; j < var.ref_len; j++) {
                if (var.ref_bases[j] >= 4) {
                    if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Invalid ref base: %s %" PRId64  " %d\n", chrom, var.pos+j, var.ref_bases[j]);
                    skip = 1; break;
                }
            }
            if (skip) continue;
            for (int j = 0; j < var.n_alt_allele; j++) {
                for (int k = 0; k < var.alt_len[j]; k++) {
                    if (var.alt_bases[j][k] >= 4) {
                        if (LONGCALLD_VERBOSE >= 2) fprintf(stderr, "Invalid alt base: %s %" PRId64  " %d\n", chrom, var.pos, var.alt_bases[j][k]);
                        skip = 1; break;
                    }
                }
                if (skip) break;
            }
            if (skip) continue;
        }
        fprintf(out_vcf, "%s\t%" PRId64 "\t.\t", chrom, var.pos);
        // ref bases
        for (int j = 0; j < var.ref_len; j++) fprintf(out_vcf, "%c", "ACGTN"[var.ref_bases[j]]);
        // alt bases
        fprintf(out_vcf, "\t");
        for (int j = 0; j < var.n_alt_allele; j++) {
            for (int k = 0; k < var.alt_len[j]; k++) {
                fprintf(out_vcf, "%c", "ACGTN"[var.alt_bases[j][k]]);
            }
            if (j < var.n_alt_allele - 1) fprintf(out_vcf, ",");
        }
        int is_sv=0, k = 0;;
        char SVLEN[1024] = "SVLEN=", tmp[1024];
        char SVTYPE[1024] = "SVTYPE=";
        for (int i = 0; i < var.n_alt_allele; i++) {
            if (abs(var.alt_len[i] - var.ref_len) >= 50) { 
                if (k > 0) {
                    strcat(SVLEN, ",");
                    strcat(SVTYPE, ",");
                }
                is_sv=1;
                sprintf(tmp, "%d", var.alt_len[i] - var.ref_len);
                strcat(SVLEN, tmp);
                sprintf(tmp, "%s", var.alt_len[i] > var.ref_len ? "INS" : "DEL");
                strcat(SVTYPE, tmp);
                k++;
            }
        }
        // QUAL, FILTER, INFO
        fprintf(out_vcf, "\t%d\tPASS\tEND=%" PRId64 "", var.QUAL, var.pos + var.ref_len - 1);
        if (is_sv) fprintf(out_vcf, ";%s;%s\t", SVTYPE, SVLEN);
        else fprintf(out_vcf, "\t");
        fprintf(out_vcf, "GT:PS\t%d|%d:%" PRId64 "\n", var.GT[0], var.GT[1], var.PS);
    }
    return ret;
}