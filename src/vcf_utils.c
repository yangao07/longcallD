#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "main.h"
#include "bam_utils.h"
#include "utils.h"
#include "htslib/bgzf.h"
#include "htslib/hfile.h"
#include "collect_var.h"
#include "call_var_main.h"
#include "htslib/vcf.h"
#include "htslib/sam.h"

extern int LONGCALLD_VERBOSE;

void write_vcf_header(bam_hdr_t *hdr, struct call_var_opt_t *opt) {
    htsFile *out_vcf = opt->out_vcf;
    char *sample_name = opt->sample_name;
    bcf_hdr_t *vcf_hdr = bcf_hdr_init("w");
    if (!vcf_hdr) _err_error_exit("Could not allocate VCF header.\n");
    // File format
    // bcf_hdr_append(vcf_hdr, "##fileformat=VCFv4.3");

    // Get current date
    time_t t = time(NULL);
    struct tm *tm = localtime(&t);
    char date[11];
    strftime(date, sizeof(date), "%Y%m%d", tm);
    char date_str[50];
    snprintf(date_str, sizeof(date_str), "##fileDate=%s", date);
    bcf_hdr_append(vcf_hdr, date_str);

    // Source information
    char source_str[100];
    snprintf(source_str, sizeof(source_str), "##source=%s version=%s", PROG, LONGCALLD_VERSION);
    bcf_hdr_append(vcf_hdr, source_str);

    // Command line
    char *cmd_str = (char*)malloc(strlen(CMD) + 6);
    snprintf(cmd_str, strlen(CMD) + 6, "##CL=%s", CMD);
    bcf_hdr_append(vcf_hdr, cmd_str);
    free(cmd_str);

    // Reference sequence information
    for (int i = 0; i < hdr->n_targets; i++) {
        char contig_str[256];
        snprintf(contig_str, sizeof(contig_str), "##contig=<ID=%s,length=%d>", hdr->target_name[i], hdr->target_len[i]);
        bcf_hdr_append(vcf_hdr, contig_str);
    }

    // FILTER fields
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=PASS,Description=\"All filters passed\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=LowQual,Description=\"Low quality variant\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=RefCall,Description=\"Reference call\">");
    bcf_hdr_append(vcf_hdr, "##FILTER=<ID=NoCall,Description=\"Site has depth=0 resulting in no call\">");

    // INFO fields
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SOMATIC,Number=0,Type=Flag,Description=\"Somatic/mosaic variant\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position on CHROM\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of structural variation\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of structural variation\">");
    // TSD info
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSD,Number=1,Type=String,Description=\"Target site duplication sequence\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDLEN,Number=1,Type=Integer,Description=\"Length of target site duplication\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=POLYALEN,Number=1,Type=Integer,Description=\"Length of polyA/T sequence\">");
    // there will be 2 TSDs for DEL
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=MEI,Number=0,Type=Flag,Description=\"Mobile element insertion\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDPOS1,Number=1,Type=Integer,Description=\"Start position of first target site duplication on CHROM\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=TSDPOS2,Number=1,Type=Integer,Description=\"Start position of second target site duplication on CHROM\">");
    bcf_hdr_append(vcf_hdr, "##INFO=<ID=REPNAME,Number=1,Type=String,Description=\"Repeat name\">");
    // bcf_hdr_append(vcf_hdr, "##INFO=<ID=REPCLSFAM,Number=1,Type=String,Description=\"Repeat class/family\">");

    // FORMAT fields
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Conditional genotype quality\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Total read depth\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Read depth for each allele\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Phred-scaled genotype likelihoods rounded to the closest integer\">");
    bcf_hdr_append(vcf_hdr, "##FORMAT=<ID=PS,Number=1,Type=Integer,Description=\"Phase set\">");

    // Add sample
    bcf_hdr_add_sample(vcf_hdr, sample_name);
    bcf_hdr_add_sample(vcf_hdr, NULL); // Finalize sample list

    // Write header to file
    if (bcf_hdr_write(out_vcf, vcf_hdr) < 0) {
        fprintf(stderr, "Error: Could not write VCF header.\n");
    }
    opt->vcf_hdr = vcf_hdr;
}

int write_var_to_vcf(var_t *vars, const struct call_var_opt_t *opt, char *chrom) {
    htsFile *out_vcf = opt->out_vcf;  // htsFile pointer
    int n_vars = vars->n;
    int n_output_vars = 0;
    int buf_m = 50000; char *buffer = (char*)malloc(buf_m * sizeof(char));  // Buffer for output string
    int len = 0;        // Length of the output string
    int min_sv_len = opt->min_sv_len;

    for (int var_i = 0; var_i < n_vars; var_i++) {
        var1_t var = vars->vars[var_i];
        if (var.n_alt_allele == 0) continue;
        if (var.DP < opt->min_dp) continue;
        if (opt->out_somatic && var.is_somatic) {
            if (var.AD[1] < opt->min_somatic_te_dp) continue;
            else if (var.AD[1] < opt->min_alt_dp && var.tsd_len <= 0) continue; 
        } else {
            if (var.AD[1] < opt->min_alt_dp) continue;
        }

        // Validate bases
        if (opt->out_amb_base == 0) {
            uint8_t skip = 0;
            for (int j = 0; j < var.ref_len; j++) {
                if (var.ref_bases[j] >= 4) {
                    if (LONGCALLD_VERBOSE >= 2) 
                        fprintf(stderr, "Invalid ref base: %s %" PRId64 " %d\n", chrom, var.pos + j, var.ref_bases[j]);
                    skip = 1; break;
                }
            }
            if (skip) continue;

            for (int j = 0; j < var.n_alt_allele; j++) {
                for (int k = 0; k < var.alt_len[j]; k++) {
                    if (var.alt_bases[j][k] >= 4) {
                        if (LONGCALLD_VERBOSE >= 2) 
                            fprintf(stderr, "Invalid alt base: %s %" PRId64 " %d\n", chrom, var.pos, var.alt_bases[j][k]);
                        skip = 1; break;
                    }
                }
                if (skip) break;
            }
            if (skip) continue;
        }

        // Write to buffer
        len = snprintf(buffer, buf_m, "%s\t%" PRId64 "\t.\t", chrom, var.pos);

        // Write REF
        int base_len = var.ref_len;
        for (int j = 0; j < var.n_alt_allele; j++) base_len += var.alt_len[j];
        if (base_len > buf_m) {
            buf_m = base_len + 10000;
            buffer = (char*)realloc(buffer, buf_m * sizeof(char));
        }
        for (int j = 0; j < var.ref_len; j++) 
            len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.ref_bases[j]]);

        // Write ALT
        len += snprintf(buffer + len, buf_m - len, "\t");
        for (int j = 0; j < var.n_alt_allele; j++) {
            for (int k = 0; k < var.alt_len[j]; k++) 
                len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.alt_bases[j][k]]);
            if (j < var.n_alt_allele - 1) len += snprintf(buffer + len, buf_m - len, ",");
        }

        // Structural Variant (SV) annotation
        int is_sv = 0, k = 0;
        char SVLEN[1024] = "SVLEN=", tmp[1024];
        char SVTYPE[1024] = "SVTYPE=";
        for (int i = 0; i < var.n_alt_allele; i++) { // assert(var.n_alt_allele == 1)
            if (abs(var.alt_len[i] - var.ref_len) >= min_sv_len) { 
                if (k > 0) {
                    strcat(SVLEN, ",");
                    strcat(SVTYPE, ",");
                }
                is_sv = 1;
                sprintf(tmp, "%d", var.alt_len[i] - var.ref_len);
                strcat(SVLEN, tmp);
                sprintf(tmp, "%s", var.alt_len[i] > var.ref_len ? "INS" : "DEL");
                strcat(SVTYPE, tmp);
                k++;
            }
        }
        
        // Write QUAL, FILTER, INFO
        len += snprintf(buffer + len, buf_m - len, "\t%d\tPASS\t", var.QUAL);
        if (var.is_somatic) len += snprintf(buffer + len, buf_m - len, "SOMATIC;");
        if (var.te_seq_i >= 0) len += snprintf(buffer + len, buf_m - len, "MEI;");
        len += snprintf(buffer + len, buf_m - len, "END=%" PRId64 "", var.pos + var.ref_len - 1);
        if (is_sv) { 
            len += snprintf(buffer + len, buf_m - len, ";%s;%s", SVTYPE, SVLEN);
            if (var.tsd_len > 0) {
                len += snprintf(buffer + len, buf_m - len, ";TSD=");
                for (int i = 0; i < var.tsd_len; ++i) len += snprintf(buffer + len, buf_m - len, "%c", "ACGTN"[var.tsd_seq[i]]);
                len += snprintf(buffer + len, buf_m - len, ";TSDLEN=%d;POLYALEN=%d;TSDPOS1=%" PRId64 "", var.tsd_len, var.polya_len, var.tsd_pos1);
                if (var.tsd_pos2 > 0) len += snprintf(buffer + len, buf_m - len, ";TSDPOS2=%" PRId64 "", var.tsd_pos2);
            }
            if (var.te_seq_i >= 0) len += snprintf(buffer + len, buf_m - len, ";REPNAME=%c%s", "+-"[var.te_is_rev], opt->te_seq_names[var.te_seq_i]);
        }
        len += snprintf(buffer + len, buf_m - len, "\t");

        // Write FORMAT and Genotype Data
        int gt1 = var.GT[0], gt2 = var.GT[1];
        int is_hom = gt1 == gt2; int gt_seperator = '|';
        if (var.PS == 0) {
            gt_seperator = '/';
            if (gt1 > gt2) { int tmp = gt1; gt1 = gt2; gt2 = tmp;  }
        }
        if (is_hom || var.PS == 0) 
            len += snprintf(buffer + len, buf_m - len, "GT:DP:AD:GQ\t%d%c%d:%d:", gt1, gt_seperator, gt2, var.DP);
        else 
            len += snprintf(buffer + len, buf_m - len, "GT:DP:AD:GQ:PS\t%d%c%d:%d:", gt1, gt_seperator, gt2, var.DP);

        for (int j = 0; j < 1 + var.n_alt_allele; j++) {
            if (j > 0) len += snprintf(buffer + len, buf_m - len, ",");
            len += snprintf(buffer + len, buf_m - len, "%d", var.AD[j]);
        }

        if (is_hom || var.PS == 0) 
            len += snprintf(buffer + len, buf_m - len, ":%d\n", var.GQ);
        else 
            len += snprintf(buffer + len, buf_m - len, ":%d:%" PRId64 "\n", var.GQ, var.PS);

        // Write to htsFile
        if (out_vcf->format.compression!=no_compression) {
            if (bgzf_write(out_vcf->fp.bgzf, buffer, len) < 0) {
                _err_error_exit("Could not write to VCF file.\n");
            }
        } else {
            if (hwrite(out_vcf->fp.hfile, buffer, len) < 0) {
                _err_error_exit("Could not write to VCF file.\n");
            }
        }
        n_output_vars++;
    }
    free(buffer);
    return n_output_vars;
}
