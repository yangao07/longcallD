#!/bin/bash
src=/homes2/yangao/data/HG002/run_rtg_eval.sh

q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.0929.vcf.gz
t_vcf=/homes2/yangao/data/HG002/V4.2.1/chr11_10M_50M.het.phased.vcf.gz
t_bed=/homes2/yangao/data/HG002/V4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
out=/homes2/yangao/data/HG002/chr11/v4.2.1_10M_50M_het_eval2
rm -rf $out 2> /dev/null
bash $src $t_vcf $t_bed $q_vcf $out
# 
# q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.snps.vcf.gz
# t_vcf=/homes2/yangao/data/HG002/V4.2.1/chr11_10M_50M.het.snps.phased.vcf.gz
# t_bed=/homes2/yangao/data/HG002/V4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
# out=/homes2/yangao/data/HG002/chr11/v4.2.1_10M_50M_het_snps_eval
# rm -rf $out 2> /dev/null
# bash $src $t_vcf $t_bed $q_vcf $out
# 
# q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.indels.vcf.gz
# t_vcf=/homes2/yangao/data/HG002/V4.2.1/chr11_10M_50M.het.indels.phased.vcf.gz
# t_bed=/homes2/yangao/data/HG002/V4.2.1/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
# out=/homes2/yangao/data/HG002/chr11/v4.2.1_10M_50M_het_indels_eval
# rm -rf $out 2> /dev/null
# bash $src $t_vcf $t_bed $q_vcf $out

q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.vcf.gz
t_vcf=/homes2/yangao/data/HG002/T2T_v1.1/chr11_10M_50M.het.vcf.gz
t_bed=/homes2/yangao/data/HG002/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
out=/homes2/yangao/data/HG002/chr11/T2T_v1.1_10M_50M_het_eval
#rm -rf $out 2> /dev/null
#echo "bash $src $t_vcf $t_bed $q_vcf $out"
#bash $src $t_vcf $t_bed $q_vcf $out

q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.snps.vcf.gz
t_vcf=/homes2/yangao/data/HG002/T2T_v1.1/chr11_10M_50M.het.snps.vcf.gz
t_bed=/homes2/yangao/data/HG002/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
out=/homes2/yangao/data/HG002/chr11/T2T_v1.1_10M_50M_het_snps_eval
#rm -rf $out 2> /dev/null
#echo "bash $src $t_vcf $t_bed $q_vcf $out"
#bash $src $t_vcf $t_bed $q_vcf $out

q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc.indels.vcf.gz
t_vcf=/homes2/yangao/data/HG002/T2T_v1.1/chr11_10M_50M.het.indels.vcf.gz
t_bed=/homes2/yangao/data/HG002/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
out=/homes2/yangao/data/HG002/chr11/T2T_v1.1_10M_50M_het_indels_eval
#rm -rf $out 2> /dev/null
#echo "bash $src $t_vcf $t_bed $q_vcf $out"
#bash $src $t_vcf $t_bed $q_vcf $out

# # whatshap
# out=V4.2.1_whatshap.out
# q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc2.vcf.gz
# t_vcf=/homes2/yangao/data/HG002/V4.2.1/chr11_10M_50M.het.phased.vcf.gz
# whatshap compare $q_vcf $t_vcf > $out
# 
# out=DeepVar_whatshap.out
# q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc2.vcf.gz
# t_vcf=/homes2/yangao/data/HG002/DeepVariant/HG002.GRCh38.deepvariant.chr11_10M_50M.hetm2.phased.vcf.gz
# whatshap compare $q_vcf $t_vcf > $out
# 
# out=T2T_v1.1_whatshap.out
# q_vcf=/homes2/yangao/data/HG002/chr11/chr11_10M_50M.lc2.vcf.gz
# t_vcf=/homes2/yangao/data/HG002/T2T_v1.1/chr11_10M_50M.het.vcf.gz
# whatshap compare $q_vcf $t_vcf > $out
