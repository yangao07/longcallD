#!/bin/bash
vcfdist=/homes2/yangao/software/miniconda3/envs/vcfdist/bin/vcfdist

if [ $# -ne 2 ]; then
    echo "Usage: $0 <in_vcf> <out_dir>"
    exit 1
fi

in_vcf=$1
out_dir=$2

bench_vcf=/homes2/yangao/data/HG002/T2T_v1.1/chr11_10M_50M.small.vcf.gz
bench_bed=/homes2/yangao/data/HG002/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_chr11_smvar.benchmark.bed
ref_fa=/homes2/yangao/data/HG002/chr11/chr11.fa

echo "$vcfdist $in_vcf $bench_vcf $ref_fa -b $bench_bed -p $out_dir -v 0"
$vcfdist $in_vcf $bench_vcf $ref_fa -b $bench_bed -p $out_dir -v 0

