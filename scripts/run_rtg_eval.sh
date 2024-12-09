#!/bin/bash
rtg=/homes2/yangao/software/rtg-tools-3.12.1/rtg

if [[ $# -ne 4 ]]; then
    echo "Usage $0 <truth.vcf.gz> <high-conf.bed> <query.vcf.gz> <output_prefix>"
    echo "
Note: vcf should be indexed with 'bcftools index -t'
"
    exit 1
fi

sdf=/homes2/yangao/data/HG002/chr11/chr11_sdf

truth_vcf=$1
truth_bed=$2
query_vcf=$3
out_pre=$4

echo "$rtg vcfeval -b $truth_vcf --bed-regions $truth_bed -c $query_vcf -o $out_pre -t $sdf"
$rtg vcfeval -b $truth_vcf --bed-regions $truth_bed -c $query_vcf -o $out_pre -t $sdf
