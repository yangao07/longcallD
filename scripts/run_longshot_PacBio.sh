#!/bin/bash
in_bam=/homes2/yangao/data/HG002/01bam/pacbio/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
# TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed
THREADS=8
out_dir=/homes2/yangao/data/HG002/longshot
out_vcf=${out_dir}/longshot_output.vcf

ls=~/.conda/envs/sniffles/bin/longshot
for chrom in $(seq 1 22) X Y; do
    out_vcf=${out_dir}/longshot_output.chr${chrom}.vcf
    # echo "${ls}  --bam $in_bam --ref $ref_fa  -r chr$chrom --out $out_vcf"
done

# concat VCFs
input_vcfs=""
for chrom in $(seq 1 22) X Y; do
    input_vcfs="${input_vcfs} ${out_dir}/longshot_output.chr${chrom}.vcf"
done

all_vcf=${out_dir}/longshot_output.all.vcf.gz
echo "bcftools concat -Oz -o $all_vcf $input_vcfs"
# bcftools concat -Oz -o $all_vcf $input_vcfs
echo "bcftools index -t $all_vcf"
