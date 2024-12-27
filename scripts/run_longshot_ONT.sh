#!/bin/bash
in_bams=(
/hlilab/21data/Nanopore/HG002/rep2-2023.05/PAO83395.pass.cram
/hlilab/yangao/data/HG002/01bam/ont/PAO89685.pass.cram
/hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO83395.pass.cram
/hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO89685.pass.cram
)

ref_fa=/hlilab/11ref/hs38.fa
out_dirs=(
/hlilab/yangao/data/HG002/longshot/ont/PAO83395/hac
/hlilab/yangao/data/HG002/longshot/ont/PAO89685/hac
/hlilab/yangao/data/HG002/longshot/ont/PAO83395/sup
/hlilab/yangao/data/HG002/longshot/ont/PAO89685/sup
)

ref_fa=/hlilab/11ref/hs38.fa

# TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed

ls=~/.conda/envs/sniffles/bin/longshot
for((i=0;i<${#in_bams[@]};i++)); do
    in_bam=${in_bams[$i]}
    out_dir=${out_dirs[$i]}
    mkdir -p ${out_dir} 2> /dev/null

    for chrom in $(seq 1 22) X Y; do
        out_vcf=${out_dir}/longshot_output.chr${chrom}.vcf
        echo "${ls}  --bam $in_bam --ref $ref_fa --strand_bias_pvalue_cutoff 0.01 -r chr$chrom --out $out_vcf"
    done

    # concat VCFs
    input_vcfs=""
    for chrom in $(seq 1 22) X Y; do
        input_vcfs="${input_vcfs} ${out_dir}/longshot_output.chr${chrom}.vcf"
    done

    all_vcf=${out_dir}/longshot_output.all.vcf.gz
    # echo "bcftools concat -Oz -o $all_vcf $input_vcfs"
    # echo "bcftools index -t $all_vcf"
done

