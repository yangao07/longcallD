#!/bin/bash
in_bams=(
# /homes2/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
/homes2/yangao/data/HG002/01bam/pacbio/HG002.m84011_220902_175841_s1.GRCh38.bam
)
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed
THREADS=16 # 4G per thread
out_dirs=(
# /homes2/yangao/data/HG002/sniffles/pacbio
/homes2/yangao/data/HG002/sniffles/pacbio_sniffles_paper
)

sample_name=/homes2/yangao/data/HG002/sniffles/reheader_sample_name.txt

for((i=0;i<${#in_bams[@]};i++)); do
    in_bam=${in_bams[$i]}
    out_dir=${out_dirs[$i]}
    mkdir -p ${out_dir}
    out_vcf=${out_dir}/sniffles_withTRF.vcf
    rm $out_vcf*
    phased_vcf=${out_dir}/sniffles_withTRF_hiphased.vcf.gz

    snf=~/.conda/envs/sniffles/bin/sniffles
    ${snf}  -i $in_bam \
            --reference $ref_fa \
            --threads ${THREADS} \
            --tandem-repeats $TRF_bed \
            -v $out_vcf

    rm ${out_vcf}.gz.*
    bcftools reheader -s $sample_name $out_vcf | bcftools view -Oz > ${out_vcf}.gz
    bcftools index -t ${out_vcf}.gz
    phased_vcf=${out_dir}/sniffles_withTRF_hiphased.vcf.gz
    rm ${phased_vcf}*
    hiphase --bam $in_bam --reference $ref_fa --ignore-read-groups --vcf ${out_vcf}.gz --output-vcf $phased_vcf
done