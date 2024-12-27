#!/bin/bash
in_bams=(
# /homes2/yangao/data/HG002/01bam/pacbio/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3_1_22_XY.bam
/homes2/yangao/data/HG002/01bam/pacbio/HG002.m84011_220902_175841_s1.GRCh38.bam
)
out_dirs=(
# /homes2/yangao/data/HG002/cuteSV/pacbio
/homes2/yangao/data/HG002/cuteSV/pacbio_sawfish_paper
)
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
samp=HG002
TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed
THREADS=16

cutesv=~/.conda/envs/sniffles/bin/cuteSV

for((i=0;i<${#in_bams[@]};i++))
do
    in_bam=${in_bams[i]}
    out_dir=${out_dirs[i]}
    out_vcf=${out_dir}/cutesv.vcf
    mkdir -p $out_dir
    rm ${out_vcf}*
    ${cutesv} $in_bam $ref_fa $out_vcf $out_dir \
        --max_cluster_bias_INS		1000 \
        --diff_ratio_merging_INS	0.9  \
        --max_cluster_bias_DEL	1000     \
        --diff_ratio_merging_DEL	0.5  \
        --genotype \
        --sample ${samp} \
        --threads ${THREADS}

    bs.sh $out_vcf
    phased_vcf=${out_dir}/cutesv_hiphased.vcf.gz
    rm ${phased_vcf}*
    hiphase --bam $in_bam --reference $ref_fa --ignore-read-groups --vcf ${out_vcf}.gz --output-vcf $phased_vcf
done