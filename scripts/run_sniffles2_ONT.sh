#!/bin/bash
in_bams=(
/homes2/yangao/data/HG002/01bam/ont/PAO83395.pass.cram
/homes2/yangao/data/HG002/01bam/ont/PAO89685.pass.cram
/homes2/yangao/data/HG002/01bam/ont/sup_acc/PAO83395.pass.cram
/homes2/yangao/data/HG002/01bam/ont/sup_acc/PAO89685.pass.cram
)
out_dirs=(
/homes2/yangao/data/HG002/sniffles/ont/PAO83395/hac
/homes2/yangao/data/HG002/sniffles/ont/PAO89685/hac
/homes2/yangao/data/HG002/sniffles/ont/PAO83395/sup
/homes2/yangao/data/HG002/sniffles/ont/PAO89685/sup
)
sample_name=/homes2/yangao/data/HG002/sniffles/reheader_sample_name.txt
ref_fa=/hlilab/11ref/hs38.fa
TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed
THREADS=16
out_vcf=${out_dir}/sniffles_withTRF.vcf

snf=~/.conda/envs/sniffles/bin/sniffles

for((i=0;i<${#in_bams[@]};i++))
do
    in_bam=${in_bams[i]}
    out_dir=${out_dirs[i]}
    out_vcf=${out_dir}/sniffles_withTRF.vcf
    mkdir -p $out_dir
#     ${snf}  -i $in_bam \
#             --reference $ref_fa \
#             --threads ${THREADS} \
#             --tandem-repeats $TRF_bed \
#             -v $out_vcf

    rm ${out_vcf}.gz.*
    bcftools reheader -s $sample_name $out_vcf | bcftools view -Oz > ${out_vcf}.gz
    bcftools index -t ${out_vcf}.gz
    # phased_vcf=${out_dir}/sniffles_withTRF_hiphased.vcf.gz
    # hiphase --bam $in_bam --reference $ref_fa --ignore-read-groups --vcf ${out_vcf}.gz --output-vcf $phased_vcf
done