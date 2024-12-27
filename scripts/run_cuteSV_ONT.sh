#!/bin/bash
cutesv=~/.conda/envs/sniffles/bin/cuteSV

in_bams=(
/homes2/yangao/data/HG002/01bam/ont/PAO83395.pass.cram
/homes2/yangao/data/HG002/01bam/ont/PAO89685.pass.cram
/homes2/yangao/data/HG002/01bam/ont/sup_acc/PAO83395.pass.cram
/homes2/yangao/data/HG002/01bam/ont/sup_acc/PAO89685.pass.cram
)
out_dirs=(
/homes2/yangao/data/HG002/cuteSV/ont/PAO83395/hac
/homes2/yangao/data/HG002/cuteSV/ont/PAO89685/hac
/homes2/yangao/data/HG002/cuteSV/ont/PAO83395/sup
/homes2/yangao/data/HG002/cuteSV/ont/PAO89685/sup
)
# sample_name=/homes2/yangao/data/HG002/sniffles/reheader_sample_name.txt
# TRF_bed=/homes2/yangao/data/TRF_bed/human_GRCh38_no_alt_analysis_set.trf.bed
ref_fa=/hlilab/11ref/hs38.fa
samp=HG002
THREADS=16
inc_bed=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed

for((i=0;i<${#in_bams[@]};i++))
do
    in_bam=${in_bams[i]}
    out_dir=${out_dirs[i]}
    out_vcf=${out_dir}/cutesv.vcf
    mkdir -p $out_dir
    rm -rf ${out_dir}/*
    ${cutesv} $in_bam $ref_fa $out_vcf $out_dir \
        --max_cluster_bias_INS 100 \
        --diff_ratio_merging_INS 0.3 \
        --max_cluster_bias_DEL 100     \
        --diff_ratio_merging_DEL 0.3 \
        --genotype \
        --sample ${samp} -include_bed ${inc_bed} \
        --threads ${THREADS}

    bs.sh $out_vcf
    # phased_vcf=${out_dir}/cutesv_hiphased.vcf.gz
    # rm ${phased_vcf}*
    # hiphase --bam $in_bam --reference $ref_fa --ignore-read-groups --vcf ${out_vcf}.gz --output-vcf $phased_vcf
done