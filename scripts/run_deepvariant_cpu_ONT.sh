#!/bin/bash
deepvariant_sif=/hlilab/yangao/software/deepvariant/deepvariant_1.8.0-cpu.sif
version=1_8_0

ref_fa=/hlilab/11ref/hs38.fa
ref_dir=/hlilab/11ref/
THREADS=8 # 4G per thread
model=ONT_R104
in_bams=(
# /hlilab/21data/Nanopore/HG002/rep2-2023.05/PAO83395.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/PAO89685.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO83395.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO89685.pass.cram
# /hlilab/yangao/data/HG002/deepvariant/ont/case_study_github/HG002_R1041_Duplex_all_Dorado_v0.1.1_400bps_pass_2_GRCh38.chr20.bam
)
OUTPUT_DIRS=(
/hlilab/yangao/data/HG002/deepvariant/ont/PAO83395/hac
/hlilab/yangao/data/HG002/deepvariant/ont/PAO89685/hac
/hlilab/yangao/data/HG002/deepvariant/ont/PAO83395/sup
/hlilab/yangao/data/HG002/deepvariant/ont/PAO89685/sup
/hlilab/yangao/data/HG002/deepvariant/ont/case_study_github/
)

module load singularity
# 1st run: deepvariant

INPUT_DIR=/hlilab/21data/Nanopore/HG002/rep2-2023.05/
INPUT_DIR2=/hlilab/yangao/data/HG002/01bam/ont/
INPUT_DIR3=/hlilab/yangao/data/HG002/01bam/ont/sup_acc/
INPUT_DIR4=/hlilab/yangao/data/HG002/deepvariant/ont/case_study_github/
ulimit -u 10000
for((i=0;i<${#in_bams[@]};i++)); do
    in_bam=${in_bams[$i]}
    OUTPUT_DIR=${OUTPUT_DIRS[$i]}
    mkdir -p ${OUTPUT_DIR} 2> /dev/null
    output1_vcf=${OUTPUT_DIR}/output1.vcf.gz
    singularity run \
        -B ${INPUT_DIR},${INPUT_DIR2},${INPUT_DIR3},${ref_dir},${OUTPUT_DIR} \
        ${deepvariant_sif} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type ${model} \
        --ref ${ref_fa} \
        --reads ${in_bam} \
        --output_vcf ${output1_vcf} \
        --num_shards ${THREADS}

    # 2nd run: whatshap
    phased_vcf=${OUTPUT_DIR}/output1.phased.vcf.gz
    haplotagged_bam=${OUTPUT_DIR}/haplotagged.bam

    whatshap phase \
            --output ${phased_vcf} \
            --reference ${ref_fa} \
            --ignore-read-groups \
            ${output1_vcf} \
            ${in_bam}
    tabix -p vcf ${phased_vcf}

    whatshap haplotag \
          --output ${haplotagged_bam} \
          --reference ${ref_fa} \
          --ignore-read-groups \
          --skip-missing-contigs \
          ${phased_vcf} \
          ${in_bam}
    samtools index ${haplotagged_bam}

    # 3rd run: deepvariant
    output2_vcf=${OUTPUT_DIR}/output2.vcf.gz
    singularity run \
        -B ${INPUT_DIR},${INPUT_DIR2},${INPUT_DIR3},${ref_dir},${OUTPUT_DIR} \
        ${deepvariant_sif} \
        /opt/deepvariant/bin/run_deepvariant \
        --model_type ${model} \
        --ref ${ref_fa} \
        --reads ${haplotagged_bam} \
        --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true,channel_list=haplotype" \
        --output_vcf ${output2_vcf} \
        --num_shards ${THREADS}

    # 4th run: whatshap
    phased_vcf2=${OUTPUT_DIR}/output2.phased.vcf.gz
    whatshap phase \
            --output ${phased_vcf2} \
            --reference ${ref_fa} \
            --ignore-read-groups \
            ${output2_vcf} \
            ${haplotagged_bam}
    bcftools index -t ${phased_vcf2}
done