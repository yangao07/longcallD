#!/bin/bash
# Usage: $0 step (1/2/3)
# if [[ $# -ne 1 ]]; then
  # echo "Usage: $0 <step>"
  # exit 1
# fi
# step=$1
deepvariant_sif=/hlilab/yangao/software/deepvariant/deepvariant_1.8.0-cpu.sif
version=1_8_0

THREADS=16 # 16G per thread
model=PACBIO
INPUT_DIR=/hlilab/yangao/data/HG002/deepvariant/test
OUTPUT_DIR=/hlilab/yangao/data/HG002/deepvariant/test
in_bam=${INPUT_DIR}/HG003.GRCh38.chr20.pFDA_truthv2.bam
ref_fa=${INPUT_DIR}/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna
mkdir -p ${OUTPUT_DIR} 2> /dev/null
module load singularity
# 1st run: deepvariant

output1_vcf=${OUTPUT_DIR}/output1.vcf.gz
ulimit -u 10000
singularity exec \
    -B ${INPUT_DIR} \
    ${deepvariant_sif} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${ref_fa} \
    --reads ${in_bam} \
    --output_vcf ${output1_vcf} \
    --regions chr20 \
    --num_shards ${THREADS}

# 2nd run: whatshap

phased_vcf=${OUTPUT_DIR}/output1.phased.vcf.gz
haplotagged_bam=${OUTPUT_DIR}/haplotagged.bam
whatshap phase \
        --output ${phased_vcf} \
        --reference ${ref_fa} \
        --chromosome chr20 \
        ${output1_vcf} \
        ${in_bam}
tabix -p vcf ${phased_vcf}
whatshap haplotag \
      --output ${haplotagged_bam} \
      --reference ${ref_fa} \
      ${phased_vcf} \
      ${in_bam}
samtools index ${haplotagged_bam}

# 3rd run: deepvariant
ulimit -u 10000
output2_vcf=${OUTPUT_DIR}/output2.vcf.gz
singularity run \
    -B ${INPUT_DIR},${OUTPUT_DIR} \
    ${deepvariant_sif} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${ref_fa} \
    --reads ${haplotagged_bam} \
    --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true,channel_list=haplotype" \
    --output_vcf ${output2_vcf} \
    --regions chr20 \
    --num_shards ${THREADS}

# 4th run: whatshap
phased_vcf2=${OUTPUT_DIR}/output2.phased.vcf.gz
whatshap phase \
        --output ${phased_vcf2} \
        --reference ${ref_fa} \
        ${output2_vcf} \
        ${haplotagged_bam}