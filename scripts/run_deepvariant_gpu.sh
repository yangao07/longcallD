#!/bin/bash
# Usage: $0 step (1/2/3)
if [[ $# -ne 1 ]]; then
  echo "Usage: $0 <step>"
  exit 1
fi
step=$1
deepvariant_sif=/hlilab/yangao/software/deepvariant/deepvariant_1.8.0-gpu.sif
version=1_8_0

THREADS=8
# region="chr11:10000000-50000000"
# region=chr11
model=PACBIO
in_bam=/hlilab/yangao/data/HG002/chr11/chr11_10M_50M.bam #HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3_1_22_XY.bam
ref_fa=/hlilab/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
INPUT_DIR=/hlilab/yangao/data/HG002
OUTPUT_DIR=/hlilab/yangao/data/HG002/deepvariant/v${version}
mkdir -p ${OUTPUT_DIR} 2> /dev/null
mkdir -p ${OUTPUT_DIR}_cpu 2> /dev/null
# 1st run: deepvariant

output1_vcf=${OUTPUT_DIR}/output1.vcf.gz
if [[ $step -eq 1 ]]; then
  singularity run --nv \
      -B ${INPUT_DIR},${OUTPUT_DIR} \
      ${deepvariant_sif} \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type ${model} \
      --ref ${ref_fa} \
      --reads ${in_bam} \
      --output_vcf ${output1_vcf} \
      --num_shards ${THREADS}
fi

# 2nd run: whatshap

phased_vcf=${OUTPUT_DIR}/output1.phased.vcf.gz
haplotagged_bam=${OUTPUT_DIR}/chr11_10M_50M.haplotagged.bam
if [[ $step -eq 2 ]]; then
  whatshap phase \
          --output ${phased_vcf} \
          --reference ${ref_fa} \
          ${output1_vcf} \
          ${in_bam}
  whatshap haplotag \
        --output ${haplotagged_bam} \
        --reference ${ref_fa} \
        ${phased_vcf} \
        ${in_bam}

  samtools index ${haplotagged_bam}
fi

# 3rd run: deepvariant
output2_vcf=${OUTPUT_DIR}_cpu/output2.vcf.gz

if [[ $step -eq 3 ]]; then
  singularity run --nv \
      -B ${INPUT_DIR},${OUTPUT_DIR}_cpu \
      ${deepvariant_sif} \
      /opt/deepvariant/bin/run_deepvariant \
      --model_type ${model} \
      --ref ${ref_fa} \
      --reads ${haplotagged_bam} \
      --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true" \
      --output_vcf ${output_vcf} \
      --num_shards ${THREADS}
fi