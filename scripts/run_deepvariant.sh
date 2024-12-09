#!/bin/bash
deepvariant_sif=/homes2/yangao/programs/deep_variant/deepvariant_1.6.1.sif

THREADS=16
model=PACBIO
in_bam=/hlilab/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
ref_fa=/hlilab/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
INPUT_DIR=/hlilab/yangao/data
OUTPUT_DIR=/homes2/yangao/programs/longcallD/scripts
if [[ $# -ne 0 ]]; then
  if [[ $# -ne 2 ]]; then
    echo "Usage: $0 <in_bam> <ref_fa>"
    exit 1
  fi
  in_bam=$1
  ref_fa=$2
  # OUTPUT_DIR=$3
  # INPUT_DIR=$(dirname ${in_bam})
fi


mkdir -p ${OUTPUT_DIR}/deepvariant 2> /dev/null
# 1st run: deepvariant
output1_vcf=${OUTPUT_DIR}/deepvariant/output1.vcf.gz
singularity exec \
    -B ${INPUT_DIR},${OUTPUT_DIR} \
    ${deepvariant_sif} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${ref_fa} \
    --reads ${in_bam} \
    --output_vcf ${output1_vcf} \
    # --regions chr20 \
    --num_shards ${THREADS}

# 2nd run: whatshap
phased_vcf=${OUTPUT_DIR}/deepvariant/deepvariant1.phased.vcf.gz
whatshap phase \
        --output ${phased_vcf} \
        --reference ${ref_fa} \
        # --chromosome chr20 \
        ${output1_vcf} \
        ${in_bam}

haplotagged_bam=${OUTPUT_DIR}/deepvariant/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.haplotagged.bam
whatshap haplotag \
        --output ${haplotagged_bam} \
        --reference ${ref_fa} \
        ${phased_vcf} \
        ${in_bam}

samtools index ${haplotagged_bam}

# 3rd run: deepvariant
output2_vcf=${OUTPUT_DIR}/deepvariant/output2.vcf.gz
singularity exec \
    -B ${INPUT_DIR},${OUTPUT_DIR} \
    ${deepvariant_sif} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type {model} \
    --ref ${ref_fa} \
    --reads ${haplotagged_bam} \
    --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true" \
    # --regions chr20 \
    --output_vcf ${output_vcf} \
    --num_shards ${THREADS}