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
in_bam=/hlilab/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3_1_22_XY.bam
ref_fa=/hlilab/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
INPUT_DIR=/hlilab/yangao/data/HG002
OUTPUT_DIR=/hlilab/yangao/data/HG002/deepvariant/v${version}_cpu_HG002
mkdir -p ${OUTPUT_DIR} 2> /dev/null
module load singularity
# 1st run: deepvariant

output1_vcf=${OUTPUT_DIR}/output1.vcf.gz
# singularity run \
#     -B ${INPUT_DIR},${OUTPUT_DIR} \
#     ${deepvariant_sif} \
#     /opt/deepvariant/bin/run_deepvariant \
#     --model_type ${model} \
#     --ref ${ref_fa} \
#     --reads ${in_bam} \
#     --output_vcf ${output1_vcf} \
#     --num_shards ${THREADS}

# 2nd run: whatshap

phased_vcf=${OUTPUT_DIR}/output1.phased.vcf.gz
haplotagged_bam=${OUTPUT_DIR}/haplotagged.bam
# whatshap phase \
#         --output ${phased_vcf} \
#         --reference ${ref_fa} \
#         --ignore-read-groups \
#         ${output1_vcf} \
#         ${in_bam}
# tabix -p vcf ${phased_vcf}
# whatshap haplotag \
#       --output ${haplotagged_bam} \
#       --reference ${ref_fa} \
#       --ignore-read-groups \
#       --skip-missing-contigs \
#       ${phased_vcf} \
#       ${in_bam}
# samtools index ${haplotagged_bam}

# 3rd run: deepvariant
output2_vcf=${OUTPUT_DIR}/output3.vcf.gz
ulimit -u 10000
singularity run \
    -B ${INPUT_DIR},${OUTPUT_DIR} \
    ${deepvariant_sif} \
    /opt/deepvariant/bin/run_deepvariant \
    --model_type ${model} \
    --ref ${ref_fa} \
    --reads ${haplotagged_bam} \
    --make_examples_extra_args="sort_by_haplotypes=true,parse_sam_aux_fields=true,channel_list=haplotype" \
    --output_vcf ${output2_vcf} \
    --num_shards ${THREADS}

# 4th run: whatshap
phased_vcf2=${OUTPUT_DIR}/output3.phased.vcf.gz
whatshap phase \
        --output ${phased_vcf2} \
        --reference ${ref_fa} \
        ${output2_vcf} \
        ${haplotagged_bam}