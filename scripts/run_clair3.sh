#!/bin/bash
# run clair3 with singularity
clair3_sif=/home/gaoy1/software/clair3-sing/clair3_latest.sif
THREADS=8
MODEL_NAME=hifi_revio
PLATFORM=hifi

in_bam=~/hg002/chr1/tmp.bam # chr22/chr22.bam
ref_fa=~/hg002/chr1/chr1.fa #chr22/chr22.fa
INPUT_DIR=~/hg002/chr11 #22
OUTPUT_DIR=~/hg002/chr11/clair3_output #22/clair3_output

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


mkdir -p ${OUTPUT_DIR}
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  ${clair3_sif}                 \
  /opt/bin/run_clair3.sh        \
  -b ${in_bam}                  \
  -f ${ref_fa}                  \
  -t ${THREADS}                 \
  -p ${PLATFORM}                \
  -m /opt/models/${MODEL_NAME}  \
  --use_whatshap_for_final_output_phasing \
  --use_whatshap_for_final_output_haplotagging \
  -o ${OUTPUT_DIR}