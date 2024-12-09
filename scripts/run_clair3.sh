#!/bin/bash
# run clair3 with singularity
clair3_sif=/homes2/yangao/software/clair3/clair3_v1_0_10.sif # Jul28/2024, v1.0.10
THREADS=16
MODEL_NAME=hifi_revio
PLATFORM=hifi

in_bam=/homes2/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
ref_fa=/homes2/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
INPUT_DIR=/homes2/yangao/data/HG002
OUTPUT_DIR=/homes2/yangao/data/HG002/clair3 #22/clair3_output

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

# mkdir -p ${OUTPUT_DIR} 2> /dev/null
singularity exec \
  -B ${INPUT_DIR},${OUTPUT_DIR} \
  ${clair3_sif}                 \
  /opt/bin/run_clair3.sh        \
  -b ${in_bam}                  \
  -f ${ref_fa}                  \
  -t ${THREADS}                 \
  -p ${PLATFORM}                \
  -m /opt/models/${MODEL_NAME}  \
  --use_longphase_for_intermediate_phasing \
  --use_longphase_for_final_output_phasing \
  -o ${OUTPUT_DIR}
