#!/bin/bash
# run clair3 with singularity
clair3_sif=/homes2/yangao/software/clair3/clair3_v1_0_10.sif # Jul28/2024, v1.0.10
THREADS=32 # 2G per thread
PLATFORM=ont

in_bams=(
# /hlilab/21data/Nanopore/HG002/rep2-2023.05/PAO83395.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/PAO89685.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO83395.pass.cram
# /hlilab/yangao/data/HG002/01bam/ont/sup_acc/PAO89685.pass.cram
)
MODEL_NAMES=(
r1041_e82_400bps_hac_v500
r1041_e82_400bps_hac_v500
r1041_e82_400bps_sup_v500
r1041_e82_400bps_sup_v500
)
INPUT_DIR=/hlilab/21data/Nanopore/HG002/rep2-2023.05/
INPUT_DIR2=/hlilab/yangao/data/HG002/01bam/ont/
INPUT_DIR3=/hlilab/yangao/data/HG002/01bam/ont/sup_acc/
ref_fa=/hlilab/11ref/hs38.fa
ref_dir=/hlilab/11ref
OUTPUT_DIRS=(
/hlilab/yangao/data/HG002/clair3/ont/PAO83395/hac
/hlilab/yangao/data/HG002/clair3/ont/PAO89685/hac
/hlilab/yangao/data/HG002/clair3/ont/PAO83395/sup
/hlilab/yangao/data/HG002/clair3/ont/PAO89685/sup
)

clair3_dir=/hlilab/yangao/software/clair3

module load singularity
for((i=0;i<${#in_bams[@]};i++)); do
    in_bam=${in_bams[$i]}
    OUTPUT_DIR=${OUTPUT_DIRS[$i]}
    MODEL_NAME=${MODEL_NAMES[$i]}
    mkdir -p ${OUTPUT_DIR} 2> /dev/null
    singularity exec \
      -B ${INPUT_DIR},${INPUT_DIR2},${INPUT_DIR3},${ref_dir},${clair3_dir},${OUTPUT_DIR} \
      ${clair3_sif}                   \
      /opt/bin/run_clair3.sh          \
      -b ${in_bam}                    \
      -f ${ref_fa}                    \
      -t ${THREADS}                   \
      -p ${PLATFORM}                  \
      -m ${clair3_dir}/${MODEL_NAME}  \
      --use_longphase_for_intermediate_phasing \
      --use_longphase_for_final_output_phasing \
      -o ${OUTPUT_DIR}
done