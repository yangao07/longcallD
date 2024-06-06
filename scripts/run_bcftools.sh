#!/bin/bash
# run bcftools to call variants
in_bam=/home/gaoy1/sdata/variant_calling/HG002/GRCh38/chr22/chr22.bam
ref_fa=/home/gaoy1/sdata/variant_calling/HG002/GRCh38/chr22/chr22.fa
THREADS=8
OUTPUT_DIR=/home/gaoy1/sdata/variant_calling/HG002/GRCh38/chr22/bcftools_output

mkdir -p $OUTPUT_DIR
echo "bcftools mpileup -f $ref_fa $in_bam | bcftools call -mv -Ob -o ${OUTPUT_DIR}/calls.bcf"
bcftools mpileup -f $ref_fa $in_bam | bcftools call -mv -Ob -o ${OUTPUT_DIR}/calls.bcf
