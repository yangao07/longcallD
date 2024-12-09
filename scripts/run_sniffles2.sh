#!/bin/bash
in_bam=/homes2/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
ref_fa=/homes2/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta.gz
THREADS=8
out_dir=/homes2/yangao/programs/longcallD/scripts
out_vcf=${out_dir}/sniffles_out/sniffles_noTRF.vcf

snf=~/.conda/envs/sniffles/bin/sniffles
${snf}  -i $in_bam \
        --reference $ref_fa \
        --threads ${THREADS} \
        -v $out_vcf