#!/bin/bash
THREADS=16 #16G per thread
in_bams=(
# /homes2/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
/homes2/yangao/data/HG002/01bam/pacbio/HG002.m84011_220902_175841_s1.GRCh38.bam
)
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta
out_dirs=(
# /homes2/yangao/data/HG002/sawfish
/homes2/yangao/data/HG002/sawfish/pacbio_sawfish_paper
)

for i in ${!in_bams[@]}; do
    in_bam=${in_bams[$i]}
    out_dir=${out_dirs[$i]}
    mkdir -p ${out_dir}
    dis_dir=${out_dir}/discover
    call_dir=${out_dir}/call
    rm -rf ${dis_dir} ${call_dir}
    sawfish discover --threads ${THREADS} \
                     --ref ${ref_fa} \
                     --bam ${in_bam} \
                     --output-dir ${dis_dir}

    sawfish joint-call --threads ${THREADS} \
                       --sample ${dis_dir} \
                       --output-dir ${call_dir}
done