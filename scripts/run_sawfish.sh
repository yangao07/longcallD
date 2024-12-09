#!/bin/bash
THREADS=16
in_bam=/homes2/yangao/data/HG002/HG002_PacBio-HiFi-Revio_20231031_48x_GRCh38-GIABv3.bam
ref_fa=/homes2/yangao/programs/longcallD/scripts/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

dis_dir=/homes2/yangao/programs/longcallD/scripts/sawfish/discover
call_dir=/homes2/yangao/programs/longcallD/scripts/sawfish/call

# rm -rf ${dis_dir} ${call_dir}
# sawfish discover --threads ${THREADS} \
#                  --ref ${ref_fa} \
#                  --bam ${in_bam} \
#                  --output-dir ${dis_dir}

sawfish joint-call --threads ${THREADS} \
                   --sample ${dis_dir} \
                   --output-dir ${call_dir}