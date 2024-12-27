#!/bin/bash
vcfdist=/homes2/yangao/software/miniconda3/envs/vcfdist/bin/vcfdist

S=0
N=3
threads=16
in_vcfs=(
/homes2/yangao/data/HG002/lc.1209.all.sorted.SV.vcf.gz
/homes2/yangao/data/HG002/sniffles/sniffles_withTRF_fixed.hiphased.pass.vcf.gz
/homes2/yangao/data/HG002/sawfish/call/genotyped.sv.pass.vcf.gz
)

out_dirs=(
/homes2/yangao/data/HG002/lc.1209.all.SV_vdfdist_sawfish/
/homes2/yangao/data/HG002/sniffles/sniffles_withTRF_pass_vcfdist_sawfish/
/homes2/yangao/data/HG002/sawfish/call/sawfish_pass_vcfdist_sawfish/
)

bench_vcf=/homes2/yangao/data/HG002/T2T_v1.1/sawfish_bench/GRCh38_HG2-T2TQ100-V1.0.SV.vcf.gz
bench_bed=/homes2/yangao/data/HG002/T2T_v1.1/sawfish_bench/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark_1_22.bed
ref_fa=/homes2/yangao/data/HG002/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta #chr11/chr11.fa

max_sv_len=10000
for i in $(seq $S $((N-1))); do
    in_vcf=${in_vcfs[$i]}
    out_dir=${out_dirs[$i]}
    echo "$vcfdist $in_vcf $bench_vcf $ref_fa --largest-variant $max_sv_len -b $bench_bed -p $out_dir"
done




