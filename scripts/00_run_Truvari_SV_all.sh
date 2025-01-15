#!/bin/bash
#truvari=~/.conda/envs/sniffles/bin/truvari

S=0
N=1
in_vcfs=(
/homes2/yangao/data/HG002/lc_1227/old.SV.vcf.gz
/homes2/yangao/data/HG002/lc/lc.0114.SV.vcf.gz
/homes2/yangao/mydata/lc.1227.SV.vcf.gz
/homes2/yangao/data/HG002/sniffles/pacbio/sniffles_withTRF_fixed.hiphased.vcf.gz
/homes2/yangao/data/HG002/sawfish/call/genotyped.sv.vcf.gz
/homes2/yangao/data/HG002/cuteSV/pacbio/cutesv_hiphased.vcf.gz
)

out_dirs=(
/homes2/yangao/data/HG002/lc_1227/old.SV_truvari
/homes2/yangao/data/HG002/lc/lc.0114.SV_truvari
# /homes2/yangao/mydata/lc.1209.SV_truvari
/homes2/yangao/mydata/lc.1227.SV_truvari
/homes2/yangao/data/HG002/sniffles/pacbio/sniffles_withTRF_truvari
/homes2/yangao/data/HG002/sawfish/call/sawfish_truvari
/homes2/yangao/data/HG002/cuteSV/pacbio/cutesv_truvari
)

bench_vcf=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1.SV.vcf.gz
bench_bed=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_stvar.benchmark.bed
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta

threads=16
for((i=$S;i<$N;i++)); do
    in_vcf=${in_vcfs[$i]}
    # echo "bcf_extract_precise.sh $in_vcf"
    # bcf_extract_precise.sh $in_vcf
    # in_vcf=${in_vcf%.vcf.gz}.precise.vcf.gz
    out_dir=${out_dirs[$i]}
    # echo "truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf  --passonly --pick ac --dup-to-ins"
    # truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf  --passonly --pick ac --dup-to-ins
    echo "truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --pick ac"
    truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --pick ac
    echo "truvari refine -m '--auto --thread ${threads}' -f $ref_fa --regions ${out_dir}/candidate.refine.bed --recount --use-region-coords --use-original-vcfs --align mafft ${out_dir}"
    truvari refine -m '--auto --thread 16' -f $ref_fa --regions ${out_dir}/candidate.refine.bed --recount --use-region-coords --use-original-vcfs --align mafft ${out_dir} 
done
