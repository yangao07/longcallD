#!/bin/bash
truvari=truvari


S=0
N=1
threads=16
in_vcfs=(
/homes2/yangao/data/HG002/longCallD/lc.1209.all.sorted.SV.vcf.gz
# /homes2/yangao/mydata/lc.1227.SV.vcf.gz
/homes2/yangao/data/HG002/sniffles/pacbio/sniffles_withTRF_fixed.hiphased.vcf.gz
/homes2/yangao/data/HG002/sawfish/call/genotyped.sv.vcf.gz
/homes2/yangao/data/HG002/sniffles/pacbio_sniffles_paper/sniffles_withTRF_hiphased.vcf.gz
/homes2/yangao/data/HG002/sawfish/pacbio_sawfish_paper/call/genotyped.sv.vcf.gz
)

out_dirs=(
# /homes2/yangao/data/HG002/longCallD/lc.1209.all.SV_truvari_sawfish
/homes2/yangao/mydata/lc.1209_SV_truvari_sawfish
# /homes2/yangao/mydata/lc.1227_SV_truvari_sawfish
/homes2/yangao/data/HG002/sniffles/pacbio/sniffles_withTRF_truvari_sawfish
/homes2/yangao/data/HG002/sawfish/call/sawfish_truvari_sawfish
/homes2/yangao/data/HG002/sniffles/pacbio_sniffles_paper/sniffles_withTRF_truvari_sawfish
/homes2/yangao/data/HG002/sawfish/pacbio_sawfish_paper/call/sawfish_truvari_sawfish
)

bench_vcf=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/sawfish_bench/GRCh38_HG2-T2TQ100-V1.0.SV.vcf.gz
bench_bed=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/sawfish_bench/GRCh38_HG2-T2TQ100-V1.0_stvar.benchmark_1_22.bed
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta #chr11/chr11.fa

for i in $(seq $S $((N-1))); do
    in_vcf=${in_vcfs[$i]}
    out_dir=${out_dirs[$i]}
    if [[ $i -ne 0 ]]; then
    echo "truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --passonly --pick ac --dup-to-ins"
    $truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --passonly --pick ac --dup-to-ins
    else
        echo "truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --pick ac"
        truvari bench -f $ref_fa --includebed $bench_bed -o $out_dir -b $bench_vcf -c $in_vcf --pick ac
    fi

    echo "truvari refine -m '--auto --thread ${threads}' -f $ref_fa --regions ${out_dir}/candidate.refine.bed --recount --use-region-coords --use-original-vcfs --align mafft ${out_dir} "
    truvari refine -m '--auto --thread 16' -f $ref_fa --regions ${out_dir}/candidate.refine.bed --recount --use-region-coords --use-original-vcfs --align mafft ${out_dir}
    # echo "truvari ga4gh --input ${out_dir} --output ${out_dir}/combined_result --with-refine"
    # truvari ga4gh --input ${out_dir} --output ${out_dir}/combined_result --with-refine
done
