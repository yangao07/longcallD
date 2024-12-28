#!/bin/bash
py=~/software/miniconda3/envs/hap.py/bin/python2
happy=~/software/miniconda3/envs/hap.py/bin/hap.py

S=0
N=1
in_vcfs=(
# /homes2/yangao/data/HG002/lc.1209.all.sorted.small.vcf.gz
/homes2/yangao/mydata/lc.1227.small.vcf.gz
/homes2/yangao/data/HG002/clair3/phased_merge_output.small.vcf.gz
/homes2/yangao/data/HG002/deepvariant/HG002.GRCh38.deepvariant.phased.small.vcf.gz
/homes2/yangao/data/HG002/deepvariant/v1_8_0_cpu_HG002/output3.small.vcf.gz
/homes2/yangao/data/HG002/longshot/pacbio/longshot_output.all.vcf.gz
# ONT
/homes2/yangao/data/HG002/clair3/ont/PAO83395/hac/phased_merge_output.vcf.gz
/homes2/yangao/data/HG002/clair3/ont/PAO83395/sup/phased_merge_output.vcf.gz
/homes2/yangao/data/HG002/clair3/ont/PAO89685/hac/phased_merge_output.vcf.gz
/homes2/yangao/data/HG002/clair3/ont/PAO89685/sup/phased_merge_output.vcf.gz
/homes2/yangao/data/HG002/deepvariant/ont/PAO83395/hac/output2.phased.vcf.gz
/homes2/yangao/data/HG002/deepvariant/ont/PAO83395/sup/output2.phased.vcf.gz
/homes2/yangao/data/HG002/deepvariant/ont/PAO89685/hac/output2.phased.vcf.gz
/homes2/yangao/data/HG002/deepvariant/ont/PAO89685/sup/output2.phased.vcf.gz
)

out_dirs=(
# /homes2/yangao/data/HG002/lc.1209.happy
/homes2/yangao/mydata/lc.1227.happy
/homes2/yangao/data/HG002/clair3/happy_small
/homes2/yangao/data/HG002/deepvariant/happy_small
/homes2/yangao/data/HG002/deepvariant/v1_8_0_cpu_HG002/output3.small_happy
/homes2/yangao/data/HG002/longshot/pacbio/happy_small
# ONT
/homes2/yangao/data/HG002/clair3/ont/PAO83395/hac/happy
/homes2/yangao/data/HG002/clair3/ont/PAO83395/sup/happy
/homes2/yangao/data/HG002/clair3/ont/PAO89685/hac/happy
/homes2/yangao/data/HG002/clair3/ont/PAO89685/sup/happy
/homes2/yangao/data/HG002/deepvariant/ont/PAO83395/hac/happy
/homes2/yangao/data/HG002/deepvariant/ont/PAO83395/sup/happy
/homes2/yangao/data/HG002/deepvariant/ont/PAO89685/hac/happy
/homes2/yangao/data/HG002/deepvariant/ont/PAO89685/sup/happy
)

n_threads=4

bench_vcf=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1.small.vcf.gz
bench_bed=/homes2/yangao/data/HG002/02benchmark_truth_set/T2T_v1.1/GRCh38_HG2-T2TQ100-V1.1_smvar.benchmark.bed
ref_fa=/homes2/yangao/data/HG002/00ref/GRCh38_GIABv3_no_alt_analysis_set_maskedGRC_decoys_MAP2K3_KMT2C_KCNJ18.fasta


for i in $(seq $S $((N-1))); do
    in_vcf=${in_vcfs[$i]}
    out_dir=${out_dirs[$i]}

    echo "$py $happy --pass-only --threads $n_threads -r $ref_fa -f $bench_bed -o $out_dir $bench_vcf $in_vcf"
    $py $happy --pass-only --threads $n_threads -r $ref_fa -f $bench_bed -o $out_dir $bench_vcf $in_vcf
done
