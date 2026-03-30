#!/bin/sh
set -eu

repo_dir=$(CDPATH= cd -- "$(dirname "$0")/.." && pwd)
tmp_dir=$(mktemp -d "${TMPDIR:-/tmp}/longcalld-trusted-input.XXXXXX")
trap 'rm -rf "$tmp_dir"' EXIT HUP INT TERM

for tool in bgzip tabix; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "Missing required tool: $tool" >&2
        exit 1
    fi
done

phased_vcf_gz="$tmp_dir/HG002_chr11_hifi_test.phased.vcf.gz"
baseline_vcf="$tmp_dir/baseline.vcf"
trusted_vcf="$tmp_dir/trusted.vcf"
baseline_sites="$tmp_dir/baseline.sites"
trusted_sites="$tmp_dir/trusted.sites"

bgzip -c "$repo_dir/test_data/HG002_chr11_hifi_test.phased.vcf" > "$phased_vcf_gz"
tabix -f -p vcf "$phased_vcf_gz"

"$repo_dir/bin/longcallD" call --hifi -V 1 \
    -o "$baseline_vcf" \
    "$repo_dir/test_data/chr11_2M.fa" \
    "$repo_dir/test_data/HG002_chr11_hifi_test.bam" \
    chr11:1300000-1301000

"$repo_dir/bin/longcallD" call --hifi -V 1 \
    --trust-input-phasing \
    --input-phased-vcf "$phased_vcf_gz" \
    -o "$trusted_vcf" \
    "$repo_dir/test_data/chr11_2M.fa" \
    "$repo_dir/test_data/HG002_chr11_hifi_test.bam" \
    chr11:1300000-1301000

if [ "$(grep -c '^##fileformat=' "$trusted_vcf")" -ne 1 ]; then
    echo "Expected exactly one VCF fileformat header line" >&2
    exit 1
fi

awk 'BEGIN{n=0}!/^#/{n++}END{exit(n==2?0:1)}' "$trusted_vcf" || {
    echo "Expected exactly 2 trusted-mode variants" >&2
    exit 1
}

awk -F '\t' '!/^#/ { print $1 "\t" $2 "\t" $4 "\t" $5 }' "$baseline_vcf" > "$baseline_sites"
awk -F '\t' '!/^#/ { print $1 "\t" $2 "\t" $4 "\t" $5 }' "$trusted_vcf" > "$trusted_sites"
cmp -s "$baseline_sites" "$trusted_sites" || {
    echo "Trusted input phasing changed the called variant sites" >&2
    exit 1
}

awk -F '\t' '
BEGIN {
    ok = 1;
    seen_1300719 = 0;
    seen_1300805 = 0;
}
/^#/ { next }
{
    split($10, fmt, ":");
    gt = fmt[1];
    ps = fmt[5];
    if ($2 == 1300719) {
        seen_1300719 = 1;
        if ($4 != "C" || $5 != "T" || gt != "0|1" || ps != "204100") ok = 0;
    } else if ($2 == 1300805) {
        seen_1300805 = 1;
        if ($4 != "T" || $5 != "G" || gt != "0|1" || ps != "204100") ok = 0;
    } else {
        ok = 0;
    }
}
END {
    exit(ok && seen_1300719 && seen_1300805 ? 0 : 1);
}
' "$trusted_vcf" || {
    echo "Trusted input phasing did not preserve the expected GT/PS assignments" >&2
    exit 1
}

echo "trusted-input-phasing test passed"
