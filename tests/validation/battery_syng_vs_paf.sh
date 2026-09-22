#!/usr/bin/env bash
#
# Battery comparison of syng-based vs PAF-based impg queries across
# many random query regions, aggregating stats so we can see the
# distribution of agreement (not just one region's results).
#
# Usage:
#   battery_syng_vs_paf.sh <syng_prefix> <paf_path> <agc_path> <n_regions> [chromosome]
#
# Emits a TSV to stdout with one row per query region:
#   region  syng_rows  paf_rows  common  syng_only  paf_only \
#     mean_ds  max_abs_ds  mean_de  max_abs_de  pct_start_5bp  pct_end_5bp

set -euo pipefail

syng_prefix="${1:?syng index prefix}"
paf_path="${2:?paf path}"
agc_path="${3:?agc path}"
n_regions="${4:?number of regions}"
chromosome="${5:-S288C#0#chrIV}"

# Determine chromosome length from the .fai sitting next to the agc's
# bgzip fasta. Fall back to 1500000 if not found.
fai_candidates=(
    "$(dirname "$agc_path")/$(basename "${agc_path%.agc}")_bgzip.fa.gz.fai"
    "$(dirname "$agc_path")/yeast235_bgzip.fa.gz.fai"
)
chr_len=0
for fai in "${fai_candidates[@]}"; do
    if [[ -f "$fai" ]]; then
        chr_len=$(awk -v c="$chromosome" '$1 == c { print $2; exit }' "$fai")
        break
    fi
done
if [[ "$chr_len" -le 0 ]]; then
    chr_len=1500000
    echo "warning: couldn't find .fai for $chromosome; using default length $chr_len" >&2
fi

# Generate region list: uniform random starts, fixed 2-5kb widths.
regions=$(awk -v n="$n_regions" -v L="$chr_len" -v c="$chromosome" \
    'BEGIN {
        srand(42)
        for (i = 0; i < n; i++) {
            w = 2000 + int(rand() * 3000)          # 2000..5000bp
            s = int(rand() * (L - w - 10000))      # keep 10kb margin
            s = (s < 0) ? 0 : s
            print c ":" s "-" (s + w)
        }
    }' | sort -u)

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT

# Header
printf "region\tsyng_rows\tpaf_rows\tcommon\tsyng_only\tpaf_only\tmean_ds\tmax_abs_ds\tmean_de\tmax_abs_de\tpct_start_5bp\tpct_end_5bp\n"

for region in $regions; do
    syng_bed="$tmpdir/$(echo "$region" | tr '/:' '__').syng.bed"
    paf_bed="$tmpdir/$(echo "$region" | tr '/:' '__').paf.bed"

    impg query -a "$syng_prefix" --sequence-files "$agc_path" \
        -r "$region" -d 10000 -o bed \
        > "$syng_bed" 2>/dev/null || { echo "syng query failed: $region" >&2; continue; }

    impg query -a "$paf_path" --sequence-files "$agc_path" \
        -r "$region" -d 10000 -o bed -x \
        2>/dev/null \
        | awk -F'\t' '{print $1"\t"$2"\t"$3"\t"$6}' \
        > "$paf_bed" || { echo "paf query failed: $region" >&2; continue; }

    syng_n=$(wc -l < "$syng_bed")
    paf_n=$(wc -l < "$paf_bed")

    # Reuse the same per-path comparison logic as compare_syng_vs_paf.sh
    # but capture summary stats only.
    common=$(join -t$'\t' -1 1 -2 1 \
        <(awk -F'\t' '!seen[$1]++ {print}' "$syng_bed" | sort) \
        <(awk -F'\t' '!seen[$1]++ {print}' "$paf_bed" | sort) \
        | wc -l)
    syng_only=$(comm -23 \
        <(awk -F'\t' '{print $1}' "$syng_bed" | sort -u) \
        <(awk -F'\t' '{print $1}' "$paf_bed" | sort -u) \
        | grep -c . || true)
    paf_only=$(comm -13 \
        <(awk -F'\t' '{print $1}' "$syng_bed" | sort -u) \
        <(awk -F'\t' '{print $1}' "$paf_bed" | sort -u) \
        | grep -c . || true)

    stats=$(join -t$'\t' -1 1 -2 1 \
        <(awk -F'\t' '!seen[$1]++ {print}' "$syng_bed" | sort) \
        <(awk -F'\t' '!seen[$1]++ {print}' "$paf_bed" | sort) \
        | awk -F'\t' '
            BEGIN { n = 0 }
            {
                ds = $2 - $5
                de = $3 - $6
                n++
                sds += ds; sde += de
                ads = ds < 0 ? -ds : ds
                ade = de < 0 ? -de : de
                if (ads > max_ads) max_ads = ads
                if (ade > max_ade) max_ade = ade
                if (ads <= 5) w5s++
                if (ade <= 5) w5e++
            }
            END {
                if (n == 0) { print "0\t0\t0\t0\t0\t0"; exit }
                printf "%.1f\t%d\t%.1f\t%d\t%.1f\t%.1f\n",
                    sds/n, max_ads+0, sde/n, max_ade+0,
                    100.0*w5s/n, 100.0*w5e/n
            }')

    printf "%s\t%d\t%d\t%d\t%d\t%d\t%s\n" \
        "$region" "$syng_n" "$paf_n" "$common" "$syng_only" "$paf_only" "$stats"
done
