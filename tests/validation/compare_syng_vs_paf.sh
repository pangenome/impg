#!/usr/bin/env bash
#
# Compare syng-based and PAF-based impg queries on the same region.
#
# Both commands are expected to produce bed-format output with columns:
#   genome  start  end  strand
#
# Usage:
#   compare_syng_vs_paf.sh <syng_bed> <paf_bed>
#
# Exits 0 regardless — this is a diagnostic tool, not a test.

set -euo pipefail

syng_bed="${1:?syng bed path}"
paf_bed="${2:?paf bed path}"

# Normalize: strip trailing whitespace, drop blank lines, sort by (genome, start).
norm() {
    awk -F'\t' 'NF>=3 { printf "%s\t%s\t%s\t%s\n", $1, $2, $3, (NF>=4 ? $4 : "+") }' "$1" \
        | sort -k1,1 -k2,2n
}

syng_norm=$(mktemp)
paf_norm=$(mktemp)
trap 'rm -f "$syng_norm" "$paf_norm"' EXIT

norm "$syng_bed" > "$syng_norm"
norm "$paf_bed" > "$paf_norm"

echo "=== Row counts ==="
printf "  syng: %d\n" "$(wc -l < "$syng_norm")"
printf "  paf:  %d\n" "$(wc -l < "$paf_norm")"

echo ""
echo "=== Coverage ==="
syng_keys=$(awk -F'\t' '{print $1}' "$syng_norm" | sort -u)
paf_keys=$(awk -F'\t' '{print $1}' "$paf_norm" | sort -u)
n_syng_only=$(comm -23 <(echo "$syng_keys") <(echo "$paf_keys") | grep -c . || true)
n_paf_only=$(comm -13 <(echo "$syng_keys") <(echo "$paf_keys") | grep -c . || true)
n_both=$(comm -12 <(echo "$syng_keys") <(echo "$paf_keys") | grep -c . || true)
printf "  syng-only paths: %d\n" "$n_syng_only"
printf "  paf-only  paths: %d\n" "$n_paf_only"
printf "  both:            %d\n" "$n_both"

if [[ "$n_syng_only" -gt 0 ]]; then
    echo ""
    echo "=== syng-only paths (first 5) ==="
    comm -23 <(echo "$syng_keys") <(echo "$paf_keys") | head -5
fi
if [[ "$n_paf_only" -gt 0 ]]; then
    echo ""
    echo "=== paf-only paths (first 5) ==="
    comm -13 <(echo "$syng_keys") <(echo "$paf_keys") | head -5
fi

echo ""
echo "=== Boundary deltas (syng - paf) on common paths ==="
# For paths that appear in both, join on name and report start/end deltas.
# When multiple rows per path exist (paralogs / inversions), take the first.
join -t$'\t' -1 1 -2 1 \
    <(awk -F'\t' '!seen[$1]++ {print}' "$syng_norm") \
    <(awk -F'\t' '!seen[$1]++ {print}' "$paf_norm") \
    | awk -F'\t' 'BEGIN{n=0}
        {
            ds = $2 - $5
            de = $3 - $6
            print $1, ds, de, $4, $7
            n++
            sumds += ds; sumde += de
            abs_ds = ds < 0 ? -ds : ds
            abs_de = de < 0 ? -de : de
            if (abs_ds > max_abs_ds) max_abs_ds = abs_ds
            if (abs_de > max_abs_de) max_abs_de = abs_de
            if (abs_ds <= 5) within5_s++
            if (abs_de <= 5) within5_e++
            if (abs_ds <= 20) within20_s++
            if (abs_de <= 20) within20_e++
            if ($4 != $7) strand_mismatch++
        }
        END {
            if (n == 0) {
                print "(no common paths)"
                exit
            }
            printf "\n"
            printf "  N common rows (first per path): %d\n", n
            printf "  mean start delta: %.1fbp\n", sumds/n
            printf "  mean end   delta: %.1fbp\n", sumde/n
            printf "  max |start delta|: %dbp\n", max_abs_ds
            printf "  max |end   delta|: %dbp\n", max_abs_de
            printf "  |start delta| <= 5bp: %d/%d (%.1f%%)\n", within5_s, n, 100.0*within5_s/n
            printf "  |end   delta| <= 5bp: %d/%d (%.1f%%)\n", within5_e, n, 100.0*within5_e/n
            printf "  |start delta| <= 20bp: %d/%d (%.1f%%)\n", within20_s, n, 100.0*within20_s/n
            printf "  |end   delta| <= 20bp: %d/%d (%.1f%%)\n", within20_e, n, 100.0*within20_e/n
            printf "  strand mismatches: %d\n", strand_mismatch+0
        }' OFS=$'\t'

echo ""
echo "=== Largest boundary deltas (top 10 outliers by |start delta|) ==="
join -t$'\t' -1 1 -2 1 \
    <(awk -F'\t' '!seen[$1]++ {print}' "$syng_norm") \
    <(awk -F'\t' '!seen[$1]++ {print}' "$paf_norm") \
    | awk -F'\t' '{
        ds = $2 - $5; de = $3 - $6
        abs_ds = ds < 0 ? -ds : ds
        printf "%d\t%s\tsyng=%s..%s(%s)\tpaf=%s..%s(%s)\n", abs_ds, $1, $2, $3, $4, $5, $6, $7
    }' OFS=$'\t' \
    | sort -k1,1 -rn \
    | head -10
