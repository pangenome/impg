#!/bin/bash
set -euo pipefail

# Demo: compare full vs partitioned pggb on a single region.
#
#   Full:        query -o gfa --gfa-engine pggb --force-large-region
#   Partitioned: query -o gfa --gfa-engine pggb:10000
#
# Validates that both GFAs contain exactly the same set of input sequences.
#
# Usage:
#   cd /path/to/impg   (repo root)
#   bash scripts/demo-partitioned.sh

IMPG="target/release/impg"
AGC="data/hprcv2/HPRC_r2_assemblies_0.6.1.agc"
OUTDIR="demo-partitioned-results"
THREADS=12
TIME="/usr/bin/time"
ODGI=$(command -v odgi 2>/dev/null || true)

REGION="CHM13#0#chr6:29000000-30000000"
PARTITION_SIZE=10000
BATCH_BYTES="50M"

mkdir -p "$OUTDIR"

PERF_TSV="$OUTDIR/performance.tsv"
VAL_TSV="$OUTDIR/validation.tsv"

printf "Command\tWall_s\tMem_MB\tSegments\tLinks\tPaths\tAvgSegBp\tStatus\n" > "$PERF_TSV"
printf "Comparison\tResult\tDetail\n" > "$VAL_TSV"

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

run_timed() {
    local log_file="$1"; shift
    if "$TIME" -v "$@" 2>"$log_file"; then
        local status=0
    else
        local status=$?
    fi
    local wall peak_kb
    wall=$(grep "Elapsed (wall clock)" "$log_file" | sed 's/.*: //')
    peak_kb=$(grep "Maximum resident set size" "$log_file" | awk '{print $NF}')

    local secs
    if [[ "$wall" == *:*:* ]]; then
        secs=$(echo "$wall" | awk -F: '{printf "%.1f", $1*3600+$2*60+$3}')
    else
        secs=$(echo "$wall" | awk -F: '{printf "%.1f", $1*60+$2}')
    fi
    local peak_mb
    peak_mb=$(awk "BEGIN{printf \"%.0f\", $peak_kb/1024}")

    echo "${secs} ${peak_mb} ${status}"
    return "$status"
}

gfa_stats() {
    local f=$1
    local segs links paths avg
    segs=$(grep -c '^S' "$f" || true)
    links=$(grep -c '^L' "$f" || true)
    paths=$(grep -c '^P' "$f" || true)
    avg=$(awk -F'\t' '/^S/{sum+=length($3); n++} END{if(n>0) printf "%.0f", sum/n; else print 0}' "$f")
    echo "$segs $links $paths $avg"
}

record() {
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "${4:-}" "${5:-}" "${6:-}" "${7:-}" "$8" >> "$PERF_TSV"
}

record_validation() {
    printf "%s\t%s\t%s\n" "$1" "$2" "${3:-}" >> "$VAL_TSV"
}

# Check that a GFA's path names match the expected sorted sequence list.
# GFA paths are named "seqname:start-end"; strip coords and deduplicate before comparing.
check_paths_vs_input() {
    local label="$1" gfa="$2" expected_sorted="$3"

    if [ ! -s "$gfa" ]; then
        echo "  SKIP: $label input-path-check — GFA missing"
        record_validation "$label input-path-check" "SKIP" "GFA missing"
        return
    fi

    local graph_paths
    if [ -n "$ODGI" ]; then
        graph_paths=$("$ODGI" paths -i "$gfa" -L 2>/dev/null | sed 's/:[0-9]*-[0-9]*$//' | sort -u)
    else
        graph_paths=$(grep '^P' "$gfa" | cut -f2 | sed 's/:[0-9]*-[0-9]*$//' | sort -u)
    fi

    local n_expected n_graph diff_out diff_lines missing extra
    n_expected=$(printf '%s\n' "$expected_sorted" | grep -c . || true)
    n_graph=$(printf '%s\n' "$graph_paths" | grep -c . || true)
    diff_out=$(diff <(printf '%s\n' "$expected_sorted") <(printf '%s\n' "$graph_paths") || true)
    diff_lines=$(printf '%s\n' "$diff_out" | grep -c . || true)
    missing=$(printf '%s\n' "$diff_out" | grep '^<' | sed 's/^< //' | head -5 | tr '\n' ' ' || true)
    extra=$(printf '%s\n'   "$diff_out" | grep '^>' | sed 's/^> //' | head -5 | tr '\n' ' ' || true)

    echo "  $label input-path-check: expected=$n_expected graph=$n_graph diff_lines=$diff_lines"
    [ -n "$missing" ] && echo "    missing from graph: $missing"
    [ -n "$extra"   ] && echo "    extra in graph:     $extra"

    if [ "$diff_lines" -eq 0 ]; then
        record_validation "$label input-path-check" "PASS" "${n_graph} paths match input"
    else
        record_validation "$label input-path-check" "FAIL" \
            "expected=${n_expected} graph=${n_graph} diff_lines=${diff_lines}"
    fi
}

compare_gfa() {
    local gfa_a="$1" gfa_b="$2"

    read -r s_a l_a p_a avg_a <<< "$(gfa_stats "$gfa_a")"
    read -r s_b l_b p_b avg_b <<< "$(gfa_stats "$gfa_b")"

    local detail_a="${s_a}S/${l_a}L/${p_a}P/avg${avg_a}bp"
    local detail_b="${s_b}S/${l_b}L/${p_b}P/avg${avg_b}bp"

    local paths_a paths_b path_diff
    paths_a=$(grep '^P' "$gfa_a" | cut -f2 | sed 's/:[0-9]*-[0-9]*$//' | sort -u)
    paths_b=$(grep '^P' "$gfa_b" | cut -f2 | sed 's/:[0-9]*-[0-9]*$//' | sort -u)
    path_diff=$(diff <(echo "$paths_a") <(echo "$paths_b") | wc -l || true)

    echo "  full:        $detail_a"
    echo "  partitioned: $detail_b"
    echo "  path coverage diff lines: $path_diff"

    if [ "$path_diff" -eq 0 ]; then
        record_validation "full vs partitioned" "PATHS_MATCH" "full=$detail_a partitioned=$detail_b"
    else
        record_validation "full vs partitioned" "PATHS_DIFFER" "full=$detail_a partitioned=$detail_b path_diff=$path_diff"
    fi
}

# -------------------------------------------------------------------------
# Setup: extract BED to get expected sequences
# -------------------------------------------------------------------------

PREFIX="$OUTDIR/run"
BED_FILE="${PREFIX}.bed"

echo "Region: $REGION  partition_size: $PARTITION_SIZE"
echo ""
echo "Extracting BED for expected sequences ..."
$IMPG query --alignment-list data/hprcv2/tpas/tpa-list.txt --sequence-files "$AGC" \
    -r "$REGION" -o bed --force-large-region \
    -O "${PREFIX}" -t "$THREADS" -v 0 2>/dev/null || true

if [ -s "$BED_FILE" ]; then
    EXPECTED_SORTED=$(cut -f1 "$BED_FILE" | sort -u)
    N_EXPECTED=$(printf '%s\n' "$EXPECTED_SORTED" | grep -c . || true)
    echo "  $N_EXPECTED unique sequences in BED"
else
    EXPECTED_SORTED=""
    echo "  WARNING: BED extraction failed"
fi

# -------------------------------------------------------------------------
# Run 1: full pggb (baseline)
# -------------------------------------------------------------------------

FULL_GFA="${PREFIX}.full.pggb.gfa"
echo ""
echo "--- Run 1: full pggb --force-large-region --batch-bytes ${BATCH_BYTES} ---"
metrics=$(run_timed "${PREFIX}.full.pggb.log" $IMPG query \
    --alignment-list data/hprcv2/tpas/tpa-list.txt --sequence-files "$AGC" \
    -r "$REGION" -o gfa --gfa-engine pggb \
    --force-large-region --batch-bytes "$BATCH_BYTES" \
    -O "${PREFIX}.full.pggb" -t "$THREADS" -v 1) || true
read -r wall mem st <<< "$metrics"
if [ "$st" -eq 0 ] && [ -s "$FULL_GFA" ]; then
    read -r s l p avg <<< "$(gfa_stats "$FULL_GFA")"
    echo "  OK  wall=${wall}s  mem=${mem}MB  ${s}S/${l}L/${p}P/avg${avg}bp"
    record "full.pggb" "$wall" "$mem" "$s" "$l" "$p" "$avg" "OK"
    check_paths_vs_input "full.pggb" "$FULL_GFA" "$EXPECTED_SORTED"
else
    echo "  FAIL"
    record "full.pggb" "$wall" "$mem" "" "" "" "" "FAIL"
fi

# -------------------------------------------------------------------------
# Run 2: partitioned pggb
# -------------------------------------------------------------------------

PART_GFA="${PREFIX}.part.pggb.gfa"
echo ""
echo "--- Run 2: partitioned pggb:${PARTITION_SIZE} ---"
metrics=$(run_timed "${PREFIX}.part.pggb.log" $IMPG query \
    --alignment-list data/hprcv2/tpas/tpa-list.txt --sequence-files "$AGC" \
    -r "$REGION" -o gfa --gfa-engine "pggb:${PARTITION_SIZE}" \
    -O "${PREFIX}.part.pggb" -t "$THREADS" -v 1) || true
read -r wall mem st <<< "$metrics"
if [ "$st" -eq 0 ] && [ -s "$PART_GFA" ]; then
    read -r s l p avg <<< "$(gfa_stats "$PART_GFA")"
    echo "  OK  wall=${wall}s  mem=${mem}MB  ${s}S/${l}L/${p}P/avg${avg}bp"
    record "part.pggb:${PARTITION_SIZE}" "$wall" "$mem" "$s" "$l" "$p" "$avg" "OK"
    check_paths_vs_input "part.pggb:${PARTITION_SIZE}" "$PART_GFA" "$EXPECTED_SORTED"
else
    echo "  FAIL"
    record "part.pggb:${PARTITION_SIZE}" "$wall" "$mem" "" "" "" "" "FAIL"
fi

# -------------------------------------------------------------------------
# Compare
# -------------------------------------------------------------------------

echo ""
echo "--- Comparison ---"
if [ -s "$FULL_GFA" ] && [ -s "$PART_GFA" ]; then
    compare_gfa "$FULL_GFA" "$PART_GFA"
else
    echo "  SKIP: one or both GFAs missing"
fi

# -------------------------------------------------------------------------
# Summary
# -------------------------------------------------------------------------

echo ""
echo "================================================================="
echo "PERFORMANCE SUMMARY"
echo "================================================================="
printf "%-26s  %7s  %7s  %6s  %6s  %6s  %7s  %s\n" \
       "Command" "Wall(s)" "Mem(MB)" "S" "L" "P" "AvgSeg" "Status"
printf "%-26s  %7s  %7s  %6s  %6s  %6s  %7s  %s\n" \
       "--------------------------" "-------" "-------" "------" "------" "------" "-------" "------"
tail -n +2 "$PERF_TSV" | while IFS=$'\t' read -r cmd wall mem s l p avg status; do
    printf "%-26s  %7s  %7s  %6s  %6s  %6s  %7s  %s\n" \
        "$cmd" "$wall" "$mem" "$s" "$l" "$p" "$avg" "$status"
done

echo ""
echo "================================================================="
echo "VALIDATION SUMMARY"
echo "================================================================="
tail -n +2 "$VAL_TSV" | while IFS=$'\t' read -r comparison result detail; do
    printf "%-40s  %-14s  %s\n" "$comparison" "$result" "$detail"
done

echo ""
echo "TSV files: $PERF_TSV  $VAL_TSV"
