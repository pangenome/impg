#!/bin/bash
set -euo pipefail

# Demo: validate --partition-size flag across engines, region sizes, and window sizes.
#
# For each (region, engine, partition-size) triple, compares:
#   - Baseline: query -o gfa (no partitioning, uses --force-large-region)
#   - Partitioned: query -o gfa --partition-size <ps>
#
# Also tests the `graph` command with --partition-size.
#
# Regions:
#   - 20k: CHM13#0#chr6:29000000-29020000
#   - 50k: CHM13#0#chr6:29000000-29050000
#
# Partition sizes: 5000, 10000
# Engines: poa (skipped when partition-size > 5000), seqwish, pggb
#
# Usage:
#   cd data/human-pangenome-tpas
#   bash ../../scripts/demo-partitioned.sh

IMPG="../../target/release/impg"
AGC="../HPRC_r2_assemblies_0.6.1.agc"
OUTDIR="../../demo-partitioned-results"
THREADS=12
TIME="/usr/bin/time"

mkdir -p "$OUTDIR"

PERF_TSV="$OUTDIR/performance.tsv"
VAL_TSV="$OUTDIR/validation.tsv"

printf "Region\tCommand\tWall_s\tMem_MB\tSegments\tLinks\tPaths\tAvgSegBp\tStatus\n" > "$PERF_TSV"
printf "Region\tEngine\tComparison\tResult\tDetail\n" > "$VAL_TSV"

# -------------------------------------------------------------------------
# Region definitions
# -------------------------------------------------------------------------
REGIONS=(
    "CHM13#0#chr6:29000000-29020000"
    "CHM13#0#chr6:29000000-29050000"
)

PARTITION_SIZES=("5000" "10000")

ENGINES=("poa" "seqwish" "pggb")

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

make_label() {
    local region=$1
    local bare chr pos_start pos_end size_bp size_label
    bare=${region##*#}
    chr=${bare%%:*}
    pos_start=${bare#*:}; pos_start=${pos_start%-*}
    pos_end=${bare#*-}
    size_bp=$((pos_end - pos_start))

    if (( size_bp >= 1000000 )); then
        size_label="$((size_bp / 1000000))M"
    elif (( size_bp >= 1000 )); then
        size_label="$((size_bp / 1000))k"
    else
        size_label="${size_bp}bp"
    fi

    local pos_label
    if (( pos_start >= 1000000000 )); then
        pos_label="$((pos_start / 1000000000))G"
    elif (( pos_start >= 1000000 )); then
        pos_label="$((pos_start / 1000000))M"
    elif (( pos_start >= 1000 )); then
        pos_label="$((pos_start / 1000))k"
    else
        pos_label="$pos_start"
    fi

    echo "${chr}_${pos_label}_${size_label}"
}

# -------------------------------------------------------------------------
# Summary accumulators
# -------------------------------------------------------------------------
declare -a SUM_LABEL SUM_CMD SUM_WALL SUM_MEM SUM_STATUS SUM_S SUM_L SUM_P SUM_AVG
IDX=0

record() {
    SUM_LABEL[$IDX]="$1"
    SUM_CMD[$IDX]="$2"
    SUM_WALL[$IDX]="$3"
    SUM_MEM[$IDX]="$4"
    SUM_STATUS[$IDX]="$5"
    SUM_S[$IDX]="${6:-}"
    SUM_L[$IDX]="${7:-}"
    SUM_P[$IDX]="${8:-}"
    SUM_AVG[$IDX]="${9:-}"
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" "${6:-}" "${7:-}" "${8:-}" "${9:-}" "$5" >> "$PERF_TSV"
    IDX=$((IDX + 1))
}

declare -a VAL_LABEL VAL_ENGINE VAL_PAIR VAL_RESULT VAL_DETAIL
VIDX=0

record_validation() {
    VAL_LABEL[$VIDX]="$1"
    VAL_ENGINE[$VIDX]="$2"
    VAL_PAIR[$VIDX]="$3"
    VAL_RESULT[$VIDX]="$4"
    VAL_DETAIL[$VIDX]="${5:-}"
    printf "%s\t%s\t%s\t%s\t%s\n" "$1" "$2" "$3" "$4" "${5:-}" >> "$VAL_TSV"
    VIDX=$((VIDX + 1))
}

compare_gfa() {
    local label=$1 engine=$2 pair_name=$3 gfa_a=$4 gfa_b=$5

    if [ ! -s "$gfa_a" ] || [ ! -s "$gfa_b" ]; then
        local reason=""
        if [ ! -s "$gfa_a" ] && [ ! -s "$gfa_b" ]; then
            reason="both missing"
        elif [ ! -s "$gfa_a" ]; then
            reason="${pair_name%% vs *} missing"
        else
            reason="${pair_name##* vs } missing"
        fi
        echo "    SKIP: $engine $pair_name — $reason"
        record_validation "$label" "$engine" "$pair_name" "SKIP" "$reason"
        return
    fi

    # Compare GFA stats (segments, links, paths)
    read -r s_a l_a p_a avg_a <<< "$(gfa_stats "$gfa_a")"
    read -r s_b l_b p_b avg_b <<< "$(gfa_stats "$gfa_b")"

    local detail_a="${s_a}S/${l_a}L/${p_a}P/avg${avg_a}"
    local detail_b="${s_b}S/${l_b}L/${p_b}P/avg${avg_b}"

    # Check path coverage (same paths represented)
    local paths_a paths_b
    paths_a=$(grep '^P' "$gfa_a" | cut -f2 | sort)
    paths_b=$(grep '^P' "$gfa_b" | cut -f2 | sort)

    local path_diff
    path_diff=$(diff <(echo "$paths_a") <(echo "$paths_b") | wc -l || true)

    echo "    $engine $pair_name:"
    echo "      baseline:    $detail_a"
    echo "      partitioned: $detail_b"
    echo "      path coverage diff lines: $path_diff"

    if [ "$p_a" = "$p_b" ] && [ "$path_diff" -eq 0 ]; then
        record_validation "$label" "$engine" "$pair_name" "PATHS_MATCH" "baseline=$detail_a partitioned=$detail_b"
    else
        record_validation "$label" "$engine" "$pair_name" "PATHS_DIFFER" "baseline=$detail_a partitioned=$detail_b path_diff=$path_diff"
    fi
}

# -------------------------------------------------------------------------
# Part 1: query command tests
# -------------------------------------------------------------------------

echo "================================================================="
echo "Part 1: query -o gfa — baseline vs --partition-size"
echo "================================================================="

for i in "${!REGIONS[@]}"; do
    REGION="${REGIONS[$i]}"
    LABEL=$(make_label "$REGION")

    bare=${REGION##*#}
    pos_start=${bare#*:}; pos_start=${pos_start%-*}
    pos_end=${bare#*-}
    size_bp=$((pos_end - pos_start))

    echo ""
    echo "-----------------------------------------------------------------"
    echo "Region: $REGION ($LABEL, ${size_bp}bp)"
    echo "-----------------------------------------------------------------"

    PREFIX="$OUTDIR/$LABEL"

    # Run baseline once per (region, engine) — skip POA (regions too large)
    for ENGINE in "${ENGINES[@]}"; do
        if [ "$ENGINE" = "poa" ]; then
            echo ""
            echo "  [SKIP] poa baseline — region too large for POA"
            continue
        fi

        BASELINE_GFA="${PREFIX}.baseline.${ENGINE}.gfa"
        echo ""
        echo "  --- engine=$ENGINE (baseline) ---"
        echo "  [baseline] query -o gfa --gfa-engine $ENGINE --force-large-region"
        metrics=$(run_timed "${PREFIX}.baseline.${ENGINE}.log" $IMPG query \
            --alignment-list tpa-list.txt --sequence-files "$AGC" \
            -r "$REGION" -o gfa --gfa-engine "$ENGINE" \
            --force-large-region \
            -O "${PREFIX}.baseline.${ENGINE}" -t "$THREADS" -v 1) || true
        read -r wall mem st <<< "$metrics"
        if [ "$st" -eq 0 ] && [ -s "$BASELINE_GFA" ]; then
            read -r s l p avg <<< "$(gfa_stats "$BASELINE_GFA")"
            record "$LABEL" "baseline.$ENGINE" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
        else
            record "$LABEL" "baseline.$ENGINE" "$wall" "$mem" "FAIL"
        fi
    done

    # Run partitioned for each (partition-size, engine), compare to baseline
    for PS in "${PARTITION_SIZES[@]}"; do
        echo ""
        echo "  ========= partition-size=$PS ========="

        for ENGINE in "${ENGINES[@]}"; do
            # Skip POA when partition-size > 5000
            if [ "$ENGINE" = "poa" ] && [ "$PS" -gt 5000 ]; then
                echo "  [SKIP] poa engine with partition-size=$PS (> 5000bp)"
                continue
            fi

            BASELINE_GFA="${PREFIX}.baseline.${ENGINE}.gfa"
            PART_GFA="${PREFIX}.ps${PS}.${ENGINE}.gfa"
            echo ""
            echo "  --- engine=$ENGINE, ps=$PS ---"
            echo "  [partitioned] query -o gfa --gfa-engine $ENGINE --partition-size $PS"
            metrics=$(run_timed "${PREFIX}.ps${PS}.${ENGINE}.log" $IMPG query \
                --alignment-list tpa-list.txt --sequence-files "$AGC" \
                -r "$REGION" -o gfa --gfa-engine "$ENGINE" \
                --partition-size "$PS" \
                -O "${PREFIX}.ps${PS}.${ENGINE}" -t "$THREADS" -v 1) || true
            read -r wall mem st <<< "$metrics"
            if [ "$st" -eq 0 ] && [ -s "$PART_GFA" ]; then
                read -r s l p avg <<< "$(gfa_stats "$PART_GFA")"
                record "$LABEL" "ps${PS}.$ENGINE" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
            else
                record "$LABEL" "ps${PS}.$ENGINE" "$wall" "$mem" "FAIL"
            fi

            # Compare to baseline
            compare_gfa "$LABEL" "$ENGINE" "baseline vs ps${PS}" "$BASELINE_GFA" "$PART_GFA"
        done
    done
done

# -------------------------------------------------------------------------
# Part 2: graph command test
# -------------------------------------------------------------------------

echo ""
echo "================================================================="
echo "Part 2: graph command — baseline vs --partition-size"
echo "================================================================="

# Use the 20k region for graph tests
GRAPH_REGION="${REGIONS[0]}"
GRAPH_LABEL=$(make_label "$GRAPH_REGION")
GRAPH_PREFIX="$OUTDIR/$GRAPH_LABEL"

echo "  Extracting FASTA for graph tests ..."
metrics=$(run_timed "${GRAPH_PREFIX}.fasta.log" $IMPG query \
    --alignment-list tpa-list.txt --sequence-files "$AGC" \
    -r "$GRAPH_REGION" -o fasta --force-large-region \
    -O "${GRAPH_PREFIX}" -t "$THREADS" -v 1) || true
read -r wall mem st <<< "$metrics"

if [ "$st" -eq 0 ] && [ -s "${GRAPH_PREFIX}.fa" ]; then
    record "$GRAPH_LABEL" "fasta" "$wall" "$mem" "OK"

    for ENGINE in "seqwish" "pggb"; do
        # Baseline graph (once per engine)
        BASELINE_GFA="${GRAPH_PREFIX}.graph.baseline.${ENGINE}.gfa"
        echo ""
        echo "  --- graph engine=$ENGINE (baseline) ---"
        echo "  [baseline] graph --gfa-engine $ENGINE"
        metrics=$(run_timed "${GRAPH_PREFIX}.graph.baseline.${ENGINE}.log" $IMPG graph \
            --sequence-files "${GRAPH_PREFIX}.fa" \
            -g "$BASELINE_GFA" \
            --gfa-engine "$ENGINE" -t "$THREADS" -v 1) || true
        read -r wall mem st <<< "$metrics"
        if [ "$st" -eq 0 ] && [ -s "$BASELINE_GFA" ]; then
            read -r s l p avg <<< "$(gfa_stats "$BASELINE_GFA")"
            record "$GRAPH_LABEL" "graph.baseline.$ENGINE" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
        else
            record "$GRAPH_LABEL" "graph.baseline.$ENGINE" "$wall" "$mem" "FAIL"
        fi

        # Partitioned graph for each partition size
        for PS in "${PARTITION_SIZES[@]}"; do
            PART_GFA="${GRAPH_PREFIX}.graph.ps${PS}.${ENGINE}.gfa"
            echo "  [partitioned] graph --gfa-engine $ENGINE --partition-size $PS"
            metrics=$(run_timed "${GRAPH_PREFIX}.graph.ps${PS}.${ENGINE}.log" $IMPG graph \
                --sequence-files "${GRAPH_PREFIX}.fa" \
                -g "$PART_GFA" \
                --gfa-engine "$ENGINE" --partition-size "$PS" -t "$THREADS" -v 1) || true
            read -r wall mem st <<< "$metrics"
            if [ "$st" -eq 0 ] && [ -s "$PART_GFA" ]; then
                read -r s l p avg <<< "$(gfa_stats "$PART_GFA")"
                record "$GRAPH_LABEL" "graph.ps${PS}.$ENGINE" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
            else
                record "$GRAPH_LABEL" "graph.ps${PS}.$ENGINE" "$wall" "$mem" "FAIL"
            fi

            compare_gfa "$GRAPH_LABEL" "$ENGINE" "graph baseline vs ps${PS}" "$BASELINE_GFA" "$PART_GFA"
        done
    done
else
    echo "  FASTA extraction failed — skipping graph tests"
    record "$GRAPH_LABEL" "fasta" "$wall" "$mem" "FAIL"
fi

# -------------------------------------------------------------------------
# Part 3: Error case validation
# -------------------------------------------------------------------------

echo ""
echo "================================================================="
echo "Part 3: Error case validation"
echo "================================================================="

# partition-size too small
echo -n "  partition-size < 1000: "
if $IMPG query --alignment-list tpa-list.txt --sequence-files "$AGC" \
    -r "${REGIONS[0]}" -o gfa --partition-size 500 -t "$THREADS" -v 0 2>&1 | grep -q "must be at least 1000"; then
    echo "PASS (correctly rejected)"
else
    echo "FAIL (should have been rejected)"
fi

# partition-size + separate-files
echo -n "  partition-size + separate-files: "
if $IMPG partition --alignment-list tpa-list.txt --sequence-files "$AGC" \
    -w 10000 --partition-size 10000 --separate-files -o gfa -t "$THREADS" -v 0 2>&1 | grep -q "mutually exclusive"; then
    echo "PASS (correctly rejected)"
else
    echo "FAIL (should have been rejected)"
fi

# similarity + partition-size
echo -n "  similarity + partition-size: "
if $IMPG similarity --alignment-list tpa-list.txt --sequence-files "$AGC" \
    -r "${REGIONS[0]}" --partition-size 10000 -t "$THREADS" -v 0 2>&1 | grep -q "not yet supported"; then
    echo "PASS (correctly rejected)"
else
    echo "FAIL (should have been rejected)"
fi

# -------------------------------------------------------------------------
# Performance summary
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "PERFORMANCE SUMMARY"
echo "================================================================="
printf "%-18s  %-30s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "Region" "Command" "Wall(s)" "Mem(MB)" "S" "L" "P" "AvgSeg" "Status"
printf "%-18s  %-30s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "------------------" "------------------------------" "-------" "-------" "------" "------" "------" "------" "------"

for ((i=0; i<IDX; i++)); do
    printf "%-18s  %-30s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
        "${SUM_LABEL[$i]}" \
        "${SUM_CMD[$i]}" \
        "${SUM_WALL[$i]}" \
        "${SUM_MEM[$i]}" \
        "${SUM_S[$i]:--}" \
        "${SUM_L[$i]:--}" \
        "${SUM_P[$i]:--}" \
        "${SUM_AVG[$i]:--}" \
        "${SUM_STATUS[$i]}"
done

# -------------------------------------------------------------------------
# Validation summary
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "VALIDATION SUMMARY"
echo "================================================================="
printf "%-18s  %-10s  %-30s  %-14s  %s\n" \
       "Region" "Engine" "Comparison" "Result" "Detail"
printf "%-18s  %-10s  %-30s  %-14s  %s\n" \
       "------------------" "----------" "------------------------------" "--------------" "------"

for ((i=0; i<VIDX; i++)); do
    printf "%-18s  %-10s  %-30s  %-14s  %s\n" \
        "${VAL_LABEL[$i]}" \
        "${VAL_ENGINE[$i]}" \
        "${VAL_PAIR[$i]}" \
        "${VAL_RESULT[$i]}" \
        "${VAL_DETAIL[$i]}"
done

echo ""
echo "TSV files saved:"
echo "  Performance: $PERF_TSV"
echo "  Validation:  $VAL_TSV"
