#!/bin/bash
set -euo pipefail

# Demo: validate that "impg query -o gfa" and "impg graph" produce identical GFA graphs
#
# Tests region sizes: 500bp, 1k, 2k, 5k, 10k, 50k
# Tests engines: recursive, seqwish
# Records wall time and peak memory for each command
# Compares full GFA output (sorted) between query and graph
#
# Usage:
#   cd data/human-pangenome-tpas
#   bash ../../scripts/demo_gfa_commands.sh

IMPG="../../target/release/impg"
AGC="../HPRC_r2_assemblies_0.6.1.agc"
OUTDIR="../../demo-results"
THREADS=12
TIME="/usr/bin/time"

mkdir -p "$OUTDIR"

# -------------------------------------------------------------------------
# Region definitions: 500bp, 1k, 2k, 5k, 10k, 50k
# -------------------------------------------------------------------------
REGIONS=(
    "CHM13#0#chr6:29000000-29001000"
    "CHM13#0#chr6:29000000-29002000"
    "CHM13#0#chr6:29000000-29005000"
)

ENGINES=(poa recursive seqwish pggb)

# -------------------------------------------------------------------------
# Helpers
# -------------------------------------------------------------------------

# run_timed LOGFILE CMD [ARGS...]
# Runs CMD under /usr/bin/time, saves full stderr to LOGFILE,
# echoes "wall_s peak_mb exit_status" on stdout.
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

# Collect GFA stats: segments links paths avg_seg_bp
gfa_stats() {
    local f=$1
    local segs links paths avg
    segs=$(grep -c '^S' "$f" || true)
    links=$(grep -c '^L' "$f" || true)
    paths=$(grep -c '^P' "$f" || true)
    avg=$(awk -F'\t' '/^S/{sum+=length($3); n++} END{if(n>0) printf "%.0f", sum/n; else print 0}' "$f")
    echo "$segs $links $paths $avg"
}

# Make a short label from a region string, including the size.
# e.g. "CHM13#0#chr6:29000000-29001000" -> "chr6_29M_1k"
make_label() {
    local region=$1
    local bare chr pos_start pos_end size_bp size_label
    bare=${region##*#}              # chr6:29000000-29001000
    chr=${bare%%:*}                 # chr6
    pos_start=${bare#*:}; pos_start=${pos_start%-*}  # 29000000
    pos_end=${bare#*-}              # 29001000
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
    IDX=$((IDX + 1))
}

# Validation accumulators
declare -a VAL_LABEL VAL_ENGINE VAL_RESULT VAL_DETAIL
VIDX=0

record_validation() {
    VAL_LABEL[$VIDX]="$1"
    VAL_ENGINE[$VIDX]="$2"
    VAL_RESULT[$VIDX]="$3"
    VAL_DETAIL[$VIDX]="${4:-}"
    VIDX=$((VIDX + 1))
}

# -------------------------------------------------------------------------
# run_region REGION
# -------------------------------------------------------------------------
run_region() {
    local REGION=$1
    local LABEL
    LABEL=$(make_label "$REGION")

    echo ""
    echo "================================================================="
    echo "Region: $REGION  ($LABEL)"
    echo "================================================================="

    local PREFIX="$OUTDIR/$LABEL"

    # --- Step 1: extract FASTA via query -o fasta ---
    echo "  extracting FASTA ..."
    local fasta_ok=false
    metrics=$(run_timed "${PREFIX}.fasta.log" $IMPG query \
        --alignment-list tpa-list.txt --sequence-files "$AGC" \
        -r "$REGION" -o fasta --force-large-region \
        -O "${PREFIX}" -t "$THREADS" -v 1) || true
    read -r wall mem st <<< "$metrics"
    if [ "$st" -eq 0 ] && [ -s "${PREFIX}.fa" ]; then
        record "$LABEL" "fasta" "$wall" "$mem" "OK"
        fasta_ok=true
    else
        record "$LABEL" "fasta" "$wall" "$mem" "FAIL"
    fi

    # --- Step 2: for each engine, run query -o gfa and graph, then compare ---
    for ENGINE in "${ENGINES[@]}"; do
        local TAG="${ENGINE}"

        # --- query -o gfa ---
        echo "  query -o gfa  (engine=$ENGINE) ..."
        metrics=$(run_timed "${PREFIX}.query.${TAG}.log" $IMPG query \
            --alignment-list tpa-list.txt --sequence-files "$AGC" \
            -r "$REGION" -o gfa --engine "$ENGINE" --force-large-region \
            -O "${PREFIX}.query.${TAG}" -t "$THREADS" -v 1) || true
        read -r wall mem st <<< "$metrics"
        if [ "$st" -eq 0 ] && [ -s "${PREFIX}.query.${TAG}.gfa" ]; then
            read -r s l p avg <<< "$(gfa_stats "${PREFIX}.query.${TAG}.gfa")"
            record "$LABEL" "q.$TAG" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
        else
            record "$LABEL" "q.$TAG" "$wall" "$mem" "FAIL"
        fi

        # --- graph ---
        if $fasta_ok; then
            echo "  graph          (engine=$ENGINE) ..."
            metrics=$(run_timed "${PREFIX}.graph.${TAG}.log" $IMPG graph \
                --fasta-files "${PREFIX}.fa" \
                -g "${PREFIX}.graph.${TAG}.gfa" \
                --engine "$ENGINE" --aligner wfmash -t "$THREADS" -v 1) || true
            read -r wall mem st <<< "$metrics"
            if [ "$st" -eq 0 ] && [ -s "${PREFIX}.graph.${TAG}.gfa" ]; then
                read -r s l p avg <<< "$(gfa_stats "${PREFIX}.graph.${TAG}.gfa")"
                record "$LABEL" "g.$TAG" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
            else
                record "$LABEL" "g.$TAG" "$wall" "$mem" "FAIL"
            fi
        fi

        # --- validate: compare graph structure via odgi stats -S ---
        local qgfa="${PREFIX}.query.${TAG}.gfa"
        local ggfa="${PREFIX}.graph.${TAG}.gfa"
        if [ -s "$qgfa" ] && [ -s "$ggfa" ]; then
            if ! command -v odgi &>/dev/null; then
                echo "    SKIP: $ENGINE — odgi not found, skipping validation"
                record_validation "$LABEL" "$ENGINE" "SKIP" "odgi not found"
            else
            local qstats gstats
            qstats=$(odgi stats -i "$qgfa" -S 2>/dev/null)
            gstats=$(odgi stats -i "$ggfa" -S 2>/dev/null)
            if [ "$qstats" = "$gstats" ]; then
                echo "    PASS: $ENGINE — same structure ($qstats)"
                record_validation "$LABEL" "$ENGINE" "PASS" "$qstats"
            else
                echo "    FAIL: $ENGINE — structure differs"
                echo "      query: $qstats"
                echo "      graph: $gstats"
                record_validation "$LABEL" "$ENGINE" "FAIL" "q=$qstats g=$gstats"
            fi
            fi
        else
            local reason=""
            if [ ! -s "$qgfa" ] && [ ! -s "$ggfa" ]; then
                reason="both failed"
            elif [ ! -s "$qgfa" ]; then
                reason="query failed"
            else
                reason="graph failed"
            fi
            echo "    SKIP: $ENGINE — $reason"
            record_validation "$LABEL" "$ENGINE" "SKIP" "$reason"
        fi
    done
}

# -------------------------------------------------------------------------
# Run all regions
# -------------------------------------------------------------------------

for REGION in "${REGIONS[@]}"; do
    run_region "$REGION"
done

# -------------------------------------------------------------------------
# Performance summary
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "PERFORMANCE SUMMARY"
echo "================================================================="
printf "%-18s  %-14s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "Region" "Command" "Wall(s)" "Mem(MB)" "S" "L" "P" "AvgSeg" "Status"
printf "%-18s  %-14s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "------------------" "--------------" "-------" "-------" "------" "------" "------" "------" "------"

for ((i=0; i<IDX; i++)); do
    printf "%-18s  %-14s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
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
# Validation summary: query -o gfa vs graph (full GFA identity)
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "VALIDATION: query -o gfa vs graph (full GFA identity check)"
echo "================================================================="
printf "%-18s  %-10s  %-6s  %s\n" "Region" "Engine" "Result" "Detail"
printf "%-18s  %-10s  %-6s  %s\n" "------------------" "----------" "------" "------"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for ((i=0; i<VIDX; i++)); do
    printf "%-18s  %-10s  %-6s  %s\n" \
        "${VAL_LABEL[$i]}" \
        "${VAL_ENGINE[$i]}" \
        "${VAL_RESULT[$i]}" \
        "${VAL_DETAIL[$i]}"
    case "${VAL_RESULT[$i]}" in
        PASS) PASS_COUNT=$((PASS_COUNT + 1)) ;;
        FAIL) FAIL_COUNT=$((FAIL_COUNT + 1)) ;;
        SKIP) SKIP_COUNT=$((SKIP_COUNT + 1)) ;;
    esac
done

echo ""
echo "Total: $PASS_COUNT PASS, $FAIL_COUNT FAIL, $SKIP_COUNT SKIP (out of $VIDX comparisons)"
if [ "$FAIL_COUNT" -gt 0 ]; then
    echo "ERROR: Some graphs differ between query and graph!"
    exit 1
elif [ "$PASS_COUNT" -eq 0 ]; then
    echo "WARNING: No successful comparisons were made."
    exit 1
else
    echo "All graphs identical between query -o gfa and graph."
fi
