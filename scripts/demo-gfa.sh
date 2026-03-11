#!/bin/bash
set -euo pipefail

# Demo: validate that "impg query -o gfa", "impg graph", and "impg align + impg graph --input-paf"
# produce equivalent GFA graphs across all sparsification strategies.
#
# Tests region sizes: 1k, 2k, 10k
# Tests engines: seqwish, pggb (all sizes) + poa (1k, 2k only)
# Tests sparsification: none, random:0.5, giant:0.95, tree:3:1:0.1, wfmash:auto
# Three command paths per (region, engine, strategy):
#   1. query -o gfa --sparsify
#   2. graph --sparsify
#   3. align --sparsify → PAF, then graph --input-paf
# Compares graph structure via odgi stats -S across all three paths
# Saves performance TSV and validation TSV to OUTDIR for later analysis.
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

# Output TSV files for performance and validation data
PERF_TSV="$OUTDIR/performance.tsv"
VAL_TSV="$OUTDIR/validation.tsv"

printf "Region\tCommand\tWall_s\tMem_MB\tSegments\tLinks\tPaths\tAvgSegBp\tStatus\n" > "$PERF_TSV"
printf "Region\tEngine\tStrategy\tComparison\tResult\tDetail\n" > "$VAL_TSV"

# -------------------------------------------------------------------------
# Region definitions: 1k, 2k, 10k
# -------------------------------------------------------------------------
REGIONS=(
    "CHM13#0#chr6:29000000-29001000"
    "CHM13#0#chr6:29000000-29002000"
    "CHM13#0#chr6:29000000-29010000"
)

# Engines for a given region size (bp).
# poa is only viable for small regions (≤2k); seqwish+pggb for all sizes.
engines_for_size() {
    local size_bp=$1
    if (( size_bp <= 2000 )); then
        echo "seqwish pggb poa"
    else
        echo "seqwish pggb"
    fi
}

STRATEGIES=(
    "none"
    "random:0.5"
    "giant:0.95"
    "tree:3:1:0.1"
    "wfmash:auto"
)

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

# Convert strategy string to filesystem-safe tag: replace : with _
# e.g. "random:0.5" -> "random_0.5", "wfmash:auto" -> "wfmash_auto"
strategy_tag() {
    echo "${1//:/_}"
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
    # Append to TSV
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" "${6:-}" "${7:-}" "${8:-}" "${9:-}" "$5" >> "$PERF_TSV"
    IDX=$((IDX + 1))
}

# Validation accumulators
declare -a VAL_LABEL VAL_ENGINE VAL_STRATEGY VAL_PAIR VAL_RESULT VAL_DETAIL
VIDX=0

record_validation() {
    VAL_LABEL[$VIDX]="$1"
    VAL_ENGINE[$VIDX]="$2"
    VAL_STRATEGY[$VIDX]="$3"
    VAL_PAIR[$VIDX]="$4"
    VAL_RESULT[$VIDX]="$5"
    VAL_DETAIL[$VIDX]="${6:-}"
    # Append to TSV
    printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$1" "$2" "$3" "$4" "$5" "${6:-}" >> "$VAL_TSV"
    VIDX=$((VIDX + 1))
}

# Compare two GFA files via odgi stats -S, record result
compare_gfa() {
    local label=$1 engine=$2 stag=$3 pair_name=$4 gfa_a=$5 gfa_b=$6

    if [ ! -s "$gfa_a" ] || [ ! -s "$gfa_b" ]; then
        local reason=""
        if [ ! -s "$gfa_a" ] && [ ! -s "$gfa_b" ]; then
            reason="both missing"
        elif [ ! -s "$gfa_a" ]; then
            reason="${pair_name%% vs *} missing"
        else
            reason="${pair_name##* vs } missing"
        fi
        echo "    SKIP: $engine/$stag $pair_name — $reason"
        record_validation "$label" "$engine" "$stag" "$pair_name" "SKIP" "$reason"
        return
    fi

    if ! command -v odgi &>/dev/null; then
        echo "    SKIP: $engine/$stag $pair_name — odgi not found"
        record_validation "$label" "$engine" "$stag" "$pair_name" "SKIP" "odgi not found"
        return
    fi

    local stats_a stats_b
    stats_a=$(odgi stats -i "$gfa_a" -S 2>/dev/null)
    stats_b=$(odgi stats -i "$gfa_b" -S 2>/dev/null)
    if [ "$stats_a" = "$stats_b" ]; then
        echo "    PASS: $engine/$stag $pair_name ($stats_a)"
        record_validation "$label" "$engine" "$stag" "$pair_name" "PASS" "$stats_a"
    else
        echo "    FAIL: $engine/$stag $pair_name"
        echo "      ${pair_name%% vs *}: $stats_a"
        echo "      ${pair_name##* vs }: $stats_b"
        record_validation "$label" "$engine" "$stag" "$pair_name" "FAIL" "${pair_name%% vs *}=$stats_a ${pair_name##* vs }=$stats_b"
    fi
}

# -------------------------------------------------------------------------
# run_region REGION
# -------------------------------------------------------------------------
run_region() {
    local REGION=$1
    local LABEL
    LABEL=$(make_label "$REGION")

    # Compute region size and select engines
    local bare=${REGION##*#}
    local pos_start=${bare#*:}; pos_start=${pos_start%-*}
    local pos_end=${bare#*-}
    local size_bp=$((pos_end - pos_start))
    local ENGINES_STR
    ENGINES_STR=$(engines_for_size "$size_bp")
    read -ra REGION_ENGINES <<< "$ENGINES_STR"

    echo ""
    echo "================================================================="
    echo "Region: $REGION  ($LABEL)  engines: ${REGION_ENGINES[*]}"
    echo "================================================================="

    local PREFIX="$OUTDIR/$LABEL"

    # --- Step 1: extract FASTA via query -o fasta (once per region) ---
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

    # --- Step 2: for each engine × strategy, run all three paths ---
    # poa engine does not run alignment, so --sparsify is invalid; only test "none".
    for ENGINE in "${REGION_ENGINES[@]}"; do
        local ENGINE_STRATEGIES=("${STRATEGIES[@]}")
        if [ "$ENGINE" = "poa" ]; then
            ENGINE_STRATEGIES=("none")
        fi
        for STRATEGY in "${ENGINE_STRATEGIES[@]}"; do
            local STAG
            STAG=$(strategy_tag "$STRATEGY")
            local TAG="${ENGINE}.${STAG}"

            echo ""
            echo "  --- engine=$ENGINE  sparsify=$STRATEGY ($TAG) ---"

            # Build sparsify flag (omit for "none")
            local SPARSIFY_FLAG=()
            if [ "$STRATEGY" != "none" ]; then
                SPARSIFY_FLAG=(--sparsify "$STRATEGY")
            fi

            # --- Path 1: query -o gfa ---
            echo "  [1/3] query -o gfa ..."
            local qgfa="${PREFIX}.query.${TAG}.gfa"
            metrics=$(run_timed "${PREFIX}.query.${TAG}.log" $IMPG query \
                --alignment-list tpa-list.txt --sequence-files "$AGC" \
                -r "$REGION" -o gfa --gfa-engine "$ENGINE" "${SPARSIFY_FLAG[@]}" \
                --force-large-region \
                -O "${PREFIX}.query.${TAG}" -t "$THREADS" -v 1) || true
            read -r wall mem st <<< "$metrics"
            if [ "$st" -eq 0 ] && [ -s "$qgfa" ]; then
                read -r s l p avg <<< "$(gfa_stats "$qgfa")"
                record "$LABEL" "q.$TAG" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
            else
                record "$LABEL" "q.$TAG" "$wall" "$mem" "FAIL"
            fi

            # --- Path 2: graph (aligns internally) ---
            local ggfa="${PREFIX}.graph.${TAG}.gfa"
            if $fasta_ok; then
                echo "  [2/3] graph ..."
                metrics=$(run_timed "${PREFIX}.graph.${TAG}.log" $IMPG graph \
                    --sequence-files "${PREFIX}.fa" \
                    -g "$ggfa" \
                    --gfa-engine "$ENGINE" --aligner wfmash "${SPARSIFY_FLAG[@]}" \
                    -t "$THREADS" -v 1) || true
                read -r wall mem st <<< "$metrics"
                if [ "$st" -eq 0 ] && [ -s "$ggfa" ]; then
                    read -r s l p avg <<< "$(gfa_stats "$ggfa")"
                    record "$LABEL" "g.$TAG" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
                else
                    record "$LABEL" "g.$TAG" "$wall" "$mem" "FAIL"
                fi
            fi

            # --- Path 3: align → PAF, then graph --input-paf ---
            local align_dir="${PREFIX}.align.${TAG}"
            local align_paf="${align_dir}/alignments.paf"
            local agfa="${PREFIX}.align_graph.${TAG}.gfa"
            if $fasta_ok; then
                # Step A: align with sparsification → PAF
                echo "  [3/3] align → graph ..."
                metrics=$(run_timed "${PREFIX}.align.${TAG}.log" $IMPG align \
                    --sequence-files "${PREFIX}.fa" \
                    -o "$align_dir" --format paf \
                    --aligner wfmash "${SPARSIFY_FLAG[@]}" \
                    -t "$THREADS" -v 1) || true
                read -r wall mem st <<< "$metrics"
                if [ "$st" -eq 0 ] && [ -s "$align_paf" ]; then
                    record "$LABEL" "a.$TAG" "$wall" "$mem" "OK"

                    # Step B: build graph from pre-computed PAF (no sparsification)
                    metrics=$(run_timed "${PREFIX}.align_graph.${TAG}.log" $IMPG graph \
                        --sequence-files "${PREFIX}.fa" \
                        -a "$align_paf" \
                        -g "$agfa" \
                        --gfa-engine "$ENGINE" -t "$THREADS" -v 1) || true
                    read -r wall mem st <<< "$metrics"
                    if [ "$st" -eq 0 ] && [ -s "$agfa" ]; then
                        read -r s l p avg <<< "$(gfa_stats "$agfa")"
                        record "$LABEL" "ag.$TAG" "$wall" "$mem" "OK" "$s" "$l" "$p" "$avg"
                    else
                        record "$LABEL" "ag.$TAG" "$wall" "$mem" "FAIL"
                    fi
                else
                    record "$LABEL" "a.$TAG" "$wall" "$mem" "FAIL"
                fi
            fi

            # --- Validate: compare all three paths pairwise ---
            compare_gfa "$LABEL" "$ENGINE" "$STAG" "query vs graph" "$qgfa" "$ggfa"
            compare_gfa "$LABEL" "$ENGINE" "$STAG" "query vs align+graph" "$qgfa" "$agfa"
            compare_gfa "$LABEL" "$ENGINE" "$STAG" "graph vs align+graph" "$ggfa" "$agfa"
        done
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
printf "%-18s  %-28s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "Region" "Command" "Wall(s)" "Mem(MB)" "S" "L" "P" "AvgSeg" "Status"
printf "%-18s  %-28s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
       "------------------" "----------------------------" "-------" "-------" "------" "------" "------" "------" "------"

for ((i=0; i<IDX; i++)); do
    printf "%-18s  %-28s  %7s  %7s  %6s  %6s  %6s  %6s  %s\n" \
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
echo "VALIDATION: cross-path graph structure comparison (odgi stats -S)"
echo "================================================================="
printf "%-18s  %-10s  %-14s  %-22s  %-6s  %s\n" \
       "Region" "Engine" "Strategy" "Comparison" "Result" "Detail"
printf "%-18s  %-10s  %-14s  %-22s  %-6s  %s\n" \
       "------------------" "----------" "--------------" "----------------------" "------" "------"

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for ((i=0; i<VIDX; i++)); do
    printf "%-18s  %-10s  %-14s  %-22s  %-6s  %s\n" \
        "${VAL_LABEL[$i]}" \
        "${VAL_ENGINE[$i]}" \
        "${VAL_STRATEGY[$i]}" \
        "${VAL_PAIR[$i]}" \
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
    echo "ERROR: Some graphs differ across command paths!"
    exit 1
elif [ "$PASS_COUNT" -eq 0 ]; then
    echo "WARNING: No successful comparisons were made."
    exit 1
else
    echo "All graphs equivalent across query / graph / align+graph paths."
fi

echo ""
echo "TSV files saved:"
echo "  Performance: $PERF_TSV"
echo "  Validation:  $VAL_TSV"
