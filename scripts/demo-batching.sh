#!/bin/bash
set -euo pipefail

# Demo: validate that --batch-bytes produces equivalent alignments to unbatched runs.
#
# Tests aligners: wfmash, fastga
# For each aligner, compares unbatched vs batched alignment output (impg align),
# then builds graphs from both PAFs to verify they're healthy.
#
# Batch sizes are per-aligner to trigger multi-batch partitioning given the cost
# models (wfmash: 500MB + 20 bytes/bp memory, fastga: 100MB + 12 bytes/bp disk).
#
# Uses TPA alignment lists and AGC archive as input (same as demo_gfa_commands.sh).
#
# Usage:
#   cd data/human-pangenome-tpas
#   bash ../../scripts/demo-batching.sh

IMPG="../../target/release/impg"
AGC="../HPRC_r2_assemblies_0.6.1.agc"
OUTDIR="../../demo-batching-results"
THREADS=12
TIME="/usr/bin/time"

mkdir -p "$OUTDIR"

# Region to test — 50k on chr6 MHC
REGION="CHM13#0#chr6:29000000-29050000"
ENGINE="seqwish"

# Batch sizes per aligner, chosen to force multiple batches with the extracted data.
# fastga cost model: 100MB overhead + 12 bytes/bp
declare -A BATCH_SIZES
BATCH_SIZES[fastga]="120M"

ALIGNERS=("fastga")

# -------------------------------------------------------------------------
# Helpers (from demo_gfa_commands.sh)
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

# Extract number of batches from a log file
# Looks for "[batch] Partitioned ... into N batches" or "no batching needed"
extract_num_batches() {
    local logfile=$1
    if grep -q "no batching needed" "$logfile" 2>/dev/null; then
        echo "1"
    elif grep -q "Partitioned.*into.*batches" "$logfile" 2>/dev/null; then
        grep "Partitioned.*into.*batches" "$logfile" | grep -oP 'into \K\d+' | head -1
    else
        echo "-"
    fi
}

# Measure peak disk usage of a directory from a monitor log (bytes\tpath format)
peak_disk_mb() {
    local logfile=$1
    if [ -s "$logfile" ]; then
        local peak_bytes
        peak_bytes=$(sort -k1 -n "$logfile" | tail -1 | awk '{print $1}')
        awk "BEGIN{printf \"%.0f\", $peak_bytes/1048576}"
    else
        echo "-"
    fi
}

# -------------------------------------------------------------------------
# Summary accumulators
# -------------------------------------------------------------------------
declare -a SUM_ALIGNER SUM_BATCH SUM_WALL SUM_MEM SUM_DISK SUM_STATUS SUM_LINES SUM_NBATCH
IDX=0

record() {
    SUM_ALIGNER[$IDX]="$1"
    SUM_BATCH[$IDX]="$2"
    SUM_WALL[$IDX]="$3"
    SUM_MEM[$IDX]="$4"
    SUM_DISK[$IDX]="$5"
    SUM_STATUS[$IDX]="$6"
    SUM_LINES[$IDX]="${7:-}"
    SUM_NBATCH[$IDX]="${8:-}"
    IDX=$((IDX + 1))
}

# Validation
declare -a VAL_ALIGNER VAL_PAIR VAL_RESULT VAL_DETAIL
VIDX=0

record_validation() {
    VAL_ALIGNER[$VIDX]="$1"
    VAL_PAIR[$VIDX]="$2"
    VAL_RESULT[$VIDX]="$3"
    VAL_DETAIL[$VIDX]="${4:-}"
    VIDX=$((VIDX + 1))
}

# -------------------------------------------------------------------------
# Main
# -------------------------------------------------------------------------

LABEL=$(make_label "$REGION")
PREFIX="$OUTDIR/$LABEL"

echo "================================================================="
echo "Batching demo: $REGION ($LABEL)"
echo "Engine: $ENGINE   Threads: $THREADS"
echo "================================================================="

# Step 1: extract FASTA
echo ""
echo "--- Extracting FASTA ---"
metrics=$(run_timed "${PREFIX}.fasta.log" $IMPG query \
    --alignment-list tpa-list.txt --sequence-files "$AGC" \
    -r "$REGION" -o fasta --force-large-region \
    -O "${PREFIX}" -t "$THREADS" -v 1) || true
read -r wall mem st <<< "$metrics"
if [ "$st" -ne 0 ] || [ ! -s "${PREFIX}.fa" ]; then
    echo "FATAL: FASTA extraction failed (exit=$st)"
    exit 1
fi
NUM_SEQS=$(grep -c '^>' "${PREFIX}.fa")
TOTAL_BP=$(grep -v '^>' "${PREFIX}.fa" | wc -c)
echo "  FASTA extracted: ${PREFIX}.fa (${NUM_SEQS} seqs, ${TOTAL_BP} bp, ${wall}s, ${mem}MB)"

# Step 2: for each aligner, run unbatched and batched alignments
declare -A PAF_FILES

for ALIGNER in "${ALIGNERS[@]}"; do
    BATCH_SIZE="${BATCH_SIZES[$ALIGNER]}"

    for MODE in "none" "$BATCH_SIZE"; do
        TAG="${ALIGNER}.batch_${MODE}"
        ALIGN_DIR="${PREFIX}.align.${TAG}"
        PAF="${ALIGN_DIR}/alignments.paf"
        LOG="${PREFIX}.align.${TAG}.log"
        DISK_LOG="${PREFIX}.align.${TAG}.disk.log"
        TMPDIR_ALIGN="/tmp/impg_batch_demo_${TAG}_$$"
        PAF_FILES["${ALIGNER}:${MODE}"]="$PAF"

        rm -rf "$ALIGN_DIR" "$TMPDIR_ALIGN"
        mkdir -p "$TMPDIR_ALIGN"

        BATCH_FLAG=()
        if [ "$MODE" != "none" ]; then
            BATCH_FLAG=(--batch-bytes "$MODE")
        fi

        echo ""
        echo "--- impg align: $ALIGNER / batch=$MODE ---"

        # Monitor disk usage of temp dir
        (while true; do du -sb "$TMPDIR_ALIGN" 2>/dev/null; sleep 0.2; done) > "$DISK_LOG" 2>/dev/null &
        MON_PID=$!

        metrics=$(run_timed "$LOG" $IMPG align \
            --sequence-files "${PREFIX}.fa" \
            -o "$ALIGN_DIR" --format paf \
            --aligner "$ALIGNER" "${BATCH_FLAG[@]}" \
            --temp-dir "$TMPDIR_ALIGN" \
            -t "$THREADS" -v 1) || true
        read -r wall mem st <<< "$metrics"

        kill $MON_PID 2>/dev/null || true
        wait $MON_PID 2>/dev/null || true

        nbatch=$(extract_num_batches "$LOG")
        disk=$(peak_disk_mb "$DISK_LOG")

        if [ "$st" -eq 0 ] && [ -s "$PAF" ]; then
            lines=$(wc -l < "$PAF")
            echo "  OK: ${lines} alignments (${wall}s, ${mem}MB RSS, ${disk}MB disk, batches=${nbatch})"
            record "$ALIGNER" "$MODE" "$wall" "$mem" "$disk" "OK" "$lines" "$nbatch"
        else
            echo "  FAIL (exit=$st)"
            record "$ALIGNER" "$MODE" "$wall" "$mem" "$disk" "FAIL" "" "$nbatch"
        fi

        rm -rf "$TMPDIR_ALIGN"
    done
done

# Step 3: validate PAF equivalence (core columns only, ignoring chain/scaffold tags)
echo ""
echo "================================================================="
echo "VALIDATION: PAF alignment equivalence (core columns 1-12)"
echo "================================================================="

for ALIGNER in "${ALIGNERS[@]}"; do
    BATCH_SIZE="${BATCH_SIZES[$ALIGNER]}"
    PAF_NONE="${PAF_FILES["${ALIGNER}:none"]}"
    PAF_BATCHED="${PAF_FILES["${ALIGNER}:${BATCH_SIZE}"]}"

    if [ ! -s "$PAF_NONE" ] || [ ! -s "$PAF_BATCHED" ]; then
        echo "  SKIP: $ALIGNER — missing PAF"
        record_validation "$ALIGNER" "PAF none vs batch=$BATCH_SIZE" "SKIP" "missing PAF"
        continue
    fi

    ndiff=$(diff <(cut -f1-12 "$PAF_NONE" | sort) <(cut -f1-12 "$PAF_BATCHED" | sort) | wc -l)
    if [ "$ndiff" -eq 0 ]; then
        echo "  PASS: $ALIGNER PAF none vs batch=$BATCH_SIZE (0 differences)"
        record_validation "$ALIGNER" "PAF none vs batch=$BATCH_SIZE" "PASS" "0 diff lines"
    else
        echo "  FAIL: $ALIGNER PAF none vs batch=$BATCH_SIZE ($ndiff diff lines)"
        record_validation "$ALIGNER" "PAF none vs batch=$BATCH_SIZE" "FAIL" "$ndiff diff lines"
    fi
done

# Step 4: build graphs from all PAFs and compare
echo ""
echo "================================================================="
echo "VALIDATION: Graph equivalence (build graph from each PAF)"
echo "================================================================="

declare -A GFA_FILES

for ALIGNER in "${ALIGNERS[@]}"; do
    BATCH_SIZE="${BATCH_SIZES[$ALIGNER]}"

    for MODE in "none" "$BATCH_SIZE"; do
        TAG="${ALIGNER}.batch_${MODE}"
        PAF="${PAF_FILES["${ALIGNER}:${MODE}"]}"
        GFA="${PREFIX}.graph.${TAG}.gfa"
        GFA_FILES["${ALIGNER}:${MODE}"]="$GFA"

        if [ ! -s "$PAF" ]; then
            echo "  SKIP: $TAG — no PAF"
            continue
        fi

        $IMPG graph \
            --sequence-files "${PREFIX}.fa" \
            -a "$PAF" -g "$GFA" \
            --gfa-engine "$ENGINE" -t "$THREADS" -v 1 2>/dev/null

        if [ -s "$GFA" ]; then
            read -r s l p avg <<< "$(gfa_stats "$GFA")"
            echo "  $TAG: ${s}S ${l}L ${p}P avg=${avg}bp"
        else
            echo "  $TAG: FAIL (no GFA produced)"
        fi
    done
done

# Compare GFA stats between batched and unbatched for each aligner
echo ""
for ALIGNER in "${ALIGNERS[@]}"; do
    BATCH_SIZE="${BATCH_SIZES[$ALIGNER]}"
    GFA_NONE="${GFA_FILES["${ALIGNER}:none"]}"
    GFA_BATCHED="${GFA_FILES["${ALIGNER}:${BATCH_SIZE}"]}"

    if [ ! -s "$GFA_NONE" ] || [ ! -s "$GFA_BATCHED" ]; then
        echo "  SKIP: $ALIGNER GFA — missing"
        record_validation "$ALIGNER" "GFA none vs batch=$BATCH_SIZE" "SKIP" "missing GFA"
        continue
    fi

    if command -v odgi &>/dev/null; then
        stats_a=$(odgi stats -i "$GFA_NONE" -S 2>/dev/null)
        stats_b=$(odgi stats -i "$GFA_BATCHED" -S 2>/dev/null)
        if [ "$stats_a" = "$stats_b" ]; then
            echo "  PASS: $ALIGNER GFA none vs batch=$BATCH_SIZE ($stats_a)"
            record_validation "$ALIGNER" "GFA none vs batch=$BATCH_SIZE" "PASS" "$stats_a"
        else
            echo "  FAIL: $ALIGNER GFA none vs batch=$BATCH_SIZE"
            echo "    none:    $stats_a"
            echo "    batched: $stats_b"
            record_validation "$ALIGNER" "GFA none vs batch=$BATCH_SIZE" "FAIL" "none=$stats_a batched=$stats_b"
        fi
    else
        # Fallback: compare line counts
        s_a=$(grep -c '^S' "$GFA_NONE")
        s_b=$(grep -c '^S' "$GFA_BATCHED")
        l_a=$(grep -c '^L' "$GFA_NONE")
        l_b=$(grep -c '^L' "$GFA_BATCHED")
        p_a=$(grep -c '^P' "$GFA_NONE")
        p_b=$(grep -c '^P' "$GFA_BATCHED")
        if [ "$s_a" = "$s_b" ] && [ "$l_a" = "$l_b" ] && [ "$p_a" = "$p_b" ]; then
            echo "  PASS: $ALIGNER GFA none vs batch=$BATCH_SIZE (${s_a}S ${l_a}L ${p_a}P)"
            record_validation "$ALIGNER" "GFA none vs batch=$BATCH_SIZE" "PASS" "${s_a}S ${l_a}L ${p_a}P"
        else
            echo "  FAIL: $ALIGNER GFA none vs batch=$BATCH_SIZE"
            echo "    none:    ${s_a}S ${l_a}L ${p_a}P"
            echo "    batched: ${s_b}S ${l_b}L ${p_b}P"
            record_validation "$ALIGNER" "GFA none vs batch=$BATCH_SIZE" "FAIL" "none=${s_a}S/${l_a}L/${p_a}P batched=${s_b}S/${l_b}L/${p_b}P"
        fi
    fi
done

# -------------------------------------------------------------------------
# Performance summary
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "PERFORMANCE SUMMARY"
echo "================================================================="
printf "%-10s  %-10s  %7s  %7s  %7s  %7s  %9s  %s\n" \
       "Aligner" "Batch" "Wall(s)" "RSS(MB)" "Disk(MB)" "Batches" "PAF_lines" "Status"
printf "%-10s  %-10s  %7s  %7s  %7s  %7s  %9s  %s\n" \
       "----------" "----------" "-------" "-------" "-------" "-------" "---------" "------"

for ((i=0; i<IDX; i++)); do
    printf "%-10s  %-10s  %7s  %7s  %7s  %7s  %9s  %s\n" \
        "${SUM_ALIGNER[$i]}" \
        "${SUM_BATCH[$i]}" \
        "${SUM_WALL[$i]}" \
        "${SUM_MEM[$i]}" \
        "${SUM_DISK[$i]:--}" \
        "${SUM_NBATCH[$i]:--}" \
        "${SUM_LINES[$i]:--}" \
        "${SUM_STATUS[$i]}"
done

# -------------------------------------------------------------------------
# Validation summary
# -------------------------------------------------------------------------
echo ""
echo "================================================================="
echo "VALIDATION SUMMARY"
echo "================================================================="

PASS_COUNT=0
FAIL_COUNT=0
SKIP_COUNT=0

for ((i=0; i<VIDX; i++)); do
    printf "%-10s  %-30s  %-6s  %s\n" \
        "${VAL_ALIGNER[$i]}" \
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
    echo "ERROR: Some batched results differ from unbatched!"
    exit 1
elif [ "$PASS_COUNT" -eq 0 ]; then
    echo "WARNING: No successful comparisons were made."
    exit 1
else
    echo "All batched results equivalent to unbatched."
fi
