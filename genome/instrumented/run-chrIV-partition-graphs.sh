#!/usr/bin/env bash
# The chrIV partition-embedded graph build + the alignment/walk dumps
# (chrIV end-to-end slice 1). Same recipe as the chrMT/chrI build:
# partition_graph_export in RAW mode, frequency mask disabled, AGC gap
# splicing, global-syng node-id continuity; serial; receipt markers,
# external RSS poller and the 64GiB guard. Receipts beside the
# chrMT/chrI set in the validation dir's partition-graphs/.
#
# Prerequisite: partition-expressibility-census-chrIV.py requests
# (writes the build list + both request files).
set -uo pipefail
PREFIX=/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng
AGC=/home/erikg/yeast/yeast235.agc
ROOT=/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results
OUT=/home/erikg/yeast/genome-balanced-diploid-validation-20260930/partition-graphs
BIN=/home/erikg/impg-genome-inference/target/release/examples/partition_graph_export
DUMP=/home/erikg/impg-genome-inference/target/release/examples/syng_homology_dump
GUARD_KB=$((64*1024*1024))  # the 64GiB guard

# run_phase <name> <stdout-file> <command...>: markers, poller, guard.
# The command's stderr (the progress log) lands in <name>.log; stdout
# lands in <stdout-file> (the export tool writes files itself, so its
# stdout file is just the log tail; the dump tools write receipts to
# stdout).
run_phase() {
  local name="$1"; local outfile="$2"; shift 2
  local log="$OUT/$name.log" rss="$OUT/$name.rss"
  rm -f "$OUT/$name.done" "$OUT/$name.exit" "$OUT/$name.wall"
  date -u +%FT%TZ > "$OUT/$name.start"
  local start; start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$@" > "$outfile" 2> "$log" &
  local pid; pid=$!
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$rss"
      if [ "${rss_kb:-0}" -gt "$GUARD_KB" ]; then
        echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB > ${GUARD_KB}kB" >> "$rss"
        kill -9 "$pid" 2>/dev/null
        break
      fi
      sleep 5
    done ) & local poller; poller=$!
  local status; status=0; wait "$pid" || status=$?
  kill "$poller" 2>/dev/null; wait "$poller" 2>/dev/null
  echo "$status" > "$OUT/$name.exit"
  echo "$(( $(date +%s) - start ))" > "$OUT/$name.wall"
  if [ "$status" -eq 0 ]; then touch "$OUT/$name.done"; fi
  return "$status"
}

mapfile -t IDS < "$OUT/partition-graph-chrIV-build-list.txt"
echo "build set: ${#IDS[@]} partitions"

run_phase build-chrIV-partition-graphs "$OUT/build-chrIV-partition-graphs.stdout" \
  "$BIN" export "$PREFIX" "$AGC" "$ROOT" "$OUT" "${IDS[@]}" || exit 1

run_phase dumps-chrIV-homology "$OUT/partition-graph-chrIV-homology.jsonl" \
  "$DUMP" homology "$PREFIX" \
  "$OUT/partition-graph-chrIV-homology-requests.tsv" \
  || exit 1

run_phase dumps-chrIV-walks "$OUT/partition-graph-chrIV-walks.jsonl" \
  "$DUMP" walks "$PREFIX" \
  "$OUT/partition-graph-chrIV-walk-requests.tsv" \
  || exit 1

date -u +%FT%TZ > "$OUT/dumps-chrIV.done"
echo "chrIV graphs + dumps complete"
