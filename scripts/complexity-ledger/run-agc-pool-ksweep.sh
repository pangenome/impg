#!/usr/bin/env bash
# THE AGC DECOMPRESSOR POOL K-SWEEP - the memory/wall trade, measured.
#
# THE BUG THIS PRICES. `AgcIndex` used to open one AGC `Decompressor`
# per ambient thread slot; a decompressor holds the archive's own
# decompression state in RAM, ~the archive's size (3.70 GB/slot against
# the 3.30 GB HPRC panel). The routes derivation therefore peaked at
# 711 GB on a 192-thread run. The pool is now bounded by
# `IMPG_AGC_DECOMPRESSOR_SLOTS` (default 8).
#
# THIS SWEEP answers the only open question: what does bounding cost in
# WALL on the one step that ran at 192 threads? Fix RAYON_NUM_THREADS
# (default 192, the committed form) and vary K alone. Two measurements
# per point, both of record:
#   (1) wall + core-s + peak RSS (the trade), and
#   (2) THE IDENTITY GATE - graph.json byte-identical across every K,
#       because sharing a decompressor must change no byte (the call
#       consumes graph.json's lane inventory; the per-lane files are
#       sealed artifacts checked by the same digest).
#
# `IMPG_AGC_DECOMPRESSOR_SLOTS=<threads>` restores the committed
# behaviour exactly, so the K=192 point IS the before receipt - no code
# revert, no second binary.
#
# Marker-gated and resumable per K. The workdir self-drains (< budget);
# each point's out-dir is removed after its receipt is mined.
set -uo pipefail
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh

SWEEP=${SWEEP:-/mnt/nvme2n1/erikg/agc-pool-ksweep}
RCP=$SWEEP/ksweep-receipt.txt
CENSUS=${CENSUS:-$EE/cut-census-L1/cut-census.jsonl}
# SMALL FIRST: K=8 is the shipped default, and the committed behaviour
# (one decompressor per thread) is the LAST point, so the safe receipts
# land before the memory-hungry one. K=1 is the serial floor - the
# committed serial form measured ~7.5h - so it is NOT in the default
# list; run it explicitly (KS="1") if the floor is ever needed.
KS=${KS:-"8 16 32 64 192"}
THREADS=${RAYON_NUM_THREADS:-192}
mkdir -p "$SWEEP"

[ -f "$CENSUS" ] || { echo "no census at $CENSUS (the routing-bind input)" >&2; exit 1; }
[ -x "$IMPG_NEW" ] || { echo "no binary at $IMPG_NEW" >&2; exit 1; }
[ -f "$RCP" ] || {
  : > "$RCP"
  echo "# THE AGC DECOMPRESSOR POOL K-SWEEP (R=$(date -u +%FT%TZ))" >> "$RCP"
  echo "# binary=$IMPG_NEW  panel=$PANEL  census=$CENSUS  threads=$THREADS" >> "$RCP"
  echo "# K wall_s core_s peak_rss_gb graph_sha256 graph_bytes" >> "$RCP"
}

peak_rss_gb() { # <pid>: sample the process tree's peak RSS, print GB
  local pid=$1 peak=0 rss
  while kill -0 "$pid" 2>/dev/null; do
    rss=$(ps -o rss= -g "$(ps -o pgid= -p "$pid" 2>/dev/null | tr -d ' ')" 2>/dev/null \
          | awk '{s+=$1} END{print s+0}')
    [ "${rss:-0}" -gt "$peak" ] && peak=$rss
    sleep 5
  done
  awk -v kb="$peak" 'BEGIN{printf "%.2f", kb/1048576}'
}

for K in $KS; do
  M=$SWEEP/K$K.done
  [ -f "$M" ] && { echo "K=$K done (marker)"; continue; }
  OUT=$SWEEP/K$K
  rm -rf "$OUT"
  echo ":: K=$K (threads=$THREADS) starting"
  t0=$(date +%s)
  cpu=$SWEEP/K$K.cpu
  ( TIMEFORMAT='%U %S'; time (
      IMPG_AGC_DECOMPRESSOR_SLOTS=$K RAYON_NUM_THREADS=$THREADS \
        "$IMPG_NEW" genome-infer derive-routes \
          --panel "$PANEL" --sources "$AGC" --read-lengths 148 \
          --census "$CENSUS" --out-dir "$OUT"
    ) ) > "$SWEEP/K$K.log" 2> "$cpu" &
  dp=$!
  peak=$(peak_rss_gb "$dp")
  wait "$dp"; st=$?
  wall=$(( $(date +%s) - t0 ))
  if [ "$st" -ne 0 ]; then
    echo "K=$K FAILED exit=$st wall=${wall}s - log: $SWEEP/K$K.log" | tee -a "$RCP"
    continue
  fi
  # (THE LAST NUMERIC LINE, not every line: the derivation writes its own
  # stderr lines into this capture (e.g. the LEGACY-PANEL OPT-IN notice),
  # and a bare awk over all lines concatenated its per-line prints into
  # garbage like "0.00.054221.9". The time report is the final line and is
  # the only one that is two floats.)
  core=$(awk '/^[0-9]+\.[0-9]+ [0-9]+\.[0-9]+$/{printf "%.1f", $1+$2; found=1} END{if(!found) printf "NA"}' "$cpu" 2>/dev/null)
  sha=$(sha256sum "$OUT/graph.json" | cut -d' ' -f1)
  bytes=$(stat -c %s "$OUT/graph.json")
  echo "K=$K ${wall} ${core} ${peak} ${sha} ${bytes}" >> "$RCP"
  echo ":: K=$K wall=${wall}s core=${core} peak=${peak}GB graph=${bytes}B"
  # THE IDENTITY GATE: every K must emit the same graph.json as the
  # first completed point (the pool shares decompressors; it may not
  # change a byte).
  ref=$(grep -m1 -E '^[0-9]+ [0-9]+ ' "$RCP" | awk '{print $5}')
  if [ -n "$ref" ] && [ "$ref" != "$sha" ]; then
    echo "K=$K IDENTITY GATE FAILED: graph.json sha256 $sha != reference $ref" | tee -a "$RCP"
    exit 1
  fi
  # the workdir self-drains: the receipt is the product, the artifact
  # is a transient (re-derivable; the campaign's own cheap rerun)
  rm -rf "$OUT"
  touch "$M"
  # THE SWEEP'S OWN DISK GATE (NOT the campaign workdir's budget: $EE is
  # a different scratch root and sits at its own bound, which would stop
  # an unrelated run for an unrelated reason). Each point's artifact is
  # small and removed above; this only guards the box's headroom.
  free_g=$(df -BG --output=avail "$SWEEP" | tail -1 | tr -dc '0-9')
  [ "${free_g:-0}" -ge 20 ] || { echo "DISK STOP after K=$K (${free_g}G free)" >> "$RCP"; exit 1; }
done

echo "=== ksweep receipt ($RCP) ==="
cat "$RCP"
echo
echo "READ IT AS: peak_rss_gb should fall with K to a knee; wall_s should"
echo "rise toward the K=1 floor. Pick the knee where wall is still inside"
echo "the once-per-panel budget (~25 min at the committed form)."
