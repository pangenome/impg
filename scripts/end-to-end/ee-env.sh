#!/usr/bin/env bash
# THE END-TO-END RECEIPT - the shared environment (the three-step
# product at the current stack: panel-prep once -> map per sample ->
# call-region per territory, zero product flags).
#
# THE RUN (the receipt's own scope): HG002 35x vs the ONE committed
# global panel (HPRC_r2_assemblies_0.6.1.syng, metadata v1 - the
# substrate acknowledgment IMPG_PANEL_LEGACY_V1_OK=1, the committed
# campaign's own env), the full accepted territory universe (5,844,
# the runtime-verification split ledger), all committed defaults.
#
# STREAM DISCIPLINE (docs/syng-gmem-bwt/stream-discipline.md, standing
# law): receipts computed at write time; answers kept (calls.jsonl +
# sequences); transients reaped in-pass by the runner's own exit path
# (trap-based, success OR failure OR kill); hard workdir budget
# (<=150G per scratch root) enforced inside every runner; df watched.
# Marker-gated everywhere; a successor RESUMES, never restarts.
set -uo pipefail

EE=${EE:-/home/erikg/hprcv2/end-to-end-receipt}
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
CRIDER=/home/erikg/hprcv2/wgs35x-fleet/c-rider
HPRC=/home/erikg/hprcv2

# THE BINARIES. The current stack (HEAD, built with the clean-PATH env
# into CARGO_TARGET_DIR=target-frontier; NEVER target/release/impg):
# every three-step product stage (panel-prep, map, call-region,
# derive-routes, the cuts) runs on it. The SUBSTRATE instruments (the
# export-profile prep chain: syng build, partition beds, graph export,
# homology walks) stay pinned to the committed target/release binaries
# (unchanged since Oct 19, the instruments the committed territory
# universe's own receipts were built with) - the honest era split,
# stated on every wall.
IMPG_NEW=/home/erikg/impg-genome-inference/target-frontier/release/impg
SP_NEW=/home/erikg/impg-genome-inference/target-frontier/release/examples/sample_projection
IMPG_SUB=/home/erikg/impg-genome-inference/target/release/impg
EX_SUB=/home/erikg/impg-genome-inference/target/release/examples
# the committed projection instrument (unchanged since the committed
# 20.79min receipt; the map step's byte-identity gate demands it)
SP_COMMITTED=/home/erikg/impg-genome-inference/target/release/examples/sample_projection

PANEL=$HPRC/HPRC_r2_assemblies_0.6.1.syng
AGC=$HPRC/HPRC_r2_assemblies_0.6.1.agc
SHARDS=$HPRC/hg002-projection/full
BAM=$HPRC/HG002.GRCh38.35x.bam
L1_AXIS=$HPRC/generic-sampling/L1/reference-axis.json
L1_BEDS=$HPRC/writeup-codification/three-step/partitions/L1/beds
L1_COMPONENT="HG002#1#chr1:4632437-4883593"

# THE WORKDIR BUDGET (the stream law's own bound, per scratch root).
EE_BUDGET_GB=${EE_BUDGET_GB:-150}
# THE WORKER GUARD: 64GiB per call/cut worker (the campaign's own
# guard-as-oracle doctrine; every guard event is a named fail-closed
# census entry, never a silent loss).
GUARD_GIB=${GUARD_GIB:-64}
# the census-producing cut and the routes derive are panel-side
# one-time steps on the whole panel: wider honest guards, stated on
# their receipts (the cut path runs unguarded by design in the
# committed receipts - here it is bounded and the bound is receipted).
ONCE_GUARD_GIB=${ONCE_GUARD_GIB:-96}

export IMPG_PANEL_LEGACY_V1_OK=1
export IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1

mkdir -p "$EE"/{ledger,receipts,logs,prep,partitions,map,reads,calls,chunks,scratch}

# ---- the shared receipts-at-write-time helpers ------------------------

ee_now() { date -u +%FT%TZ; }
ee_log() { echo "$(ee_now) $*" >> "$EE/logs/driver.log"; }

# budget guard: fail the caller when the workdir exceeds the bound
ee_budget_check() {
  local gb
  gb=$(du -sBG "$EE" 2>/dev/null | awk '{ gsub(/G/,""); print int($1) }')
  if [ "${gb:-0}" -ge "$EE_BUDGET_GB" ]; then
    echo "WORKDIR BUDGET EXCEEDED: ${gb}G >= ${EE_BUDGET_GB}G at $EE" >&2
    return 1
  fi
}

# THE ONE STEP RUNNER (marker-gated; wall + core-s receipts at write
# time; resumable): ee_step <dir> <name> <cmd...>
#   -> <dir>/<name>.{done,exit,wall,cpu,log,err}
ee_step() {
  local dir=$1 name=$2; shift 2
  ee_step_out "$dir" "$name" "$dir/$name.log" "$@"
}

# the same, with the command's stdout going to <outfile> instead of
# the step log (for the stream-to-file instruments)
ee_step_out() {
  local dir=$1 name=$2 outfile=$3; shift 3
  local t0 cpu st
  [ -f "$dir/$name.done" ] && return 0
  mkdir -p "$dir"
  t0=$(date +%s)
  cpu=$( { TIMEFORMAT='%U %S'; time "$@" > "$outfile" 2> "$dir/$name.err"; } 2>&1 )
  st=$?
  echo "$(( $(date +%s) - t0 ))" > "$dir/$name.wall"
  echo "$cpu" > "$dir/$name.cpu"
  if [ "$st" -eq 0 ]; then
    touch "$dir/$name.done"
  else
    echo "$st" > "$dir/$name.exit"
  fi
  return "$st"
}

# the RSS guard poller: <pid> <guard_gib> <receipt-file> [label]
# kills the pid when its VmRSS crosses the guard; the event is
# receipted (the fail-closed census's own entry). Prints the poller pid.
ee_guard() {
  local pid=$1 gib=$2 receipt=$3 label=${4:-worker}
  local kb=$((gib*1024*1024))
  (
    while kill -0 "$pid" 2>/dev/null; do
      local rss
      rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || { sleep 10; continue; }
      if [ "${rss:-0}" -gt "$kb" ]; then
        echo "$(ee_now) GUARD kill $label pid=$pid rss=${rss}kB > ${gib}GiB" >> "$receipt"
        kill -9 "$pid" 2>/dev/null
        break
      fi
      sleep 10
    done
  ) &
  echo $!
}

# the guard poller over a whole invocation tree (every direct
# descendant worker over the guard is killed - the per-worker
# doctrine at the fleet shape; the parent then fails at its own
# first-failure discipline and the leg is receipted fail-closed).
ee_guard_tree() {
  local root=$1 gib=$2 receipt=$3 label=${4:-fleet}
  local kb=$((gib*1024*1024))
  (
    while kill -0 "$root" 2>/dev/null; do
      local pid rss
      for pid in $(ps -eo pid,ppid --no-headers | awk -v r="$root" '$2==r{print $1}'); do
        rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
        if [ "${rss:-0}" -gt "$kb" ]; then
          echo "$(ee_now) GUARD kill $label pid=$pid rss=${rss}kB > ${gib}GiB" >> "$receipt"
          kill -9 "$pid" 2>/dev/null
        fi
      done
      sleep 10
    done
  ) &
  echo $!
}

# peak-RSS tracker: <pid> <receipt-file> [label] - records the max
ee_rss_track() {
  local pid=$1 receipt=$2 label=${3:-worker}
  (
    local peak=0 rss
    while kill -0 "$pid" 2>/dev/null; do
      rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || { sleep 10; continue; }
      [ "${rss:-0}" -gt "$peak" ] && peak=$rss
      sleep 10
    done
    echo "$(ee_now) RSSPEAK $label pid=$pid peak=${peak}kB" >> "$receipt"
  ) &
  echo $!
}
