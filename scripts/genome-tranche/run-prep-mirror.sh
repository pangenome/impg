#!/usr/bin/env bash
# THE GENOME TRANCHE - the per-territory substrate MIRROR + REGEN (the
# campaign's own prep chain, run-prep-one.sh VERBATIM, with RV pointed at
# the tranche's own work dir so the campaign's prep/ ledgers are never
# touched): the definitional beds/ledgers (locus.bed, locus.refined.bed,
# axis-path.txt, split.info.json) are COPIED from the runtime
# verification's prep/<T> (the close-out reap kept exactly these), then
# the marker-gated chain re-runs extract -> syng -> beds -> export ->
# walks -> maps in the tranche dir (the ~22s-on-demand regen class; the
# beds segdup tail is the honest risk, walls receipted per step).
# usage: run-prep-mirror.sh <TID>   (env: WT, THREADS)
set -uo pipefail
WT=${WT:-/home/erikg/hprcv2/global-panel-call-mode/genome-tranche1}
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
THREADS=${THREADS:-4}
T=$1
P=$WT/prep/$T
mkdir -p "$P"
for f in locus.bed locus.refined.bed axis-path.txt split.info.json locus.members.bed; do
  if [ ! -f "$P/$f" ] && [ -f "$RV/prep/$T/$f" ]; then cp "$RV/prep/$T/$f" "$P/$f"; fi
done
# the campaign's chain, marker-gated, env-driven (RV=$WT: every path it
# touches lives in the tranche dir; the AGC + the c-rider are absolute)
RV=$WT THREADS=$THREADS bash "$WT/run-prep-one.sh" "$T"
