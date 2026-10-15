#!/usr/bin/env bash
# THE DOSAGE-SURFACE FLEET: the copy-number layer over the chain
# layer's two molecules for every component of the validation sample
# (the chrMT/chrI gates already stand), then the per-component gate.
# Marker-idempotent per component:
#   dosage-<C>/                     the dosage run (dosage.jsonl +
#                                   the test-mode artifact)
#   dosage-gate-<C>.txt/.done      the gate (the thin-layer proof +
#                                   the QC arithmetic + the agreement
#                                   table)
# Usage: run-dosage-fleet.sh [COMP ...]  (default: the 15 remaining)
set -uo pipefail
cd /home/erikg/impg-genome-inference || exit 1
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
PANEL=/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng
IMPG=/home/erikg/impg-genome-inference/target/release/impg
COMPONENTS="${*:-chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI}"

truth_qv_of() {
  # chrMT's truth-QV artifact lives at the dedicated truthqv run of
  # record (the gate component's own demonstration run).
  case "$1" in
    chrMT) echo "$D/cli-product-chrMT-truthqv/calls.jsonl.truth-qv.jsonl" ;;
    *)     echo "$D/cli-product-$1/calls.jsonl.truth-qv.jsonl" ;;
  esac
}
folds_of() {
  case "$1" in
    chrMT) echo "$D/cli-product-chrMT-truthqv/instrument-receipt.jsonl" ;;
    *)     echo "$D/cli-product-$1/instrument-receipt.jsonl" ;;
  esac
}

for C in $COMPONENTS; do
  echo "=== $C $(date -u +%FT%TZ)"
  if [ ! -f "$D/dosage-$C.done" ]; then
    rm -rf "$D/dosage-$C"
    stdbuf -oL -eL "$IMPG" genome-infer emit-dosage \
      --molecules "$D/chain-molecules-$C/molecules.jsonl" \
      --panel "$PANEL" \
      --partition-graphs "$D/locality-graphs-$C" \
      --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
      --truth-qv-file "$(truth_qv_of "$C")" \
      --folds "$(folds_of "$C")" \
      --out-dir "$D/dosage-$C" \
      > "$D/dosage-$C.log" 2>&1
    echo $? > "$D/dosage-$C.exit"
    [ "$(cat "$D/dosage-$C.exit")" = 0 ] \
      || { echo "$C dosage run FAILED"; exit 1; }
    touch "$D/dosage-$C.done"
  fi
  echo "=== $C dosage done $(date -u +%FT%TZ)"
  if [ ! -f "$D/dosage-gate-$C.done" ]; then
    python3 -B genome/instrumented/dosage-gate.py "$C" \
      "$D/chain-molecules-$C/molecules.jsonl" \
      "$D/dosage-$C" \
      "$D/locality-graphs-$C" \
      "$(truth_qv_of "$C")" \
      > "$D/dosage-gate-$C.txt" 2>&1
    echo $? > "$D/dosage-gate-$C.exit"
    [ "$(cat "$D/dosage-gate-$C.exit")" = 0 ] \
      || { echo "$C gate FAILED"; exit 1; }
    touch "$D/dosage-gate-$C.done"
  fi
  echo "=== $C gate done $(date -u +%FT%TZ)"
done
touch "$D/dosage-fleet.done"
echo "=== FLEET COMPLETE $(date -u +%FT%TZ)"
