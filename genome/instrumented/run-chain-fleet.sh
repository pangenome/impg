#!/usr/bin/env bash
# THE CHAIN-MOLECULES FLEET: the chain/phasing layer over the CLI
# product's closed calls for every component of the validation
# sample (the chrMT/chrI gates already stand), then the per-component
# gate. Marker-idempotent per component:
#   chain-molecules-<C>/            the chain run (molecules.jsonl +
#                                  the test-mode artifact)
#   chain-gate-<C>.txt/.done       the gate (the thin-layer proof +
#                                  the switch table)
# Usage: run-chain-fleet.sh [COMP ...]  (default: the 15 remaining)
set -uo pipefail
cd /home/erikg/impg-genome-inference || exit 1
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
PANEL=/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng
ROUTES=/home/erikg/yeast/genome-panel-routes-rebuild-v1/routes
IMPG=/home/erikg/impg-genome-inference/target/release/impg
COMPONENTS="${*:-chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI}"

for C in $COMPONENTS; do
  echo "=== $C $(date -u +%FT%TZ)"
  if [ ! -f "$D/chain-molecules-$C.done" ]; then
    rm -rf "$D/chain-molecules-$C"
    stdbuf -oL -eL "$IMPG" genome-infer chain-molecules \
      --calls "$D/cli-product-$C/calls.jsonl" \
      --panel "$PANEL" --routes "$ROUTES" \
      --partition-graphs "$D/locality-graphs-$C" \
      --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
      --truth-qv-file "$D/cli-product-$C/calls.jsonl.truth-qv.jsonl" \
      --out-dir "$D/chain-molecules-$C" \
      > "$D/chain-molecules-$C.log" 2>&1
    echo $? > "$D/chain-molecules-$C.exit"
    [ "$(cat "$D/chain-molecules-$C.exit")" = 0 ] \
      || { echo "$C chain run FAILED"; exit 1; }
    touch "$D/chain-molecules-$C.done"
  fi
  echo "=== $C chain done $(date -u +%FT%TZ)"
  if [ ! -f "$D/chain-gate-$C.done" ]; then
    python3 -B genome/instrumented/chain-molecules-gate.py "$C" \
      "$D/cli-product-$C/calls.jsonl" \
      "$D/chain-molecules-$C" \
      "$D/cli-product-$C/calls.jsonl.truth-qv.jsonl" \
      > "$D/chain-gate-$C.txt" 2>&1
    echo $? > "$D/chain-gate-$C.exit"
    [ "$(cat "$D/chain-gate-$C.exit")" = 0 ] \
      || { echo "$C gate FAILED"; exit 1; }
    touch "$D/chain-gate-$C.done"
  fi
  echo "=== $C gate done $(date -u +%FT%TZ)"
done
touch "$D/chain-fleet.done"
echo "=== FLEET COMPLETE $(date -u +%FT%TZ)"
