#!/usr/bin/env bash
# THE CLI-PRODUCT FLEET (the owner's product swap, stage 4): the
# per-locus diplotype-calling CLI run end-to-end over the full
# diploid validation sample (17 components), one component at a time,
# marker-idempotent per step. Per component:
#   1. the census (regenerated through the committed multicensus
#      runner when the receipt is absent - the validation dir was
#      cleaned of the census scratch; the emission is truth-free);
#   2. the CLI product run (the default truth-free calls.jsonl plus
#      the --truth-qv TEST MODE with the validation sample's truth
#      haplotypes as input - the assessment side lives ONLY in the
#      separate calls.jsonl.truth-qv.jsonl artifact);
#   3. the per-component gate vs the SURVIVING committed QV receipts
#      of record (the called folds, the truth ranks, the sequence QV
#      - every field), written to cli-product-gate-<C>.txt.
# The final aggregate lands at cli-product-gate-AGGREGATE.txt (the
# before/after receipt vs the committed assessment numbers).
# Usage: run-cli-product-fleet.sh [COMP ...]   (default: all 15 fleet
# components; chrMT/chrI already gated in the stage-1/2 sessions)
set -uo pipefail
cd /home/erikg/impg-genome-inference || exit 1
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930
PANEL=/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng
ROUTES=/home/erikg/yeast/genome-panel-routes-rebuild-v1/routes
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
IMPG=/home/erikg/impg-genome-inference/target/release/impg
COMPONENTS="${*:-chrII chrIII chrIV chrV chrVI chrVII chrVIII chrIX chrX chrXI chrXII chrXIII chrXIV chrXV chrXVI}"

for C in $COMPONENTS; do
  echo "=== $C $(date -u +%FT%TZ)"
  # step 1: the census (the committed runner re-creates the scratch dir)
  if [ ! -s "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" ]; then
    rm -f "$D/run-cosine-multicensus-pilot-$C.done"
    bash genome/instrumented/run-cosine-multicensus.sh "$C" \
      > "$D/cli-fleet-census-$C.log" 2>&1
    [ -s "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" ] \
      || { echo "$C census FAILED"; exit 1; }
  fi
  echo "=== $C census ready $(date -u +%FT%TZ)"
  # step 2: the CLI product run (the loci from the reference axis)
  if [ ! -f "$D/cli-product-$C.done" ]; then
    read -r NLOCI <<< "$(python3 - "$C" << 'EOF'
import json, sys
c = sys.argv[1]
axis = json.load(open(
    "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
))["intervals"]
print(sum(1 for iv in axis if iv["component"] == f"S288C#0#{c}"))
EOF
)"
    LOCI="$(python3 -c "print(','.join(str(i) for i in range($NLOCI)))")"
    CONTIG="$C"
    rm -rf "$D/cli-product-$C"
    ( cd "$D" && "$IMPG" genome-infer call-diplotypes \
        --panel "$PANEL" --routes "$ROUTES" \
        --partition-graphs "$D/locality-graphs-$C" \
        --component "S288C#0#$C" \
        --census "$D/cosine-diagnostic-scratch/$C/cosine-multi-census.jsonl" \
        --reads "$D/reads.fastq.gz" \
        --derive-cache "$D/realign-derive-cache-chrMT.bin" \
        --loci "$LOCI" \
        --truth-qv "S288C#0#$C" "SK1#0#$CONTIG" \
        --qv-biwfa "$QVBIWFA" \
        --out-dir "$D/cli-product-$C" ) \
      > "$D/cli-product-$C.log" 2>&1
    echo $? > "$D/cli-product-$C.exit"
    [ "$(cat "$D/cli-product-$C.exit")" = 0 ] || { echo "$C CLI run FAILED"; exit 1; }
    touch "$D/cli-product-$C.done"
  fi
  echo "=== $C CLI run done $(date -u +%FT%TZ)"
  # step 3: the gate vs the surviving committed QV receipts
  if [ ! -f "$D/cli-product-gate-$C.done" ]; then
    python3 - "$C" << 'EOF' > "$D/cli-product-gate-$C.txt" 2>&1 \
      || { echo "$C gate FAILED"; exit 1; }
import json, sys
D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
c = sys.argv[1]
qv = {}
for line in open(f"{D}/realign-sequence-qv-locality-{c}.jsonl"):
    r = json.loads(line); qv[r["locus"]] = r
test = {}
for line in open(f"{D}/cli-product-{c}/calls.jsonl.truth-qv.jsonl"):
    r = json.loads(line); test[r["locus"]] = r
FIELDS = ("truth_rank", "rank1", "assignment", "likegt_assignment",
          "assignment_disagreement", "pairs", "chosen_pairs", "total_edits",
          "total_columns", "mismatches", "gap_columns", "error", "identity",
          "qv", "perfect")
diffs = []
for L, r in sorted(qv.items()):
    t = test.get(L)
    if t is None or not t.get("truth_pair_expressible"):
        diffs.append((L, "expressibility", r["truth_rank"], None)); continue
    for k in FIELDS:
        if r.get(k) != t.get(k):
            diffs.append((L, k, repr(r.get(k))[:60], repr(t.get(k))[:60]))
    for field in ("called_folds", "truth_folds"):
        if [sorted(f["strains"]) for f in r[field]] != [sorted(f["strains"]) for f in t[field]]:
            diffs.append((L, field, "strains", ""))
n_expressible = sum(1 for t in test.values() if t.get("truth_pair_expressible"))
n_rank1 = sum(1 for t in test.values() if t.get("rank1"))
print(f"component {c}: loci {len(test)}, expressible {n_expressible}, rank-1 {n_rank1}")
print(f"gate vs committed QV receipts ({len(qv)} expressible loci): "
      + ("FIELD-IDENTICAL" if not diffs else f"{len(diffs)} DIFFS"))
for d in diffs[:20]:
    print("  diff:", d)
EOF
    [ "$(grep -c 'FIELD-IDENTICAL' "$D/cli-product-gate-$C.txt")" -ge 1 ] \
      || { echo "$C gate FAILED (diffs present)"; exit 1; }
    touch "$D/cli-product-gate-$C.done"
  fi
  echo "=== $C GATE CLOSED $(date -u +%FT%TZ)"
done
echo "FLEET COMPLETE $(date -u +%FT%TZ)"
