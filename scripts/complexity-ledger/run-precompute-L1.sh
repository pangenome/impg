#!/usr/bin/env bash
# THE PANEL-SIDE PER-LOCUS PRECOMPUTE RUNNER AT THE STANDARD TARGET
# (SS10.6.23 slice 1's answer gates at L1, the five-locality chain's
# committed 10kb recipe — run-memo-L1.sh's own substrate shape VERBATIM,
# rebuilt transiently in bounded scratch from the READ-ONLY source assets
# at /home/erikg/hprcv2/generic-sampling/L1), marker-gated, budget-guarded
# (150G), receipts at write time, reaped by the worker's checkpoint
# (`reap` subcommand).
#
# THE LEGS: precompute (the ONCE-PER-PANEL price — the precompute step's
# own wall + asset size, receipted beside the asset it builds), call-off
# (the new binary's LIVE path, IMPG_PANEL_PRECOMPUTE=off — the
# toggle-restore gate + the before-timing), call-on (the SERVE — the
# after-timing), BOTH against the committed memo receipts' answer digests
# (calls.jsonl AND truth-qv sha256) and the committed fold-set digests,
# parsed from docs/syng-gmem-bwt/memoization-receipts/call-base.receipt.txt
# (the committed record of record).
# usage: run-precompute-L1.sh [reap]
set -uo pipefail
G=/home/erikg/hprcv2/generic-sampling
WORK=/mnt/nvme2n1/erikg/precompute-L1
KEEP=/mnt/nvme2n1/erikg/precompute-L1-keep
WORKDIR_BUDGET_GB=150
REPO=/home/erikg/impg-genome-inference
# the precompute increment's binaries: PRE = the new build (the serve +
# the toggle); BASE = the committed HEAD rebuild (aaefecc) the substrate
# stages pin to (the one-era discipline).
IMPG_PRE=$REPO/target-frontier/release/impg-precompute-wip
IMPG_NEW=$REPO/target-frontier/release/impg-baseline-aaefecc
IMPG_BASE=$REPO/target-frontier/release/impg-baseline-aaefecc
COMMITTED_RECEIPT=$REPO/docs/syng-gmem-bwt/memoization-receipts/call-base.receipt.txt
IMPG_SUB=/home/erikg/impg-genome-inference/target/release/impg
EXAMPLES=/home/erikg/impg-genome-inference/target/release/examples
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
P0=$G/L1
READS="$P0/hg002.L1.300x.all.fq.gz"
P=$WORK
W=10000

do_reap() {
  echo "$(date -u +%FT%TZ) reap: deleting the transient scratch $WORK" | tee -a "$KEEP/reap.log"
  rm -rf "$WORK"
  echo "$(date -u +%FT%TZ) scratch reaped; keep dir: $KEEP ($(du -sh "$KEEP" 2>/dev/null | cut -f1))" | tee -a "$KEEP/reap.log"
}
if [ "${1:-}" = "reap" ]; then do_reap; exit 0; fi
mkdir -p "$WORK" "$KEEP"

budget_check() {
  local gb; gb=$(( $(du -sk "$WORK" 2>/dev/null | cut -f1) / 1048576 ))
  echo "$(date -u +%FT%TZ) workdir ${gb}G of ${WORKDIR_BUDGET_GB}G budget" >> "$WORK/budget.log"
  [ "$gb" -lt "$WORKDIR_BUDGET_GB" ] || { echo "BUDGET EXCEEDED ${gb}G" | tee -a "$WORK/budget.log"; exit 3; }
}
guard_poll() { local pid=$1 rssf=$2 cap_kb=${3:-$((64*1024*1024))}
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt "$cap_kb" ]; then echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$rssf"; kill -9 "$pid" 2>/dev/null; break; fi
      sleep 15
    done ) &
}
runjob() { local name=$1 cap_kb=$2; shift 2
  if [ -f "$WORK/$name.done" ]; then echo "=== $name done (marker)"; return 0; fi
  budget_check
  local start pid st; start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$@" > "$WORK/$name.log" 2> "$WORK/$name.err" &
  pid=$!; guard_poll "$pid" "$WORK/$name.rss" "$cap_kb"; local poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$WORK/$name.exit"; echo "$(( $(date +%s) - start ))" > "$WORK/$name.wall"
  [ "$st" -eq 0 ] && touch "$WORK/$name.done" || { echo "=== $name FAILED ($st)"; return 1; }
}

BEDS="$P/partition-w$((W/1000))k-d$((W/10000))k/results"; ROUTES="$P/routes"; LGRAPHS="$P/locality-graphs"
RAX="$P/reference-axis.json"
CENSUS="$P/cosine-multi-census.jsonl"; DCACHE="$P/realign-derive-cache.bin"
PANEL="$P/locus-pos64.syng"
AXIS_PATH=$(awk -F'\t' '$1 ~ /^HG002#1#/ {print $1":"$2"-"$3; exit}' "$P0/locus.refined.bed")
TRUTH2_PATH=$(awk -F'\t' '$1 ~ /^HG002#2#/ {print $1":"$2"-"$3; exit}' "$P0/locus.refined.bed")

runjob s0-syng $((64*1024*1024)) "$IMPG_NEW" syng -f "$P0/locus.fa" --syncmer-length 63 --smer-length 8 \
  --syncmer-seed 7 --position-sample-rate 64 -o "$PANEL" || exit 1
mkdir -p "$BEDS"
printf '%s\n' "$AXIS_PATH" > "$P/partition-starting-sequences.txt"
runjob s1-partition $((64*1024*1024)) "$IMPG_NEW" partition \
  -a "$PANEL" -w $W -d $((W/10)) \
  --starting-sequences-file "$P/partition-starting-sequences.txt" \
  --selection-mode haplotype \
  --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
  --no-rehome-singletons \
  --sequence-files "$P0/locus.fa" \
  -o bed --separate-files --output-folder "$BEDS" \
  -t 16 || exit 1
NP=$(ls "$BEDS" | grep -c '\.bed$')
runjob s2-catalog $((64*1024*1024)) "$IMPG_NEW" genome-infer build-catalog \
  --bed-dir "$BEDS" --panel "$PANEL" --out-dir "$P/catalog" || exit 1
runjob s2b-routes $((64*1024*1024)) "$IMPG_NEW" genome-infer build-panel-routes \
  --panel "$PANEL" --catalog "$P/catalog/catalog.json" \
  --sources "$P0/locus.fa" --read-lengths 148 --out-dir "$ROUTES" || exit 1
rm -f "$P/catalog/catalog.json"
IDS=$(seq 0 $((NP-1)))
mkdir -p "$P/export-graphs"
runjob s3-export $((64*1024*1024)) "$EXAMPLES/partition_graph_export" export "$PANEL" "$P0/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
AXIS_PIDS=$(for i in $(seq 0 $((NP-1))); do grep -qF "$AXIS_PATH" "$BEDS/partition$i.bed" 2>/dev/null && echo $i; done | tr '\n' ' ')
if [ ! -f "$P/s4-pggb.done" ]; then
  budget_check
  start=$(date +%s)
  LOCUS_DIR="$P0" BED_ROOT="$BEDS" PGGB_OUT="$P/pggb-graphs" \
    bash "$G/run-pggb-builds-generic.sh" $AXIS_PIDS > "$P/s4-pggb.log" 2>&1
  bad=0; for i in $AXIS_PIDS; do
    [ "$(cat "$P/pggb-graphs/partition$i/graph.gfa.exit" 2>/dev/null)" = "0" ] || bad=1; done
  [ "$bad" = 0 ] && { touch "$P/s4-pggb.done"; echo "$(( $(date +%s) - start ))" > "$P/s4-pggb.wall"; } \
    || { echo "=== pggb FAILED"; exit 1; }
fi
for i in $AXIS_PIDS; do cp "$BEDS/partition$i.bed" "$P/pggb-graphs/partition$i/members.bed" 2>/dev/null; done
mkdir -p "$P/pggb-projection"
runjob s5-projection $((64*1024*1024)) "$EXAMPLES/pggb_projection" --panel "$PANEL" --agc "$P0/locus.fa" \
  --export-graphs "$P/export-graphs" --pggb-dir "$P/pggb-graphs" \
  --out-dir "$P/pggb-projection" --partitions "$(echo $AXIS_PIDS | tr ' ' ',')" || exit 1
if [ ! -f "$P/s6-maps.done" ]; then
  budget_check
  start=$(date +%s)
  python3 "$G/locality-maps.py" "$AXIS_PATH" "$BEDS" "$P/export-graphs" \
    "$P/pggb-graphs" "$P/pggb-projection" "$LGRAPHS" \
    "$RAX" "$P" "$P/locality-domain-structure.txt" \
    > "$P/s6-maps.log" 2> "$P/s6-maps.err" \
    && { touch "$P/s6-maps.done"; echo "$(( $(date +%s) - start ))" > "$P/s6-maps.wall"; } \
    || { echo "=== maps FAILED"; exit 1; }
fi
NW=$(ls "$LGRAPHS" | grep -c 'window[0-9]*.gfa.map.json')
echo "=== L1: $NW windows"
if [ ! -f "$P/s7-census.done" ]; then
  budget_check
  start=$(date +%s)
  IMPG_FORCE_DIPLOID_DIAGNOSTIC=1 IMPG_MULTIPLICITY_VARIANT=gentle-beta \
  IMPG_COSINE_READ_MATCHED=1 \
  IMPG_COSINE_MULTI_CENSUS="$CENSUS" \
  IMPG_COSINE_DOMAIN_PILOT=1 \
    taskset -c 0-63 nice -n 10 \
    stdbuf -oL -eL "$EXAMPLES/panel_route_mem_routed" \
    --panel "$PANEL" --routes "$ROUTES" \
    --sample "$P0/sample/sample.membwt" \
    --reads "$READS" \
    --axis "$RAX" \
    --bed-directory "$BEDS" \
    --component "$AXIS_PATH" \
    --depth 137.3 --read-length 148 \
    --routing-universe component --spine --spine-phasing \
    --reference-pair "truth:$P0/truth-routes/HG002-1.json:$P0/truth-routes/HG002-2.json" \
    > "$P/s7-census.log" 2> "$P/s7-census.err" &
  cpid=$!; guard_poll "$cpid" "$P/s7-census.rss" $((64*1024*1024)); cpoll=$!
  cst=0; wait "$cpid" || cst=$?; kill "$cpoll" 2>/dev/null; wait "$cpoll" 2>/dev/null
  echo "$(( $(date +%s) - start ))" > "$P/s7-census.wall"
  [ "$cst" -eq 0 ] && touch "$P/s7-census.done" || { echo "=== census FAILED ($cst)"; exit 1; }
fi
if [ ! -f "$P/s8-derive-cache.done" ]; then
  budget_check
  python3 -c "
import sys
for i in range($NW): print('{\"locus\": %d, \"records\": []}' % i)" > "$P/ingredients-stub.jsonl"
  start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$EXAMPLES/partition_realign_likelihood" \
    --panel "$PANEL" --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --ingredients "$P/ingredients-stub.jsonl" \
    --reads "$READS" \
    --derive-cache "$DCACHE" \
    --read-length 148 \
    --out "$P/realign-stub-out.jsonl" \
    > "$P/s8-derive-cache.log" 2> "$P/s8-derive-cache.err" &
  pid=$!; guard_poll "$pid" "$P/s8-derive-cache.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$(( $(date +%s) - start ))" > "$P/s8-derive-cache.wall"
  if [ "$st" -eq 0 ] || { [ -s "$DCACHE" ] && grep -q 'the read hull differs from the placed path segment' "$P/s8-derive-cache.err"; }; then
    touch "$P/s8-derive-cache.done"
  else echo "=== derive cache FAILED ($st)"; exit 1; fi
fi
LOCI=$(python3 -c "
import json
a=json.load(open('$RAX'))
n=sum(1 for iv in a['intervals'] if iv['component'] == '$AXIS_PATH')
print(','.join(str(i) for i in range(n)))")
echo "=== CLI loci: $LOCI"

run_call() { local tag=$1 bin=$2 pre_env=${3:-auto}
  if [ -f "$P/call-$tag.done" ]; then echo "=== call-$tag done (marker)"; return 0; fi
  budget_check
  rm -rf "$P/call-$tag"
  local start pid st
  start=$(date +%s)
  local pre_env_args=()
  [ "$pre_env" = "off" ] && pre_env_args+=("IMPG_PANEL_PRECOMPUTE=off")
  [ "$pre_env" = "on" ] && pre_env_args+=("IMPG_PANEL_PRECOMPUTE=on")
  env RAYON_NUM_THREADS=1 IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 IMPG_PANEL_LEGACY_V1_OK=1 \
    IMPG_LOCUS_SCORING_CENSUS=1 IMPG_LOCUS_PHASE_CENSUS=1 IMPG_FOLD_SET_DIGEST=1 \
    "${pre_env_args[@]}" \
    nice -n 10 stdbuf -oL -eL "$bin" genome-infer call-diplotypes \
    --panel "$PANEL" \
    --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    --census "$CENSUS" \
    --reads "$READS" \
    --derive-cache "$DCACHE" \
    --loci "$LOCI" \
    --truth-qv "$AXIS_PATH" "$TRUTH2_PATH" \
    --far-bin-prune=false \
    --qv-biwfa "$QVBIWFA" \
    --out-dir "$P/call-$tag" \
    > "$P/call-$tag.log" 2> "$P/call-$tag.err" &
  pid=$!; guard_poll "$pid" "$P/call-$tag.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/call-$tag.exit"; echo "$(( $(date +%s) - start ))" > "$P/call-$tag.wall"
  {
    echo "# call-$tag receipt ($(date -u +%FT%TZ)); binary $bin; loci $LOCI"
    echo "# wall_s $(cat "$P/call-$tag.wall") exit $(cat "$P/call-$tag.exit")"
    if [ -f "$P/call-$tag/calls.jsonl" ]; then
      echo "calls_sha256 $(sha256sum "$P/call-$tag/calls.jsonl" | cut -d' ' -f1) $(stat -c%s "$P/call-$tag/calls.jsonl") bytes"
    fi
    if [ -f "$P/call-$tag/calls.jsonl.truth-qv.jsonl" ]; then
      echo "truthqv_sha256 $(sha256sum "$P/call-$tag/calls.jsonl.truth-qv.jsonl" | cut -d' ' -f1) $(stat -c%s "$P/call-$tag/calls.jsonl.truth-qv.jsonl") bytes"
    fi
    grep -h "panel-precompute" "$P/call-$tag.err" || true
    grep -h "phase inputs" "$P/call-$tag.err" || true
    grep -h "phase binding" "$P/call-$tag.err" || true
    grep -h "phase rows" "$P/call-$tag.err" || true
    grep -h "locus .* folds.*s\\]$" "$P/call-$tag.err" || true
    grep -h "phase census" "$P/call-$tag.err" || true
    grep -h "alignment-connected folds" "$P/call-$tag.err" || true
    grep -h "scoring census" "$P/call-$tag.err" || true
    grep -h "locus .* done:" "$P/call-$tag.err" || true
    python3 - "$P/call-$tag" <<'EOF'
import json, sys, os
FNV_OFFSET = 0xCBF29CE484222325
FNV_PRIME = 0x100000001B3
def fnv_digest(folds):
    d = FNV_OFFSET
    def feed(b):
        nonlocal d
        for byte in b:
            d ^= byte
            d = (d * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
        
    for fold in folds:
        feed(fold["length"].to_bytes(8, "little"))
        for m in fold["members"]:
            feed(m["path_name"].encode())
            feed(m["start"].to_bytes(8, "little"))
            feed(m["end"].to_bytes(8, "little"))
    return f"{d:016x}"
d = sys.argv[1]
path = os.path.join(d, "instrument-receipt.jsonl")
if os.path.exists(path):
    for line in open(path, errors="replace"):
        try: r = json.loads(line)
        except Exception: continue
        if not isinstance(r, dict) or "locus" not in r: continue
        fsd = r.get("fold_set_digest")
        if fsd:
            print(f"fold_set_digest locus {r['locus']} folds {fsd['folds']} {fsd['digest']} (in-process)")
        else:
            fi = r.get("fold_identities")
            if fi:
                print(f"fold_set_digest locus {r['locus']} folds {len(fi)} {fnv_digest(fi)} (offline-mirror)")
        pc = r.get("phase_census")
        if pc: print(f"phase_census locus {r['locus']} " + " ".join(
            f"{k}={v}" for k, v in sorted(pc.items()) if k != "note"))
        w = r.get("walls")
        if w: print(f"walls locus {r['locus']} " + " ".join(f"{k}={v}" for k, v in sorted(w.items())))
EOF
  } > "$WORK/call-$tag.receipt.txt" 2>&1
  cp -a "$WORK/call-$tag.receipt.txt" "$KEEP/"
  [ "$st" -eq 0 ] && [ -s "$P/call-$tag/calls.jsonl" ] && touch "$P/call-$tag.done" && return 0
  echo "=== call-$tag FAILED ($st)"; return 1
}

# ---- THE PRECOMPUTE STEP (the once-per-panel price, receipted): the
# row context + the fold layouts, keyed on content into the locality
# graphs' own panel-precompute-loci section.
if [ ! -f "$P/precompute.done" ]; then
  budget_check
  start=$(date +%s)
  env IMPG_PANEL_LEGACY_V1_OK=1 RAYON_NUM_THREADS=1 \
    nice -n 10 stdbuf -oL -eL "$IMPG_PRE" genome-infer precompute-panel-loci \
    --panel "$PANEL" --routes "$ROUTES" \
    --partition-graphs "$LGRAPHS" \
    --component "$AXIS_PATH" \
    > "$P/precompute.log" 2> "$P/precompute.err" &
  pid=$!; guard_poll "$pid" "$P/precompute.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/precompute.exit"; echo "$(( $(date +%s) - start ))" > "$P/precompute.wall"
  {
    echo "# precompute receipt ($(date -u +%FT%TZ)); binary $IMPG_PRE; component $AXIS_PATH"
    echo "# wall_s $(cat "$P/precompute.wall") exit $(cat "$P/precompute.exit")"
    cat "$P/precompute.err"
    echo "asset_bytes $(du -sb "$LGRAPHS/panel-precompute-loci" 2>/dev/null | cut -f1)"
  } > "$WORK/precompute.receipt.txt" 2>&1
  cp -a "$WORK/precompute.receipt.txt" "$KEEP/"
  [ "$st" -eq 0 ] && touch "$P/precompute.done" || { echo "=== precompute FAILED ($st)"; exit 1; }
fi

# THE PRECOMPUTE INCREMENT'S LEGS at the standard target: the LIVE path
# (the toggle-restore gate + the before-timing), the SERVE (the
# after-timing), BOTH against the committed memo receipts' digests.
run_call precompute-off "$IMPG_PRE" off || exit 1
run_call precompute-on "$IMPG_PRE" on || exit 1

# ---- THE GATES: answers + fold-set digests vs the committed receipts
{
  committed_calls=$(grep -m1 '^calls_sha256' "$COMMITTED_RECEIPT" | awk '{print $2}')
  committed_truth=$(grep -m1 '^truthqv_sha256' "$COMMITTED_RECEIPT" | awk '{print $2}')
  pass=1
  for tag in precompute-off precompute-on; do
    c=$(grep -m1 '^calls_sha256' "$WORK/call-$tag.receipt.txt" | awk '{print $2}')
    t=$(grep -m1 '^truthqv_sha256' "$WORK/call-$tag.receipt.txt" | awk '{print $2}')
    [ "$c" = "$committed_calls" ] && echo "GATE calls $tag PASS ($c)" || { echo "GATE calls $tag FAIL ($c != $committed_calls)"; pass=0; }
    [ "$t" = "$committed_truth" ] && echo "GATE truth-qv $tag PASS ($t)" || { echo "GATE truth-qv $tag FAIL ($t != $committed_truth)"; pass=0; }
  done
  for tag in precompute-off precompute-on; do
    for locus_digest in $(grep '^fold_set_digest' "$COMMITTED_RECEIPT" | awk '{print $3":"$6}'); do
      locus=${locus_digest%%:*}; digest=${locus_digest##*:}
      got=$(grep -m1 "fold_set_digest locus $locus " "$WORK/call-$tag.receipt.txt" | awk '{print $6}')
      [ "$got" = "$digest" ] || { echo "GATE fold digest $tag locus $locus FAIL ($got != $digest)"; pass=0; }
    done
    n_committed=$(grep -c '^fold_set_digest' "$COMMITTED_RECEIPT")
    n_got=$(grep -c '^fold_set_digest' "$WORK/call-$tag.receipt.txt")
    [ "$n_got" = "$n_committed" ] && echo "GATE fold digests $tag PASS ($n_got/$n_committed identical)" || { echo "GATE fold digests $tag FAIL ($n_got vs $n_committed)"; pass=0; }
  done
  [ "$pass" = 1 ] && echo "=== L1 PRECOMPUTE GATES: ALL PASS" || echo "=== L1 PRECOMPUTE GATES: FAILURES (above)"
} | tee "$KEEP/gates.txt"
echo "=== L1 PRECOMPUTE GATES COMPLETE - receipts in $KEEP; reap with: $0 reap"
