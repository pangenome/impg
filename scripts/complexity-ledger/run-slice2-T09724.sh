#!/usr/bin/env bash
# THE BIND-SCAFFOLDING SLICE RUNNER (SS10.6.25 slice 2's gates +
# timing at the fold-heavy tranche-2 member T09724 — the committed
# global-frame callC form, the run-memo-T09724.sh substrate shape VERBATIM:
# extract -> syng -> partition w10k d1k -> export -> walks -> maps -> reads
# -> rvcut -> cutA -> routes), marker-gated and resumable, budget-guarded
# (150G), receipts at write time, reaped by the worker's checkpoint
# (`reap`).
#
# THE LEGS: precompute (the ONCE-PER-PANEL price at the GLOBAL frame),
# call-off (the new binary's LIVE path, IMPG_PANEL_PRECOMPUTE=off — the
# toggle-restore gate + the before-timing), call-on (the SERVE), BOTH
# against the committed memo receipts' answer digests + fold-set digests
# (docs/syng-gmem-bwt/memoization-receipts/call-base-all.receipt.txt).
# usage: run-slice2-T09724.sh [reap]
#
# THE LEGS: base-l9/memo-l9 (the committed form's class-scan term before/after
# at the fold-heavy locus), base-all/memo-all (the whole territory: the
# calls/truth-qv sha256 answer-digest gates + the per-locus memo hit-rate
# table under the phase census), memo-off-all (the toggle-restore gate: the
# memo binary with IMPG_CLASS_SCAN_MEMO=off must equal the direct form byte
# for byte), memo-premerge-all (the pre-merge anatomy form — the unmerged
# bit-equal rows are where the memo's classes concentrate; digest vs the
# committed mergefix pre-merge receipts).
set -uo pipefail
T=T09724
HPRC=/home/erikg/hprcv2
MEMBERS=$HPRC/panel/tranche2/members/$T
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
WORK=/mnt/nvme2n1/erikg/slice2-$T
KEEP=/mnt/nvme2n1/erikg/slice2-$T-keep
WORKDIR_BUDGET_GB=150
REPO=/home/erikg/impg-genome-inference
# the precompute increment's binaries: PRE = the new build (the serve +
# the toggle); the substrate stages pin to the committed binaries (the
# one-era discipline)
IMPG_PRE=$REPO/target-frontier/release/impg-slice2-wip
COMMITTED_RECEIPT=$REPO/docs/syng-gmem-bwt/memoization-receipts/call-base-all.receipt.txt
IMPG_SUB=/home/erikg/impg-genome-inference/target/release/impg
EXAMPLES=/home/erikg/impg-genome-inference/target/release/examples
AGC_BIN=/home/erikg/agc/bin/agc
AGC=$HPRC/HPRC_r2_assemblies_0.6.1.agc
GSYNG=$HPRC/HPRC_r2_assemblies_0.6.1.syng
SHARDS=$HPRC/hg002-projection/full
BAM=$HPRC/HG002.GRCh38.35x.bam
QVBIWFA=/home/erikg/impg-genome-inference/genome/instrumented/qv-biwfa/target/release/qv-biwfa
AXIS=$(cat "$MEMBERS/axis-path.txt")
TRUTH2=$(awk -F'\t' -v t=$T '$1==t {print $3}' $HPRC/panel/tranche2/truth2.tsv)
P=$WORK/prep

do_reap() {
  echo "$(date -u +%FT%TZ) reap: deleting the transient scratch $WORK" | tee -a "$KEEP/reap.log"
  rm -rf "$WORK"
  echo "$(date -u +%FT%TZ) scratch reaped; keep dir: $KEEP ($(du -sh "$KEEP" 2>/dev/null | cut -f1))" | tee -a "$KEEP/reap.log"
}
if [ "${1:-}" = "reap" ]; then do_reap; exit 0; fi

mkdir -p "$WORK" "$KEEP" "$P"

budget_check() {
  local gb
  gb=$(( $(du -sk "$WORK" 2>/dev/null | cut -f1) / 1048576 ))
  echo "$(date -u +%FT%TZ) workdir ${gb}G of ${WORKDIR_BUDGET_GB}G budget" >> "$WORK/budget.log"
  if [ "$gb" -ge "$WORKDIR_BUDGET_GB" ]; then
    echo "$(date -u +%FT%TZ) BUDGET EXCEEDED: ${gb}G >= ${WORKDIR_BUDGET_GB}G - stopping fail-closed" | tee -a "$WORK/budget.log"
    exit 3
  fi
}

guard_poll() { local pid=$1 rssf=$2 cap_kb=${3:-$((300*1024*1024))}
  ( while kill -0 "$pid" 2>/dev/null; do
      rss_kb=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null) || continue
      echo "$(date -u +%FT%TZ) $pid ${rss_kb:-0}" >> "$rssf"
      if [ "${rss_kb:-0}" -gt "$cap_kb" ]; then echo "$(date -u +%FT%TZ) GUARD kill ${rss_kb}kB" >> "$rssf"; kill -9 "$pid" 2>/dev/null; break; fi
      sleep 15
    done ) &
}

stage() { local name=$1 cap_kb=$2; shift 2
  if [ -f "$P/$name.done" ]; then echo "=== $name done (marker)"; return 0; fi
  budget_check
  local start pid st; start=$(date +%s)
  nice -n 10 stdbuf -oL -eL "$@" > "$P/$name.log" 2> "$P/$name.err" &
  pid=$!; guard_poll "$pid" "$P/$name.rss" "$cap_kb"; local poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/$name.exit"; echo "$(( $(date +%s) - start ))" > "$P/$name.wall"
  [ "$st" -eq 0 ] && touch "$P/$name.done" || { echo "=== $name FAILED ($st)"; return 1; }
  return 0
}

# ---- stage extract: locus.fa from the KEPT extract specs (the committed
# AGC getctg recipe; the kept extract.specs.txt IS the committed artifact)
if [ ! -f "$P/extract.done" ]; then
  budget_check
  python3 - "$MEMBERS" "$P" "$AGC_BIN" "$AGC" <<'EOF'
import subprocess, sys, os
members_dir, P, agc_bin, agc = sys.argv[1:5]
specs = [l.strip() for l in open(f"{members_dir}/extract.specs.txt") if l.strip()]
wanted = []
for spec in specs:
    path, span = spec.split("@", 1)
    lo, hi = span.split(":")[-1].split("-")
    # the spec span is 1-based inclusive; the bed row is 0-based half-open
    wanted.append((path, int(lo) - 1, int(hi)))
raw = f"{P}/extract.raw.fa"
with open(f"{P}/extract.specs.txt", "w") as f:
    f.write("\n".join(specs) + "\n")
subprocess.run(f"xargs -a {P}/extract.specs.txt {agc_bin} getctg {agc} > {raw}",
              shell=True, check=True, timeout=3600)
recs, buf = [], []
for line in open(raw):
    if line.startswith(">"):
        if buf: recs.append("".join(buf))
        buf = []
    elif line.strip(): buf.append(line.strip())
if buf: recs.append("".join(buf))
assert len(recs) == len(wanted), f"EXTRACTION COUNT MISMATCH {len(recs)} vs {len(wanted)}"
total = 0
with open(f"{P}/locus.fa", "w") as fa, open(f"{P}/locus.members.bed", "w") as mb:
    for (path, s, e), seq in zip(wanted, recs):
        assert len(seq) == e - s, f"{path}:{s}-{e}: len {len(seq)} != bed {e - s}"
        total += len(seq)
        fa.write(f">{path}:{s}-{e}\n")
        for i in range(0, len(seq), 80): fa.write(seq[i:i+80] + "\n")
        mb.write(f"{path}\t{s}\t{e}\t+\n")
os.remove(raw)
print(f"extract: {len(wanted)} members, {total} bp")
EOF
  st=$?
  echo "$st" > "$P/extract.exit"
  [ "$st" -eq 0 ] && touch "$P/extract.done" || { echo "=== extract FAILED"; exit 1; }
fi

# ---- stage syng: the locality panel (committed k63 s8 seed7 pos64)
stage syng $((64*1024*1024)) "$IMPG_SUB" syng -f "$P/locus.fa" --syncmer-length 63 \
  --smer-length 8 --syncmer-seed 7 --position-sample-rate 64 -o "$P/locus-pos64.syng" || exit 1

# ---- stage beds: the committed w10k d1k haplotype partition
BEDS="$P/partition-w10k-d1k/results"
if [ ! -f "$P/beds.done" ]; then
  mkdir -p "$BEDS"
  printf '%s\n' "$AXIS" > "$P/starting-sequences.txt"
  stage beds $((64*1024*1024)) "$IMPG_SUB" partition \
    -a "$P/locus-pos64.syng" -w 10000 -d 1000 \
    --starting-sequences-file "$P/starting-sequences.txt" \
    --selection-mode haplotype \
    --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5 \
    --no-rehome-singletons \
    --sequence-files "$P/locus.fa" \
    -o bed --separate-files --output-folder "$BEDS" \
    -t 8 || exit 1
fi
NP=$(ls "$BEDS" | grep -c '\.bed$'); echo "=== $NP partitions"

# ---- stage export + walks
IDS=$(seq 0 $((NP-1)))
mkdir -p "$P/export-graphs"
stage export $((64*1024*1024)) "$EXAMPLES/partition_graph_export" export "$P/locus-pos64.syng" \
  "$P/locus.fa" "$BEDS" "$P/export-graphs" $IDS || exit 1
if [ ! -f "$P/walks.done" ]; then
  budget_check
  python3 "$RV/rv-mk-walk-requests.py" "$BEDS" "$P/locus-pos64.syng.names" \
    "$P/walk-requests.tsv" || exit 1
  start=$(date +%s)
  nice -n 10 "$EXAMPLES/syng_homology_dump" walks "$P/locus-pos64.syng" "$P/walk-requests.tsv" \
    > "$P/export-walks.jsonl" 2> "$P/walks.err" \
    && { touch "$P/walks.done"; echo "$(( $(date +%s) - start ))" > "$P/walks.wall"; } \
    || { echo 1 > "$P/walks.exit"; exit 1; }
fi

# ---- stage maps: the locality graphs (K=63) + the reference axis + the
# domain structure (the KEPT member smalls are the free determinism gates)
if [ ! -f "$P/maps.done" ]; then
  budget_check
  start=$(date +%s)
  python3 /home/erikg/hprcv2/wgs35x-fleet/c-rider/locality-maps-export.py \
    "$AXIS" "$BEDS" "$P/export-graphs" "$P/export-walks.jsonl" \
    "$P/locality-graphs" "$P/reference-axis.json" "$P" \
    "$P/locality-domain-structure.txt" 63 \
    > "$P/maps.log" 2> "$P/maps.err" \
    && { touch "$P/maps.done"; echo "$(( $(date +%s) - start ))" > "$P/maps.wall"; } \
    || { echo "=== maps FAILED"; exit 1; }
  cmp -s "$P/reference-axis.json" "$MEMBERS/reference-axis.json" \
    && echo "reference-axis IDENTICAL to the kept member small" >> "$P/maps.log" \
    || echo "reference-axis DIFFERS from the kept member small" >> "$P/maps.log"
  cmp -s "$P/locality-domain-structure.txt" "$MEMBERS/locality-domain-structure.txt" \
    && echo "domain-structure IDENTICAL to the kept member small" >> "$P/maps.log" \
    || echo "domain-structure DIFFERS from the kept member small" >> "$P/maps.log"
fi
NW=$(ls "$P/locality-graphs" | grep -c 'window[0-9]*.gfa.map.json')
LOCI=$(seq -s, 0 $((NW-1)))
echo "=== $T: $NW windows, loci $LOCI"

# ---- stage reads: the BAM fetch (primaries in window +/-50kb, flat fastq)
if [ ! -f "$P/reads.done" ]; then
  budget_check
  start=$(date +%s)
  samtools view -u -F 0x904 "$BAM" "chr20:59950000-60300000" \
    | samtools fastq - 2> "$P/reads-fq.err" \
    | gzip -1 > "$WORK/reads.fq.gz"
  st=${PIPESTATUS[0]}
  echo "$st" > "$P/reads.exit"; echo "$(( $(date +%s) - start ))" > "$P/reads.wall"
  [ "$st" -eq 0 ] && touch "$P/reads.done" || { echo "=== reads FAILED ($st)"; exit 1; }
fi

# ---- stage rvcut: the restriction basis (the committed kept-read set);
# the census/cache are reaped in-pass (only the fasta + stats ride forward)
if [ ! -f "$P/rvcut.done" ]; then
  budget_check
  mkdir -p "$WORK/rvcut/$T"
  printf '{"name": "%s", "panel": "%s/locus-pos64.syng", "axis": "%s/reference-axis.json", "bed_directory": "%s/partition-w10k-d1k/results", "component": "%s", "sources": "%s/locus.fa", "out_dir": "%s/rvcut/%s"}\n' \
    "$T" "$P" "$P" "$P" "$AXIS" "$P" "$WORK" "$T" > "$WORK/rvcut-manifest.json"
  start=$(date +%s)
  RAYON_NUM_THREADS=32 IMPG_CUT_LOCAL_DERIVE=1 \
    nice -n 10 stdbuf -oL -eL "$EXAMPLES/sample_projection" cut-batch \
    --global-panel "$GSYNG" --agc "$AGC" --shards-dir "$SHARDS" \
    --manifest "$WORK/rvcut-manifest.json" > "$P/rvcut.log" 2> "$P/rvcut.err" &
  pid=$!; guard_poll "$pid" "$P/rvcut.rss" $((64*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/rvcut.exit"; echo "$(( $(date +%s) - start ))" > "$P/rvcut.wall"
  if [ "$st" -eq 0 ] && [ -f "$WORK/rvcut/$T/cut-stats.json" ]; then
    touch "$P/rvcut.done"
    rm -f "$WORK/rvcut/$T/cut-census.jsonl" "$WORK/rvcut/$T/cut-cache.bin"
  else echo "=== rvcut FAILED ($st)"; exit 1; fi
fi

# ---- stage cutA: the restricted global-frame census (the committed cut)
if [ ! -f "$P/cutA.done" ]; then
  budget_check
  start=$(date +%s)
  IMPG_CUT_READSET_FASTA="$WORK/rvcut/$T/cut-reads.fasta" \
  RAYON_NUM_THREADS=32 nice -n 10 stdbuf -oL -eL "$EXAMPLES/sample_projection" cut --frame global \
    --record-origin 1 \
    --global-panel "$GSYNG" \
    --axis "$P/reference-axis.json" --bed-directory "$BEDS" --component "$AXIS" \
    --agc "$AGC" --shards-dir "$SHARDS" \
    --out-dir "$WORK/cutA" \
    > "$P/cutA.log" 2> "$P/cutA.err" &
  pid=$!; guard_poll "$pid" "$P/cutA.rss" $((120*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/cutA.exit"; echo "$(( $(date +%s) - start ))" > "$P/cutA.wall"
  { [ "$st" -eq 0 ] && [ -f "$WORK/cutA/cut-stats.json" ]; } && touch "$P/cutA.done" \
    || { echo "=== cutA FAILED ($st)"; exit 1; }
  # the restriction basis' consumer is done: reap the fasta in-pass
  rm -f "$WORK/rvcut/$T/cut-reads.fasta"
fi

# ---- stage routes: the whole-panel routes inventory (the re-derivation the
# reap receipt names at ~25min once-per-panel)
if [ ! -f "$P/routes.done" ]; then
  budget_check
  rm -rf "$WORK/routes"
  IMPG_PANEL_LEGACY_V1_OK=1 RAYON_NUM_THREADS=32 stage routes $((700*1024*1024)) "$IMPG_SUB" genome-infer derive-routes \
    --panel "$GSYNG" --sources "$AGC" --read-lengths 148 \
    --census "$WORK/cutA/cut-census.jsonl" --out-dir "$WORK/routes" || exit 1
fi

# ---- the calls: the committed callC assessment form VERBATIM (merge ON,
# serial RAYON=1 = core-seconds, --far-bin-prune, truth-qv + allwave), the
# phase-census + fold-set-digest + scoring-census instruments ON.
run_call() { local tag=$1 bin=$2 loci=$3 merge_env=${4:-1} pre_env=${5:-auto}
  if [ -f "$P/call-$tag.done" ]; then echo "=== call-$tag done (marker)"; return 0; fi
  budget_check
  rm -rf "$WORK/call-$tag"
  local start pid st
  start=$(date +%s)
  local merge_env_args=()
  [ "$merge_env" = "1" ] && merge_env_args+=("IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1")
  local pre_env_args=()
  [ "$pre_env" = "off" ] && pre_env_args+=("IMPG_PANEL_PRECOMPUTE=off")
  [ "$pre_env" = "on" ] && pre_env_args+=("IMPG_PANEL_PRECOMPUTE=on")
  env RAYON_NUM_THREADS=1 IMPG_PANEL_LEGACY_V1_OK=1 \
    "${merge_env_args[@]}" "${pre_env_args[@]}" \
    IMPG_TRUTH_QV_ALLWAVE=1 IMPG_TRUTH_QV_ALLWAVE_BIN="$QVBIWFA" \
    IMPG_LOCUS_SCORING_CENSUS=1 IMPG_LOCUS_PHASE_CENSUS=1 IMPG_FOLD_SET_DIGEST=1 \
    nice -n 10 stdbuf -oL -eL "$bin" genome-infer call-diplotypes \
    --panel "$GSYNG" \
    --routes "$WORK/routes" \
    --partition-graphs "$P/locality-graphs" \
    --component "$AXIS" \
    --census "$WORK/cutA/cut-census.jsonl" \
    --reads "$WORK/reads.fq.gz" \
    --derive-cache "$WORK/cutA/cut-cache.bin" \
    --loci "$loci" \
    --truth-qv "$AXIS" "$TRUTH2" \
    --far-bin-prune \
    --qv-biwfa "$QVBIWFA" \
    --global-frame --record-origin 1 \
    --out-dir "$WORK/call-$tag" \
    > "$P/call-$tag.log" 2> "$P/call-$tag.err" &
  pid=$!; guard_poll "$pid" "$P/call-$tag.rss" $((300*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/call-$tag.exit"; echo "$(( $(date +%s) - start ))" > "$P/call-$tag.wall"
  # THE WRITE-TIME RECEIPTS (counters + digests, never raw dumps)
  {
    echo "# call-$tag receipt ($(date -u +%FT%TZ)); binary $bin; loci $loci"
    echo "# wall_s $(cat "$P/call-$tag.wall") exit $(cat "$P/call-$tag.exit")"
    if [ -f "$WORK/call-$tag/calls.jsonl" ]; then
      echo "calls_sha256 $(sha256sum "$WORK/call-$tag/calls.jsonl" | cut -d' ' -f1) $(stat -c%s "$WORK/call-$tag/calls.jsonl") bytes"
    fi
    if [ -f "$WORK/call-$tag/calls.jsonl.truth-qv.jsonl" ]; then
      echo "truthqv_sha256 $(sha256sum "$WORK/call-$tag/calls.jsonl.truth-qv.jsonl" | cut -d' ' -f1) $(stat -c%s "$WORK/call-$tag/calls.jsonl.truth-qv.jsonl") bytes"
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
    grep -h "phase summary" "$P/call-$tag.err" || true
    python3 - "$WORK/call-$tag" <<'EOF'
import json, sys, os
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
                FNV_OFFSET = 0xCBF29CE484222325
                FNV_PRIME = 0x100000001B3
                d0 = FNV_OFFSET
                for fold in fi:
                    for byte in fold["length"].to_bytes(8, "little"): d0 ^= byte; d0 = (d0 * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
                    for m in fold["members"]:
                        for byte in m["path_name"].encode(): d0 ^= byte; d0 = (d0 * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
                        for byte in m["start"].to_bytes(8, "little"): d0 ^= byte; d0 = (d0 * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
                        for byte in m["end"].to_bytes(8, "little"): d0 ^= byte; d0 = (d0 * FNV_PRIME) & 0xFFFFFFFFFFFFFFFF
                print(f"fold_set_digest locus {r['locus']} folds {len(fi)} {d0:016x} (offline-mirror)")
        pc = r.get("phase_census")
        if pc: print(f"phase_census locus {r['locus']} " + " ".join(
            f"{k}={v}" for k, v in sorted(pc.items()) if k != "note"))
        w = r.get("walls")
        if w: print(f"walls locus {r['locus']} " + " ".join(f"{k}={v}" for k, v in sorted(w.items())))
EOF
  } > "$WORK/call-$tag.receipt.txt" 2>&1
  cp -a "$WORK/call-$tag.receipt.txt" "$KEEP/"
  [ "$st" -eq 0 ] && [ -s "$WORK/call-$tag/calls.jsonl" ] && touch "$P/call-$tag.done" \
    && return 0
  echo "=== call-$tag FAILED ($st)"; return 1
}

# ---- THE PRECOMPUTE STEP (the once-per-panel price at the GLOBAL
# frame, receipted): the row context + the fold layouts, keyed on
# content into the locality graphs' own panel-precompute-loci section.
if [ ! -f "$P/precompute.done" ]; then
  budget_check
  start=$(date +%s)
  env IMPG_PANEL_LEGACY_V1_OK=1 RAYON_NUM_THREADS=1 \
    nice -n 10 stdbuf -oL -eL "$IMPG_PRE" genome-infer precompute-panel-loci \
    --panel "$GSYNG" --routes "$WORK/routes" \
    --partition-graphs "$P/locality-graphs" \
    --component "$AXIS" \
    --global-frame --record-origin 1 \
    > "$P/precompute.log" 2> "$P/precompute.err" &
  pid=$!; guard_poll "$pid" "$P/precompute.rss" $((700*1024*1024)); poll=$!
  st=0; wait "$pid" || st=$?; kill "$poll" 2>/dev/null; wait "$poll" 2>/dev/null
  echo "$st" > "$P/precompute.exit"; echo "$(( $(date +%s) - start ))" > "$P/precompute.wall"
  {
    echo "# precompute receipt ($(date -u +%FT%TZ)); binary $IMPG_PRE; component $AXIS; frame global origin 1"
    echo "# wall_s $(cat "$P/precompute.wall") exit $(cat "$P/precompute.exit")"
    cat "$P/precompute.err"
    echo "asset_bytes $(du -sb "$P/locality-graphs/panel-precompute-loci" 2>/dev/null | cut -f1)"
  } > "$WORK/precompute.receipt.txt" 2>&1
  cp -a "$WORK/precompute.receipt.txt" "$KEEP/"
  [ "$st" -eq 0 ] && touch "$P/precompute.done" || { echo "=== precompute FAILED ($st)"; exit 1; }
fi

# THE PRECOMPUTE INCREMENT'S LEGS: the LIVE path (the toggle-restore
# gate + the before-timing) and the SERVE (the after-timing), BOTH
# against the committed receipts' answer digests + fold digests.
run_call slice2-off-all "$IMPG_PRE" "$LOCI" 1 off || exit 1
run_call slice2-on-all "$IMPG_PRE" "$LOCI" 1 on || exit 1

# ---- THE GATES: answers + fold-set digests vs the committed receipts
{
  committed_calls=$(grep -m1 '^calls_sha256' "$COMMITTED_RECEIPT" | awk '{print $2}')
  committed_truth=$(grep -m1 '^truthqv_sha256' "$COMMITTED_RECEIPT" | awk '{print $2}')
  pass=1
  for tag in slice2-off-all slice2-on-all; do
    c=$(grep -m1 '^calls_sha256' "$WORK/call-$tag.receipt.txt" | awk '{print $2}')
    t=$(grep -m1 '^truthqv_sha256' "$WORK/call-$tag.receipt.txt" | awk '{print $2}')
    [ "$c" = "$committed_calls" ] && echo "GATE calls $tag PASS ($c)" || { echo "GATE calls $tag FAIL ($c != $committed_calls)"; pass=0; }
    [ "$t" = "$committed_truth" ] && echo "GATE truth-qv $tag PASS ($t)" || { echo "GATE truth-qv $tag FAIL ($t != $committed_truth)"; pass=0; }
  done
  for tag in slice2-off-all slice2-on-all; do
    for locus_digest in $(grep '^fold_set_digest' "$COMMITTED_RECEIPT" | awk '{print $3":"$6}'); do
      locus=${locus_digest%%:*}; digest=${locus_digest##*:}
      got=$(grep -m1 "fold_set_digest locus $locus " "$WORK/call-$tag.receipt.txt" | awk '{print $6}')
      [ "$got" = "$digest" ] || { echo "GATE fold digest $tag locus $locus FAIL ($got != $digest)"; pass=0; }
    done
    n_committed=$(grep -c '^fold_set_digest' "$COMMITTED_RECEIPT")
    n_got=$(grep -c '^fold_set_digest' "$WORK/call-$tag.receipt.txt")
    [ "$n_got" = "$n_committed" ] && echo "GATE fold digests $tag PASS ($n_got/$n_committed identical)" || { echo "GATE fold digests $tag FAIL ($n_got vs $n_committed)"; pass=0; }
  done
  [ "$pass" = 1 ] && echo "=== T09724 SLICE2 GATES: ALL PASS" || echo "=== T09724 SLICE2 GATES: FAILURES (above)"
} | tee "$KEEP/gates-slice2-T09724.txt"
echo "=== T09724 SLICE2 GATES COMPLETE - receipts in $KEEP; reap with: $0 reap"
