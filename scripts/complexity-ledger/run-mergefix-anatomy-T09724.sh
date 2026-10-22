#!/usr/bin/env bash
# THE MERGE-FIX ANATOMY RUNNER (the fold-merge OPEN term's own instrument run,
# under THE STREAM DISCIPLINE): ONE tranche-2 S_EXP territory (T09724, the
# T03226-class shape — far-bin share 43.46% vs T03226's 43.4%, the fold-heavy
# l9 at 809.8s serial in the committed tranche-2 receipt) REBUILT transiently
# in a bounded scratch dir from the KEPT definitional smalls
# (/home/erikg/hprcv2/panel/tranche2/members/T09724/) + the surviving panel
# assets (the global syng/AGC, the HG002 BAM, the projection shards), the
# committed tranche recipe VERBATIM in stage shape (extract -> syng ->
# partition w10k d1k -> export -> walks -> maps -> reads -> rvcut -> cutA ->
# derive-routes -> the callC assessment form), marker-gated and resumable.
#
# THE DISCIPLINE'S OWN RULES, ENFORCED HERE:
#   - HARD WORKDIR BUDGET 150G (WORKDIR_BUDGET_GB): checked before every
#     stage; a stage that would cross it stops and says so (fail-closed).
#   - RECEIPTS AT WRITE TIME: every stage's .wall/.exit/.rss receipts + the
#     calls' phase-census/digest/answer-digest receipt lines are computed as
#     they run; no raw dump is banked for later analysis.
#   - TRANSIENT FASTAS REAPED IN-PASS: extract.raw.fa (the extraction), the
#     rvcut census/cache (the restriction cut) — deleted the moment their
#     consumer completes, inside the stage that produced them.
#   - RAW INTERMEDIATES BY NAMED EXCEPTION (rule 6): the scratch's
#     cut-census.jsonl + cut-cache.bin + reads.fq.gz + substrate (consumer:
#     the anatomy/gate calls below; budget: 150G; REAP TIME: the runner's
#     terminal state — the `reap` subcommand, run by the worker at its
#     checkpoint, success OR failure). Only the receipts + the product
#     answers (<out>/calls.jsonl, the truth-qv sidecar — KB-scale) survive.
# usage: run-anatomy.sh [reap]   (no arg: run/resume every pending stage)
set -uo pipefail
T=T09724
HPRC=/home/erikg/hprcv2
MEMBERS=$HPRC/panel/tranche2/members/$T
RV=/home/erikg/hprcv2/wgs35x-fleet/runtime-verification
WORK=/mnt/nvme2n1/erikg/mergefix-$T
KEEP=/mnt/nvme2n1/erikg/mergefix-$T-keep
WORKDIR_BUDGET_GB=150
IMPG_NEW=/home/erikg/impg-genome-inference/target-frontier/release/impg-instr-prefix-a67a
# the fix binary: set IMPG_FIX when the fix build exists (the fix legs are
# the after-gates; they run only with it)
IMPG_FIX=${IMPG_FIX:-}
IMPG_BASE=/home/erikg/impg-genome-inference/target-frontier/release/impg-baseline-9204f9e
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
  IMPG_PANEL_LEGACY_V1_OK=1 RAYON_NUM_THREADS=32 stage routes $((700*1024*1024)) "$IMPG_NEW" genome-infer derive-routes \
    --panel "$GSYNG" --sources "$AGC" --read-lengths 148 \
    --census "$WORK/cutA/cut-census.jsonl" --out-dir "$WORK/routes" || exit 1
fi

# ---- the calls: the committed callC assessment form VERBATIM (merge ON,
# serial RAYON=1 = core-seconds, --far-bin-prune, truth-qv + allwave), the
# phase-census + fold-set-digest + scoring-census instruments ON.
run_call() { local tag=$1 bin=$2 loci=$3 merge_env=${4:-1}
  if [ -f "$P/call-$tag.done" ]; then echo "=== call-$tag done (marker)"; return 0; fi
  budget_check
  rm -rf "$WORK/call-$tag"
  local start pid st
  start=$(date +%s)
  local merge_env_args=()
  [ "$merge_env" = "1" ] && merge_env_args+=("IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1")
  env RAYON_NUM_THREADS=1 IMPG_PANEL_LEGACY_V1_OK=1 \
    "${merge_env_args[@]}" \
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

# the legs, in order: the fold-heavy locus 9 alone first (the shape
# receipt) at the baseline binary AND the instrumented binary (the
# instrumentation's own answer-digest gate: identical answers), then the
# whole territory at both (the corpus-share before-number + the per-locus
# phase-census anatomy)
run_call base-l9 "$IMPG_BASE" 9 || exit 1
run_call instr-l9 "$IMPG_NEW" 9 || exit 1
run_call base-all "$IMPG_BASE" "$LOCI" || exit 1
run_call instr-all "$IMPG_NEW" "$LOCI" || exit 1
# THE PRE-MERGE ANATOMY LEG (the ledger's own l15 form: no merge env) —
# the fold-heavy shape's uncollapsed tie sets are where the per-pair
# constant bites; the phase census decomposes it per locus
run_call instr-premerge-all "$IMPG_NEW" "$LOCI" 0 || exit 1
echo "=== BASELINE ANATOMY COMPLETE - receipts in $KEEP"
# THE FIX LEGS (the after-gates): the same territory, the same forms, the
# fix binary - answers byte-identical (digest gate), walls the after-numbers
if [ -z "$IMPG_FIX" ]; then
  echo "=== IMPG_FIX unset: the fix legs are SKIPPED (set IMPG_FIX=<the fix binary> and relaunch)"
  exit 0
fi
run_call fix-all "$IMPG_FIX" "$LOCI" || exit 1
run_call fix-premerge-all "$IMPG_FIX" "$LOCI" 0 || exit 1
echo "=== FIX LEGS COMPLETE - receipts in $KEEP; reap with: $0 reap"
