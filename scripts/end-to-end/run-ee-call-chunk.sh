#!/usr/bin/env bash
# THE END-TO-END RECEIPT, STEP 3 - THE CALL CHUNK (the fleet form,
# all committed defaults, ZERO product flags): one call-region
# invocation over a 32-territory chunk of the universe. The parent
# loads the ONE global panel and stays the load owner (fork-shared
# COW); the width-32 child builds its own AGC index, batch-cuts the
# chunk at the RESTRICTED census profile (the basis auto-derived from
# each territory's map-step cut-reads.fasta), and runs the chunk's 32
# calls concurrently at the T16 rung, with the FAR-BIN PRUNE default
# ON and the record origin auto-derived from the registry. The 64GiB
# guard over every worker (the child is the worker; a guard kill is a
# named fail-closed census entry - the leg re-runs on resume).
#
# THE CHUNK ROOT: the call step's own --map interface is per-invocation
# (one reads path, one root); the chunk root assembles the per-
# territory map legs' own cut outputs (symlinks - the batch-
# composition invariance receipt is the warrant: single == batch
# member) + the chunk's own reads fetch (the committed +-50kb rule
# per territory, -M multi-region: a SUPERSET per member, harmless to
# the census-id-keyed reads lookup) + the manifest in the map step's
# own schema. usage: run-ee-call-chunk.sh <chunk-id>
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh
CID=$1
LIST=$EE/chunks/chunk-$CID.list
L=$EE/ledger/chunk-$CID
RCP=$EE/receipts/call-chunk-$CID.txt
mkdir -p "$L"
[ -s "$LIST" ] || { echo "no chunk list $LIST" >&2; exit 1; }
mapfile -t TS < "$LIST"
[ -f "$L/call.done" ] && exit 0
# resumable: all members already terminal -> done
pending=()
for T in "${TS[@]}"; do
  [ -f "$EE/ledger/$T/call.done" ] || pending+=("$T")
done
[ "${#pending[@]}" -eq 0 ] && { touch "$L/call.done"; exit 0; }
[ "${#pending[@]}" -ne "${#TS[@]}" ] && \
  echo "$(ee_now) chunk $CID partial resume: ${#pending[@]} of ${#TS[@]} pending" >> "$RCP"
TS=("${pending[@]}")
ee_budget_check || { echo budget >> "$RCP"; exit 1; }

for T in "${TS[@]}"; do
  [ -f "$EE/ledger/$T/map.done" ] || { echo "$T map not done" >> "$RCP"; exit 1; }
done

# ---- the chunk reads fetch (the committed rule, multi-region, -M) ----
REGIONS=""
for T in "${TS[@]}"; do
  read -r CHROM START END < <(awk -v t="$T" '$4==t{print $1, $2, $3; exit}' "$RV/bed4-genome.tsv")
  REGIONS="$REGIONS $CHROM:$((START-50000))-$((END+50000))"
done
if [ ! -f "$L/reads.done" ]; then
  ee_step_out "$L" reads "$EE/reads/chunk-$CID.fq.gz" \
    bash -c 'samtools view -u -M -F 0x904 "$1" $2 | samtools fastq - 2>"$3" | gzip -1' \
    _ "$BAM" "$REGIONS" "$L/reads.samtools.err" || exit 1
  touch "$L/reads.done"
fi

# ---- the chunk map root (the call step's own interface) -----------
ROOT=$EE/scratch/callroot-$CID
rm -rf "$ROOT"; mkdir -p "$ROOT"
for T in "${TS[@]}"; do
  ln -s "$EE/map/$T/$T" "$ROOT/$T"
done
python3 - "$ROOT" "$EE/reads/chunk-$CID.fq.gz" "${TS[@]}" <<'PY'
import json, sys
root, reads, ts = sys.argv[1], sys.argv[2], sys.argv[3:]
json.dump({"version": 1,
           "panel": "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng",
           "shards_dir": "/home/erikg/hprcv2/hg002-projection/full",
           "reads": reads,
           "partitions": list(ts)}, open(f"{root}/manifest.json", "w"))
PY

# ---- THE CALL (zero product flags; the fleet form's own defaults) --
: > "$RCP"
echo "$(ee_now) chunk $CID: ${#TS[@]} territories, the fleet form (COW parent load + width-32 child, T16, restricted, farbin-prune default)" >> "$RCP"
OUT=$EE/calls/chunk-$CID
rm -rf "$OUT"
t0=$(date +%s)
cpu=$L/call.cpu
# resumable attempt naming: the call step's out-dir discipline refuses
# an existing dir, so partial-resume attempts take fresh suffixes; the
# per-territory markers stay authoritative
ATT=$(ls -d "$EE/calls"/chunk-$CID.a* 2>/dev/null | wc -l)
OUT=$EE/calls/chunk-$CID.a$ATT
( TIMEFORMAT='%U %S'; time (
    "$IMPG_NEW" genome-infer call-region \
      --map "$ROOT" --partitions "$EE/partitions" \
      --region "$(IFS=,; echo "${TS[*]}")" \
      --out-dir "$OUT"
  ) ) > "$L/call.log" 2> "$cpu" &
cpid=$!
gp=$(ee_guard_tree "$cpid" "$GUARD_GIB" "$RCP" "chunk-$CID-worker")
tp=$(ee_rss_track "$cpid" "$RCP" "chunk-$CID-parent")
wait "$cpid"; st=$?
kill "$gp" "$tp" 2>/dev/null; wait "$gp" "$tp" 2>/dev/null
wall=$(( $(date +%s) - t0 ))
echo "$(ee_now) chunk $CID call-region exit=$st wall=${wall}s cpu='$(cat $cpu 2>/dev/null)'" >> "$RCP"
echo "guard events: $(grep -c GUARD "$RCP" || true)" >> "$RCP"
[ "$st" -ne 0 ] && { echo "$st" > "$L/call.exit"; exit "$st"; }

# ---- per-territory markers + receipts at write time ---------------
OUT_SAVE=$OUT
gzip -9 -c "$L/call.log" > "$EE/receipts/call-chunk-$CID.log.gz"
for T in "${TS[@]}"; do
  TD=$OUT/$T
  [ -f "$TD/calls.jsonl" ] || { echo "$(ee_now) chunk $CID: $T has no calls.jsonl" >> "$RCP"; continue; }
  touch "$EE/ledger/$T/call.done"
  # keep the receipts smalls beside the answers
  for f in cut-stats.json; do
    [ -f "$OUT/$T-cut/$f" ] && cp "$OUT/$T-cut/$f" "$EE/ledger/$T/restricted-cut-stats.json"
    [ -f "$EE/map/$T/$T/$f" ] && cp "$EE/map/$T/$T/$f" "$EE/ledger/$T/map-cut-stats.json"
  done
  # THE INSTRUMENT-RECEIPT MINE (the stream law: counters + walls
  # computed at write time; the raw per-record/per-ingredient streams
  # are the reaped class - the tranche-2 lesson, ~600MB/territory at
  # L1-class vs a 35KB calls.jsonl)
  python3 - "$TD" "$EE/ledger/$T" "$CID" >> "$RCP" 2>&1 <<'PY'
import json, os, re, sys
td, led, cid = sys.argv[1], sys.argv[2], sys.argv[3]
keep = re.compile(r"count|seconds|_s$|mult|pairs|prune|served|fallback|records|loci|locus$|folds|classes|memo|merge|qual_entries|scan")
def walk(o, path, out):
    if isinstance(o, dict):
        for k, v in o.items(): walk(v, f"{path}.{k}", out)
    elif isinstance(o, list):
        out[f"{path}#len"] = len(o)
    elif isinstance(o, (int, float, bool)) or o is None:
        if keep.search(path): out[path] = o
mined = {"chunk": cid,
         "calls_bytes": os.path.getsize(f"{td}/calls.jsonl"),
         "loci": sum(1 for _ in open(f"{td}/calls.jsonl"))}
ir = f"{td}/instrument-receipt.jsonl"
if os.path.exists(ir):
    per_locus = []
    for line in open(ir):
        try: rec = json.loads(line)
        except Exception: continue
        small = {}; walk(rec, "", small)
        per_locus.append(small)
    mined["instrument"] = per_locus
open(f"{led}/leg.json", "a").write(json.dumps(mined) + "\n")
print(f"   {os.path.basename(led)}: calls {mined['calls_bytes']}B loci {mined['loci']} "
      f"instrument mined {len(mined.get('instrument', []))} locus lines")
PY
  # the raw instrument streams are the reaped class (the answers +
  # the mined counters + the chunk's gzipped log are the receipts)
  rm -f "$TD"/instrument-receipt.* 2>/dev/null
  cp "$TD/manifest.json" "$EE/ledger/$T/call-manifest.json" 2>/dev/null || true
done
OUT=$OUT_SAVE

# ---- THE JIT REAPER (the stream law: transients reaped in-pass at
# the leg's completion; the answers + receipts are the keep set) ----
freed=$(du -sBG "$ROOT" "$EE/reads/chunk-$CID.fq.gz" 2>/dev/null | awk '{s+=$1} END{print s+0}')
rm -rf "$ROOT" "$EE/reads/chunk-$CID.fq.gz"
for T in "${TS[@]}"; do
  [ -f "$EE/ledger/$T/call.done" ] || continue
  p=$EE/prep/$T
  # the substrate + the precompute asset (the named exception's own
  # reap condition: the territory's call is complete)
  rm -rf "$p/locus.fa" "$p/extract.raw.fa" "$p"/locus-pos64.syng* \
        "$p/partition-w10k-d1k" "$p/export-graphs" "$p/export-walks.jsonl" \
        "$p/walk-requests.tsv" "$p/locality-graphs" 2>/dev/null
  rm -rf "$EE/map/$T" "$EE/reads/$T.fq.gz" "$OUT/$T-cut" 2>/dev/null
done
touch "$L/call.done"
echo "$(ee_now) chunk $CID COMPLETE: ${#TS[@]} territories, wall ${wall}s; reaped ~${freed}G of transients (answers + receipts kept)" >> "$RCP"
echo "chunk $CID done wall=${wall}s"
