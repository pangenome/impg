#!/usr/bin/env bash
# rv-close-reap.sh - THE CAMPAIGN-CLOSE REAP (the owner's standing disk ruling:
# "never keep huge amounts of intermediate data"; the campaign STOPPED BY OWNER
# before all-terminal, so the built-ahead substrate and the cut evidence that
# the dead wave will never consume are regenerable bulk, not receipts).
# THE CLOSE CLASSES (both idempotent, both measured before/after):
#   CLOSE-A: the per-territory REGENERABLE substrate (locus.fa, the syng
#     slices, export-graphs, export-walks.jsonl, walk-requests.tsv,
#     locality-graphs - the JIT reaper proved regen ~22s/territory, the
#     T03194 full-chain re-entry receipt) on ACCEPTED territories that will
#     never be called (no call terminal) - plus the same files' residue on
#     terminal/rejected territories (idempotent residue sweep, expected ~0 -
#     the wave's own in-lane reap took them at call time).  Producer markers
#     (extract/syng/export/walks/walkreq/maps/substrate .done) DROP with the
#     reaped outputs, per the JIT reaper's idempotence doctrine.
#   CLOSE-B: the map/ CUT-EVIDENCE BULK (cut-cache.bin, cut-census.jsonl,
#     cut-reads.fasta) for territories without a live call lane - the cut is
#     deterministic and byte-identity-gated (regression-L1), the wave is dead,
#     the call will never read them.  map/<T>/<T>/cut-stats.json KEPT (the
#     cut receipt rv-rollup.py reads).
# NEVER TOUCHED: the input BAM, hg002-projection/full (outside this tree),
# every definitional bed/ledger (locus.refined.bed, locus.members.bed,
# locus.bed, split.*, extract.specs.txt, axis-path.txt, reference-axis.json,
# partition-w10k-d1k/, *.wall, *.err, *.log, markers of terminal steps),
# receipts/, calls/, cuts/, query/, partitions/, reads/, wave-logs/,
# poison.tsv, the handoff, the datasheet.
# LIVE-LANE EXCLUSIONS: never touch a territory with a live call lane, a
# live cut-batch manifest entry, or a live run-prep-one.
RV=${RV:-/home/erikg/hprcv2/wgs35x-fleet/runtime-verification}
LOG=$RV/close-reap.log
cd "$RV" || exit 1
log(){ echo "$(printf '%(%FT%TZ)T' -1) $*" >> "$LOG"; }
dfgb(){ df -BG --output=avail / | tail -1 | tr -dc '0-9'; }

log "=== campaign-close reap start (stopped-by-owner close; the honest tails resolved) ==="
DF0=$(dfgb); PREP0=$(du -sm --apparent-size prep 2>/dev/null | cut -f1); MAP0=$(du -sm --apparent-size map 2>/dev/null | cut -f1)
log "before: df ${DF0}G avail, prep ${PREP0}M apparent, map ${MAP0}M apparent"

# 1. live-lane exclusions (never touch a territory a live process works on)
EXCL=" "
for p in $(pgrep -f "one_territory" 2>/dev/null); do
  [ "$p" = "$$" ] && continue
  T=$(tr '\0' '\n' < /proc/$p/cmdline 2>/dev/null | tail -1)
  [ -n "$T" ] && EXCL="$EXCL$T "
done
for p in $(pgrep -f "sample_projection cut-batch" 2>/dev/null); do
  M=$(tr '\0' '\n' < /proc/$p/cmdline 2>/dev/null | awk '/manifest/{getline; print}')
  [ -n "$M" ] && while read -r t; do EXCL="$EXCL$t "; done < <(sed -n 's/.*"name": "\([^"]*\)".*/\1/p' "$M" 2>/dev/null)
done
for p in $(pgrep -f "run-prep-one" 2>/dev/null); do
  T=$(tr '\0' '\n' < /proc/$p/cmdline 2>/dev/null | tail -1)
  [ -n "$T" ] && EXCL="$EXCL$T "
done
log "live-lane exclusions:${EXCL}"
excl(){ case "$EXCL" in *" $1 "*) return 0;; *) return 1;; esac; }

# 2. the never-called set: accepted minus call-terminal (calls.jsonl or call.exit)
> /tmp/close-called.lst
ls -d receipts/*/ 2>/dev/null | sed 's|receipts/||; s|/||' >> /tmp/close-called.lst
ls -d prep/*/call.exit 2>/dev/null | sed 's|prep/||; s|/.*||' >> /tmp/close-called.lst
sort -u -o /tmp/close-called.lst /tmp/close-called.lst
NEVER=$(comm -23 <(sort prep/accepted.txt) /tmp/close-called.lst)
log "never-called accepted territories: $(wc -w <<< "$NEVER") (callok+callfail terminal: $(wc -l < /tmp/close-called.lst))"

# 3. CLOSE-A: the regenerable substrate, never-called first (the bulk), then
#    the idempotent residue sweep over every prep/T* dir.
REGEN="locus.fa locus.fa.fai export-walks.jsonl walk-requests.tsv export-graphs locality-graphs"
> /tmp/close-A.lst; AMK=""
for T in $NEVER; do
  excl "$T" && continue
  P=prep/$T; [ -d "$P" ] || continue
  for f in $REGEN; do [ -e "$P/$f" ] && echo "$P/$f" >> /tmp/close-A.lst; done
  for f in "$P"/locus-pos64.syng.*; do [ -e "$f" ] && echo "$f" >> /tmp/close-A.lst; done
  AMK="$AMK $T"
done
# residue sweep (idempotent): the same regenerable names on every other prep dir
for P in prep/T*; do
  T=$(basename "$P")
  case " $AMK " in *" $T "*) continue;; esac
  excl "$T" && continue
  for f in $REGEN; do [ -e "$P/$f" ] && echo "$P/$f" >> /tmp/close-A.lst; done
  for f in "$P"/locus-pos64.syng.*; do [ -e "$f" ] && echo "$f" >> /tmp/close-A.lst; done
done
mb_of(){ [ -s "$1" ] && xargs -r -n 256 du -sm --apparent-size < "$1" 2>/dev/null | awk '{s+=$1} END{print s+0}' || echo 0; }
A_MB=$(mb_of /tmp/close-A.lst); A_N=$(wc -l < /tmp/close-A.lst)
xargs -r -a /tmp/close-A.lst rm -rf --
for T in $AMK; do rm -f "prep/$T/extract.done" "prep/$T/syng.done" "prep/$T/export.done" "prep/$T/walks.done" "prep/$T/walkreq.done" "prep/$T/maps.done" "prep/$T/substrate.done"; done
log "CLOSE-A reaped: $A_N files/dirs, ${A_MB}M regenerable substrate (locus.fa/syng slices/export maps/walks; the ~22s regen receipt) on $(wc -w <<< "$AMK") never-called territories + residue sweep; producer markers dropped (idempotence); beds/ledgers/axis/partition beds KEPT"

# 4. CLOSE-B: the map/ cut-evidence bulk (cut-stats.json KEPT - the rollup's receipt)
> /tmp/close-B.lst; BMK=""
for D in map/T*/*; do
  [ -d "$D" ] || continue
  T=$(basename "$(dirname "$D")")
  excl "$T" && continue
  for f in cut-cache.bin cut-census.jsonl cut-reads.fasta; do [ -e "$D/$f" ] && echo "$D/$f" >> /tmp/close-B.lst; done
  BMK="$BMK $T"
done
B_MB=$(mb_of /tmp/close-B.lst); B_N=$(wc -l < /tmp/close-B.lst)
xargs -r -a /tmp/close-B.lst rm -rf --
log "CLOSE-B reaped: $B_N files, ${B_MB}M cut-evidence bulk (cache/census/reads - the deterministic byte-identity-gated cut, regenerable) on $(wc -w <<< "$BMK") territories; map/<T>/<T>/cut-stats.json KEPT (the cut receipt)"

DF1=$(dfgb); PREP1=$(du -sm --apparent-size prep 2>/dev/null | cut -f1); MAP1=$(du -sm --apparent-size map 2>/dev/null | cut -f1)
log "after: df ${DF1}G avail, prep ${PREP1}M apparent, map ${MAP1}M apparent - freed prep $((PREP0-PREP1))M + map $((MAP0-MAP1))M apparent (df ${DF0}G -> ${DF1}G)"
log "=== campaign-close reap done ==="
