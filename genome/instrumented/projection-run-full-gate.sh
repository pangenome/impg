#!/usr/bin/env bash
# THE FULL-SAMPLE IDENTITY GATE (the census source is the artifact):
# (1) BYTE-COMPARE the dedup'd occurrence stream: run_dedup emits each
#     hash partition SORTED (BTreeMap), so occurrences.bin is a
#     canonical function of the pattern multiset alone - byte identity is
#     the strongest form (same records, same occurrences).
# (2) per-shard COMMUTATIVE PACKET-MULTISET HASH over BOTH full shard
#     sets (shard-hash.py): the 100% identity standard. Index-window
#     sampling is NOT valid full-vs-full: the packet order within a
#     shard is the send-in-par scheduling's own nondeterminism layered
#     over arrival-batch stream interleaving, so equal-index windows
#     hold different batch memberships (measured: reads/bases equal,
#     records +-17 in 409,600 - the batch-edge jitter class, receipts in
#     gates/full.windows.*.txt); the whole-file multiset cancels it.
# (3) projection-stats counters: reads/bases/records/read-lengths equal.
set -uo pipefail
B=/home/erikg/hprcv2/proj-bandwidth
FULL=/home/erikg/hprcv2/hg002-projection/full
OPT="$B/full-optimized"
G="$B/gates"; mkdir -p "$G"
F="$G/full.committed.txt"
: > "$F"
ok=0
echo "== (1) occurrences.bin byte-compare (canonical: partitions sorted)" >> "$F"
if cmp "$FULL/occurrences.bin" "$OPT/occurrences.bin" >> "$F" 2>&1; then
  echo "   occurrences.bin BYTE-IDENTICAL" >> "$F"
else
  echo "   occurrences.bin DIFFERS" >> "$F"; ok=1
fi
echo "== (3) projection-stats counters" >> "$F"
python3 - "$FULL/projection-stats.json" "$OPT/projection-stats.json" >> "$F" 2>&1 <<'PY'
import json, sys
a = json.load(open(sys.argv[1])); b = json.load(open(sys.argv[2]))
for k in ("reads", "bases", "reads_with_mems", "records", "read_lengths", "per_shard_reads"):
    if a[k] != b[k]:
        print(f"   MISMATCH {k}: {a[k]} vs {b[k]}"); sys.exit(1)
print(f"   counters identical (reads {a['reads']}, records {a['records']})")
PY
[ $? -ne 0 ] && ok=1
echo "== (2) per-shard commutative packet-multiset hash (the 100% standard)" >> "$F"
ls "$OPT"/shard.*.bin.gz "$FULL"/shard.*.bin.gz \
  | xargs -P 24 -n 1 nice -n 19 python3 "$B/shard-hash.py" \
  > "$G/full-shard-hashes.gate.txt" 2>&1
python3 - "$G/full-shard-hashes.gate.txt" >> "$F" 2>&1 <<'PY'
import re, sys
from collections import defaultdict
d = defaultdict(dict)
for line in open(sys.argv[1]):
    t = line.split()
    m = re.search(r'shard\.(fifo_\d+)\.bin\.gz', t[0])
    side = 'optimized' if t[0].startswith('/home/erikg/hprcv2/proj-bandwidth') else 'committed'
    d[m.group(1)][side] = (int(t[2]), int(t[4]), t[-1])
bad = [s for s, v in d.items() if len(v) != 2 or v['committed'] != v['optimized']]
print(f"   shards {len(d)}: " + ("ALL IDENTICAL (packets + records + 128-bit multiset sum per shard)"
      if not bad and len(d) == 19 else f"MISMATCH {bad}"))
sys.exit(1 if bad or len(d) != 19 else 0)
PY
[ $? -ne 0 ] && ok=1
[ "$ok" -eq 0 ] && echo "FULL-SAMPLE GATE PASS" >> "$F" || echo "FULL-SAMPLE GATE FAIL" >> "$F"
echo "full gate exit $ok"; exit $ok
