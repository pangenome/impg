#!/usr/bin/env python3
"""THE WHOLE-GENOME RUNTIME VERIFICATION - the refined-bed split (the old
query/split-refined-beds.py was disk-cleaned with the query dir; this is the
re-implementation from the observed committed conventions): the query stage's
shard output beds (rows: path, start, end, TID, ., strand) split per window
into prep/<TID>/locus.refined.bed (path, start, end, strand) + the territory's
locus.bed row (the CHM13#0 bridge target, from grid/bridge-accepted.bed), and
the ACCEPTANCE LEDGER: a territory is callable iff its refined bed carries
EXACTLY ONE HG002#1 axis row (the axis component path = that row's
path:start-end). The runtime receipt needs no truth; HG002#2 rows are recorded
as info (the accuracy-grade set additionally required them - the delta named
in the datasheet). Rejected classes named: no-axis-row / multi-axis-row /
window-missing-from-query.
usage: rv-split-refined-beds.py <RV> <grid-bridge-accepted.bed> <out.bed>...
"""
import os, sys

RV, BRIDGE = sys.argv[1:3]
outs = sys.argv[3:]
PREP = f"{RV}/prep"
os.makedirs(PREP, exist_ok=True)

# the bridge rows: TID, chrom, grch38S, grch38E, CHM13 target, tS, tE, alen, matches
bridge = {}
for line in open(BRIDGE):
    a = line.rstrip("\n").split("\t")
    if len(a) >= 7:
        bridge[a[0]] = a

seen = {}
for path in outs:
    for line in open(path):
        a = line.rstrip("\n").split("\t")
        if len(a) < 6 or not a[0]:
            continue
        seen.setdefault(a[3], []).append((a[0], int(a[1]), int(a[2]), a[5]))

accepted, rejected = [], []
stats = {"windows": 0, "rows": 0, "hg2_rows": 0}
for tid in sorted(bridge):
    stats["windows"] += 1
    rows = seen.get(tid)
    b = bridge[tid]
    d = f"{PREP}/{tid}"
    os.makedirs(d, exist_ok=True)
    with open(f"{d}/locus.bed", "w") as f:
        f.write(f"{b[4]}\t{b[5]}\t{b[6]}\t{tid}\n")
    if not rows:
        rejected.append((tid, "window-missing-from-query", 0, 0))
        open(f"{d}/split.rejected", "w").write("window-missing-from-query\n")
        continue
    stats["rows"] += len(rows)
    slack = int(b[6]) - int(b[5])
    axis = [r for r in rows
            if r[0] == f"HG002#1#{b[1]}"
            and r[1] < int(b[6]) + slack and r[2] > int(b[5]) - slack]
    hg2 = sum(1 for r in rows if r[0].startswith("HG002#2#"))
    stats["hg2_rows"] += hg2
    with open(f"{d}/locus.refined.bed", "w") as f:
        for p, s, e, strand in rows:
            f.write(f"{p}\t{s}\t{e}\t{strand}\n")
    with open(f"{d}/split.info.json", "w") as f:
        import json
        json.dump({"rows": len(rows), "axis_rows": len(axis),
                   "hg002_2_rows": hg2}, f)
    if len(axis) == 0:
        rejected.append((tid, "no-same-chrom-axis-row", len(rows), hg2))
        open(f"{d}/split.rejected", "w").write("no-same-chrom-axis-row\n")
    elif len(axis) > 1:
        rejected.append((tid, "multi-same-chrom-axis-row", len(rows), hg2))
        open(f"{d}/split.rejected", "w").write("multi-same-chrom-axis-row\n")
    else:
        accepted.append(tid)
        with open(f"{d}/axis-path.txt", "w") as f:
            f.write(f"{axis[0][0]}:{axis[0][1]}-{axis[0][2]}\n")

with open(f"{PREP}/accepted.txt", "w") as f:
    f.write("\n".join(accepted) + ("\n" if accepted else ""))
with open(f"{PREP}/rejected.tsv", "w") as f:
    for tid, cls, rows, hg2 in rejected:
        f.write(f"{tid}\t{cls}\t{rows}\t{hg2}\n")
import json
with open(f"{PREP}/split-stats.json", "w") as f:
    json.dump(stats, f, indent=1)
print(f"windows {stats['windows']} accepted {len(accepted)} rejected {len(rejected)} "
      f"rows {stats['rows']}")
