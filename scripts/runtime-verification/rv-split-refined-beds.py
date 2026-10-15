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
in the datasheet). Rejected classes named (the bridge rule, same as
rv-split-incremental.py): no-axis-bridge / ambiguous-axis-bridge /
no-image-row / multi-image-row / transitive-expansion /
window-missing-from-query. Requires axis-bridges.paf (the drift-correct
axis - the old fixed-slack rule is gone).
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

# THE AXIS BRIDGE (the same rule as rv-split-incremental.py - the final
# reconciliation split must not carry a different axis rule than the
# incremental one): the window's UNIQUE primary same-chrom asm5 row against
# HG002#1#chr1-22 (rv-axis-bridge.py, axis-bridges.paf); the axis ROW is then
# the refined bed's same-contig row with MAXIMUM OVERLAP vs the bridge
# interval. The committed HG002#1 assembly's coordinate offset DRIFTS along
# the chromosome (+113kb at 3.3M -> +3.8Mb past 150M on chr1) - the old
# fixed +/-250kb slack rule wrongly rejected every territory past the drift
# point; the bridge handles the drift, never a slack.
axis_bridge = {}
paf = f"{RV}/axis-bridges.paf"
if os.path.exists(paf):
    for line in open(paf):
        a = line.rstrip("\n").split("\t")
        if len(a) < 11:
            continue
        axis_bridge.setdefault(a[0], []).append((a[5], int(a[7]), int(a[8])))
else:
    sys.exit("axis-bridges.paf missing - the final split requires the bridge")

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
    ab = [iv for iv in axis_bridge.get(tid, []) if iv[0] == f"HG002#1#{b[1]}"]
    axis = []
    if len(ab) == 1:
        tS, tE = ab[0][1], ab[0][2]
        cand = [(min(x[2], tE) - max(x[1], tS), x) for x in rows
                if x[0] == f"HG002#1#{b[1]}"]
        cand = [c for c in cand if c[0] > 0]
        if cand:
            best = max(c[0] for c in cand)
            axis = [c[1] for c in cand if c[0] == best]
    hg2 = sum(1 for r in rows if r[0].startswith("HG002#2#"))
    stats["hg2_rows"] += hg2
    with open(f"{d}/locus.refined.bed", "w") as f:
        for p, s, e, strand in rows:
            f.write(f"{p}\t{s}\t{e}\t{strand}\n")
    with open(f"{d}/split.info.json", "w") as f:
        import json
        json.dump({"rows": len(rows), "axis_rows": len(axis),
                   "hg002_2_rows": hg2}, f)
    npaths = len({x[0] for x in rows})
    if npaths > 2000 or len(rows) > 2 * npaths:
        rejected.append((tid, f"transitive-expansion rows {len(rows)} paths {npaths}", len(rows), hg2))
        open(f"{d}/split.rejected", "w").write(f"transitive-expansion rows {len(rows)} paths {npaths}\n")
    elif len(axis) == 1:
        accepted.append(tid)
        with open(f"{d}/axis-path.txt", "w") as f:
            f.write(f"{axis[0][0]}:{axis[0][1]}-{axis[0][2]}\n")
    else:
        if len(ab) == 0:
            cls = "no-axis-bridge"
        elif len(ab) > 1:
            cls = "ambiguous-axis-bridge"
        elif len(axis) == 0:
            cls = "no-image-row"
        else:
            cls = "multi-image-row"
        rejected.append((tid, cls, len(rows), hg2))
        open(f"{d}/split.rejected", "w").write(cls + "\n")

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
