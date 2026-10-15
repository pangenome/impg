#!/usr/bin/env python3
"""THE INCREMENTAL SPLIT (the pipeline overlap lever): the query stage's
per-window completion lines (query/shard-*.err: 'query complete: <TID>') are
the boundary signal; the out beds are line-buffered, so a completed window's
rows are durable. Every invocation: find newly-completed TIDs, extract their
rows from the live shard outputs, write prep/<TID>/{locus.refined.bed,
locus.bed, axis-path.txt | split.rejected}, and APPEND to prep/accepted.txt
(one line per TID, idempotent via the per-TID split.done marker). The final
full split (rv-split-refined-beds.py, over the completed shards) remains the
closing reconciliation.
usage: rv-split-incremental.py <RV>
"""
import json, os, re, sys

RV = sys.argv[1]
Q = f"{RV}/query"
PREP = f"{RV}/prep"
os.makedirs(PREP, exist_ok=True)

BRIDGE = f"{RV}/../grid/bridge-accepted.bed"
bridge = {}
for line in open(BRIDGE):
    a = line.rstrip("\n").split("\t")
    if len(a) >= 7:
        bridge[a[0]] = a

# 1. the completed set from the per-window log lines
completed = set()
for err in sorted(os.listdir(Q)):
    if not (err.startswith("shard-") and err.endswith(".err")):
        continue
    for line in open(f"{Q}/{err}", errors="replace"):
        m = re.search(r"query complete: (\S+?) \(", line)
        if m and m.group(1) in bridge:
            completed.add(m.group(1))

# 2. the already-split set
done = set()
donef = f"{PREP}/split-done.txt"
if os.path.exists(donef):
    done = {l.strip() for l in open(donef) if l.strip()}
new = sorted(completed - done)
if not new:
    print(f"incremental split: no new completed windows (completed {len(completed)}, split {len(done)})")
    sys.exit(0)
want = set(new)

# 3. one pass over the live shard outputs for the new windows' rows
rows = {t: [] for t in new}
for err in sorted(os.listdir(Q)):
    if not (err.startswith("shard-") and err.endswith(".out.bed")):
        continue
    p = f"{Q}/{err}"
    if not os.path.exists(p):
        continue
    with open(p, errors="replace") as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 6 and a[3] in want:
                rows[a[3]].append((a[0], int(a[1]), int(a[2]), a[5]))

# 4. write per territory + the acceptance ledger (append-only)
accepted_new, rejected_new = [], []
with open(donef, "a") as df:
    for tid in new:
        b = bridge[tid]
        d = f"{PREP}/{tid}"
        os.makedirs(d, exist_ok=True)
        r = rows[tid]
        with open(f"{d}/locus.bed", "w") as f:
            f.write(f"{b[4]}\t{b[5]}\t{b[6]}\t{tid}\n")
        if not r:
            with open(f"{d}/split.rejected", "w") as f:
                f.write("window-missing-from-query\n")
            rejected_new.append((tid, "window-missing-from-query"))
        else:
            # THE AXIS ROW RULE (deterministic, parameter-free): the axis is
            # the HG002#1 row on the SAME-NAMED contig as the bridge chrom,
            # within one territory-span slack of the bridge target (the
            # cross-chromosomal rows are segdup hits, never the axis). 0 or
            # >1 same-chrom candidates = the fail-closed classes.
            slack = int(b[6]) - int(b[5])
            axis = [x for x in r
                    if x[0] == f"HG002#1#{b[1]}"
                    and x[1] < int(b[6]) + slack and x[2] > int(b[5]) - slack]
            with open(f"{d}/locus.refined.bed", "w") as f:
                for pth, s, e, strand in r:
                    f.write(f"{pth}\t{s}\t{e}\t{strand}\n")
            with open(f"{d}/split.info.json", "w") as f:
                json.dump({"rows": len(r), "axis_rows": len(axis),
                           "hg002_2_rows": sum(1 for x in r if x[0].startswith("HG002#2#"))}, f)
            if len(axis) == 1:
                with open(f"{d}/axis-path.txt", "w") as f:
                    f.write(f"{axis[0][0]}:{axis[0][1]}-{axis[0][2]}\n")
                accepted_new.append(tid)
            else:
                cls = ("no-same-chrom-axis-row" if len(axis) == 0
                       else "multi-same-chrom-axis-row")
                with open(f"{d}/split.rejected", "w") as f:
                    f.write(cls + "\n")
                rejected_new.append((tid, cls))
        df.write(tid + "\n")
if accepted_new:
    with open(f"{PREP}/accepted.txt", "a") as f:
        for t in accepted_new:
            f.write(t + "\n")
if rejected_new:
    with open(f"{PREP}/rejected-incremental.tsv", "a") as f:
        for t, cls in rejected_new:
            f.write(f"{t}\t{cls}\n")
print(f"incremental split: completed {len(completed)} newly split {len(new)} "
      f"accepted +{len(accepted_new)} rejected +{len(rejected_new)}")
