#!/usr/bin/env python3
"""GATE 1, part 2: THE SEAM-EXTENSION CENSUS over the chrMT/chrI
localities - the measurement of the seam rejoin under the alignment-
induced locality domain.

THE SEAM REJOIN (exact, derived): a universe path's fold material
extends along the path's own partition-BED tiling (abutting/
overlapping rows - the panel's proven coordinate-continuous spine)
from its locality-component row, including every chained row that
carries at least one interned syng node SHARED with the window's
axis-row walk (the export graph's own P-line spelling, signed global
syncmer ids, gap splices excluded). The partition boundary severed
one continuous material; the shared-node coverage names exactly the
severed piece - the ortholog material the window's territory
projects onto that path. The chain stops at rows carrying no shared
node; a disjoint repeat copy is never chained in.

The walk nodes of a tiling row are testable only where that row's
partition has an export graph (the committed build-set closure,
1,414 partitions - truth-side closed by construction: every SK1
holder partition is built). Rows in unbuilt partitions are counted
honestly as UNTESTABLE and bound the chain.

Per window the census reports: the universe paths, how many extend,
the hull bp added, the truth-second (SK1) verdict (does the truth's
ortholog coverage need rows beyond the window's partition - the
severed-truth class), and the untestable-boundary count.

Receipt-side only. Usage: locality-seam-census.py chrI
"""
import json
import os
import sys
from collections import defaultdict

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
PGGB = f"{D}/pggb-partition-graphs-chrMT-chrI"
MAPS = f"{D}/pggb-projected-graphs-chrMT-chrI"
EXPORT = f"{D}/partition-graphs"
BEDDIR = "/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results"
NUM_SYNCMER_NODES = 10024605
TRUTH_SECOND = "SK1#0#"


class UnionFind:
    def __init__(self):
        self.p = {}

    def find(self, x):
        self.p.setdefault(x, x)
        while self.p[x] != x:
            self.p[x] = self.p[self.p[x]]
            x = self.p[x]
        return x

    def union(self, a, b):
        ra, rb = self.find(a), self.find(b)
        if ra != rb:
            self.p[ra] = rb


def read_gfa_paths(path):
    """P lines only: row name -> interned syng node abs ids."""
    paths = {}
    with open(path) as f:
        for line in f:
            if line[0] != "P":
                continue
            fields = line.rstrip("\n").split("\t")
            nodes = set()
            for s in fields[2].split(","):
                s = s.strip()
                if not s:
                    continue
                ident = s.rstrip("+-")
                try:
                    v = int(ident)
                except ValueError:
                    continue
                if v <= NUM_SYNCMER_NODES:
                    nodes.add(abs(v))
            paths[fields[1]] = nodes
    return paths


def read_gfa_links_components(path):
    """union-find over S/L, plus the paths' segment lists."""
    uf = UnionFind()
    paths = {}
    with open(path) as f:
        for line in f:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "L":
                uf.union(fields[1], fields[3])
            elif fields[0] == "P":
                steps = [s.strip().rstrip("+-") for s in fields[2].split(",")
                         if s.strip()]
                paths[fields[1]] = steps
                for a, b in zip(steps, steps[1:]):
                    uf.union(a, b)
    return uf, paths


class ExportCache:
    def __init__(self):
        self.built = set()
        for fn in os.listdir(EXPORT):
            if fn.endswith(".gfa"):
                self.built.add(int(fn.replace("partition", "")
                                   .replace(".gfa", "")))
        self.cache = {}

    def nodes(self, partition):
        if partition not in self.cache:
            self.cache[partition] = read_gfa_paths(
                f"{EXPORT}/partition{partition}.gfa")
        return self.cache[partition]


def row_name(path_name, start, end):
    return f"{path_name}:{start}-{end}"


def main():
    comp = sys.argv[1] if len(sys.argv) > 1 else "chrI"
    axis_path = f"S288C#0#{comp}"
    out_path = f"{D}/locality-seam-census-{comp}.txt"

    # the windows (axis rows ranked by start)
    axis_rows = []
    for fn in sorted(os.listdir(MAPS)):
        if not fn.endswith(".gfa.map.json"):
            continue
        m = json.load(open(os.path.join(MAPS, fn)))
        for mem in m["members"]:
            if mem["path_name"] == axis_path:
                axis_rows.append((mem["start"], mem["end"], m["partition"]))
    axis_rows.sort()

    # the panel tiling per path, per partition (the BEDs)
    print("# loading the partition BEDs (19,421 partitions)...",
          file=sys.stderr)
    tiling = defaultdict(list)
    for fn in os.listdir(BEDDIR):
        if not fn.endswith(".bed"):
            continue
        try:
            pid = int(fn.replace("partition", "").replace(".bed", ""))
        except ValueError:
            continue
        with open(os.path.join(BEDDIR, fn)) as f:
            for line in f:
                a = line.rstrip("\n").split("\t")
                if len(a) < 3:
                    continue
                tiling[a[0]].append((pid, int(a[1]), int(a[2])))
    for p in tiling:
        tiling[p].sort(key=lambda t: (t[1], t[2]))

    ec = ExportCache()
    lines = []
    tot_paths = tot_extended = 0
    tot_added_bp = 0
    for w, (s, e, part) in enumerate(axis_rows):
        uf, gpaths = read_gfa_links_components(
            f"{PGGB}/partition{part}/graph.gfa")
        aname = row_name(axis_path, s, e)
        first_seg = gpaths[aname][0]
        acomp = uf.find(first_seg)
        universe = set()
        for rn, steps in gpaths.items():
            if any(uf.find(seg) == acomp for seg in steps):
                universe.add(rn)
        # the axis walk's node set
        axis_nodes = ec.nodes(part).get(row_name(axis_path, s, e), set())
        # per universe path: chain from the component row
        ext_paths = 0
        added_bp = 0
        untestable = 0
        sk1_verdict = "no-extension"
        for rn in sorted(universe):
            path_name = rn.rsplit(":", 1)[0]
            if path_name == axis_path:
                continue  # the axis fold is the window's own row
            rows = tiling.get(path_name)
            if not rows:
                continue
            # index the component rows of this path in this partition
            comp_rows = [i for i, (pid, rs, re) in enumerate(rows)
                         if pid == part and row_name(path_name, rs, re) in universe]
            if not comp_rows:
                continue
            # BFS along the tiling from each component row; include a
            # chained row iff it carries a shared axis-walk node
            included = set(comp_rows)
            frontier = list(comp_rows)
            while frontier:
                i = frontier.pop()
                _, cs, ce = rows[i]
                for j in range(len(rows)):
                    if j in included:
                        continue
                    pid2, rs2, re2 = rows[j]
                    # adjacency: the intervals overlap or abut
                    adj = rs2 <= ce and re2 >= cs
                    if not adj:
                        continue
                    if pid2 == part:
                        # same-partition rows are already component
                        # candidates; only universe rows chain on
                        if row_name(path_name, rs2, re2) in universe:
                            included.add(j)
                            frontier.append(j)
                        continue
                    if pid2 not in ec.built:
                        untestable += 1
                        continue
                    nodes = ec.nodes(pid2).get(
                        row_name(path_name, rs2, re2))
                    if nodes is None:
                        untestable += 1
                        continue
                    if nodes & axis_nodes:
                        included.add(j)
                        frontier.append(j)
            # the added material: chained rows beyond the component
            # rows' own span
            comp_span = (min(rows[i][1] for i in comp_rows),
                         max(rows[i][2] for i in comp_rows))
            ext_rows = [rows[i] for i in included
                        if i not in set(comp_rows)]
            shared_ext = [r for r in ext_rows
                          if r[0] in ec.built and ec.nodes(r[0]).get(
                              row_name(path_name, r[1], r[2]),
                              set()) & axis_nodes]
            if shared_ext:
                ext_paths += 1
                lo = min(comp_span[0], min(r[1] for r in shared_ext))
                hi = max(comp_span[1], max(r[2] for r in shared_ext))
                added_bp += (hi - lo) - (comp_span[1] - comp_span[0])
                if path_name.startswith(TRUTH_SECOND):
                    sk1_verdict = (
                        f"EXTENDED +{(hi - lo) - (comp_span[1] - comp_span[0])}bp "
                        f"over {len(shared_ext)} severed rows")
        tot_paths += len(universe)
        tot_extended += ext_paths
        tot_added_bp += added_bp
        lines.append(
            f"window {w}\tpartition {part}\tuniverse_paths {len(universe)}\t"
            f"extended {ext_paths}\tadded_bp {added_bp}\t"
            f"untestable_boundary {untestable}\tsk1 {sk1_verdict}")
        print(lines[-1], file=sys.stderr)
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"# universe paths total {tot_paths}, extended {tot_extended}, "
          f"added bp total {tot_added_bp}")
    print(f"# written {out_path}")


main()
