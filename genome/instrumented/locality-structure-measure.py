#!/usr/bin/env python3
"""GATE 1 of the alignment-induced locality domain: MEASURE the new
locality structure over the built chrMT/chrI pggb graphs, against the
partition domains of the committed normalized-rule receipts.

THE LOCALITY DEFINITION UNDER MEASUREMENT (exact, derived, no
thresholds):
  * the locality of a window is the alignment-connected structure of
    its axis partition's pggb graph that contains the window's axis
    row - the connected component (undirected, over the graph's own
    segment/link structure) whose segments the axis row's pggb path
    traverses;
  * the locus's universe is the connected material: the member rows
    whose pggb path traverses the locality component (partition
    membership over-joins disconnected material - those rows leave);
  * the seam rejoin: a universe path's fold material extends along
    the path's own partition-BED tiling (abutting/overlapping rows,
    the panel's proven coordinate-continuous spine) to the hull of
    the rows that carry interned syng nodes SHARED with the window's
    axis-row walk - the partition boundary severed the same
    continuous material, and the shared-node coverage names exactly
    the severed piece (the ortholog material the window's territory
    projects onto that path).

This script measures, per window of chrMT/chrI:
  (1) the partition's pggb component census (component count, the
      largest, the axis component's row share, rows leaving the
      universe, rows straddling components);
  (2) the expressibility pre-check: the truth-second path's row in
      the axis component? (both homologs present in the locality);
  (3) the seam-extension census for every universe path: the tiling
      rows outside the axis partition that carry axis-walk shared
      nodes, and the resulting fold hull sizes (the per-locus
      universe size change);
  (4) at the residual (non-rank-1) loci: the winner folds' component
      membership (would over-join separation remove the winner?).

Receipt-side only. Usage: locality-structure-measure.py chrI
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
NAMES = "/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names"
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


class Gfa:
    """Minimal GFA reader: segments (name -> length), links, paths."""

    def __init__(self, path):
        self.path = path
        self.segs = {}
        self.paths = {}
        links = []
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    self.segs[fields[1]] = len(fields[2])
                elif fields[0] == "L":
                    links.append((fields[1], fields[3]))
                elif fields[0] == "P":
                    self.paths[fields[1]] = [
                        s.strip().rstrip("+-") for s in fields[2].split(",")
                        if s.strip()
                    ]
        self.uf = UnionFind()
        for a, b in links:
            self.uf.union(a, b)
        # a path's own traversal also asserts adjacency through its
        # segments: consecutive steps of a P line are connected
        # material by construction (the path spells a real sequence)
        for steps in self.paths.values():
            for a, b in zip(steps, steps[1:]):
                self.uf.union(a, b)

    def row_components(self, row_name):
        out = set()
        for seg in self.paths.get(row_name, []):
            out.add(self.uf.find(seg))
        return out


def row_name(path_name, start, end):
    return f"{path_name}:{start}-{end}"


def export_walk_nodes(partition, row_name_):
    """The row's interned syng node ids from the export partition GFA's
    P line (signed global syncmer ids, gap splices excluded)."""
    gfa = EXPORT_CACHE[partition]
    out = set()
    for ident in gfa.paths.get(row_name_, []):
        if int(ident) <= NUM_SYNCMER_NODES:
            out.add(abs(int(ident)))
    return out


def load_bed_tiling(paths_of_interest):
    """The panel's partition BED rows for the requested paths, per
    partition: path -> list of (partition, start, end)."""
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
                if len(a) < 3 or a[0] not in paths_of_interest:
                    continue
                s, e = int(a[1]), int(a[2])
                tiling[a[0]].append((pid, s, e))
    for p in tiling:
        tiling[p].sort(key=lambda t: (t[1], t[2]))
    return tiling


EXPORT_CACHE = {}


def main():
    comp = sys.argv[1] if len(sys.argv) > 1 else "chrI"
    axis_path = f"S288C#0#{comp}"
    receipt = f"{D}/realign-exhaustive-territory-{comp}.jsonl"
    out_path = f"{D}/locality-structure-{comp}.txt"

    # the windows: the axis rows ranked by start (the scorer's own
    # derivation), one per window, from the projected maps
    axis_rows = []
    for fn in sorted(os.listdir(MAPS)):
        if not fn.endswith(".gfa.map.json"):
            continue
        m = json.load(open(os.path.join(MAPS, fn)))
        for mem in m["members"]:
            if mem["path_name"] == axis_path:
                axis_rows.append((mem["start"], mem["end"], m["partition"]))
    axis_rows.sort()
    print(f"# locality-structure measurement, component {comp}")
    print(f"# windows: {len(axis_rows)}")

    # the receipt's per-locus facts (truth folds, winners, ranks)
    receipt_facts = {}
    for line in open(receipt):
        d = json.loads(line)
        fi = d.get("fold_identities") or []
        facts = {"rank": d.get("truth_rank"),
                 "expressible": d.get("truth_pair_expressible"),
                 "winner_rows": set(), "truth_rows": set()}
        for i in d.get("best_fold_indices") or []:
            for m in fi[i]["members"]:
                facts["winner_rows"].add((m["path_name"], m["start"], m["end"]))
        for i in d.get("truth_folds") or []:
            for m in fi[i]["members"]:
                facts["truth_rows"].add((m["path_name"], m["start"], m["end"]))
        receipt_facts[d["locus"]] = facts

    lines = []
    tot_left = 0
    tot_universe = 0
    tot_members = 0
    for w, (s, e, part) in enumerate(axis_rows):
        gfa = Gfa(f"{PGGB}/partition{part}/graph.gfa")
        EXPORT_CACHE[part] = Gfa(f"{EXPORT}/partition{part}.gfa")
        aname = row_name(axis_path, s, e)
        acomps = gfa.row_components(aname)
        assert len(acomps) >= 1, f"axis row {aname} absent from the pggb graph"
        # the locality component = the component holding the axis
        # row's FIRST segment (the walk's own spine; a row whose
        # material straddles components belongs where its walk
        # starts - stated rule)
        first_seg = gfa.paths[aname][0]
        acomp = gfa.uf.find(first_seg)
        # component census
        comp_rows = defaultdict(set)
        straddle = 0
        for rn in gfa.paths:
            cs = gfa.row_components(rn)
            for c in cs:
                comp_rows[c].add(rn)
            if len(cs) > 1:
                straddle += 1
        universe_rows = sorted(comp_rows[acomp])
        n_members = len(gfa.paths)
        n_left = n_members - len(universe_rows)
        tot_left += n_left
        tot_universe += len(universe_rows)
        tot_members += n_members
        # the truth-second row in this partition
        sk1_rows = [rn for rn in gfa.paths if rn.startswith(TRUTH_SECOND)]
        sk1_in = [rn for rn in sk1_rows if acomp in gfa.row_components(rn)]
        # the winner rows at this locus (if any) - in the locality?
        facts = receipt_facts.get(w)
        win_in = win_out = 0
        if facts:
            for (p, ps, pe) in facts["winner_rows"]:
                rn = row_name(p, ps, pe)
                cs = gfa.row_components(rn) if rn in gfa.paths else set()
                if acomp in cs:
                    win_in += 1
                else:
                    win_out += 1
        lines.append(
            f"window {w}\tpartition {part}\tmembers {n_members}\t"
            f"components {len(comp_rows)}\tlargest {max(len(v) for v in comp_rows.values())}\t"
            f"axis_component_rows {len(universe_rows)}\tleaving {n_left}\t"
            f"straddling {straddle}\t"
            f"sk1_rows {len(sk1_rows)}/in {len(sk1_in)}\t"
            f"rank {facts['rank'] if facts else '-'}\t"
            f"winner_rows in/out {win_in}/{win_out}"
        )
    with open(out_path, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"# members total {tot_members}, universe total {tot_universe}, "
          f"leaving total {tot_left}")
    print("\n".join(lines))
    print(f"# written {out_path}")


main()
