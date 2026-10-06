#!/usr/bin/env python3
"""THE ALIGNMENT-INDUCED LOCALITY DOMAIN - the construction tool (the
owner's go: the candidate domains come from ALIGNMENT-INDUCED
LOCALITIES, not partition memberships).

THE LOCALITY (exact, derived, no thresholds, no tuning constants):
  * the locality of a window is the alignment-connected structure of
    its axis partition's pggb graph that contains the window's axis
    row - the connected component (undirected, over the graph's own
    segment/link structure, a path's own consecutive steps included)
    whose segments the axis row's pggb path traverses. Partition
    membership over-joins disconnected material: those rows leave
    the locus's universe (the partition-541 class).
  * THE FOLD MATERIAL on each universe path = the path-interval hull
    of the MAXIMAL MONOTONE AXIS-CORRESPONDENCE through the path's
    locality rows - the order-preserving (path position, axis
    position) pairs pooled from BOTH committed structures:
      (i) the pggb shared segments (the within-partition real
          alignment: segments traversed by both the axis row's and
          the path's pggb path - divergent pockets sit INSIDE the
          aligned span as variant bubbles, and material beyond the
          alignment - the diverged edge pockets the committed
          edge-refinement could not claim - is EXCLUDED, repairing
          the chrI-L13 winner-hull class), and
     (ii) the interned syng nodes SHARED with the axis row's walk
          (the cross-partition seam coverage: the panel's proven
          coordinate-continuous tiling chains the path's rows across
          partition boundaries, and a chained row joins the
          correspondence only when its shared nodes EXTEND the chain
          monotonically at its ends - the severed ortholog material
          rejoins (the tiled-elsewhere class), interior repeat
          coincidences do not).
    The fold is the path's own material over that hull (the window
    coverage of the correspondence's ends, +K+W at the top); the
    axis path's fold is the window's axis row itself (the
    territory). A path with NO correspondence pair (connected
    through other rows, nothing projecting onto the territory)
    keeps its full component rows - the committed fold, decided by
    the standing image rule.
  * THE RECORDS re-derive per locality: the territory-touch
    convention on the new locality extents - an occurrence votes at
    the locus iff its interval overlaps a locality fold interval on
    the occurrence's own path (the scorer's locality mode).

Node positions come from the committed pggb projection sidecars
(per-row absolute-bp walks, the spell-equality-proven syng
coordinates); segment positions from the pggb GFA's own 0M
cumulative spellings (the pggb paths spell the rows exactly).
Rows in partitions without a pggb build are UNTESTABLE for the seam
chain - the committed build-set closure bounds the seam repair
(truth-side closed by construction), counted honestly in the
structure receipt.

Emit per window: locality-graphs-<C>/window<w>.gfa.map.json (the
scorer's own PartitionMap schema, partition id = the window id) +
the structure receipt (the gate-1 measurement). Receipt-side only;
assessment-side. Usage:
  locality-domain-maps.py chrI [--measure-only]
"""
import json
import os
import sys
from collections import defaultdict

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
PGGB = f"{D}/pggb-partition-graphs-chrMT-chrI"
PROJ = f"{PGGB}/projection"
MAPS = f"{D}/pggb-projected-graphs-chrMT-chrI"
BEDDIR = "/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results"
K, W = 63, 8  # the syng's own syncmer window parameters
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


class PggbGfa:
    """The per-partition pggb build: segments, 0M links (positions are
    cumulative), P paths, union-find components."""

    def __init__(self, path):
        self.segs = {}
        self.paths = {}
        uf = UnionFind()
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    self.segs[fields[1]] = len(fields[2])
                elif fields[0] == "L":
                    uf.union(fields[1], fields[3])
                elif fields[0] == "P":
                    self.paths[fields[1]] = [
                        s.strip().rstrip("+-") for s in fields[2].split(",")
                        if s.strip()
                    ]
        for steps in self.paths.values():
            for a, b in zip(steps, steps[1:]):
                uf.union(a, b)
        self.uf = uf
        self._pos = {}

    def positions(self, row_name):
        """(segment -> cumulative start bp on the row's own spelling;
        the pggb paths spell the rows exactly, offset 0)."""
        if row_name not in self._pos:
            pos = {}
            acc = 0
            for seg in self.paths[row_name]:
                pos[seg] = acc
                acc += self.segs[seg]
            self._pos[row_name] = pos
        return self._pos[row_name]

    def row_components(self, row_name):
        return {self.uf.find(s) for s in self.paths[row_name]}


def row_name(path_name, start, end):
    return f"{path_name}:{start}-{end}"


def lis_core(pairs):
    """The maximal order-preserving (monotone in both coordinates)
    chain through the correspondence pairs - the colinear alignment
    core. Sorted by axis position, the longest strictly increasing
    subsequence on the path position; repeat-coincidence pairs that
    violate the colinear structure drop out of the hull."""
    if not pairs:
        return []
    order = sorted(range(len(pairs)), key=lambda i: (pairs[i][1], pairs[i][0]))
    tails_pb = []   # per chain length: the smallest path bp
    tails_idx = []  # per chain length: the pair index achieving it
    preds = [-1] * len(pairs)
    for idx in order:
        pb = pairs[idx][0]
        lo, hi = 0, len(tails_pb)
        while lo < hi:
            mid = (lo + hi) // 2
            if tails_pb[mid] < pb:
                lo = mid + 1
            else:
                hi = mid
        preds[idx] = tails_idx[lo - 1] if lo > 0 else -1
        if lo == len(tails_pb):
            tails_pb.append(pb)
            tails_idx.append(idx)
        else:
            tails_pb[lo] = pb
            tails_idx[lo] = idx
    chain = []
    cur = tails_idx[-1]
    while cur >= 0:
        chain.append(pairs[cur])
        cur = preds[cur]
    chain.reverse()
    return chain


class LocalityBuilder:
    def __init__(self):
        self.pggb_cache = {}
        self.walk_cache = {}
        self.pggb_built = set()
        for d in os.listdir(PGGB):
            if d.startswith("partition") and os.path.isdir(f"{PGGB}/{d}"):
                self.pggb_built.add(int(d.replace("partition", "")))
        # the panel tiling per path: path -> [(partition, start, end)]
        self.tiling = defaultdict(list)
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
                    self.tiling[a[0]].append((pid, int(a[1]), int(a[2])))
        for p in self.tiling:
            self.tiling[p].sort(key=lambda t: (t[1], t[2]))

    def pggb(self, partition):
        if partition not in self.pggb_cache:
            self.pggb_cache[partition] = PggbGfa(
                f"{PGGB}/partition{partition}/graph.gfa")
        return self.pggb_cache[partition]

    def walks(self, partition):
        """row_name -> {abs node id: absolute bp} from the projection
        sidecar (the spell-equality-proven syng coordinates)."""
        if partition not in self.walk_cache:
            table = {}
            with open(f"{PROJ}/partition{partition}.pggb.projection.jsonl") as f:
                for line in f:
                    r = json.loads(line)
                    rn = r["p"]
                    table[rn] = {}
                    for node, bp in r["w"]:
                        n = abs(int(node))
                        if n not in table[rn]:
                            table[rn][n] = int(bp)
            self.walk_cache[partition] = table
        return self.walk_cache[partition]

    def build_window(self, axis_path, s, e, part, axis_partitions):
        """One window's locality: (members, structure)."""
        pg = self.pggb(part)
        aname = row_name(axis_path, s, e)
        if aname not in pg.paths:
            raise SystemExit(f"axis row {aname} absent from the pggb graph")
        first_seg = pg.paths[aname][0]
        acomp = pg.uf.find(first_seg)
        universe = [rn for rn, steps in pg.paths.items()
                    if any(pg.uf.find(seg) == acomp for seg in steps)]
        axis_pos = pg.positions(aname)
        axis_segset = set(axis_pos)
        axis_node_bp = self.walks(part)[aname]

        members = [{"path_name": axis_path, "start": s, "end": e}]
        structure = {
            "partition": part,
            "universe_rows": len(universe),
            "members": 0,
            "extended_paths": 0,
            "truncated_paths": 0,
            "untestable_rows": 0,
        }
        # (Folds are PER ROW, one fold per universe row, each extended
        # along its own tiling seam chain - NOT merged per path: the
        # committed domain folds identical copies (same sequence, same
        # relative walk) into ONE candidate, and a per-path hull of a
        # path's several rows both breaks that coalescing (an identical
        # twin row merged with a longer row of the same path becomes a
        # separate rival carrying the longer material - the chrMT-L2
        # regression class) and makes the truth's own first-homolog
        # support a rival. The fold interval = the exact extent hull
        # of the correspondence pairs (segment pairs carry their
        # segment length, node pairs the window extent K+W), so
        # identical rows yield identical intervals and coalesce.)
        for rn in sorted(universe):
            path_name, se = rn.rsplit(":", 1)
            rs, re = (int(x) for x in se.split("-"))
            if path_name == axis_path:
                continue  # the axis fold is the axis row verbatim
            rpos = pg.positions(rn)
            # (i) the row's own pggb shared segments (the colinear
            # LIS core; repeat-coincidence segments drop out)
            raw = []
            for seg, pp in rpos.items():
                if seg in axis_segset:
                    raw.append((rs + pp, s + axis_pos[seg],
                                pg.segs[seg]))
            core = lis_core([(p, a) for p, a, _ in raw])
            extent = {p: x for p, a, x in raw}
            pairs = [(p, a, extent[p]) for p, a in core]
            # (ii) the row's own cross-partition seam chain
            tiling = [(pid, ts, te) for (pid, ts, te)
                      in self.tiling.get(path_name, [])
                      if pid != part and pid in axis_partitions
                      and ts < re and te > rs]
            admitted = set()
            changed = True
            while changed and pairs:
                changed = False
                lo_p = min(p for p, a, x in pairs)
                hi_p = max(p + x for p, a, x in pairs)
                lo_a = min(a for p, a, x in pairs)
                hi_a = max(a for p, a, x in pairs)
                for (pid, ts, te) in tiling:
                    key = (pid, ts, te)
                    if key in admitted:
                        continue
                    if not (ts <= hi_p and te >= lo_p):
                        continue
                    if pid not in self.pggb_built:
                        continue
                    nodes = self.walks(pid).get(
                        row_name(path_name, ts, te))
                    if not nodes:
                        continue
                    shared = [(bp, axis_node_bp[n], K + W)
                              for n, bp in nodes.items()
                              if n in axis_node_bp]
                    if not shared:
                        continue
                    grows = any(
                        (ab < lo_a and pb < lo_p)
                        or (ab > hi_a and pb > hi_p)
                        for pb, ab, _ in shared)
                    if not grows:
                        continue
                    admitted.add(key)
                    pairs.extend(
                        (pb, ab, K + W) for pb, ab, _ in shared
                        if (ab < lo_a and pb < lo_p)
                        or (ab > hi_a and pb > hi_p))
                    changed = True
            if pairs:
                lo = min(p for p, a, x in pairs)
                hi = max(p + x for p, a, x in pairs)
                structure["untestable_rows"] += sum(
                    1 for (pid, ts, te) in tiling
                    if pid not in self.pggb_built
                    and ts <= hi and te >= lo
                    and (pid, ts, te) not in admitted)
                if lo < rs or hi > re:
                    structure["extended_paths"] += 1
                if lo > rs or hi < re:
                    structure["truncated_paths"] += 1
                members.append({"path_name": path_name,
                                "start": int(lo), "end": int(hi)})
            else:
                # no correspondence: the committed fold (the full row,
                # the standing image rule decides)
                members.append({"path_name": path_name,
                                "start": rs, "end": re})
        structure["members"] = len(members)
        return members, structure


def comp_start(rn):
    return int(rn.rsplit(":", 1)[1].split("-")[0])


def key_not_admitted(admitted, pid, ts, te):
    return (pid, ts, te) not in admitted


def main():
    comp = sys.argv[1] if len(sys.argv) > 1 else "chrI"
    measure_only = "--measure-only" in sys.argv
    axis_path = f"S288C#0#{comp}"
    outdir = f"{D}/locality-graphs-{comp}"
    if not measure_only:
        os.makedirs(outdir, exist_ok=True)

    axis_rows = []
    for fn in sorted(os.listdir(MAPS)):
        if not fn.endswith(".gfa.map.json"):
            continue
        m = json.load(open(os.path.join(MAPS, fn)))
        for mem in m["members"]:
            if mem["path_name"] == axis_path:
                axis_rows.append((mem["start"], mem["end"], m["partition"]))
    axis_rows.sort()

    lb = LocalityBuilder()
    axis_partitions = {p for (_, _, p) in axis_rows}
    lines = [f"# the alignment-induced locality domain, component {comp}",
             f"# windows {len(axis_rows)}"]
    for w, (s, e, part) in enumerate(axis_rows):
        members, st = lb.build_window(axis_path, s, e, part, axis_partitions)
        sk1 = [m for m in members
               if m["path_name"] == f"SK1#0#{comp}"]
        sk1_desc = ""
        if len(sk1) == 1:
            sk1_desc = (f"\tsk1 [{sk1[0]['start']},{sk1[0]['end']})"
                        f" len {sk1[0]['end'] - sk1[0]['start']}")
        elif len(sk1) > 1:
            sk1_desc = f"\tsk1 {len(sk1)} folds"
        lines.append(
            f"window {w}\tpartition {part}\taxis [{s},{e})\t"
            f"universe_rows {st['universe_rows']}\t"
            f"members {st['members']}\t"
            f"extended {st['extended_paths']}\ttruncated "
            f"{st['truncated_paths']}\tuntestable {st['untestable_rows']}"
            + sk1_desc
        )
        print(lines[-1], file=sys.stderr)
        if not measure_only:
            with open(f"{outdir}/window{w}.gfa.map.json", "w") as f:
                json.dump({
                    "partition": w,
                    "construction": {
                        "material": "the alignment-induced locality domain: "
                        "the pggb component containing the axis row; the "
                        "fold material = the maximal monotone axis "
                        "correspondence (pggb shared segments within the "
                        "partition, interned-node seam chain across the "
                        "panel tiling)",
                        "axis_row": [s, e],
                        "axis_partition": part,
                    },
                    "members": members,
                }, f)
    out = f"{D}/locality-domain-structure-{comp}.txt"
    with open(out, "w") as f:
        f.write("\n".join(lines) + "\n")
    print(f"# written {out}" + ("" if measure_only else f" + {outdir}"))


main()
