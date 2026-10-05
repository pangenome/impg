#!/usr/bin/env python3
"""THE CALLED-VS-TRUTH SEQUENCE QV (the owner's correction of the
log-gap "QV-like" pattern — the likegt sequence-QV pattern, NOT the
model-internal log-gap conversion).

Per expressible locus of the closed components (chrIX EXCLUDED: its
exhaustive run is still in flight): align the CALLED diplotype against
the ACTUAL diplotype.

  (1) THE MATERIAL: the called class's material is the winner fold
      pair's spelled sequences (best_fold_indices, the two homolog
      candidates); the truth pair's spelled sequences are the truth
      folds' (truth_folds: the folds carrying the S288C and SK1 truth
      rows). A fold's spelled sequence is re-derived TWO ways and both
      are gated per locus: (a) the PANEL path window — yeast235.txt
      carries the panel as one concatenated sequence with the names
      TSV, so the member row [start, end) window is a direct fetch; the
      fold criterion (all members share the identical row sequence) is
      re-verified against every member's own panel window; (b) the
      PARTITION GFA — the member row's P-line spelling (L-overlap
      trimmed, gap segments included, the committed checker phase-2
      derivation) must CONTAIN the panel window at a front-overhang
      offset. The ingredients sidecars are never opened.

  (2) THE ASSIGNMENT (the owner's yardstick): the two called homologs
      are matched to the two truth homologs by best assignment WITHOUT
      replacement — both orders scored, better taken. Each order is
      scored by the assignment-wide per-base error rate
      (sum of the two pairs' edits) / (sum of the two pairs'
      alignment columns); the lower rate wins; ties broken to the
      lower total edit count, then to the direct order. likegt's own
      rule (the higher MEAN pair identity) is also computed and any
      disagreement is reported, never hidden.

  (3) THE ALIGNMENT: biWFA (lib_wfa2, the workspace's own pinned
      revision, the same machinery likegt uses) — gap-affine End2End,
      likegt's penalties (match 0 / mismatch 4 / gap-open 6 /
      gap-extend 2), Medium memory, no heuristic, Alignment scope.
      The CIGAR is walked against both sequences byte-by-byte.

  (4) THE CONVENTION (the standard Phred form, exact): a mismatch
      column counts one error; a gap column (either side) counts one
      error — indels enter the rate as individual gap columns; the
      End2End span leaves NO unaligned overhang: every base of every
      sequence is in exactly one column, so the denominator is the
      total number of alignment columns over the chosen assignment's
      two pairs (the union of the assignment's aligned material).
      Per-base error = (mismatches + gap columns) / (total columns).
      QV = -10 * log10(error), the likegt floor-and-cap: an error at
      or below 1e-6 (a perfect call, zero edits) is reported as QV 60.
      identity = 1 - error.

Receipt-side only: reads the committed receipts and partition GFAs,
no instrument input, no thresholds, no product change; chrIX excluded.
Usage:
  realign-sequence-qv.py --component C   (per-component pass; writes
        realign-sequence-qv-C.jsonl beside the receipts)
  realign-sequence-qv.py --tables       (the aggregate; reads every
        per-locus jsonl, writes realign-sequence-qv-tables.txt)
"""
import json
import math
import os
import subprocess
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"
PANEL = "/home/erikg/yeast/yeast235.txt"
PANEL_NAMES = "/home/erikg/yeast/yeast235.txt.names.tsv"
HERE = os.path.dirname(os.path.abspath(__file__))
HELPER = f"{HERE}/qv-biwfa/target/release/qv-biwfa"
QV_CAP = 60.0          # the likegt convention: perfect call = QV 60
ERROR_FLOOR = 1e-6     # the likegt floor under the cap

# The 16 closed components (chrIX EXCLUDED — its exhaustive run is
# still in flight; the aggregate states this).
COMPONENTS = [
    "chrMT", "chrI", "chrIV", "chrII", "chrIII", "chrV", "chrVI",
    "chrVII", "chrVIII", "chrX", "chrXI", "chrXII", "chrXIII",
    "chrXIV", "chrXV", "chrXVI",
]


def complement(base):
    return {"A": "T", "C": "G", "G": "C", "T": "A"}.get(base, base)


def revcomp(seq):
    return "".join(complement(b) for b in reversed(seq))


# ----------------------------------------------------- the panel fetch
class Panel:
    """The panel as one concatenated sequence + the names TSV."""

    def __init__(self):
        self.names = {}
        with open(PANEL_NAMES) as f:
            for line in f:
                name, off, length = line.split()
                self.names[name] = (int(off), int(length))
        self.handle = open(PANEL)

    def window(self, path_name, start, end):
        offset, length = self.names[path_name]
        if end > length:
            raise KeyError(
                f"{path_name}: window [{start},{end}) exceeds path length {length}"
            )
        self.handle.seek(offset + start)
        return self.handle.read(end - start)


# ------------------------------------------------- the partition GFAs
class Gfa:
    """One partition GFA: S segments, P lines, L overlaps (the
    committed checker phase-2 derivation, verbatim in spirit)."""

    def __init__(self, partition):
        self.partition = partition
        self.segs = {}
        self.paths = {}
        self.links = {}
        with open(f"{GRAPHS}/partition{partition}.gfa") as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    self.segs[fields[1]] = fields[2]
                elif fields[0] == "P":
                    self.paths[fields[1]] = fields[2].split(",")
                elif fields[0] == "L":
                    self.links[(fields[1], fields[2], fields[3], fields[4])] = int(
                        fields[5].rstrip("M")
                    )

    def spelled(self, row_name):
        """The row's spelled sequence (L-overlap trimmed, gap segments
        included, '-' steps reverse complements)."""
        content = {}
        parts = []
        prev = None
        for step in self.paths[row_name]:
            name, sign = step[:-1], step[-1]
            if (name, sign) not in content:
                content[(name, sign)] = (
                    self.segs[name] if sign == "+" else revcomp(self.segs[name])
                )
            if prev is not None:
                overlap = self.links[(prev[0], prev[1], name, sign)]
                parts.append(content[(name, sign)][overlap:])
            else:
                parts.append(content[(name, sign)])
            prev = (name, sign)
        return "".join(parts)


GFAS = {}


def gfa(partition):
    if partition not in GFAS:
        GFAS[partition] = Gfa(partition)
    return GFAS[partition]


# ------------------------------------------------------ the biWFA helper
class AlignFailure(Exception):
    pass


class BiWfa:
    def __init__(self):
        self.proc = subprocess.Popen(
            [HELPER],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            text=True,
        )
        self.count = 0

    def align(self, a, b):
        self.proc.stdin.write(f"{a}\t{b}\n")
        self.proc.stdin.flush()
        resp = self.proc.stdout.readline().strip()
        self.count += 1
        if resp.startswith("FAIL"):
            raise AlignFailure(resp)
        matches, mismatches, ins, dels = (int(x) for x in resp.split("\t"))
        return {
            "matches": matches,
            "mismatches": mismatches,
            "ins": ins,
            "dels": dels,
            "columns": matches + mismatches + ins + dels,
            "edits": mismatches + ins + dels,
        }

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()


def qv_of(error):
    """QV = -10*log10(error) with the likegt floor-and-cap."""
    return -10.0 * math.log10(max(error, ERROR_FLOOR)) if error > 0 else QV_CAP


# ------------------------------------------------------- census classes
def census_classes(component):
    """The per-locus expressibility class from the component's census
    receipt (the fleet tables' class source; IN-AXIS / NEIGHBOR (the
    tiled-elsewhere class) / FOREIGN / CONTIG-END)."""
    path = f"{GRAPHS}/partition-graph-{component}-census.jsonl"
    classes = {}
    if not os.path.exists(path):
        return classes
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            if r.get("question") != "a":
                continue
            if r["verdict"] == "IN-AXIS-PARTITION":
                cls = "in-axis"
            elif r.get("absent_class") == "CONTIG-END-LENGTH-POLYMORPHISM":
                cls = "contig-end"
            elif r.get("split_class") == "NEIGHBOR-AXIS":
                cls = "neighbor"
            else:
                cls = "foreign"
            classes[r["locus"]] = cls
    return classes


def fold_strains(fold):
    return sorted({m["path_name"].split("#")[0] for m in fold["members"]})


# ------------------------------------------------------ the per-locus QV
def locus_qv(d, panel, aligner, gates):
    """The called-vs-truth per-locus record, or None when the locus
    is not truth-pair-expressible. Every used fold is gated BOTH ways:
    the panel-window fold criterion (all members identical) and the
    GFA containment (the checker phase-2 derivation)."""
    if not d["truth_pair_expressible"]:
        return None
    called = list(d["best_fold_indices"])
    truth = list(d["truth_folds"])
    if len(called) != 2 or len(truth) != 2:
        raise AlignFailure(
            f"locus {d['locus']}: unexpected fold-pair shape "
            f"called={called} truth={truth}"
        )
    partition = d["partition"]
    g = gfa(partition)
    seqs = {}
    for fi in sorted(set(called + truth)):
        fold = d["fold_identities"][fi]
        if fold["length"] <= 0:
            raise AlignFailure(f"locus {d['locus']} fold {fi}: empty fold")
        first = fold["members"][0]
        seq = panel.window(first["path_name"], first["start"], first["end"])
        if len(seq) != fold["length"]:
            raise AlignFailure(
                f"locus {d['locus']} fold {fi}: panel window length "
                f"{len(seq)} != fold length {fold['length']}"
            )
        # (a) the fold criterion re-verified: every member's own panel
        # window is the identical row sequence.
        for m in fold["members"]:
            other = panel.window(m["path_name"], m["start"], m["end"])
            if len(other) != fold["length"] or other != seq:
                raise AlignFailure(
                    f"locus {d['locus']} fold {fi}: member "
                    f"{m['path_name']}:{m['start']}-{m['end']} panel "
                    "window differs from the fold's row sequence"
                )
        # (b) the GFA containment: the panel window appears in the
        # member row's P-line spelling at a front-overhang offset.
        spelled = g.spelled(f"{first['path_name']}:{first['start']}-{first['end']}")
        offset = spelled.find(seq)
        if offset < 0 or spelled[offset:offset + len(seq)] != seq:
            raise AlignFailure(
                f"locus {d['locus']} fold {fi}: panel window absent from "
                "the GFA spelling of "
                f"{first['path_name']}:{first['start']}-{first['end']}"
            )
        gates["folds_verified"] += 1
        seqs[fi] = seq
    # the 2x2 alignment matrix (deduplicated per locus)
    pair_cache = {}

    def pair(ci, ti):
        key = (called[ci], truth[ti])
        if key not in pair_cache:
            pair_cache[key] = aligner.align(seqs[called[ci]], seqs[truth[ti]])
        return pair_cache[key]

    p00, p11 = pair(0, 0), pair(1, 1)
    p01, p10 = pair(0, 1), pair(1, 0)
    orders = {
        "direct": (p00, p11),
        "crossed": (p01, p10),
    }
    scores = {}
    for name, (pa, pb) in orders.items():
        edits = pa["edits"] + pb["edits"]
        columns = pa["columns"] + pb["columns"]
        error = (edits / columns) if columns else 0.0
        scores[name] = (error, edits, columns)
    # the assignment WITHOUT replacement: both orders scored by the
    # assignment-wide per-base error rate, the lower taken; ties to
    # the lower total edits, then the direct order.
    chosen = min(
        ("direct", "crossed"),
        key=lambda n: (scores[n][0], scores[n][1], 0 if n == "direct" else 1),
    )
    # likegt's own rule (the higher mean pair identity) — disagreement
    # is reported, never hidden.
    mean_id = {
        name: sum(
            1.0 - (p["edits"] / p["columns"] if p["columns"] else 0.0)
            for p in orders[name]
        ) / 2.0
        for name in orders
    }
    likegt_choice = "direct" if mean_id["direct"] >= mean_id["crossed"] else "crossed"
    error, edits, columns = scores[chosen]
    pa, pb = orders[chosen]
    return {
        "locus": d["locus"],
        "partition": partition,
        "truth_rank": d["truth_rank"],
        "rank1": d["truth_rank"] == 1,
        "called_folds": [
            {"index": fi, "strains": fold_strains(d["fold_identities"][fi])}
            for fi in called
        ],
        "truth_folds": [
            {"index": fi, "strains": fold_strains(d["fold_identities"][fi])}
            for fi in truth
        ],
        "assignment": chosen,
        "likegt_assignment": likegt_choice,
        "assignment_disagreement": chosen != likegt_choice,
        "pairs": {
            "called0_truth0": p00,
            "called1_truth1": p11,
            "called0_truth1": p01,
            "called1_truth0": p10,
        },
        "chosen_pairs": [pa, pb],
        "total_edits": edits,
        "total_columns": columns,
        "mismatches": pa["mismatches"] + pb["mismatches"],
        "gap_columns": pa["ins"] + pa["dels"] + pb["ins"] + pb["dels"],
        "error": error,
        "identity": 1.0 - error,
        "qv": qv_of(error),
        "perfect": edits == 0,
    }


def run_component(component):
    receipt_path = f"{D}/realign-exhaustive-{component}.jsonl"
    if not os.path.exists(receipt_path):
        print(f"{component}: NO RECEIPT, skipped", flush=True)
        return False
    out_path = f"{D}/realign-sequence-qv-{component}.jsonl"
    panel = Panel()
    aligner = BiWfa()
    gates = {"folds_verified": 0}
    n_loci = n_expressible = 0
    classes = census_classes(component)
    with open(receipt_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            d = json.loads(line)
            n_loci += 1
            rec = None
            try:
                rec = locus_qv(d, panel, aligner, gates)
            except AlignFailure as e:
                print(f"{component} locus {d['locus']}: ALIGNMENT GATE FAIL: {e}",
                      flush=True)
                raise
            if rec is None:
                continue
            n_expressible += 1
            rec["class"] = classes.get(rec["locus"], "pilot (no census receipt)")
            rec["component"] = component
            fout.write(json.dumps(rec) + "\n")
    aligner.close()
    print(
        f"{component}: {n_loci} loci, {n_expressible} truth-pair-expressible, "
        f"{gates['folds_verified']} folds dual-gated (panel fold criterion + "
        f"GFA containment), {aligner.count} biWFA alignments; wrote {out_path}",
        flush=True,
    )
    return True


# --------------------------------------------------------- the tables
def nearest_rank(values, p):
    """The nearest-rank percentile: sorted ascending, the value at
    index int(round(p/100 * (n-1))) — deterministic, stated."""
    s = sorted(values)
    if not s:
        return None
    idx = int(round((p / 100.0) * (len(s) - 1)))
    return s[idx]


def median(values):
    s = sorted(values)
    n = len(s)
    if not n:
        return None
    if n % 2:
        return s[n // 2]
    return (s[n // 2 - 1] + s[n // 2]) / 2.0


def frac(values, threshold):
    if not values:
        return None
    return sum(1 for v in values if v >= threshold) / len(values)


def run_tables():
    out_lines = []
    say = lambda t="": (print(t, flush=True), out_lines.append(t))
    records = []
    for component in COMPONENTS:
        path = f"{D}/realign-sequence-qv-{component}.jsonl"
        if not os.path.exists(path):
            say(f"WARNING: {component} per-locus QV receipt absent")
            continue
        with open(path) as f:
            for line in f:
                records.append(json.loads(line))
    say("== THE CALLED-VS-TRUTH SEQUENCE QV — THE AGGREGATE OF RECORD")
    say(f"   components: {len(COMPONENTS)} closed (chrIX EXCLUDED — its exhaustive")
    say("   run is still in flight); loci: every truth-pair-expressible locus")
    say(f"   records: {len(records)}")
    say("")
    say("THE CONVENTION (the honest statement):")
    say("  * material: called = the winner fold pair's spelled sequences;")
    say("    truth = the truth folds' (the S288C/SK1 rows). Every used fold")
    say("    dual-gated: all members' panel windows identical (the fold")
    say("    criterion) AND the panel window contained in the member row's")
    say("    partition-GFA P-line spelling (the checker phase-2 derivation).")
    say("  * assignment: both orders scored by the assignment-wide per-base")
    say("    error (summed edits / summed columns), the lower taken (the")
    say("    owner's injective no-replacement yardstick); likegt's mean-identity")
    say("    rule computed beside it, any disagreement reported below.")
    say("  * alignment: biWFA gap-affine End2End (likegt's penalties: match 0,")
    say("    mismatch 4, gap-open 6, gap-extend 2), CIGAR walked byte-by-byte.")
    say("  * per-base error = (mismatch columns + gap columns) / (total")
    say("    alignment columns over the chosen assignment's two pairs).")
    say("    Indels enter as individual gap columns (either side); the End2End")
    say("    span leaves NO unaligned overhang — every base of every sequence")
    say("    sits in exactly one column, so the denominator is the union of")
    say("    the assignment's aligned material. QV = -10*log10(error), the")
    say("    likegt floor-and-cap: error <= 1e-6 (a perfect call, zero edits)")
    say("    is reported as QV 60. identity = 1 - error.")
    say("  * percentiles: nearest-rank on the sorted values; median the")
    say("    standard middle (mean of the two middles at even n).")
    say("")
    qvs = [r["qv"] for r in records]
    perfect = sum(1 for r in records if r["perfect"])
    say("== (1) THE GENOME-WIDE QV DISTRIBUTION (all expressible loci)")
    say(f"   n={len(qvs)}  perfect calls (QV 60, zero edits): {perfect}")
    if qvs:
        say(f"   median {median(qvs):.2f}  p10 {nearest_rank(qvs,10):.2f}  "
            f"p90 {nearest_rank(qvs,90):.2f}  min {min(qvs):.2f}  max {max(qvs):.2f}")
        say(f"   QV>=40: {frac(qvs,40)*100:.2f}%   QV>=30: {frac(qvs,30)*100:.2f}%   "
            f"QV>=20: {frac(qvs,20)*100:.2f}%")
    say("")
    say("== (2) THE PER-COMPONENT TABLE")
    say("component  n_expressible  rank1  non-rank1  medianQV  p10QV  p90QV  "
        "minQV  >=40%  >=30%  >=20%")
    for component in COMPONENTS:
        rs = [r for r in records if r["component"] == component]
        if not rs:
            continue
        cq = [r["qv"] for r in rs]
        n1 = sum(1 for r in rs if r["rank1"])
        say(f"{component:<10} {len(rs):<14} {n1:<6} {len(rs)-n1:<10} "
            f"{median(cq):<8.2f} {nearest_rank(cq,10):<6.2f} "
            f"{nearest_rank(cq,90):<6.2f} {min(cq):<6.2f} "
            f"{frac(cq,40)*100:<6.2f} {frac(cq,30)*100:<6.2f} {frac(cq,20)*100:.2f}")
    say("")
    say("== (3) THE REFRACTIVE TABLE: rank-1 calls vs non-rank-1 calls")
    r1 = [r for r in records if r["rank1"]]
    nr1 = [r for r in records if not r["rank1"]]
    for label, rs in (("RANK-1 (the called class IS the truth class)", r1),
                      ("NON-RANK-1 (a different class called)", nr1)):
        q = [r["qv"] for r in rs]
        e = [r["error"] for r in rs]
        say(f"   {label}: n={len(rs)}")
        if rs:
            say(f"     QV: median {median(q):.2f}  p10 {nearest_rank(q,10):.2f}  "
                f"p90 {nearest_rank(q,90):.2f}  min {min(q):.2f}")
            say(f"     per-base error: median {median(e):.6g}  "
                f"p90 {nearest_rank(e,90):.6g}  max {max(e):.6g}")
            say(f"     identity >=99.9%: {sum(1 for r in rs if r['identity']>=0.999)}/{len(rs)}"
                f" ({sum(1 for r in rs if r['identity']>=0.999)/len(rs)*100:.2f}%)   "
                f">=99%: {sum(1 for r in rs if r['identity']>=0.99)}/{len(rs)}   "
                f">=95%: {sum(1 for r in rs if r['identity']>=0.95)}/{len(rs)}")
            say(f"     perfect (zero edits): {sum(1 for r in rs if r['perfect'])}")
    say("")
    say("   THE OWNER'S 'EFFECTIVELY THE SAME' CLAIM, MEASURED:")
    if nr1:
        for thr in (0.999, 0.99, 0.95):
            n = sum(1 for r in nr1 if r["identity"] >= thr)
            say(f"     non-rank-1 calls >= {thr*100:.1f}% identical to the truth "
                f"diplotype: {n}/{len(nr1)} ({n/len(nr1)*100:.2f}%)")
        say(f"     median identity of non-rank-1 calls: {median([r['identity'] for r in nr1])*100:.4f}%")
    say("")
    say("== (4) THE NON-RANK-1 CALLS BY CENSUS CLASS")
    by_class = {}
    for r in nr1:
        by_class.setdefault(r["class"], []).append(r)
    say("class        n  medianQV  p10QV  minQV  medianIdentity%")
    for cls in sorted(by_class):
        rs = by_class[cls]
        q = [r["qv"] for r in rs]
        say(f"{cls:<12} {len(rs):<3} {median(q):<9.2f} {nearest_rank(q,10):<7.2f} "
            f"{min(q):<7.2f} {median([r['identity'] for r in rs])*100:.4f}")
    say("")
    say("== (5) THE WORST CALLS (the non-rank-1 loci by QV ascending, named)")
    say("component  locus  QV  error  mismatches  gapcols  called-folds  truth-folds  class  assignment")
    for r in sorted(nr1, key=lambda r: r["qv"])[:25]:
        called = " | ".join("/".join(f["strains"][:3]) +
                            ("+" if len(f["strains"]) > 3 else "")
                            for f in r["called_folds"])
        truth = " | ".join("/".join(f["strains"][:3]) +
                           ("+" if len(f["strains"]) > 3 else "")
                           for f in r["truth_folds"])
        say(f"{r['component'] or '?':<10} L{r['locus']:<5} {r['qv']:.2f}  "
            f"{r['error']:.6g}  {r['mismatches']:<10} {r['gap_columns']:<7} "
            f"{called:<24} {truth:<24} {r['class']:<8} {r['assignment']}")
    say("")
    dis = [r for r in records if r["assignment_disagreement"]]
    say(f"== (6) THE ASSIGNMENT-RULE CROSS-CHECK: likegt's mean-identity rule "
        f"disagrees with the error-rate rule at {len(dis)} of {len(records)} loci")
    for r in dis[:10]:
        say(f"   {r['component']} L{r['locus']}: error-rate {r['assignment']} vs "
            f"likegt {r['likegt_assignment']}")
    say("")
    say("(the perfect-call convention: zero-edit loci carry QV 60 exactly;")
    say(" their per-pair CIGARs are full match columns, verified by the same")
    say(" biWFA walk — no fast path is taken)")
    with open(f"{D}/realign-sequence-qv-tables.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"wrote {D}/realign-sequence-qv-tables.txt", flush=True)


def main():
    args = sys.argv[1:]
    if "--tables" in args:
        run_tables()
        return
    for i, a in enumerate(args):
        if a == "--component":
            run_component(args[i + 1])
            return
    for component in COMPONENTS:
        run_component(component)


if __name__ == "__main__":
    main()
