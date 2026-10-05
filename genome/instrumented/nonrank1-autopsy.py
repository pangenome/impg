#!/usr/bin/env python3
"""THE CAUSAL AUTOPSY OF THE 562 NON-RANK-1 CALLS (the owner's demand:
WHY are they called with gap-dominated ~10% structural divergence).

Measurement only — NO FIXES IMPLEMENTED. Receipt-side: streams the
committed main per-locus receipts (realign-exhaustive-<C>.jsonl),
joins the committed sequence-QV receipts (realign-sequence-qv-<C>.
jsonl, the QV stage of record), re-aligns the CHOSEN assignment's two
pairs under the same biWFA machinery with the helper's --cigar mode
(default protocol byte-identical; the counts must equal the stored QV
stage's pairs exactly — a built-in cross-check), and reads the census
receipts for the row-geometry facts. The ingredients sidecars are
never opened. chrIX EXCLUDED (its exhaustive run is still in flight).

THE HYPOTHESIS UNDER TEST (stated as hypothesis — measured, not
assumed): EXTENT ASYMMETRY — partition/window boundaries fragment the
truth's in-domain rows while rival rows span longer across the
seams, collecting orphaned read mass the fragmented truth cannot
place (the unplaced-mass advantage), so the likelihood correctly
prefers the longer explainer, and the QV comparison over mismatched
extents produces end-gaps.

PER NON-RANK-1 EXPRESSIBLE LOCUS, measured:

  (1) THE EXTENT COMPARISON — the called folds' covered intervals and
      lengths (fold_identities members + length) vs the truth folds';
      which is longer and by how much (called_len_sum - truth_len_sum);
      the per-end overhangs in bases FROM THE ALIGNMENT (leading /
      trailing gap columns per side = the called-vs-truth extent
      overhang at each end; the aligned overlap is the columns
      between them).

  (2) THE ALIGNMENT GAP ANATOMY — from the walked biWFA CIGAR of each
      chosen pair: the leading and trailing gap-column runs (END-GAPS,
      the extent overhang) vs the interior gap columns (TRUE MATERIAL
      SUBSTITUTION) vs the mismatch columns, per locus in bases; the
      interior gap runs; the end-gap share of all edits.

  (3) THE LIKELIHOOD DECOMPOSITION — from the receipt's own
      log_gap_decomposition (the instrument's exact semantics): the
      BOTH-PLACED EVIDENCE TERM (both_placed_evidence_gap = the
      winner-minus-truth summed over units BOTH pairs place; negative
      = truth-favored) and the UNPLACED-DIFFERENTIAL TERM (residual =
      log_gap - both_placed_evidence_gap, exactly the winner-minus-
      truth over the units NOT both place: the winner's gain on units
      only it places minus the truth's gain on units only it places;
      units neither places score at the same derived ELSEWHERE branch
      E for both and cancel exactly). Per locus: who wins each term;
      the unplaced masses and unit counts; the E-scale of the mass
      asymmetry (mass_asym * E, the stated magnitude scale of the
      unplaced term, NOT its exact LL — the exact LL of the
      one-side-placed units is not in the receipt).

  (4) THE CENSUS CLASS and the row-geometry facts — the census
      verdict (the fleet-table convention, re-derived from the census
      receipt), the window, the axis partition, and the SK1 window
      ortholog's tiling placement partitions (the strongest-anchors
      hit's placement list — the seam-split / partition-boundary cut
      points named per locus); the truth folds' member rows (the
      in-domain row set), member counts, and whether the truth folds
      carry the SK1 rows.

  (5)+(6) THE AGGREGATE VERDICT AND THE FIX MAPPING — the --tables
      pass: the signature cross-tab, the primary-cause assignment
      (stated precedence), the per-cause QV profiles, and the fix
      mapping with the measured expected-conversion counts (a locus
      converts under a mass-re-attribution class of repair iff the
      both-placed evidence already favors the truth: after full
      neutralization of the unplaced-differential term the shared
      evidence decides, and it decides for the truth exactly when
      both_placed_evidence_gap < 0).

Standing rules: truth assessment-side only; no thresholds (the
classification conventions are stated, like the census's 50bp rule);
scoreboard machinery unmodified; no selection swap, no PR push; the
102GB receipts streamed, never slurped; ingredients sidecars
untouchable.

Usage:
  nonrank1-autopsy.py --component C   (per-component pass; writes
        realign-nonrank1-autopsy-C.jsonl beside the receipts)
  nonrank1-autopsy.py --tables        (the aggregate; reads every
        per-locus jsonl, writes realign-nonrank1-autopsy-tables.txt)
"""
import json
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"

# the committed QV stage of record, imported for its Panel / Gfa /
# census machinery (the same dual-gated material derivation)
import importlib.util

_spec = importlib.util.spec_from_file_location(
    "qvmod", os.path.join(HERE, "realign-sequence-qv.py")
)
qvmod = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(qvmod)

Panel = qvmod.Panel
gfa = qvmod.gfa
census_classes = qvmod.census_classes
fold_strains = qvmod.fold_strains
AlignFailure = qvmod.AlignFailure
COMPONENTS = qvmod.COMPONENTS
HELPER = qvmod.HELPER


class CigarAligner:
    """The biWFA helper in --cigar mode: counts + the walked CIGAR op
    string (the default-mode protocol is unchanged; the counts must
    equal the committed QV stage's stored pair counts — asserted)."""

    def __init__(self):
        self.proc = subprocess.Popen(
            [HELPER, "--cigar"],
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
        fields = resp.split("\t")
        matches, mismatches, ins, dels = (int(x) for x in fields[:4])
        return {
            "matches": matches,
            "mismatches": mismatches,
            "ins": ins,
            "dels": dels,
            "columns": matches + mismatches + ins + dels,
            "edits": mismatches + ins + dels,
            "cigar": fields[4],
        }

    def close(self):
        self.proc.stdin.close()
        self.proc.wait()


MATCH_OPS = ("M", "=", "X")


def cigar_anatomy(pair):
    """The end-gap vs interior-gap anatomy of one walked CIGAR.

    Orientation: the pair was aligned as (called, truth) — a 'D' col
    consumes the CALLED sequence (the called material extends beyond
    the truth's covered extent at that end), an 'I' col consumes the
    TRUTH sequence. The leading run = the gap columns before the
    first match/mismatch column; the trailing run = after the last.
    """
    cigar = pair["cigar"]
    n = len(cigar)
    first = next((i for i, c in enumerate(cigar) if c in MATCH_OPS), None)
    if first is None:
        # no aligned column at all: the entire CIGAR is end-gap
        return {
            "called_overhang_5p": cigar.count("D"),
            "truth_overhang_5p": cigar.count("I"),
            "called_overhang_3p": 0,
            "truth_overhang_3p": 0,
            "interior_gap_cols": 0,
            "interior_gap_runs": 0,
            "aligned_overlap_cols": 0,
        }
    last = max(i for i, c in enumerate(cigar) if c in MATCH_OPS)
    lead, trail = cigar[:first], cigar[last + 1:]
    interior = cigar[first:last + 1]
    runs = 0
    in_run = False
    for c in interior:
        if c in ("I", "D"):
            if not in_run:
                runs += 1
                in_run = True
        else:
            in_run = False
    return {
        "called_overhang_5p": lead.count("D"),
        "truth_overhang_5p": lead.count("I"),
        "called_overhang_3p": trail.count("D"),
        "truth_overhang_3p": trail.count("I"),
        "interior_gap_cols": interior.count("I") + interior.count("D"),
        "interior_gap_runs": runs,
        "aligned_overlap_cols": len(interior),
    }


def fold_extent(fold, max_rows=8):
    """A fold's covered rows and length (the extent facts)."""
    members = fold["members"]
    rows = [
        {"path": m["path_name"], "start": m["start"], "end": m["end"]}
        for m in members
    ]
    return {
        "length": fold["length"],
        "n_members": len(members),
        "rows": rows[:max_rows],
        "paths": sorted({m["path_name"].split("#")[0] for m in members}),
    }


def census_geometry(classes, component):
    """The per-locus row-geometry facts from the census receipt: the
    window, the axis partition, the verdict detail, and the SK1
    window ortholog's tiling placement (the strongest-anchors hit's
    placement partitions — the seam-split / cut points)."""
    path = f"{GRAPHS}/partition-graph-{component}-census.jsonl"
    geo = {}
    if not os.path.exists(path):
        return geo
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            if r.get("question") != "a":
                continue
            hits = r.get("sk1_window_hits") or []
            best = None
            if hits:
                best = max(hits, key=lambda h: h.get("anchors", 0))
            placement = (best or {}).get("placement") or []
            geo[r["locus"]] = {
                "window": r.get("window"),
                "axis_partition": r.get("axis_partition"),
                "windowed_kind": r.get("windowed_kind"),
                "ortholog_hit_anchors": (best or {}).get("anchors"),
                "ortholog_partitions": sorted(
                    {p["partition"] for p in placement}
                ),
                "ortholog_partition_count": len(
                    {p["partition"] for p in placement}
                ),
                "ortholog_hit_interval": (best or {}).get("interval"),
            }
    return geo


def load_qv_receipts(component):
    """The committed QV stage's per-locus records, by locus (small;
    the join source for the stored assignment/pairs/QV)."""
    path = f"{D}/realign-sequence-qv-{component}.jsonl"
    out = {}
    if not os.path.exists(path):
        return out
    with open(path) as f:
        for line in f:
            r = json.loads(line)
            out[r["locus"]] = r
    return out


def locus_autopsy(d, panel, aligner, qv_rec, geo_rec, gates):
    """The per-locus autopsy record, or None when the locus is not an
    expressible non-rank-1 locus."""
    if not d["truth_pair_expressible"]:
        return None
    rank = d["truth_rank"]
    if rank is None or rank == 1:
        return None
    if qv_rec is None:
        raise AlignFailure(f"locus {d['locus']}: no QV-stage record")
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
        # the GFA containment gate (the checker phase-2 derivation);
        # the full-member fold-criterion gate is the committed QV
        # stage's, cross-checked below via the stored pair counts.
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

    # the 2x2 matrix + the assignment, EXACTLY the QV stage's rule; the
    # counts must equal the stored pairs (the cross-check).
    pair_cache = {}

    def pair(ci, ti):
        key = (called[ci], truth[ti])
        if key not in pair_cache:
            pair_cache[key] = aligner.align(seqs[called[ci]], seqs[truth[ti]])
        return pair_cache[key]

    p00, p11 = pair(0, 0), pair(1, 1)
    p01, p10 = pair(0, 1), pair(1, 0)
    for name, stored in (
        ("called0_truth0", p00), ("called1_truth1", p11),
        ("called0_truth1", p01), ("called1_truth0", p10),
    ):
        s = qv_rec["pairs"][name]
        for k in ("matches", "mismatches", "ins", "dels"):
            if stored[k] != s[k]:
                raise AlignFailure(
                    f"locus {d['locus']}: pair {name} count {k} "
                    f"{stored[k]} != stored QV {s[k]}"
                )
    orders = {"direct": (p00, p11), "crossed": (p01, p10)}
    scores = {}
    for name, (pa, pb) in orders.items():
        edits = pa["edits"] + pb["edits"]
        columns = pa["columns"] + pb["columns"]
        error = (edits / columns) if columns else 0.0
        scores[name] = (error, edits, columns)
    chosen = min(
        ("direct", "crossed"),
        key=lambda n: (scores[n][0], scores[n][1], 0 if n == "direct" else 1),
    )
    if chosen != qv_rec["assignment"]:
        raise AlignFailure(
            f"locus {d['locus']}: assignment {chosen} != stored "
            f"{qv_rec['assignment']}"
        )

    # (2) the anatomy of the chosen assignment's two pairs
    anat = [cigar_anatomy(p) for p in orders[chosen]]
    called_overhang = sum(
        a["called_overhang_5p"] + a["called_overhang_3p"] for a in anat
    )
    truth_overhang = sum(
        a["truth_overhang_5p"] + a["truth_overhang_3p"] for a in anat
    )
    end_gap_cols = called_overhang + truth_overhang
    interior_gap_cols = sum(a["interior_gap_cols"] for a in anat)
    mismatches = sum(p["mismatches"] for p in orders[chosen])
    total_edits = end_gap_cols + interior_gap_cols + mismatches
    if total_edits != qv_rec["total_edits"]:
        raise AlignFailure(
            f"locus {d['locus']}: anatomy edits {total_edits} != stored "
            f"QV {qv_rec['total_edits']}"
        )

    # (1) the extent comparison
    called_extents = [fold_extent(d["fold_identities"][fi]) for fi in called]
    truth_extents = [fold_extent(d["fold_identities"][fi]) for fi in truth]
    called_len_sum = sum(e["length"] for e in called_extents)
    truth_len_sum = sum(e["length"] for e in truth_extents)

    # (3) the likelihood decomposition (the receipt's own fields)
    lg = d["log_gap_decomposition"]
    log_gap = d["log_gap"]
    both_placed = lg["both_placed_evidence_gap"]
    residual = log_gap - both_placed
    mass_asym = lg["truth_unplaced_mass"] - lg["winner_unplaced_mass"]
    e_elsewhere = d["scoring"]["elsewhere_log_prob"]

    rec = {
        "component": d.get("component"),
        "locus": d["locus"],
        "partition": partition,
        "class": qv_rec["class"],
        "truth_rank": rank,
        "log_gap": log_gap,
        "truth_log_likelihood": d["truth_log_likelihood"],
        "best_log_likelihood": d["best_log_likelihood"],
        "qv": qv_rec["qv"],
        "identity": qv_rec["identity"],
        "error": qv_rec["error"],
        # (1) the extent comparison
        "called_folds": [
            {"index": fi, **e} for fi, e in zip(called, called_extents)
        ],
        "truth_folds": [
            {"index": fi, **e} for fi, e in zip(truth, truth_extents)
        ],
        "called_len_sum": called_len_sum,
        "truth_len_sum": truth_len_sum,
        "len_delta": called_len_sum - truth_len_sum,
        "called_longer": called_len_sum > truth_len_sum,
        # (2) the alignment gap anatomy, in bases
        "anatomy_pairs": [
            {k: v for k, v in a.items()} for a in anat
        ],
        "end_gap_cols": end_gap_cols,
        "interior_gap_cols": interior_gap_cols,
        "interior_gap_runs": sum(a["interior_gap_runs"] for a in anat),
        "mismatches": mismatches,
        "total_edits": total_edits,
        "called_overhang_cols": called_overhang,
        "truth_overhang_cols": truth_overhang,
        "aligned_overlap_cols": sum(a["aligned_overlap_cols"] for a in anat),
        "end_gap_share_of_edits": (
            end_gap_cols / total_edits if total_edits else 0.0
        ),
        "net_end_gap_len_delta": called_overhang - truth_overhang,
        # (3) the likelihood decomposition
        "both_placed_evidence_gap": both_placed,
        "both_placed_term_winner": (
            "truth" if both_placed < 0
            else "winner" if both_placed > 0
            else "tie"
        ),
        "residual_unplaced_gap": residual,
        "unplaced_term_winner": (
            "winner" if residual > 0
            else "truth" if residual < 0
            else "tie"
        ),
        "truth_unplaced_mass": lg["truth_unplaced_mass"],
        "truth_unplaced_units": lg["truth_unplaced_units"],
        "winner_unplaced_mass": lg["winner_unplaced_mass"],
        "winner_unplaced_units": lg["winner_unplaced_units"],
        "unplaced_mass_asym": mass_asym,
        "unplaced_units_asym": (
            lg["truth_unplaced_units"] - lg["winner_unplaced_units"]
        ),
        "elsewhere_log_prob": e_elsewhere,
        "e_scale_mass_asym": mass_asym * e_elsewhere,
    }
    # (4) the census class + row geometry
    if geo_rec is not None:
        for k, v in geo_rec.items():
            rec[f"census_{k}"] = v
    # the truth folds' SK1 membership (the in-domain row-set fact)
    rec["truth_folds_carry_sk1"] = any(
        any(m["path_name"].startswith("SK1") for m in
            d["fold_identities"][fi]["members"])
        for fi in truth
    )
    # the signatures (measured; the conventions stated in --tables)
    rec["s1_unplaced_mass_advantage"] = mass_asym > 0
    rec["s2_both_placed_truth_favored"] = both_placed < 0
    rec["s3_end_gap_dominated"] = end_gap_cols > interior_gap_cols
    rec["s4_interior_gap_dominated"] = interior_gap_cols > end_gap_cols
    return rec


def run_component(component):
    receipt_path = f"{D}/realign-exhaustive-{component}.jsonl"
    if not os.path.exists(receipt_path):
        print(f"{component}: NO RECEIPT, skipped", flush=True)
        return False
    out_path = f"{D}/realign-nonrank1-autopsy-{component}.jsonl"
    panel = Panel()
    aligner = CigarAligner()
    gates = {"folds_verified": 0}
    qv_by_locus = load_qv_receipts(component)
    classes = census_classes(component)
    geo = census_geometry(classes, component)
    n_loci = n_autopsy = 0
    with open(receipt_path) as fin, open(out_path, "w") as fout:
        for line in fin:
            d = json.loads(line)
            n_loci += 1
            qv_rec = qv_by_locus.get(d["locus"])
            geo_rec = geo.get(d["locus"])
            rec = locus_autopsy(d, panel, aligner, qv_rec, geo_rec, gates)
            if rec is None:
                continue
            rec["component"] = component
            n_autopsy += 1
            fout.write(json.dumps(rec) + "\n")
    aligner.close()
    print(
        f"{component}: {n_loci} loci streamed, {n_autopsy} non-rank-1 "
        f"expressible autopsied, {gates['folds_verified']} folds gated, "
        f"{aligner.count} biWFA alignments; wrote {out_path}",
        flush=True,
    )
    return True


# --------------------------------------------------------- the tables
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


def primary_cause(r):
    """The stated precedence (most specific structural class first):

      1. contig-end      — the census class absent/contig-end (SK1's
                           contig ends short of the window);
      2. extent-asymmetry — the hypothesis's full measured shape: the
                           winner holds the unplaced-mass advantage,
                           the truth wins the both-placed evidence,
                           and the anatomy is end-gap-dominated;
      3. repeat-domain substitution (unplaced-advantage, interior
                           anatomy) — the winner holds the unplaced
                           advantage and the truth wins the shared
                           evidence, but the divergence is interior
                           (true material substitution: the called
                           fold is a diverged copy, not an extent
                           mismatch);
      4. rival-dominant  — the winner wins the both-placed evidence
                           too (no shared-evidence rescue exists);
      5. shared-evidence rival — no unplaced-mass advantage at all;
                           the rival wins purely on jointly placed
                           reads;
      6. boundary        — the rare exact ties / pure-mismatch shapes
                           the data names.
    """
    if r["class"].startswith("absent/"):
        return "contig-end"
    s1 = r["s1_unplaced_mass_advantage"]
    s2 = r["s2_both_placed_truth_favored"]
    s3 = r["s3_end_gap_dominated"]
    s4 = r["s4_interior_gap_dominated"]
    if s1 and s2 and s3:
        return "extent-asymmetry"
    if s1 and s2 and s4:
        return "repeat-domain substitution (unplaced-adv, interior)"
    if s1 and s2:
        return "boundary"
    if s1:
        return "rival-dominant (winner wins both-placed too)"
    return "shared-evidence rival (no unplaced advantage)"


def fix_lever(r, cause):
    """The ordered fix mapping (stated precedence, most specific
    structural disease first; a locus is counted by ONE lever)."""
    if cause == "contig-end":
        return "contig-end handling"
    cls = r["class"]
    if cls == "tiled-elsewhere/neighbor":
        return "partition-spine seam repair"
    if cause.startswith("repeat-domain"):
        return "repeat-domain repair"
    if cause == "extent-asymmetry":
        return "extent normalization of the comparison domain"
    # rival-dominant / shared-evidence / boundary
    return "no mass-re-attribution lever (rival wins shared evidence too)"


def run_tables():
    out_lines = []
    say = lambda t="": (print(t, flush=True), out_lines.append(t))
    records = []
    for component in COMPONENTS:
        path = f"{D}/realign-nonrank1-autopsy-{component}.jsonl"
        if not os.path.exists(path):
            say(f"WARNING: {component} autopsy receipt absent")
            continue
        with open(path) as f:
            for line in f:
                records.append(json.loads(line))
    say("== THE CAUSAL AUTOPSY OF THE NON-RANK-1 CALLS — THE AGGREGATE OF RECORD")
    say(f"   components: {len(COMPONENTS)} closed (chrIX EXCLUDED — its exhaustive")
    say("   run is still in flight); cohort: every truth-pair-expressible locus")
    say("   whose truth rank is not 1 (the QV stage's non-rank-1 set)")
    say(f"   records: {len(records)}")
    say("")
    say("THE CONVENTIONS (the honest statement):")
    say("  * material/alignment: the QV stage's own (biWFA End2End, likegt")
    say("    penalties, the chosen injective assignment); this stage re-aligns")
    say("    the chosen pairs with the helper's --cigar mode and ASSERTS the")
    say("    counts equal the committed QV receipts pair-for-pair (every locus")
    say("    passed; the assignment re-derived identically).")
    say("  * anatomy: the leading/trailing gap-column runs of the walked CIGAR")
    say("    are END-GAPS (the extent overhang: 'D' cols consume the CALLED")
    say("    sequence, 'I' cols the TRUTH); the gap columns between the first")
    say("    and last match/mismatch column are INTERIOR GAPS (true material")
    say("    substitution); mismatches are counted separately, in bases.")
    say("  * likelihood terms: both_placed_evidence_gap is the receipt's own")
    say("    winner-minus-truth over the units BOTH pairs place (negative =")
    say("    truth-favored); the residual_unplaced_gap = log_gap - that term")
    say("    is EXACTLY the winner-minus-truth over the units not both place")
    say("    (units neither places score at the same derived E and cancel);")
    say("    e_scale_mass_asym = mass_asym * E is the stated MAGNITUDE SCALE")
    say("    of the unplaced term, not its exact LL.")
    say("  * signatures: S1 unplaced-mass advantage (truth cannot place more")
    say("    mass than the winner); S2 both-placed truth-favored; S3/S4 the")
    say("    end-gap-dominated / interior-gap-dominated anatomy. The primary")
    say("    cause assigns each locus ONE cause under the stated precedence")
    say("    (contig-end class first, then the extent-asymmetry shape, then")
    say("    the interior-substitution shape, then the rival-dominant and")
    say("    shared-evidence shapes, then the boundary shapes).")
    say("")

    say("== (1) THE EXTENT MEASUREMENT (the hypothesis's domain)")
    longer = sum(1 for r in records if r["called_longer"])
    shorter = sum(1 for r in records if r["len_delta"] < 0)
    equal = len(records) - longer - shorter
    say(f"   called diplotype LONGER than truth: {longer}/{len(records)} "
        f"({longer/len(records)*100:.1f}%)   shorter: {shorter}   equal: {equal}")
    say(f"   median len_delta (called-truth, bp): "
        f"{median([r['len_delta'] for r in records]):.0f}")
    say(f"   median called overhang (end-gap cols on the called side): "
        f"{median([r['called_overhang_cols'] for r in records]):.0f}")
    say(f"   median truth overhang (end-gap cols on the truth side): "
        f"{median([r['truth_overhang_cols'] for r in records]):.0f}")
    say(f"   the end-gap share of all edits: median "
        f"{median([r['end_gap_share_of_edits'] for r in records])*100:.1f}%   "
        f"loci end-gap-dominated (S3): "
        f"{sum(1 for r in records if r['s3_end_gap_dominated'])}   "
        f"interior-dominated (S4): "
        f"{sum(1 for r in records if r['s4_interior_gap_dominated'])}")
    say(f"   the net length delta explained by end gaps vs interior gaps "
        f"(medians): {median([r['net_end_gap_len_delta'] for r in records]):.0f} "
        f"vs {median([r['interior_gap_cols'] for r in records]):.0f} bp")
    say("")

    say("== (2) THE LIKELIHOOD TERM WINNERS (who wins which term)")
    say("   the both-placed evidence term: "
        f"truth-favored {sum(1 for r in records if r['s2_both_placed_truth_favored'])}"
        f"   winner-favored "
        f"{sum(1 for r in records if r['both_placed_term_winner']=='winner')}")
    say("   the unplaced-differential term: "
        f"winner-favored {sum(1 for r in records if r['unplaced_term_winner']=='winner')}"
        f"   truth-favored {sum(1 for r in records if r['unplaced_term_winner']=='truth')}")
    say(f"   the unplaced-mass asymmetry (truth mass the winner also fails to place subtracted): "
        f"median {median([r['unplaced_mass_asym'] for r in records]):.0f} reads")
    say(f"   positive (S1: the truth cannot place more mass) at "
        f"{sum(1 for r in records if r['s1_unplaced_mass_advantage'])}/{len(records)} loci")
    say("")

    say("== (3) THE SIGNATURE CROSS-TAB (the hypothesis's shape measured)")
    say("S1(unplaced-adv)  S2(shared->truth)  S3(end-dominated)  n  medianQV  "
        "median len_delta")
    from collections import defaultdict
    sig = defaultdict(list)
    for r in records:
        key = (r["s1_unplaced_mass_advantage"],
               r["s2_both_placed_truth_favored"],
               r["s3_end_gap_dominated"])
        sig[key].append(r)
    for key in sorted(sig, reverse=True):
        rs = sig[key]
        say(f"{str(key[0]):<16} {str(key[1]):<17} {str(key[2]):<17} "
            f"{len(rs):<4} {median([r['qv'] for r in rs]):<9.2f} "
            f"{median([r['len_delta'] for r in rs]):.0f}")
    say("   THE HYPOTHESIS'S FULL SHAPE (S1 & S2 & S3): "
        f"{len(sig.get((True, True, True), []))} of {len(records)} "
        f"({len(sig.get((True, True, True), []))/len(records)*100:.1f}%)")
    say("")

    say("== (4) THE PRIMARY CAUSE TABLE (the causal decomposition, stated precedence)")
    causes = defaultdict(list)
    for r in records:
        causes[primary_cause(r)].append(r)
    say("cause  n  frac  medianQV  medianIdentity%  medianEndGapShare%  "
        "medianMassAsym  medianLenDelta")
    for cause in sorted(causes, key=lambda c: -len(causes[c])):
        rs = causes[cause]
        say(f"{cause:<52} {len(rs):<4} {len(rs)/len(records)*100:<5.1f} "
            f"{median([r['qv'] for r in rs]):<9.2f} "
            f"{median([r['identity'] for r in rs])*100:<16.2f} "
            f"{median([r['end_gap_share_of_edits'] for r in rs])*100:<19.1f} "
            f"{median([r['unplaced_mass_asym'] for r in rs]):<15.0f} "
            f"{median([r['len_delta'] for r in rs]):.0f}")
    say("")

    say("== (5) THE CAUSE x CENSUS-CLASS CROSS-TAB")
    cc = defaultdict(int)
    for r in records:
        cc[(primary_cause(r), r["class"])] += 1
    classes_seen = sorted({r["class"] for r in records})
    say("cause \\ class: " + "  ".join(classes_seen))
    for cause in sorted(causes, key=lambda c: -len(causes[c])):
        say(f"{cause:<52} " + "  ".join(
            str(cc.get((cause, cls), 0)) for cls in classes_seen))
    say("")

    say("== (6) THE ROW-GEOMETRY FACTS (the truth's row set, per cause)")
    for cause in sorted(causes, key=lambda c: -len(causes[c])):
        rs = causes[cause]
        with_geo = [r for r in rs if r.get("census_axis_partition") is not None]
        if not with_geo:
            say(f"{cause}: pilot component (no census receipt)")
            continue
        multi = [r for r in with_geo
                 if r.get("census_ortholog_partition_count", 0) > 1]
        say(f"{cause}: SK1 window ortholog tiles >1 partition at "
            f"{len(multi)}/{len(with_geo)}; median tiling partitions "
            f"{median([r['census_ortholog_partition_count'] for r in with_geo]):.0f}; "
            f"truth folds carry SK1 rows at "
            f"{sum(1 for r in rs if r['truth_folds_carry_sk1'])}/{len(rs)}")
    say("")

    say("== (7) THE FIX MAPPING (NO FIXES IMPLEMENTED — the measured expected conversions)")
    say("   THE CONVERSION CRITERION (stated): a mass-re-attribution class of")
    say("   repair (extent normalization / seam repair / repeat-domain repair)")
    say("   converts a locus to truth rank-1 iff, with the unplaced-differential")
    say("   term neutralized, the both-placed evidence decides for the truth:")
    say("   both_placed_evidence_gap < 0. Loci with both_placed >= 0 are NOT")
    say("   convertible by mass re-attribution alone (the rival also wins the")
    say("   jointly placed reads) — named, not hidden.")
    convertible = [r for r in records if r["s2_both_placed_truth_favored"]]
    say(f"   total convertible under full neutralization: {len(convertible)}/{len(records)}")
    say("")
    say("   lever (stated precedence: the most specific structural disease first)  "
        "addressed  convertible  not-convertible")
    lever_counts = defaultdict(lambda: [0, 0])
    for r in records:
        cause = primary_cause(r)
        lever = fix_lever(r, cause)
        if r["s2_both_placed_truth_favored"]:
            lever_counts[lever][0] += 1
        else:
            lever_counts[lever][1] += 1
    order = sorted(lever_counts, key=lambda l: -lever_counts[l][0])
    for lever in order:
        conv, notconv = lever_counts[lever]
        say(f"   {lever:<66} {conv + notconv:<10} {conv:<12} {notconv}")
    say("")
    say("   THE HONEST ORDERING: repair the structural diseases in the order")
    for lever in order:
        if lever_counts[lever][0]:
            say(f"     {lever_counts[lever][0]:>4} expected conversions  <- {lever}")
    say("")

    say("== (8) THE PER-LOCUS NAMED EXTREMES (the largest end-gap and")
    say("   interior-gap loci, and the largest unplaced-mass asymmetries)")
    for label, key in (("end-gap cols", "end_gap_cols"),
                       ("interior-gap cols", "interior_gap_cols"),
                       ("unplaced-mass asym", "unplaced_mass_asym")):
        say(f"   by {label}:")
        for r in sorted(records, key=lambda r: -r[key])[:5]:
            say(f"     {r['component']} L{r['locus']}: {r[key]:.0f} "
                f"(QV {r['qv']:.2f}, len_delta {r['len_delta']}, "
                f"both_placed {r['both_placed_evidence_gap']:.1f}, "
                f"class {r['class']})")
    say("")
    say("(measurement only — no fixes implemented; the per-locus receipts")
    say(" realign-nonrank1-autopsy-<C>.jsonl carry the full per-locus extent,")
    say(" anatomy, decomposition and row-geometry records)")
    with open(f"{D}/realign-nonrank1-autopsy-tables.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"wrote {D}/realign-nonrank1-autopsy-tables.txt", flush=True)


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
