#!/usr/bin/env python3
"""territory-conversion-count.py -- the stage-2 fleet deliverable table:
the extent-normalization convention's conversion count vs the causal
autopsy's prediction, per component and per autopsy cause.

For every component: stream the COMMITTED default-convention receipts
and the TERRITORY-NORMALIZED receipts side by side; per locus state
the before/after truth rank and log gap (CONVERT / HOLD / REGRESS /
moved / unchanged); join the non-rank-1-before loci with the causal
autopsy's per-locus receipts and re-derive the primary cause under the
autopsy's stated precedence (contig-end class first, then the extent
shape S1+S2+end-dominated, then the interior shape S1+S2+interior-
dominated, then rival-dominant S1+not-S2, then the shared-evidence
rival, then the boundary shapes); aggregate the conversions by cause
with the autopsy's 307-loci extent-class expectation beside them.

Usage: territory-conversion-count.py [--out FILE]
"""
import json
import os
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
COMPONENTS = [
    "chrMT", "chrI", "chrIV", "chrII", "chrIII", "chrV", "chrVI",
    "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII",
    "chrXIII", "chrXIV", "chrXV", "chrXVI",
]


def autopsy_cause(record):
    """The autopsy's stated precedence, re-derived from the per-locus
    receipt's own signature fields."""
    if record is None:
        return "pilot-ungated"  # chrMT/chrI carry no census class
    cls = record.get("class") or ""
    if "contig-end" in cls:
        return "contig-end"
    s1 = record.get("s1_unplaced_mass_advantage")
    s2 = record.get("s2_both_placed_truth_favored")
    s3 = record.get("s3_end_gap_dominated")
    s4 = record.get("s4_interior_gap_dominated")
    if s1 and s2 and s3:
        return "extent"
    if s1 and s2 and s4:
        return "repeat-domain"
    if s1 and not s2:
        return "rival-dominant"
    if not s1 and not s2:
        return "shared-evidence"
    return "boundary"


def main():
    out = None
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    lines = []
    rows = []

    def say(text):
        print(text, flush=True)
        lines.append(text)

    say("== THE EXTENT-NORMALIZATION CONVERSION TABLE (the territory convention vs the committed convention)")
    say("component  loci  expr  rank1-before  rank1-after  CONVERT  REGRESS  moved")
    total = {
        "loci": 0, "expr": 0, "rank1_b": 0, "rank1_a": 0,
        "convert": 0, "regress": 0, "moved": 0,
    }
    conv_by_cause = {}
    extent_by_component = {}
    regressions = []
    pending = []
    for c in COMPONENTS:
        committed_path = f"{D}/realign-exhaustive-{c}.jsonl"
        territory_path = f"{D}/realign-exhaustive-territory-{c}.jsonl"
        if not os.path.exists(territory_path):
            pending.append(c)
            say(f"{c:<10} PENDING (no territory receipt)")
            continue
        autopsy = {}
        autopsy_path = f"{D}/realign-nonrank1-autopsy-{c}.jsonl"
        if os.path.exists(autopsy_path):
            for line in open(autopsy_path):
                r = json.loads(line)
                autopsy[r["locus"]] = r
        stats = dict(loci=0, expr=0, rank1_b=0, rank1_a=0, convert=0, regress=0, moved=0)
        comp_rows = []
        with open(committed_path) as fb, open(territory_path) as ft:
            for bline, tline in zip(fb, ft):
                b, t = json.loads(bline), json.loads(tline)
                assert b["locus"] == t["locus"], f"{c}: locus order mismatch"
                stats["loci"] += 1
                if not b["truth_pair_expressible"]:
                    continue
                stats["expr"] += 1
                br, tr = b["truth_rank"], t["truth_rank"]
                stats["rank1_b"] += br == 1
                stats["rank1_a"] += tr == 1
                verdict = "unchanged"
                cause = None
                tie = None
                if br == 1 and tr == 1:
                    verdict = "HOLD"
                    if t.get("called_class_count", 1) > 1:
                        tie = "hold-tied"
                elif br == 1:
                    verdict = "REGRESS"
                    regressions.append((c, b["locus"]))
                elif tr == 1:
                    verdict = "CONVERT"
                    cause = autopsy_cause(autopsy.get(b["locus"]))
                    conv_by_cause[cause] = conv_by_cause.get(cause, 0) + 1
                    # THE TIE HONESTY: under the normalized convention a
                    # converted locus can be a UNIQUE win (the winner IS
                    # the truth pair) or a likelihood TIE (the truth is in
                    # the called set at the max LL, but called[0] — the
                    # committed winner-by-enumeration-order — is a
                    # likelihood-identical rival whose extra row material
                    # is OUT-OF-IMAGE by construction; the likelihood
                    # cannot distinguish the pair on the territory).
                    if t["best_fold_indices"] == t["truth_folds"]:
                        tie = (
                            "unique"
                            if t.get("called_class_count", 1) == 1
                            else "winner-is-truth-tied"
                        )
                    elif t.get("truth_in_called_set"):
                        tie = "truth-tied-not-first"
                    else:
                        tie = "NOT-IN-CALLED-SET"
                elif br != tr:
                    verdict = "moved"
                comp_rows.append(
                    {
                        "component": c,
                        "locus": b["locus"],
                        "rank_before": br,
                        "rank_after": tr,
                        "gap_before": b["log_gap"],
                        "gap_after": t["log_gap"],
                        "verdict": verdict,
                        "autopsy_cause": cause,
                        "convert_kind": tie,
                        "called_class_count": t.get("called_class_count"),
                    }
                )
                if verdict in ("CONVERT", "moved"):
                    stats["moved" if verdict == "moved" else "convert"] += 1
                stats["regress"] += verdict == "REGRESS"
                stats["moved"] += verdict == "moved"
        rows.extend(comp_rows)
        extent_conv = sum(
            1
            for r in comp_rows
            if r["verdict"] == "CONVERT" and r["autopsy_cause"] == "extent"
        )
        extent_by_component[c] = extent_conv
        for key in total:
            total[key] += stats[key]
        say(
            f"{c:<10} {stats['loci']:<5} {stats['expr']:<4} {stats['rank1_b']:<13} "
            f"{stats['rank1_a']:<12} {stats['convert']:<8} {stats['regress']:<8} {stats['moved']}"
        )
    say("")
    say(f"THE NEW GENOME AGGREGATE (closed components): loci {total['loci']}, "
        f"expressible {total['expr']}, truth rank-1 {total['rank1_b']} -> {total['rank1_a']} "
        f"({100.0 * total['rank1_a'] / max(1, total['expr']):.1f}% of expressible)")
    say(f"CONVERT {total['convert']}  REGRESS {total['regress']}  moved {total['moved']}")
    if regressions:
        say(f"REGRESSIONS (design failures): {regressions}")
    say("")
    say("== THE CONVERSIONS BY THE AUTOPSY'S CAUSE (the join with the causal autopsy receipts)")
    for cause in sorted(conv_by_cause, key=lambda k: -conv_by_cause[k]):
        say(f"   {cause:<20} {conv_by_cause[cause]}")
    kinds = {}
    for r in rows:
        if r["verdict"] == "CONVERT":
            kinds[r["convert_kind"]] = kinds.get(r["convert_kind"], 0) + 1
    say(f"   THE CONVERTS' TIE ANATOMY: {kinds} "
        "(unique = the winner IS the truth pair; truth-tied-not-first = the truth "
        "is in the called set at the max LL but called[0] is a "
        "likelihood-identical rival whose extra material is out-of-image)")
    say(f"   the extent-class conversions: {sum(extent_by_component.values())} "
        f"(the autopsy's measured extent class: 307 loci, 274 expected "
        f"conversions under full neutralization)")
    say(f"   extent conversions by component: {extent_by_component}")
    if pending:
        say(f"PENDING components: {pending}")
    if out:
        with open(out, "w") as f:
            for r in rows:
                f.write(json.dumps(r) + "\n")
        say(f"per-locus rows -> {out}")
    with open(f"{D}/territory-conversion-tables.txt", "w") as f:
        f.write("\n".join(lines) + "\n")
    say(f"tables -> {D}/territory-conversion-tables.txt")


if __name__ == "__main__":
    main()
