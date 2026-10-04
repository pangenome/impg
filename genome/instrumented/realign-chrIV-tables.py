#!/usr/bin/env python3
"""The chrIV end-to-end tables (slice 2 of the chrIV stage):

  (1) THE DOMAIN-SCALING TABLE — per locus: the partition's member
      rows, the post-fold candidate count, the unordered class-pair
      count, the unit count, the record count, the factorized
      placements — beside the chrI profile (per-locus means, the
      largest loci), with the linear-or-superlinear verdict computed
      from the numbers (chrIV is ~6.7x chrI in window count and bp).
  (2) THE TRUTH-RANK TABLE — per locus: the slice-1 expressibility
      class (the census receipt of record), the instrument's truth
      expressibility, truth rank, log gap, the winner identity (the
      best folds' member-strain sets), QUAL, the tied-class count;
      the aggregate rank-1 counts by class, the named groups.
  (3) THE WALL/MEMORY LADDER — the serial and 4-wide run markers'
      external walls and RSS poller peaks, the receipts' in-process
      phase walls (arithmetic shown: every phase sum checked against
      the whole-instrument wall), beside the chrI 4-wide profile and
      the pre-run extrapolations.

Receipt-side only: reads the committed receipts and run markers, no
instrument input, no thresholds (the class statements come from the
committed slice-1 census receipt; the ~6.7x scale factors are measured
bp/window ratios, stated as such).

Usage: realign-chrIV-tables.py [--skip-tables]   (prints to stdout and
writes realign-chrIV-tables.txt beside the receipts)
"""
import json
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"
CHRIV_CENSUS = f"{GRAPHS}/partition-graph-chrIV-census.jsonl"

out_lines = []


def say(text=""):
    print(text, flush=True)
    out_lines.append(text)


def load_receipt(prefix, component):
    per = {}
    for line in open(f"{D}/{prefix}.jsonl"):
        d = json.loads(line)
        per[d["locus"]] = d
    return per


def member_rows_by_partition():
    import glob
    rows = {}
    for path in glob.glob(f"{GRAPHS}/*.gfa.map.json"):
        m = json.load(open(path))
        rows[m["partition"]] = m["counts"]["members"]
    return rows


def winner_identity(d):
    """The winner class's fold member strain sets, compact."""
    names = []
    for index in d["best_fold_indices"]:
        fold = d["fold_identities"][index]
        strains = sorted({m["path_name"].split("#")[0] for m in fold["members"]})
        names.append("/".join(strains[:4]) + (f"+{len(strains)-4}" if len(strains) > 4 else ""))
    return names


def main():
    skip_tables = "--skip-tables" in sys.argv

    # ------------------------------------------------ the receipts
    iv = load_receipt("realign-exhaustive-chrIV", "chrIV")
    iv_serial = load_receipt("realign-par-serial-chrIV", "chrIV")
    ic = load_receipt("realign-exhaustive-chrI", "chrI")
    members = member_rows_by_partition()

    census = {}
    for line in open(CHRIV_CENSUS):
        r = json.loads(line)
        if r.get("question") == "a":
            census[r["locus"]] = r["verdict"]

    # ------------------------------------------ (1) the domain scaling
    say("== (1) THE DOMAIN-SCALING TABLE (chrIV, all 167 loci; per-locus domain sizes)")
    say("locus  partition  member_rows  folds  class_pairs  units  records  placements")
    tot_rows = tot_folds = tot_pairs = tot_units = tot_records = tot_place = 0
    per_locus = []
    for locus in sorted(iv):
        d = iv[locus]
        rows = members[d["partition"]]
        pairs = d["class_count"]
        place = d["validation"]["factorized_checked"]
        per_locus.append(
            (locus, d["partition"], rows, d["folds"], pairs, d["unit_count"],
             d["record_count"], place)
        )
        tot_rows += rows
        tot_folds += d["folds"]
        tot_pairs += pairs
        tot_units += d["unit_count"]
        tot_records += d["record_count"]
        tot_place += place
    for (locus, partition, rows, folds, pairs, units, records, place) in per_locus:
        say(f"L{locus:<5} {partition:<9} {rows:<12} {folds:<6} {pairs:<12} "
            f"{units:<6} {records:<8} {place}")
    say(f"TOTALS: member rows {tot_rows}, folds {tot_folds}, class pairs {tot_pairs}, "
        f"units {tot_units}, records {tot_records}, factorized placements {tot_place}")

    ic_pairs = sum(d["class_count"] for d in ic.values())
    ic_units = sum(d["unit_count"] for d in ic.values())
    ic_records = sum(d["record_count"] for d in ic.values())
    ic_place = sum(d["validation"]["factorized_checked"] for d in ic.values())
    n_iv = len(iv)
    n_ic = len(ic)
    say("")
    say("THE CHRIV vs CHRI PROFILE (the scaling verdict, measured):")
    say(f"  windows:                {n_iv} vs {n_ic}  ({n_iv/n_ic:.2f}x)")
    say(f"  class pairs (post-fold): {tot_pairs} vs {ic_pairs}  ({tot_pairs/ic_pairs:.2f}x)")
    say(f"  units:                   {tot_units} vs {ic_units}  ({tot_units/ic_units:.2f}x)")
    say(f"  records:                {tot_records} vs {ic_records}  ({tot_records/ic_records:.2f}x)")
    say(f"  placements:             {tot_place} vs {ic_place}  ({tot_place/ic_place:.2f}x)")
    say(f"  per-locus class pairs:   {tot_pairs/n_iv:.0f} vs {ic_pairs/n_ic:.0f}")
    say(f"  per-locus units:         {tot_units/n_iv:.0f} vs {ic_units/n_ic:.0f}")
    biggest = sorted(per_locus, key=lambda t: -t[4])[:10]
    say("  THE LARGEST LOCI (class pairs; the repeat/subtelomeric neighborhoods):")
    for (locus, partition, rows, folds, pairs, units, records, place) in biggest:
        kind = census.get(locus, "?")
        say(f"    L{locus}: partition {partition}, member rows {rows}, folds {folds}, "
            f"class pairs {pairs}, units {units}, placements {place} [{kind}]")
    if not skip_tables:
        pass  # tables already emitted above

    # ------------------------------------------- (2) the truth-rank table
    say("")
    say("== (2) THE TRUTH-RANK TABLE (all 167 loci; classes from the slice-1 census receipt)")
    say("locus  class                expressible  rank   log_gap          winner  QUAL    tied")
    agg = {}
    rank1 = []
    expressible = []
    for locus in sorted(iv):
        d = iv[locus]
        cls = census.get(locus, "?")
        ex = d["truth_pair_expressible"]
        rank = d["truth_rank"]
        gap = d["log_gap"]
        qual = "unbounded" if d["qual"] is None else f"{d['qual']:.2f}"
        win = "|".join(winner_identity(d))
        tied = d["truth_tied_classes"] if d["truth_tied_classes"] is not None else "-"
        say(f"L{locus:<5} {cls:<20} {str(ex):<12} {str(rank):<6} "
            f"{('%.2f' % gap) if gap is not None else 'n/a':<16} {win:<28} {qual:<8} {tied}")
        agg.setdefault((cls, ex), []).append(locus)
        if ex:
            expressible.append(locus)
            if rank == 1:
                rank1.append(locus)
    say("")
    say("THE AGGREGATE (expressibility x class, with truth rank-1 counts):")
    for (cls, ex), loci in sorted(agg.items()):
        r1 = [l for l in loci if iv[l]["truth_rank"] == 1]
        say(f"  {cls} expressible={ex}: {len(loci)} loci, truth rank-1 {len(r1)}")
    say(f"  TOTAL: {len(expressible)} expressible of {n_iv}; truth rank-1 {len(rank1)}: {rank1}")

    # ------------------------------------ (3) the wall/memory ladder
    say("")
    say("== (3) THE WALL/MEMORY LADDER (external walls + poller RSS peaks; in-process phase sums)")

    def run_stats(tag):
        exit_code = int(open(f"{D}/run-{tag}.exit").read().strip())
        wall = int(open(f"{D}/run-{tag}.wall").read().strip())
        peak = 0
        for line in open(f"{D}/run-{tag}.rss"):
            peak = max(peak, int(line.split()[2]))
        return exit_code, wall, peak

    for tag, label in (
        ("realignpar-serial-chrIV-chrIV", "chrIV serial (the identity base)"),
        ("realignexhaustive-chrIV-chrIV", "chrIV 4-wide (the exhaustive run of record)"),
    ):
        try:
            exit_code, wall, peak = run_stats(tag)
        except FileNotFoundError:
            say(f"  {label}: markers missing (not yet run)")
            continue
        say(f"  {label}: exit {exit_code}, external wall {wall}s, poller RSS peak "
            f"{peak/1048576:.2f} GB")
        # in-process phase sums (the serial receipt's per-locus walls)
        src = iv_serial if "serial" in tag else iv
        phases = {}
        for d in src.values():
            for key, value in d["walls"].items():
                if key.endswith("_seconds"):
                    phases[key] = phases.get(key, 0.0) + value
        say("    in-process phase sums (per-locus, from the receipts):")
        for key in sorted(phases):
            say(f"      {key}: {phases[key]:.1f}s")
        say(f"      loci-phase total: {sum(v for k, v in phases.items() if k != 'locus_seconds'):.1f}s; "
            f"locus_seconds sum {phases.get('locus_seconds', 0):.1f}s")

    # the chrI 4-wide reference
    try:
        exit_code, wall, peak = run_stats("realignpar-4-chrI")
        say(f"  chrI 4-wide reference: exit {exit_code}, wall {wall}s, peak {peak/1048576:.2f} GB")
    except FileNotFoundError:
        say("  chrI 4-wide reference: markers missing")

    with open(f"{D}/realign-chrIV-tables.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    say("")
    say(f"(tables written to {D}/realign-chrIV-tables.txt)")


if __name__ == "__main__":
    main()
