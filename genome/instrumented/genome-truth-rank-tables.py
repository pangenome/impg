#!/usr/bin/env python3
"""THE WHOLE-GENOME TRUTH-RANK TABLE (the owner's go: all 17 components
of the balanced-diploid validation under the marginal realignment
instrument).

Per component (chrMT, chrI, chrIV are the closed committed records;
the 14 fleet components run the chrIV recipe end-to-end):
  * the per-locus truth expressibility + truth rank-1 counts, by the
    component census receipt's expressibility class where one exists;
  * the component aggregate: truth rank-1 / expressible / loci;
  * the measured walls with the arithmetic shown (the multicensus run,
    the anchor projection, the serial identity base, the 4-wide
    exhaustive run of record) and the external poller RSS peaks;
  * pending components reported as pending (the script is re-runnable
    mid-fleet; the done-markers carry the state).

Receipt-side only; no thresholds; no instrument inputs; the class
statements come from the committed census receipts.
Usage: genome-truth-rank-tables.py   (prints and writes
genome-truth-rank-tables.txt at the validation dir)
"""
import json
import os
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20230930"
D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{D}/partition-graphs"
# (Stage 2: --receipt-prefix P points the table at an alternative
# receipt family - the territory-normalized receipts; --out-suffix S
# keeps the committed aggregate of record untouched.)
import sys as _sys
RECEIPT_PREFIX = "realign-exhaustive"
OUT_SUFFIX = ""
if "--receipt-prefix" in _sys.argv:
    RECEIPT_PREFIX = _sys.argv[_sys.argv.index("--receipt-prefix") + 1]
if "--out-suffix" in _sys.argv:
    OUT_SUFFIX = _sys.argv[_sys.argv.index("--out-suffix") + 1]

COMPONENTS = [
    "chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
    "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII",
    "chrXIV", "chrXV", "chrXVI",
]

out_lines = []


def say(text=""):
    print(text, flush=True)
    out_lines.append(text)


def marker(tag, field):
    """A run marker's value, or None when absent."""
    path = f"{D}/run-{tag}.{field}"
    if not os.path.exists(path):
        return None
    value = open(path).read().strip()
    return value if value else None


def poller_peak(tag):
    path = f"{D}/run-{tag}.rss"
    if not os.path.exists(path):
        return None
    peak = 0
    for line in open(path):
        fields = line.split()
        if len(fields) >= 3 and fields[2].isdigit():
            peak = max(peak, int(fields[2]))
    return peak or None


def component_walls(c):
    """The measured walls + RSS peaks, from the run markers of record
    (the tag names differ between the committed pilots and the fleet:
    named explicitly, never guessed)."""
    census_tag = (
        f"cosine-multicensus-pilot-{c}"
        if c in ("chrMT", "chrI")
        else f"multicensus-pilot-{c}"
    )
    def _first_existing(field, *tags):
        for t in tags:
            v = marker(t, field)
            if v is not None:
                return t
        return None

    anchor_primary = _first_existing(
        "exit", f"anchorproj-{c}", f"anchorproj-{c}-{c}"
    )
    anchor_tag = anchor_primary or f"anchorproj-{c}"
    serial_tag = {  # chrMT's identity base is the committed slice-E
        # serial run's own marker (the phase-3 ladder re-baselines came
        # later); chrI's is the phase-3 serial re-baseline marker
        "chrMT": "realignexhaustive-chrMT-chrMT",
    }.get(c, f"realignpar-serial-{c}-{c}")
    wide_tag = (
        f"realignpar-4-{c}-{c}" if c in ("chrMT", "chrI")
        else f"realignexhaustive-{c}-{c}"
    )
    phases = []
    for label, tag in (
        ("census", census_tag),
        ("anchor", anchor_tag),
        ("serial", serial_tag),
        ("4-wide", wide_tag),
    ):
        phases.append((label, tag, marker(tag, "wall"), poller_peak(tag)))
    return phases


def main():
    say("== THE WHOLE-GENOME TRUTH-RANK TABLE (17 components, the balanced-diploid validation)")
    say("component  loci  expressible  truth-rank1  census-in-axis  neighbor  foreign  contig-end")
    total_loci = total_expressible = total_rank1 = 0
    per_component = {}
    pending = []
    for c in COMPONENTS:
        receipt = f"{D}/{RECEIPT_PREFIX}-{c}.jsonl"
        if not os.path.exists(receipt):
            pending.append(c)
            say(f"{c:<10} PENDING (no exhaustive receipt)")
            continue
        loci = expressible = rank1 = 0
        for line in open(receipt):
            d = json.loads(line)
            loci += 1
            if d["truth_pair_expressible"]:
                expressible += 1
                if d["truth_rank"] == 1:
                    rank1 += 1
        # the census classes (the 15 components with census receipts;
        # chrMT/chrI's class records live in the committed slice-E docs)
        classes = {}
        census = f"{GRAPHS}/partition-graph-{c}-census.jsonl"
        if os.path.exists(census):
            counts = {"IN-AXIS-PARTITION": 0, "NEIGHBOR": 0, "FOREIGN": 0, "CONTIG-END": 0}
            for line in open(census):
                r = json.loads(line)
                if r.get("question") != "a":
                    continue
                if r["verdict"] == "IN-AXIS-PARTITION":
                    counts["IN-AXIS-PARTITION"] += 1
                elif r["absent_class"] == "CONTIG-END-LENGTH-POLYMORPHISM":
                    counts["CONTIG-END"] += 1
                elif r["split_class"] == "NEIGHBOR-AXIS":
                    counts["NEIGHBOR"] += 1
                else:
                    counts["FOREIGN"] += 1
            classes = counts
        per_component[c] = (loci, expressible, rank1, classes)
        total_loci += loci
        total_expressible += expressible
        total_rank1 += rank1
        cls = classes or {}
        say(
            f"{c:<10} {loci:<5} {expressible:<12} {rank1:<13} "
            f"{cls.get('IN-AXIS-PARTITION', '-'): <14} {cls.get('NEIGHBOR', '-'):<9} "
            f"{cls.get('FOREIGN', '-'):<8} {cls.get('CONTIG-END', '-')}"
        )
    say("")
    say(f"THE AGGREGATE (closed components): loci {total_loci}, expressible "
        f"{total_expressible}, truth rank-1 {total_rank1}")
    if pending:
        say(f"   PENDING: {pending}")
    say("")
    say("== THE MEASURED WALLS (external run markers; RSS = the external poller peak)")
    say("component  census  anchor  serial  4-wide  |  census_GB  anchor_GB  serial_GB  4wide_GB")
    wall_sums = {"census": 0, "anchor": 0, "serial": 0, "4-wide": 0}
    have_wall = True
    for c in COMPONENTS:
        phases = component_walls(c)
        row = [c]
        missing = False
        for label, tag, wall, peak in phases:
            if wall is None:
                row.append("-")
                missing = True
            else:
                row.append(f"{wall}s")
                if not missing:
                    wall_sums[label] += int(wall)
            row_peaks = None
        peaks = []
        for label, tag, wall, peak in phases:
            peaks.append("-" if peak is None else f"{peak/1048576:.1f}")
        if missing:
            say(f"{c:<10} " + "  ".join(row[1:]) + "   (markers pending)")
            have_wall = False
        else:
            say(f"{c:<10} {row[1]:<8} {row[2]:<8} {row[3]:<8} {row[4]:<8} |  "
                + "  ".join(f"{p:>9}" for p in peaks))
    say("")
    if have_wall:
        total = sum(wall_sums.values())
        say("THE SUMMED WALLS (arithmetic shown, closed components):")
        say(f"  census {wall_sums['census']}s + anchor {wall_sums['anchor']}s + "
            f"serial {wall_sums['serial']}s + 4-wide {wall_sums['4-wide']}s = {total}s "
            f"({total/60:.1f} min)")
        say(f"  the fleet-critical path (census + anchor + 4-wide, the serial base a gate "
            f"artifact): {wall_sums['census'] + wall_sums['anchor'] + wall_sums['4-wide']}s")
    else:
        say("  (wall sums deferred: markers still pending above)")
    say("")
    say("== THE HONEST WHOLE-GENOME VERDICT (vs the grounded chrIV-based "
        "~20-min-at-4-wide extrapolation)")
    if pending:
        say(f"  MID-FLEET: {17 - len(pending)} of 17 components closed; the verdict below "
            f"covers the closed set; rerun when the done-markers land.")
    if total_loci:
        closed_mb = sum(
            per_component[c][0] for c in per_component
        )
        say(f"  closed components: {len(per_component)} of 17, {total_loci} loci, "
            f"truth rank-1 {total_rank1}/{total_expressible} expressible "
            f"({100.0 * total_rank1 / max(1, total_expressible):.1f}% of expressible)")
    with open(f"{D}/genome-truth-rank-tables{OUT_SUFFIX}.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    say(f"(written to {D}/genome-truth-rank-tables{OUT_SUFFIX}.txt)")


if __name__ == "__main__":
    main()
