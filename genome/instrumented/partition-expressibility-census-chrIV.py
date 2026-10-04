#!/usr/bin/env python3
"""THE EXPRESSIBILITY CENSUS AT CHRIV SCALE (the chrIV end-to-end stage,
slice 1: survey, partition graphs, expressibility — the owner's 2026-11
ruling extended from the chrMT/chrI build to the first 1.5 Mb
chromosome; BUILD AND MEASURE only, no inference rewiring, no scoring).

The census questions at chrIV scale:
  (a) THE 805-CLASS TEST over the balanced-diploid validation's chrIV
      component (167 loci, truth S288C+SK1, truth_pair_in_candidate_domain
      0/167 in the windowed frame): per locus, is the SK1 homolog's
      material PRESENT BY ALIGNMENT in the partition structure — the
      ORTHOLOG statement (the balanced receipt's own SK1 truth interval,
      at the 98 assessable loci) and the POSITIONAL statement (the
      window's coordinates on the SK1 truth contig — a coordinate
      statement, at every locus) — and WHERE does it sit: in the
      window's own axis partition, or in FOREIGN/REPEAT partitions (the
      partition-spine disease class, holders named)? The contig-end
      length-polymorphism class (SK1#0#chrIV ends at 1,486,921 bp, 80 kb
      before S288C#0#chrIV's 1,566,853) is named as its own class, the
      chrMT L13 precedent;
  (c) THE NEAR-TWIN CONDENSATION: chrIV's near-twin route pairs —
      flagged by the stated structural rule below — as shared-node
      paths + variant pockets across the built graphs;
  (d) THE TRUTH PAIR'S JOINT SPELLABILITY: per locus, is BOTH homologs'
      material jointly spellable in the partition graphs by alignment
      (S288C row in the axis partition by construction; SK1 by its
      ortholog/positional placement)?

THE NEAR-TWIN FLAG RULE (stated, structural, pre-scoring — the way the
chrI BTE#3/#4 twins were flagged, derived from the panel's own partition
placement, not from any scoring receipt): a pair of chrIV-family panel
paths is a NEAR-TWIN ROUTE PAIR iff EVERY member row of one path has an
interval-close counterpart (both ends within 50 bp — the stated
coordinate-agreement window, the census's own convention, like the
homology padding) in the same partition of the other path, i.e. the
interval-close coverage of the shorter path's tiled bp is 100%. At
chrIV this flags exactly two pairs (the full >= 0.3 spectrum is emitted
in the survey receipt): AAA#0#chrIV/SGDref#0#chrIV (100% EXACT-identical
rows — the identical-through-graph class, the chrI 5397/5545 analogue)
and CLL#0#chrIV/CLL#1#chrIV (100% interval-close, 0% exact, offsets
within +-16 bp — the near-identical same-strain copy class, the chrMT
5599/5616 analogue); the measured runner-up ALI#0/ALI#1 sits at 72.5%
and the next pair at 48.8% — a natural gap, no tuning.

Inputs (receipts on disk; nothing is re-derived from truth information
the panel does not already carry):
  * the partition membership BEDs (the completed impg partition run,
    commit 295bca9, 19,421 partitions, 100% of panel bp);
  * the chrIV partition GFAs + .map.json sidecars (partition_graph_export,
    the same recipe: raw mode, frequency mask disabled, AGC gap splicing,
    global-syng node-id continuity);
  * the balanced-diploid validation's chrIV receipt
    (genotype-distance-balanced-chrIV.json — the windowed-frame status:
    assessable/bracketed, the SK1 ortholog intervals);
  * the syng homology receipts (query_region_with_anchors from every
    window extent, padding 120) and the walk receipts (the panel's own
    path walks for every truth/twin member row) — both produced by
    syng_homology_dump from the panel index.

Modes:
  requests — derive the build list + write
              partition-graph-chrIV-homology-requests.tsv and
              partition-graph-chrIV-walk-requests.tsv (for the dump runs);
  census   — compute the census, write
              partition-graph-chrIV-census.jsonl and print the tables.

Verdict rules (stated; no thresholds enter anything):
  * a locus's SK1 HOMOLOG statement = the ORTHOLOG interval (the
    balanced receipt's own SK1 truth interval — the alignment-truth
    statement) where the locus is assessable, else the POSITIONAL region
    (the window's coordinates on the SK1 truth contig);
  * PRESENT BY ALIGNMENT = the statement interval is covered by SK1
    member rows of the partition structure (the partition BEDs are
    themselves alignment-induced);
  * IN-AXIS-PARTITION / TILED-ELSEWHERE / PARTIAL-ELSEWHERE / ABSENT
    over that interval; a locus is in the PARTITION-SPINE DISEASE CLASS
    when the SK1 verdict is TILED/PARTIAL-ELSEWHERE (the material is in
    the partition structure but in foreign/repeat partitions, holders
    named), and in the CONTIG-END LENGTH-POLYMORPHISM class when the
    window lies beyond the SK1 truth contig's end.

Assessment-side only. No product file is read or written. The census is
a pure measurement of the panel's own alignment structure.
"""

import collections
import json
import os
import sys

DATA = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
GRAPHS = f"{DATA}/partition-graphs"
BEDS = "/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results"
AXIS_FILE = (
    "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
)
NAMES = "/home/erikg/yeast/syng-k63-s8-seed7-acgt-only/yeast235.syng.names"
BALANCED = f"{DATA}/genotype-distance-balanced-chrIV.json"
PADDING = 120  # the syng probe's own existing choice
INTERVAL_CLOSE = 50  # the stated coordinate-agreement window (bp)
SPECTRUM_FLOOR = 0.3  # the reported spectrum floor (a reporting bound, not a flag rule)
TRUTH0_NAME = "S288C#0#chrIV"
TRUTH1_NAME = "SK1#0#chrIV"

BUILD_LIST = f"{GRAPHS}/partition-graph-chrIV-build-list.txt"
HOMOLOGY_REQUESTS = f"{GRAPHS}/partition-graph-chrIV-homology-requests.tsv"
WALK_REQUESTS = f"{GRAPHS}/partition-graph-chrIV-walk-requests.tsv"
HOMOLOGY_RECEIPT = f"{GRAPHS}/partition-graph-chrIV-homology.jsonl"
WALK_RECEIPT = f"{GRAPHS}/partition-graph-chrIV-walks.jsonl"
CENSUS_RECEIPT = f"{GRAPHS}/partition-graph-chrIV-census.jsonl"


def load_names():
    """id -> (name, length); the names file is the panel's own statement
    of the path lengths (used for the contig-end class)."""
    names = {}
    with open(NAMES) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            names[int(fields[0])] = (fields[1], int(fields[2]))
    return names


def scan_beds():
    """The panel's own alignment-induced membership, in full."""
    by_partition = {}
    by_path = collections.defaultdict(list)
    for entry in sorted(os.listdir(BEDS)):
        if not (entry.startswith("partition") and entry.endswith(".bed")):
            continue
        partition = int(entry[len("partition") : -len(".bed")])
        rows = []
        with open(os.path.join(BEDS, entry)) as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                rows.append((fields[0], int(fields[1]), int(fields[2])))
        by_partition[partition] = rows
        for name, start, end in rows:
            by_path[name].append((partition, start, end))
    return by_partition, by_path


def load_axis():
    """The 167 chrIV windows of the balanced-diploid validation (the
    axis file's S288C#0#chrIV intervals, one per S288C member row)."""
    axis = []
    for interval in json.load(open(AXIS_FILE))["intervals"]:
        if interval["component"] == TRUTH0_NAME:
            axis.append(
                (interval["start"], interval["end"], int(interval["group"][len("partition") :]))
            )
    axis.sort()
    return axis


def load_balanced():
    """The windowed-frame status per locus (the committed balanced
    chrIV receipt): assessable/bracketed + the SK1 ortholog intervals."""
    document = json.load(open(BALANCED))
    rows = {record["locus"]: record for record in document["rows"]}
    return document, rows


def interval_close_stats(rows_by_path, name_a, name_b):
    """Interval-close coverage of a's tiled bp by b's rows in the same
    partitions (both directions reported by the caller)."""
    by_part = collections.defaultdict(list)
    for partition, start, end in rows_by_path.get(name_b, ()):
        by_part[partition].append((start, end))
    matched_bp = exact_bp = total_bp = 0
    for partition, start, end in rows_by_path.get(name_a, ()):
        span = end - start
        total_bp += span
        for start2, end2 in by_part.get(partition, ()):
            if abs(start2 - start) <= INTERVAL_CLOSE and abs(end2 - end) <= INTERVAL_CLOSE:
                matched_bp += span
                if start2 == start and end2 == end:
                    exact_bp += span
                break
    return matched_bp, exact_bp, total_bp


def derive_near_twins(by_path, names):
    """THE FLAG RULE (stated above): full interval-close coverage of the
    shorter path's tiled bp, over the chrIV-family panel paths. Returns
    (flagged pairs as name tuples, the >= SPECTRUM_FLOOR spectrum, the
    chrIV-family path list)."""
    family = sorted(name for name in by_path if name.rsplit("#", 1)[1] == "chrIV")
    spectrum = []
    for i, name_a in enumerate(family):
        for name_b in family[i + 1:]:
            bp_ab, _, total_a = interval_close_stats(by_path, name_a, name_b)
            bp_ba, _, total_b = interval_close_stats(by_path, name_b, name_a)
            if not total_a or not total_b:
                continue
            frac = min(bp_ab / total_a, bp_ba / total_b)
            if frac >= SPECTRUM_FLOOR:
                _, exact_ab, _ = interval_close_stats(by_path, name_a, name_b)
                spectrum.append(
                    {
                        "pair": [name_a, name_b],
                        "interval_close_fraction": round(frac, 4),
                        "exact_fraction": round(exact_ab / total_a, 4),
                    }
                )
    spectrum.sort(key=lambda entry: -entry["interval_close_fraction"])
    flagged = [
        tuple(entry["pair"])
        for entry in spectrum
        if entry["interval_close_fraction"] >= 1.0
    ]
    return flagged, spectrum, family


def build_set(axis, by_path, twin_names):
    """Axis partitions + the truth-side SK1 holders + the near-twin
    holders — the chrIV graph build set (derived, not hardcoded)."""
    partitions = {partition for _, _, partition in axis}
    for name in [TRUTH1_NAME] + list(twin_names):
        partitions |= {p for p, _, _ in by_path.get(name, ())}
    return sorted(partitions)


def placement(rows_by_path, name, start, end):
    """Partitions holding member rows of `name` overlapping [start,end)."""
    hits = []
    for partition, row_start, row_end in rows_by_path.get(name, ()):
        if row_start < end and start < row_end:
            hits.append((partition, row_start, row_end))
    return sorted(hits)


def union_coverage(intervals, start, end):
    """Total bp of [start,end) covered by the union of the intervals."""
    covered = 0
    reach = start
    for row_start, row_end in sorted(intervals):
        if row_end <= reach:
            continue
        covered += min(row_end, end) - max(row_start, reach)
        reach = max(reach, row_end)
        if reach >= end:
            break
    return covered


def statement_verdict(rows, partition, start, end):
    """The stated verdict classes over a statement interval."""
    covered = union_coverage([(s, e) for _, s, e in rows], start, end)
    if any(p == partition for p, _, _ in rows):
        return "IN-AXIS-PARTITION", covered
    if not rows:
        return "ABSENT", covered
    if covered >= end - start:
        return "TILED-ELSEWHERE", covered
    return "PARTIAL-ELSEWHERE", covered


def write_requests():
    names = load_names()
    ids = {name: path for path, (name, _) in names.items()}
    axis = load_axis()
    by_partition, by_path = scan_beds()
    flagged, spectrum, family = derive_near_twins(by_path, names)
    twin_names = sorted({name for pair in flagged for name in pair})
    build = build_set(axis, by_path, twin_names)

    with open(BUILD_LIST, "w") as handle:
        handle.write("\n".join(str(p) for p in build) + "\n")

    truth0, truth1 = ids[TRUTH0_NAME], ids[TRUTH1_NAME]
    homology_lines = []
    for locus, (start, end, partition) in enumerate(axis):
        targets = {truth1}
        for name in twin_names:
            if any(p == partition for p, _, _ in by_path.get(name, ())):
                targets.add(ids[name])
        homology_lines.append(
            f"{truth0}\t{locus}\t{start}\t{end}\t{PADDING}\t"
            + ",".join(str(t) for t in sorted(targets))
        )

    walk_names = {TRUTH0_NAME, TRUTH1_NAME} | set(twin_names)
    walk_lines = []
    for partition in build:
        for name, start, end in by_partition[partition]:
            if name in walk_names:
                walk_lines.append(f"{ids[name]}\t{start}\t{end}\tpartition{partition}")

    with open(HOMOLOGY_REQUESTS, "w") as handle:
        handle.write("\n".join(homology_lines) + "\n")
    with open(WALK_REQUESTS, "w") as handle:
        handle.write("\n".join(walk_lines) + "\n")
    print(f"build set: {len(build)} partitions -> {BUILD_LIST}")
    print(f"near-twin flag: {['/'.join(pair) for pair in flagged]}")
    print(f"wrote {len(homology_lines)} homology requests -> {HOMOLOGY_REQUESTS}")
    print(f"wrote {len(walk_lines)} walk requests -> {WALK_REQUESTS}")


def load_homology():
    hits = collections.defaultdict(list)
    with open(HOMOLOGY_RECEIPT) as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("request_done"):
                continue
            hits[(record["query_path"], record["window"], record["target_path"])].append(record)
    return hits


def load_walks():
    walks = {}
    with open(WALK_RECEIPT) as handle:
        for line in handle:
            record = json.loads(line)
            walks[(record["path"], record["start"], record["end"], record["tag"])] = record["steps"]
    return walks


def load_gfa_paths(partition):
    paths = {}
    with open(f"{GRAPHS}/partition{partition}.gfa") as handle:
        for line in handle:
            if line.startswith("P\t"):
                fields = line.rstrip("\n").split("\t")
                paths[fields[1]] = [
                    (int(step[:-1]), step[-1] == "+") for step in fields[2].split(",")
                ]
    return paths


def anchor_count(hit):
    return hit["anchors"] if isinstance(hit["anchors"], int) else len(hit["anchors"])


def run_census():
    names = load_names()
    ids = {name: path for path, (name, _) in names.items()}
    truth0, truth1 = ids[TRUTH0_NAME], ids[TRUTH1_NAME]
    sk1_len = names[truth1][1]
    axis = load_axis()
    by_partition, by_path = scan_beds()
    balanced, balanced_rows = load_balanced()
    flagged, spectrum, family = derive_near_twins(by_path, names)
    twin_names = sorted({name for pair in flagged for name in pair})
    build = build_set(axis, by_path, twin_names)
    hits = load_homology()
    walks = load_walks()

    census = []

    # ------------------------------------------------------------------
    # (s) the survey record (the stage's question 1, receipt-side)
    # ------------------------------------------------------------------
    sk1_rows = by_path.get(TRUTH1_NAME, ())
    axis_partitions = {p for _, _, p in axis}
    sk1_axis_rows = [r for r in sk1_rows if r[0] in axis_partitions]
    sk1_holder_partitions = sorted({p for p, _, _ in sk1_rows} - axis_partitions)
    twin_holder_partitions = sorted(
        {p for name in twin_names for p, _, _ in by_path.get(name, ())}
        - axis_partitions
        - set(sk1_holder_partitions)
    )
    census.append(
        {
            "question": "s",
            "component": "chrIV",
            "partitions_genome_wide": len(by_partition),
            "chriv_family_paths": len(family),
            "chriv_family_rows": sum(len(by_path.get(n, ())) for n in family),
            "chriv_family_partitions": len(
                {p for n in family for p, _, _ in by_path.get(n, ())}
            ),
            "loci": len(axis),
            "axis_partitions": len(axis_partitions),
            "s288c_rows": len(by_path.get(TRUTH0_NAME, ())),
            "s288c_len_bp": names[truth0][1],
            "sk1_rows": len(sk1_rows),
            "sk1_len_bp": sk1_len,
            "sk1_rows_in_axis_partitions": len(sk1_axis_rows),
            "sk1_holder_partitions": sk1_holder_partitions,
            "twin_pairs": [list(pair) for pair in flagged],
            "twin_holder_partitions": twin_holder_partitions,
            "near_twin_flag_rule": (
                "a pair of chrIV-family paths is flagged iff every member row of one "
                "path has an interval-close counterpart (both ends within "
                f"{INTERVAL_CLOSE} bp) in the same partition of the other path "
                "(full interval-close coverage of the shorter path's tiled bp)"
            ),
            "interval_close_spectrum_ge_floor": spectrum,
            "build_set": build,
        }
    )

    # ------------------------------------------------------------------
    # (a) the 805-class test: all 167 loci (the windowed frame's
    #     truth_pair_in_candidate_domain is 0/167 — the whole chromosome
    #     is the 805 class at chrIV scale)
    # ------------------------------------------------------------------
    assert all(
        not record["truth_pair_in_candidate_domain"]
        for record in balanced["rows"]
    ), "the balanced chrIV receipt must carry truth_pair_in_candidate_domain=false everywhere"
    axis_partitions = {partition for _, _, partition in axis}
    disease = 0
    disease_foreign = 0
    contig_end = 0
    in_axis = 0
    for locus, (start, end, partition) in enumerate(axis):
        balanced_record = balanced_rows[locus]
        assessable = balanced_record["assessable"]
        window_hits = hits.get((truth0, locus, truth1), [])
        record = {
            "question": "a",
            "component": "chrIV",
            "locus": locus,
            "window": [start, end],
            "axis_partition": partition,
            "windowed_kind": (
                "assessable"
                if assessable
                else f"bracketed:{balanced_record.get('bracket_reason')}"
            ),
            "truth_pair_in_candidate_domain": balanced_record[
                "truth_pair_in_candidate_domain"
            ],
            "sk1_window_hits": [
                {
                    "interval": [hit["start"], hit["end"]],
                    "strand": hit["strand"],
                    "anchors": anchor_count(hit),
                    "placement": [
                        {"partition": p, "row": [s, e]}
                        for p, s, e in placement(
                            by_path, TRUTH1_NAME, hit["start"], hit["end"]
                        )
                    ],
                }
                for hit in window_hits
            ],
        }
        # the POSITIONAL statement (a coordinate statement): the window's
        # coordinates on the SK1 truth contig
        positional_rows = placement(by_path, TRUTH1_NAME, start, end)
        positional_verdict, positional_covered = statement_verdict(
            positional_rows, partition, start, end
        )
        record["positional"] = {
            "rows": [
                {"partition": p, "row": [s, e]} for p, s, e in positional_rows
            ],
            "bp_covered_by_rows": positional_covered,
            "window_bp": end - start,
            "verdict": positional_verdict,
        }
        # the ORTHOLOG statement (the balanced receipt's own SK1 truth
        # interval — the alignment-truth statement) where assessable
        if assessable:
            o_start, o_end = balanced_record["SK1_truth_interval"]
            o_rows = placement(by_path, TRUTH1_NAME, o_start, o_end)
            o_verdict, o_covered = statement_verdict(o_rows, partition, o_start, o_end)
            record["ortholog"] = {
                "interval": [o_start, o_end],
                "strand": balanced_record["SK1_truth_strand"],
                "rows": [{"partition": p, "row": [s, e]} for p, s, e in o_rows],
                "bp_covered_by_rows": o_covered,
                "bp": o_end - o_start,
                "verdict": o_verdict,
            }
            verdict = o_verdict
        else:
            record["ortholog"] = None
            verdict = positional_verdict
        # the absent class: beyond the SK1 truth contig's end is a length
        # polymorphism, not a placement gap (the chrMT L13 precedent)
        if verdict == "ABSENT" or (
            verdict == "PARTIAL-ELSEWHERE" and start >= sk1_len
        ):
            absent_class = "CONTIG-END-LENGTH-POLYMORPHISM"
        elif verdict == "ABSENT":
            absent_class = "NO-ROWS"
        else:
            absent_class = "NONE"
        # the split class (the partition-spine disease decomposition):
        # NEIGHBOR-AXIS = every holder is an axis partition (a coordinate
        # offset class -- the material tiles in an adjacent axis
        # partition); FOREIGN-REPEAT = at least one holder sits outside
        # the axis partition set (the disease class, holders named)
        statement_rows = (
            (record["ortholog"] or {}).get("rows")
            if assessable
            else record["positional"]["rows"]
        )
        holders = {row["partition"] for row in statement_rows}
        if verdict in ("TILED-ELSEWHERE", "PARTIAL-ELSEWHERE"):
            split_class = (
                "NEIGHBOR-AXIS" if holders <= axis_partitions
                else "FOREIGN-REPEAT"
            )
        else:
            split_class = "NONE"
        if verdict == "IN-AXIS-PARTITION":
            in_axis += 1
        elif absent_class == "CONTIG-END-LENGTH-POLYMORPHISM":
            contig_end += 1
        else:
            disease += 1
            if split_class == "FOREIGN-REPEAT":
                disease_foreign += 1
        record["verdict"] = verdict
        record["absent_class"] = absent_class
        record["split_class"] = split_class
        positional_hits = [
            hit for hit in window_hits if hit["start"] < end and start < hit["end"]
        ]
        record["window_query"] = {
            "hits": len(window_hits),
            "positional_hit_bp": sum(
                min(hit["end"], end) - max(hit["start"], start)
                for hit in positional_hits
            ),
            "best_hit": None
            if not window_hits
            else {
                "interval": [
                    max(window_hits, key=anchor_count)["start"],
                    max(window_hits, key=anchor_count)["end"],
                ],
                "anchors": max(anchor_count(hit) for hit in window_hits),
            },
        }
        census.append(record)

    census.append(
        {
            "question": "a-summary",
            "component": "chrIV",
            "total_loci": len(axis),
            "in_axis_partition": in_axis,
            "split_elsewhere": disease,
            "split_neighbor_axis": disease - disease_foreign,
            "split_foreign_repeat": disease_foreign,
            "contig_end_length_polymorphism": contig_end,
        }
    )

    # ------------------------------------------------------------------
    # (c) the near-twin condensation (the census recipe, per flagged pair)
    # ------------------------------------------------------------------
    for pair in flagged:
        left_name, right_name = pair
        left, right = ids[left_name], ids[right_name]
        record = {
            "question": "c",
            "component": "chrIV",
            "routes": [left_name, right_name],
            "route_ids": [left, right],
            "partitions": {},
        }
        membership = {}
        for name in (left_name, right_name):
            membership[name] = [
                {"partition": p, "row": [s, e]}
                for p, s, e in sorted(by_path.get(name, ()))
            ]
        record["membership"] = membership
        left_rows = by_path.get(left_name, ())
        right_rows = by_path.get(right_name, ())
        common_partitions = sorted(
            {p for p, _, _ in left_rows} & {p for p, _, _ in right_rows}
        )
        record["common_partitions"] = common_partitions
        measured_partitions = [p for p in common_partitions if p in build]
        record["unmeasured_common_partitions"] = [
            p for p in common_partitions if p not in build
        ]
        left_nodes_total, right_nodes_total, shared_nodes_total = set(), set(), set()
        for partition in measured_partitions:
            paths = load_gfa_paths(partition)
            left_nodes = {
                abs(node)
                for name, steps in paths.items()
                if name.startswith(left_name + ":")
                for node, _ in steps
            }
            right_nodes = {
                abs(node)
                for name, steps in paths.items()
                if name.startswith(right_name + ":")
                for node, _ in steps
            }
            shared = left_nodes & right_nodes
            left_nodes_total |= left_nodes
            right_nodes_total |= right_nodes
            shared_nodes_total |= shared
            record["partitions"][str(partition)] = {
                "left_nodes": len(left_nodes),
                "right_nodes": len(right_nodes),
                "shared_nodes": len(shared),
                "left_only": sorted(left_nodes - right_nodes),
                "right_only": sorted(right_nodes - left_nodes),
            }
        record["totals"] = {
            "left_nodes": len(left_nodes_total),
            "right_nodes": len(right_nodes_total),
            "shared_nodes": len(shared_nodes_total),
        }
        # variant pockets in bp from the panel's own walks (the built
        # partitions only), the census recipe
        pockets = {"left": [], "right": []}
        for side, name, path_id, other_total in (
            ("left", left_name, left, right_nodes_total),
            ("right", right_name, right, left_nodes_total),
        ):
            for partition in common_partitions:
                for row_p, start, end in by_path.get(name, ()):
                    if row_p != partition:
                        continue
                    steps = walks.get((path_id, start, end, f"partition{partition}"))
                    if steps is None:
                        continue
                    run = None
                    for bp, node in steps:
                        if abs(node) not in other_total:
                            if run is None:
                                run = [bp, bp]
                            else:
                                run[1] = bp
                        else:
                            if run is not None:
                                pockets[side].append([run[0], run[1]])
                                run = None
                    if run is not None:
                        pockets[side].append([run[0], run[1]])
            pockets[side].sort()
        record["variant_pockets_bp"] = pockets
        census.append(record)

    # ------------------------------------------------------------------
    # (d) the truth pair's joint spellability at every locus
    # ------------------------------------------------------------------
    joint_in_axis = 0
    for locus, (start, end, partition) in enumerate(axis):
        balanced_record = balanced_rows[locus]
        assessable = balanced_record["assessable"]
        window_hits = hits.get((truth0, locus, truth1), [])
        s288c_rows = placement(by_path, TRUTH0_NAME, start, end)
        record = {
            "question": "d",
            "component": "chrIV",
            "locus": locus,
            "window": [start, end],
            "axis_partition": partition,
            "windowed_kind": (
                "assessable"
                if assessable
                else f"bracketed:{balanced_record.get('bracket_reason')}"
            ),
            "s288c_rows": [
                {"partition": p, "row": [s, e]} for p, s, e in s288c_rows
            ],
        }
        positional_rows = placement(by_path, TRUTH1_NAME, start, end)
        positional_verdict, positional_covered = statement_verdict(
            positional_rows, partition, start, end
        )
        sk1 = {
            "statement": "ortholog" if assessable else "positional",
            "window_hits": len(window_hits),
            "positional_rows": [
                {"partition": p, "row": [s, e]} for p, s, e in positional_rows
            ],
            "positional_bp_covered_by_rows": positional_covered,
            "positional_verdict": positional_verdict,
        }
        if assessable:
            o_start, o_end = balanced_record["SK1_truth_interval"]
            o_rows = placement(by_path, TRUTH1_NAME, o_start, o_end)
            o_verdict, o_covered = statement_verdict(o_rows, partition, o_start, o_end)
            sk1.update(
                {
                    "ortholog_interval": [o_start, o_end],
                    "ortholog_rows": [
                        {"partition": p, "row": [s, e]} for p, s, e in o_rows
                    ],
                    "ortholog_bp_covered_by_rows": o_covered,
                    "verdict": o_verdict,
                }
            )
            verdict = o_verdict
            statement_rows = o_rows
        else:
            verdict = positional_verdict
            statement_rows = positional_rows
        sk1["verdict"] = verdict
        # the joint verdict: S288C's window row is the axis row by
        # construction; the pair is jointly spellable in one partition
        # graph iff the SK1 statement verdict is IN-AXIS-PARTITION
        if verdict == "IN-AXIS-PARTITION":
            joint = "JOINT-IN-AXIS-PARTITION"
            joint_in_axis += 1
        elif verdict in ("TILED-ELSEWHERE", "PARTIAL-ELSEWHERE"):
            holders = {p for p, _, _ in statement_rows}
            joint = (
                "SK1-SPLIT-NEIGHBOR-AXIS-PARTITIONS"
                if holders <= axis_partitions
                else "SK1-SPLIT-FOREIGN-PARTITIONS"
            )
        else:
            joint = "SK1-ABSENT-CONTIG-END"
        foreign_holders = sorted(
            {p for p, _, _ in statement_rows} - {partition}
        )
        record["sk1"] = sk1
        record["joint_verdict"] = joint
        record["sk1_foreign_holders"] = foreign_holders
        census.append(record)

    with open(CENSUS_RECEIPT, "w") as handle:
        for record in census:
            handle.write(json.dumps(record) + "\n")
    print(f"wrote {len(census)} census records -> {CENSUS_RECEIPT}")

    # ------------------------------------------------------------------
    # headline tables
    # ------------------------------------------------------------------
    survey = next(r for r in census if r["question"] == "s")
    print("\n== (s) the chrIV survey")
    print(
        f"  {survey['partitions_genome_wide']} partitions genome-wide; "
        f"chrIV-family: {survey['chriv_family_paths']} paths / "
        f"{survey['chriv_family_rows']} rows / {survey['chriv_family_partitions']} partitions"
    )
    print(
        f"  the balanced validation's chrIV component: {survey['loci']} loci over "
        f"{survey['axis_partitions']} axis partitions; S288C#0#chrIV "
        f"{survey['s288c_rows']} rows / {survey['s288c_len_bp']} bp; SK1#0#chrIV "
        f"{survey['sk1_rows']} rows / {survey['sk1_len_bp']} bp "
        f"({survey['sk1_rows_in_axis_partitions']} rows in axis partitions; "
        f"{len(survey['sk1_holder_partitions'])} SK1 holder partitions: "
        f"{survey['sk1_holder_partitions']})"
    )
    print(
        f"  near-twin routes flagged: {['/'.join(p) for p in survey['twin_pairs']]}; "
        f"{len(survey['twin_holder_partitions'])} twin-only holder partitions: "
        f"{survey['twin_holder_partitions']}; build set {len(survey['build_set'])} partitions"
    )
    print(f"  interval-close spectrum >= {SPECTRUM_FLOOR}:")
    for entry in survey["interval_close_spectrum_ge_floor"]:
        print(
            f"    {'/'.join(entry['pair'])}: close {entry['interval_close_fraction']:.3f} "
            f"exact {entry['exact_fraction']:.3f}"
        )

    summary = next(r for r in census if r["question"] == "a-summary")
    print("\n== (a) the 805-class test at chrIV scale (all 167 loci)")
    print(
        f"  IN-AXIS-PARTITION {summary['in_axis_partition']} | "
        f"SPLIT-ELSEWHERE {summary['split_elsewhere']} "
        f"(NEIGHBOR-AXIS {summary['split_neighbor_axis']} | "
        f"FOREIGN-REPEAT {summary['split_foreign_repeat']}) | "
        f"CONTIG-END-LENGTH-POLYMORPHISM {summary['contig_end_length_polymorphism']}"
    )
    for record in census:
        if record.get("question") != "a":
            continue
        if record["split_class"] != "FOREIGN-REPEAT" and record["absent_class"] != "CONTIG-END-LENGTH-POLYMORPHISM":
            continue
        where = ", ".join(
            f"partition{row['partition']}{row['row']}"
            for row in (
                record["ortholog"]["rows"] if record["ortholog"] else record["positional"]["rows"]
            )
        )
        print(
            f"  L{record['locus']} ({record['windowed_kind']}, window {record['window']}, "
            f"axis p{record['axis_partition']}): {record['verdict']}"
            f" | {where}"
            f" | window-query hits {record['window_query']['hits']}"
        )

    print("\n== (c) the near-twin condensation")
    for record in census:
        if record.get("question") != "c":
            continue
        totals = record["totals"]
        left_pct = 100.0 * totals["shared_nodes"] / max(1, totals["left_nodes"])
        right_pct = 100.0 * totals["shared_nodes"] / max(1, totals["right_nodes"])
        print(
            f"  {'/'.join(record['routes'])}: shared {totals['shared_nodes']}/"
            f"{totals['left_nodes']} ({left_pct:.1f}%) and {totals['shared_nodes']}/"
            f"{totals['right_nodes']} ({right_pct:.1f}%) across "
            f"{len(record['common_partitions'])} common partitions "
            f"({len(record['unmeasured_common_partitions'])} unbuilt)"
        )
        print(
            f"    variant pockets (bp runs from the panel walks): "
            f"left-only {len(record['variant_pockets_bp']['left'])} runs, "
            f"right-only {len(record['variant_pockets_bp']['right'])} runs"
        )

    print("\n== (d) the truth pair's joint spellability (per locus)")
    counts = collections.Counter()
    for record in census:
        if record.get("question") == "d":
            counts[record["joint_verdict"]] += 1
    print(f"  verdict counts: {dict(counts)}; JOINT-IN-AXIS {joint_in_axis}/{len(axis)}")
    for record in census:
        if record.get("question") != "d":
            continue
        if record["joint_verdict"] in (
            "JOINT-IN-AXIS-PARTITION",
            "SK1-SPLIT-NEIGHBOR-AXIS-PARTITIONS",
        ):
            continue
        print(
            f"  L{record['locus']} ({record['windowed_kind']}, window {record['window']}, "
            f"axis p{record['axis_partition']}): {record['joint_verdict']} "
            f"| SK1 {record['sk1']['verdict']} | holders {record['sk1_foreign_holders']}"
        )


def main():
    if len(sys.argv) < 2 or sys.argv[1] not in ("requests", "census"):
        print(__doc__)
        sys.exit(1)
    if sys.argv[1] == "requests":
        write_requests()
    else:
        run_census()


if __name__ == "__main__":
    main()
