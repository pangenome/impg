#!/usr/bin/env python3
"""Independent receipt checker for THE CHRIV PARTITION-EMBEDDED GRAPH
BUILD AND THE EXPRESSIBILITY CENSUS AT CHRIV SCALE (chrIV end-to-end
slice 1; the chrMT/chrI checker's recipe extended to the first 1.5 Mb
chromosome).

Validates, from the receipts on disk and nothing else (every
re-derivation below is an independent implementation — the checker does
not import or execute the census code):

  phase 1 — run markers: the graph build and the alignment/walk dumps
            exited 0 with done markers, RSS peaks under the 64GiB
            guard, and all build-set GFAs and mapping sidecars exist;
  phase 2 — GFA structure: S/L/P counts match the mapping sidecars, gap
            segment ids are contiguous from n_syncmer_nodes+1, every S
            id is in range, every P step is an emitted S node, and the
            P names are exactly the partition's member rows;
  phase 3 — request re-derivation: the build list and both request
            files (homology, walks) are re-derived from the axis, the
            partition BEDs, the balanced chrIV receipt, and the stated
            near-twin flag rule (full interval-close coverage, both
            ends within 50 bp), and must match the committed files
            exactly;
  phase 4 — receipt validation: every homology request has a completion
            marker with matching extents, every emitted hit targets a
            requested path, and every walk request has exactly one
            walk receipt line;
  phase 5 — SPELL EQUALITY (the id-continuity proof): for every walk
            receipt, the corresponding GFA P line spells EXACTLY the
            panel's own walk — same signed global syncmer node ids in
            the same order (gap segments, ids > n_syncmer_nodes, are
            skipped): the partition graph's coordinates ARE the global
            syng interning, measured, not assumed;
  phase 6 — verdict re-derivation: the survey counts, the (a) 805-class
            verdicts (ortholog + positional statements, the absent
            classes), the (c) near-twin flag spectrum + shared-node
            structure + variant-pocket counts, and the (d) per-locus
            joint-spellability verdicts are re-derived from the raw
            receipts and must match the census receipt record by
            record.

Assessment-side only. No product file is read or written.
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
PADDING = 120
INTERVAL_CLOSE = 50
SPECTRUM_FLOOR = 0.3
GUARD_KB = 64 * 1024 * 1024
TRUTH0_NAME = "S288C#0#chrIV"
TRUTH1_NAME = "SK1#0#chrIV"

failures = []


def check(condition, message):
    if condition:
        print(f"ok: {message}")
    else:
        failures.append(message)
        print(f"FAIL: {message}")


def load_names():
    names = {}
    with open(NAMES) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            names[int(fields[0])] = (fields[1], int(fields[2]))
    return names


def scan_beds():
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
                if len(fields) >= 3:
                    rows.append((fields[0], int(fields[1]), int(fields[2])))
        by_partition[partition] = rows
        for name, start, end in rows:
            by_path[name].append((partition, start, end))
    return by_partition, by_path


def load_axis():
    axis = []
    for interval in json.load(open(AXIS_FILE))["intervals"]:
        if interval["component"] == TRUTH0_NAME:
            axis.append(
                (interval["start"], interval["end"], int(interval["group"][len("partition") :]))
            )
    axis.sort()
    return axis


def coverage(rows, start, end):
    """Independent union-coverage of [start,end) by the rows."""
    covered = 0
    last = start
    for row_start, row_end in sorted((s, e) for _, s, e in rows):
        if row_end <= last:
            continue
        covered += min(row_end, end) - max(row_start, last)
        last = max(last, row_end)
        if last >= end:
            break
    return covered


def overlapping(rows_by_path, name, start, end):
    return sorted(
        (p, s, e) for p, s, e in rows_by_path.get(name, ()) if s < end and start < e
    )


def close_bp(rows_by_path, name_a, name_b):
    """Independent interval-close coverage of a's tiled bp by b's rows
    in the same partitions (both ends within INTERVAL_CLOSE bp)."""
    by_part = collections.defaultdict(list)
    for partition, start, end in rows_by_path.get(name_b, ()):
        by_part[partition].append((start, end))
    matched = exact = total = 0
    for partition, start, end in rows_by_path.get(name_a, ()):
        total += end - start
        for start2, end2 in by_part.get(partition, ()):
            if abs(start2 - start) <= INTERVAL_CLOSE and abs(end2 - end) <= INTERVAL_CLOSE:
                matched += end - start
                if start2 == start and end2 == end:
                    exact += end - start
                break
    return matched, exact, total


def derive_twins(by_path):
    """Independent re-derivation of the near-twin flag rule and the
    reported spectrum (the checker's own implementation)."""
    family = sorted(name for name in by_path if name.rsplit("#", 1)[1] == "chrIV")
    spectrum = []
    for i, name_a in enumerate(family):
        for name_b in family[i + 1:]:
            ab, _, total_a = close_bp(by_path, name_a, name_b)
            ba, _, total_b = close_bp(by_path, name_b, name_a)
            if not total_a or not total_b:
                continue
            frac = min(ab / total_a, ba / total_b)
            if frac >= SPECTRUM_FLOOR:
                _, exact_ab, _ = close_bp(by_path, name_a, name_b)
                spectrum.append(
                    {
                        "pair": [name_a, name_b],
                        "interval_close_fraction": round(frac, 4),
                        "exact_fraction": round(exact_ab / total_a, 4),
                    }
                )
    spectrum.sort(key=lambda entry: -entry["interval_close_fraction"])
    flagged = [
        tuple(entry["pair"]) for entry in spectrum if entry["interval_close_fraction"] >= 1.0
    ]
    return flagged, spectrum, family


def derive_build(axis, by_path, twin_names):
    partitions = {p for _, _, p in axis}
    for name in [TRUTH1_NAME] + list(twin_names):
        partitions |= {p for p, _, _ in by_path.get(name, ())}
    return sorted(partitions)


def verdict_of(rows, partition, start, end):
    covered = coverage(rows, start, end)
    if any(p == partition for p, _, _ in rows):
        return "IN-AXIS-PARTITION", covered
    if not rows:
        return "ABSENT", covered
    if covered >= end - start:
        return "TILED-ELSEWHERE", covered
    return "PARTIAL-ELSEWHERE", covered


def phase1_markers(build):
    print("== phase 1: run markers")
    for name in (
        "build-chrIV-partition-graphs",
        "dumps-chrIV-homology",
        "dumps-chrIV-walks",
    ):
        check(os.path.exists(f"{GRAPHS}/{name}.done"), f"marker {name}.done exists")
        code = open(f"{GRAPHS}/{name}.exit").read().strip()
        check(code == "0", f"{name}.exit == 0")
        check(os.path.exists(f"{GRAPHS}/{name}.log"), f"log {name}.log exists")
        wall = open(f"{GRAPHS}/{name}.wall").read().strip()
        print(f"   {name}: wall {wall}s")
    check(os.path.exists(f"{GRAPHS}/dumps-chrIV.done"), "marker dumps-chrIV.done exists")
    peaks = {}
    for name in (
        "build-chrIV-partition-graphs",
        "dumps-chrIV-homology",
        "dumps-chrIV-walks",
    ):
        peak = 0
        path = f"{GRAPHS}/{name}.rss"
        if os.path.exists(path):
            with open(path) as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) >= 3 and fields[2].isdigit():
                        peak = max(peak, int(fields[2]))
        peaks[name] = peak
        check(peak < GUARD_KB, f"{name}: RSS peak {peak}kB under the 64GiB guard")
    missing = [
        p for p in build
        if not os.path.exists(f"{GRAPHS}/partition{p}.gfa")
        or not os.path.exists(f"{GRAPHS}/partition{p}.gfa.map.json")
    ]
    check(not missing, f"all {len(build)} GFAs + maps present (missing: {missing})")
    return peaks


def phase2_gfa_structure(by_partition, build):
    print("== phase 2: GFA structure")
    n_nodes = json.load(open(f"{GRAPHS}/partition{build[0]}.gfa.map.json"))["n_syncmer_nodes"]
    check(
        all(
            json.load(open(f"{GRAPHS}/partition{p}.gfa.map.json"))["n_syncmer_nodes"] == n_nodes
            for p in build
        ),
        "n_syncmer_nodes identical across all maps (the global interning)",
    )
    problems = []
    for partition in build:
        mapping = json.load(open(f"{GRAPHS}/partition{partition}.gfa.map.json"))
        segments = {}
        path_lines = {}
        links = 0
        with open(f"{GRAPHS}/partition{partition}.gfa") as handle:
            for line in handle:
                if line.startswith("S\t"):
                    fields = line.split("\t")
                    segments[int(fields[1])] = fields[2].strip()
                elif line.startswith("L\t"):
                    links += 1
                elif line.startswith("P\t"):
                    fields = line.rstrip("\n").split("\t")
                    path_lines[fields[1]] = [
                        (int(step[:-1]), step[-1]) for step in fields[2].split(",")
                    ]
        counts = mapping["counts"]
        expected_names = {
            f"{name}:{start}-{end}" for name, start, end in by_partition[partition]
        }
        gap_ids = sorted(i for i in segments if i > n_nodes)
        if len(segments) != counts["segments"]:
            problems.append(f"p{partition}: S count")
        if links != counts["links"]:
            problems.append(f"p{partition}: L count")
        if len(path_lines) != counts["members"]:
            problems.append(f"p{partition}: P count != members")
        if set(path_lines) != expected_names:
            problems.append(f"p{partition}: P names != member rows")
        if len(gap_ids) != counts["gap_segments"]:
            problems.append(f"p{partition}: gap count")
        if gap_ids and gap_ids != list(range(n_nodes + 1, n_nodes + 1 + len(gap_ids))):
            problems.append(f"p{partition}: gap ids not contiguous")
        for name, steps in path_lines.items():
            if any(abs(node) not in segments for node, _ in steps):
                problems.append(f"p{partition}: P {name} steps not all S")
        if any(not (0 < i <= n_nodes + counts["gap_segments"]) for i in segments):
            problems.append(f"p{partition}: S id out of range")
    check(not problems, f"GFA structure clean on all {len(build)} (problems: {problems[:6]})")
    return n_nodes


def phase3_requests(names, axis, by_path, build):
    print("== phase 3: request re-derivation")
    ids = {name: path for path, (name, _) in names.items()}
    flagged, spectrum, _ = derive_twins(by_path)
    twin_names = sorted({name for pair in flagged for name in pair})
    check(
        [len(flagged), [list(p) for p in flagged]] == [2, [["AAA#0#chrIV", "SGDref#0#chrIV"], ["CLL#0#chrIV", "CLL#1#chrIV"]]],
        f"near-twin flag re-derived: {['/'.join(p) for p in flagged]}",
    )
    derived_build = derive_build(axis, by_path, twin_names)
    check(
        derived_build == build,
        f"build list re-derived exactly ({len(build)} partitions)",
    )
    homology_lines = []
    for locus, (start, end, partition) in enumerate(axis):
        targets = {ids[TRUTH1_NAME]}
        for name in twin_names:
            if any(p == partition for p, _, _ in by_path.get(name, ())):
                targets.add(ids[name])
        homology_lines.append(
            f"{ids[TRUTH0_NAME]}\t{locus}\t{start}\t{end}\t{PADDING}\t"
            + ",".join(str(t) for t in sorted(targets))
        )
    walk_names = {TRUTH0_NAME, TRUTH1_NAME} | set(twin_names)
    walk_lines = []
    for partition in build:
        with open(f"{BEDS}/partition{partition}.bed") as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 3:
                    continue
                name, start, end = fields[0], int(fields[1]), int(fields[2])
                if name in walk_names:
                    walk_lines.append(f"{ids[name]}\t{start}\t{end}\tpartition{partition}")
    committed_build = [
        int(line)
        for line in open(f"{GRAPHS}/partition-graph-chrIV-build-list.txt").read().split()
    ]
    committed_homology = open(f"{GRAPHS}/partition-graph-chrIV-homology-requests.tsv").read()
    committed_walks = open(f"{GRAPHS}/partition-graph-chrIV-walk-requests.tsv").read()
    check(
        committed_build == derived_build,
        f"committed build list == re-derived ({len(committed_build)} partitions)",
    )
    check(
        committed_homology == "\n".join(homology_lines) + "\n",
        f"homology requests re-derived exactly ({len(homology_lines)} lines)",
    )
    check(
        committed_walks == "\n".join(walk_lines) + "\n",
        f"walk requests re-derived exactly ({len(walk_lines)} lines)",
    )
    return homology_lines, walk_lines, flagged, spectrum


def phase4_receipts(homology_lines, walk_lines):
    print("== phase 4: receipt validation")
    markers = {}
    hits = collections.defaultdict(list)
    with open(f"{GRAPHS}/partition-graph-chrIV-homology.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("request_done"):
                markers[(record["query_path"], record["window"], record["qs"], record["qe"])] = record["hits"]
            else:
                hits[(record["query_path"], record["window"], record["target_path"])].append(record)
    marker_problems = []
    for request in homology_lines:
        query, window, qs, qe, padding, targets = request.split("\t")
        if (int(query), int(window), int(qs), int(qe)) not in markers:
            marker_problems.append(request)
    check(not marker_problems, f"every homology request has a completion marker ({marker_problems[:2]})")
    target_problems = []
    for (query, window, target), _ in hits.items():
        request = next(
            (r for r in homology_lines
             if r.split("\t")[0] == str(query) and r.split("\t")[1] == str(window)),
            None,
        )
        if request is None or str(target) not in request.split("\t")[5].split(","):
            target_problems.append((query, window, target))
    check(not target_problems, f"every emitted hit targets a requested path ({target_problems[:3]})")
    walks = {}
    with open(f"{GRAPHS}/partition-graph-chrIV-walks.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            walks[(record["path"], record["start"], record["end"], record["tag"])] = record["steps"]
    walk_keys = {
        (int(f[0]), int(f[1]), int(f[2]), f[3])
        for f in (line.split("\t") for line in walk_lines)
    }
    check(set(walks) == walk_keys, f"walk receipts == walk requests exactly ({len(walks)} walks)")
    return hits, walks


def phase5_spells(walks, n_nodes, names):
    print("== phase 5: SPELL EQUALITY (id continuity)")
    # index the P lines of every built partition once
    problems = []
    spelled = 0
    by_partition_steps = collections.defaultdict(dict)
    wanted = collections.defaultdict(set)
    for (path, start, end, tag) in walks:
        partition = int(tag[len("partition") :])
        wanted[partition].add(f"{names[path][0]}:{start}-{end}")
    for partition, names_wanted in wanted.items():
        with open(f"{GRAPHS}/partition{partition}.gfa") as handle:
            for line in handle:
                if not line.startswith("P\t"):
                    continue
                fields = line.rstrip("\n").split("\t")
                if fields[1] in names_wanted:
                    by_partition_steps[partition][fields[1]] = [
                        (int(step[:-1]), step[-1] == "+") for step in fields[2].split(",")
                    ]
    for (path, start, end, tag), steps in walks.items():
        partition = int(tag[len("partition") :])
        wanted_name = f"{names[path][0]}:{start}-{end}"
        gfa_steps = by_partition_steps[partition].get(wanted_name)
        if gfa_steps is None:
            problems.append(f"{tag} path {path} [{start},{end}): no P line")
            continue
        syncmer_steps = [(node, forward) for node, forward in gfa_steps if abs(node) <= n_nodes]
        walk_steps = [(abs(node), node > 0) for _, node in steps]
        if syncmer_steps != walk_steps:
            problems.append(
                f"{tag} path {path} [{start},{end}): spell mismatch "
                f"({len(syncmer_steps)} gfa vs {len(walk_steps)} walk steps)"
            )
        else:
            spelled += 1
    check(not problems, f"every walk is spelled exactly by its GFA P line ({problems[:4]})")
    check(spelled == len(walks), f"all {len(walks)} walks verified ({spelled})")
    return spelled


def phase6_verdicts(names, axis, by_path, build, hits, flagged, spectrum):
    print("== phase 6: verdict re-derivation")
    ids = {name: path for path, (name, _) in names.items()}
    sk1_len = names[ids[TRUTH1_NAME]][1]
    balanced = json.load(open(BALANCED))
    balanced_rows = {record["locus"]: record for record in balanced["rows"]}
    census = {}
    with open(f"{GRAPHS}/partition-graph-chrIV-census.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            key = record.get("question")
            if key == "a":
                census[("a", record["locus"])] = record
            elif key == "c":
                census[("c", tuple(record["routes"]))] = record
            elif key == "d":
                census[("d", record["locus"])] = record
            elif key == "s":
                census[("s",)] = record
            elif key == "a-summary":
                census[("a-summary",)] = record

    # survey record
    problems = []
    survey = census[("s",)]
    family = sorted(name for name in by_path if name.rsplit("#", 1)[1] == "chrIV")
    check(
        survey["partitions_genome_wide"] == 19421,
        f"survey: {survey['partitions_genome_wide']} partitions genome-wide",
    )
    check(
        survey["chriv_family_paths"] == len(family)
        and survey["chriv_family_rows"] == sum(len(by_path.get(n, ())) for n in family),
        "survey: chrIV-family path/row counts re-derived",
    )
    axis_partitions = {p for _, _, p in axis}
    sk1_rows = by_path.get(TRUTH1_NAME, ())
    check(
        survey["loci"] == len(axis) == 167
        and survey["axis_partitions"] == len(axis_partitions)
        and survey["s288c_rows"] == len(by_path.get(TRUTH0_NAME, ()))
        and survey["sk1_rows"] == len(sk1_rows)
        and survey["sk1_len_bp"] == sk1_len
        and survey["sk1_rows_in_axis_partitions"]
        == len([r for r in sk1_rows if r[0] in axis_partitions])
        and survey["sk1_holder_partitions"]
        == sorted({p for p, _, _ in sk1_rows} - axis_partitions),
        "survey: locus/axis/truth-side counts re-derived",
    )
    twin_names = sorted({name for pair in flagged for name in pair})
    check(
        [list(p) for p in flagged] == survey["twin_pairs"]
        and spectrum == survey["interval_close_spectrum_ge_floor"],
        "survey: the near-twin flag + spectrum re-derived exactly",
    )
    check(
        survey["build_set"] == build,
        "survey: build set re-derived",
    )

    # (a) per-locus verdicts
    axis_partitions = {p for _, _, p in axis}
    in_axis = split_foreign = contig_end = split_total = 0
    for locus, (start, end, partition) in enumerate(axis):
        record = census[("a", locus)]
        balanced_record = balanced_rows[locus]
        assessable = balanced_record["assessable"]
        positional_rows = overlapping(by_path, TRUTH1_NAME, start, end)
        positional_verdict, positional_covered = verdict_of(positional_rows, partition, start, end)
        if assessable:
            o_start, o_end = balanced_record["SK1_truth_interval"]
            o_rows = overlapping(by_path, TRUTH1_NAME, o_start, o_end)
            verdict, covered = verdict_of(o_rows, partition, o_start, o_end)
            statement_rows = o_rows
        else:
            verdict, covered = positional_verdict, positional_covered
            statement_rows = positional_rows
        if verdict == "ABSENT" or (verdict == "PARTIAL-ELSEWHERE" and start >= sk1_len):
            absent_class = "CONTIG-END-LENGTH-POLYMORPHISM"
        elif verdict == "ABSENT":
            absent_class = "NO-ROWS"
        else:
            absent_class = "NONE"
        holders = {p for p, _, _ in statement_rows}
        if verdict in ("TILED-ELSEWHERE", "PARTIAL-ELSEWHERE"):
            split_class = "NEIGHBOR-AXIS" if holders <= axis_partitions else "FOREIGN-REPEAT"
        else:
            split_class = "NONE"
        if record["verdict"] != verdict or record["absent_class"] != absent_class:
            problems.append(
                f"(a) L{locus}: {record['verdict']}/{record['absent_class']} != {verdict}/{absent_class}"
            )
        if record["split_class"] != split_class:
            problems.append(f"(a) L{locus}: split class {record['split_class']} != {split_class}")
        if record["positional"]["verdict"] != positional_verdict:
            problems.append(f"(a) L{locus}: positional verdict mismatch")
        if record["positional"]["bp_covered_by_rows"] != positional_covered:
            problems.append(f"(a) L{locus}: positional coverage mismatch")
        if assessable:
            if record["ortholog"] is None:
                problems.append(f"(a) L{locus}: ortholog statement missing")
            elif record["ortholog"]["verdict"] != verdict or record["ortholog"]["bp_covered_by_rows"] != covered:
                problems.append(f"(a) L{locus}: ortholog verdict/coverage mismatch")
        if verdict == "IN-AXIS-PARTITION":
            in_axis += 1
        elif absent_class == "CONTIG-END-LENGTH-POLYMORPHISM":
            contig_end += 1
        else:
            split_total += 1
            if split_class == "FOREIGN-REPEAT":
                split_foreign += 1
    check(not problems, f"(a) per-locus verdicts match the census receipt ({problems[:4]})")
    summary = census[("a-summary",)]
    check(
        summary["in_axis_partition"] == in_axis
        and summary["split_elsewhere"] == split_total
        and summary["split_neighbor_axis"] == split_total - split_foreign
        and summary["split_foreign_repeat"] == split_foreign
        and summary["contig_end_length_polymorphism"] == contig_end,
        f"(a) headline re-derived: IN-AXIS {in_axis} / SPLIT {split_total} "
        f"(NEIGHBOR-AXIS {split_total - split_foreign} / FOREIGN-REPEAT {split_foreign}) "
        f"/ CONTIG-END {contig_end} of {len(axis)}",
    )

    # (c) the near-twin condensation
    problems = []
    for pair in flagged:
        left_name, right_name = pair
        left_rows = by_path.get(left_name, ())
        right_rows = by_path.get(right_name, ())
        common = sorted({p for p, _, _ in left_rows} & {p for p, _, _ in right_rows})
        measured = [p for p in common if p in build]
        record = census[("c", pair)]
        if sorted(record["common_partitions"]) != common:
            problems.append(f"(c) {pair}: common partition lists differ")
        left_total, right_total, shared_total = set(), set(), set()
        for partition in measured:
            path_lines = {}
            with open(f"{GRAPHS}/partition{partition}.gfa") as handle:
                for line in handle:
                    if line.startswith("P\t"):
                        fields = line.rstrip("\n").split("\t")
                        path_lines[fields[1]] = [
                            (int(step[:-1]), step[-1] == "+") for step in fields[2].split(",")
                        ]
            left_nodes = {
                abs(node)
                for name, steps in path_lines.items()
                if name.startswith(left_name + ":")
                for node, _ in steps
            }
            right_nodes = {
                abs(node)
                for name, steps in path_lines.items()
                if name.startswith(right_name + ":")
                for node, _ in steps
            }
            left_total |= left_nodes
            right_total |= right_nodes
            shared_total |= left_nodes & right_nodes
        if record["totals"] != {
            "left_nodes": len(left_total),
            "right_nodes": len(right_total),
            "shared_nodes": len(shared_total),
        }:
            problems.append(f"(c) {pair}: totals differ")
        # variant pockets re-derivation: bp runs on nodes the other twin
        # does not share anywhere in the built common partitions —
        # re-derive the pocket counts from the walks receipt
        walks = {}
        with open(f"{GRAPHS}/partition-graph-chrIV-walks.jsonl") as handle:
            for line in handle:
                w = json.loads(line)
                walks[(w["path"], w["start"], w["end"], w["tag"])] = w["steps"]
        pockets = {"left": 0, "right": 0}
        for side, name, path_id, other_total in (
            ("left", left_name, ids[left_name], right_total),
            ("right", right_name, ids[right_name], left_total),
        ):
            for partition in common:
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
                                pockets[side] += 1
                                run = None
                    if run is not None:
                        pockets[side] += 1
        if (
            len(record["variant_pockets_bp"]["left"]) != pockets["left"]
            or len(record["variant_pockets_bp"]["right"]) != pockets["right"]
        ):
            problems.append(
                f"(c) {pair}: pocket run counts differ "
                f"({len(record['variant_pockets_bp']['left'])}/{len(record['variant_pockets_bp']['right'])} "
                f"vs {pockets['left']}/{pockets['right']})"
            )
    check(not problems, f"(c) near-twin condensation totals match ({problems[:3]})")

    # (d) per-locus joint spellability
    problems = []
    joint_in_axis = 0
    for locus, (start, end, partition) in enumerate(axis):
        record = census[("d", locus)]
        balanced_record = balanced_rows[locus]
        assessable = balanced_record["assessable"]
        positional_rows = overlapping(by_path, TRUTH1_NAME, start, end)
        positional_verdict, positional_covered = verdict_of(positional_rows, partition, start, end)
        if assessable:
            o_start, o_end = balanced_record["SK1_truth_interval"]
            o_rows = overlapping(by_path, TRUTH1_NAME, o_start, o_end)
            verdict = verdict_of(o_rows, partition, o_start, o_end)[0]
            statement_rows = o_rows
        else:
            verdict = positional_verdict
            statement_rows = positional_rows
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
        foreign_holders = sorted({p for p, _, _ in statement_rows} - {partition})
        if (
            record["joint_verdict"] != joint
            or record["sk1_foreign_holders"] != foreign_holders
            or record["sk1"]["verdict"] != verdict
        ):
            problems.append(
                f"(d) L{locus}: {record['joint_verdict']}/{record['sk1']['verdict']} != {joint}/{verdict}"
            )
    check(not problems, f"(d) joint-spellability verdicts match ({problems[:4]})")
    print(f"   (d) JOINT-IN-AXIS re-derived: {joint_in_axis}/{len(axis)}")


def main():
    names = load_names()
    by_partition, by_path = scan_beds()
    axis = load_axis()
    flagged, spectrum, _ = derive_twins(by_path)
    twin_names = sorted({name for pair in flagged for name in pair})
    build = derive_build(axis, by_path, twin_names)
    phase1_markers(build)
    n_nodes = phase2_gfa_structure(by_partition, build)
    homology_lines, walk_lines, flagged, spectrum = phase3_requests(names, axis, by_path, build)
    hits, walks = phase4_receipts(homology_lines, walk_lines)
    spelled = phase5_spells(walks, n_nodes, names)
    phase6_verdicts(names, axis, by_path, build, hits, flagged, spectrum)
    print("\n== summary")
    if failures:
        print(f"{len(failures)} FAILURES:")
        for failure in failures:
            print(f"  {failure}")
        sys.exit(1)
    print("ALL PHASES PASS")


if __name__ == "__main__":
    main()
