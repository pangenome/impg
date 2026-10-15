#!/usr/bin/env python3
"""Independent receipt checker for THE PARTITION-EMBEDDED GRAPH BUILD and
THE EXPRESSIBILITY CENSUS (the owner's ruling, 2026-11-06).

Validates, from the receipts on disk and nothing else (every re-derivation
below is an independent implementation — the checker does not import or
execute the census code):

  phase 1 — run markers: the graph build and the alignment/walk dumps
            exited 0 with done markers; all 48 GFAs and mapping sidecars
            exist;
  phase 2 — GFA structure: S/L/P counts match the mapping sidecars, gap
            segment ids are contiguous from n_syncmer_nodes+1, every S id
            is in range, every P step is an emitted S node, and the P
            names are exactly the partition's member rows
            (`name:start-end`);
  phase 3 — request re-derivation: both request files (homology, walks)
            are re-derived from the axis, the partition BEDs, the
            committed likelihood receipts, and the multi-census receipts,
            and must match the committed request files exactly;
  phase 4 — receipt validation: every homology request has a completion
            marker with matching extents, every emitted hit targets a
            requested path, and every walk request has exactly one walk
            receipt line;
  phase 5 — SPELL EQUALITY (the id-continuity proof): for every walk
            receipt, the corresponding GFA P line spells EXACTLY the
            panel's own walk — same signed global syncmer node ids in the
            same order (gap segments, ids > n_syncmer_nodes, are skipped):
            the partition graph's coordinates ARE the global syng
            interning, measured, not assumed;
  phase 6 — verdict re-derivation: the census's (a) 805-class verdicts,
            (b) pocket placement verdicts, (c) near-twin shared-node
            structure, and (d) per-locus truth-spell verdicts are
            re-derived from the raw receipts and must match the census
            receipt record by record.

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
SCRATCH = f"{DATA}/cosine-diagnostic-scratch"
PADDING = 120
COMPONENTS = ("chrMT", "chrI")
TRUTH_PATHS = {"chrI": (9579, 9613), "chrMT": (9580, 9614)}
TRUTH_NAMES = {
    "chrI": ("S288C#0#chrI", "SK1#0#chrI"),
    "chrMT": ("S288C#0#chrMT", "SK1#0#chrMT"),
}
TWINS = {"chrMT": (5599, 5616), "chrI": (5397, 5545)}
TAIL_LOCI = {"chrI": [7, 8, 13, 15, 16, 18]}
BUILT = sorted(
    [p for p in range(21)]
    + [14900, 6382, 14901, 14902, 15582, 15968, 14903, 14904,
       14905, 14906, 15586, 16122, 16242, 16360]
    + [16239, 15966, 15972]
    + [102, 16771, 14610, 3134, 18304, 14346, 14482, 14606, 3093, 1899]
    + [15970, 16121]
)

failures = []


def check(condition, message):
    if condition:
        print(f"ok: {message}")
    else:
        failures.append(message)
        print(f"FAIL: {message}")


def coverage(rows, start, end):
    """Independent union-coverage of [start,end) by the rows (a checker-local
    implementation; the census has its own)."""
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


def load_names():
    names = {}
    with open(NAMES) as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            names[int(fields[0])] = fields[1]
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
    axis = {}
    for interval in json.load(open(AXIS_FILE))["intervals"]:
        component = interval["component"].rsplit("#", 1)[1]
        if component in COMPONENTS:
            axis.setdefault(component, []).append(
                (interval["start"], interval["end"], int(interval["group"][9:]))
            )
    return axis


def derive_pocket_carriers(locus, names):
    """Independent re-derivation of the pocket carriers at a chrI tail
    window: off-truth occurrences sharing covered nodes with the truth
    pair's window occurrences, on contigs outside the chrI family."""
    truth = set(TRUTH_PATHS["chrI"])
    carriers = collections.Counter()
    with open(f"{SCRATCH}/chrI/cosine-multi-census.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            truth_nodes = set()
            for occurrence in record["occurrences"]:
                if occurrence["path"] in truth and locus in occurrence["partitions"]:
                    for interval in occurrence["intervals"]:
                        truth_nodes.update(interval)
            if not truth_nodes:
                continue
            for occurrence in record["occurrences"]:
                if occurrence["path"] in truth:
                    continue
                contig = names[occurrence["path"]].rsplit("#", 1)[1]
                if contig.startswith("chrI"):
                    continue
                nodes = set()
                for interval in occurrence["intervals"]:
                    nodes.update(interval)
                if nodes & truth_nodes:
                    carriers[names[occurrence["path"]]] += 1
    return carriers


def overlapping(rows_by_path, name, start, end):
    return sorted(
        (p, s, e)
        for p, s, e in rows_by_path.get(name, ())
        if s < end and start < e
    )


def phase1_markers():
    print("== phase 1: run markers")
    for marker in ("build-partition-graphs.done", "dumps.done"):
        check(os.path.exists(f"{GRAPHS}/{marker}"), f"marker {marker} exists")
    for receipt in ("build-partition-graphs.log", "dumps.log"):
        check(os.path.exists(f"{GRAPHS}/{receipt}"), f"log {receipt} exists")
    for exit_file in ("partition-graph-homology.exit", "partition-graph-walks.exit"):
        code = open(f"{GRAPHS}/{exit_file}").read().strip()
        check(code == "0", f"{exit_file} == 0")
    missing = [
        p for p in BUILT
        if not os.path.exists(f"{GRAPHS}/partition{p}.gfa")
        or not os.path.exists(f"{GRAPHS}/partition{p}.gfa.map.json")
    ]
    check(not missing, f"all {len(BUILT)} GFAs + maps present (missing: {missing})")


def phase2_gfa_structure(by_partition):
    print("== phase 2: GFA structure")
    n_nodes = json.load(open(f"{GRAPHS}/partition{BUILT[0]}.gfa.map.json"))["n_syncmer_nodes"]
    check(
        all(
            json.load(open(f"{GRAPHS}/partition{p}.gfa.map.json"))["n_syncmer_nodes"] == n_nodes
            for p in BUILT
        ),
        "n_syncmer_nodes identical across all maps (the global interning)",
    )
    problems = []
    for partition in BUILT:
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
        if gap_ids and gap_ids != list(
            range(n_nodes + 1, n_nodes + 1 + len(gap_ids))
        ):
            problems.append(f"p{partition}: gap ids not contiguous")
        for name, steps in path_lines.items():
            if any(abs(node) not in segments for node, _ in steps):
                problems.append(f"p{partition}: P {name} steps not all S")
        if any(not (0 < i <= n_nodes + counts["gap_segments"]) for i in segments):
            problems.append(f"p{partition}: S id out of range")
    check(not problems, f"GFA structure clean on all {len(BUILT)} (problems: {problems[:6]})")
    return n_nodes


def phase3_requests(names, axis, by_partition, by_path):
    print("== phase 3: request re-derivation")
    ids = {name: path for path, name in names.items()}
    homology_lines = []
    walk_lines = []
    for component in COMPONENTS:
        truth0, truth1 = TRUTH_PATHS[component]
        for locus, (start, end, partition) in enumerate(axis[component]):
            targets = {truth1}
            if component == "chrI" and locus in TAIL_LOCI["chrI"]:
                targets.update(ids[n] for n in derive_pocket_carriers(locus, names))
            for twin in TWINS[component]:
                if any(p == partition for p, _, _ in by_path.get(names[twin], ())):
                    targets.add(twin)
            homology_lines.append(
                f"{truth0}\t{locus}\t{start}\t{end}\t{PADDING}\t"
                + ",".join(str(t) for t in sorted(targets))
            )
    walk_names = {n for c in COMPONENTS for n in TRUTH_NAMES[c]} | {
        names[t] for c in COMPONENTS for t in TWINS[c]
    }
    for partition in BUILT:
        for name, start, end in by_partition[partition]:
            if name in walk_names:
                walk_lines.append(f"{ids[name]}\t{start}\t{end}\tpartition{partition}")
    committed_homology = open(f"{GRAPHS}/partition-graph-homology-requests.tsv").read()
    committed_walks = open(f"{GRAPHS}/partition-graph-walk-requests.tsv").read()
    check(
        committed_homology == "\n".join(homology_lines) + "\n",
        f"homology requests re-derived exactly ({len(homology_lines)} lines)",
    )
    check(
        committed_walks == "\n".join(walk_lines) + "\n",
        f"walk requests re-derived exactly ({len(walk_lines)} lines)",
    )
    return homology_lines, walk_lines


def phase4_receipts(homology_lines, walk_lines):
    print("== phase 4: receipt validation")
    markers = {}
    hits = collections.defaultdict(list)
    with open(f"{GRAPHS}/partition-graph-homology.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            if record.get("request_done"):
                markers[(record["query_path"], record["window"], record["qs"], record["qe"])] = record["hits"]
            else:
                hits[(record["query_path"], record["window"], record["target_path"])].append(record)
    marker_problems = []
    for request in homology_lines:
        query, window, qs, qe, padding, targets = request.split("\t")
        key = (int(query), int(window), int(qs), int(qe))
        if key not in markers:
            marker_problems.append(f"missing marker {key}")
    check(not marker_problems, f"every homology request has a completion marker ({marker_problems[:3]})")
    target_problems = []
    for (query, window, target), target_hits in hits.items():
        request = next(
            (r for r in homology_lines
             if r.split("\t")[0] == str(query) and r.split("\t")[1] == str(window)),
            None,
        )
        if request is None or str(target) not in request.split("\t")[5].split(","):
            target_problems.append((query, window, target))
    check(not target_problems, f"every emitted hit targets a requested path ({target_problems[:3]})")
    walks = {}
    with open(f"{GRAPHS}/partition-graph-walks.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            walks[(record["path"], record["start"], record["end"], record["tag"])] = record["steps"]
    walk_keys = {
        (int(f[0]), int(f[1]), int(f[2]), f[3]) for f in
        (line.split("\t") for line in walk_lines)
    }
    check(
        set(walks) == walk_keys,
        f"walk receipts == walk requests exactly ({len(walks)} walks)",
    )
    return hits, walks


def phase5_spells(walks, n_nodes, names):
    print("== phase 5: SPELL EQUALITY (id continuity)")
    problems = []
    spelled = 0
    for (path, start, end, tag), steps in walks.items():
        partition = int(tag[len("partition") :])
        wanted = f"{names[path]}:{start}-{end}"
        gfa_steps = None
        with open(f"{GRAPHS}/partition{partition}.gfa") as handle:
            for line in handle:
                if line.startswith(f"P\t{wanted}\t"):
                    fields = line.rstrip("\n").split("\t")
                    gfa_steps = [
                        (int(step[:-1]), step[-1] == "+") for step in fields[2].split(",")
                    ]
                    break
        if gfa_steps is None:
            problems.append(f"{tag} path {path} [{start},{end}): no P line for {wanted}")
            continue
        # skip gap segments (ids beyond the syncmer interning)
        syncmer_steps = [(node, forward) for node, forward in gfa_steps if abs(node) <= n_nodes]
        # the walk emits (bp, node) with the node's storage sign; the GFA
        # step sign is forward when the id is positive
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


def phase6_verdicts(names, axis, by_partition, by_path, hits):
    print("== phase 6: verdict re-derivation")
    census = {}
    with open(f"{GRAPHS}/partition-graph-census.jsonl") as handle:
        for line in handle:
            record = json.loads(line)
            key = record.get("question")
            if key == "a":
                census[("a", record["component"], record["locus"])] = record
            elif key == "b":
                census[("b", record["locus"])] = record
            elif key == "c":
                census[("c", record["component"])] = record
            elif key == "d":
                census[("d", record["component"], record["locus"])] = record
            elif key == "a-summary":
                census[("a-summary",)] = record

    # (a)
    problems = []
    present = 0
    non_expressible = []
    status_kinds = {}
    for component in COMPONENTS:
        with open(f"{DATA}/cosine-graph-likelihood-readmatched-{component}.jsonl") as handle:
            for line in handle:
                record = json.loads(line)
                pieces = record["truth_piece_presence"]
                if record["truth_pair_expressible"]:
                    kind = "expressible"
                elif all(pieces):
                    kind = "bracketed"
                else:
                    kind = "non-expressible"
                status_kinds[(component, record["locus"])] = kind
    for component in COMPONENTS:
        truth0, truth1 = TRUTH_PATHS[component]
        for locus, (start, end, partition) in enumerate(axis[component]):
            if status_kinds[(component, locus)] != "non-expressible":
                continue
            non_expressible.append((component, locus))
            rows = overlapping(by_path, TRUTH_NAMES[component][1], start, end)
            covered = coverage(rows, start, end)
            if any(p == partition for p, _, _ in rows):
                verdict = "IN-AXIS-PARTITION"
            elif covered >= end - start:
                verdict = "TILED-ELSEWHERE"
            elif rows:
                verdict = "PARTIAL-ELSEWHERE"
            else:
                verdict = "ABSENT"
            record = census.get(("a", component, locus))
            if record is None:
                problems.append(f"(a) {component} L{locus}: no census record")
                continue
            if record["verdict"] != verdict:
                problems.append(f"(a) {component} L{locus}: {record['verdict']} != {verdict}")
            elif record["positional_bp_covered_by_rows"] != covered:
                problems.append(f"(a) {component} L{locus}: coverage mismatch")
            else:
                if verdict != "ABSENT":
                    present += 1
    summary = census[("a-summary",)]
    check(
        summary["present_by_alignment"] == present
        and len(summary["non_expressible_loci"]) == len(non_expressible),
        f"(a) headline re-derived: {present}/{len(non_expressible)} present by alignment",
    )
    check(not problems, f"(a) verdicts match the census receipt ({problems[:4]})")

    # (b)
    problems = []
    ids = {name: path for path, name in names.items()}
    truth0 = TRUTH_PATHS["chrI"][0]
    for locus in TAIL_LOCI["chrI"]:
        _, _, partition = axis["chrI"][locus]
        carriers = derive_pocket_carriers(locus, names)
        record = census[("b", locus)]
        census_carriers = {p["carrier"]: p for p in record["pockets"]}
        if set(census_carriers) != set(carriers):
            problems.append(f"(b) L{locus}: carrier sets differ")
        for carrier, _ in carriers.items():
            pocket = census_carriers.get(carrier)
            if pocket is None:
                continue
            window_hits = hits.get((truth0, locus, ids[carrier]), [])
            if window_hits:
                start, end = window_hits[0]["start"], window_hits[0]["end"]
                rows = overlapping(by_path, carrier, start, end)
                verdict = (
                    "IN-BY-ALIGNMENT" if any(p == partition for p, _, _ in rows)
                    else "OUT-ALIGNED-ELSEWHERE" if rows
                    else "OUT-UNPLACED"
                )
            else:
                verdict = "NO-ALIGNMENT-INTO-WINDOW"
            if pocket["verdict"] != verdict:
                problems.append(f"(b) L{locus} {carrier}: {pocket['verdict']} != {verdict}")
    check(not problems, f"(b) pocket verdicts match the census receipt ({problems[:4]})")

    # (c)
    problems = []
    for component in COMPONENTS:
        left, right = TWINS[component]
        left_name, right_name = names[left], names[right]
        left_rows = by_path.get(left_name, ())
        right_rows = by_path.get(right_name, ())
        common = sorted({p for p, _, _ in left_rows} & {p for p, _, _ in right_rows})
        # the shared-node structure is measured on the BUILT graphs only
        # (the census states the same restriction)
        measured = [p for p in common if p in BUILT]
        record = census[("c", component)]
        if sorted(record["common_partitions"]) != common:
            problems.append(f"(c) {component}: common partition lists differ")
        if record["unmeasured_common_partitions"] != [
            p for p in common if p not in BUILT
        ]:
            problems.append(f"(c) {component}: unmeasured list differs")
        left_total, right_total, shared_total = set(), set(), set()
        for partition in measured:
            segments, path_lines = None, {}
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
            problems.append(f"(c) {component}: totals differ")
    check(not problems, f"(c) near-twin shared-node totals match ({problems})")

    # (d)
    problems = []
    for component in COMPONENTS:
        truth1_name = TRUTH_NAMES[component][1]
        for locus, (start, end, partition) in enumerate(axis[component]):
            record = census[("d", component, locus)]
            rows = overlapping(by_path, truth1_name, start, end)
            covered = coverage(rows, start, end)
            if any(p == partition for p, _, _ in rows):
                verdict = "IN-AXIS-PARTITION"
            elif covered >= end - start:
                verdict = "TILED-ELSEWHERE"
            elif rows:
                verdict = "PARTIAL-ELSEWHERE"
            else:
                verdict = "ABSENT"
            if record["sk1"]["verdict"] != verdict:
                problems.append(
                    f"(d) {component} L{locus}: {record['sk1']['verdict']} != {verdict}"
                )
            elif record["sk1"]["positional_bp_covered_by_rows"] != covered:
                problems.append(f"(d) {component} L{locus}: coverage mismatch")
    check(not problems, f"(d) truth-spell verdicts match the census receipt ({problems[:4]})")


def main():
    names = load_names()
    by_partition, by_path = scan_beds()
    axis = load_axis()
    phase1_markers()
    n_nodes = phase2_gfa_structure(by_partition)
    homology_lines, walk_lines = phase3_requests(names, axis, by_partition, by_path)
    hits, walks = phase4_receipts(homology_lines, walk_lines)
    spelled = phase5_spells(walks, n_nodes, names)
    phase6_verdicts(names, axis, by_partition, by_path, hits)
    print("\n== summary")
    if failures:
        print(f"{len(failures)} FAILURES:")
        for failure in failures:
            print(f"  {failure}")
        sys.exit(1)
    print("ALL PHASES PASS")


if __name__ == "__main__":
    main()
