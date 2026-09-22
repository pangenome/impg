#!/usr/bin/env python3
"""Audit C4 crush POASTA replacement cycle diagnostics.

The script consumes the debug tree emitted by IMPG_CRUSH_DEBUG_DIR and writes
small TSV reports for:

* per-GFA segment/link/path counts and directed SCC/self-loop counts,
* exact path-sequence preservation for the final laced graph,
* per-applied-candidate range-preservation flags, and
* matches from applied POASTA replacements back to raw graph_to_gfa dumps.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


sys.setrecursionlimit(1_000_000)
ORIENTS = {"+", "-"}
WALK_STEP_RE = re.compile(r"([><])([^><]+)")
RC = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")


@dataclass(frozen=True)
class Link:
    from_segment: str
    from_orient: str
    to_segment: str
    to_orient: str


@dataclass
class GfaData:
    path: Path
    segments: dict[str, bytes]
    links: list[Link]
    paths: dict[str, list[tuple[str, str]]]
    record_counts: dict[str, int]


@dataclass
class CycleSummary:
    label: str
    kind: str
    path: Path
    segments: int
    links: int
    paths: int
    walks: int
    path_steps: int
    segment_bp: int
    directed_sccs_gt1: int
    directed_self_loop_edges: int
    segment_self_loop_links: int
    largest_directed_scc: int
    cyclic_oriented_nodes: int
    path_digest: str
    path_bp_total: int


def rev_orient(orient: str) -> str:
    return "-" if orient == "+" else "+"


def oriented_sequence(sequence: bytes, orient: str) -> bytes:
    if orient == "+":
        return sequence
    return sequence.translate(RC)[::-1]


def parse_path_steps(raw_steps: str) -> list[tuple[str, str]]:
    if raw_steps in {"", "*"}:
        return []
    steps: list[tuple[str, str]] = []
    for token in raw_steps.split(","):
        if not token:
            continue
        orient = token[-1]
        if orient not in ORIENTS:
            raise ValueError(f"invalid P-line step orientation in {token!r}")
        steps.append((token[:-1], orient))
    return steps


def parse_walk_steps(raw_walk: str) -> list[tuple[str, str]]:
    if raw_walk in {"", "*"}:
        return []
    steps: list[tuple[str, str]] = []
    for match in WALK_STEP_RE.finditer(raw_walk):
        orient = "+" if match.group(1) == ">" else "-"
        steps.append((match.group(2), orient))
    return steps


def walk_name(fields: list[str]) -> str:
    if len(fields) < 7:
        return "malformed_walk"
    sample, hap, seqid, start, end = fields[1:6]
    return f"{sample}#{hap}#{seqid}:{start}-{end}"


def parse_gfa(path: Path) -> GfaData:
    segments: dict[str, bytes] = {}
    links: list[Link] = []
    paths: dict[str, list[tuple[str, str]]] = {}
    record_counts: dict[str, int] = defaultdict(int)
    duplicate_counts: dict[str, int] = defaultdict(int)

    with path.open("rb") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            raw_line = raw_line.rstrip(b"\n")
            if not raw_line:
                continue
            fields = raw_line.split(b"\t")
            record = fields[0].decode("ascii", errors="replace")
            record_counts[record] += 1
            if record == "S" and len(fields) >= 3:
                segments[fields[1].decode()] = fields[2]
            elif record == "L" and len(fields) >= 5:
                links.append(
                    Link(
                        fields[1].decode(),
                        fields[2].decode(),
                        fields[3].decode(),
                        fields[4].decode(),
                    )
                )
            elif record == "P" and len(fields) >= 3:
                name = fields[1].decode()
                duplicate_counts[name] += 1
                if duplicate_counts[name] > 1:
                    name = f"{name}#{duplicate_counts[name]}"
                paths[name] = parse_path_steps(fields[2].decode())
            elif record == "W" and len(fields) >= 7:
                text_fields = [field.decode() for field in fields]
                name = walk_name(text_fields)
                duplicate_counts[name] += 1
                if duplicate_counts[name] > 1:
                    name = f"{name}#{duplicate_counts[name]}"
                paths[name] = parse_walk_steps(text_fields[6])
            elif record in {"P", "W"}:
                raise ValueError(f"{path}:{line_number}: malformed {record} record")

    return GfaData(path=path, segments=segments, links=links, paths=paths, record_counts=dict(record_counts))


def path_sequence(data: GfaData, steps: list[tuple[str, str]]) -> bytes:
    parts: list[bytes] = []
    for segment, orient in steps:
        try:
            sequence = data.segments[segment]
        except KeyError as exc:
            raise KeyError(f"{data.path}: path references missing segment {segment!r}") from exc
        parts.append(oriented_sequence(sequence, orient))
    return b"".join(parts)


def path_sequences(data: GfaData) -> dict[str, bytes]:
    return {name: path_sequence(data, steps) for name, steps in data.paths.items()}


def path_digest_and_bp(data: GfaData) -> tuple[str, int]:
    digest = hashlib.sha256()
    total_bp = 0
    for name, sequence in sorted(path_sequences(data).items()):
        total_bp += len(sequence)
        digest.update(name.encode())
        digest.update(b"\0")
        digest.update(hashlib.sha256(sequence).digest())
        digest.update(b"\0")
    return digest.hexdigest(), total_bp


def directed_edges(data: GfaData) -> list[tuple[tuple[str, str], tuple[str, str], Link]]:
    edges: list[tuple[tuple[str, str], tuple[str, str], Link]] = []
    for link in data.links:
        if link.from_orient not in ORIENTS or link.to_orient not in ORIENTS:
            continue
        forward = ((link.from_segment, link.from_orient), (link.to_segment, link.to_orient), link)
        reverse = (
            (link.to_segment, rev_orient(link.to_orient)),
            (link.from_segment, rev_orient(link.from_orient)),
            link,
        )
        edges.append(forward)
        edges.append(reverse)
    return edges


def tarjan_scc(nodes: Iterable[tuple[str, str]], adjacency: dict[tuple[str, str], list[tuple[str, str]]]) -> list[list[tuple[str, str]]]:
    index = 0
    stack: list[tuple[str, str]] = []
    on_stack: set[tuple[str, str]] = set()
    indices: dict[tuple[str, str], int] = {}
    lowlinks: dict[tuple[str, str], int] = {}
    components: list[list[tuple[str, str]]] = []

    def visit(node: tuple[str, str]) -> None:
        nonlocal index
        indices[node] = index
        lowlinks[node] = index
        index += 1
        stack.append(node)
        on_stack.add(node)

        for target in adjacency.get(node, []):
            if target not in indices:
                visit(target)
                lowlinks[node] = min(lowlinks[node], lowlinks[target])
            elif target in on_stack:
                lowlinks[node] = min(lowlinks[node], indices[target])

        if lowlinks[node] == indices[node]:
            component: list[tuple[str, str]] = []
            while True:
                item = stack.pop()
                on_stack.remove(item)
                component.append(item)
                if item == node:
                    break
            components.append(component)

    for node in nodes:
        if node not in indices:
            visit(node)
    return components


def summarize_cycles(label: str, kind: str, path: Path) -> CycleSummary:
    data = parse_gfa(path)
    adjacency: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    nodes: set[tuple[str, str]] = set()
    self_loop_edges = 0
    segment_self_loop_links = 0

    for segment in data.segments:
        nodes.add((segment, "+"))
        nodes.add((segment, "-"))

    for source, target, link in directed_edges(data):
        nodes.add(source)
        nodes.add(target)
        adjacency[source].append(target)
        if source == target:
            self_loop_edges += 1
        if link.from_segment == link.to_segment:
            segment_self_loop_links += 1

    components = tarjan_scc(nodes, adjacency)
    sccs_gt1 = [component for component in components if len(component) > 1]
    self_loop_nodes = {source for source, target, _ in directed_edges(data) if source == target}
    cyclic_nodes = set(self_loop_nodes)
    for component in sccs_gt1:
        cyclic_nodes.update(component)
    path_digest, path_bp_total = path_digest_and_bp(data)

    return CycleSummary(
        label=label,
        kind=kind,
        path=path,
        segments=len(data.segments),
        links=len(data.links),
        paths=data.record_counts.get("P", 0),
        walks=data.record_counts.get("W", 0),
        path_steps=sum(len(steps) for steps in data.paths.values()),
        segment_bp=sum(len(sequence) for sequence in data.segments.values()),
        directed_sccs_gt1=len(sccs_gt1),
        directed_self_loop_edges=self_loop_edges,
        segment_self_loop_links=segment_self_loop_links,
        largest_directed_scc=max((len(component) for component in sccs_gt1), default=0),
        cyclic_oriented_nodes=len(cyclic_nodes),
        path_digest=path_digest,
        path_bp_total=path_bp_total,
    )


def compare_paths(reference: Path, observed: Path) -> dict[str, object]:
    ref_sequences = path_sequences(parse_gfa(reference))
    obs_sequences = path_sequences(parse_gfa(observed))
    ref_names = set(ref_sequences)
    obs_names = set(obs_sequences)
    common = ref_names & obs_names
    mismatched = [name for name in sorted(common) if ref_sequences[name] != obs_sequences[name]]
    return {
        "reference": reference,
        "observed": observed,
        "reference_paths": len(ref_sequences),
        "observed_paths": len(obs_sequences),
        "common_paths": len(common),
        "missing_paths": len(ref_names - obs_names),
        "extra_paths": len(obs_names - ref_names),
        "mismatched_paths": len(mismatched),
        "first_mismatches": ",".join(mismatched[:10]),
    }


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def find_frontier_dirs(debug_dir: Path) -> list[Path]:
    return sorted(path for path in debug_dir.glob("applied_frontier_*") if path.is_dir())


def candidate_rows(frontier: Path) -> list[dict[str, str]]:
    metadata = frontier / "applied-candidates.tsv"
    if not metadata.exists():
        return []
    rows = read_tsv(metadata)
    for row in rows:
        row["_frontier"] = frontier.name
        row["_frontier_path"] = str(frontier)
    return rows


def replacement_path(frontier: Path, row: dict[str, str]) -> Path:
    rel = row.get("replacement_gfa") or f"{row['label']}/replacement.gfa"
    return frontier / rel


def collect_summaries(debug_dir: Path, final_gfa: Path | None) -> tuple[list[CycleSummary], list[dict[str, str]]]:
    summaries: list[CycleSummary] = []
    candidates: list[dict[str, str]] = []
    for frontier in find_frontier_dirs(debug_dir):
        final_laced = frontier / "final_laced.gfa"
        if final_laced.exists():
            summaries.append(summarize_cycles(f"{frontier.name}/final_laced", "final_laced_debug", final_laced))
        for row in candidate_rows(frontier):
            path = replacement_path(frontier, row)
            if not path.exists():
                continue
            candidates.append(row)
            summaries.append(summarize_cycles(row["label"], "applied_replacement", path))

    if final_gfa is not None and final_gfa.exists():
        summaries.append(summarize_cycles("final_output", "final_output", final_gfa))

    for build_dir in sorted(debug_dir.glob("replacement_*_poasta")):
        raw = build_dir / "raw.graph_to_gfa.gfa"
        exact = build_dir / "exact.normalized.gfa"
        if raw.exists():
            summaries.append(summarize_cycles(build_dir.name, "poasta_raw_graph_to_gfa", raw))
        if exact.exists():
            summaries.append(summarize_cycles(build_dir.name, "poasta_exact_normalized", exact))
    return summaries, candidates


def cycle_summary_rows(summaries: list[CycleSummary], root: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for summary in summaries:
        try:
            rel = summary.path.relative_to(root)
        except ValueError:
            rel = summary.path
        rows.append(
            {
                "label": summary.label,
                "kind": summary.kind,
                "gfa": rel,
                "segments": summary.segments,
                "links": summary.links,
                "paths": summary.paths,
                "walks": summary.walks,
                "path_steps": summary.path_steps,
                "segment_bp": summary.segment_bp,
                "path_bp_total": summary.path_bp_total,
                "directed_sccs_gt1": summary.directed_sccs_gt1,
                "directed_self_loop_edges": summary.directed_self_loop_edges,
                "segment_self_loop_links": summary.segment_self_loop_links,
                "largest_directed_scc": summary.largest_directed_scc,
                "cyclic_oriented_nodes": summary.cyclic_oriented_nodes,
                "path_digest": summary.path_digest,
            }
        )
    return rows


def candidate_preservation_rows(debug_dir: Path) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for lengths_path in sorted(debug_dir.glob("applied_frontier_*/*/path-lengths.tsv")):
        label = lengths_path.parent.name
        lengths = read_tsv(lengths_path)
        failed = [
            row
            for row in lengths
            if row.get("preserved") != "true" or row.get("observed_bp", "").startswith("error:") or row.get("observed_bp") == "missing"
        ]
        rows.append(
            {
                "label": label,
                "paths_checked": len(lengths),
                "paths_failed": len(failed),
                "first_failed_paths": ",".join(row.get("replacement_path", "") for row in failed[:10]),
                "path_lengths_tsv": lengths_path.relative_to(debug_dir),
            }
        )
    return rows


def poasta_build_match_rows(debug_dir: Path, candidates: list[dict[str, str]]) -> list[dict[str, object]]:
    digest_to_builds: dict[str, list[tuple[Path, CycleSummary, CycleSummary | None]]] = defaultdict(list)
    for build_dir in sorted(debug_dir.glob("replacement_*_poasta")):
        exact = build_dir / "exact.normalized.gfa"
        raw = build_dir / "raw.graph_to_gfa.gfa"
        if not exact.exists():
            continue
        exact_summary = summarize_cycles(build_dir.name, "poasta_exact_normalized", exact)
        raw_summary = summarize_cycles(build_dir.name, "poasta_raw_graph_to_gfa", raw) if raw.exists() else None
        digest_to_builds[exact_summary.path_digest].append((build_dir, exact_summary, raw_summary))

    rows: list[dict[str, object]] = []
    for row in candidates:
        if row.get("method") != "Poasta":
            continue
        frontier = Path(row["_frontier_path"])
        applied = replacement_path(frontier, row)
        if not applied.exists():
            continue
        applied_summary = summarize_cycles(row["label"], "applied_replacement", applied)
        matches = digest_to_builds.get(applied_summary.path_digest, [])
        if not matches:
            rows.append(
                {
                    "label": row["label"],
                    "matched_build_dir": "",
                    "matched": "false",
                    "applied_directed_sccs_gt1": applied_summary.directed_sccs_gt1,
                    "raw_directed_sccs_gt1": "",
                    "exact_directed_sccs_gt1": "",
                    "applied_self_loops": applied_summary.directed_self_loop_edges,
                    "raw_self_loops": "",
                    "exact_self_loops": "",
                    "raw_gfa": "",
                    "exact_gfa": "",
                }
            )
            continue
        for build_dir, exact_summary, raw_summary in matches:
            rows.append(
                {
                    "label": row["label"],
                    "matched_build_dir": build_dir.name,
                    "matched": "true",
                    "applied_directed_sccs_gt1": applied_summary.directed_sccs_gt1,
                    "raw_directed_sccs_gt1": raw_summary.directed_sccs_gt1 if raw_summary else "",
                    "exact_directed_sccs_gt1": exact_summary.directed_sccs_gt1,
                    "applied_self_loops": applied_summary.directed_self_loop_edges,
                    "raw_self_loops": raw_summary.directed_self_loop_edges if raw_summary else "",
                    "exact_self_loops": exact_summary.directed_self_loop_edges,
                    "raw_gfa": (build_dir / "raw.graph_to_gfa.gfa").relative_to(debug_dir),
                    "exact_gfa": (build_dir / "exact.normalized.gfa").relative_to(debug_dir),
                }
            )
    return rows


def max_numeric_segment_id(data: GfaData) -> int:
    max_id = 0
    for segment in data.segments:
        try:
            max_id = max(max_id, int(segment))
        except ValueError:
            continue
    return max_id


def candidate_id_ranges(input_gfa: Path, candidates: list[dict[str, str]]) -> list[tuple[int, int, str, str]]:
    next_id = max_numeric_segment_id(parse_gfa(input_gfa)) + 1
    ranges: list[tuple[int, int, str, str]] = []
    for row in sorted(candidates, key=lambda item: int(item.get("rank", "0") or "0")):
        segment_count = int(row.get("replacement_segments", "0") or "0")
        if segment_count == 0:
            continue
        start = next_id
        end = next_id + segment_count - 1
        ranges.append((start, end, row.get("label", ""), row.get("method", "")))
        next_id = end + 1
    return ranges


def segment_origin(segment: str, ranges: list[tuple[int, int, str, str]], input_max_id: int) -> tuple[str, str]:
    try:
        segment_id = int(segment)
    except ValueError:
        return ("non_numeric", segment)
    if segment_id <= input_max_id:
        return ("input_seed", "input_seed")
    for start, end, label, method in ranges:
        if start <= segment_id <= end:
            return (method or "replacement", label)
    return ("unmapped_replacement", "unmapped_replacement")


def final_cycle_origin_rows(input_gfa: Path, final_gfa: Path, candidates: list[dict[str, str]]) -> list[dict[str, object]]:
    input_data = parse_gfa(input_gfa)
    final_data = parse_gfa(final_gfa)
    input_max_id = max_numeric_segment_id(input_data)
    ranges = candidate_id_ranges(input_gfa, candidates)

    adjacency: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
    nodes = {(segment, orient) for segment in final_data.segments for orient in ORIENTS}
    self_edges: list[tuple[tuple[str, str], Link]] = []
    for source, target, link in directed_edges(final_data):
        nodes.add(source)
        nodes.add(target)
        adjacency[source].append(target)
        if source == target:
            self_edges.append((source, link))

    summary: dict[tuple[str, str, str], dict[str, object]] = {}

    def add_cycle(cycle_type: str, component: list[tuple[str, str]]) -> None:
        class_names: set[str] = set()
        labels: set[str] = set()
        for segment, _orient in component:
            class_name, label = segment_origin(segment, ranges, input_max_id)
            class_names.add(class_name)
            labels.add(label)
        origin_classes = ",".join(sorted(class_names))
        origin_labels = ",".join(sorted(labels))
        key = (cycle_type, origin_classes, origin_labels)
        example = ",".join(f"{segment}{orient}" for segment, orient in sorted(component)[:30])
        if key not in summary:
            summary[key] = {
                "cycle_type": cycle_type,
                "origin_classes": origin_classes,
                "origin_labels": origin_labels,
                "count": 0,
                "max_size": 0,
                "example_nodes": example,
            }
        summary[key]["count"] = int(summary[key]["count"]) + 1
        summary[key]["max_size"] = max(int(summary[key]["max_size"]), len(component))

    for source, _link in self_edges:
        add_cycle("directed_self_loop", [source])
    for component in tarjan_scc(nodes, adjacency):
        if len(component) > 1:
            add_cycle("directed_scc", component)

    return sorted(
        summary.values(),
        key=lambda row: (str(row["cycle_type"]), str(row["origin_classes"]), str(row["origin_labels"])),
    )


def write_cycle_examples(debug_dir: Path, output_dir: Path, max_examples: int) -> None:
    rows: list[dict[str, object]] = []
    gfas = [path for path in debug_dir.glob("applied_frontier_*/*/replacement.gfa")]
    gfas.extend(path for path in debug_dir.glob("applied_frontier_*/final_laced.gfa"))
    gfas.extend(path for path in debug_dir.glob("replacement_*_poasta/raw.graph_to_gfa.gfa"))
    gfas.extend(path for path in debug_dir.glob("replacement_*_poasta/exact.normalized.gfa"))
    for gfa in sorted(gfas):
        data = parse_gfa(gfa)
        adjacency: dict[tuple[str, str], list[tuple[str, str]]] = defaultdict(list)
        nodes = {(segment, orient) for segment in data.segments for orient in ORIENTS}
        self_edges: list[tuple[tuple[str, str], tuple[str, str], Link]] = []
        for source, target, link in directed_edges(data):
            adjacency[source].append(target)
            nodes.add(source)
            nodes.add(target)
            if source == target:
                self_edges.append((source, target, link))
        components = [component for component in tarjan_scc(nodes, adjacency) if len(component) > 1]
        label = gfa.parent.name if gfa.name == "replacement.gfa" else gfa.parent.name + "/" + gfa.name
        for idx, (source, _, link) in enumerate(self_edges[:max_examples], start=1):
            rows.append(
                {
                    "gfa": gfa.relative_to(debug_dir),
                    "label": label,
                    "example_type": "directed_self_loop",
                    "index": idx,
                    "size": 1,
                    "nodes": f"{source[0]}{source[1]}",
                    "link": f"{link.from_segment}{link.from_orient}->{link.to_segment}{link.to_orient}",
                }
            )
        for idx, component in enumerate(components[:max_examples], start=1):
            rows.append(
                {
                    "gfa": gfa.relative_to(debug_dir),
                    "label": label,
                    "example_type": "directed_scc",
                    "index": idx,
                    "size": len(component),
                    "nodes": ",".join(f"{segment}{orient}" for segment, orient in sorted(component)[:50]),
                    "link": "",
                }
            )
    write_tsv(
        output_dir / "cycle-examples.tsv",
        rows,
        ["gfa", "label", "example_type", "index", "size", "nodes", "link"],
    )


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--debug-dir", type=Path, required=True)
    parser.add_argument("--input-gfa", type=Path, required=True)
    parser.add_argument("--final-gfa", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--max-examples", type=int, default=5)
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)

    summaries, candidates = collect_summaries(args.debug_dir, args.final_gfa)
    write_tsv(
        args.output_dir / "gfa-cycle-summary.tsv",
        cycle_summary_rows(summaries, args.debug_dir),
        [
            "label",
            "kind",
            "gfa",
            "segments",
            "links",
            "paths",
            "walks",
            "path_steps",
            "segment_bp",
            "path_bp_total",
            "directed_sccs_gt1",
            "directed_self_loop_edges",
            "segment_self_loop_links",
            "largest_directed_scc",
            "cyclic_oriented_nodes",
            "path_digest",
        ],
    )

    final_debug = None
    frontier_dirs = find_frontier_dirs(args.debug_dir)
    if frontier_dirs:
        candidate_final = frontier_dirs[-1] / "final_laced.gfa"
        if candidate_final.exists():
            final_debug = candidate_final
    preservation_rows = [compare_paths(args.input_gfa, args.final_gfa)]
    if final_debug is not None and final_debug != args.final_gfa:
        preservation_rows.append(compare_paths(args.input_gfa, final_debug))
    write_tsv(
        args.output_dir / "path-preservation.tsv",
        preservation_rows,
        [
            "reference",
            "observed",
            "reference_paths",
            "observed_paths",
            "common_paths",
            "missing_paths",
            "extra_paths",
            "mismatched_paths",
            "first_mismatches",
        ],
    )

    write_tsv(
        args.output_dir / "candidate-path-preservation.tsv",
        candidate_preservation_rows(args.debug_dir),
        ["label", "paths_checked", "paths_failed", "first_failed_paths", "path_lengths_tsv"],
    )

    write_tsv(
        args.output_dir / "poasta-build-matches.tsv",
        poasta_build_match_rows(args.debug_dir, candidates),
        [
            "label",
            "matched_build_dir",
            "matched",
            "applied_directed_sccs_gt1",
            "raw_directed_sccs_gt1",
            "exact_directed_sccs_gt1",
            "applied_self_loops",
            "raw_self_loops",
            "exact_self_loops",
            "raw_gfa",
            "exact_gfa",
        ],
    )

    write_tsv(
        args.output_dir / "final-cycle-origins.tsv",
        final_cycle_origin_rows(args.input_gfa, args.final_gfa, candidates),
        ["cycle_type", "origin_classes", "origin_labels", "count", "max_size", "example_nodes"],
    )

    write_cycle_examples(args.debug_dir, args.output_dir, args.max_examples)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
