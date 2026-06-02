#!/usr/bin/env python3
"""
Report whether a GFA looks like a focused, clean local variation graph.

The checks are deliberately topology/layout based:

- one dominant connected component
- no excessive long-range links in the current GFA segment order
- no excessive long jumps between consecutive path steps
- path starts/ends mostly converge on common head/tail nodes
- tips are mostly explained as path endpoints

This is not a biological truth check. It is a deterministic smoke test for
obvious graph construction failures such as repeat-glue, scattered long-range
links, fragmented path starts, unsupported tips, and disconnected debris.
"""

from __future__ import annotations

import argparse
import collections
import dataclasses
import math
import re
import sys
from pathlib import Path


STEP_RE = re.compile(r"(.+)([+-])$")
W_STEP_RE = re.compile(r"([<>])([^<>]+)")


@dataclasses.dataclass
class ParsedGfa:
    order: dict[str, int]
    seq_counts: collections.Counter[str]
    edges: list[tuple[str, str]]
    paths: list[list[str]]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("gfa", nargs="+", help="GFA file(s) to check")
    parser.add_argument("--tsv", type=Path, help="Write metrics TSV")
    parser.add_argument("--quiet", action="store_true", help="Only print TSV rows")
    parser.add_argument(
        "--max-link-jump-p99",
        type=int,
        default=5_000,
        help="Fail if p99 link span in segment-order units exceeds N",
    )
    parser.add_argument(
        "--max-path-jump-p99",
        type=int,
        default=5_000,
        help="Fail if p99 consecutive path-step span in segment-order units exceeds N",
    )
    parser.add_argument(
        "--max-link-jump-frac",
        type=float,
        default=0.25,
        help="Fail if max link span exceeds this fraction of segment count",
    )
    parser.add_argument(
        "--min-largest-component-frac",
        type=float,
        default=0.98,
        help="Fail if the largest undirected component contains less than this fraction of nodes",
    )
    parser.add_argument(
        "--max-components",
        type=int,
        default=1,
        help="Fail if connected component count exceeds N",
    )
    parser.add_argument(
        "--max-internal-tips",
        type=int,
        default=0,
        help="Fail if more than N degree-0/1 nodes are not path endpoints",
    )
    parser.add_argument(
        "--min-common-start-frac",
        type=float,
        default=0.70,
        help="Fail if fewer than this fraction of paths share the most common start node",
    )
    parser.add_argument(
        "--min-common-end-frac",
        type=float,
        default=0.70,
        help="Fail if fewer than this fraction of paths share the most common end node",
    )
    parser.add_argument(
        "--warn-duplicate-sequence-frac",
        type=float,
        default=0.10,
        help="Warn if more than this fraction of segment sequences are duplicated",
    )
    parser.add_argument(
        "--repeat-context-max-minor",
        type=int,
        default=2,
        help="Count a repeated-node context as loop-like if non-dominant context support is at most N",
    )
    parser.add_argument(
        "--repeat-context-min-dominance",
        type=float,
        default=0.80,
        help="Count a repeated-node context as loop-like if the dominant context support fraction is at least F",
    )
    return parser.parse_args()


def split_p_steps(raw: str) -> list[str]:
    if raw == "*":
        return []
    nodes: list[str] = []
    for token in raw.split(","):
        match = STEP_RE.fullmatch(token)
        if match:
            nodes.append(match.group(1))
    return nodes


def split_w_steps(raw: str) -> list[str]:
    if raw == "*":
        return []
    return [match.group(2) for match in W_STEP_RE.finditer(raw)]


def parse_gfa(path: Path) -> ParsedGfa:
    order: dict[str, int] = {}
    seq_counts: collections.Counter[str] = collections.Counter()
    edges: list[tuple[str, str]] = []
    paths: list[list[str]] = []

    with path.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            tag = fields[0]
            if tag == "S" and len(fields) >= 3:
                node = fields[1]
                if node not in order:
                    order[node] = len(order)
                seq_counts[fields[2]] += 1
            elif tag == "L" and len(fields) >= 5:
                edges.append((fields[1], fields[3]))
            elif tag == "P" and len(fields) >= 3:
                steps = split_p_steps(fields[2])
                if steps:
                    paths.append(steps)
            elif tag == "W" and len(fields) >= 7:
                steps = split_w_steps(fields[6])
                if steps:
                    paths.append(steps)

    return ParsedGfa(order=order, seq_counts=seq_counts, edges=edges, paths=paths)


def quantile(values: list[int], frac: float) -> int:
    if not values:
        return 0
    idx = math.ceil(frac * len(values)) - 1
    idx = min(max(idx, 0), len(values) - 1)
    return values[idx]


def component_sizes(nodes: list[str], edges: list[tuple[str, str]]) -> list[int]:
    parent = {node: node for node in nodes}
    size = {node: 1 for node in nodes}

    def find(x: str) -> str:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    def union(a: str, b: str) -> None:
        if a not in parent or b not in parent:
            return
        ra = find(a)
        rb = find(b)
        if ra == rb:
            return
        if size[ra] < size[rb]:
            ra, rb = rb, ra
        parent[rb] = ra
        size[ra] += size[rb]

    for a, b in edges:
        union(a, b)

    counts: collections.Counter[str] = collections.Counter(find(node) for node in nodes)
    return sorted(counts.values(), reverse=True)


def jump_stats(order: dict[str, int], pairs: list[tuple[str, str]]) -> tuple[int, int, int]:
    jumps = sorted(abs(order[a] - order[b]) for a, b in pairs if a in order and b in order)
    return quantile(jumps, 0.95), quantile(jumps, 0.99), (jumps[-1] if jumps else 0)


def path_jump_pairs(paths: list[list[str]]) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    for steps in paths:
        pairs.extend(zip(steps, steps[1:]))
    return pairs


def endpoint_stats(paths: list[list[str]]) -> tuple[collections.Counter[str], collections.Counter[str]]:
    starts: collections.Counter[str] = collections.Counter()
    ends: collections.Counter[str] = collections.Counter()
    for steps in paths:
        if steps:
            starts[steps[0]] += 1
            ends[steps[-1]] += 1
    return starts, ends


def top_fraction(counter: collections.Counter[str], total: int) -> tuple[str, int, float]:
    if total == 0 or not counter:
        return "", 0, 0.0
    node, count = counter.most_common(1)[0]
    return node, count, count / total


def local_repeat_context_stats(
    paths: list[list[str]],
    max_minor: int,
    min_dominance: float,
) -> tuple[int, int]:
    """Detect near-single-copy nodes with rare repeated contexts in one path.

    This is the signature of single-syncmer repeat loops: a node is usually seen
    in one local context, but a small number of paths contain the same node a
    second time with another predecessor/successor pair.
    """
    if max_minor <= 0 or not paths:
        return 0, 0

    context_counts: dict[str, collections.Counter[tuple[str, str]]] = collections.defaultdict(
        collections.Counter
    )
    for steps in paths:
        for idx, node in enumerate(steps):
            left = steps[idx - 1] if idx else "^"
            right = steps[idx + 1] if idx + 1 < len(steps) else "$"
            context_counts[node][(left, right)] += 1

    dominant: dict[str, tuple[str, str]] = {}
    for node, counts in context_counts.items():
        if len(counts) <= 1:
            continue
        ranked = counts.most_common()
        if len(ranked) > 1 and ranked[0][1] == ranked[1][1]:
            continue
        total = sum(counts.values())
        major_context, major_count = ranked[0]
        minor = total - major_count
        if minor <= max_minor and major_count / total >= min_dominance:
            dominant[node] = major_context

    clone_nodes: set[str] = set()
    clone_occurrences = 0
    for steps in paths:
        path_counts = collections.Counter(steps)
        for idx, node in enumerate(steps):
            if path_counts[node] <= 1 or node not in dominant:
                continue
            left = steps[idx - 1] if idx else "^"
            right = steps[idx + 1] if idx + 1 < len(steps) else "$"
            if (left, right) != dominant[node]:
                clone_nodes.add(node)
                clone_occurrences += 1

    return len(clone_nodes), clone_occurrences


def evaluate(path: Path, args: argparse.Namespace) -> dict[str, object]:
    gfa = parse_gfa(path)
    nodes = list(gfa.order)
    n_nodes = len(nodes)
    degree: collections.Counter[str] = collections.Counter()
    for a, b in gfa.edges:
        degree[a] += 1
        degree[b] += 1

    comps = component_sizes(nodes, gfa.edges)
    largest_component = comps[0] if comps else 0
    largest_component_frac = largest_component / n_nodes if n_nodes else 0.0

    starts, ends = endpoint_stats(gfa.paths)
    endpoint_nodes = set(starts) | set(ends)
    tips = [node for node in nodes if degree[node] <= 1]
    internal_tips = [node for node in tips if node not in endpoint_nodes]

    link_p95, link_p99, link_max = jump_stats(gfa.order, gfa.edges)
    path_p95, path_p99, path_max = jump_stats(gfa.order, path_jump_pairs(gfa.paths))

    start_node, start_count, start_frac = top_fraction(starts, len(gfa.paths))
    end_node, end_count, end_frac = top_fraction(ends, len(gfa.paths))

    duplicate_sequences = sum(1 for count in gfa.seq_counts.values() if count > 1)
    max_duplicate_count = max(gfa.seq_counts.values(), default=0)
    duplicate_sequence_frac = duplicate_sequences / n_nodes if n_nodes else 0.0
    local_repeat_context_nodes, local_repeat_context_occurrences = local_repeat_context_stats(
        gfa.paths,
        args.repeat_context_max_minor,
        args.repeat_context_min_dominance,
    )

    failures: list[str] = []
    warnings: list[str] = []
    if len(comps) > args.max_components:
        failures.append(f"components>{args.max_components}")
    if largest_component_frac < args.min_largest_component_frac:
        failures.append("largest_component_frac")
    if len(internal_tips) > args.max_internal_tips:
        failures.append(f"internal_tips>{args.max_internal_tips}")
    if start_frac < args.min_common_start_frac:
        failures.append("common_start_frac")
    if end_frac < args.min_common_end_frac:
        failures.append("common_end_frac")
    if link_p99 > args.max_link_jump_p99:
        failures.append("link_jump_p99")
    if path_p99 > args.max_path_jump_p99:
        failures.append("path_jump_p99")
    if n_nodes and link_max > int(n_nodes * args.max_link_jump_frac):
        failures.append("link_jump_max_frac")
    if duplicate_sequence_frac > args.warn_duplicate_sequence_frac:
        warnings.append("duplicate_sequence_frac")

    return {
        "file": str(path),
        "status": "FAIL" if failures else "PASS",
        "failures": ",".join(failures),
        "warnings": ",".join(warnings),
        "segments": n_nodes,
        "links": len(gfa.edges),
        "paths": len(gfa.paths),
        "components": len(comps),
        "largest_component_frac": f"{largest_component_frac:.6f}",
        "tips": len(tips),
        "internal_tips": len(internal_tips),
        "common_start": start_node,
        "common_start_frac": f"{start_frac:.6f}",
        "common_start_count": start_count,
        "common_end": end_node,
        "common_end_frac": f"{end_frac:.6f}",
        "common_end_count": end_count,
        "link_jump_p95": link_p95,
        "link_jump_p99": link_p99,
        "link_jump_max": link_max,
        "path_jump_p95": path_p95,
        "path_jump_p99": path_p99,
        "path_jump_max": path_max,
        "duplicate_sequences": duplicate_sequences,
        "duplicate_sequence_frac": f"{duplicate_sequence_frac:.6f}",
        "max_duplicate_count": max_duplicate_count,
        "local_repeat_context_nodes": local_repeat_context_nodes,
        "local_repeat_context_occurrences": local_repeat_context_occurrences,
    }


FIELDS = [
    "file",
    "status",
    "failures",
    "warnings",
    "segments",
    "links",
    "paths",
    "components",
    "largest_component_frac",
    "tips",
    "internal_tips",
    "common_start",
    "common_start_frac",
    "common_start_count",
    "common_end",
    "common_end_frac",
    "common_end_count",
    "link_jump_p95",
    "link_jump_p99",
    "link_jump_max",
    "path_jump_p95",
    "path_jump_p99",
    "path_jump_max",
    "duplicate_sequences",
    "duplicate_sequence_frac",
    "max_duplicate_count",
    "local_repeat_context_nodes",
    "local_repeat_context_occurrences",
]


def row_to_tsv(row: dict[str, object]) -> str:
    return "\t".join(str(row[field]) for field in FIELDS)


def main() -> int:
    args = parse_args()
    rows = [evaluate(Path(raw), args) for raw in args.gfa]
    lines = ["\t".join(FIELDS), *(row_to_tsv(row) for row in rows)]
    output = "\n".join(lines) + "\n"
    if args.tsv:
        args.tsv.parent.mkdir(parents=True, exist_ok=True)
        args.tsv.write_text(output)
    if not args.quiet or not args.tsv:
        print(output, end="")
    return 1 if any(row["status"] == "FAIL" for row in rows) else 0


if __name__ == "__main__":
    raise SystemExit(main())
