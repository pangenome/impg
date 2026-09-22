#!/usr/bin/env python3
"""Diagnose C4 residual underalignment motifs and run local polish variants.

This is an experiment driver for task `diagnose-residual-underaligned`.
It intentionally writes heavy outputs outside the repository and keeps only
the reusable driver and final report under version control.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


DEFAULT_GFA = Path(
    "/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.poa2kb.gfa"
)
DEFAULT_OUT = Path("/home/erikg/impg/data/c4_residual_motif_polish_20260605T150000Z")
DEFAULT_COMPARE = Path("target/release/examples/compare_gfa_paths")
PUBLIC_BASE = "https://hypervolu.me/~erik/impg/"
UPLOAD_TARGET = "erik@hypervolu.me:www/impg/"
DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass(frozen=True)
class Step:
    segment: str
    strand: str

    @property
    def token(self) -> str:
        return f"{self.segment}{self.strand}"


@dataclass
class PathRecord:
    name: str
    steps: list[Step]
    prefix_bp: list[int]


@dataclass
class Gfa:
    path: Path
    segments: dict[str, str]
    order: dict[str, int]
    order_bp: dict[str, int]
    links: list[tuple[str, str, str, str]]
    paths: list[PathRecord]
    path_by_name: dict[str, PathRecord]
    coverage_visits: Counter[str]
    coverage_paths: dict[str, set[str]]
    transition_support: Counter[tuple[str, str]]
    link_set: set[tuple[str, str, str, str]]
    self_loop_nodes: set[str]


@dataclass
class Candidate:
    chunk_id: str
    kind: str
    example_path: str
    core_start_step: int
    core_end_step: int
    left_anchor: Step
    right_anchor: Step
    core_steps: list[Step]
    involved_paths: list[str]
    core_bp: int
    max_order_jump: int
    max_white_space_bp: int
    link_support: int
    repeat_max_run: int = 0
    score: int = 0


@dataclass
class Window:
    path_name: str
    start_step: int
    end_step: int
    steps: list[Step]
    sequence: str

    @property
    def traversal_bp(self) -> int:
        return len(self.sequence)


@dataclass
class VariantResult:
    chunk_id: str
    variant: str
    method: str
    min_match: str
    gfa: Path
    exit_status: int
    compare_status: str
    compare_stdout: Path | None
    compare_stderr: Path | None
    segments: int = 0
    links: int = 0
    paths: int = 0
    path_steps: int = 0
    segment_bp: int = 0
    singleton_bp: int = 0
    self_loop_count: int = 0
    visual_improved: str = ""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gfa", type=Path, default=DEFAULT_GFA)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--compare-bin", type=Path, default=DEFAULT_COMPARE)
    parser.add_argument("--threads", type=int, default=8)
    parser.add_argument("--max-offshoots", type=int, default=5)
    parser.add_argument("--max-selfloops", type=int, default=3)
    parser.add_argument("--test-chunks", type=int, default=5)
    parser.add_argument("--flank-steps", type=int, default=1)
    parser.add_argument("--max-anchor-steps", type=int, default=120)
    parser.add_argument("--max-window-bp", type=int, default=10_000)
    parser.add_argument("--min-window-paths", type=int, default=20)
    parser.add_argument("--max-sparse-paths", type=int, default=5)
    parser.add_argument("--min-flank-paths", type=int, default=100)
    parser.add_argument("--min-order-jump", type=int, default=1000)
    parser.add_argument("--min-white-space-bp", type=int, default=10_000)
    parser.add_argument("--ks", default="79,127,191,311")
    parser.add_argument("--impg", default="impg")
    parser.add_argument("--gfasort", default="gfasort")
    parser.add_argument("--gfalook", default="gfalook")
    parser.add_argument("--abpoa-bin", default="abpoa")
    parser.add_argument("--upload-target", default=UPLOAD_TARGET)
    parser.add_argument("--public-base", default=PUBLIC_BASE)
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--no-render", action="store_true")
    parser.add_argument("--no-upload", action="store_true")
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def token_to_step(token: str) -> Step:
    token = token.strip()
    if not token:
        raise ValueError("empty GFA path token")
    strand = token[-1]
    if strand not in "+-":
        raise ValueError(f"path token lacks strand: {token}")
    return Step(token[:-1], strand)


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def step_sequence(gfa: Gfa, step: Step) -> str:
    seq = gfa.segments[step.segment]
    return seq if step.strand == "+" else revcomp(seq)


def path_subsequence(gfa: Gfa, steps: list[Step]) -> str:
    return "".join(step_sequence(gfa, step) for step in steps)


def parse_gfa(path: Path) -> Gfa:
    segments: dict[str, str] = {}
    order: dict[str, int] = {}
    order_bp: dict[str, int] = {}
    links: list[tuple[str, str, str, str]] = []
    paths: list[PathRecord] = []
    cumulative_bp = 0

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S":
                name = fields[1]
                seq = fields[2]
                order[name] = len(order)
                order_bp[name] = cumulative_bp
                segments[name] = seq
                cumulative_bp += len(seq)
            elif fields[0] == "L":
                links.append((fields[1], fields[2], fields[3], fields[4]))
            elif fields[0] == "P":
                steps = [token_to_step(token) for token in fields[2].split(",") if token]
                prefix = [0]
                total = 0
                for step in steps:
                    total += len(segments[step.segment])
                    prefix.append(total)
                paths.append(PathRecord(fields[1], steps, prefix))

    coverage_visits: Counter[str] = Counter()
    coverage_paths: dict[str, set[str]] = defaultdict(set)
    transition_support: Counter[tuple[str, str]] = Counter()
    for path_record in paths:
        for step in path_record.steps:
            coverage_visits[step.segment] += 1
            coverage_paths[step.segment].add(path_record.name)
        for left, right in zip(path_record.steps, path_record.steps[1:]):
            transition_support[(left.token, right.token)] += 1

    link_set = set(links)
    self_loop_nodes = {src for src, _so, dst, _do in links if src == dst}
    return Gfa(
        path=path,
        segments=segments,
        order=order,
        order_bp=order_bp,
        links=links,
        paths=paths,
        path_by_name={path_record.name: path_record for path_record in paths},
        coverage_visits=coverage_visits,
        coverage_paths=dict(coverage_paths),
        transition_support=transition_support,
        link_set=link_set,
        self_loop_nodes=self_loop_nodes,
    )


def segment_path_count(gfa: Gfa, segment: str) -> int:
    return len(gfa.coverage_paths.get(segment, ()))


def order_jump(gfa: Gfa, left: Step, right: Step) -> int:
    return abs(gfa.order[left.segment] - gfa.order[right.segment])


def white_space_bp(gfa: Gfa, left: Step, right: Step) -> int:
    left_pos = gfa.order_bp[left.segment]
    right_pos = gfa.order_bp[right.segment]
    left_len = len(gfa.segments[left.segment])
    right_len = len(gfa.segments[right.segment])
    lo_end = min(left_pos + left_len, right_pos + right_len)
    hi_start = max(left_pos, right_pos)
    return max(0, hi_start - lo_end)


def core_involved_paths(gfa: Gfa, core_steps: Iterable[Step]) -> list[str]:
    names: set[str] = set()
    for step in core_steps:
        names.update(gfa.coverage_paths.get(step.segment, set()))
    return sorted(names)


def discover_offshoots(gfa: Gfa, args: argparse.Namespace) -> list[Candidate]:
    candidates: list[Candidate] = []
    seen_cores: set[tuple[str, ...]] = set()
    for path_record in gfa.paths:
        steps = path_record.steps
        i = 1
        while i < len(steps) - 1:
            if segment_path_count(gfa, steps[i].segment) > args.max_sparse_paths:
                i += 1
                continue
            start = i
            while (
                i < len(steps) - 1
                and segment_path_count(gfa, steps[i].segment) <= args.max_sparse_paths
            ):
                i += 1
            end = i - 1
            left = steps[start - 1]
            right = steps[end + 1]
            if (
                segment_path_count(gfa, left.segment) < args.min_flank_paths
                or segment_path_count(gfa, right.segment) < args.min_flank_paths
            ):
                continue
            core = steps[start : end + 1]
            core_key = tuple(step.token for step in core)
            if core_key in seen_cores:
                continue
            seen_cores.add(core_key)
            jumps = [
                order_jump(gfa, steps[start - 1], steps[start]),
                order_jump(gfa, steps[end], steps[end + 1]),
            ]
            spaces = [
                white_space_bp(gfa, steps[start - 1], steps[start]),
                white_space_bp(gfa, steps[end], steps[end + 1]),
            ]
            if max(jumps) < args.min_order_jump and max(spaces) < args.min_white_space_bp:
                continue
            involved = core_involved_paths(gfa, core)
            if len(involved) > args.max_sparse_paths:
                continue
            core_bp = sum(len(gfa.segments[step.segment]) for step in core)
            support = max(
                gfa.transition_support[(steps[start - 1].token, steps[start].token)],
                gfa.transition_support[(steps[end].token, steps[end + 1].token)],
            )
            chunk_id = f"off_{len(candidates) + 1:03d}_{core[0].segment}_{start}"
            score = max(jumps) * 1_000_000 + max(spaces) * 10 + core_bp
            candidates.append(
                Candidate(
                    chunk_id=chunk_id,
                    kind="singleton_offshoot",
                    example_path=path_record.name,
                    core_start_step=start,
                    core_end_step=end,
                    left_anchor=left,
                    right_anchor=right,
                    core_steps=core,
                    involved_paths=involved,
                    core_bp=core_bp,
                    max_order_jump=max(jumps),
                    max_white_space_bp=max(spaces),
                    link_support=support,
                    score=score,
                )
            )
    candidates.sort(key=lambda c: (-c.score, c.example_path, c.core_start_step))
    for rank, candidate in enumerate(candidates, 1):
        candidate.chunk_id = f"off_{rank:03d}_{candidate.core_steps[0].segment}_{candidate.core_start_step}"
    return candidates


def discover_selfloops(gfa: Gfa, args: argparse.Namespace) -> list[Candidate]:
    best: dict[str, tuple[int, PathRecord, int, int]] = {}
    for path_record in gfa.paths:
        steps = path_record.steps
        i = 0
        while i < len(steps):
            j = i + 1
            while j < len(steps) and steps[j].token == steps[i].token:
                j += 1
            run_len = j - i
            segment = steps[i].segment
            if segment in gfa.self_loop_nodes and run_len >= 2 and i > 0 and j < len(steps):
                prior = best.get(segment)
                if prior is None or run_len > prior[0]:
                    best[segment] = (run_len, path_record, i, j - 1)
            i = j

    candidates: list[Candidate] = []
    for segment, (run_len, path_record, start, end) in best.items():
        steps = path_record.steps
        left = steps[start - 1]
        right = steps[end + 1]
        if (
            segment_path_count(gfa, left.segment) < args.min_flank_paths
            or segment_path_count(gfa, right.segment) < args.min_flank_paths
        ):
            continue
        core = steps[start : end + 1]
        involved = core_involved_paths(gfa, core)
        core_bp = sum(len(gfa.segments[step.segment]) for step in core)
        jumps = [
            order_jump(gfa, steps[start - 1], steps[start]),
            order_jump(gfa, steps[end], steps[end + 1]),
        ]
        spaces = [
            white_space_bp(gfa, steps[start - 1], steps[start]),
            white_space_bp(gfa, steps[end], steps[end + 1]),
        ]
        chunk_id = f"loop_{len(candidates) + 1:03d}_{segment}_{start}"
        score = run_len * 10_000_000 + gfa.coverage_visits[segment] * 100 + max(jumps)
        candidates.append(
            Candidate(
                chunk_id=chunk_id,
                kind="self_loop_repeat",
                example_path=path_record.name,
                core_start_step=start,
                core_end_step=end,
                left_anchor=left,
                right_anchor=right,
                core_steps=core,
                involved_paths=involved,
                core_bp=core_bp,
                max_order_jump=max(jumps),
                max_white_space_bp=max(spaces),
                link_support=gfa.coverage_visits[segment],
                repeat_max_run=run_len,
                score=score,
            )
        )
    candidates.sort(key=lambda c: (-c.score, c.example_path, c.core_start_step))
    for rank, candidate in enumerate(candidates, 1):
        candidate.chunk_id = f"loop_{rank:03d}_{candidate.core_steps[0].segment}_{candidate.core_start_step}"
    return candidates


def find_anchor_pair(path_record: PathRecord, left: Step, right: Step, max_steps: int) -> tuple[int, int] | None:
    left_indices = [i for i, step in enumerate(path_record.steps) if step == left]
    if not left_indices:
        return None
    right_indices = [i for i, step in enumerate(path_record.steps) if step == right]
    if not right_indices:
        return None
    best: tuple[int, int] | None = None
    best_span = sys.maxsize
    for left_i in left_indices:
        for right_i in right_indices:
            if right_i <= left_i:
                continue
            span = right_i - left_i + 1
            if span <= max_steps and span < best_span:
                best = (left_i, right_i)
                best_span = span
    return best


def extract_windows(gfa: Gfa, candidate: Candidate, args: argparse.Namespace) -> list[Window]:
    windows: list[Window] = []
    for path_record in gfa.paths:
        pair = find_anchor_pair(
            path_record,
            candidate.left_anchor,
            candidate.right_anchor,
            args.max_anchor_steps,
        )
        if pair is None:
            continue
        start, end = pair
        steps = path_record.steps[start : end + 1]
        seq = path_subsequence(gfa, steps)
        if not seq or len(seq) > args.max_window_bp:
            continue
        windows.append(Window(path_record.name, start, end, steps, seq))
    windows.sort(key=lambda w: (w.path_name, w.start_step, w.end_step))
    return windows


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def fasta_escape(name: str) -> str:
    # GFA path names in this data set have no whitespace. Guard anyway so
    # seqwish/impg receives names that can round-trip as FASTA headers.
    return re.sub(r"\s+", "_", name)


def wrap_sequence(seq: str, width: int = 80) -> str:
    return "\n".join(seq[i : i + width] for i in range(0, len(seq), width))


def materialize_window_gfa(gfa: Gfa, candidate_dir: Path, candidate: Candidate, windows: list[Window]) -> tuple[Path, Path]:
    fasta = candidate_dir / f"{candidate.chunk_id}.windows.fa"
    local_gfa = candidate_dir / f"{candidate.chunk_id}.input.gfa"
    used_segments: set[str] = set()
    used_edges: set[tuple[str, str, str, str]] = set()
    for window in windows:
        used_segments.update(step.segment for step in window.steps)
        for left, right in zip(window.steps, window.steps[1:]):
            edge = (left.segment, left.strand, right.segment, right.strand)
            if edge in gfa.link_set:
                used_edges.add(edge)

    with fasta.open("w", encoding="utf-8") as handle:
        for window in windows:
            handle.write(f">{fasta_escape(window.path_name)}\n{wrap_sequence(window.sequence)}\n")

    with local_gfa.open("w", encoding="utf-8") as handle:
        handle.write("H\tVN:Z:1.0\n")
        for segment in sorted(used_segments, key=lambda s: gfa.order[s]):
            handle.write(f"S\t{segment}\t{gfa.segments[segment]}\n")
        for src, so, dst, do in sorted(
            used_edges,
            key=lambda e: (gfa.order[e[0]], gfa.order[e[2]], e[0], e[1], e[2], e[3]),
        ):
            handle.write(f"L\t{src}\t{so}\t{dst}\t{do}\t0M\n")
        for window in windows:
            tokens = ",".join(step.token for step in window.steps)
            handle.write(f"P\t{window.path_name}\t{tokens}\t*\n")
    return fasta, local_gfa


def gfa_metrics(path: Path) -> dict[str, int]:
    if not path.exists():
        return {
            "segments": 0,
            "links": 0,
            "paths": 0,
            "path_steps": 0,
            "segment_bp": 0,
            "singleton_bp": 0,
            "self_loop_count": 0,
        }
    gfa = parse_gfa(path)
    singleton_bp = sum(
        len(seq)
        for segment, seq in gfa.segments.items()
        if gfa.coverage_visits.get(segment, 0) == 1
    )
    return {
        "segments": len(gfa.segments),
        "links": len(gfa.links),
        "paths": len(gfa.paths),
        "path_steps": sum(len(path_record.steps) for path_record in gfa.paths),
        "segment_bp": sum(len(seq) for seq in gfa.segments.values()),
        "singleton_bp": singleton_bp,
        "self_loop_count": len(gfa.self_loop_nodes),
    }


def shell_join(command: Iterable[object]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def run_command(
    command: list[object],
    stdout: Path,
    stderr: Path,
    *,
    cwd: Path | None = None,
    env: dict[str, str] | None = None,
) -> int:
    stdout.parent.mkdir(parents=True, exist_ok=True)
    stderr.parent.mkdir(parents=True, exist_ok=True)
    with stdout.open("wb") as out, stderr.open("wb") as err:
        proc = subprocess.run([str(part) for part in command], cwd=cwd, env=env, stdout=out, stderr=err)
    return proc.returncode


def compare_paths(compare_bin: Path, expected: Path, observed: Path, out_prefix: Path) -> tuple[str, Path, Path]:
    stdout = out_prefix.parent / f"{out_prefix.name}.compare.stdout.log"
    stderr = out_prefix.parent / f"{out_prefix.name}.compare.stderr.log"
    if not compare_bin.exists():
        return "missing_compare_bin", stdout, stderr
    rc = run_command([compare_bin, expected, observed], stdout, stderr)
    return ("pass" if rc == 0 else f"fail:{rc}"), stdout, stderr


def run_seqwish_variant(
    args: argparse.Namespace,
    candidate_dir: Path,
    chunk_id: str,
    fasta: Path,
    local_input: Path,
    k: int,
) -> VariantResult:
    label = f"seqwish_k{k}"
    out_gfa = candidate_dir / f"{chunk_id}.{label}.gfa"
    logs = candidate_dir / "logs"
    debug_dir = candidate_dir / f"debug_{label}"
    temp_dir = Path("/tmp") / f"c4motif_{chunk_id}_{label}"
    command = [
        args.impg,
        "graph",
        "--sequence-files",
        fasta,
        "--gfa-engine",
        "seqwish",
        "--fastga",
        "--num-mappings",
        "1:many",
        "--scaffold-filter",
        "1:many",
        "--scaffold-jump",
        "0",
        "--min-match-len",
        str(k),
        "--temp-dir",
        temp_dir,
        "--debug-dir",
        debug_dir,
        "-g",
        out_gfa,
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    rc = 0
    if args.force or not out_gfa.exists():
        rc = run_command(command, logs / f"{label}.stdout.log", logs / f"{label}.stderr.log")
    status, cmp_out, cmp_err = compare_paths(args.compare_bin, local_input, out_gfa, candidate_dir / f"{chunk_id}.{label}")
    metrics = gfa_metrics(out_gfa)
    return VariantResult(
        chunk_id=chunk_id,
        variant=label,
        method="seqwish",
        min_match=str(k),
        gfa=out_gfa,
        exit_status=rc,
        compare_status=status,
        compare_stdout=cmp_out,
        compare_stderr=cmp_err,
        **metrics,
    )


def run_poa_variant(
    args: argparse.Namespace,
    candidate_dir: Path,
    chunk_id: str,
    fasta: Path,
    local_input: Path,
) -> VariantResult:
    label = "poa"
    out_gfa = candidate_dir / f"{chunk_id}.{label}.gfa"
    logs = candidate_dir / "logs"
    command = [
        args.impg,
        "graph",
        "--sequence-files",
        fasta,
        "--gfa-engine",
        "poa",
        "-g",
        out_gfa,
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    rc = 0
    if args.force or not out_gfa.exists():
        rc = run_command(command, logs / f"{label}.stdout.log", logs / f"{label}.stderr.log")
    status, cmp_out, cmp_err = compare_paths(args.compare_bin, local_input, out_gfa, candidate_dir / f"{chunk_id}.{label}")
    metrics = gfa_metrics(out_gfa)
    return VariantResult(
        chunk_id=chunk_id,
        variant=label,
        method="poa",
        min_match=".",
        gfa=out_gfa,
        exit_status=rc,
        compare_status=status,
        compare_stdout=cmp_out,
        compare_stderr=cmp_err,
        **metrics,
    )


def run_abpoa_variant(
    args: argparse.Namespace,
    candidate_dir: Path,
    chunk_id: str,
    local_input: Path,
) -> VariantResult:
    label = "abpoa_crush"
    out_gfa = candidate_dir / f"{chunk_id}.{label}.gfa"
    logs = candidate_dir / "logs"
    command = [
        args.impg,
        "crush",
        "-g",
        local_input,
        "-o",
        out_gfa,
        "--method",
        "abpoa",
        "--abpoa-bin",
        args.abpoa_bin,
        "--max-iterations",
        "3",
        "--max-traversal-len",
        "10k",
        "--max-median-traversal-len",
        "10k",
        "--max-total-sequence",
        "10m",
        "--max-traversals",
        "10000",
        "--polish-rounds",
        "0",
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    rc = 0
    if args.force or not out_gfa.exists():
        rc = run_command(command, logs / f"{label}.stdout.log", logs / f"{label}.stderr.log")
    status, cmp_out, cmp_err = compare_paths(args.compare_bin, local_input, out_gfa, candidate_dir / f"{chunk_id}.{label}")
    metrics = gfa_metrics(out_gfa)
    return VariantResult(
        chunk_id=chunk_id,
        variant=label,
        method="abpoa_crush",
        min_match=".",
        gfa=out_gfa,
        exit_status=rc,
        compare_status=status,
        compare_stdout=cmp_out,
        compare_stderr=cmp_err,
        **metrics,
    )


def score_metrics(result: VariantResult) -> tuple[int, int, int, int]:
    return (result.singleton_bp, result.self_loop_count, result.segment_bp, result.segments)


def visually_improves(input_result: VariantResult, result: VariantResult) -> bool:
    if result.compare_status != "pass":
        return False
    if result.singleton_bp > input_result.singleton_bp:
        return False
    if result.segment_bp > input_result.segment_bp * 2:
        return False
    return (
        result.singleton_bp < input_result.singleton_bp
        or result.self_loop_count < input_result.self_loop_count
        or result.segment_bp < input_result.segment_bp
    )


def render_pair(args: argparse.Namespace, out_dir: Path, input_gfa: Path, best: VariantResult) -> dict[str, str]:
    render_dir = out_dir / "renders"
    render_dir.mkdir(parents=True, exist_ok=True)
    rows: dict[str, str] = {}
    for label, gfa_path in [("before", input_gfa), ("after", best.gfa)]:
        sorted_gfa = render_dir / f"{best.chunk_id}.{label}.Ygs.gfa"
        png = render_dir / f"{best.chunk_id}.{label}.Ygs.gfalook-m.png"
        if args.force or not sorted_gfa.exists():
            run_command(
                [args.gfasort, "-i", gfa_path, "-o", sorted_gfa, "-p", "Ygs", "-t", str(args.threads), "-v", "1"],
                render_dir / f"{best.chunk_id}.{label}.gfasort.stdout.log",
                render_dir / f"{best.chunk_id}.{label}.gfasort.stderr.log",
            )
        if args.force or not png.exists():
            run_command(
                [args.gfalook, "-i", sorted_gfa, "-o", png, "-m", "-x", "1800", "-y", "1000", "-a", "3", "-t", str(args.threads), "-v", "1"],
                render_dir / f"{best.chunk_id}.{label}.gfalook.stdout.log",
                render_dir / f"{best.chunk_id}.{label}.gfalook.stderr.log",
            )
        rows[label] = str(png)
    return rows


def upload_pngs(args: argparse.Namespace, pngs: dict[str, str]) -> dict[str, str]:
    if args.no_upload:
        return {}
    uploaded: dict[str, str] = {}
    for label, png in pngs.items():
        path = Path(png)
        remote_name = path.name
        rc = subprocess.run(["scp", str(path), f"{args.upload_target}{remote_name}"]).returncode
        if rc == 0:
            uploaded[label] = f"{args.public_base}{remote_name}"
    return uploaded


def write_command_manifest(args: argparse.Namespace, out_dir: Path) -> None:
    with (out_dir / "commands.md").open("w", encoding="utf-8") as handle:
        handle.write("# C4 motif-local polish commands\n\n")
        handle.write("Driver invocation:\n\n```bash\n")
        handle.write(shell_join([sys.executable, *sys.argv]) + "\n")
        handle.write("```\n\n")
        handle.write("Tool paths:\n\n")
        for name in ["impg", "gfasort", "gfalook", "abpoa-bin"]:
            value = getattr(args, name.replace("-", "_"))
            handle.write(f"- {name}: `{value}`\n")


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    write_command_manifest(args, args.out_dir)
    gfa = parse_gfa(args.gfa)
    offshoots = discover_offshoots(gfa, args)
    selfloops = discover_selfloops(gfa, args)

    discovered_rows: list[dict[str, object]] = []
    for source, candidates in [("offshoot", offshoots), ("selfloop", selfloops)]:
        for rank, candidate in enumerate(candidates, 1):
            discovered_rows.append(
                {
                    "source": source,
                    "rank": rank,
                    "chunk_id": candidate.chunk_id,
                    "kind": candidate.kind,
                    "example_path": candidate.example_path,
                    "core_start_step": candidate.core_start_step,
                    "core_end_step": candidate.core_end_step,
                    "left_anchor": candidate.left_anchor.token,
                    "right_anchor": candidate.right_anchor.token,
                    "core_steps": ",".join(step.token for step in candidate.core_steps),
                    "core_bp": candidate.core_bp,
                    "involved_path_count": len(candidate.involved_paths),
                    "involved_paths": ";".join(candidate.involved_paths),
                    "max_order_jump": candidate.max_order_jump,
                    "max_white_space_bp": candidate.max_white_space_bp,
                    "link_support": candidate.link_support,
                    "repeat_max_run": candidate.repeat_max_run,
                    "score": candidate.score,
                }
            )
    discovered_fields = [
        "source",
        "rank",
        "chunk_id",
        "kind",
        "example_path",
        "core_start_step",
        "core_end_step",
        "left_anchor",
        "right_anchor",
        "core_steps",
        "core_bp",
        "involved_path_count",
        "involved_paths",
        "max_order_jump",
        "max_white_space_bp",
        "link_support",
        "repeat_max_run",
        "score",
    ]
    write_tsv(args.out_dir / "discovered_chunks.tsv", discovered_rows, discovered_fields)

    selected_pool = offshoots[: args.max_offshoots] + selfloops[: args.max_selfloops]
    selected: list[tuple[Candidate, list[Window]]] = []
    window_rows: list[dict[str, object]] = []
    for candidate in selected_pool:
        windows = extract_windows(gfa, candidate, args)
        if len(windows) < args.min_window_paths:
            continue
        selected.append((candidate, windows))
        if len(selected) >= args.test_chunks:
            break
    # If top offshoots filled the list, force at least one self-loop candidate in
    # the tested representative set so direct repeat motifs are exercised.
    if not any(candidate.kind == "self_loop_repeat" for candidate, _ in selected):
        for candidate in selfloops:
            windows = extract_windows(gfa, candidate, args)
            if len(windows) >= args.min_window_paths:
                if len(selected) >= args.test_chunks:
                    selected[-1] = (candidate, windows)
                else:
                    selected.append((candidate, windows))
                break

    all_results: list[VariantResult] = []
    selected_rows: list[dict[str, object]] = []
    best_result: VariantResult | None = None
    best_input: Path | None = None
    ks = [int(k.strip()) for k in args.ks.split(",") if k.strip()]
    for candidate, windows in selected:
        candidate_dir = args.out_dir / "chunks" / candidate.chunk_id
        candidate_dir.mkdir(parents=True, exist_ok=True)
        fasta, local_input = materialize_window_gfa(gfa, candidate_dir, candidate, windows)
        input_metrics = gfa_metrics(local_input)
        max_traversal_bp = max((window.traversal_bp for window in windows), default=0)
        median_traversal_bp = sorted(window.traversal_bp for window in windows)[len(windows) // 2]
        selected_rows.append(
            {
                "chunk_id": candidate.chunk_id,
                "kind": candidate.kind,
                "example_path": candidate.example_path,
                "core_start_step": candidate.core_start_step,
                "core_end_step": candidate.core_end_step,
                "left_anchor": candidate.left_anchor.token,
                "right_anchor": candidate.right_anchor.token,
                "core_steps": ",".join(step.token for step in candidate.core_steps),
                "core_bp": candidate.core_bp,
                "involved_path_count": len(candidate.involved_paths),
                "involved_paths": ";".join(candidate.involved_paths),
                "window_path_count": len(windows),
                "min_traversal_bp": min(window.traversal_bp for window in windows),
                "median_traversal_bp": median_traversal_bp,
                "max_traversal_bp": max_traversal_bp,
                "input_gfa": local_input,
                "fasta": fasta,
                **{f"input_{key}": value for key, value in input_metrics.items()},
            }
        )
        for window in windows:
            window_rows.append(
                {
                    "chunk_id": candidate.chunk_id,
                    "path_name": window.path_name,
                    "start_step": window.start_step,
                    "end_step": window.end_step,
                    "step_count": len(window.steps),
                    "traversal_bp": window.traversal_bp,
                    "start_segment": window.steps[0].token,
                    "end_segment": window.steps[-1].token,
                }
            )
        if args.no_run:
            continue
        results: list[VariantResult] = []
        for k in ks:
            results.append(run_seqwish_variant(args, candidate_dir, candidate.chunk_id, fasta, local_input, k))
        if max_traversal_bp <= 2_000:
            results.append(run_poa_variant(args, candidate_dir, candidate.chunk_id, fasta, local_input))
        if max_traversal_bp <= 10_000:
            results.append(run_abpoa_variant(args, candidate_dir, candidate.chunk_id, local_input))
        input_result = VariantResult(
            chunk_id=candidate.chunk_id,
            variant="input",
            method="input",
            min_match=".",
            gfa=local_input,
            exit_status=0,
            compare_status="self",
            compare_stdout=None,
            compare_stderr=None,
            **input_metrics,
        )
        passing = [result for result in results if result.exit_status == 0 and result.compare_status == "pass"]
        if passing:
            best = min(passing, key=score_metrics)
            for result in results:
                result.visual_improved = "yes" if visually_improves(input_result, result) else "no"
            improving = [result for result in passing if visually_improves(input_result, result)]
            if improving:
                best = min(improving, key=score_metrics)
                if best_result is None or score_metrics(best) < score_metrics(best_result):
                    best_result = best
                    best_input = local_input
        all_results.extend([input_result, *results])

    selected_fields = [
        "chunk_id",
        "kind",
        "example_path",
        "core_start_step",
        "core_end_step",
        "left_anchor",
        "right_anchor",
        "core_steps",
        "core_bp",
        "involved_path_count",
        "involved_paths",
        "window_path_count",
        "min_traversal_bp",
        "median_traversal_bp",
        "max_traversal_bp",
        "input_gfa",
        "fasta",
        "input_segments",
        "input_links",
        "input_paths",
        "input_path_steps",
        "input_segment_bp",
        "input_singleton_bp",
        "input_self_loop_count",
    ]
    write_tsv(args.out_dir / "selected_chunks.tsv", selected_rows, selected_fields)
    write_tsv(
        args.out_dir / "windows.tsv",
        window_rows,
        [
            "chunk_id",
            "path_name",
            "start_step",
            "end_step",
            "step_count",
            "traversal_bp",
            "start_segment",
            "end_segment",
        ],
    )

    result_rows = [
        {
            "chunk_id": result.chunk_id,
            "variant": result.variant,
            "method": result.method,
            "min_match": result.min_match,
            "gfa": result.gfa,
            "exit_status": result.exit_status,
            "compare_status": result.compare_status,
            "compare_stdout": result.compare_stdout or "",
            "compare_stderr": result.compare_stderr or "",
            "segments": result.segments,
            "links": result.links,
            "paths": result.paths,
            "path_steps": result.path_steps,
            "segment_bp": result.segment_bp,
            "singleton_bp": result.singleton_bp,
            "self_loop_count": result.self_loop_count,
            "visual_improved": result.visual_improved,
        }
        for result in all_results
    ]
    write_tsv(
        args.out_dir / "variant_results.tsv",
        result_rows,
        [
            "chunk_id",
            "variant",
            "method",
            "min_match",
            "gfa",
            "exit_status",
            "compare_status",
            "compare_stdout",
            "compare_stderr",
            "segments",
            "links",
            "paths",
            "path_steps",
            "segment_bp",
            "singleton_bp",
            "self_loop_count",
            "visual_improved",
        ],
    )

    if best_result is not None and best_input is not None and not args.no_render:
        pngs = render_pair(args, args.out_dir, best_input, best_result)
        uploads = upload_pngs(args, pngs)
        upload_rows = [
            {"label": label, "png": pngs.get(label, ""), "url": uploads.get(label, "")}
            for label in sorted(pngs)
        ]
        write_tsv(args.out_dir / "uploaded_pngs.tsv", upload_rows, ["label", "png", "url"])

    print(f"wrote {args.out_dir}")
    print(f"discovered offshoots={len(offshoots)} selfloops={len(selfloops)} selected={len(selected)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
