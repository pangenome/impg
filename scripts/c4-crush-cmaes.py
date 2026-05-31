#!/usr/bin/env python3
"""
CMA-ES optimizer for C4 syng/crush graph construction.

The optimizer is deliberately external to impg: every candidate is run as an
ordinary impg command, graph-report is run afterwards, and only this wrapper
marks a candidate invalid or assigns an objective value.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import json
import math
import os
import random
import re
import shlex
import shutil
import signal
import subprocess
import sys
import threading
import time
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Any, Callable, Iterable


try:
    import cma  # type: ignore

    HAVE_CMA = True
    CMA_IMPORT_ERROR = ""
except Exception as exc:  # pragma: no cover - exercised only without cma
    cma = None
    HAVE_CMA = False
    CMA_IMPORT_ERROR = str(exc)


METHOD_FAMILIES = (
    "auto",
    "sweepga",
    "chain-povu",
    "iterative-multi-level",
    "coverage-multi-bubble",
)

INVALID_OBJECTIVE = 1.0e12
PATH_STEP_RE = re.compile(r"(.+)([+-])$")
WALK_STEP_RE = re.compile(r"([<>])([^<>]+)")
DEFAULT_FREQUENCY_BIN_LABELS = (
    "0-0.01",
    "0.01-0.05",
    "0.05-0.1",
    "0.1-0.25",
    "0.25-0.5",
    "0.5-0.75",
    "0.75-1",
    "1-1.5",
    "1.5-2",
    "2-3",
    ">=3",
)
DEFAULT_LENGTH_BIN_LABELS = (
    "0-1",
    "1-2",
    "2-4",
    "4-8",
    "8-16",
    "16-32",
    "32-64",
    "64-128",
    "128-256",
    "256-512",
    "512-1024",
    "1024-2048",
    "2048-4096",
    "4096-8192",
    "8192-16384",
    "16384-32768",
    "32768-65536",
    ">=65536",
)


@dataclasses.dataclass(frozen=True)
class Dimension:
    name: str
    description: str
    decode: Callable[[float], Any]
    range_doc: str

    def value(self, raw: float) -> Any:
        return self.decode(clamp01(raw))


@dataclasses.dataclass
class TimedCommandResult:
    command: list[str]
    returncode: int
    timeout: bool
    wall_seconds: float
    stdout: str
    stderr: str
    time_file: str
    max_rss_kb: int | None
    time_metrics: dict[str, str]


@dataclasses.dataclass
class TrialPlan:
    trial_id: int
    generation: int
    candidate_index: int
    vector: list[float]
    decoded: dict[str, Any]
    trial_dir: Path
    primary_command: list[str]
    graph_report_command: list[str]
    output_gfa: Path
    graph_report: Path
    path_validation: Path


@dataclasses.dataclass
class TrialResult:
    trial_id: int
    generation: int
    candidate_index: int
    vector: list[float]
    decoded: dict[str, Any]
    valid: bool
    invalid_reasons: list[str]
    raw_objective: float | None
    tell_objective: float
    metrics: dict[str, Any]
    objective_terms: dict[str, float]
    output_gfa: str
    graph_report: str
    path_validation: str
    primary: dict[str, Any]
    graph_report_run: dict[str, Any] | None
    optimizer: str
    started_at: str
    finished_at: str


@dataclasses.dataclass(frozen=True)
class HistogramBin:
    label: str
    lower: float
    upper: float | None
    center: float


@dataclasses.dataclass(frozen=True)
class GraphShape:
    path_count: int
    segment_count: int
    total_bp: int
    path_steps: int
    bp_weighted_depth_mean: float
    singleton_bp: int
    frequency_distribution: tuple[float, ...]
    length_distribution: tuple[float, ...]


@dataclasses.dataclass(frozen=True)
class TargetProfile:
    source_gfa: str | None
    source_freq_bins: str | None
    source_summary: str | None
    freq_graph: str | None
    frequency_bins: tuple[HistogramBin, ...]
    frequency_distribution: tuple[float, ...]
    length_bins: tuple[HistogramBin, ...]
    length_distribution: tuple[float, ...] | None
    total_bp: float | None
    path_count: int | None
    bp_weighted_depth_mean: float | None
    singleton_bp: float | None


class Optimizer:
    name = "fallback-random"

    def ask(self, count: int) -> list[list[float]]:
        raise NotImplementedError

    def tell(self, vectors: list[list[float]], objectives: list[float]) -> None:
        raise NotImplementedError

    def replay_generation(self, vectors: list[list[float]], objectives: list[float]) -> None:
        self.tell(vectors, objectives)


class CmaOptimizer(Optimizer):
    name = "cma-es"

    def __init__(self, dimensions: int, population_size: int, seed: int) -> None:
        options = {
            "bounds": [0.0, 1.0],
            "popsize": population_size,
            "seed": seed,
            "verb_disp": 0,
            "verb_log": 0,
            "verbose": -9,
        }
        self.es = cma.CMAEvolutionStrategy([0.5] * dimensions, 0.22, options)

    def ask(self, count: int) -> list[list[float]]:
        return [list(map(float, vector)) for vector in self.es.ask(count)]

    def tell(self, vectors: list[list[float]], objectives: list[float]) -> None:
        self.es.tell(vectors, objectives)

    def replay_generation(self, vectors: list[list[float]], objectives: list[float]) -> None:
        self.es.ask(len(vectors))
        self.es.tell(vectors, objectives)


class RandomOptimizer(Optimizer):
    name = "fallback-random"

    def __init__(self, dimensions: int, seed: int) -> None:
        self.dimensions = dimensions
        self.rng = random.Random(seed)

    def ask(self, count: int) -> list[list[float]]:
        return [[self.rng.random() for _ in range(self.dimensions)] for _ in range(count)]

    def tell(self, vectors: list[list[float]], objectives: list[float]) -> None:
        del vectors, objectives


def clamp01(value: float) -> float:
    if math.isnan(value):
        return 0.5
    return min(1.0, max(0.0, float(value)))


def linear_float(name: str, low: float, high: float, description: str) -> Dimension:
    def decode(raw: float) -> float:
        return low + clamp01(raw) * (high - low)

    return Dimension(name, description, decode, f"linear float [{low}, {high}]")


def log_float(name: str, low: float, high: float, description: str) -> Dimension:
    log_low = math.log(low)
    log_high = math.log(high)

    def decode(raw: float) -> float:
        return math.exp(log_low + clamp01(raw) * (log_high - log_low))

    return Dimension(name, description, decode, f"log float [{low}, {high}]")


def linear_int(name: str, low: int, high: int, description: str) -> Dimension:
    def decode(raw: float) -> int:
        return int(round(low + clamp01(raw) * (high - low)))

    return Dimension(name, description, decode, f"linear integer [{low}, {high}]")


def log_int(name: str, low: int, high: int, description: str) -> Dimension:
    log_low = math.log(low)
    log_high = math.log(high)

    def decode(raw: float) -> int:
        return int(round(math.exp(log_low + clamp01(raw) * (log_high - log_low))))

    return Dimension(name, description, decode, f"log integer [{low}, {high}]")


def odd_linear_int(name: str, low: int, high: int, description: str) -> Dimension:
    if low % 2 == 0:
        low += 1
    if high % 2 == 0:
        high -= 1

    def decode(raw: float) -> int:
        value = int(round(low + clamp01(raw) * (high - low)))
        if value % 2 == 0:
            value += 1 if value < high else -1
        return max(low, min(high, value))

    return Dimension(name, description, decode, f"linear odd integer [{low}, {high}]")


def decode_dimensions(dimensions: list[Dimension], vector: Iterable[float]) -> dict[str, Any]:
    return {dim.name: dim.value(raw) for dim, raw in zip(dimensions, vector)}


def build_dimensions(args: argparse.Namespace) -> list[Dimension]:
    dims: list[Dimension] = []
    if args.mode == "full-query":
        dims.extend(
            [
                log_float(
                    "mask_top_fraction",
                    0.00005,
                    0.003,
                    "syng mask top-frequency drop fraction",
                ),
                log_int("mask_max_occ", 100, 2500, "syng mask maximum local occurrence count"),
                linear_int("mask_min_run", 2, 8, "minimum consecutive shared syncmer run"),
                odd_linear_int(
                    "mask_sequence_k",
                    63,
                    311,
                    "exact repeated sequence support floor for weak syncmers",
                ),
            ]
        )
        if args.syng_engine == "syng-local":
            dims.extend(
                [
                    odd_linear_int("syng_k", 63, 311, "regional syng-local syncmer k"),
                    linear_int("syng_s", 8, 31, "regional syng-local syncmer s"),
                ]
            )

    dims.extend(
        [
            linear_int("max_rounds", 1, 4, "crush frontier replacement rounds"),
            log_int("min_traversal_len", 1, 20_000, "minimum longest traversal length"),
            log_int("max_traversal_len", 1_000, 120_000, "maximum direct traversal length"),
            log_int("max_median_traversal_len", 200, 30_000, "maximum median direct traversal length"),
            log_int("max_traversals", 64, 25_000, "maximum traversals per direct candidate"),
            linear_int("polish_rounds", 1, 4, "small-tangle polish rounds"),
            log_int("polish_max_traversal_len", 500, 50_000, "maximum polish traversal length"),
            log_int(
                "polish_max_median_traversal_len",
                100,
                10_000,
                "maximum median polish traversal length",
            ),
            log_int("replacement_flank_bp", 1, 1_000, "bubble flank context bp"),
        ]
    )

    if args.method_family == "auto":
        dims.extend(
            [
                log_int("auto_poasta_max_traversal_len", 1_000, 60_000, "auto POASTA median cutoff"),
            ]
        )
    elif args.method_family == "sweepga":
        dims.extend(
            [
                linear_int("k_nearest", 1, 8, "pair-sampling nearest neighbors"),
                linear_int("k_farthest", 0, 4, "pair-sampling farthest neighbors"),
                linear_int("pair_trees", 1, 4, "independent pair-sampling trees"),
                log_float("random_fraction", 0.0001, 0.05, "deterministic random pair fraction"),
                linear_int("mash_k", 11, 21, "mash k-mer size"),
                odd_linear_int("seqwish_k", 51, 501, "replacement seqwish exact-match floor"),
                log_int("min_match_length", 31, 501, "replacement local exact-run floor"),
                log_int("replacement_min_map_length", 50, 2_000, "minimum replacement mapping length"),
                linear_float("replacement_min_identity", 0.0, 0.99, "minimum replacement identity"),
                log_int("max_pair_alignments", 500, 100_000, "maximum selected pair alignments"),
                log_int("max_replacement_paf_bytes", 1_000_000, 512_000_000, "maximum replacement PAF bytes"),
                linear_int("sweepga_kmer_frequency", 0, 100, "SweepGA/FastGA k-mer frequency"),
            ]
        )
    elif args.method_family == "chain-povu":
        dims.append(log_int("chain_target_bp", 2_000, 100_000, "target subtree block span"))
    elif args.method_family in {"iterative-multi-level", "coverage-multi-bubble"}:
        dims.extend(
            [
                log_int("window_target_bp", 2_000, 120_000, "multi-bubble window target span"),
                linear_int("max_window_sites", 2, 20, "maximum sites merged into one window"),
                log_int("candidate_limit", 16, 768, "maximum generated candidates per round"),
                linear_int("min_objective_delta", 0, 10, "diagnostic local objective floor"),
            ]
        )
    else:  # pragma: no cover - argparse choices prevent this
        raise ValueError(f"unsupported method family: {args.method_family}")
    return dims


def normalize_decoded(decoded: dict[str, Any], args: argparse.Namespace) -> dict[str, Any]:
    out = dict(decoded)
    out["max_median_traversal_len"] = min(out["max_median_traversal_len"], out["max_traversal_len"])
    out["polish_max_median_traversal_len"] = min(
        out["polish_max_median_traversal_len"], out["polish_max_traversal_len"]
    )
    if args.mode == "full-query" and args.syng_engine == "syng-local":
        out["syng_s"] = max(1, min(int(out["syng_s"]), (int(out["syng_k"]) - 1) // 2))
    return out


def crush_stage_params(decoded: dict[str, Any], args: argparse.Namespace) -> list[tuple[str, Any]]:
    params: list[tuple[str, Any]] = [
        ("method", args.method_family),
        ("max-rounds", decoded["max_rounds"]),
        ("min-traversal-len", decoded["min_traversal_len"]),
        ("max-traversal-len", decoded["max_traversal_len"]),
        ("max-median-traversal-len", decoded["max_median_traversal_len"]),
        ("max-traversals", decoded["max_traversals"]),
        ("polish-rounds", decoded["polish_rounds"]),
        ("polish-max-traversal-len", decoded["polish_max_traversal_len"]),
        ("polish-max-median-traversal-len", decoded["polish_max_median_traversal_len"]),
        ("replacement-flank-bp", decoded["replacement_flank_bp"]),
    ]
    if args.method_family == "auto":
        params.extend(
            [
                ("auto-poasta-max-len", decoded["auto_poasta_max_traversal_len"]),
            ]
        )
    elif args.method_family == "sweepga":
        params.extend(
            [
                ("k-nearest", decoded["k_nearest"]),
                ("k-farthest", decoded["k_farthest"]),
                ("pair-trees", decoded["pair_trees"]),
                ("random-fraction", f"{decoded['random_fraction']:.6g}"),
                ("mash-k", decoded["mash_k"]),
                ("seqwish-k", decoded["seqwish_k"]),
                ("min-match-length", decoded["min_match_length"]),
                ("replacement-min-map-length", decoded["replacement_min_map_length"]),
                ("replacement-min-identity", f"{decoded['replacement_min_identity']:.6g}"),
                ("max-pair-alignments", decoded["max_pair_alignments"]),
                ("max-replacement-paf-bytes", decoded["max_replacement_paf_bytes"]),
                ("kmer-frequency", decoded["sweepga_kmer_frequency"]),
            ]
        )
    elif args.method_family == "chain-povu":
        params.append(("chain-target-bp", decoded["chain_target_bp"]))
    elif args.method_family == "iterative-multi-level":
        params.extend(multi_level_params(decoded, args, objective="size", repeat_aware=False))
    elif args.method_family == "coverage-multi-bubble":
        params.extend(multi_level_params(decoded, args, objective="coverage", repeat_aware=True))
    return params


def multi_level_params(
    decoded: dict[str, Any], args: argparse.Namespace, objective: str, repeat_aware: bool
) -> list[tuple[str, Any]]:
    return [
        ("window-mode", args.multi_level_window_mode),
        ("window-target-bp", decoded["window_target_bp"]),
        ("max-window-sites", decoded["max_window_sites"]),
        ("candidate-limit", decoded["candidate_limit"]),
        ("min-objective-delta", decoded["min_objective_delta"]),
        ("objective", objective),
        ("repeat-aware-boundaries", str(repeat_aware).lower()),
    ]


def build_engine_string(decoded: dict[str, Any], args: argparse.Namespace) -> str:
    stages: list[str] = [args.syng_engine]
    if args.syng_engine == "syng-local":
        stages[0] = (
            f"syng-local:blunt,k={decoded['syng_k']},s={decoded['syng_s']},seed={args.syng_seed}"
        )
    else:
        stages[0] = "syng:blunt"
    stages.append(
        "mask,"
        f"top={decoded['mask_top_fraction']:.8g},"
        f"max-occ={decoded['mask_max_occ']},"
        f"min-run={decoded['mask_min_run']},"
        f"sequence-k={decoded['mask_sequence_k']}"
    )
    crush_params = ",".join(f"{key}={value}" for key, value in crush_stage_params(decoded, args))
    stages.append(f"crush,{crush_params}")
    if args.engine_suffix:
        stages.append(args.engine_suffix)
    return "gfa:" + ":".join(stages)


def make_optimizer(dimensions: int, args: argparse.Namespace) -> Optimizer:
    if HAVE_CMA:
        return CmaOptimizer(dimensions, args.population_size, args.seed)
    return RandomOptimizer(dimensions, args.seed)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Run ask/tell CMA-ES over C4 syng/crush graph-building parameters.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--mode", choices=["crush-only", "full-query"], required=True)
    parser.add_argument("--method-family", choices=METHOD_FAMILIES, default="auto")
    parser.add_argument("--max-trials", type=int, default=40)
    parser.add_argument("--population-size", type=int, default=8)
    parser.add_argument("--jobs", type=int, default=1)
    parser.add_argument("--threads", type=int, default=4)
    parser.add_argument("--timeout", type=int, default=1800, help="Per external command timeout in seconds")
    parser.add_argument("--seed", type=int, default=7)
    parser.add_argument("--study-dir", type=Path, required=True)
    parser.add_argument("--impg-bin", default=default_impg_bin())
    parser.add_argument("--dry-run", action="store_true", help="Print one decoded trial and commands without running impg")
    parser.add_argument("--no-resume", action="store_true", help="Ignore an existing trial log for optimizer replay")

    parser.add_argument("--input-gfa", type=Path, help="Input GFA for --mode crush-only")

    parser.add_argument("--index", help="syng index prefix/path for --mode full-query")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--target-range", help="impg query -r value for --mode full-query")
    group.add_argument("--target-bed", type=Path, help="impg query -b value; intended for one-row BED studies")
    parser.add_argument("--sequence-files", nargs="+", default=[], help="FASTA/AGC files passed to impg query")
    parser.add_argument("--sequence-list", type=Path, help="Sequence-list file passed to impg query")
    parser.add_argument("--expected-path-count", type=int, help="Required output path count for --mode full-query")
    parser.add_argument("--syng-engine", choices=["syng", "syng-local"], default="syng-local")
    parser.add_argument("--syng-seed", type=int, default=7)
    parser.add_argument("--merge-distance", default="100k")
    parser.add_argument("--engine-suffix", default="sort,pipeline=Ygs")
    parser.add_argument(
        "--multi-level-window-mode",
        choices=["sibling", "sliding", "outward", "combined"],
        default="combined",
    )

    parser.add_argument("--objective-bp-scale", type=float, help="Override total-bp objective normalization")
    parser.add_argument(
        "--target-gfa",
        type=Path,
        help="PGGB/control GFA used to derive target path-depth frequency, node length, and total-bp profiles",
    )
    parser.add_argument(
        "--target-freq-bins",
        type=Path,
        help="TSV with freq_bin and bp_frac columns; bin labels define target frequency histogram bins",
    )
    parser.add_argument(
        "--target-summary",
        type=Path,
        help="TSV summary for target metadata such as paths, segment_bp, bp_weighted_depth_mean, and singleton_bp",
    )
    parser.add_argument(
        "--target-frequency-weight",
        type=float,
        default=1000.0,
        help="Weight applied to target frequency TV+ordered-distance score",
    )
    parser.add_argument(
        "--target-length-weight",
        type=float,
        default=100.0,
        help="Weight applied to bp-weighted node length distribution distance",
    )
    parser.add_argument(
        "--target-total-bp-weight",
        type=float,
        default=25.0,
        help="Weight applied to excess total segment bp ratio relative to the target",
    )
    parser.add_argument(
        "--target-whitespace-weight",
        type=float,
        default=1.0,
        help="Weight applied to legacy path-whitespace penalties when target-shape scoring is enabled",
    )
    parser.add_argument("--report-top", type=int, default=20)
    parser.add_argument("--query-extra-arg", action="append", default=[], help="Extra single token for impg query")
    parser.add_argument("--crush-extra-arg", action="append", default=[], help="Extra single token for impg crush")
    args = parser.parse_args()
    validate_args(args)
    return args


def default_impg_bin() -> str:
    local = Path("target/release/impg")
    if local.exists():
        return str(local)
    found = shutil.which("impg")
    return found if found else str(local)


def validate_args(args: argparse.Namespace) -> None:
    if args.max_trials < 1:
        raise SystemExit("--max-trials must be >= 1")
    if args.population_size < 2:
        raise SystemExit("--population-size must be >= 2")
    if args.jobs < 1:
        raise SystemExit("--jobs must be >= 1")
    if args.threads < 1:
        raise SystemExit("--threads must be >= 1")
    if args.timeout < 1:
        raise SystemExit("--timeout must be >= 1")
    if args.report_top < 1:
        raise SystemExit("--report-top must be >= 1")
    for name in (
        "target_frequency_weight",
        "target_length_weight",
        "target_total_bp_weight",
        "target_whitespace_weight",
    ):
        if getattr(args, name) < 0.0:
            raise SystemExit(f"--{name.replace('_', '-')} must be >= 0")
    for name in ("target_gfa", "target_freq_bins", "target_summary"):
        path = getattr(args, name)
        if path and not path.exists():
            raise SystemExit(f"--{name.replace('_', '-')} does not exist: {path}")
    if args.sequence_files and args.sequence_list:
        raise SystemExit("--sequence-files and --sequence-list are mutually exclusive")
    if args.mode == "crush-only":
        if not args.input_gfa:
            raise SystemExit("--input-gfa is required for --mode crush-only")
    else:
        if not args.index:
            raise SystemExit("--index is required for --mode full-query")
        if not args.target_range and not args.target_bed:
            raise SystemExit("--target-range or --target-bed is required for --mode full-query")
        if args.expected_path_count is None:
            raise SystemExit("--expected-path-count is required for --mode full-query")
        if not args.sequence_files and not args.sequence_list:
            raise SystemExit("--sequence-files or --sequence-list is required for --mode full-query")


def utc_now() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def ensure_study(args: argparse.Namespace, dimensions: list[Dimension], optimizer: Optimizer) -> None:
    args.study_dir.mkdir(parents=True, exist_ok=True)
    meta_path = args.study_dir / "study.json"
    meta = {
        "created_or_updated_at": utc_now(),
        "mode": args.mode,
        "method_family": args.method_family,
        "optimizer": optimizer.name,
        "cma_available": HAVE_CMA,
        "cma_import_error": "" if HAVE_CMA else CMA_IMPORT_ERROR,
        "seed": args.seed,
        "population_size": args.population_size,
        "dimensions": [
            {
                "name": dim.name,
                "description": dim.description,
                "range": dim.range_doc,
            }
            for dim in dimensions
        ],
        "objective": objective_description(args),
        "target_profile": target_profile_json(getattr(args, "target_profile", None)),
    }
    meta_path.write_text(json.dumps(meta, indent=2, sort_keys=True) + "\n")


def objective_description(args: argparse.Namespace) -> dict[str, Any]:
    if getattr(args, "target_profile", None) is not None:
        return {
            "formula": (
                "target_frequency_weight*(target_frequency_tv + target_frequency_emd) + "
                "target_length_weight*target_length_distance + "
                "target_total_bp_weight*max(0,total_segment_bp/target_total_bp - 1) + "
                "target_whitespace_weight*legacy_path_whitespace_penalty + "
                "duplicate_sequence_frac + local_repeat_context_nodes/segments"
            ),
            "target_frequency_weight": args.target_frequency_weight,
            "target_length_weight": args.target_length_weight,
            "target_total_bp_weight": args.target_total_bp_weight,
            "target_whitespace_weight": args.target_whitespace_weight,
            "invalid_trials": INVALID_OBJECTIVE,
        }
    return {
        "formula": (
            "total_segment_bp/bp_scale + 2*(singleton_bp/total_segment_bp) + "
            "10*(path_white_space_bp_p99/total_segment_bp) + "
            "5*(path_white_space_bp_max/total_segment_bp) + "
            "(path_white_space_bridges_ge_threshold/paths) + duplicate_sequence_frac + "
            "(local_repeat_context_nodes/segments) - "
            "(node_coverage_bp_weighted_mean/paths)"
        ),
        "invalid_trials": INVALID_OBJECTIVE,
    }


def load_completed_trials(log_path: Path) -> list[dict[str, Any]]:
    if not log_path.exists():
        return []
    rows: list[dict[str, Any]] = []
    with log_path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            try:
                row = json.loads(line)
            except json.JSONDecodeError:
                continue
            if row.get("event") == "trial":
                rows.append(row)
    return rows


def replay_completed_generations(
    optimizer: Optimizer, completed: list[dict[str, Any]], population_size: int
) -> int:
    by_generation: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in completed:
        by_generation[int(row.get("generation", -1))].append(row)

    next_generation = 0
    for generation in sorted(by_generation):
        rows = sorted(by_generation[generation], key=lambda r: int(r.get("candidate_index", 0)))
        if generation != next_generation:
            break
        if len(rows) < population_size:
            break
        rows = rows[:population_size]
        vectors = [list(map(float, row["vector"])) for row in rows]
        objectives = [float(row["tell_objective"]) for row in rows]
        optimizer.replay_generation(vectors, objectives)
        next_generation += 1
    return next_generation


def next_trial_id(completed: list[dict[str, Any]]) -> int:
    if not completed:
        return 1
    return max(int(row.get("trial_id", 0)) for row in completed) + 1


def build_trial_plan(
    trial_id: int,
    generation: int,
    candidate_index: int,
    vector: list[float],
    decoded: dict[str, Any],
    args: argparse.Namespace,
) -> TrialPlan:
    trial_dir = args.study_dir / f"trial-{trial_id:06d}"
    output_gfa = trial_dir / "output.gfa"
    report_tsv = trial_dir / "graph-report.tsv"
    path_validation = trial_dir / "path-validation.json"
    if args.mode == "crush-only":
        primary = build_crush_command(args, decoded, output_gfa)
    else:
        primary = build_query_command(args, decoded, trial_dir)
    graph_report = [
        args.impg_bin,
        "graph-report",
        "-g",
        str(output_gfa),
        "-o",
        str(report_tsv),
        "--format",
        "tsv",
        "--top",
        str(args.report_top),
    ]
    return TrialPlan(
        trial_id=trial_id,
        generation=generation,
        candidate_index=candidate_index,
        vector=vector,
        decoded=decoded,
        trial_dir=trial_dir,
        primary_command=primary,
        graph_report_command=graph_report,
        output_gfa=output_gfa,
        graph_report=report_tsv,
        path_validation=path_validation,
    )


def build_crush_command(args: argparse.Namespace, decoded: dict[str, Any], output_gfa: Path) -> list[str]:
    cmd = [
        args.impg_bin,
        "crush",
        "-g",
        str(args.input_gfa),
        "-o",
        str(output_gfa),
        "-t",
        str(args.threads),
    ]
    for key, value in crush_stage_params(decoded, args):
        cli_key = "sweepga-kmer-frequency" if key == "kmer-frequency" else key
        cmd.extend([f"--{cli_key}", str(value)])
    cmd.extend(args.crush_extra_arg)
    return cmd


def build_query_command(args: argparse.Namespace, decoded: dict[str, Any], trial_dir: Path) -> list[str]:
    engine = build_engine_string(decoded, args)
    prefix = trial_dir / "output"
    cmd = [args.impg_bin, "query", "-a", str(args.index), "-d", str(args.merge_distance)]
    if args.target_range:
        cmd.extend(["-r", args.target_range])
    else:
        cmd.extend(["-b", str(args.target_bed)])
    if args.sequence_files:
        cmd.append("--sequence-files")
        cmd.extend(map(str, args.sequence_files))
    if args.sequence_list:
        cmd.extend(["--sequence-list", str(args.sequence_list)])
    cmd.extend(["-o", engine, "-O", str(prefix), "-t", str(args.threads)])
    cmd.extend(args.query_extra_arg)
    return cmd


def command_script(plan: TrialPlan) -> str:
    return (
        "#!/bin/sh\n"
        "set -eu\n\n"
        "# Primary impg command\n"
        f"{shlex.join(plan.primary_command)}\n\n"
        "# External graph-report scoring command\n"
        f"{shlex.join(plan.graph_report_command)}\n"
    )


def run_trial(plan: TrialPlan, args: argparse.Namespace, objective_scale: float | None, optimizer_name: str) -> TrialResult:
    started_at = utc_now()
    plan.trial_dir.mkdir(parents=True, exist_ok=True)
    (plan.trial_dir / "decoded.json").write_text(json.dumps(plan.decoded, indent=2, sort_keys=True) + "\n")
    command_path = plan.trial_dir / "command.sh"
    command_path.write_text(command_script(plan))
    command_path.chmod(0o755)

    invalid_reasons: list[str] = []
    graph_report_run: TimedCommandResult | None = None
    metrics: dict[str, Any] = {}
    objective_terms: dict[str, float] = {}
    raw_objective: float | None = None

    primary = run_timed(
        plan.primary_command,
        plan.trial_dir / "stdout.log",
        plan.trial_dir / "stderr.log",
        plan.trial_dir / "time.txt",
        args.timeout,
    )
    if primary.timeout:
        invalid_reasons.append("primary_timeout")
    if primary.returncode != 0:
        invalid_reasons.append(f"primary_exit_{primary.returncode}")

    if args.mode == "full-query":
        discovered_outputs = discover_query_output_gfas(plan.trial_dir)
        if len(discovered_outputs) > 1:
            invalid_reasons.append(f"multiple_output_gfa:{len(discovered_outputs)}")
        if discovered_outputs:
            discovered = discovered_outputs[-1]
            if discovered != plan.output_gfa:
                shutil.copyfile(discovered, plan.output_gfa)
        else:
            invalid_reasons.append("missing_output_gfa")
    elif not plan.output_gfa.exists():
        invalid_reasons.append("missing_output_gfa")

    if plan.output_gfa.exists():
        graph_report_run = run_timed(
            plan.graph_report_command,
            plan.trial_dir / "graph-report.stdout.log",
            plan.trial_dir / "graph-report.stderr.log",
            plan.trial_dir / "graph-report.time.txt",
            args.timeout,
        )
        if graph_report_run.timeout:
            invalid_reasons.append("graph_report_timeout")
        if graph_report_run.returncode != 0:
            invalid_reasons.append(f"graph_report_exit_{graph_report_run.returncode}")
        if plan.graph_report.exists() and graph_report_run.returncode == 0:
            try:
                metrics = parse_report_tsv(plan.graph_report)
            except Exception as exc:
                invalid_reasons.append(f"graph_report_parse_error:{exc}")
    else:
        graph_report_run = None

    path_validation = validate_paths(plan, args, metrics)
    plan.path_validation.write_text(json.dumps(path_validation, indent=2, sort_keys=True) + "\n")
    if not path_validation.get("ok", False):
        invalid_reasons.extend(path_validation.get("reasons", ["path_validation_failed"]))

    if metrics and not invalid_reasons:
        try:
            raw_objective, objective_terms = compute_objective(metrics, args, objective_scale, plan.output_gfa)
        except Exception as exc:
            invalid_reasons.append(f"objective_error:{exc}")

    tell_objective = raw_objective if raw_objective is not None else INVALID_OBJECTIVE + plan.trial_id
    result = TrialResult(
        trial_id=plan.trial_id,
        generation=plan.generation,
        candidate_index=plan.candidate_index,
        vector=plan.vector,
        decoded=plan.decoded,
        valid=not invalid_reasons,
        invalid_reasons=invalid_reasons,
        raw_objective=raw_objective,
        tell_objective=tell_objective,
        metrics=metrics,
        objective_terms=objective_terms,
        output_gfa=str(plan.output_gfa),
        graph_report=str(plan.graph_report),
        path_validation=str(plan.path_validation),
        primary=timed_result_json(primary),
        graph_report_run=timed_result_json(graph_report_run) if graph_report_run else None,
        optimizer=optimizer_name,
        started_at=started_at,
        finished_at=utc_now(),
    )
    (plan.trial_dir / "trial.json").write_text(
        json.dumps(dataclasses.asdict(result), indent=2, sort_keys=True) + "\n"
    )
    return result


def run_timed(
    command: list[str], stdout_path: Path, stderr_path: Path, time_path: Path, timeout: int
) -> TimedCommandResult:
    time_bin = shutil.which("time")
    timed_command = command
    if time_bin:
        timed_command = [time_bin, "-v", "-o", str(time_path), *command]
    start = time.monotonic()
    timeout_hit = False
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        try:
            proc = subprocess.Popen(
                timed_command,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
            )
        except OSError as exc:
            stderr.write(f"failed to start command: {exc}\n".encode())
            wall = time.monotonic() - start
            return TimedCommandResult(
                command=command,
                returncode=127,
                timeout=False,
                wall_seconds=wall,
                stdout=str(stdout_path),
                stderr=str(stderr_path),
                time_file=str(time_path),
                max_rss_kb=None,
                time_metrics={},
            )
        try:
            returncode = proc.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            timeout_hit = True
            terminate_process_group(proc)
            returncode = proc.returncode if proc.returncode is not None else -signal.SIGKILL
    wall = time.monotonic() - start
    metrics = parse_time_file(time_path)
    return TimedCommandResult(
        command=command,
        returncode=returncode,
        timeout=timeout_hit,
        wall_seconds=wall,
        stdout=str(stdout_path),
        stderr=str(stderr_path),
        time_file=str(time_path),
        max_rss_kb=parse_int(metrics.get("Maximum resident set size (kbytes)")),
        time_metrics=metrics,
    )


def terminate_process_group(proc: subprocess.Popen[bytes]) -> None:
    try:
        os.killpg(proc.pid, signal.SIGTERM)
    except ProcessLookupError:
        return
    try:
        proc.wait(timeout=5)
    except subprocess.TimeoutExpired:
        try:
            os.killpg(proc.pid, signal.SIGKILL)
        except ProcessLookupError:
            pass
        proc.wait()


def parse_time_file(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    metrics: dict[str, str] = {}
    with path.open(errors="replace") as handle:
        for line in handle:
            if ":" not in line:
                continue
            key, value = line.split(":", 1)
            metrics[key.strip()] = value.strip()
    return metrics


def timed_result_json(result: TimedCommandResult) -> dict[str, Any]:
    return {
        "command": result.command,
        "returncode": result.returncode,
        "timeout": result.timeout,
        "wall_seconds": result.wall_seconds,
        "stdout": result.stdout,
        "stderr": result.stderr,
        "time_file": result.time_file,
        "max_rss_kb": result.max_rss_kb,
        "time_metrics": result.time_metrics,
    }


def discover_query_output_gfas(trial_dir: Path) -> list[Path]:
    exact = trial_dir / "output.gfa"
    if exact.exists():
        return [exact]
    candidates = sorted(trial_dir.rglob("*.gfa"), key=lambda p: (p.stat().st_mtime, str(p)))
    return candidates


def parse_report_tsv(path: Path) -> dict[str, Any]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        row = next(reader)
    parsed: dict[str, Any] = {}
    for key, value in row.items():
        if value is None:
            parsed[key] = value
        elif re.fullmatch(r"-?\d+", value):
            parsed[key] = int(value)
        else:
            try:
                parsed[key] = float(value)
            except ValueError:
                parsed[key] = value
    return parsed


def load_target_profile(args: argparse.Namespace) -> TargetProfile | None:
    if not (args.target_gfa or args.target_freq_bins or args.target_summary):
        return None

    frequency_bins = parse_histogram_bins(DEFAULT_FREQUENCY_BIN_LABELS)
    frequency_distribution: tuple[float, ...] | None = None
    freq_graph: str | None = None
    if args.target_freq_bins:
        frequency_bins, frequency_distribution, freq_graph = read_frequency_bins_tsv(args.target_freq_bins)

    length_bins = parse_histogram_bins(DEFAULT_LENGTH_BIN_LABELS)
    summary = read_target_summary_tsv(args.target_summary) if args.target_summary else {}
    target_shape: GraphShape | None = None
    if args.target_gfa:
        target_shape = graph_shape_from_gfa(args.target_gfa, frequency_bins, length_bins)
        frequency_distribution = target_shape.frequency_distribution
    elif args.target_length_weight > 0.0:
        raise SystemExit("--target-gfa is required when --target-length-weight is greater than zero")

    if frequency_distribution is None:
        raise SystemExit("--target-gfa or --target-freq-bins is required for target-shape scoring")

    total_bp = (
        float(target_shape.total_bp)
        if target_shape is not None
        else optional_positive_float(summary.get("segment_bp") or summary.get("total_segment_bp"))
    )
    path_count = (
        target_shape.path_count
        if target_shape is not None
        else optional_int(summary.get("paths") or summary.get("path_count"))
    )
    bp_weighted_depth_mean = (
        target_shape.bp_weighted_depth_mean
        if target_shape is not None
        else optional_positive_float(
            summary.get("bp_weighted_depth_mean") or summary.get("node_coverage_bp_weighted_mean")
        )
    )
    singleton_bp = (
        float(target_shape.singleton_bp)
        if target_shape is not None
        else optional_positive_float(summary.get("singleton_bp"))
    )
    length_distribution = target_shape.length_distribution if target_shape is not None else None

    if args.target_total_bp_weight > 0.0 and not total_bp:
        raise SystemExit("--target-gfa or --target-summary with segment_bp is required for total-bp target")

    return TargetProfile(
        source_gfa=str(args.target_gfa) if args.target_gfa else None,
        source_freq_bins=str(args.target_freq_bins) if args.target_freq_bins else None,
        source_summary=str(args.target_summary) if args.target_summary else None,
        freq_graph=freq_graph,
        frequency_bins=frequency_bins,
        frequency_distribution=frequency_distribution,
        length_bins=length_bins,
        length_distribution=length_distribution,
        total_bp=total_bp,
        path_count=path_count,
        bp_weighted_depth_mean=bp_weighted_depth_mean,
        singleton_bp=singleton_bp,
    )


def target_profile_json(target: TargetProfile | None) -> dict[str, Any] | None:
    if target is None:
        return None
    return {
        "source_gfa": target.source_gfa,
        "source_freq_bins": target.source_freq_bins,
        "source_summary": target.source_summary,
        "freq_graph": target.freq_graph,
        "frequency_bins": [bin_json(item) for item in target.frequency_bins],
        "frequency_distribution": list(target.frequency_distribution),
        "length_bins": [bin_json(item) for item in target.length_bins],
        "length_distribution": list(target.length_distribution) if target.length_distribution else None,
        "total_bp": target.total_bp,
        "path_count": target.path_count,
        "bp_weighted_depth_mean": target.bp_weighted_depth_mean,
        "singleton_bp": target.singleton_bp,
    }


def bin_json(item: HistogramBin) -> dict[str, Any]:
    return {
        "label": item.label,
        "lower": item.lower,
        "upper": item.upper,
        "center": item.center,
    }


def parse_histogram_bins(labels: Iterable[str]) -> tuple[HistogramBin, ...]:
    parsed: list[tuple[str, float, float | None]] = []
    for raw_label in labels:
        label = str(raw_label).strip()
        if not label:
            continue
        if label.startswith(">="):
            parsed.append((label, float(label[2:]), None))
            continue
        if "-" not in label:
            raise ValueError(f"malformed histogram bin label: {label}")
        low, high = label.split("-", 1)
        parsed.append((label, float(low), float(high)))

    bins: list[HistogramBin] = []
    previous_width = 1.0
    for label, lower, upper in parsed:
        if upper is None:
            center = lower + previous_width / 2.0
        else:
            center = (lower + upper) / 2.0
            previous_width = max(upper - lower, 1.0e-12)
        bins.append(HistogramBin(label=label, lower=lower, upper=upper, center=center))
    if not bins:
        raise ValueError("empty histogram bin list")
    return tuple(bins)


def read_frequency_bins_tsv(path: Path) -> tuple[tuple[HistogramBin, ...], tuple[float, ...], str | None]:
    labels: list[str] = []
    fractions: list[float] = []
    selected_graph: str | None = None
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            graph = row.get("graph") or ""
            if selected_graph is None:
                selected_graph = graph
            if graph != selected_graph:
                continue
            label = row.get("freq_bin")
            if label is None:
                raise ValueError(f"{path} is missing a freq_bin column")
            labels.append(label)
            fractions.append(positive_float(row.get("bp_frac"), 0.0))
    if not labels:
        raise ValueError(f"{path} did not contain any target frequency bins")
    return parse_histogram_bins(labels), tuple(normalize_distribution(fractions)), selected_graph


def read_target_summary_tsv(path: Path) -> dict[str, Any]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        try:
            row = next(reader)
        except StopIteration:
            return {}
    return dict(row)


def graph_shape_from_gfa(
    path: Path, frequency_bins: tuple[HistogramBin, ...], length_bins: tuple[HistogramBin, ...]
) -> GraphShape:
    segment_lengths: dict[str, int] = {}
    occurrence_depth: dict[str, int] = defaultdict(int)
    path_count = 0
    path_steps = 0

    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == "S":
                segment_lengths[fields[1]] = gfa_segment_length(fields)
            elif len(fields) >= 3 and fields[0] == "P":
                path_count += 1
                steps = parse_p_steps(fields[2])
                path_steps += len(steps)
                for node, _orientation in steps:
                    occurrence_depth[node] += 1
            elif len(fields) >= 7 and fields[0] == "W":
                path_count += 1
                steps = parse_w_steps(fields[6])
                path_steps += len(steps)
                for node, _orientation in steps:
                    occurrence_depth[node] += 1

    total_bp = sum(segment_lengths.values())
    path_denominator = max(1, path_count)
    frequency_bp = [0.0] * len(frequency_bins)
    length_bp = [0.0] * len(length_bins)
    singleton_bp = 0
    weighted_depth_sum = 0.0

    for node, length in segment_lengths.items():
        depth = occurrence_depth.get(node, 0)
        weighted_depth_sum += depth * length
        if depth == 1:
            singleton_bp += length
        frequency = depth / path_denominator
        frequency_bp[histogram_bin_index(frequency_bins, frequency)] += length
        length_bp[histogram_bin_index(length_bins, float(length))] += length

    bp_weighted_depth_mean = weighted_depth_sum / total_bp if total_bp > 0 else 0.0
    return GraphShape(
        path_count=path_count,
        segment_count=len(segment_lengths),
        total_bp=total_bp,
        path_steps=path_steps,
        bp_weighted_depth_mean=bp_weighted_depth_mean,
        singleton_bp=singleton_bp,
        frequency_distribution=tuple(normalize_distribution(frequency_bp)),
        length_distribution=tuple(normalize_distribution(length_bp)),
    )


def gfa_segment_length(fields: list[str]) -> int:
    sequence = fields[2]
    if sequence != "*":
        return len(sequence)
    for tag in fields[3:]:
        parts = tag.split(":", 2)
        if len(parts) == 3 and parts[0] == "LN":
            try:
                return max(0, int(parts[2]))
            except ValueError:
                return 0
    return 0


def histogram_bin_index(bins: tuple[HistogramBin, ...], value: float) -> int:
    if value < bins[0].lower:
        return 0
    for idx, item in enumerate(bins):
        if item.upper is None:
            if value >= item.lower:
                return idx
        elif item.lower <= value < item.upper:
            return idx
    return len(bins) - 1


def normalize_distribution(values: Iterable[float]) -> list[float]:
    items = [max(0.0, float(value)) for value in values]
    total = sum(items)
    if total <= 0.0:
        return [0.0 for _ in items]
    return [value / total for value in items]


def total_variation(left: Iterable[float], right: Iterable[float]) -> float:
    return 0.5 * sum(abs(a - b) for a, b in zip(left, right))


def ordered_distribution_distance(
    left: Iterable[float],
    right: Iterable[float],
    bins: tuple[HistogramBin, ...],
    transform: Callable[[float], float] | None = None,
) -> float:
    left_items = list(left)
    right_items = list(right)
    if len(left_items) != len(right_items):
        raise ValueError("distribution lengths differ")
    if len(left_items) <= 1:
        return total_variation(left_items, right_items)

    centers = [item.center for item in bins]
    if transform is not None:
        centers = [transform(center) for center in centers]
    span = max(centers) - min(centers)
    if span <= 0.0:
        span = 1.0

    cdf_delta = 0.0
    distance = 0.0
    for idx in range(len(left_items) - 1):
        cdf_delta += left_items[idx] - right_items[idx]
        distance += abs(cdf_delta) * max(0.0, centers[idx + 1] - centers[idx])
    return distance / span


def compute_objective(
    metrics: dict[str, Any], args: argparse.Namespace, objective_scale: float | None, output_gfa: Path
) -> tuple[float, dict[str, float]]:
    legacy_terms = compute_legacy_objective_terms(metrics, args, objective_scale)
    target_profile = getattr(args, "target_profile", None)
    if target_profile is None:
        return sum(legacy_terms.values()), legacy_terms

    target_terms = compute_target_objective_terms(output_gfa, metrics, args, target_profile, legacy_terms)
    return target_terms["raw_objective_weighted_total"], {**legacy_terms, **target_terms}


def compute_legacy_objective_terms(
    metrics: dict[str, Any], args: argparse.Namespace, objective_scale: float | None
) -> dict[str, float]:
    total_bp = positive_float(metrics.get("total_segment_bp"), 1.0)
    paths = positive_float(metrics.get("paths"), max(1, args.expected_path_count or 1))
    segments = positive_float(metrics.get("segments"), 1.0)
    bp_scale = objective_scale if objective_scale and objective_scale > 0 else 1_000_000.0

    terms = {
        "total_segment_bp_scaled": positive_float(metrics.get("total_segment_bp"), 0.0) / bp_scale,
        "singleton_bp_fraction_x2": 2.0 * positive_float(metrics.get("singleton_bp"), 0.0) / total_bp,
        "path_white_space_bp_p99_fraction_x10": 10.0
        * positive_float(metrics.get("path_white_space_bp_p99"), 0.0)
        / total_bp,
        "path_white_space_bp_max_fraction_x5": 5.0
        * positive_float(metrics.get("path_white_space_bp_max"), 0.0)
        / total_bp,
        "long_white_space_bridges_per_path": positive_float(
            metrics.get("path_white_space_bridges_ge_threshold"), 0.0
        )
        / paths,
        "duplicate_sequence_frac": positive_float(metrics.get("duplicate_sequence_frac"), 0.0),
        "local_repeat_context_nodes_fraction": positive_float(
            metrics.get("local_repeat_context_nodes"), 0.0
        )
        / segments,
        "coverage_reward": -positive_float(metrics.get("node_coverage_bp_weighted_mean"), 0.0)
        / paths,
    }
    return terms


def compute_target_objective_terms(
    output_gfa: Path,
    metrics: dict[str, Any],
    args: argparse.Namespace,
    target: TargetProfile,
    legacy_terms: dict[str, float],
) -> dict[str, float]:
    trial_shape = graph_shape_from_gfa(output_gfa, target.frequency_bins, target.length_bins)
    frequency_tv = total_variation(
        trial_shape.frequency_distribution,
        target.frequency_distribution,
    )
    frequency_emd = ordered_distribution_distance(
        trial_shape.frequency_distribution,
        target.frequency_distribution,
        target.frequency_bins,
    )
    frequency_score = frequency_tv + frequency_emd
    frequency_weighted = args.target_frequency_weight * frequency_score

    if target.length_distribution is None:
        length_tv = 0.0
        length_emd = 0.0
        length_distance = 0.0
    else:
        length_tv = total_variation(trial_shape.length_distribution, target.length_distribution)
        length_emd = ordered_distribution_distance(
            trial_shape.length_distribution,
            target.length_distribution,
            target.length_bins,
            transform=lambda value: math.log2(value + 1.0),
        )
        length_distance = length_tv + length_emd
    length_weighted = args.target_length_weight * length_distance

    target_total_bp = target.total_bp or 0.0
    trial_total_bp = positive_float(metrics.get("total_segment_bp"), float(trial_shape.total_bp))
    if target_total_bp > 0.0:
        total_bp_ratio = trial_total_bp / target_total_bp
        total_bp_excess = max(0.0, trial_total_bp - target_total_bp)
        total_bp_excess_ratio = max(0.0, total_bp_ratio - 1.0)
    else:
        total_bp_ratio = 0.0
        total_bp_excess = 0.0
        total_bp_excess_ratio = 0.0
    total_bp_weighted = args.target_total_bp_weight * total_bp_excess_ratio

    whitespace_penalty = (
        legacy_terms["path_white_space_bp_p99_fraction_x10"]
        + legacy_terms["path_white_space_bp_max_fraction_x5"]
        + legacy_terms["long_white_space_bridges_per_path"]
    )
    whitespace_weighted = args.target_whitespace_weight * whitespace_penalty
    residual_graph_penalty = (
        legacy_terms["duplicate_sequence_frac"] + legacy_terms["local_repeat_context_nodes_fraction"]
    )
    raw_objective = (
        frequency_weighted
        + length_weighted
        + total_bp_weighted
        + whitespace_weighted
        + residual_graph_penalty
    )

    return {
        "target_frequency_tv": frequency_tv,
        "target_frequency_emd": frequency_emd,
        "target_frequency_score": frequency_score,
        "target_frequency_weighted": frequency_weighted,
        "target_length_tv": length_tv,
        "target_length_emd": length_emd,
        "target_length_distance": length_distance,
        "target_length_weighted": length_weighted,
        "target_total_bp": target_total_bp,
        "target_total_bp_ratio": total_bp_ratio,
        "target_total_bp_excess": total_bp_excess,
        "target_total_bp_excess_ratio": total_bp_excess_ratio,
        "target_total_bp_weighted": total_bp_weighted,
        "target_whitespace_penalty": whitespace_penalty,
        "target_whitespace_weighted": whitespace_weighted,
        "target_residual_graph_penalty_weighted": residual_graph_penalty,
        "target_trial_path_count": float(trial_shape.path_count),
        "target_trial_bp_weighted_depth_mean": trial_shape.bp_weighted_depth_mean,
        "target_bp_weighted_depth_mean": target.bp_weighted_depth_mean or 0.0,
        "target_trial_singleton_bp": float(trial_shape.singleton_bp),
        "target_singleton_bp": target.singleton_bp or 0.0,
        "raw_objective_weighted_total": raw_objective,
    }


def positive_float(value: Any, default: float) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    if not math.isfinite(out):
        return default
    return max(0.0, out)


def optional_positive_float(value: Any) -> float | None:
    if value is None:
        return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(out) or out < 0.0:
        return None
    return out


def optional_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(str(value).replace(",", ""))
    except ValueError:
        return None


def parse_int(value: Any) -> int | None:
    if value is None:
        return None
    try:
        return int(str(value).replace(",", ""))
    except ValueError:
        return None


def validate_paths(plan: TrialPlan, args: argparse.Namespace, metrics: dict[str, Any]) -> dict[str, Any]:
    if not plan.output_gfa.exists():
        return {"ok": False, "reasons": ["missing_output_gfa"]}
    try:
        observed = path_sequences(plan.output_gfa)
    except Exception as exc:
        return {"ok": False, "reasons": [f"observed_path_parse_error:{exc}"]}

    if args.mode == "full-query":
        expected = args.expected_path_count
        observed_count = len(observed)
        report_count = int(metrics.get("paths", observed_count)) if metrics else observed_count
        ok = observed_count == expected and report_count == expected
        reasons = []
        if observed_count != expected:
            reasons.append(f"expected_path_count_mismatch:{observed_count}!={expected}")
        if report_count != expected:
            reasons.append(f"report_path_count_mismatch:{report_count}!={expected}")
        return {
            "ok": ok,
            "mode": args.mode,
            "expected_paths": expected,
            "observed_paths": observed_count,
            "report_paths": report_count,
            "reasons": reasons,
        }

    try:
        expected_paths = path_sequences(args.input_gfa)
    except Exception as exc:
        return {"ok": False, "reasons": [f"input_path_parse_error:{exc}"]}
    missing = sorted(set(expected_paths) - set(observed))
    extra = sorted(set(observed) - set(expected_paths))
    mismatches = sorted(
        name for name in set(expected_paths).intersection(observed) if expected_paths[name] != observed[name]
    )
    reasons = []
    if missing:
        reasons.append(f"missing_paths:{len(missing)}")
    if extra:
        reasons.append(f"extra_paths:{len(extra)}")
    if mismatches:
        reasons.append(f"spelling_mismatches:{len(mismatches)}")
    return {
        "ok": not reasons,
        "mode": args.mode,
        "expected_paths": len(expected_paths),
        "observed_paths": len(observed),
        "missing_paths": len(missing),
        "extra_paths": len(extra),
        "spelling_mismatches": len(mismatches),
        "first_missing": missing[0] if missing else None,
        "first_extra": extra[0] if extra else None,
        "first_spelling_mismatch": describe_first_mismatch(
            mismatches[0], expected_paths[mismatches[0]], observed[mismatches[0]]
        )
        if mismatches
        else None,
        "reasons": reasons,
    }


def path_sequences(path: Path) -> dict[str, str]:
    segments: dict[str, str] = {}
    paths: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "S" and len(fields) >= 3:
                segments[fields[1]] = fields[2]
            elif fields[0] == "P" and len(fields) >= 3:
                paths[fields[1]] = spell_steps(parse_p_steps(fields[2]), segments)
            elif fields[0] == "W" and len(fields) >= 7:
                paths[w_line_name(fields)] = spell_steps(parse_w_steps(fields[6]), segments)
    return paths


def parse_p_steps(raw: str) -> list[tuple[str, str]]:
    if raw == "*":
        return []
    steps: list[tuple[str, str]] = []
    for token in raw.split(","):
        match = PATH_STEP_RE.fullmatch(token)
        if not match:
            raise ValueError(f"malformed P step: {token}")
        steps.append((match.group(1), match.group(2)))
    return steps


def parse_w_steps(raw: str) -> list[tuple[str, str]]:
    if raw == "*":
        return []
    return [(match.group(2), "+" if match.group(1) == ">" else "-") for match in WALK_STEP_RE.finditer(raw)]


def spell_steps(steps: list[tuple[str, str]], segments: dict[str, str]) -> str:
    seqs: list[str] = []
    for node, orientation in steps:
        if node not in segments:
            raise ValueError(f"path references unknown segment: {node}")
        seq = segments[node]
        seqs.append(seq if orientation == "+" else reverse_complement(seq))
    return "".join(seqs)


def reverse_complement(seq: str) -> str:
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def w_line_name(fields: list[str]) -> str:
    sample = fields[1] if len(fields) > 1 else "*"
    hap = fields[2] if len(fields) > 2 else "*"
    seq = fields[3] if len(fields) > 3 else "*"
    start = fields[4] if len(fields) > 4 else "*"
    end = fields[5] if len(fields) > 5 else "*"
    return f"{sample}#{hap}#{seq}:{start}-{end}"


def describe_first_mismatch(name: str, expected: str, observed: str) -> dict[str, Any]:
    first = next(
        (idx for idx, pair in enumerate(zip(expected, observed)) if pair[0].upper() != pair[1].upper()),
        min(len(expected), len(observed)),
    )
    return {
        "path": name,
        "offset": first,
        "expected_len": len(expected),
        "observed_len": len(observed),
        "expected_preview": expected[max(0, first - 40) : first + 40],
        "observed_preview": observed[max(0, first - 40) : first + 40],
    }


def append_trial_log(log_path: Path, result: TrialResult, lock: threading.Lock) -> None:
    row = dataclasses.asdict(result)
    row["event"] = "trial"
    with lock:
        with log_path.open("a") as handle:
            handle.write(json.dumps(row, sort_keys=True) + "\n")


def dry_run(args: argparse.Namespace, dimensions: list[Dimension]) -> int:
    vector = [0.5] * len(dimensions)
    decoded = normalize_decoded(decode_dimensions(dimensions, vector), args)
    plan = build_trial_plan(1, 0, 0, vector, decoded, args)
    print(
        json.dumps(
            {
                "mode": args.mode,
                "method_family": args.method_family,
                "dimensions": [
                    {
                        "name": dim.name,
                        "description": dim.description,
                        "range": dim.range_doc,
                    }
                    for dim in dimensions
                ],
                "decoded": decoded,
                "target_profile": target_profile_json(getattr(args, "target_profile", None)),
                "primary_command": plan.primary_command,
                "graph_report_command": plan.graph_report_command,
                "command_sh": command_script(plan),
            },
            indent=2,
            sort_keys=True,
        )
    )
    return 0


def objective_scale_for_study(args: argparse.Namespace) -> float | None:
    if args.objective_bp_scale:
        return args.objective_bp_scale
    if args.mode != "crush-only" or not args.input_gfa:
        return None
    try:
        total_bp = sum(len(seq) for seq in read_segments(args.input_gfa).values() if seq != "*")
    except Exception:
        return None
    return float(total_bp) if total_bp > 0 else None


def read_segments(path: Path) -> dict[str, str]:
    segments: dict[str, str] = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == "S":
                segments[fields[1]] = fields[2]
    return segments


def main() -> int:
    args = parse_args()
    args.target_profile = load_target_profile(args)
    dimensions = build_dimensions(args)
    if args.dry_run:
        return dry_run(args, dimensions)

    optimizer = make_optimizer(len(dimensions), args)
    ensure_study(args, dimensions, optimizer)
    log_path = args.study_dir / "trials.jsonl"
    completed = [] if args.no_resume else load_completed_trials(log_path)
    start_generation = 0 if args.no_resume else replay_completed_generations(
        optimizer, completed, args.population_size
    )
    trial_id = 1 if args.no_resume else next_trial_id(completed)
    completed_count = 0 if args.no_resume else len(completed)
    if completed_count >= args.max_trials:
        print(f"study already has {completed_count} trial(s), max-trials={args.max_trials}", file=sys.stderr)
        return 0

    objective_scale = objective_scale_for_study(args)
    log_lock = threading.Lock()
    generation = start_generation

    print(
        f"optimizer={optimizer.name} dimensions={len(dimensions)} "
        f"completed={completed_count} next_generation={generation}",
        file=sys.stderr,
    )
    if not HAVE_CMA:
        print(f"warning: cma import failed, using deterministic random fallback: {CMA_IMPORT_ERROR}", file=sys.stderr)

    while completed_count < args.max_trials:
        remaining = args.max_trials - completed_count
        batch_size = min(args.population_size, remaining)
        vectors = optimizer.ask(batch_size)
        plans: list[TrialPlan] = []
        for candidate_index, vector in enumerate(vectors):
            decoded = normalize_decoded(decode_dimensions(dimensions, vector), args)
            plans.append(build_trial_plan(trial_id, generation, candidate_index, vector, decoded, args))
            trial_id += 1

        results: list[TrialResult] = []
        with ThreadPoolExecutor(max_workers=args.jobs) as executor:
            future_to_plan = {
                executor.submit(run_trial, plan, args, objective_scale, optimizer.name): plan for plan in plans
            }
            for future in as_completed(future_to_plan):
                result = future.result()
                append_trial_log(log_path, result, log_lock)
                results.append(result)
                status = "valid" if result.valid else "invalid:" + ",".join(result.invalid_reasons)
                print(
                    f"trial={result.trial_id:06d} generation={result.generation} "
                    f"candidate={result.candidate_index} objective={result.tell_objective:.6g} {status}",
                    file=sys.stderr,
                )

        results.sort(key=lambda r: r.candidate_index)
        if len(results) == args.population_size:
            optimizer.tell([r.vector for r in results], [r.tell_objective for r in results])
        else:
            print(
                f"final partial generation has {len(results)} candidate(s); "
                "not telling CMA-ES because no further ask is needed",
                file=sys.stderr,
            )
        completed_count += len(results)
        generation += 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
