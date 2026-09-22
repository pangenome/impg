#!/usr/bin/env python3
"""Run and summarize the C4 run-aware high-frequency mask crush sweep."""

from __future__ import annotations

import argparse
import csv
import json
import math
import os
import re
import shlex
import subprocess
import sys
import time
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable


INDEX = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng"
AGC = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc"
REGION = "GRCh38#0#chr6:31891045-32123783"
MERGE_DISTANCE = "50k"
TARGET_GFA = (
    "/home/erikg/impg/data/local_syng_parameter_sweep_20260529T193452Z/"
    "pggb_control/graphs/C4_GRCh38_53kb.gfa"
)
DEFAULT_OUT_DIR = "/home/erikg/impg/data/c4_highfreq_mask_crush_sweep"
EXPECTED_PATHS = 465

FREQUENCY_BINS = (
    ("0-0.01", 0.0, 0.01, 0.005),
    ("0.01-0.05", 0.01, 0.05, 0.03),
    ("0.05-0.1", 0.05, 0.1, 0.075),
    ("0.1-0.25", 0.1, 0.25, 0.175),
    ("0.25-0.5", 0.25, 0.5, 0.375),
    ("0.5-0.75", 0.5, 0.75, 0.625),
    ("0.75-1", 0.75, 1.0, 0.875),
    ("1-1.5", 1.0, 1.5, 1.25),
    ("1.5-2", 1.5, 2.0, 1.75),
    ("2-3", 2.0, 3.0, 2.5),
    (">=3", 3.0, None, 3.5),
)

LENGTH_BINS = (
    ("0-1", 0.0, 1.0, 0.5),
    ("1-2", 1.0, 2.0, 1.5),
    ("2-4", 2.0, 4.0, 3.0),
    ("4-8", 4.0, 8.0, 6.0),
    ("8-16", 8.0, 16.0, 12.0),
    ("16-32", 16.0, 32.0, 24.0),
    ("32-64", 32.0, 64.0, 48.0),
    ("64-128", 64.0, 128.0, 96.0),
    ("128-256", 128.0, 256.0, 192.0),
    ("256-512", 256.0, 512.0, 384.0),
    ("512-1024", 512.0, 1024.0, 768.0),
    ("1024-2048", 1024.0, 2048.0, 1536.0),
    ("2048-4096", 2048.0, 4096.0, 3072.0),
    ("4096-8192", 4096.0, 8192.0, 6144.0),
    ("8192-16384", 8192.0, 16384.0, 12288.0),
    ("16384-32768", 16384.0, 32768.0, 24576.0),
    ("32768-65536", 32768.0, 65536.0, 49152.0),
    (">=65536", 65536.0, None, 98304.0),
)


@dataclass(frozen=True)
class SweepConfig:
    name: str
    top: float
    max_occ: int
    freq_run: int
    freq_span: int
    note: str
    mask_extra: tuple[str, ...] = ()
    crush_extra: tuple[str, ...] = ()

    def mask_params(self) -> str:
        parts = [f"top={self.top:g}", f"freq-run={self.freq_run}", f"freq-span={self.freq_span}"]
        if self.max_occ > 0:
            parts.append(f"max-occ={self.max_occ}")
        parts.extend(self.mask_extra)
        parts.append("freq-run-aware=true")
        return ",".join(parts)


CONFIGS = (
    SweepConfig("default_top0005_run10_span0", 0.0005, 0, 10, 0, "current run-aware default"),
    SweepConfig("top001_run3_span0", 0.001, 0, 3, 0, "stronger top, permissive run rescue"),
    SweepConfig("top001_run5_span0", 0.001, 0, 5, 0, "stronger top, mid run rescue"),
    SweepConfig("top001_run10_span0", 0.001, 0, 10, 0, "stronger top, default run rescue"),
    SweepConfig("top001_run16_span0", 0.001, 0, 16, 0, "stronger top, strict run rescue"),
    SweepConfig("top001_run10_span500", 0.001, 0, 10, 500, "default run rescue plus 500 bp span"),
    SweepConfig("top001_run10_span1k", 0.001, 0, 10, 1000, "default run rescue plus 1 kb span"),
    SweepConfig("top002_run10_span1k", 0.002, 0, 10, 1000, "2x top with 1 kb span rescue"),
    SweepConfig("maxocc1000_run10_span1k", 0.0, 1000, 10, 1000, "absolute occurrence cap"),
    SweepConfig(
        "unified_top001_run10_span1k",
        0.001,
        0,
        10,
        1000,
        "unified explicit-HF plus spectrum-glue run/span rescue",
    ),
    SweepConfig(
        "unified_maxocc1000_run10_span1k",
        0.0,
        1000,
        10,
        1000,
        "unified absolute cap plus spectrum-glue run/span rescue",
    ),
    SweepConfig(
        "top001_run10_span1k_sharedrun10_seq1k",
        0.001,
        0,
        10,
        1000,
        "apply run/span rescue to generic shared/scaffold-context mask",
        ("min-run=10", "sequence-k=1k"),
    ),
    SweepConfig(
        "top001_run10_span1k_noskip_polish50k",
        0.001,
        0,
        10,
        1000,
        "disable candidate-size skip and raise POASTA polish caps",
        (),
        (
            "max-traversal-len=0",
            "max-median-traversal-len=0",
            "polish-max-traversal-len=50k",
            "polish-max-median-traversal-len=10k",
        ),
    ),
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path(DEFAULT_OUT_DIR))
    parser.add_argument("--impg", default="impg")
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--timeout", type=int, default=0, help="Per-command timeout in seconds; 0 disables")
    parser.add_argument(
        "--configs",
        help="Comma-separated config names to run/summarize; defaults to the full compact matrix",
    )
    parser.add_argument("--skip-runs", action="store_true")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--render-top", type=int, default=3)
    parser.add_argument("--no-upload", action="store_true")
    parser.add_argument("--upload-host", default="erik@hypervolu.me")
    parser.add_argument("--upload-dir", default="www/impg")
    return parser.parse_args()


def selected_configs(args: argparse.Namespace) -> tuple[SweepConfig, ...]:
    if not args.configs:
        return CONFIGS
    wanted = [item.strip() for item in args.configs.split(",") if item.strip()]
    by_name = {config.name: config for config in CONFIGS}
    missing = [name for name in wanted if name not in by_name]
    if missing:
        raise SystemExit(f"unknown --configs value(s): {', '.join(missing)}")
    return tuple(by_name[name] for name in wanted)


def run_logged(
    cmd: list[str],
    stdout_path: Path,
    stderr_path: Path,
    timeout: int,
    env: dict[str, str] | None = None,
) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    with stdout_path.open("w") as stdout, stderr_path.open("w") as stderr:
        try:
            completed = subprocess.run(
                cmd,
                stdout=stdout,
                stderr=stderr,
                text=True,
                timeout=timeout or None,
                env=env,
            )
            return completed.returncode
        except subprocess.TimeoutExpired:
            stderr.write(f"\nTIMEOUT after {timeout}s\n")
            return 124


def timed_command(cmd: list[str]) -> list[str]:
    return ["/usr/bin/time", "-v", *cmd]


def engine_string(config: SweepConfig) -> str:
    crush_params = [
        "method=auto",
        "auto-2tier=true",
        "auto-poasta-max-len=10k",
        "max-traversal-len=10k",
        "max-median-traversal-len=1k",
        "min-traversal-len=5k",
        "max-rounds=until-done",
        "seqwish-k=311",
        "max-pair-alignments=0",
        "max-paf-bytes=0",
        "no-filter=true",
        "polish-method=poasta",
        "polish-rounds=until-done",
        "polish-max-traversal-len=10k",
        "polish-max-median-traversal-len=1k",
    ]
    crush_params.extend(config.crush_extra)
    return f"gfa:syng:mask,{config.mask_params()}:crush,{','.join(override_params(crush_params))}:nosort"


def override_params(params: list[str]) -> list[str]:
    order: list[str] = []
    by_key: dict[str, str] = {}
    for param in params:
        key = param.split("=", 1)[0]
        if key not in by_key:
            order.append(key)
        by_key[key] = param
    return [by_key[key] for key in order]


def query_command(args: argparse.Namespace, config: SweepConfig, prefix: Path) -> list[str]:
    return [
        args.impg,
        "query",
        "-t",
        str(args.threads),
        "-a",
        INDEX,
        "-r",
        REGION,
        "--sequence-files",
        AGC,
        "-d",
        MERGE_DISTANCE,
        "-o",
        engine_string(config),
        "-O",
        str(prefix),
        "-v",
        "1",
    ]


def graph_report_command(args: argparse.Namespace, gfa: Path, output: Path) -> list[str]:
    return [
        args.impg,
        "graph-report",
        "-g",
        str(gfa),
        "-o",
        str(output),
        "--format",
        "tsv",
        "--top",
        "20",
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]


def sort_command(args: argparse.Namespace, in_gfa: Path, out_gfa: Path) -> list[str]:
    return ["gfasort", "-i", str(in_gfa), "-o", str(out_gfa), "-p", "Ygs", "-t", str(args.threads)]


def read_report_tsv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return dict(rows[0]) if rows else {}


def parse_time_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    for line in path.read_text(errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            match = re.search(r"\):\s*(.+)$", line)
            values["wall"] = match.group(1).strip() if match else line.rsplit(":", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            values["max_rss_kb"] = line.rsplit(":", 1)[1].strip()
        elif "Exit status" in line:
            values["exit_status"] = line.rsplit(":", 1)[1].strip()
    return values


def parse_mask_log(path: Path) -> dict[str, str]:
    text = path.read_text(errors="replace") if path.exists() else ""
    out = {
        "hf_selected_nodes": "",
        "hf_total_nodes": "",
        "hf_local_max_occ": "",
        "hf_selected_min_occ": "",
        "hf_selected_max_occ": "",
        "hf_occurrences": "",
        "hf_rescued_occurrences": "",
        "hf_run_supported": "",
        "hf_sequence_supported": "",
        "hf_private_split_occurrences": "",
        "spectrum_glue_nodes": "",
        "spectrum_glue_occurrences": "",
        "spectrum_glue_rescued_occurrences": "",
        "spectrum_glue_run_supported": "",
        "spectrum_glue_sequence_supported": "",
        "spectrum_glue_private_split_occurrences": "",
        "spectrum_glue_scaffold_candidates": "",
        "spectrum_glue_dense_signatures": "",
        "scaffold_weak_occurrences": "",
        "scaffold_split_occurrences": "",
        "scaffold_min_run": "",
        "scaffold_supported": "",
        "scaffold_candidates": "",
        "scaffold_dense_signatures": "",
        "scaffold_sequence_k": "",
        "scaffold_sequence_supported": "",
    }
    selected_re = re.search(
        r"selected\s+(\d+)\s+/\s+(\d+)\s+high-frequency local syncmer node\(s\).*?"
        r"local max=(\d+), selected min/max=(\d+)/(\d+)",
        text,
    )
    if selected_re:
        out.update(
            {
                "hf_selected_nodes": selected_re.group(1),
                "hf_total_nodes": selected_re.group(2),
                "hf_local_max_occ": selected_re.group(3),
                "hf_selected_min_occ": selected_re.group(4),
                "hf_selected_max_occ": selected_re.group(5),
            }
        )
    occurrence_re = re.search(
        r"high-frequency occurrence mask selected\s+(\d+)\s+node\(s\) covering\s+(\d+)\s+"
        r"local occurrence\(s\); rescued/supported\s+(\d+)\s+occurrence\(s\)\s+"
        r"\(run-supported=(\d+), sequence-supported=(\d+).*?\), privately split\s+(\d+)\s+",
        text,
    )
    if occurrence_re:
        out.update(
            {
                "hf_selected_nodes": occurrence_re.group(1),
                "hf_occurrences": occurrence_re.group(2),
                "hf_rescued_occurrences": occurrence_re.group(3),
                "hf_run_supported": occurrence_re.group(4),
                "hf_sequence_supported": occurrence_re.group(5),
                "hf_private_split_occurrences": occurrence_re.group(6),
            }
        )
    if "high-frequency occurrence mask disabled" in text:
        out["hf_occurrence_mask"] = "disabled"
    elif "high-frequency occurrence mask selected" in text:
        out["hf_occurrence_mask"] = "run-aware"
    else:
        out["hf_occurrence_mask"] = ""
    spectrum_re = re.search(
        r"spectrum-selected\s+(\d+)\s+dispersed high-copy scaffold-glue node\(s\) covering\s+"
        r"(\d+)\s+local occurrence\(s\)",
        text,
    )
    if spectrum_re:
        out["spectrum_glue_nodes"] = spectrum_re.group(1)
        out["spectrum_glue_occurrences"] = spectrum_re.group(2)
    scaffold_re = re.search(
        r"scaffold-context evaluated\s+\d+\s+local shared syncmer occurrence\(s\).*?"
        r"\(weak=(\d+), split=(\d+), min_run=(\d+), spectrum-glue-filtered-nodes=(\d+), "
        r"spectrum-glue-filtered-occurrences=(\d+), scaffold-supported=(\d+), "
        r"scaffold-candidates=(\d+), dense-signatures=(\d+), sequence_k=(\d+), "
        r"sequence-supported=(\d+)(?:, [^)]*)?\)",
        text,
    )
    if scaffold_re:
        out.update(
            {
                "scaffold_weak_occurrences": scaffold_re.group(1),
                "scaffold_split_occurrences": scaffold_re.group(2),
                "scaffold_min_run": scaffold_re.group(3),
                "spectrum_glue_nodes": scaffold_re.group(4),
                "spectrum_glue_occurrences": scaffold_re.group(5),
                "scaffold_supported": scaffold_re.group(6),
                "scaffold_candidates": scaffold_re.group(7),
                "scaffold_dense_signatures": scaffold_re.group(8),
                "scaffold_sequence_k": scaffold_re.group(9),
                "scaffold_sequence_supported": scaffold_re.group(10),
            }
        )
    spectrum_policy_re = re.search(
        r"spectrum-glue-rescued=(\d+), spectrum-glue-run-supported=(\d+), "
        r"spectrum-glue-sequence-supported=(\d+), spectrum-glue-private-split=(\d+), "
        r"spectrum-glue-scaffold-candidates=(\d+), spectrum-glue-dense-signatures=(\d+)",
        text,
    )
    if spectrum_policy_re:
        out.update(
            {
                "spectrum_glue_rescued_occurrences": spectrum_policy_re.group(1),
                "spectrum_glue_run_supported": spectrum_policy_re.group(2),
                "spectrum_glue_sequence_supported": spectrum_policy_re.group(3),
                "spectrum_glue_private_split_occurrences": spectrum_policy_re.group(4),
                "spectrum_glue_scaffold_candidates": spectrum_policy_re.group(5),
                "spectrum_glue_dense_signatures": spectrum_policy_re.group(6),
            }
        )
    return out


def gfa_segment_length(fields: list[str]) -> int:
    if len(fields) >= 3 and fields[2] != "*":
        return len(fields[2])
    for tag in fields[3:]:
        parts = tag.split(":", 2)
        if len(parts) == 3 and parts[0] == "LN":
            return int(parts[2])
    return 0


def bin_index(bins: tuple[tuple[str, float, float | None, float], ...], value: float) -> int:
    for idx, (_label, lower, upper, _center) in enumerate(bins):
        if upper is None:
            if value >= lower:
                return idx
        elif lower <= value < upper:
            return idx
    return len(bins) - 1


def normalize(values: Iterable[float]) -> tuple[float, ...]:
    items = [max(0.0, float(value)) for value in values]
    total = sum(items)
    if total <= 0.0:
        return tuple(0.0 for _ in items)
    return tuple(value / total for value in items)


def parse_p_steps(value: str) -> list[tuple[str, str]]:
    steps = []
    if value == "*":
        return steps
    for raw in value.split(","):
        if not raw:
            continue
        steps.append((raw[:-1], raw[-1]))
    return steps


def parse_w_steps(value: str) -> list[tuple[str, str]]:
    out = []
    idx = 0
    while idx < len(value):
        orient = value[idx]
        idx += 1
        start = idx
        while idx < len(value) and value[idx] not in "><":
            idx += 1
        out.append((value[start:idx], "+" if orient == ">" else "-"))
    return out


def graph_shape(path: Path) -> dict[str, object]:
    lengths: dict[str, int] = {}
    depths: dict[str, int] = defaultdict(int)
    path_count = 0
    path_steps = 0
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3 and fields[0] == "S":
                lengths[fields[1]] = gfa_segment_length(fields)
            elif len(fields) >= 3 and fields[0] == "P":
                path_count += 1
                steps = parse_p_steps(fields[2])
                path_steps += len(steps)
                for node, _orient in steps:
                    depths[node] += 1
            elif len(fields) >= 7 and fields[0] == "W":
                path_count += 1
                steps = parse_w_steps(fields[6])
                path_steps += len(steps)
                for node, _orient in steps:
                    depths[node] += 1
    total_bp = sum(lengths.values())
    freq_bp = [0.0 for _ in FREQUENCY_BINS]
    length_bp = [0.0 for _ in LENGTH_BINS]
    singleton_bp = 0
    weighted_depth = 0.0
    denom = max(1, path_count)
    for node, length in lengths.items():
        depth = depths.get(node, 0)
        weighted_depth += depth * length
        if depth == 1:
            singleton_bp += length
        freq_bp[bin_index(FREQUENCY_BINS, depth / denom)] += length
        length_bp[bin_index(LENGTH_BINS, float(length))] += length
    return {
        "paths": path_count,
        "segments": len(lengths),
        "path_steps": path_steps,
        "total_segment_bp": total_bp,
        "bp_weighted_depth": weighted_depth / total_bp if total_bp else 0.0,
        "singleton_bp": singleton_bp,
        "frequency_distribution": normalize(freq_bp),
        "length_distribution": normalize(length_bp),
    }


def total_variation(left: Iterable[float], right: Iterable[float]) -> float:
    return 0.5 * sum(abs(a - b) for a, b in zip(left, right))


def ordered_distance(
    left: Iterable[float],
    right: Iterable[float],
    bins: tuple[tuple[str, float, float | None, float], ...],
    log_centers: bool = False,
) -> float:
    left_items = list(left)
    right_items = list(right)
    centers = [math.log2(center + 1.0) if log_centers else center for *_rest, center in bins]
    span = max(centers) - min(centers)
    if span <= 0.0:
        span = 1.0
    cdf_delta = 0.0
    distance = 0.0
    for idx in range(len(left_items) - 1):
        cdf_delta += left_items[idx] - right_items[idx]
        distance += abs(cdf_delta) * max(0.0, centers[idx + 1] - centers[idx])
    return distance / span


def target_metrics(shape: dict[str, object], target: dict[str, object], report: dict[str, str]) -> dict[str, str]:
    freq_dist = shape["frequency_distribution"]
    target_freq = target["frequency_distribution"]
    length_dist = shape["length_distribution"]
    target_length = target["length_distribution"]
    freq_tv = total_variation(freq_dist, target_freq)
    freq_emd = ordered_distance(freq_dist, target_freq, FREQUENCY_BINS)
    length_tv = total_variation(length_dist, target_length)
    length_emd = ordered_distance(length_dist, target_length, LENGTH_BINS, log_centers=True)
    total_bp = float(report.get("total_segment_bp") or shape["total_segment_bp"] or 0.0)
    target_bp = float(target["total_segment_bp"])
    bp_ratio = total_bp / target_bp if target_bp else 0.0
    ws_p99 = positive_float(report.get("path_white_space_bp_p99"), 0.0)
    ws_max = positive_float(report.get("path_white_space_bp_max"), 0.0)
    long_ws = positive_float(report.get("path_white_space_bridges_ge_threshold"), 0.0)
    paths = max(1.0, positive_float(report.get("paths"), float(shape["paths"])))
    whitespace_penalty = (10.0 * ws_p99 / max(1.0, total_bp)) + (5.0 * ws_max / max(1.0, total_bp)) + (
        long_ws / paths
    )
    score = (
        1000.0 * (freq_tv + freq_emd)
        + 100.0 * (length_tv + length_emd)
        + 25.0 * max(0.0, bp_ratio - 1.0)
        + whitespace_penalty
    )
    return {
        "target_frequency_tv": f"{freq_tv:.6f}",
        "target_frequency_emd": f"{freq_emd:.6f}",
        "target_length_distance": f"{(length_tv + length_emd):.6f}",
        "target_total_bp_ratio": f"{bp_ratio:.6f}",
        "target_shape_score": f"{score:.3f}",
    }


def positive_float(value: object, default: float) -> float:
    try:
        out = float(str(value).replace(",", ""))
    except (TypeError, ValueError):
        return default
    return out if math.isfinite(out) and out >= 0.0 else default


def summarize_run(args: argparse.Namespace, config: SweepConfig, target: dict[str, object]) -> dict[str, str]:
    run_dir = args.out_dir / config.name
    gfa = run_dir / "run.nosort.gfa"
    sorted_gfa = run_dir / "run.Ygs.gfa"
    report = run_dir / "graph-report.tsv"
    metrics = read_report_tsv(report)
    shape = graph_shape(gfa) if gfa.exists() else {}
    mask = parse_mask_log(run_dir / "query.stderr.log")
    time_values = parse_time_log(run_dir / "query.stderr.log")
    sort_time = parse_time_log(run_dir / "sort.stderr.log")
    report_time = parse_time_log(run_dir / "graph-report.stderr.log")
    target_values = target_metrics(shape, target, metrics) if shape and metrics else {}
    paths = str(metrics.get("paths") or shape.get("paths") or "")
    path_count_ok = str(paths == str(EXPECTED_PATHS)).lower()
    row = {
        "name": config.name,
        "note": config.note,
        "top": f"{config.top:g}",
        "max_occ": str(config.max_occ or ""),
        "freq_run": str(config.freq_run),
        "freq_span": str(config.freq_span),
        "gfa": str(gfa),
        "sorted_gfa": str(sorted_gfa),
        "exit_status": time_values.get("exit_status", ""),
        "paths": paths,
        "path_count_ok": path_count_ok,
        "segments": str(metrics.get("segments") or shape.get("segments") or ""),
        "links": str(metrics.get("links") or ""),
        "total_segment_bp": str(metrics.get("total_segment_bp") or shape.get("total_segment_bp") or ""),
        "bp_weighted_depth": str(
            metrics.get("node_coverage_bp_weighted_mean") or shape.get("bp_weighted_depth") or ""
        ),
        "singleton_bp": str(metrics.get("singleton_bp") or shape.get("singleton_bp") or ""),
        "path_white_space_bp_p99": str(metrics.get("path_white_space_bp_p99") or ""),
        "path_white_space_bp_max": str(metrics.get("path_white_space_bp_max") or ""),
        "path_white_space_bridges_ge_threshold": str(
            metrics.get("path_white_space_bridges_ge_threshold") or ""
        ),
        "wall": time_values.get("wall", ""),
        "max_rss_kb": time_values.get("max_rss_kb", ""),
        "sort_wall": sort_time.get("wall", ""),
        "sort_max_rss_kb": sort_time.get("max_rss_kb", ""),
        "graph_report_wall": report_time.get("wall", ""),
        "graph_report_max_rss_kb": report_time.get("max_rss_kb", ""),
    }
    row.update(mask)
    row.update(target_values)
    return row


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def run_sweep(args: argparse.Namespace) -> None:
    args.out_dir.mkdir(parents=True, exist_ok=True)
    configs = selected_configs(args)
    target = graph_shape(Path(TARGET_GFA))
    write_json(args.out_dir / "target-shape.json", clean_shape_for_json(target))
    target_report = args.out_dir / "pggb-target.graph-report.tsv"
    if args.force or not target_report.exists():
        cmd = graph_report_command(args, Path(TARGET_GFA), target_report)
        run_logged(
            timed_command(cmd),
            args.out_dir / "pggb-target.graph-report.stdout.log",
            args.out_dir / "pggb-target.graph-report.stderr.log",
            args.timeout,
        )
    for config in configs:
        run_dir = args.out_dir / config.name
        run_dir.mkdir(parents=True, exist_ok=True)
        prefix = run_dir / "run.nosort"
        nosort_gfa = run_dir / "run.nosort.gfa"
        command = query_command(args, config, prefix)
        (run_dir / "engine.txt").write_text(engine_string(config) + "\n")
        (run_dir / "command.sh").write_text("#!/bin/sh\nset -eu\n" + shlex.join(command) + "\n")
        os.chmod(run_dir / "command.sh", 0o755)
        if not args.skip_runs and (args.force or not nosort_gfa.exists()):
            print(f"[sweep] running {config.name}", flush=True)
            rc = run_logged(
                timed_command(command),
                run_dir / "query.stdout.log",
                run_dir / "query.stderr.log",
                args.timeout,
            )
            print(f"[sweep] {config.name} exit={rc}", flush=True)
        if nosort_gfa.exists() and (args.force or not (run_dir / "graph-report.tsv").exists()):
            print(f"[sweep] graph-report {config.name}", flush=True)
            rc = run_logged(
                timed_command(graph_report_command(args, nosort_gfa, run_dir / "graph-report.tsv")),
                run_dir / "graph-report.stdout.log",
                run_dir / "graph-report.stderr.log",
                args.timeout,
            )
            print(f"[sweep] graph-report {config.name} exit={rc}", flush=True)
    rows = [summarize_run(args, config, target) for config in configs]
    write_tsv(args.out_dir / "metrics.tsv", rows)
    write_json(args.out_dir / "configs.json", [config.__dict__ for config in configs])
    render_and_upload(args, rows)


def clean_shape_for_json(shape: dict[str, object]) -> dict[str, object]:
    out = dict(shape)
    out["frequency_distribution"] = list(out["frequency_distribution"])
    out["length_distribution"] = list(out["length_distribution"])
    return out


def render_and_upload(args: argparse.Namespace, rows: list[dict[str, str]]) -> None:
    valid_rows = [
        row
        for row in rows
        if row.get("path_count_ok") == "true" and row.get("target_shape_score") and Path(row["gfa"]).exists()
    ]
    valid_rows.sort(key=lambda row: float(row["target_shape_score"]))
    rendered = []
    for row in valid_rows[: max(0, args.render_top)]:
        name = row["name"]
        gfa = Path(row["gfa"])
        sorted_gfa = Path(row["sorted_gfa"])
        if args.force or not sorted_gfa.exists():
            rc = run_logged(
                timed_command(sort_command(args, gfa, sorted_gfa)),
                args.out_dir / name / "sort.stdout.log",
                args.out_dir / name / "sort.stderr.log",
                args.timeout,
            )
            if rc != 0:
                continue
        png = args.out_dir / name / f"{name}.Ygs.gfalook-m.png"
        if args.force or not png.exists():
            cmd = ["gfalook", "-i", str(sorted_gfa), "-o", str(png), "-m", "-x", "2200", "-y", "1200"]
            rc = run_logged(
                timed_command(cmd),
                args.out_dir / name / "gfalook.stdout.log",
                args.out_dir / name / "gfalook.stderr.log",
                args.timeout,
            )
            if rc != 0:
                continue
        url = ""
        if not args.no_upload:
            remote_name = f"c4-highfreq-mask-crush-sweep-{name}.png"
            remote = f"{args.upload_host}:{args.upload_dir}/{remote_name}"
            rc = subprocess.run(["scp", str(png), remote]).returncode
            if rc == 0:
                url = f"https://hypervolu.me/~erik/impg/{remote_name}"
        row["png"] = str(png)
        row["png_url"] = url
        rendered.append(row)
    write_tsv(args.out_dir / "rendered.tsv", rendered)


def main() -> int:
    args = parse_args()
    run_sweep(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
