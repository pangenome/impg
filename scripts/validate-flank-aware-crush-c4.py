#!/usr/bin/env python3
"""Run bounded flank-aware crush validation on sanity fixtures and C4.

The script is intentionally narrow for the WG validation task. It records every
command, keeps stdout/stderr under an output directory, compares exact GFA path
spellings, extracts flank/candidate diagnostics from implementation logs, and
writes summary TSV/JSON artifacts for the evaluation report.
"""

from __future__ import annotations

import argparse
import csv
import datetime as dt
import gzip
import json
import os
import re
import shlex
import subprocess
import sys
import time
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Callable, Iterable


DEFAULT_C4_INPUT = Path(
    "/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa"
)
DEFAULT_CURRENT_BEST = Path(
    "/home/erikg/impg/data/c4_spoa2k_default_selfloop_verify_20260606T065905Z/"
    "sorted/c4.k311.seed.spoa2k.default-selfloop.Ygs.gfa"
)
DEFAULT_CURRENT_BEST_REPORT = Path(
    "/home/erikg/impg/data/c4_spoa2k_default_selfloop_verify_20260606T065905Z/"
    "reports/c4.k311.seed.spoa2k.default-selfloop.Ygs.graph-report.tsv"
)
DEFAULT_CURRENT_BEST_PNG = Path(
    "/home/erikg/impg/data/c4_spoa2k_default_selfloop_verify_20260606T065905Z/"
    "renders/c4.k311.seed.spoa2k.default-selfloop.Ygs.mean-depth.png"
)
DEFAULT_CURRENT_BEST_URL = (
    "http://hypervolu.me/~erik/impg/"
    "c4-spoa2k-default-selfloop-verify-20260606T065905Z.Ygs.mean-depth.png"
)
DEFAULT_1TO1_BEST = Path(
    "/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/"
    "sorted/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.gfa"
)
DEFAULT_1TO1_BEST_REPORT = Path(
    "/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/"
    "reports/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.graph-report.tsv"
)
DEFAULT_1TO1_BEST_URL = (
    "http://hypervolu.me/~erik/impg/"
    "c4-k311-1to1-noscaffold-spoa2k-default-selfloop-20260606T070148Z.Ygs.mean-depth.png"
)
DEFAULT_PGGB = Path("/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa")
DEFAULT_PGGB_REPORT = Path(
    "/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/"
    "before/pggb_control.graph-report.tsv"
)
DEFAULT_PGGB_URL = "https://hypervolu.me/~erik/impg/c4-pggb-control.png"

DNA_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


@dataclass
class CommandResult:
    label: str
    cmd: list[str]
    stdout: Path
    stderr: Path
    returncode: int


@dataclass
class GfaRun:
    label: str
    kind: str
    input_gfa: Path
    output_gfa: Path
    sorted_gfa: Path
    flank_bp: int
    method: str
    max_iterations: int
    threads: int
    timeout: int


def parse_args() -> argparse.Namespace:
    stamp = dt.datetime.now(dt.UTC).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(f"/home/erikg/impg/data/flank_aware_crush_c4_{stamp}"),
    )
    parser.add_argument("--impg", default="/home/erikg/.cargo/bin/impg")
    parser.add_argument("--gfasort", default="/home/erikg/.cargo/bin/gfasort")
    parser.add_argument("--gfalook", default="/home/erikg/.cargo/bin/gfalook")
    parser.add_argument("--c4-input", type=Path, default=DEFAULT_C4_INPUT)
    parser.add_argument("--c4-flank-bp", type=int, default=64)
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--sanity-threads", type=int, default=4)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--skip-cargo-tests", action="store_true")
    parser.add_argument("--skip-c4", action="store_true")
    parser.add_argument("--skip-render", action="store_true")
    parser.add_argument("--no-upload", action="store_true")
    parser.add_argument("--upload-host", default="erik@hypervolu.me")
    parser.add_argument("--upload-dir", default="www/impg")
    parser.add_argument("--crush-timeout", type=int, default=900)
    parser.add_argument("--c4-crush-timeout", type=int, default=1800)
    parser.add_argument("--report-timeout", type=int, default=600)
    parser.add_argument("--render-timeout", type=int, default=300)
    return parser.parse_args()


def sh_join(cmd: Iterable[str]) -> str:
    return shlex.join([str(part) for part in cmd])


def mkdirs(root: Path) -> None:
    for child in [
        "graphs",
        "sorted",
        "reports",
        "validation",
        "renders",
        "logs",
        "sanity",
        "metadata",
    ]:
        (root / child).mkdir(parents=True, exist_ok=True)


def append_command(out_dir: Path, label: str, cmd: list[str]) -> None:
    commands = out_dir / "commands.sh"
    with commands.open("a") as handle:
        handle.write(f"\n# {label}\n")
        handle.write(sh_join(cmd) + "\n")


def run_logged(
    out_dir: Path,
    label: str,
    cmd: list[str],
    stdout: Path,
    stderr: Path,
    timeout: int,
    cwd: Path | None = None,
    allow_fail: bool = False,
    use_time: bool = True,
    stderr_filter: Callable[[str], str] | None = None,
) -> CommandResult:
    stdout.parent.mkdir(parents=True, exist_ok=True)
    stderr.parent.mkdir(parents=True, exist_ok=True)
    timed_cmd = ["/usr/bin/time", "-v", *cmd] if use_time else cmd
    append_command(out_dir, label, timed_cmd)
    with stdout.open("w") as out, stderr.open("w") as err:
        try:
            if stderr_filter is None:
                completed = subprocess.run(
                    timed_cmd,
                    stdout=out,
                    stderr=err,
                    text=True,
                    timeout=timeout or None,
                    cwd=str(cwd) if cwd else None,
                )
                rc = completed.returncode
            else:
                proc = subprocess.Popen(
                    timed_cmd,
                    stdout=out,
                    stderr=subprocess.PIPE,
                    text=True,
                    cwd=str(cwd) if cwd else None,
                    bufsize=1,
                )
                assert proc.stderr is not None
                deadline = time.monotonic() + timeout if timeout else None
                rc = 0
                while True:
                    line = proc.stderr.readline()
                    if line:
                        err.write(stderr_filter(line))
                    elif proc.poll() is not None:
                        break
                    if deadline is not None and time.monotonic() > deadline:
                        proc.kill()
                        err.write(f"\nTIMEOUT after {timeout}s\n")
                        rc = 124
                        break
                if rc != 124:
                    rc = proc.wait()
        except subprocess.TimeoutExpired:
            err.write(f"\nTIMEOUT after {timeout}s\n")
            rc = 124
    result = CommandResult(label, timed_cmd, stdout, stderr, rc)
    if rc != 0 and not allow_fail:
        raise RuntimeError(f"{label} failed with exit {rc}; see {stderr}")
    return result


def run_capture(args: argparse.Namespace, out_dir: Path, label: str, cmd: list[str]) -> str:
    append_command(out_dir, label, cmd)
    try:
        completed = subprocess.run(cmd, text=True, capture_output=True, timeout=60)
    except Exception as exc:  # noqa: BLE001
        return f"ERROR: {exc}"
    text = (completed.stdout + completed.stderr).strip()
    if completed.returncode != 0 and not text:
        return f"exit {completed.returncode}"
    return text.splitlines()[0] if text else ""


def write_tool_versions(args: argparse.Namespace, out_dir: Path) -> None:
    rows = []
    for label, cmd in [
        ("git_head", ["git", "rev-parse", "HEAD"]),
        ("git_head_short", ["git", "rev-parse", "--short", "HEAD"]),
        ("impg", [args.impg, "--version"]),
        ("gfasort", [args.gfasort, "--version"]),
        ("gfalook", [args.gfalook, "--version"]),
        ("zstd", ["zstd", "--version"]),
    ]:
        rows.append({"tool": label, "value": run_capture(args, out_dir, f"version.{label}", cmd)})
    rows.append({"tool": "utc_started", "value": dt.datetime.now(dt.UTC).isoformat()})
    write_tsv(out_dir / "metadata" / "tool_versions.tsv", rows)


def sanity_gfa_text() -> str:
    # Two independent indel sites with repeated local context. The first site is
    # at the path start, so its left flank is shortened to zero. The second site
    # reaches the path end, so its right flank is shortened to zero.
    return """H\tVN:Z:1.0
S\t1\tAAC
S\t2\tG
S\t3\tGT
S\t4\tTTA
S\t5\tAAC
S\t6\tC
S\t7\tCT
S\t8\tTTA
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t5\t+\t6\t+\t0M
L\t6\t+\t8\t+\t0M
L\t5\t+\t7\t+\t0M
L\t7\t+\t8\t+\t0M
P\tsanity_ref\t1+,2+,4+,5+,6+,8+\t*
P\tsanity_alt\t1+,3+,4+,5+,7+,8+\t*
"""


def write_sanity_fixture(out_dir: Path) -> Path:
    path = out_dir / "sanity" / "flank_required_repeated_path_end.gfa"
    path.write_text(sanity_gfa_text())
    return path


def crush_command(run: GfaRun, impg: str) -> list[str]:
    return [
        impg,
        "crush",
        "-g",
        str(run.input_gfa),
        "-o",
        str(run.output_gfa),
        "--method",
        run.method,
        "--max-iterations",
        str(run.max_iterations),
        "--max-traversal-len",
        "2000" if run.kind == "c4" else "10000",
        "--max-median-traversal-len",
        "2000" if run.kind == "c4" else "10000",
        "--max-total-sequence",
        "1m",
        "--max-traversals",
        "10k",
        "--replacement-flank-bp",
        str(run.flank_bp),
        "-t",
        str(run.threads),
        "-v",
        "1",
    ]


def sort_command(run: GfaRun, gfasort: str) -> list[str]:
    return [
        gfasort,
        "-i",
        str(run.output_gfa),
        "-o",
        str(run.sorted_gfa),
        "-p",
        "Ygs",
        "-t",
        str(run.threads),
        "-v",
        "1",
    ]


def graph_report_commands(args: argparse.Namespace, label: str, gfa: Path) -> tuple[list[str], list[str]]:
    tsv = [
        args.impg,
        "graph-report",
        "-g",
        str(gfa),
        "-o",
        str(args.out_dir / "reports" / f"{label}.graph-report.tsv"),
        "--format",
        "tsv",
        "--povu",
        "--top",
        "20",
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    js = [
        args.impg,
        "graph-report",
        "-g",
        str(gfa),
        "-o",
        str(args.out_dir / "reports" / f"{label}.graph-report.json"),
        "--format",
        "json",
        "--povu",
        "--top",
        "20",
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    return tsv, js


def render_command(args: argparse.Namespace, gfa: Path, png: Path) -> list[str]:
    return [
        args.gfalook,
        "-i",
        str(gfa),
        "-o",
        str(png),
        "-m",
        "-x",
        "3200",
        "-y",
        "1800",
        "-a",
        "3",
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]


def revcomp(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1]


def parse_gfa_sequences(path: Path) -> dict[str, str]:
    opener = gzip.open if path.suffix == ".gz" else open
    segments: dict[str, str] = {}
    paths: dict[str, str] = {}
    with opener(path, "rt") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S" and len(fields) >= 3:
                segments[fields[1]] = fields[2]
            elif fields[0] == "P" and len(fields) >= 3:
                seq = []
                if fields[2] != "*":
                    for raw in fields[2].split(","):
                        if not raw:
                            continue
                        node = raw[:-1]
                        orient = raw[-1]
                        piece = segments[node]
                        seq.append(revcomp(piece) if orient == "-" else piece)
                paths[fields[1]] = "".join(seq)
            elif fields[0] == "W" and len(fields) >= 7:
                seq = []
                walk = fields[6]
                idx = 0
                while idx < len(walk):
                    orient = walk[idx]
                    idx += 1
                    start = idx
                    while idx < len(walk) and walk[idx] not in "><":
                        idx += 1
                    node = walk[start:idx]
                    piece = segments[node]
                    seq.append(revcomp(piece) if orient == "<" else piece)
                paths[fields[3]] = "".join(seq)
    return paths


def compare_paths(expected: Path, observed: Path, label: str, out_dir: Path) -> dict[str, str]:
    expected_paths = parse_gfa_sequences(expected)
    observed_paths = parse_gfa_sequences(observed)
    missing = sorted(set(expected_paths) - set(observed_paths))
    extra = sorted(set(observed_paths) - set(expected_paths))
    mismatched = sorted(
        name
        for name in set(expected_paths) & set(observed_paths)
        if expected_paths[name] != observed_paths[name]
    )
    detail_path = out_dir / "validation" / f"{label}.path-detail.tsv"
    detail_rows = []
    for name in missing:
        detail_rows.append({"path": name, "status": "missing", "expected_len": len(expected_paths[name]), "observed_len": ""})
    for name in extra:
        detail_rows.append({"path": name, "status": "extra", "expected_len": "", "observed_len": len(observed_paths[name])})
    for name in mismatched[:100]:
        detail_rows.append(
            {
                "path": name,
                "status": "sequence_mismatch",
                "expected_len": len(expected_paths[name]),
                "observed_len": len(observed_paths[name]),
            }
        )
    if detail_rows:
        write_tsv(detail_path, detail_rows)
    else:
        detail_path.write_text("path\tstatus\texpected_len\tobserved_len\n")
    return {
        "label": label,
        "expected_gfa": str(expected),
        "observed_gfa": str(observed),
        "expected_paths": str(len(expected_paths)),
        "observed_paths": str(len(observed_paths)),
        "missing_paths": str(len(missing)),
        "extra_paths": str(len(extra)),
        "path_name_mismatches": str(len(missing) + len(extra)),
        "sequence_spelling_mismatches": str(len(mismatched)),
        "full_path_preservation": str(not missing and not extra and not mismatched).lower(),
        "detail_tsv": str(detail_path),
        "expected_path_names": ",".join(sorted(expected_paths)),
        "observed_path_names": ",".join(sorted(observed_paths)),
    }


def read_tsv_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return dict(rows[0]) if rows else {}


def write_tsv(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    if not keys:
        keys = ["status"]
        rows = [{"status": "empty"}]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=keys)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_time_log(path: Path) -> dict[str, str]:
    values: dict[str, str] = {}
    if not path.exists():
        return values
    with path.open(errors="replace") as handle:
        for line in handle:
            if "Elapsed (wall clock) time" in line:
                values["wall"] = line.rsplit("):", 1)[-1].strip()
            elif "Maximum resident set size" in line:
                values["max_rss_kb"] = line.rsplit(":", 1)[-1].strip()
            elif "Exit status" in line:
                values["exit_status"] = line.rsplit(":", 1)[-1].strip()
    return values


ROUND_RE = re.compile(
    r"crush round (\d+): (\d+) POVU site\(s\).*?, (\d+) candidate\(s\).*?, (\d+) selected"
)
NO_ELIGIBLE_RE = re.compile(r"crush round (\d+): no eligible candidates from (\d+) POVU site\(s\)")
FINAL_RE = re.compile(r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds")
SKIPPED_RE = re.compile(r"skipped (\d+) explicit .* candidate\(s\) outside direct length budgets")
PROGRESS_RE = re.compile(r"replacement build progress \d+/\d+ \(([^)]+)\)")
FLANK_OCC_RE = re.compile(
    r"orientation=(?P<orientation>forward|reverse).*?"
    r"left\(path_side=[^,]+,path_req=(?P<left_req>\d+),path_actual=(?P<left_actual>\d+),"
    r"canonical_side=[^,]+,canonical_actual=(?P<left_canon>\d+),.*?trunc=(?P<left_trunc>[^)]+)\) "
    r"right\(path_side=[^,]+,path_req=(?P<right_req>\d+),path_actual=(?P<right_actual>\d+),"
    r"canonical_side=[^,]+,canonical_actual=(?P<right_canon>\d+),.*?trunc=(?P<right_trunc>[^)]+)\) "
    r"trim\(left=(?P<trim_left>\d+),right=(?P<trim_right>\d+),"
)
TRIM_RE = re.compile(
    r"orientation=(?P<orientation>forward|reverse).*?"
    r"left_trim=(?P<trim_left>\d+) right_trim=(?P<trim_right>\d+).*?"
    r"requested_flanks=\((?P<left_req>\d+),(?P<right_req>\d+)\) "
    r"actual_path_flanks=\((?P<left_actual>\d+),(?P<right_actual>\d+)\) "
    r"canonical_flanks=\((?P<left_canon>\d+),(?P<right_canon>\d+)\)"
)


def counter_dist(values: Iterable[int]) -> str:
    counts = Counter(values)
    return ";".join(f"{key}:{counts[key]}" for key in sorted(counts))


def counter_items(counter: Counter[str]) -> str:
    return ";".join(f"{key}:{counter[key]}" for key in sorted(counter))


def parse_candidate_accounting(label: str, stderr_path: Path, flank_bp: int, method: str) -> dict[str, str]:
    selected_total = 0
    attempted_round_candidates = 0
    povu_sites = 0
    no_eligible_sites = 0
    budget_skips = 0
    progress_statuses: Counter[str] = Counter()
    left_actual: list[int] = []
    right_actual: list[int] = []
    left_canon: list[int] = []
    right_canon: list[int] = []
    trim_left: list[int] = []
    trim_right: list[int] = []
    left_trunc: Counter[str] = Counter()
    right_trunc: Counter[str] = Counter()
    orientations: Counter[str] = Counter()
    trim_boundary_errors = 0
    path_invalid_errors = 0
    resolved = 0
    bailed = 0
    candidates_seen = 0
    rounds = 0
    if not stderr_path.exists():
        line_iter: Iterable[str] = []
    else:
        line_iter = stderr_path.open(errors="replace")
    for line in line_iter:
        if match := ROUND_RE.search(line):
            povu_sites += int(match.group(2))
            attempted_round_candidates += int(match.group(3))
            selected_total += int(match.group(4))
        if match := NO_ELIGIBLE_RE.search(line):
            no_eligible_sites += int(match.group(2))
        if match := SKIPPED_RE.search(line):
            budget_skips += int(match.group(1))
        if match := PROGRESS_RE.search(line):
            progress_statuses[match.group(1)] += 1
        if "interior span" in line or "interior bounds inverted" in line:
            trim_boundary_errors += 1
        if "path-invalid" in line or "failed exact path" in line or "path-corruption" in line:
            path_invalid_errors += 1
        if match := FINAL_RE.search(line):
            resolved = int(match.group(1))
            bailed = int(match.group(2))
            candidates_seen = int(match.group(3))
            rounds = int(match.group(4))
        if "crush flank-aware trim:" in line:
            if trim := TRIM_RE.search(line):
                orientations[trim.group("orientation")] += 1
                left_req = int(trim.group("left_req"))
                right_req = int(trim.group("right_req"))
                left_val = int(trim.group("left_actual"))
                right_val = int(trim.group("right_actual"))
                left_actual.append(left_val)
                right_actual.append(right_val)
                left_canon.append(int(trim.group("left_canon")))
                right_canon.append(int(trim.group("right_canon")))
                trim_left.append(int(trim.group("trim_left")))
                trim_right.append(int(trim.group("trim_right")))
                left_trunc["shortened" if left_val < left_req else "full"] += 1
                right_trunc["shortened" if right_val < right_req else "full"] += 1
        if "crush flank-aware candidate:" in line or "crush flank-aware build:" in line:
            for occ in FLANK_OCC_RE.finditer(line):
                orientations[occ.group("orientation")] += 1
                left_actual.append(int(occ.group("left_actual")))
                right_actual.append(int(occ.group("right_actual")))
                left_canon.append(int(occ.group("left_canon")))
                right_canon.append(int(occ.group("right_canon")))
                trim_left.append(int(occ.group("trim_left")))
                trim_right.append(int(occ.group("trim_right")))
                left_trunc[occ.group("left_trunc")] += 1
                right_trunc[occ.group("right_trunc")] += 1
    if hasattr(line_iter, "close"):
        line_iter.close()  # type: ignore[attr-defined]
    if not candidates_seen:
        candidates_seen = attempted_round_candidates
    hard_rejections = bailed + path_invalid_errors + trim_boundary_errors
    return {
        "label": label,
        "method": method,
        "requested_flank_bp": str(flank_bp),
        "povu_sites_seen_round_sum": str(povu_sites),
        "no_eligible_final_round_sites": str(no_eligible_sites),
        "attempted_candidates_seen": str(candidates_seen),
        "round_candidate_sum": str(attempted_round_candidates),
        "selected_candidates": str(selected_total),
        "accepted_candidates": str(resolved),
        "bailed_candidates": str(bailed),
        "hard_path_corruption_rejections": str(hard_rejections),
        "trim_boundary_ambiguity_or_span_errors": str(trim_boundary_errors),
        "path_invalid_errors": str(path_invalid_errors),
        "budget_ineligible_candidates": str(budget_skips),
        "quality_guard_rejections": "0",
        "metric_based_filter_used": "false",
        "replacement_progress_status_samples": counter_items(progress_statuses),
        "rounds_with_replacements": str(rounds),
        "flank_occurrences_logged": str(len(left_actual)),
        "orientation_distribution": counter_items(orientations),
        "left_actual_path_bp_distribution": counter_dist(left_actual),
        "right_actual_path_bp_distribution": counter_dist(right_actual),
        "left_actual_canonical_bp_distribution": counter_dist(left_canon),
        "right_actual_canonical_bp_distribution": counter_dist(right_canon),
        "left_trim_bp_distribution": counter_dist(trim_left),
        "right_trim_bp_distribution": counter_dist(trim_right),
        "left_truncation_distribution": counter_items(left_trunc),
        "right_truncation_distribution": counter_items(right_trunc),
    }


def compact_crush_stderr_line(line: str) -> str:
    if "crush flank-aware trim:" in line:
        line = re.sub(r"candidate_id=.*? method=", "candidate_id=<omitted> method=", line)
        line = re.sub(r" path=.*? orientation=", " path=<omitted> orientation=", line)
    elif ("crush flank-aware candidate:" in line or "crush flank-aware build:" in line) and len(line) > 4000:
        line = line[:4000].rstrip("\n") + " ...<truncated flank diagnostic>\n"
    return line


def summarize_povu_json(path: Path) -> dict[str, str]:
    if not path.exists():
        return {
            "povu_max_reference_span_steps": "",
            "povu_top_reference_span_steps": "",
        }
    try:
        data = json.loads(path.read_text())
    except json.JSONDecodeError:
        return {
            "povu_max_reference_span_steps": "",
            "povu_top_reference_span_steps": "",
        }
    sites = (((data.get("povu") or {}).get("top_sites")) or [])
    spans = [int(site.get("reference_span_steps") or 0) for site in sites]
    return {
        "povu_max_reference_span_steps": str(max(spans) if spans else 0),
        "povu_top_reference_span_steps": ",".join(str(value) for value in sorted(spans, reverse=True)[:10]),
    }


def metric_row(
    label: str,
    kind: str,
    gfa: Path,
    sorted_gfa: Path,
    report_tsv: Path,
    report_json: Path,
    path_row: dict[str, str] | None,
    png: Path | None = None,
    png_url: str = "",
) -> dict[str, str]:
    report = read_tsv_row(report_tsv)
    fields = [
        "status",
        "failures",
        "warnings",
        "segments",
        "links",
        "paths",
        "path_steps",
        "total_segment_bp",
        "node_coverage_mean",
        "node_coverage_bp_weighted_mean",
        "node_coverage_p10",
        "node_coverage_median",
        "node_coverage_p90",
        "singleton_nodes",
        "singleton_bp",
        "path_white_space_bp_p99",
        "path_white_space_bp_max",
        "path_white_space_bridges_ge_threshold",
        "segment_occupancy_bp_fraction",
        "segment_white_space_bp_fraction",
        "segment_white_space_bp_total",
        "sparse_coverage_segment_bp",
        "duplicate_sequence_frac",
        "local_repeat_context_nodes",
        "local_repeat_context_occurrences",
        "direct_self_loop_edges",
        "direct_self_loop_nodes",
        "adjacent_same_node_path_steps",
        "adjacent_same_step_path_steps",
        "self_loop_repeat_runs",
        "self_loop_max_repeat_run_len",
        "self_loop_length_buckets",
        "povu_sites",
        "povu_leaf_sites",
    ]
    row = {
        "label": label,
        "kind": kind,
        "gfa": str(gfa),
        "sorted_gfa": str(sorted_gfa),
        "graph_report_tsv": str(report_tsv),
        "graph_report_json": str(report_json),
        "png": str(png or ""),
        "png_url": png_url,
    }
    for field in fields:
        row[field] = str(report.get(field, ""))
    row.update(summarize_povu_json(report_json))
    if path_row:
        for key in [
            "expected_paths",
            "observed_paths",
            "missing_paths",
            "extra_paths",
            "path_name_mismatches",
            "sequence_spelling_mismatches",
            "full_path_preservation",
            "detail_tsv",
        ]:
            row[key] = path_row.get(key, "")
    else:
        row.update(
            {
                "expected_paths": "",
                "observed_paths": "",
                "missing_paths": "",
                "extra_paths": "",
                "path_name_mismatches": "",
                "sequence_spelling_mismatches": "",
                "full_path_preservation": "",
                "detail_tsv": "",
            }
        )
    return row


def compress_gfa(out_dir: Path, label: str, gfa: Path) -> str:
    if not gfa.exists():
        return ""
    zst = Path(str(gfa) + ".zst")
    cmd = ["zstd", "-f", "-T0", str(gfa), "-o", str(zst)]
    run_logged(
        out_dir,
        f"{label}.zstd",
        cmd,
        out_dir / "logs" / f"{label}.zstd.stdout.log",
        out_dir / "logs" / f"{label}.zstd.stderr.log",
        300,
        allow_fail=True,
        use_time=True,
    )
    return str(zst) if zst.exists() else ""


def upload_png(args: argparse.Namespace, label: str, png: Path) -> str:
    if args.no_upload or not png.exists():
        return ""
    remote_name = f"flank-aware-crush-c4-{label}-{args.out_dir.name}.png"
    remote = f"{args.upload_host}:{args.upload_dir}/{remote_name}"
    result = run_logged(
        args.out_dir,
        f"{label}.upload",
        ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", str(png), remote],
        args.out_dir / "logs" / f"{label}.upload.stdout.log",
        args.out_dir / "logs" / f"{label}.upload.stderr.log",
        120,
        allow_fail=True,
        use_time=True,
    )
    if result.returncode == 0:
        return f"http://hypervolu.me/~erik/impg/{remote_name}"
    return ""


def run_cargo_sanity(args: argparse.Namespace) -> dict[str, str]:
    if args.skip_cargo_tests:
        return {"label": "cargo.flank_sanity", "status": "skipped", "log": ""}
    cmd = ["zsh", "-lc", "source ./env.sh && cargo test --lib flank -- --nocapture"]
    result = run_logged(
        args.out_dir,
        "cargo.flank_sanity",
        cmd,
        args.out_dir / "logs" / "cargo.flank_sanity.stdout.log",
        args.out_dir / "logs" / "cargo.flank_sanity.stderr.log",
        900,
        cwd=Path.cwd(),
        allow_fail=True,
        use_time=True,
    )
    stdout = result.stdout.read_text(errors="replace") if result.stdout.exists() else ""
    stderr = result.stderr.read_text(errors="replace") if result.stderr.exists() else ""
    text = stdout + "\n" + stderr
    return {
        "label": "cargo.flank_sanity",
        "status": "pass" if result.returncode == 0 else "fail",
        "exit_status": str(result.returncode),
        "log": str(result.stderr),
        "tests_covering": (
            "flank-aware repeated occurrence; path-boundary short flanks; "
            "reverse occurrence trim/lacing; global resolver flank trim; "
            "hard trim-boundary rejection; diagnostics formatting"
        ),
        "test_result_summary": next((line for line in text.splitlines() if "test result:" in line), ""),
    }


def ensure_report_for_baseline(
    args: argparse.Namespace,
    label: str,
    gfa: Path,
    known_report: Path | None,
) -> tuple[Path, Path]:
    tsv_out = args.out_dir / "reports" / f"{label}.graph-report.tsv"
    json_out = args.out_dir / "reports" / f"{label}.graph-report.json"
    if known_report and known_report.exists() and not args.force:
        tsv_out.write_text(known_report.read_text())
    else:
        tsv_cmd, _ = graph_report_commands(args, label, gfa)
        run_logged(
            args.out_dir,
            f"{label}.graph-report.tsv",
            tsv_cmd,
            args.out_dir / "logs" / f"{label}.graph-report.tsv.stdout.log",
            args.out_dir / "logs" / f"{label}.graph-report.tsv.stderr.log",
            args.report_timeout,
            allow_fail=True,
        )
    if not json_out.exists() or args.force:
        _, json_cmd = graph_report_commands(args, label, gfa)
        run_logged(
            args.out_dir,
            f"{label}.graph-report.json",
            json_cmd,
            args.out_dir / "logs" / f"{label}.graph-report.json.stdout.log",
            args.out_dir / "logs" / f"{label}.graph-report.json.stderr.log",
            args.report_timeout,
            allow_fail=True,
        )
    return tsv_out, json_out


def run_gfa_variant(args: argparse.Namespace, run: GfaRun) -> tuple[dict[str, str], dict[str, str], str]:
    if args.force or not run.output_gfa.exists():
        run_logged(
            args.out_dir,
            f"{run.label}.crush",
            crush_command(run, args.impg),
            args.out_dir / "logs" / f"{run.label}.crush.stdout.log",
            args.out_dir / "logs" / f"{run.label}.crush.stderr.log",
            run.timeout,
            stderr_filter=compact_crush_stderr_line if run.flank_bp else None,
        )
    if args.force or not run.sorted_gfa.exists():
        run_logged(
            args.out_dir,
            f"{run.label}.gfasort",
            sort_command(run, args.gfasort),
            args.out_dir / "logs" / f"{run.label}.gfasort.stdout.log",
            args.out_dir / "logs" / f"{run.label}.gfasort.stderr.log",
            args.report_timeout,
        )
    tsv_cmd, json_cmd = graph_report_commands(args, run.label, run.sorted_gfa)
    if args.force or not (args.out_dir / "reports" / f"{run.label}.graph-report.tsv").exists():
        run_logged(
            args.out_dir,
            f"{run.label}.graph-report.tsv",
            tsv_cmd,
            args.out_dir / "logs" / f"{run.label}.graph-report.tsv.stdout.log",
            args.out_dir / "logs" / f"{run.label}.graph-report.tsv.stderr.log",
            args.report_timeout,
        )
    if args.force or not (args.out_dir / "reports" / f"{run.label}.graph-report.json").exists():
        run_logged(
            args.out_dir,
            f"{run.label}.graph-report.json",
            json_cmd,
            args.out_dir / "logs" / f"{run.label}.graph-report.json.stdout.log",
            args.out_dir / "logs" / f"{run.label}.graph-report.json.stderr.log",
            args.report_timeout,
        )
    path_row = compare_paths(run.input_gfa, run.sorted_gfa, run.label, args.out_dir)
    accounting = parse_candidate_accounting(
        run.label,
        args.out_dir / "logs" / f"{run.label}.crush.stderr.log",
        run.flank_bp,
        run.method,
    )
    zst = compress_gfa(args.out_dir, run.label, run.output_gfa)
    return path_row, accounting, zst


def main() -> int:
    args = parse_args()
    mkdirs(args.out_dir)
    (args.out_dir / "commands.sh").write_text("#!/bin/sh\nset -eu\n")
    os.chmod(args.out_dir / "commands.sh", 0o755)
    write_tool_versions(args, args.out_dir)

    sanity_test_row = run_cargo_sanity(args)
    sanity_gfa = write_sanity_fixture(args.out_dir)
    sanity_runs = [
        GfaRun(
            label="sanity.flank0",
            kind="sanity",
            input_gfa=sanity_gfa,
            output_gfa=args.out_dir / "graphs" / "sanity.flank0.gfa",
            sorted_gfa=args.out_dir / "sorted" / "sanity.flank0.Ygs.gfa",
            flank_bp=0,
            method="poa",
            max_iterations=1,
            threads=args.sanity_threads,
            timeout=args.crush_timeout,
        ),
        GfaRun(
            label="sanity.flank4",
            kind="sanity",
            input_gfa=sanity_gfa,
            output_gfa=args.out_dir / "graphs" / "sanity.flank4.gfa",
            sorted_gfa=args.out_dir / "sorted" / "sanity.flank4.Ygs.gfa",
            flank_bp=4,
            method="poa",
            max_iterations=1,
            threads=args.sanity_threads,
            timeout=args.crush_timeout,
        ),
    ]
    if not args.c4_input.exists() and not args.skip_c4:
        raise FileNotFoundError(f"C4 input does not exist: {args.c4_input}")
    c4_runs = [] if args.skip_c4 else [
        GfaRun(
            label="c4.flank0",
            kind="c4",
            input_gfa=args.c4_input,
            output_gfa=args.out_dir / "graphs" / "c4.flank0.gfa",
            sorted_gfa=args.out_dir / "sorted" / "c4.flank0.Ygs.gfa",
            flank_bp=0,
            method="poa",
            max_iterations=5,
            threads=args.threads,
            timeout=args.c4_crush_timeout,
        ),
        GfaRun(
            label=f"c4.flank{args.c4_flank_bp}",
            kind="c4",
            input_gfa=args.c4_input,
            output_gfa=args.out_dir / "graphs" / f"c4.flank{args.c4_flank_bp}.gfa",
            sorted_gfa=args.out_dir / "sorted" / f"c4.flank{args.c4_flank_bp}.Ygs.gfa",
            flank_bp=args.c4_flank_bp,
            method="poa",
            max_iterations=5,
            threads=args.threads,
            timeout=args.c4_crush_timeout,
        ),
    ]

    path_rows: list[dict[str, str]] = []
    accounting_rows: list[dict[str, str]] = []
    artifact_rows: list[dict[str, str]] = []
    metric_rows: list[dict[str, str]] = []
    url_rows: list[dict[str, str]] = []

    for run in [*sanity_runs, *c4_runs]:
        path_row, accounting, zst = run_gfa_variant(args, run)
        path_rows.append(path_row)
        accounting_rows.append(accounting)
        artifact_rows.append(
            {
                "label": run.label,
                "input_gfa": str(run.input_gfa),
                "output_gfa": str(run.output_gfa),
                "sorted_gfa": str(run.sorted_gfa),
                "zst": zst,
                "crush_stderr": str(args.out_dir / "logs" / f"{run.label}.crush.stderr.log"),
                "graph_report_tsv": str(args.out_dir / "reports" / f"{run.label}.graph-report.tsv"),
                "graph_report_json": str(args.out_dir / "reports" / f"{run.label}.graph-report.json"),
            }
        )
        png = None
        url = ""
        if run.kind == "c4" and not args.skip_render:
            png = args.out_dir / "renders" / f"{run.label}.Ygs.gfalook-m.png"
            if args.force or not png.exists():
                run_logged(
                    args.out_dir,
                    f"{run.label}.gfalook",
                    render_command(args, run.sorted_gfa, png),
                    args.out_dir / "logs" / f"{run.label}.gfalook.stdout.log",
                    args.out_dir / "logs" / f"{run.label}.gfalook.stderr.log",
                    args.render_timeout,
                    allow_fail=True,
                )
            url = upload_png(args, run.label, png)
            url_rows.append({"label": run.label, "png": str(png), "url": url})
        metric_rows.append(
            metric_row(
                run.label,
                run.kind,
                run.output_gfa,
                run.sorted_gfa,
                args.out_dir / "reports" / f"{run.label}.graph-report.tsv",
                args.out_dir / "reports" / f"{run.label}.graph-report.json",
                path_row,
                png,
                url,
            )
        )

    baseline_specs = [
        ("current_best_spoa2k_default_selfloop", DEFAULT_CURRENT_BEST, DEFAULT_CURRENT_BEST_REPORT, DEFAULT_CURRENT_BEST_URL),
        ("current_best_1to1_noscaffold", DEFAULT_1TO1_BEST, DEFAULT_1TO1_BEST_REPORT, DEFAULT_1TO1_BEST_URL),
        ("pggb_control", DEFAULT_PGGB, DEFAULT_PGGB_REPORT, DEFAULT_PGGB_URL),
    ]
    for label, gfa, known_report, url in baseline_specs:
        if not gfa.exists():
            artifact_rows.append({"label": label, "missing_artifact": str(gfa)})
            continue
        report_tsv, report_json = ensure_report_for_baseline(args, label, gfa, known_report)
        metric_rows.append(
            metric_row(
                label,
                "baseline",
                gfa,
                gfa,
                report_tsv,
                report_json,
                None,
                None,
                url,
            )
        )
        url_rows.append({"label": label, "png": "", "url": url})

    write_tsv(args.out_dir / "validation" / "path_preservation.tsv", path_rows)
    write_tsv(args.out_dir / "reports" / "candidate_accounting.tsv", accounting_rows)
    write_tsv(args.out_dir / "reports" / "metrics.tsv", metric_rows)
    write_tsv(args.out_dir / "artifacts.tsv", artifact_rows)
    write_tsv(args.out_dir / "renders" / "uploaded_urls.tsv", url_rows)
    write_tsv(args.out_dir / "reports" / "sanity_tests.tsv", [sanity_test_row])

    status = {
        "out_dir": str(args.out_dir),
        "path_preservation_tsv": str(args.out_dir / "validation" / "path_preservation.tsv"),
        "candidate_accounting_tsv": str(args.out_dir / "reports" / "candidate_accounting.tsv"),
        "metrics_tsv": str(args.out_dir / "reports" / "metrics.tsv"),
        "uploaded_urls_tsv": str(args.out_dir / "renders" / "uploaded_urls.tsv"),
        "sanity_tests_tsv": str(args.out_dir / "reports" / "sanity_tests.tsv"),
    }
    (args.out_dir / "run_summary.json").write_text(json.dumps(status, indent=2, sort_keys=True) + "\n")
    print(json.dumps(status, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
