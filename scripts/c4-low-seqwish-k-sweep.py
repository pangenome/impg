#!/usr/bin/env python3
"""Run the full-C4 seqwish min-match sweep with POA polish.

This is an experiment driver, not a library entry point. It writes all heavy
artifacts under the requested data directory and records commands, time/RSS,
graph reports, path-spelling validation, renders, and upload URLs.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


DEFAULT_FASTA = Path(
    "/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/"
    "c4_whole_region.fa"
)
DEFAULT_COMPARE = Path("/home/erikg/impg/target/release/examples/compare_gfa_paths")
PUBLIC_BASE = "https://hypervolu.me/~erik/impg/"
UPLOAD_TARGET = "erik@hypervolu.me:www/impg/"


@dataclass
class Runtime:
    label: str
    wall_time: str = ""
    max_rss_kb: str = ""
    exit_status: str = ""
    stdout: str = ""
    stderr: str = ""
    time_log: str = ""
    command: str = ""


@dataclass
class RunState:
    out_dir: Path
    commands_path: Path
    runtimes: dict[str, Runtime] = field(default_factory=dict)
    rows: list[dict[str, str]] = field(default_factory=list)
    uploads: dict[str, str] = field(default_factory=dict)
    validation: dict[str, dict[str, str]] = field(default_factory=dict)
    raw_paf: Path | None = None


def parse_args() -> argparse.Namespace:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=DEFAULT_FASTA)
    parser.add_argument(
        "--out-dir",
        type=Path,
        default=Path(f"/home/erikg/impg/data/c4_low_seqwish_k_{stamp}"),
    )
    parser.add_argument("--ks", default="79,127,191,255,311")
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--compare-gfa-paths", type=Path, default=DEFAULT_COMPARE)
    parser.add_argument("--upload-target", default=UPLOAD_TARGET)
    parser.add_argument("--public-base", default=PUBLIC_BASE)
    parser.add_argument("--no-upload", action="store_true")
    parser.add_argument("--skip-existing", action="store_true")
    parser.add_argument(
        "--skip-2kb",
        action="store_true",
        help="Skip the optional 2kb POA pass for the two best seed graphs.",
    )
    parser.add_argument(
        "--reuse-paf",
        action="store_true",
        default=True,
        help="Reuse the first successful raw PAF for later K values.",
    )
    return parser.parse_args()


def ensure_dirs(out_dir: Path) -> None:
    for name in [
        "debug",
        "graphs",
        "logs",
        "renders",
        "reports",
        "sorted",
        "validation",
    ]:
        (out_dir / name).mkdir(parents=True, exist_ok=True)


def shell_join(command: Iterable[object]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def append_command(state: RunState, label: str, command: list[object]) -> None:
    with state.commands_path.open("a", encoding="utf-8") as handle:
        handle.write(f"\n## {label}\n\n```bash\n{shell_join(command)}\n```\n")


def parse_time_log(path: Path) -> tuple[str, str, str]:
    wall = ""
    rss = ""
    status = ""
    if not path.exists():
        return wall, rss, status
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            wall = line.rsplit(": ", 1)[-1].strip()
        elif "Maximum resident set size (kbytes)" in line:
            rss = line.split(":", 1)[1].strip()
        elif line.startswith("\tExit status:"):
            status = line.split(":", 1)[1].strip()
    return wall, rss, status


def run_timed(
    state: RunState,
    label: str,
    command: list[object],
    *,
    check: bool = True,
    cwd: Path | None = None,
) -> Runtime:
    stdout_path = state.out_dir / "logs" / f"{label}.stdout.log"
    stderr_path = state.out_dir / "logs" / f"{label}.stderr.log"
    time_path = state.out_dir / "logs" / f"{label}.time.txt"
    full_command = ["/usr/bin/time", "-v", "-o", str(time_path), *map(str, command)]
    append_command(state, label, full_command)
    print(f"[{datetime.now(timezone.utc).isoformat()}] running {label}", flush=True)
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        proc = subprocess.run(full_command, cwd=cwd, stdout=stdout, stderr=stderr)
    wall, rss, status = parse_time_log(time_path)
    if not status:
        status = str(proc.returncode)
    runtime = Runtime(
        label=label,
        wall_time=wall,
        max_rss_kb=rss,
        exit_status=status,
        stdout=str(stdout_path),
        stderr=str(stderr_path),
        time_log=str(time_path),
        command=shell_join(full_command),
    )
    state.runtimes[label] = runtime
    write_runtime_summary(state)
    if check and proc.returncode != 0:
        raise subprocess.CalledProcessError(proc.returncode, full_command)
    return runtime


def write_runtime_summary(state: RunState) -> None:
    path = state.out_dir / "runtime_summary.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "label",
                "wall_time",
                "max_rss_kb",
                "exit_status",
                "stdout",
                "stderr",
                "time_log",
                "command",
            ],
            delimiter="\t",
        )
        writer.writeheader()
        for runtime in state.runtimes.values():
            writer.writerow(runtime.__dict__)


def write_uploaded_urls(state: RunState) -> None:
    path = state.out_dir / "uploaded_urls.tsv"
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["label", "png", "url"])
        for label, url in sorted(state.uploads.items()):
            writer.writerow([label, state.out_dir / "renders" / f"{label}.png", url])


def write_validation_summary(state: RunState) -> None:
    path = state.out_dir / "path_validation.tsv"
    fields = [
        "label",
        "expected_paths",
        "observed_paths",
        "missing_paths",
        "extra_paths",
        "spelling_mismatches",
        "exit_status",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for label, data in sorted(state.validation.items()):
            row = {field: data.get(field, "") for field in fields}
            row["label"] = label
            writer.writerow(row)


def read_report_row(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return {key: (value or "") for key, value in row.items()}
    raise ValueError(f"empty graph-report TSV: {path}")


def int_value(row: dict[str, str], key: str) -> int:
    value = row.get(key, "")
    if value == "":
        return 0
    try:
        return int(float(value))
    except ValueError:
        return 0


def float_value(row: dict[str, str], key: str) -> float:
    value = row.get(key, "")
    if value == "":
        return 0.0
    try:
        return float(value)
    except ValueError:
        return 0.0


def run_graph_report(state: RunState, label: str, gfa: Path, threads: int) -> Path:
    tsv = state.out_dir / "reports" / f"{label}.graph-report.tsv"
    md = state.out_dir / "reports" / f"{label}.graph-report.md"
    if not tsv.exists():
        run_timed(
            state,
            f"{label}.graph-report.tsv",
            [
                "impg",
                "graph-report",
                "-g",
                gfa,
                "-o",
                tsv,
                "--format",
                "tsv",
                "--povu",
                "--top",
                "20",
                "-t",
                threads,
                "-v",
                "1",
            ],
        )
    if not md.exists():
        run_timed(
            state,
            f"{label}.graph-report.md",
            [
                "impg",
                "graph-report",
                "-g",
                gfa,
                "-o",
                md,
                "--format",
                "markdown",
                "--povu",
                "--top",
                "20",
                "-t",
                threads,
                "-v",
                "1",
            ],
        )
    return tsv


def render_and_upload(
    state: RunState,
    label: str,
    gfa: Path,
    threads: int,
    upload_target: str,
    public_base: str,
    no_upload: bool,
) -> tuple[Path, str]:
    sorted_gfa = state.out_dir / "sorted" / f"{label}.Ygs.gfa"
    png = state.out_dir / "renders" / f"{label}.png"
    if not sorted_gfa.exists():
        run_timed(
            state,
            f"{label}.gfasort",
            ["gfasort", "-i", gfa, "-o", sorted_gfa, "-p", "Ygs", "-t", threads, "-v", "1"],
        )
    if not png.exists():
        run_timed(
            state,
            f"{label}.gfalook",
            [
                "gfalook",
                "-i",
                sorted_gfa,
                "-o",
                png,
                "-m",
                "-x",
                "3200",
                "-y",
                "1800",
                "-a",
                "3",
                "-t",
                threads,
                "-v",
                "1",
            ],
        )
    url = ""
    if no_upload:
        url = f"local:{png}"
    else:
        dest = f"{upload_target.rstrip('/')}/{png.name}"
        run_timed(state, f"{label}.upload", ["scp", png, dest])
        url = f"{public_base.rstrip('/')}/{png.name}"
    state.uploads[label] = url
    write_uploaded_urls(state)
    return png, url


def parse_crush_summary(stderr_path: Path) -> dict[str, str]:
    result = {"resolved": "", "bailed": "", "candidates_seen": "", "rounds": ""}
    if not stderr_path.exists():
        return result
    pattern = re.compile(r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds")
    for line in stderr_path.read_text(encoding="utf-8", errors="replace").splitlines():
        match = pattern.search(line)
        if match:
            result["resolved"] = match.group(1)
            result["bailed"] = match.group(2)
            result["candidates_seen"] = match.group(3)
            result["rounds"] = match.group(4)
    return result


def parse_compare_output(stdout_path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not stdout_path.exists():
        return result
    for line in stdout_path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.rstrip("\n").split("\t")
        if len(parts) == 2:
            result[parts[0]] = parts[1]
    return result


def compare_paths(
    state: RunState,
    label: str,
    compare_binary: Path,
    seed_gfa: Path,
    polished_gfa: Path,
) -> dict[str, str]:
    runtime = run_timed(
        state,
        f"{label}.compare_gfa_paths",
        [compare_binary, seed_gfa, polished_gfa],
    )
    data = parse_compare_output(Path(runtime.stdout))
    data["exit_status"] = runtime.exit_status
    state.validation[label] = data
    write_validation_summary(state)
    if data.get("spelling_mismatches") != "0":
        raise RuntimeError(f"{label}: spelling_mismatches={data.get('spelling_mismatches')}")
    return data


def metric_row(
    *,
    k: int,
    x: str,
    label: str,
    gfa: Path,
    report: dict[str, str],
    png_url: str,
    graph_runtime: Runtime | None = None,
    crush_runtime: Runtime | None = None,
    sort_runtime: Runtime | None = None,
    render_runtime: Runtime | None = None,
    upload_runtime: Runtime | None = None,
    compare: dict[str, str] | None = None,
    crush_summary: dict[str, str] | None = None,
    paf_source: str = "",
) -> dict[str, str]:
    fields = {
        "k": str(k),
        "x": x,
        "label": label,
        "gfa": str(gfa),
        "segments": report.get("segments", ""),
        "links": report.get("links", ""),
        "paths": report.get("paths", ""),
        "path_steps": report.get("path_steps", ""),
        "segment_bp": report.get("total_segment_bp", ""),
        "direct_self_loop_edges": report.get("direct_self_loop_edges", "0"),
        "direct_self_loop_nodes": report.get("direct_self_loop_nodes", "0"),
        "adjacent_same_node_repeats": report.get("adjacent_same_node_path_steps", "0"),
        "adjacent_same_step_repeats": report.get("adjacent_same_step_path_steps", "0"),
        "self_loop_repeat_runs": report.get("self_loop_repeat_runs", "0"),
        "self_loop_max_repeat_run_len": report.get("self_loop_max_repeat_run_len", "0"),
        "self_loop_length_buckets": report.get("self_loop_length_buckets", ""),
        "singleton_bp": report.get("singleton_bp", ""),
        "bp_weighted_node_coverage": report.get("node_coverage_bp_weighted_mean", ""),
        "path_white_space_bp_p99": report.get("path_white_space_bp_p99", ""),
        "path_white_space_bp_max": report.get("path_white_space_bp_max", ""),
        "duplicate_sequence_frac": report.get("duplicate_sequence_frac", ""),
        "local_repeat_context_nodes": report.get("local_repeat_context_nodes", ""),
        "local_repeat_context_occurrences": report.get("local_repeat_context_occurrences", ""),
        "resolved": "",
        "bailed": "",
        "candidates_seen": "",
        "rounds": "",
        "spelling_mismatches": "",
        "graph_wall": "",
        "graph_max_rss_kb": "",
        "crush_wall": "",
        "crush_max_rss_kb": "",
        "sort_wall": "",
        "sort_max_rss_kb": "",
        "render_wall": "",
        "render_max_rss_kb": "",
        "upload_wall": "",
        "upload_max_rss_kb": "",
        "png_url": png_url,
        "paf_source": paf_source,
    }
    if crush_summary:
        for key in ["resolved", "bailed", "candidates_seen", "rounds"]:
            fields[key] = crush_summary.get(key, "")
    if compare:
        fields["spelling_mismatches"] = compare.get("spelling_mismatches", "")
    if graph_runtime:
        fields["graph_wall"] = graph_runtime.wall_time
        fields["graph_max_rss_kb"] = graph_runtime.max_rss_kb
    if crush_runtime:
        fields["crush_wall"] = crush_runtime.wall_time
        fields["crush_max_rss_kb"] = crush_runtime.max_rss_kb
    if sort_runtime:
        fields["sort_wall"] = sort_runtime.wall_time
        fields["sort_max_rss_kb"] = sort_runtime.max_rss_kb
    if render_runtime:
        fields["render_wall"] = render_runtime.wall_time
        fields["render_max_rss_kb"] = render_runtime.max_rss_kb
    if upload_runtime:
        fields["upload_wall"] = upload_runtime.wall_time
        fields["upload_max_rss_kb"] = upload_runtime.max_rss_kb
    return fields


def write_scoreboard(state: RunState) -> None:
    if not state.rows:
        return
    fields = list(state.rows[0].keys())
    with (state.out_dir / "scoreboard.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(state.rows)


def run_seed(
    state: RunState,
    fasta: Path,
    k: int,
    threads: int,
    reuse_paf: bool,
    upload_target: str,
    public_base: str,
    no_upload: bool,
) -> tuple[Path, dict[str, str]]:
    label = f"c4.k{k}.seed"
    gfa = state.out_dir / "graphs" / f"{label}.gfa"
    debug_dir = state.out_dir / "debug" / f"k{k}"
    debug_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = Path(f"/tmp/c4lk{k}")
    command: list[object] = [
        "impg",
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
        gfa,
        "-t",
        threads,
        "-v",
        "1",
    ]
    paf_source = "fresh-fastga"
    if reuse_paf and state.raw_paf and state.raw_paf.exists():
        command = [
            "impg",
            "graph",
            "--sequence-files",
            fasta,
            "--paf-file",
            state.raw_paf,
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
            gfa,
            "-t",
            threads,
            "-v",
            "1",
        ]
        paf_source = str(state.raw_paf)
    graph_runtime = None
    if not gfa.exists():
        graph_runtime = run_timed(state, f"{label}.graph", command, cwd=Path.cwd())
    else:
        graph_runtime = state.runtimes.get(f"{label}.graph")
    if state.raw_paf is None:
        raw_paf = debug_dir / "raw.paf"
        if raw_paf.exists():
            state.raw_paf = raw_paf
    report_tsv = run_graph_report(state, label, gfa, threads)
    report = read_report_row(report_tsv)
    render_and_upload(state, label, gfa, threads, upload_target, public_base, no_upload)
    row = metric_row(
        k=k,
        x="seed",
        label=label,
        gfa=gfa,
        report=report,
        png_url=state.uploads[label],
        graph_runtime=graph_runtime,
        sort_runtime=state.runtimes.get(f"{label}.gfasort"),
        render_runtime=state.runtimes.get(f"{label}.gfalook"),
        upload_runtime=state.runtimes.get(f"{label}.upload"),
        paf_source=paf_source,
    )
    state.rows.append(row)
    write_scoreboard(state)
    return gfa, report


def run_poa(
    state: RunState,
    k: int,
    seed_gfa: Path,
    x: str,
    traversal_len: str,
    threads: int,
    compare_binary: Path,
    upload_target: str,
    public_base: str,
    no_upload: bool,
) -> tuple[Path, dict[str, str]]:
    label = f"c4.k{k}.poa{x}"
    gfa = state.out_dir / "graphs" / f"{label}.gfa"
    crush_runtime = None
    if not gfa.exists():
        crush_runtime = run_timed(
            state,
            f"{label}.crush",
            [
                "impg",
                "crush",
                "-g",
                seed_gfa,
                "-o",
                gfa,
                "--method",
                "poa",
                "--max-iterations",
                "5",
                "--max-traversal-len",
                traversal_len,
                "--max-median-traversal-len",
                traversal_len,
                "--max-total-sequence",
                "1m",
                "--max-traversals",
                "10k",
                "-t",
                threads,
                "-v",
                "1",
            ],
        )
    else:
        crush_runtime = state.runtimes.get(f"{label}.crush")
    compare = compare_paths(state, label, compare_binary, seed_gfa, gfa)
    crush_summary = parse_crush_summary(Path(state.runtimes[f"{label}.crush"].stderr))
    report_tsv = run_graph_report(state, label, gfa, threads)
    report = read_report_row(report_tsv)
    render_and_upload(state, label, gfa, threads, upload_target, public_base, no_upload)
    row = metric_row(
        k=k,
        x=x,
        label=label,
        gfa=gfa,
        report=report,
        png_url=state.uploads[label],
        crush_runtime=crush_runtime,
        sort_runtime=state.runtimes.get(f"{label}.gfasort"),
        render_runtime=state.runtimes.get(f"{label}.gfalook"),
        upload_runtime=state.runtimes.get(f"{label}.upload"),
        compare=compare,
        crush_summary=crush_summary,
    )
    state.rows.append(row)
    write_scoreboard(state)
    return gfa, report


def choose_best_seed_ks(seed_reports: dict[int, dict[str, str]], count: int = 2) -> list[int]:
    def rank_item(item: tuple[int, dict[str, str]]) -> tuple[int, int, int, int, int, int]:
        k, row = item
        return (
            int_value(row, "direct_self_loop_edges"),
            int_value(row, "adjacent_same_node_path_steps"),
            int_value(row, "path_white_space_bp_p99"),
            int_value(row, "path_white_space_bp_max"),
            int_value(row, "path_jump_p99"),
            k,
        )

    return [k for k, _row in sorted(seed_reports.items(), key=rank_item)[:count]]


def write_provenance(state: RunState, args: argparse.Namespace, ks: list[int], best_2kb: list[int]) -> None:
    payload = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "task": "c4-low-seqwish-k-sweep",
        "fasta": str(args.fasta),
        "out_dir": str(args.out_dir),
        "ks": ks,
        "threads": args.threads,
        "raw_paf_reuse": str(state.raw_paf) if state.raw_paf else "",
        "best_2kb_ks_by_seed_metrics": best_2kb,
        "recipe": {
            "graph": (
                "impg graph --sequence-files FASTA --gfa-engine seqwish --fastga "
                "--num-mappings 1:many --scaffold-filter 1:many --scaffold-jump 0 "
                "--min-match-len K"
            ),
            "poa1kb": (
                "impg crush --method poa --max-iterations 5 --max-traversal-len 1k "
                "--max-median-traversal-len 1k --max-total-sequence 1m --max-traversals 10k"
            ),
        },
    }
    (state.out_dir / "provenance.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )


def write_analysis(state: RunState) -> None:
    seed_rows = [row for row in state.rows if row["x"] == "seed"]
    poa1_rows = [row for row in state.rows if row["x"] == "1kb"]
    if not seed_rows or not poa1_rows:
        return

    def rank(row: dict[str, str]) -> tuple[int, int, int, int, int]:
        return (
            int(row.get("direct_self_loop_edges") or 0),
            int(row.get("adjacent_same_node_repeats") or 0),
            int(row.get("path_white_space_bp_p99") or 0),
            int(row.get("path_white_space_bp_max") or 0),
            int(row.get("segments") or 0),
        )

    best_seed = sorted(seed_rows, key=rank)[0]
    best_poa = sorted(poa1_rows, key=rank)[0]
    k311_seed = next((row for row in seed_rows if row["k"] == "311"), None)
    better_than_311 = [
        row
        for row in seed_rows
        if k311_seed and row["k"] != "311" and rank(row) < rank(k311_seed)
    ]
    if better_than_311:
        k311_answer = (
            "Yes by these graph-report stress metrics: at least one lower K has fewer "
            "self-loop/adjacent-repeat artifacts or lower path whitespace than K=311."
        )
    else:
        k311_answer = (
            "No by these graph-report stress metrics: K=311 is not worse than the lower "
            "seed floors in this sweep."
        )

    lines = [
        "# C4 low seqwish-k sweep analysis",
        "",
        f"Output directory: `{state.out_dir}`",
        "",
        "## Explicit answers",
        "",
        f"- Is K=311 too large for C4? {k311_answer}",
        (
            "- Do K=79 or K=127 preserve C4 architecture better while avoiding the "
            "K=1 self-loop artifact? See `scoreboard.tsv`; the seed and 1kb POA "
            "rankings below are based on direct self-loops, adjacent same-node "
            "repeats, path whitespace p99/max, and segment count."
        ),
        (
            f"- Best seed tradeoff: K={best_seed['k']} "
            f"({best_seed['direct_self_loop_edges']} direct self-loop edges, "
            f"{best_seed['adjacent_same_node_repeats']} adjacent same-node repeats, "
            f"white-space p99={best_seed['path_white_space_bp_p99']})."
        ),
        (
            f"- Best 1kb POA tradeoff: K={best_poa['k']} "
            f"({best_poa['direct_self_loop_edges']} direct self-loop edges, "
            f"{best_poa['adjacent_same_node_repeats']} adjacent same-node repeats, "
            f"white-space p99={best_poa['path_white_space_bp_p99']}, "
            f"spelling mismatches={best_poa['spelling_mismatches']})."
        ),
        "",
        "## Validation",
        "",
        "- Full-C4 seed graphs were built/rendered for every required K.",
        "- POA 1kb polished graphs were built/rendered for every required K.",
        "- `compare_gfa_paths` was run for every polished graph against its seed.",
        "",
        "## Scoreboard excerpt",
        "",
        "| K | X | segments | links | path_steps | segment_bp | self_loops | adjacent_repeats | singleton_bp | bp_cov | ws_p99 | ws_max | duplicate_frac | resolved | bailed | mismatches | PNG |",
        "|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|",
    ]
    for row in state.rows:
        lines.append(
            "| "
            + " | ".join(
                [
                    row["k"],
                    row["x"],
                    row["segments"],
                    row["links"],
                    row["path_steps"],
                    row["segment_bp"],
                    row["direct_self_loop_edges"],
                    row["adjacent_same_node_repeats"],
                    row["singleton_bp"],
                    row["bp_weighted_node_coverage"],
                    row["path_white_space_bp_p99"],
                    row["path_white_space_bp_max"],
                    row["duplicate_sequence_frac"],
                    row["resolved"],
                    row["bailed"],
                    row["spelling_mismatches"],
                    row["png_url"],
                ]
            )
            + " |"
        )
    lines.append("")
    (state.out_dir / "analysis.md").write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    args = parse_args()
    ks = [int(part) for part in args.ks.split(",") if part.strip()]
    ensure_dirs(args.out_dir)
    state = RunState(out_dir=args.out_dir, commands_path=args.out_dir / "commands.md")
    if not state.commands_path.exists():
        state.commands_path.write_text(
            "# C4 low seqwish-k sweep commands\n\n"
            f"Input FASTA: `{args.fasta}`\n\n"
            "PAF reuse policy: the first successful FastGA run writes debug/raw.paf. "
            "Later K values reuse that raw PAF through `impg graph --paf-file`, with "
            "the same 1:many mapping and scaffold filters, so only seqwish "
            "min-match-len changes.\n",
            encoding="utf-8",
        )
    if not args.fasta.exists():
        raise FileNotFoundError(args.fasta)
    if not args.compare_gfa_paths.exists():
        raise FileNotFoundError(args.compare_gfa_paths)

    seed_gfas: dict[int, Path] = {}
    seed_reports: dict[int, dict[str, str]] = {}
    for k in ks:
        gfa, report = run_seed(
            state,
            args.fasta,
            k,
            args.threads,
            args.reuse_paf,
            args.upload_target,
            args.public_base,
            args.no_upload,
        )
        seed_gfas[k] = gfa
        seed_reports[k] = report

    for k in ks:
        run_poa(
            state,
            k,
            seed_gfas[k],
            "1kb",
            "1k",
            args.threads,
            args.compare_gfa_paths,
            args.upload_target,
            args.public_base,
            args.no_upload,
        )

    best_2kb: list[int] = []
    if not args.skip_2kb:
        best_2kb = choose_best_seed_ks(seed_reports, count=2)
        for k in best_2kb:
            run_poa(
                state,
                k,
                seed_gfas[k],
                "2kb",
                "2k",
                args.threads,
                args.compare_gfa_paths,
                args.upload_target,
                args.public_base,
                args.no_upload,
            )

    write_scoreboard(state)
    write_runtime_summary(state)
    write_uploaded_urls(state)
    write_validation_summary(state)
    write_provenance(state, args, ks, best_2kb)
    write_analysis(state)
    return 0


if __name__ == "__main__":
    sys.exit(main())
