#!/usr/bin/env python3
"""Run the C4 small-bubble POA/POASTA threshold sweep."""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import re
import shlex
import subprocess
import sys
from datetime import datetime, timezone


ROOT = Path("/home/erikg/impg/.wg-worktrees/agent-450/data/c4_threshold_sweep_20260603T154441Z")
SEED = Path(
    "/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/"
    "one_many_minmatch1_scaffold0.initial.gfa"
)
FASTA = Path("/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa")
IMPG = Path("/home/erikg/.cargo/bin/impg")
TIME = Path("/usr/bin/time")
EXPECTED_PATHS = 465
THREADS = 32
UPLOAD_HOST = "erik@hypervolu.me"
UPLOAD_DIR = "www/impg"

MATRIX = [
    ("poa", 100),
    ("poa", 500),
    ("poa", 1000),
    ("poa", 2000),
    ("poasta", 1000),
    ("poasta", 2500),
    ("poasta", 5000),
]

REPORT_KEYS = [
    "segments",
    "links",
    "paths",
    "path_steps",
    "total_segment_bp",
    "node_path_step_visits",
    "node_coverage_mean",
    "node_coverage_bp_weighted_mean",
    "node_coverage_p10",
    "node_coverage_median",
    "node_coverage_p90",
    "singleton_nodes",
    "singleton_bp",
    "high_coverage_threshold",
    "high_coverage_nodes",
    "high_coverage_bp",
    "link_jump_p99",
    "link_jump_max",
    "path_jump_p99",
    "path_jump_max",
    "path_white_space_bp_p99",
    "path_white_space_bp_max",
    "path_white_space_bridges_ge_threshold",
    "segment_white_space_bp_fraction",
    "segment_white_space_bp_total",
    "sparse_coverage_segment_bp",
    "path_depth_median",
    "path_depth_p95",
    "path_depth_max",
    "reused_nodes",
    "duplicate_sequence_frac",
    "local_repeat_context_nodes",
    "local_repeat_context_occurrences",
]


def log(message: str) -> None:
    ROOT.mkdir(parents=True, exist_ok=True)
    stamp = datetime.now(timezone.utc).isoformat()
    line = f"[{stamp}] {message}"
    print(line, flush=True)
    with (ROOT / "driver.log").open("a") as handle:
        handle.write(line + "\n")


def run_logged(cmd: list[str], stdout_path: Path, stderr_path: Path, cwd: Path | None = None) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    log("RUN " + shlex.join(cmd))
    with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
        proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, stdout=stdout, stderr=stderr)
    log(f"EXIT {proc.returncode} " + shlex.join(cmd))
    return proc.returncode


def timed(cmd: list[str]) -> list[str]:
    return [str(TIME), "-v", *cmd]


def revcomp(seq: bytes) -> bytes:
    table = bytes.maketrans(b"ACGTNacgtn", b"TGCANtgcan")
    return seq.translate(table)[::-1]


def parse_step_token(token: str) -> tuple[str, bool]:
    token = token.strip()
    if not token:
        raise ValueError("empty path step token")
    if token[-1] in "+-":
        return token[:-1], token[-1] == "-"
    if token[0] in "><":
        return token[1:], token[0] == "<"
    raise ValueError(f"cannot parse path step orientation: {token!r}")


def gfa_path_signatures(gfa: Path) -> dict[str, dict[str, str]]:
    segments: dict[str, bytes] = {}
    path_steps: list[tuple[str, list[tuple[str, bool]]]] = []
    with gfa.open("rb") as handle:
        for raw in handle:
            if raw.startswith(b"S\t"):
                parts = raw.rstrip(b"\n").split(b"\t")
                if len(parts) >= 3:
                    segments[parts[1].decode()] = parts[2]
            elif raw.startswith(b"P\t"):
                parts = raw.rstrip(b"\n").split(b"\t")
                if len(parts) >= 3:
                    name = parts[1].decode()
                    tokens = parts[2].decode().split(",")
                    path_steps.append((name, [parse_step_token(token) for token in tokens if token]))
    out: dict[str, dict[str, str]] = {}
    for name, steps in path_steps:
        digest = hashlib.sha256()
        length = 0
        for node, rev in steps:
            seq = segments[node]
            if rev:
                seq = revcomp(seq)
            digest.update(seq)
            length += len(seq)
        out[name] = {
            "sha256": digest.hexdigest(),
            "length": str(length),
            "steps": str(len(steps)),
        }
    return out


def fasta_signatures(fasta: Path) -> dict[str, dict[str, str]]:
    out: dict[str, dict[str, str]] = {}
    name: str | None = None
    digest: hashlib._Hash | None = None
    length = 0
    with fasta.open("rb") as handle:
        for raw in handle:
            line = raw.rstrip(b"\n\r")
            if line.startswith(b">"):
                if name is not None and digest is not None:
                    out[name] = {"sha256": digest.hexdigest(), "length": str(length)}
                name = line[1:].split()[0].decode()
                digest = hashlib.sha256()
                length = 0
            elif name is not None and digest is not None:
                seq = line.strip()
                digest.update(seq.upper())
                length += len(seq)
    if name is not None and digest is not None:
        out[name] = {"sha256": digest.hexdigest(), "length": str(length)}
    return out


def validate_paths(
    gfa: Path, seed_paths: dict[str, dict[str, str]], fasta_paths: dict[str, dict[str, str]]
) -> dict[str, str]:
    if not gfa.exists():
        return {
            "path_validation_status": "missing_gfa",
            "path_count": "0",
            "path_count_ok": "false",
            "path_names_match_seed": "false",
            "path_spelling_mismatches_vs_seed": "",
            "missing_path_names": "",
            "extra_path_names": "",
            "fasta_name_set_ok": "false" if fasta_paths else "",
            "fasta_spelling_mismatches": "",
        }
    try:
        paths = gfa_path_signatures(gfa)
    except Exception as exc:
        return {
            "path_validation_status": f"parse_error:{exc}",
            "path_count": "",
            "path_count_ok": "false",
            "path_names_match_seed": "false",
            "path_spelling_mismatches_vs_seed": "",
            "missing_path_names": "",
            "extra_path_names": "",
            "fasta_name_set_ok": "",
            "fasta_spelling_mismatches": "",
        }
    names = set(paths)
    seed_names = set(seed_paths)
    missing = sorted(seed_names - names)
    extra = sorted(names - seed_names)
    spelling = sorted(
        name
        for name in names & seed_names
        if paths[name]["sha256"] != seed_paths[name]["sha256"]
        or paths[name]["length"] != seed_paths[name]["length"]
    )
    fasta_name_set_ok = ""
    fasta_spelling = []
    if fasta_paths:
        fasta_names = set(fasta_paths)
        fasta_name_set_ok = str(names == fasta_names).lower()
        fasta_spelling = sorted(
            name
            for name in names & fasta_names
            if paths[name]["sha256"] != fasta_paths[name]["sha256"]
            or paths[name]["length"] != fasta_paths[name]["length"]
        )
    return {
        "path_validation_status": "ok",
        "path_count": str(len(paths)),
        "path_count_ok": str(len(paths) == EXPECTED_PATHS).lower(),
        "path_names_match_seed": str(not missing and not extra).lower(),
        "path_spelling_mismatches_vs_seed": str(len(spelling)),
        "missing_path_names": ",".join(missing[:20]),
        "extra_path_names": ",".join(extra[:20]),
        "fasta_name_set_ok": fasta_name_set_ok,
        "fasta_spelling_mismatches": str(len(fasta_spelling)) if fasta_paths else "",
    }


def parse_time_log(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(errors="replace").splitlines():
        if "Elapsed (wall clock) time" in line:
            match = re.search(r"\):\s*(.+)$", line)
            out["wall"] = match.group(1).strip() if match else line.rsplit(":", 1)[-1].strip()
        elif "Maximum resident set size" in line:
            out["max_rss_kb"] = line.rsplit(":", 1)[1].strip()
        elif "Exit status" in line:
            out["exit_status"] = line.rsplit(":", 1)[1].strip()
    return out


def parse_crush_log(path: Path) -> dict[str, str]:
    text = path.read_text(errors="replace") if path.exists() else ""
    rounds: dict[int, str] = {}
    no_eligible: list[str] = []
    for line in text.splitlines():
        match = re.search(r"crush round (\d+): resolved (\d+)/(\d+) replacement", line)
        if match:
            rounds[int(match.group(1))] = f"{match.group(2)}/{match.group(3)}"
            continue
        match = re.search(r"crush round (\d+): (\d+) resolved \(from sites", line)
        if match:
            rounds[int(match.group(1))] = match.group(2)
            continue
        match = re.search(r"crush iterative-multi-level round (\d+): applied (\d+) candidate", line)
        if match:
            rounds[int(match.group(1))] = match.group(2)
            continue
        match = re.search(r"crush round (\d+): no eligible candidates", line)
        if match:
            no_eligible.append(match.group(1))
    total = re.search(r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds", text)
    return {
        "resolved_sites_per_round": ",".join(f"r{r}={rounds[r]}" for r in sorted(rounds)),
        "no_eligible_rounds": ",".join(no_eligible),
        "resolved_total": total.group(1) if total else "",
        "bailed_total": total.group(2) if total else "",
        "candidates_seen": total.group(3) if total else "",
        "rounds": total.group(4) if total else "",
    }


def read_report_tsv(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return dict(rows[0]) if rows else {}


def write_tsv(path: Path, rows: list[dict[str, str]]) -> None:
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=keys)
        writer.writeheader()
        writer.writerows(rows)


def write_json(path: Path, payload: object) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n")


def crush_command(method: str, threshold: int, output: Path) -> list[str]:
    return [
        str(IMPG),
        "crush",
        "--gfa",
        str(SEED),
        "--output",
        str(output),
        "--method",
        method,
        "--max-iterations",
        "20",
        "--max-span",
        "0",
        "--max-traversal-len",
        "100k",
        "--max-median-traversal-len",
        str(threshold),
        "--max-total-sequence",
        "500m",
        "--max-traversals",
        "100k",
        "--poa-scoring",
        "1,4,6,2,26,1",
        "--threads",
        str(THREADS),
        "-v",
        "1",
    ]


def graph_report_command(gfa: Path, output: Path) -> list[str]:
    return [
        str(IMPG),
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
        str(THREADS),
        "-v",
        "1",
    ]


def main() -> int:
    ROOT.mkdir(parents=True, exist_ok=True)
    log(f"output root: {ROOT}")
    log(f"seed: {SEED}")
    log(f"impg: {IMPG}")
    seed_paths = gfa_path_signatures(SEED)
    fasta_paths = fasta_signatures(FASTA) if FASTA.exists() else {}
    seed_validation = validate_paths(SEED, seed_paths, fasta_paths)
    write_json(ROOT / "seed_path_validation.json", seed_validation)
    write_json(
        ROOT / "manifest.json",
        {
            "root": str(ROOT),
            "seed": str(SEED),
            "fasta": str(FASTA),
            "impg": str(IMPG),
            "expected_paths": EXPECTED_PATHS,
            "threads": THREADS,
            "matrix": [{"method": m, "threshold": t} for m, t in MATRIX],
            "seed_path_validation": seed_validation,
            "created_utc": datetime.now(timezone.utc).isoformat(),
        },
    )

    rows: list[dict[str, str]] = []
    validations: list[dict[str, str]] = []
    uploads: list[dict[str, str]] = []

    for method, threshold in MATRIX:
        name = f"{method}_median{threshold}"
        run_dir = ROOT / name
        run_dir.mkdir(parents=True, exist_ok=True)
        out_gfa = run_dir / f"{name}.gfa"
        sorted_gfa = run_dir / f"{name}.Ygs.gfa"
        png = run_dir / f"c4-smallbubble-{method}-median{threshold}.png"
        remote_name = f"c4-smallbubble-{method}-median{threshold}.png"
        report = run_dir / "graph-report.tsv"

        row: dict[str, str] = {
            "name": name,
            "method": method,
            "max_median_traversal_len": str(threshold),
            "gfa": str(out_gfa),
            "sorted_gfa": str(sorted_gfa),
            "png": str(png),
            "png_remote_name": remote_name,
        }

        cmd = crush_command(method, threshold, out_gfa)
        (run_dir / "command.sh").write_text("#!/bin/sh\nset -eu\n" + shlex.join(cmd) + "\n")
        os.chmod(run_dir / "command.sh", 0o755)
        rc = run_logged(timed(cmd), run_dir / "crush.stdout.log", run_dir / "crush.stderr.log")
        row["crush_rc"] = str(rc)
        row.update(parse_time_log(run_dir / "crush.stderr.log"))
        row.update(parse_crush_log(run_dir / "crush.stderr.log"))

        validation = {"name": name, "method": method, "max_median_traversal_len": str(threshold)}
        validation.update(validate_paths(out_gfa, seed_paths, fasta_paths))
        validations.append(validation)
        row.update(validation)

        if out_gfa.exists():
            rc = run_logged(
                timed(graph_report_command(out_gfa, report)),
                run_dir / "graph-report.stdout.log",
                run_dir / "graph-report.stderr.log",
            )
            row["graph_report_rc"] = str(rc)
            for key, value in parse_time_log(run_dir / "graph-report.stderr.log").items():
                row[f"graph_report_{key}"] = value
            metrics = read_report_tsv(report)
            for key in REPORT_KEYS:
                row[key] = metrics.get(key, "")

            rc = run_logged(
                timed(["gfasort", "-i", str(out_gfa), "-o", str(sorted_gfa), "-p", "Ygs", "-t", str(THREADS)]),
                run_dir / "gfasort.stdout.log",
                run_dir / "gfasort.stderr.log",
            )
            row["gfasort_rc"] = str(rc)
            for key, value in parse_time_log(run_dir / "gfasort.stderr.log").items():
                row[f"gfasort_{key}"] = value

            if sorted_gfa.exists():
                rc = run_logged(
                    timed(["gfalook", "-i", str(sorted_gfa), "-o", str(png), "-m", "-x", "2200", "-y", "1200"]),
                    run_dir / "gfalook.stdout.log",
                    run_dir / "gfalook.stderr.log",
                )
                row["gfalook_rc"] = str(rc)
                for key, value in parse_time_log(run_dir / "gfalook.stderr.log").items():
                    row[f"gfalook_{key}"] = value
                if png.exists():
                    remote = f"{UPLOAD_HOST}:{UPLOAD_DIR}/{remote_name}"
                    rc = run_logged(["scp", str(png), remote], run_dir / "scp.stdout.log", run_dir / "scp.stderr.log")
                    url = f"https://hypervolu.me/~erik/impg/{remote_name}" if rc == 0 else ""
                    row["scp_rc"] = str(rc)
                    row["png_url"] = url
                    uploads.append(
                        {
                            "name": name,
                            "png": str(png),
                            "remote_name": remote_name,
                            "scp_rc": str(rc),
                            "url": url,
                        }
                    )
                else:
                    row["gfalook_png_exists"] = "false"
            else:
                row["sorted_gfa_exists"] = "false"
        else:
            row["output_gfa_exists"] = "false"

        rows.append(row)
        write_tsv(ROOT / "metrics.tsv", rows)
        write_tsv(ROOT / "path-validation.tsv", validations)
        write_tsv(ROOT / "uploads.tsv", uploads)

    write_tsv(ROOT / "metrics.tsv", rows)
    write_tsv(ROOT / "path-validation.tsv", validations)
    write_tsv(ROOT / "uploads.tsv", uploads)
    log("sweep complete")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
