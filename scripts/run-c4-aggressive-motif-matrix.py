#!/usr/bin/env python3
"""Run the full-C4 aggressive motif-window POA/abPOA comparison.

The heavy graph outputs are intentionally written outside the repository. This
driver keeps the comparison reproducible and writes compact TSV summaries that
can be committed or referenced from an evaluation report.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import shutil
import signal
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


DEFAULT_INPUT = Path(
    "/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.poa2kb.gfa"
)
DEFAULT_OUT = Path("/home/erikg/impg/data/c4_aggressive_motif_matrix_20260605T210000Z")
DEFAULT_PUBLIC_BASE = "http://hypervolu.me/~erik/impg/"
DEFAULT_UPLOAD_TARGET = "erik@hypervolu.me:www/impg/"
GRAPH_REPORT_FIELDS = [
    "segments",
    "links",
    "paths",
    "path_steps",
    "total_segment_bp",
    "singleton_bp",
    "direct_self_loop_edges",
    "adjacent_same_node_path_steps",
    "self_loop_repeat_runs",
    "node_coverage_bp_weighted_mean",
    "path_white_space_bp_p99",
    "path_white_space_bp_max",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-gfa", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--impg", type=Path, default=Path("target/release/impg"))
    parser.add_argument(
        "--compare-bin",
        type=Path,
        default=Path("target/release/examples/compare_gfa_paths"),
    )
    parser.add_argument("--gfasort", default="gfasort")
    parser.add_argument("--gfalook", default="gfalook")
    parser.add_argument("--abpoa-bin", default="abpoa")
    parser.add_argument("--threads", type=int, default=32)
    parser.add_argument("--methods", default="poa,abpoa")
    parser.add_argument("--windows", default="10k,25k")
    parser.add_argument("--candidate-limit", type=int, default=512)
    parser.add_argument("--max-iterations", type=int, default=8)
    parser.add_argument(
        "--poa-scoring",
        default="1,4,4,2,24,1",
        help=(
            "match,mismatch,gap_open1,gap_extend1,gap_open2,gap_extend2. "
            "Default matches abPOA 1.5.x convex-gap defaults apart from "
            "making the match/mismatch values explicit."
        ),
    )
    parser.add_argument("--timeout-minutes", type=float, default=90.0)
    parser.add_argument("--public-base", default=DEFAULT_PUBLIC_BASE)
    parser.add_argument("--upload-target", default=DEFAULT_UPLOAD_TARGET)
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--no-run", action="store_true")
    parser.add_argument("--no-render", action="store_true")
    parser.add_argument("--no-upload", action="store_true")
    return parser.parse_args()


def ensure_dirs(out_dir: Path) -> dict[str, Path]:
    dirs = {
        "graphs": out_dir / "graphs",
        "reports": out_dir / "reports",
        "logs": out_dir / "logs",
        "renders": out_dir / "renders",
        "bin": out_dir / "bin",
    }
    for path in dirs.values():
        path.mkdir(parents=True, exist_ok=True)
    return dirs


def size_label(raw: str) -> str:
    return raw.lower().replace("_", "").replace(" ", "")


def sh(cmd: list[object]) -> str:
    return shlex.join(str(part) for part in cmd)


def append_command(out_dir: Path, label: str, cmd: list[object]) -> None:
    with (out_dir / "commands.sh").open("a", encoding="utf-8") as handle:
        handle.write(f"# {label}\n")
        handle.write(sh(cmd))
        handle.write("\n\n")


def run_command(
    label: str,
    cmd: list[object],
    stdout_path: Path,
    stderr_path: Path,
    time_path: Path | None,
    timeout_s: float | None,
    env: dict[str, str] | None = None,
) -> int:
    stdout_path.parent.mkdir(parents=True, exist_ok=True)
    stderr_path.parent.mkdir(parents=True, exist_ok=True)
    run_cmd = [str(part) for part in cmd]
    if time_path is not None:
        time_path.parent.mkdir(parents=True, exist_ok=True)
        run_cmd = ["/usr/bin/time", "-v", "-o", str(time_path), *run_cmd]
    with stdout_path.open("wb") as stdout, stderr_path.open("ab") as stderr:
        stderr.write(f"\n# RUN {label}\n# CMD {sh(cmd)}\n".encode())
        try:
            proc = subprocess.Popen(
                run_cmd,
                stdout=stdout,
                stderr=stderr,
                start_new_session=True,
                env=env,
            )
            try:
                return proc.wait(timeout=timeout_s)
            except subprocess.TimeoutExpired:
                os.killpg(proc.pid, signal.SIGTERM)
                try:
                    return proc.wait(timeout=15)
                except subprocess.TimeoutExpired:
                    os.killpg(proc.pid, signal.SIGKILL)
                    proc.wait()
                    stderr.write(f"\n# TIMEOUT {timeout_s:.1f}s; killed process group\n".encode())
                    return 124
        except FileNotFoundError as err:
            stderr.write(f"\n# EXEC ERROR {err}\n".encode())
            return 127


def parse_elapsed_seconds(raw: str) -> float | None:
    raw = raw.strip()
    if not raw:
        return None
    parts = raw.split(":")
    try:
        if len(parts) == 3:
            return int(parts[0]) * 3600 + int(parts[1]) * 60 + float(parts[2])
        if len(parts) == 2:
            return int(parts[0]) * 60 + float(parts[1])
        return float(raw)
    except ValueError:
        return None


def parse_time_log(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not path.exists():
        return result
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        line = line.strip()
        if line.startswith("Elapsed (wall clock) time"):
            if "):" in line:
                value = line.split("):", 1)[1].strip()
            else:
                value = line.rsplit(":", 1)[-1].strip()
            seconds = parse_elapsed_seconds(value)
            if seconds is not None:
                result["wall_s"] = f"{seconds:.2f}"
            result["wall_raw"] = value
        elif line.startswith("Maximum resident set size"):
            result["max_rss_kb"] = line.rsplit(":", 1)[-1].strip()
        elif line.startswith("Exit status"):
            result["time_exit_status"] = line.rsplit(":", 1)[-1].strip()
    return result


def read_tsv_row(path: Path) -> dict[str, str]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", errors="replace", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return dict(row)
    return {}


def parse_compare_stdout(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not path.exists():
        return result
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = line.strip().split()
        if len(parts) == 2:
            result[parts[0]] = parts[1]
    return result


def parse_crush_log(path: Path) -> dict[str, str]:
    text = path.read_text(encoding="utf-8", errors="replace") if path.exists() else ""
    result: dict[str, str] = {}
    summary = re.search(
        r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds",
        text,
    )
    if summary:
        result.update(
            {
                "resolved": summary.group(1),
                "bailed": summary.group(2),
                "candidates_seen": summary.group(3),
                "rounds": summary.group(4),
            }
        )
    generated: list[str] = []
    applied: list[str] = []
    for match in re.finditer(
        r"crush iterative-multi-level round (\d+): .*?(\d+) generated candidate\(s\), (\d+) considered",
        text,
    ):
        generated.append(f"r{match.group(1)}={match.group(2)}/{match.group(3)}")
    for match in re.finditer(
        r"crush iterative-multi-level round (\d+): applied (\d+) candidate\(s\)",
        text,
    ):
        applied.append(f"r{match.group(1)}={match.group(2)}")
    if generated:
        result["generated_by_round"] = ",".join(generated)
    if applied:
        result["applied_by_round"] = ",".join(applied)
    if "no generated candidate(s)" in text:
        result["stopped_reason"] = "stable:no-generated-candidates"
    elif "0 applicable candidate(s)" in text:
        result["stopped_reason"] = "stable:no-applicable-candidates"
    elif summary and summary.group(4) == "8":
        result["stopped_reason"] = "max-iterations"
    return result


def token_node(token: str) -> str:
    token = token.strip()
    if token and token[-1] in "+-":
        return token[:-1]
    return token


def parse_gfa_repeats(path: Path, input_segments: set[str]) -> list[dict[str, str]]:
    segments: dict[str, str] = {}
    direct_loops: Counter[str] = Counter()
    adjacent_steps: Counter[str] = Counter()
    repeat_runs: Counter[str] = Counter()
    max_run: Counter[str] = Counter()
    paths_with_repeat: dict[str, set[str]] = defaultdict(set)

    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if not fields:
                continue
            if fields[0] == "S":
                segments[fields[1]] = fields[2]
            elif fields[0] == "L" and len(fields) >= 4 and fields[1] == fields[3]:
                direct_loops[fields[1]] += 1
            elif fields[0] == "P" and len(fields) >= 3:
                path_name = fields[1]
                nodes = [token_node(token) for token in fields[2].split(",") if token]
                idx = 0
                while idx < len(nodes):
                    end = idx + 1
                    while end < len(nodes) and nodes[end] == nodes[idx]:
                        end += 1
                    run_len = end - idx
                    if run_len > 1:
                        node = nodes[idx]
                        adjacent_steps[node] += run_len - 1
                        repeat_runs[node] += 1
                        max_run[node] = max(max_run[node], run_len)
                        paths_with_repeat[node].add(path_name)
                    idx = end

    nodes = set(direct_loops) | set(adjacent_steps)
    rows: list[dict[str, str]] = []
    for node in nodes:
        seq = segments.get(node, "")
        rows.append(
            {
                "node": node,
                "len": str(len(seq)),
                "seq_prefix": seq[:48],
                "direct_self_loop_edges": str(direct_loops[node]),
                "adjacent_same_node_path_steps": str(adjacent_steps[node]),
                "self_loop_repeat_runs": str(repeat_runs[node]),
                "max_repeat_run_len": str(max_run[node]),
                "paths_with_repeat": str(len(paths_with_repeat.get(node, set()))),
                "in_input": "yes" if node in input_segments else "no",
            }
        )
    rows.sort(
        key=lambda row: (
            -int(row["direct_self_loop_edges"]),
            -int(row["adjacent_same_node_path_steps"]),
            -int(row["self_loop_repeat_runs"]),
            row["node"],
        )
    )
    return rows


def read_segment_names(path: Path) -> set[str]:
    names: set[str] = set()
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("S\t"):
                parts = line.split("\t", 3)
                if len(parts) >= 2:
                    names.add(parts[1])
    return names


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def graph_report(args: argparse.Namespace, dirs: dict[str, Path], label: str, gfa: Path) -> Path:
    out = dirs["reports"] / f"{label}.graph-report.tsv"
    if out.exists() and not args.force:
        return out
    cmd = [
        args.impg,
        "graph-report",
        "-g",
        gfa,
        "-o",
        out,
        "--format",
        "tsv",
        "--povu",
        "--top",
        "20",
        "-t",
        args.threads,
        "-v",
        "1",
    ]
    append_command(args.out_dir, f"{label}.graph-report", cmd)
    rc = run_command(
        f"{label}.graph-report",
        cmd,
        dirs["logs"] / f"{label}.graph-report.stdout.log",
        dirs["logs"] / f"{label}.graph-report.stderr.log",
        dirs["logs"] / f"{label}.graph-report.time.txt",
        args.timeout_minutes * 60,
    )
    if rc != 0:
        raise RuntimeError(f"graph-report failed for {label} with exit {rc}")
    return out


def compare_paths(args: argparse.Namespace, dirs: dict[str, Path], label: str, gfa: Path) -> int:
    cmd = [args.compare_bin, args.input_gfa, gfa]
    append_command(args.out_dir, f"{label}.compare_gfa_paths", cmd)
    return run_command(
        f"{label}.compare_gfa_paths",
        cmd,
        dirs["reports"] / f"{label}.compare_gfa_paths.stdout.txt",
        dirs["reports"] / f"{label}.compare_gfa_paths.stderr.txt",
        dirs["logs"] / f"{label}.compare_gfa_paths.time.txt",
        args.timeout_minutes * 60,
    )


def make_abpoa_wrapper(args: argparse.Namespace, dirs: dict[str, Path]) -> tuple[Path, Path]:
    real_abpoa = shutil.which(args.abpoa_bin) if not os.path.isabs(args.abpoa_bin) else args.abpoa_bin
    if not real_abpoa:
        raise FileNotFoundError(f"abPOA binary not found: {args.abpoa_bin}")
    real = Path(real_abpoa).resolve()
    log = dirs["reports"] / "abpoa.invocations.tsv"
    wrapper = dirs["bin"] / "abpoa-log-wrapper.sh"
    wrapper.write_text(
        "#!/usr/bin/env bash\n"
        "set -euo pipefail\n"
        "printf '%s\\t%s\\n' \"$(date -u +%FT%TZ)\" \"$*\" >> \"${ABPOA_INVOCATION_LOG}\"\n"
        f"exec {shlex.quote(str(real))} \"$@\"\n",
        encoding="utf-8",
    )
    wrapper.chmod(0o755)
    return wrapper, log


def run_variant(
    args: argparse.Namespace,
    dirs: dict[str, Path],
    method: str,
    window: str,
    abpoa_wrapper: Path,
    abpoa_log: Path,
    input_segments: set[str],
) -> dict[str, str]:
    label = f"c4.k311.{method}.motif{size_label(window)}.cand{args.candidate_limit}.iter{args.max_iterations}"
    gfa = dirs["graphs"] / f"{label}.gfa"
    stdout = dirs["logs"] / f"{label}.crush.stdout.log"
    stderr = dirs["logs"] / f"{label}.crush.stderr.log"
    time_log = dirs["logs"] / f"{label}.crush.time.txt"
    row: dict[str, str] = {
        "label": label,
        "method": method,
        "motif_window_bp": size_label(window),
        "candidate_limit": str(args.candidate_limit),
        "max_iterations": str(args.max_iterations),
        "poa_scoring": args.poa_scoring,
        "gfa": str(gfa),
    }
    cmd = [
        args.impg,
        "crush",
        "-g",
        args.input_gfa,
        "-o",
        gfa,
        "--method",
        method,
        "--window-mode",
        "motif",
        "--max-iterations",
        args.max_iterations,
        "--motif-max-window-bp",
        window,
        "--candidate-limit",
        args.candidate_limit,
        "--poa-scoring",
        args.poa_scoring,
        "--polish-method",
        "poa",
        "--polish-rounds",
        "until-done",
        "-t",
        args.threads,
        "-v",
        "1",
    ]
    env = os.environ.copy()
    if method == "abpoa":
        cmd.extend(["--abpoa-bin", abpoa_wrapper])
        env["ABPOA_INVOCATION_LOG"] = str(abpoa_log)
    append_command(args.out_dir, f"{label}.crush", cmd)
    if args.no_run:
        row["exit_status"] = "not-run"
    elif gfa.exists() and not args.force:
        row["exit_status"] = "existing"
    else:
        rc = run_command(
            f"{label}.crush",
            cmd,
            stdout,
            stderr,
            time_log,
            args.timeout_minutes * 60,
            env=env,
        )
        row["exit_status"] = str(rc)
    row.update(parse_time_log(time_log))
    row.update(parse_crush_log(stderr))

    if gfa.exists():
        compare_rc = compare_paths(args, dirs, label, gfa)
        row["compare_exit_status"] = str(compare_rc)
        compare = parse_compare_stdout(dirs["reports"] / f"{label}.compare_gfa_paths.stdout.txt")
        row.update({f"compare_{key}": value for key, value in compare.items()})
        row["path_preserving"] = (
            "yes"
            if compare_rc == 0
            and compare.get("missing_paths") == "0"
            and compare.get("extra_paths") == "0"
            and compare.get("spelling_mismatches") == "0"
            else "no"
        )
    else:
        row["compare_exit_status"] = "no-output"
        row["path_preserving"] = "no"

    if row.get("path_preserving") == "yes":
        report = graph_report(args, dirs, label, gfa)
        metrics = read_tsv_row(report)
        for field in GRAPH_REPORT_FIELDS:
            row[field] = metrics.get(field, "")
        repeats = parse_gfa_repeats(gfa, input_segments)
        repeat_fields = [
            "label",
            "node",
            "len",
            "seq_prefix",
            "direct_self_loop_edges",
            "adjacent_same_node_path_steps",
            "self_loop_repeat_runs",
            "max_repeat_run_len",
            "paths_with_repeat",
            "in_input",
            "log_evidence",
        ]
        log_text = stderr.read_text(encoding="utf-8", errors="replace") if stderr.exists() else ""
        repeat_rows = []
        for item in repeats[:50]:
            item = dict(item)
            item["label"] = label
            node = item["node"]
            if item["in_input"] == "no":
                item["log_evidence"] = "replacement-introduced"
            elif node in log_text:
                item["log_evidence"] = "node-id-emitted-in-crush-log"
            else:
                item["log_evidence"] = "not-emitted-in-top-candidate-log"
            repeat_rows.append(item)
        write_tsv(dirs["reports"] / f"{label}.remaining-loops.tsv", repeat_rows, repeat_fields)
    return row


def row_int(row: dict[str, str], key: str, default: int = 10**18) -> int:
    try:
        return int(float(row.get(key, "")))
    except ValueError:
        return default


def choose_best(rows: list[dict[str, str]], method: str) -> dict[str, str] | None:
    eligible = [row for row in rows if row.get("method") == method and row.get("path_preserving") == "yes"]
    if not eligible:
        return None
    return sorted(
        eligible,
        key=lambda row: (
            row_int(row, "direct_self_loop_edges"),
            row_int(row, "adjacent_same_node_path_steps"),
            row_int(row, "self_loop_repeat_runs"),
            row_int(row, "singleton_bp"),
            float(row.get("wall_s") or "1e18"),
        ),
    )[0]


def render_upload_best(
    args: argparse.Namespace,
    dirs: dict[str, Path],
    rows: list[dict[str, str]],
) -> list[dict[str, str]]:
    rendered: list[dict[str, str]] = []
    if args.no_render:
        return rendered
    for method in ["poa", "abpoa"]:
        best = choose_best(rows, method)
        if best is None:
            continue
        label = best["label"]
        gfa = Path(best["gfa"])
        sorted_gfa = dirs["graphs"] / f"{label}.Ygs.gfa"
        png = dirs["renders"] / f"{label}.Ygs.gfalook-m.png"
        sort_cmd = [args.gfasort, "-i", gfa, "-o", sorted_gfa, "-p", "Ygs", "-t", args.threads, "-v", "1"]
        look_cmd = [args.gfalook, "-i", sorted_gfa, "-o", png, "-m"]
        append_command(args.out_dir, f"{label}.gfasort", sort_cmd)
        sort_rc = run_command(
            f"{label}.gfasort",
            sort_cmd,
            dirs["logs"] / f"{label}.gfasort.stdout.log",
            dirs["logs"] / f"{label}.gfasort.stderr.log",
            dirs["logs"] / f"{label}.gfasort.time.txt",
            args.timeout_minutes * 60,
        )
        render_row = {
            "method": method,
            "label": label,
            "gfa": str(gfa),
            "sorted_gfa": str(sorted_gfa),
            "png": str(png),
            "sort_exit_status": str(sort_rc),
            "look_exit_status": "not-run",
            "upload_exit_status": "not-run",
            "url": "",
            "render_status": "ok",
        }
        if sort_rc != 0:
            render_row["render_status"] = "gfasort-failed"
            rendered.append(render_row)
            continue
        append_command(args.out_dir, f"{label}.gfalook", look_cmd)
        look_rc = run_command(
            f"{label}.gfalook",
            look_cmd,
            dirs["logs"] / f"{label}.gfalook.stdout.log",
            dirs["logs"] / f"{label}.gfalook.stderr.log",
            dirs["logs"] / f"{label}.gfalook.time.txt",
            args.timeout_minutes * 60,
        )
        render_row["look_exit_status"] = str(look_rc)
        if look_rc != 0:
            render_row["render_status"] = "gfalook-failed"
            rendered.append(render_row)
            continue
        upload_rc = "not-run"
        url = f"{args.public_base.rstrip('/')}/{png.name}"
        render_row["url"] = url
        if not args.no_upload:
            upload_cmd = ["scp", png, args.upload_target]
            append_command(args.out_dir, f"{label}.upload", upload_cmd)
            upload_rc = str(
                run_command(
                    f"{label}.upload",
                    upload_cmd,
                    dirs["logs"] / f"{label}.upload.stdout.log",
                    dirs["logs"] / f"{label}.upload.stderr.log",
                    dirs["logs"] / f"{label}.upload.time.txt",
                    args.timeout_minutes * 60,
                )
            )
        render_row["upload_exit_status"] = upload_rc
        if upload_rc not in {"0", "not-run"}:
            render_row["render_status"] = "upload-failed"
        rendered.append(render_row)
    write_tsv(
        dirs["reports"] / "rendered-best.tsv",
        rendered,
        [
            "method",
            "label",
            "gfa",
            "sorted_gfa",
            "png",
            "sort_exit_status",
            "look_exit_status",
            "upload_exit_status",
            "url",
            "render_status",
        ],
    )
    return rendered


def main() -> int:
    args = parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    dirs = ensure_dirs(args.out_dir)
    (args.out_dir / "commands.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n\n", encoding="utf-8")
    (args.out_dir / "commands.sh").chmod(0o755)

    if not args.input_gfa.exists():
        raise FileNotFoundError(args.input_gfa)
    if not args.impg.exists():
        raise FileNotFoundError(args.impg)
    if not args.compare_bin.exists():
        raise FileNotFoundError(args.compare_bin)

    input_segments = read_segment_names(args.input_gfa)
    before_report = graph_report(args, dirs, "before", args.input_gfa)
    before_repeats = parse_gfa_repeats(args.input_gfa, input_segments)
    write_tsv(
        dirs["reports"] / "before.remaining-loops.tsv",
        [dict(row, label="before", log_evidence="input") for row in before_repeats[:50]],
        [
            "label",
            "node",
            "len",
            "seq_prefix",
            "direct_self_loop_edges",
            "adjacent_same_node_path_steps",
            "self_loop_repeat_runs",
            "max_repeat_run_len",
            "paths_with_repeat",
            "in_input",
            "log_evidence",
        ],
    )

    abpoa_wrapper, abpoa_log = make_abpoa_wrapper(args, dirs)
    methods = [item.strip() for item in args.methods.split(",") if item.strip()]
    windows = [item.strip() for item in args.windows.split(",") if item.strip()]
    rows: list[dict[str, str]] = []
    for method in methods:
        if method not in {"poa", "abpoa"}:
            raise ValueError(f"unsupported method for this matrix: {method}")
        for window in windows:
            rows.append(
                run_variant(
                    args,
                    dirs,
                    method,
                    window,
                    abpoa_wrapper,
                    abpoa_log,
                    input_segments,
                )
            )
            write_summary(args, dirs, before_report, rows)
    rendered = render_upload_best(args, dirs, rows)
    write_summary(args, dirs, before_report, rows, rendered)
    return 0


def write_summary(
    args: argparse.Namespace,
    dirs: dict[str, Path],
    before_report: Path,
    rows: list[dict[str, str]],
    rendered: list[dict[str, str]] | None = None,
) -> None:
    before = read_tsv_row(before_report)
    before_row = {
        "label": "before",
        "method": "input",
        "motif_window_bp": "",
        "candidate_limit": "",
        "max_iterations": "",
        "poa_scoring": "",
        "gfa": str(args.input_gfa),
        "exit_status": "",
        "compare_exit_status": "",
        "path_preserving": "yes",
        "resolved": "",
        "bailed": "",
        "candidates_seen": "",
        "rounds": "",
        "generated_by_round": "",
        "applied_by_round": "",
        "stopped_reason": "",
        "wall_s": "",
        "max_rss_kb": "",
    }
    for field in GRAPH_REPORT_FIELDS:
        before_row[field] = before.get(field, "")

    fields = [
        "label",
        "method",
        "motif_window_bp",
        "candidate_limit",
        "max_iterations",
        "poa_scoring",
        "gfa",
        "exit_status",
        "compare_exit_status",
        "path_preserving",
        "resolved",
        "bailed",
        "candidates_seen",
        "rounds",
        "generated_by_round",
        "applied_by_round",
        "stopped_reason",
        "wall_s",
        "max_rss_kb",
        *GRAPH_REPORT_FIELDS,
        "compare_expected_paths",
        "compare_observed_paths",
        "compare_missing_paths",
        "compare_extra_paths",
        "compare_spelling_mismatches",
    ]
    write_tsv(dirs["reports"] / "matrix-summary.tsv", [before_row, *rows], fields)
    if rendered is not None:
        write_tsv(
            dirs["reports"] / "rendered-best.tsv",
            rendered,
            ["method", "label", "gfa", "sorted_gfa", "png", "upload_exit_status", "url"],
        )


if __name__ == "__main__":
    raise SystemExit(main())
