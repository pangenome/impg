#!/usr/bin/env python3
"""Run the C4 graph-quality optimization scoreboard.

The driver writes heavyweight C4 artifacts under a data directory and keeps
the committed surface to this script plus small TSV/Markdown summaries.  Each
trial records the command, runtime, exact path validation, graph-report TSV
metrics, PAF yields when the trial performs query/graph construction, and an
explicit recheck of the ten diagnosed left-edge paths.
"""

from __future__ import annotations

import argparse
import csv
import os
import re
import shlex
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Iterable


INDEX = Path("/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng")
AGC = Path("/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc")
CONTROL_ROOT = Path("/home/erikg/impg/data/c4_control_provenance_20260610T154843Z")
CONTROL_GFA = CONTROL_ROOT / "graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa"
CONTROL_SORTED = CONTROL_ROOT / "sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa"
CONTROL_SOURCE_GFA = CONTROL_ROOT / "validation/current_53kb_pggb.source.gfa"
CONTROL_PATH_SPELLINGS = CONTROL_ROOT / "debug/current_53kb_pggb_query/path_spellings.tsv"
BED = CONTROL_ROOT / "c4_53kb.bed"
DEFAULT_IMPG = Path("/home/erikg/impg/target/release/impg")
DEFAULT_COMPARE = Path("/home/erikg/impg/target/release/examples/compare_gfa_paths")
PUBLIC_BASE = "http://hypervolu.me/~erik/impg/"
UPLOAD_TARGET = "erik@hypervolu.me:www/impg/"
REGION_BASENAME = "C4_GRCh38_53kb"

SUSPECT_LEFT_PATHS = [
    "HG01952#2#CM087999.1:31995790-32075638",
    "HG01960#2#JBHIHN010000011.1:32078636-32158486",
    "HG02027#2#JBHDTY010000008.1:5188422-5268271",
    "HG02132#2#CM086908.1:31997546-32077395",
    "HG02178#2#CM089938.1:31943310-32055920",
    "HG03710#2#CM086882.1:31974825-32061063",
    "HG03874#1#CM089432.1:31915349-32027959",
    "NA18565#2#CM094342.1:31929266-32009115",
    "NA18945#2#CM101630.1:32018016-32097865",
    "NA18960#2#CM101559.1:31867170-31947019",
]

GRAPH_METRICS = [
    "paths",
    "segments",
    "links",
    "path_steps",
    "total_segment_bp",
    "node_coverage_bp_weighted_mean",
    "bp_weighted_cov",
    "singleton_bp",
    "sparse_coverage_segment_bp",
    "segment_white_space_bp_fraction",
    "segment_white_space_bp_total",
    "path_white_space_bp_p99",
    "path_white_space_bp_max",
    "path_white_space_bridges_ge_threshold",
    "direct_self_loop_edges",
    "adjacent_same_node_path_steps",
    "povu_sites",
    "povu_leaf_sites",
]

SCOREBOARD_FIELDS = [
    "trial",
    "variant_kind",
    "method_settings",
    "status",
    "exact_path_validation",
    "expected_paths",
    "observed_paths",
    "missing_paths",
    "extra_paths",
    "spelling_mismatches",
    "runtime_wall",
    "runtime_max_rss_kb",
    "runtime_exit_status",
    "raw_paf_records",
    "filtered_paf_records",
    "final_paf_records",
    "resolved_replacements",
    "bailed_replacements",
    "candidates_seen",
    "rounds",
    *GRAPH_METRICS,
    "left_suspect_fixed_count",
    "left_suspect_extended_count",
    "left_suspect_still_late_count",
    "left_suspect_missing_count",
    "gfa",
    "sorted_gfa",
    "source_gfa",
    "graph_report_tsv",
    "compare_stdout",
    "compare_stderr",
    "png_url",
    "command",
]

LEFT_FIELDS = [
    "trial",
    "suspect_path",
    "suspect_key",
    "current_start",
    "current_end",
    "current_length",
    "trial_path",
    "trial_start",
    "trial_end",
    "trial_length",
    "start_delta_vs_current",
    "first_node",
    "first_node_graph_bp",
    "left_side_status",
]


@dataclass(frozen=True)
class Trial:
    name: str
    kind: str
    settings: str
    command: tuple[str, ...]
    source_gfa: Path | None
    debug_dir: Path | None


@dataclass
class Runtime:
    wall: str = ""
    rss_kb: str = ""
    exit_status: str = ""
    command: str = ""
    stdout: Path | None = None
    stderr: Path | None = None
    time_log: Path | None = None


def parse_args() -> argparse.Namespace:
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=Path(f"/home/erikg/impg/data/c4_graph_quality_{stamp}"))
    parser.add_argument("--impg", type=Path, default=DEFAULT_IMPG)
    parser.add_argument("--compare-gfa-paths", type=Path, default=DEFAULT_COMPARE)
    parser.add_argument("--threads", type=int, default=16)
    parser.add_argument("--timeout", type=int, default=3600, help="Per heavy command timeout in seconds")
    parser.add_argument("--trials", help="Comma-separated trial names; default runs all non-control trials")
    parser.add_argument(
        "--preset",
        choices=["graph-quality", "repeat-filter-sweep"],
        default="graph-quality",
        help="Named trial set to run when --trials is not supplied",
    )
    parser.add_argument("--include-control", action="store_true", help="Also score the current impg PGGB-style control row")
    parser.add_argument("--force", action="store_true")
    parser.add_argument("--render-labels", help="Comma-separated trial names to render with gfalook -m and upload")
    parser.add_argument("--no-upload", action="store_true")
    parser.add_argument("--upload-target", default=UPLOAD_TARGET)
    parser.add_argument("--public-base", default=PUBLIC_BASE)
    return parser.parse_args()


def ensure_dirs(root: Path) -> None:
    for name in ["debug", "graphs", "logs", "reports", "sorted", "validation", "renders"]:
        (root / name).mkdir(parents=True, exist_ok=True)


def shell_join(command: Iterable[object]) -> str:
    return " ".join(shlex.quote(str(part)) for part in command)


def append_command(out_dir: Path, label: str, command: Iterable[object]) -> None:
    path = out_dir / "commands.md"
    with path.open("a", encoding="utf-8") as handle:
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
    out_dir: Path,
    label: str,
    command: list[object],
    *,
    timeout: int = 0,
    cwd: Path | None = None,
    force: bool = False,
    output_exists: Path | None = None,
) -> Runtime:
    stdout = out_dir / "logs" / f"{label}.stdout.log"
    stderr = out_dir / "logs" / f"{label}.stderr.log"
    time_log = out_dir / "logs" / f"{label}.time.txt"
    if output_exists is not None and output_exists.exists() and not force:
        wall, rss, status = parse_time_log(time_log)
        return Runtime(wall, rss, status or "0", shell_join(command), stdout, stderr, time_log)

    full_command: list[object] = ["/usr/bin/time", "-v", "-o", time_log]
    if timeout > 0:
        full_command.extend(["/usr/bin/timeout", f"{timeout}s"])
    full_command.extend(command)
    append_command(out_dir, label, full_command)
    print(f"[{datetime.now(timezone.utc).isoformat()}] {label}", flush=True)
    with stdout.open("wb") as out, stderr.open("wb") as err:
        proc = subprocess.run([str(part) for part in full_command], cwd=cwd, stdout=out, stderr=err)
    wall, rss, status = parse_time_log(time_log)
    return Runtime(wall, rss, status or str(proc.returncode), shell_join(full_command), stdout, stderr, time_log)


def read_tsv_row(path: Path) -> dict[str, str]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            return {key: value or "" for key, value in row.items()}
    return {}


def parse_compare(path: Path) -> dict[str, str]:
    out: dict[str, str] = {}
    if not path.exists():
        return out
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if "\t" in line:
            key, value = line.split("\t", 1)
            out[key] = value
    return out


def count_records(path: Path | None) -> str:
    if path is None or not path.exists():
        return ""
    count = 0
    with path.open("rb") as handle:
        for _ in handle:
            count += 1
    return str(count)


def parse_crush_summary(stderr: Path | None) -> dict[str, str]:
    if stderr is None or not stderr.exists():
        return {}
    text = stderr.read_text(encoding="utf-8", errors="replace")
    match = re.search(r"crush: (\d+) resolved, (\d+) bailed, (\d+) candidates seen across (\d+) rounds", text)
    if not match:
        return {}
    return {
        "resolved_replacements": match.group(1),
        "bailed_replacements": match.group(2),
        "candidates_seen": match.group(3),
        "rounds": match.group(4),
    }


def read_fasta(path: Path) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    name: str | None = None
    chunks: list[str] = []
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(chunks)))
                name = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        records.append((name, "".join(chunks)))
    return records


def source_gfa_from_fasta(fasta: Path, output: Path) -> None:
    records = read_fasta(fasta)
    with output.open("w", encoding="utf-8") as handle:
        handle.write("H\tVN:Z:1.0\n")
        for idx, (name, seq) in enumerate(records):
            sid = f"source_{idx}"
            handle.write(f"S\t{sid}\t{seq}\n")
            handle.write(f"P\t{name}\t{sid}+\t*\n")


def gfa_path_first_offsets(gfa: Path) -> dict[str, tuple[str, int]]:
    seg_start: dict[str, int] = {}
    offset = 0
    first: dict[str, tuple[str, int]] = {}
    with gfa.open("r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("S\t"):
                parts = line.rstrip("\n").split("\t")
                seg_start[parts[1]] = offset
                offset += len(parts[2])
            elif line.startswith("P\t"):
                parts = line.rstrip("\n").split("\t")
                steps = [step[:-1] for step in parts[2].split(",") if step]
                if steps:
                    first_node = steps[0]
                    first[parts[1]] = (first_node, seg_start.get(first_node, -1))
    return first


def parse_interval_name(path_name: str) -> tuple[str, int | None, int | None]:
    if ":" not in path_name:
        return path_name, None, None
    key, coords = path_name.rsplit(":", 1)
    match = re.fullmatch(r"(\d+)-(\d+)", coords)
    if not match:
        return path_name, None, None
    return key, int(match.group(1)), int(match.group(2))


def path_lengths(path_spellings: Path | None, source_gfa: Path | None = None) -> dict[str, int]:
    if path_spellings is not None and path_spellings.exists():
        out: dict[str, int] = {}
        with path_spellings.open("r", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            for row in reader:
                try:
                    out[row["sequence"]] = int(row["length"])
                except (KeyError, ValueError):
                    continue
        return out
    if source_gfa is not None and source_gfa.exists():
        seqs: dict[str, int] = {}
        current_path: str | None = None
        seg_len: dict[str, int] = {}
        with source_gfa.open("r", encoding="utf-8", errors="replace") as handle:
            for line in handle:
                if line.startswith("S\t"):
                    _, sid, seq, *_ = line.rstrip("\n").split("\t")
                    seg_len[sid] = len(seq)
                elif line.startswith("P\t"):
                    parts = line.rstrip("\n").split("\t")
                    current_path = parts[1]
                    steps = [step[:-1] for step in parts[2].split(",") if step]
                    seqs[current_path] = sum(seg_len.get(step, 0) for step in steps)
        return seqs
    return {}


def left_recheck(
    trial: str,
    sorted_gfa: Path,
    path_spellings_file: Path | None,
    source_gfa: Path | None,
) -> tuple[list[dict[str, str]], dict[str, str]]:
    current_lengths = path_lengths(CONTROL_PATH_SPELLINGS, CONTROL_SOURCE_GFA)
    trial_lengths = path_lengths(path_spellings_file, source_gfa)
    first_offsets = gfa_path_first_offsets(sorted_gfa) if sorted_gfa.exists() else {}

    by_key: dict[str, str] = {}
    for path_name in trial_lengths:
        key, _, _ = parse_interval_name(path_name)
        by_key.setdefault(key, path_name)

    rows: list[dict[str, str]] = []
    counts = {
        "left_suspect_fixed_count": 0,
        "left_suspect_extended_count": 0,
        "left_suspect_still_late_count": 0,
        "left_suspect_missing_count": 0,
    }

    for suspect in SUSPECT_LEFT_PATHS:
        key, current_start, current_end = parse_interval_name(suspect)
        trial_path = suspect if suspect in trial_lengths else by_key.get(key, "")
        trial_key, trial_start, trial_end = parse_interval_name(trial_path) if trial_path else (key, None, None)
        first_node, offset = first_offsets.get(trial_path, ("", -1))
        current_len = current_lengths.get(suspect, 0)
        trial_len = trial_lengths.get(trial_path, 0) if trial_path else 0
        delta = ""
        if current_start is not None and trial_start is not None:
            delta = str(trial_start - current_start)
        if not trial_path:
            status = "missing"
            counts["left_suspect_missing_count"] += 1
        elif offset >= 0 and offset < 1000:
            status = "left_side_present_in_graph_layout"
            counts["left_suspect_fixed_count"] += 1
        elif current_start is not None and trial_start is not None and trial_start < current_start:
            status = "source_interval_extended_left_but_layout_still_late"
            counts["left_suspect_extended_count"] += 1
            counts["left_suspect_still_late_count"] += 1
        else:
            status = "still_late_first_node"
            counts["left_suspect_still_late_count"] += 1
        rows.append(
            {
                "trial": trial,
                "suspect_path": suspect,
                "suspect_key": key,
                "current_start": str(current_start or ""),
                "current_end": str(current_end or ""),
                "current_length": str(current_len or ""),
                "trial_path": trial_path,
                "trial_start": str(trial_start or ""),
                "trial_end": str(trial_end or ""),
                "trial_length": str(trial_len or ""),
                "start_delta_vs_current": delta,
                "first_node": first_node,
                "first_node_graph_bp": str(offset if offset >= 0 else ""),
                "left_side_status": status,
            }
        )
    return rows, {key: str(value) for key, value in counts.items()}


def selected_trials(args: argparse.Namespace) -> list[str]:
    if args.trials:
        return [item.strip() for item in args.trials.split(",") if item.strip()]
    if args.preset == "repeat-filter-sweep":
        return [
            "q_seed_sensitive_pggb",
            "q_seed_nm11_sf11_pggb",
            "q_seed_nm1many_sf11_pggb",
            "q_seed_nm11_sf11_nojump_pggb",
            "q_seed_minid98_pggb",
            "q_seed_minid99_pggb",
        ]
    return [
        "q_seed_sensitive_pggb",
        "q_boundary_wide_pggb",
        "crush_coverage_repeat_aware",
        "crush_motif_local_right",
        "crush_outward_guarded",
        "crush_chain_povu_smooth",
    ]


def render_labels(args: argparse.Namespace) -> set[str]:
    if not args.render_labels:
        return set()
    return {item.strip() for item in args.render_labels.split(",") if item.strip()}


def query_trial(args: argparse.Namespace, name: str, extra: list[str], settings: str) -> Trial:
    debug = args.out_dir / "debug" / f"{name}_query"
    out_prefix = args.out_dir / "graphs" / name
    command = [
        str(args.impg),
        "query",
        "-a",
        str(INDEX),
        "-b",
        str(BED),
        "-d",
        "100k",
        "--sequence-files",
        str(AGC),
        "--debug-dir",
        str(debug),
        *extra,
        "-o",
        "gfa:pggb",
        "-O",
        str(out_prefix),
        "-t",
        str(args.threads),
        "-v",
        "1",
    ]
    return Trial(name, "query-pggb", settings, tuple(command), None, debug)


def crush_trial(args: argparse.Namespace, name: str, extra: list[str], settings: str) -> Trial:
    out_gfa = args.out_dir / "graphs" / f"{name}.gfa"
    command = [
        str(args.impg),
        "crush",
        "--gfa",
        str(CONTROL_GFA),
        "--output",
        str(out_gfa),
        "-t",
        str(args.threads),
        "-v",
        "1",
        *extra,
    ]
    return Trial(name, "post-control-crush", settings, tuple(command), CONTROL_SOURCE_GFA, None)


def trial_catalog(args: argparse.Namespace) -> dict[str, Trial]:
    return {
        "control_current_pggb": Trial(
            "control_current_pggb",
            "control",
            "current impg PGGB-style control; no construction rerun",
            tuple(),
            CONTROL_SOURCE_GFA,
            CONTROL_ROOT / "debug/current_53kb_pggb_query",
        ),
        "q_seed_sensitive_pggb": query_trial(
            args,
            "q_seed_sensitive_pggb",
            ["--syng-seed-drop-top-fraction", "0", "--syng-seed-walk-anchors", "3"],
            "pggb query with high-copy seed drop disabled and 3-walk seeds for left-edge sensitivity",
        ),
        "q_seed_nm11_sf11_pggb": query_trial(
            args,
            "q_seed_nm11_sf11_pggb",
            [
                "--syng-seed-drop-top-fraction",
                "0",
                "--syng-seed-walk-anchors",
                "3",
                "--num-mappings",
                "1:1",
                "--scaffold-filter",
                "1:1",
            ],
            "seed-sensitive pggb plus strict 1:1 mapping-axis and 1:1 scaffold-chain filtering",
        ),
        "q_seed_nm1many_sf11_pggb": query_trial(
            args,
            "q_seed_nm1many_sf11_pggb",
            [
                "--syng-seed-drop-top-fraction",
                "0",
                "--syng-seed-walk-anchors",
                "3",
                "--num-mappings",
                "1:many",
                "--scaffold-filter",
                "1:1",
            ],
            "seed-sensitive pggb with one best mapping per query, many per target, and strict 1:1 scaffold filtering",
        ),
        "q_seed_nm11_sf11_nojump_pggb": query_trial(
            args,
            "q_seed_nm11_sf11_nojump_pggb",
            [
                "--syng-seed-drop-top-fraction",
                "0",
                "--syng-seed-walk-anchors",
                "3",
                "--num-mappings",
                "1:1",
                "--scaffold-filter",
                "1:1",
                "--scaffold-jump",
                "0",
            ],
            "seed-sensitive pggb with 1:1 mapping/scaffold filtering and scaffold chaining disabled",
        ),
        "q_seed_minid98_pggb": query_trial(
            args,
            "q_seed_minid98_pggb",
            [
                "--syng-seed-drop-top-fraction",
                "0",
                "--syng-seed-walk-anchors",
                "3",
                "--min-aln-identity",
                "0.98",
            ],
            "seed-sensitive pggb with first-class post-alignment PAF identity filter at 98%",
        ),
        "q_seed_minid99_pggb": query_trial(
            args,
            "q_seed_minid99_pggb",
            [
                "--syng-seed-drop-top-fraction",
                "0",
                "--syng-seed-walk-anchors",
                "3",
                "--min-aln-identity",
                "0.99",
            ],
            "seed-sensitive pggb with first-class post-alignment PAF identity filter at 99%",
        ),
        "q_boundary_wide_pggb": query_trial(
            args,
            "q_boundary_wide_pggb",
            ["--syng-padding", "1000", "--syng-extension", "5000", "--syng-extend-budget", "5000"],
            "pggb query with 1kb syng padding and 5kb source-side/boundary extension",
        ),
        "crush_coverage_repeat_aware": crush_trial(
            args,
            "crush_coverage_repeat_aware",
            [
                "--method",
                "coverage-multi-bubble",
                "--objective",
                "coverage",
                "--repeat-aware-boundaries=true",
                "--window-mode",
                "combined",
                "--max-window-sites",
                "12",
                "--max-iterations",
                "2",
                "--candidate-limit",
                "192",
                "--min-match-length",
                "off",
                "--max-pair-alignments",
                "0",
                "--max-replacement-paf-bytes",
                "0",
                "--sweepga-no-filter=true",
                "--polish-rounds",
                "0",
            ],
            "occurrence-aware coverage multi-bubble windows with repeat-aware boundaries",
        ),
        "crush_motif_local_right": crush_trial(
            args,
            "crush_motif_local_right",
            [
                "--method",
                "motif-local",
                "--motif-max-sparse-paths",
                "5",
                "--motif-min-flank-paths",
                "12",
                "--motif-min-order-jump",
                "1k",
                "--motif-max-window-bp",
                "10k",
                "--max-iterations",
                "2",
                "--polish-rounds",
                "0",
            ],
            "motif-local right-repeat offshoot windows with explicit flank support",
        ),
        "crush_outward_multibubble": crush_trial(
            args,
            "crush_outward_multibubble",
            [
                "--method",
                "iterative-multi-level",
                "--window-mode",
                "outward",
                "--max-window-sites",
                "10",
                "--window-target-bp",
                "30k",
                "--candidate-limit",
                "192",
                "--max-iterations",
                "2",
                "--min-match-length",
                "off",
                "--max-pair-alignments",
                "0",
                "--max-replacement-paf-bytes",
                "0",
                "--sweepga-no-filter=true",
                "--polish-rounds",
                "0",
            ],
            "outward occurrence-window iterative multi-bubble smoothing for right condensation",
        ),
        "crush_outward_guarded": crush_trial(
            args,
            "crush_outward_guarded",
            [
                "--method",
                "iterative-multi-level",
                "--window-mode",
                "outward",
                "--max-window-sites",
                "6",
                "--window-target-bp",
                "8k",
                "--candidate-limit",
                "32",
                "--max-iterations",
                "1",
                "--max-pair-alignments",
                "10000",
                "--max-replacement-paf-bytes",
                "67108864",
                "--polish-rounds",
                "0",
            ],
            "guarded outward occurrence windows with bounded candidate count and replacement PAF size",
        ),
        "crush_chain_povu_smooth": crush_trial(
            args,
            "crush_chain_povu_smooth",
            [
                "--method",
                "chain-povu",
                "--chain-target-bp",
                "10k",
                "--polish-method",
                "smooth",
                "--polish-rounds",
                "1",
                "--max-iterations",
                "1",
            ],
            "chain-POVU smoothxg-to-POASTA occurrence block construction",
        ),
    }


def trial_gfa_path(args: argparse.Namespace, trial: Trial) -> Path:
    if trial.name == "control_current_pggb":
        return CONTROL_GFA
    if trial.kind == "query-pggb":
        return args.out_dir / "graphs" / trial.name / f"{REGION_BASENAME}.gfa"
    return args.out_dir / "graphs" / f"{trial.name}.gfa"


def run_sort(args: argparse.Namespace, trial: str, gfa: Path) -> tuple[Path, Runtime]:
    if trial == "control_current_pggb":
        return CONTROL_SORTED, Runtime(command="pre-existing control sorted GFA", exit_status="0")
    sorted_gfa = args.out_dir / "sorted" / f"{trial}.Ygs.gfa"
    runtime = run_timed(
        args.out_dir,
        f"{trial}.gfasort",
        ["gfasort", "-i", gfa, "-o", sorted_gfa, "-p", "Ygs", "-t", str(args.threads), "-v", "1"],
        timeout=args.timeout,
        force=args.force,
        output_exists=sorted_gfa,
    )
    return sorted_gfa, runtime


def run_graph_report(args: argparse.Namespace, trial: str, sorted_gfa: Path) -> tuple[Path, dict[str, str], Runtime]:
    report = args.out_dir / "reports" / f"{trial}.graph-report.tsv"
    runtime = run_timed(
        args.out_dir,
        f"{trial}.graph-report",
        [
            str(args.impg),
            "graph-report",
            "-g",
            sorted_gfa,
            "-o",
            report,
            "--format",
            "tsv",
            "--povu",
            "--top",
            "20",
            "-t",
            str(args.threads),
            "-v",
            "1",
        ],
        timeout=args.timeout,
        force=args.force,
        output_exists=report,
    )
    return report, read_tsv_row(report), runtime


def run_compare(args: argparse.Namespace, trial: str, source_gfa: Path, sorted_gfa: Path) -> tuple[dict[str, str], Runtime]:
    stdout = args.out_dir / "validation" / f"{trial}.compare_gfa_paths.stdout.log"
    runtime = run_timed(
        args.out_dir,
        f"{trial}.compare_gfa_paths",
        [args.compare_gfa_paths, source_gfa, sorted_gfa],
        timeout=args.timeout,
        force=args.force,
        output_exists=stdout,
    )
    # run_timed's generic stdout path is under logs; copy-free parse from there.
    return parse_compare(runtime.stdout or stdout), runtime


def render_and_upload(args: argparse.Namespace, trial: str, sorted_gfa: Path) -> str:
    png = args.out_dir / "renders" / f"{trial}.Ygs.gfalook-m.png"
    render = run_timed(
        args.out_dir,
        f"{trial}.gfalook",
        ["gfalook", "-i", sorted_gfa, "-o", png, "-m", "-x", "3200", "-y", "1800", "-a", "3", "-t", str(args.threads), "-v", "1"],
        timeout=args.timeout,
        force=args.force,
        output_exists=png,
    )
    if render.exit_status not in ("0", "") or not png.exists():
        return ""
    if args.no_upload:
        return ""
    remote_name = f"c4-graph-quality-{trial}-{args.out_dir.name}.Ygs.gfalook-m.png"
    remote = f"{args.upload_target.rstrip('/')}/{remote_name}"
    upload = run_timed(
        args.out_dir,
        f"{trial}.upload",
        ["scp", "-o", "BatchMode=yes", "-o", "ConnectTimeout=20", png, remote],
        timeout=120,
        force=args.force,
    )
    if upload.exit_status not in ("0", ""):
        return ""
    return f"{args.public_base.rstrip('/')}/{remote_name}"


def write_rows(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fields})


def run_trial(args: argparse.Namespace, trial: Trial, should_render: bool) -> tuple[dict[str, str], list[dict[str, str]]]:
    gfa = trial_gfa_path(args, trial)
    runtime = Runtime(command="pre-existing control graph", exit_status="0")
    status = "ok"
    source_gfa = trial.source_gfa
    path_spellings_file = trial.debug_dir / "path_spellings.tsv" if trial.debug_dir else None

    if trial.command:
        output = gfa
        runtime = run_timed(
            args.out_dir,
            f"{trial.name}.{trial.kind}",
            list(trial.command),
            timeout=args.timeout,
            force=args.force,
            output_exists=output,
        )
        if runtime.exit_status not in ("0", ""):
            status = "command_failed"
    if trial.kind == "query-pggb" and trial.debug_dir is not None:
        source_gfa = args.out_dir / "validation" / f"{trial.name}.source.gfa"
        fasta = trial.debug_dir / "combined.fa"
        if fasta.exists() and (args.force or not source_gfa.exists()):
            source_gfa_from_fasta(fasta, source_gfa)

    if not gfa.exists():
        row = {
            "trial": trial.name,
            "variant_kind": trial.kind,
            "method_settings": trial.settings,
            "status": status if status != "ok" else "missing_gfa",
            "runtime_wall": runtime.wall,
            "runtime_max_rss_kb": runtime.rss_kb,
            "runtime_exit_status": runtime.exit_status,
            "command": runtime.command,
        }
        return row, []

    sorted_gfa, _sort_runtime = run_sort(args, trial.name, gfa)
    report_path, report, _report_runtime = run_graph_report(args, trial.name, sorted_gfa)

    compare: dict[str, str] = {}
    compare_runtime = Runtime()
    if source_gfa is not None and source_gfa.exists():
        compare, compare_runtime = run_compare(args, trial.name, source_gfa, sorted_gfa)
    exact = "pass" if compare.get("spelling_mismatches") == "0" and compare.get("missing_paths") == "0" and compare.get("extra_paths") == "0" else "fail"

    left_rows, left_counts = left_recheck(trial.name, sorted_gfa, path_spellings_file, source_gfa)
    png_url = render_and_upload(args, trial.name, sorted_gfa) if should_render else ""
    crush_summary = parse_crush_summary(runtime.stderr)
    debug = trial.debug_dir
    metrics = {metric: report.get(metric, "") for metric in GRAPH_METRICS}
    metrics["bp_weighted_cov"] = report.get("node_coverage_bp_weighted_mean", metrics.get("bp_weighted_cov", ""))

    row = {
        "trial": trial.name,
        "variant_kind": trial.kind,
        "method_settings": trial.settings,
        "status": status,
        "exact_path_validation": exact,
        "expected_paths": compare.get("expected_paths", ""),
        "observed_paths": compare.get("observed_paths", ""),
        "missing_paths": compare.get("missing_paths", ""),
        "extra_paths": compare.get("extra_paths", ""),
        "spelling_mismatches": compare.get("spelling_mismatches", ""),
        "runtime_wall": runtime.wall,
        "runtime_max_rss_kb": runtime.rss_kb,
        "runtime_exit_status": runtime.exit_status,
        "raw_paf_records": count_records(debug / "raw.paf") if debug else "",
        "filtered_paf_records": count_records(debug / "filtered.paf") if debug else "",
        "final_paf_records": count_records(debug / "final.paf") if debug else "",
        "resolved_replacements": crush_summary.get("resolved_replacements", ""),
        "bailed_replacements": crush_summary.get("bailed_replacements", ""),
        "candidates_seen": crush_summary.get("candidates_seen", ""),
        "rounds": crush_summary.get("rounds", ""),
        **metrics,
        **left_counts,
        "gfa": str(gfa),
        "sorted_gfa": str(sorted_gfa),
        "source_gfa": str(source_gfa or ""),
        "graph_report_tsv": str(report_path),
        "compare_stdout": str(compare_runtime.stdout or ""),
        "compare_stderr": str(compare_runtime.stderr or ""),
        "png_url": png_url,
        "command": runtime.command,
    }
    return row, left_rows


def write_markdown_summary(out_dir: Path, scoreboard: list[dict[str, str]], left_rows: list[dict[str, str]]) -> None:
    md = out_dir / "summary.md"
    cols = [
        "trial",
        "exact_path_validation",
        "segments",
        "total_segment_bp",
        "singleton_bp",
        "sparse_coverage_segment_bp",
        "segment_white_space_bp_fraction",
        "path_white_space_bp_p99",
        "path_white_space_bp_max",
        "path_white_space_bridges_ge_threshold",
        "povu_sites",
        "left_suspect_fixed_count",
        "left_suspect_extended_count",
        "left_suspect_still_late_count",
        "png_url",
    ]
    with md.open("w", encoding="utf-8") as handle:
        handle.write("# C4 Graph-Quality Scoreboard Summary\n\n")
        handle.write(f"Run root: `{out_dir}`\n\n")
        handle.write("| " + " | ".join(cols) + " |\n")
        handle.write("|" + "|".join(["---"] * len(cols)) + "|\n")
        for row in scoreboard:
            handle.write("| " + " | ".join(row.get(col, "") for col in cols) + " |\n")
        handle.write("\n## Left-Edge Recheck\n\n")
        handle.write("The full machine-readable table is `left_edge_recheck.tsv`.\n\n")
        handle.write("| trial | suspect_path | trial_path | first_node_graph_bp | left_side_status |\n")
        handle.write("|---|---|---|---:|---|\n")
        for row in left_rows:
            handle.write(
                "| "
                + " | ".join(
                    [
                        row.get("trial", ""),
                        row.get("suspect_path", ""),
                        row.get("trial_path", ""),
                        row.get("first_node_graph_bp", ""),
                        row.get("left_side_status", ""),
                    ]
                )
                + " |\n"
            )


def main() -> int:
    args = parse_args()
    ensure_dirs(args.out_dir)
    if not args.impg.exists():
        raise FileNotFoundError(args.impg)
    if not args.compare_gfa_paths.exists():
        raise FileNotFoundError(args.compare_gfa_paths)

    wanted = selected_trials(args)
    if args.include_control:
        wanted = ["control_current_pggb", *wanted]
    catalog = trial_catalog(args)
    missing = [name for name in wanted if name not in catalog]
    if missing:
        raise SystemExit(f"unknown trial(s): {', '.join(missing)}")

    render = render_labels(args)
    scoreboard: list[dict[str, str]] = []
    left: list[dict[str, str]] = []
    for name in wanted:
        row, left_rows = run_trial(args, catalog[name], name in render)
        scoreboard.append(row)
        left.extend(left_rows)
        write_rows(args.out_dir / "scoreboard.tsv", scoreboard, SCOREBOARD_FIELDS)
        write_rows(args.out_dir / "left_edge_recheck.tsv", left, LEFT_FIELDS)
        write_markdown_summary(args.out_dir, scoreboard, left)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
