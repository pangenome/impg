#!/usr/bin/env python3
"""Run the C4 syng tail diagnosis matrix and summarize FASTA lengths."""

from __future__ import annotations

import argparse
import math
import os
import shlex
import statistics
import subprocess
import sys
import time
from pathlib import Path


DEFAULT_REGION = "GRCh38#0#chr6:31891045-32123783"
DEFAULT_SYNG = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng"
DEFAULT_AGC = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc"


VARIANTS: list[tuple[str, list[str]]] = [
    ("d50k", ["-d", "50k"]),
    ("d10k", ["-d", "10k"]),
    ("d1k", ["-d", "1k"]),
    ("d0", ["-d", "0"]),
    ("no_merge", ["--no-merge"]),
    ("d50k_frac02", ["-d", "50k", "--syng-min-chain-fraction", "0.02"]),
    ("d50k_frac05", ["-d", "50k", "--syng-min-chain-fraction", "0.05"]),
    ("d50k_frac10", ["-d", "50k", "--syng-min-chain-fraction", "0.10"]),
    ("d50k_frac20", ["-d", "50k", "--syng-min-chain-fraction", "0.20"]),
    ("d50k_frac50", ["-d", "50k", "--syng-min-chain-fraction", "0.50"]),
    ("d50k_frac80", ["-d", "50k", "--syng-min-chain-fraction", "0.80"]),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--binary", default="target/debug/impg")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--syng", default=DEFAULT_SYNG)
    parser.add_argument("--agc", default=DEFAULT_AGC)
    parser.add_argument("--region", default=DEFAULT_REGION)
    parser.add_argument("--threads", default="32")
    parser.add_argument("--variant", action="append", choices=[name for name, _ in VARIANTS])
    parser.add_argument("--skip-existing", action="store_true")
    return parser.parse_args()


def fasta_lengths(path: Path) -> list[tuple[str, int]]:
    records: list[tuple[str, int]] = []
    header: str | None = None
    length = 0
    with path.open("rt", encoding="utf-8") as handle:
        for raw in handle:
            line = raw.rstrip("\n")
            if line.startswith(">"):
                if header is not None:
                    records.append((header, length))
                header = line[1:]
                length = 0
            elif line:
                length += len(line)
    if header is not None:
        records.append((header, length))
    return records


def nearest_rank(sorted_lengths: list[int], fraction: float) -> int:
    if not sorted_lengths:
        return 0
    idx = max(0, math.ceil(fraction * len(sorted_lengths)) - 1)
    return sorted_lengths[min(idx, len(sorted_lengths) - 1)]


def summarize(records: list[tuple[str, int]]) -> dict[str, str]:
    lengths = sorted(length for _, length in records)
    if not lengths:
        return {
            "records": "0",
            "median": "0",
            "p95": "0",
            "p98": "0",
            "max": "0",
            "gt265k": "0",
            "gt285k": "0",
        }
    return {
        "records": str(len(lengths)),
        "median": str(round(statistics.median(lengths))),
        "p95": str(nearest_rank(lengths, 0.95)),
        "p98": str(nearest_rank(lengths, 0.98)),
        "max": str(lengths[-1]),
        "gt265k": str(sum(1 for length in lengths if length > 265_000)),
        "gt285k": str(sum(1 for length in lengths if length > 285_000)),
    }


def write_lengths(path: Path, records: list[tuple[str, int]]) -> None:
    with path.open("wt", encoding="utf-8") as handle:
        handle.write("record\tlength\n")
        for header, length in sorted(records, key=lambda item: (-item[1], item[0])):
            handle.write(f"{header}\t{length}\n")


def main() -> int:
    args = parse_args()
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    selected = dict(VARIANTS)
    names = args.variant or [name for name, _ in VARIANTS]
    commands_path = outdir / "commands.sh"
    metrics_path = outdir / "metrics.tsv"
    env = os.environ.copy()
    env["SYNG_PROFILE"] = "1"
    env["RAYON_NUM_THREADS"] = args.threads
    env["LD_LIBRARY_PATH"] = f"target/debug:{env.get('LD_LIBRARY_PATH', '')}".rstrip(":")

    with commands_path.open("wt", encoding="utf-8") as commands, metrics_path.open(
        "wt", encoding="utf-8"
    ) as metrics:
        metrics.write(
            "variant\tstatus\tseconds\trecords\tmedian\tp95\tp98\tmax\tgt265k\tgt285k\tfasta\n"
        )
        commands.write("#!/usr/bin/env bash\nset -euo pipefail\n")
        commands.write(f"export SYNG_PROFILE=1\nexport RAYON_NUM_THREADS={shlex.quote(args.threads)}\n")
        commands.write('export LD_LIBRARY_PATH="target/debug:${LD_LIBRARY_PATH:-}"\n\n')

        for name in names:
            variant_args = selected[name]
            prefix = outdir / name
            fasta = prefix.with_suffix(".fa")
            stdout_log = outdir / f"{name}.stdout.log"
            stderr_log = outdir / f"{name}.stderr.log"
            cmd = [
                args.binary,
                "query",
                "-a",
                args.syng,
                "-r",
                args.region,
                *variant_args,
                "--sequence-files",
                args.agc,
                "-o",
                "fasta",
                "--reverse-complement",
                "-O",
                str(prefix),
                "-t",
                args.threads,
                "-v",
                "1",
            ]
            quoted = " ".join(shlex.quote(part) for part in cmd)
            commands.write(f"{quoted} > {shlex.quote(str(stdout_log))} 2> {shlex.quote(str(stderr_log))}\n")
            commands.flush()

            if args.skip_existing and fasta.exists():
                status = 0
                seconds = 0.0
            else:
                print(f"[run] {name}: {quoted}", flush=True)
                start = time.monotonic()
                with stdout_log.open("wb") as stdout, stderr_log.open("wb") as stderr:
                    proc = subprocess.run(cmd, stdout=stdout, stderr=stderr, env=env, check=False)
                seconds = time.monotonic() - start
                status = proc.returncode

            if fasta.exists():
                records = fasta_lengths(fasta)
                write_lengths(outdir / f"{name}.lengths.tsv", records)
                summary = summarize(records)
            else:
                summary = summarize([])
            metrics.write(
                "\t".join(
                    [
                        name,
                        str(status),
                        f"{seconds:.2f}",
                        summary["records"],
                        summary["median"],
                        summary["p95"],
                        summary["p98"],
                        summary["max"],
                        summary["gt265k"],
                        summary["gt285k"],
                        str(fasta),
                    ]
                )
                + "\n"
            )
            metrics.flush()
            print(f"[done] {name}: status={status} metrics={summary}", flush=True)
            if status != 0:
                return status

    commands_path.chmod(0o755)
    return 0


if __name__ == "__main__":
    sys.exit(main())
