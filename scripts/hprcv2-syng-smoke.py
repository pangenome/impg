#!/usr/bin/env python3
"""
Run a syng query smoke panel against the local HPRCv2 index.

The panel is intentionally mixed: short windows, ordinary genic windows,
SV-rich immune/CNV loci, subtelomeric windows, and small satellite-proximal
windows. The goal is not truth-set validation; it is a broad "do not drop
everything, do not explode, do not run forever" smoke test.
"""

from __future__ import annotations

import argparse
import os
import re
import statistics
import subprocess
import sys
import time
from pathlib import Path


DEFAULT_PREFIX = "HPRC_r2_assemblies_0.6.1.syng"
DEFAULT_AGC = "HPRC_r2_assemblies_0.6.1.agc"


PANEL = [
    ("boring_chr1_100kb", "boring", "CHM13", "chr1", 100_000_000, 100_100_000),
    ("boring_chr6_100kb", "boring", "CHM13", "chr6", 50_000_000, 50_100_000),
    ("short_chr1_250bp", "short", "CHM13", "chr1", 100_000_000, 100_000_250),
    ("short_chr1_1kb", "short", "CHM13", "chr1", 100_000_000, 100_001_000),
    ("c4_grch38_53kb", "c4_rccx", "GRCh38", "chr6", 31_982_056, 32_035_418),
    ("c4_chm13_83kb", "c4_rccx", "CHM13", "chr6", 31_972_057, 32_055_418),
    ("c4_rccx_400kb", "c4_rccx", "CHM13", "chr6", 31_850_000, 32_250_000),
    ("mhc_class_i_100kb", "immune", "CHM13", "chr6", 30_000_000, 30_100_000),
    ("hla_dq_100kb", "immune", "CHM13", "chr6", 33_500_000, 33_600_000),
    ("lpa_100kb", "cnv_vntr", "CHM13", "chr6", 160_450_000, 160_550_000),
    ("amy_100kb", "cnv", "CHM13", "chr1", 103_500_000, 103_600_000),
    ("fcgr_100kb", "immune_cnv", "CHM13", "chr1", 161_450_000, 161_550_000),
    ("rhd_rhce_100kb", "immune_sv", "CHM13", "chr1", 25_200_000, 25_300_000),
    ("kir_100kb", "immune_cnv", "CHM13", "chr19", 54_500_000, 54_600_000),
    ("mapt_17q21_100kb", "inversion", "CHM13", "chr17", 45_000_000, 45_100_000),
    ("chr8p23_100kb", "inversion", "CHM13", "chr8", 12_000_000, 12_100_000),
    ("subtel_chr1_p_50kb", "subtelomere", "CHM13", "chr1", 10_000, 60_000),
    ("subtel_chr1_q_50kb", "subtelomere", "CHM13", "chr1", 248_300_000, 248_350_000),
    ("peri_cen_chr1_50kb", "pericentromere", "CHM13", "chr1", 130_000_000, 130_050_000),
    ("peri_cen_chr6_50kb", "pericentromere", "CHM13", "chr6", 60_000_000, 60_050_000),
    ("sat_chr8_10kb", "satellite", "CHM13", "chr8", 45_000_000, 45_010_000),
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hprcv2-dir", default="/home/erikg/hprcv2")
    parser.add_argument("--prefix", default=DEFAULT_PREFIX)
    parser.add_argument("--agc", default=DEFAULT_AGC)
    parser.add_argument("--impg", default=str(Path.home() / ".cargo/bin/impg"))
    parser.add_argument("--out-dir", default="/home/erikg/hprcv2/smoke/syng_query_panel")
    parser.add_argument("--timeout", type=int, default=3600)
    parser.add_argument("--merge-distance", type=int, default=10_000)
    parser.add_argument("--syng-extension", type=int, default=0)
    parser.add_argument("--syng-extend-budget", type=int, default=1_000)
    parser.add_argument("--syng-seed-walk-anchors", type=int, default=5)
    parser.add_argument("--syng-min-chain-fraction", type=float, default=0.5)
    parser.add_argument("--limit", type=int, default=0, help="Run only the first N panel rows")
    parser.add_argument("--raw", action="store_true", help="Use --syng-raw")
    return parser.parse_args()


def load_pansn_name_map(names_path: Path) -> dict[tuple[str, str], str]:
    mapping: dict[tuple[str, str], str] = {}
    with names_path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            full = parts[1]
            short = full.split(" ", 1)[0]
            pansn = short.split("#", 2)
            if len(pansn) == 3:
                sample, _hap, chrom = pansn
                mapping[(sample, chrom)] = full
    return mapping


def quantile(values: list[int], frac: float) -> int:
    if not values:
        return 0
    idx = int((len(values) - 1) * frac)
    return values[idx]


def parse_bed(output_bed: Path) -> dict[str, list[tuple[str, int]]]:
    by_label: dict[str, list[tuple[str, int]]] = {}
    with output_bed.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                by_label.setdefault("__malformed__", []).append((line.rstrip("\n"), 0))
                continue
            try:
                start = int(fields[1])
                end = int(fields[2])
            except ValueError:
                by_label.setdefault("__malformed__", []).append((line.rstrip("\n"), 0))
                continue
            by_label.setdefault(fields[3], []).append((fields[0], max(0, end - start)))
    return by_label


def parse_profiles(stderr_path: Path) -> dict[str, dict[str, str]]:
    profiles: dict[str, dict[str, str]] = {}
    current = None
    key_value = re.compile(r"([a-zA-Z_]+)=([^ )]+)")
    with stderr_path.open() as handle:
        for line in handle:
            if "Syng query:" in line:
                current = line.split("Syng query:", 1)[1].split(" (", 1)[0].strip()
                profiles.setdefault(current, {})
            elif current and ("EMITPROF" in line or "PROF chain_anchors" in line):
                bucket = profiles.setdefault(current, {})
                for key, value in key_value.findall(line):
                    bucket[key] = value
    return profiles


def main() -> int:
    args = parse_args()
    hprcv2 = Path(args.hprcv2_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    prefix = hprcv2 / args.prefix
    agc = hprcv2 / args.agc
    names_path = hprcv2 / f"{args.prefix}.names"
    if not names_path.exists():
        print(f"missing names file: {names_path}", file=sys.stderr)
        return 2

    name_map = load_pansn_name_map(names_path)
    rows = PANEL[: args.limit] if args.limit else PANEL
    panel_bed = out_dir / "panel.bed"
    panel_meta = out_dir / "panel.meta.tsv"
    with panel_bed.open("w") as bed, panel_meta.open("w") as meta:
        meta.write("label\tcategory\tsample\tchrom\tstart\tend\tquery_len\tfull_name\n")
        for label, category, sample, chrom, start, end in rows:
            full = name_map.get((sample, chrom))
            if full is None:
                raise RuntimeError(f"No syng name found for {sample} {chrom}")
            bed.write(f"{full}\t{start}\t{end}\t{label}\n")
            meta.write(
                f"{label}\t{category}\t{sample}\t{chrom}\t{start}\t{end}\t{end-start}\t{full}\n"
            )

    output_bed = out_dir / ("panel.raw.bed" if args.raw else "panel.refined.bed")
    stderr_log = out_dir / ("panel.raw.stderr" if args.raw else "panel.refined.stderr")
    summary_tsv = out_dir / ("panel.raw.summary.tsv" if args.raw else "panel.refined.summary.tsv")

    cmd = [
        args.impg,
        "query",
        "-a",
        str(prefix),
        "-b",
        str(panel_bed),
        "--sequence-files",
        str(agc),
        "-o",
        "bed",
        "-d",
        str(args.merge_distance),
        "--syng-extension",
        str(args.syng_extension),
        "--syng-extend-budget",
        str(args.syng_extend_budget),
        "--syng-seed-walk-anchors",
        str(args.syng_seed_walk_anchors),
        "--syng-min-chain-fraction",
        str(args.syng_min_chain_fraction),
    ]
    if args.raw:
        cmd.append("--syng-raw")

    env = os.environ.copy()
    env["SYNG_PROFILE"] = "1"
    env["SYNG_EMIT_PROFILE"] = "1"
    started = time.time()
    with output_bed.open("w") as stdout, stderr_log.open("w") as stderr:
        proc = subprocess.run(
            cmd,
            cwd=hprcv2,
            env=env,
            stdout=stdout,
            stderr=stderr,
            timeout=args.timeout,
            check=False,
            text=True,
        )
    elapsed = time.time() - started
    if proc.returncode != 0:
        print(f"query failed with exit code {proc.returncode}; see {stderr_log}", file=sys.stderr)
        return proc.returncode

    by_label = parse_bed(output_bed)
    profiles = parse_profiles(stderr_log)
    stdout_has_diag = "read GBWT" in output_bed.read_text(errors="replace")
    with summary_tsv.open("w") as out:
        out.write(
            "label\tcategory\tquery_len\tn\tunique_targets\tmin\tp05\tmedian\tp95\tmax\tmean"
            "\tseeds\tlocated_hits\tchains\tafter_filter\tspan_cap\tstdout_has_diag\n"
        )
        for label, category, _sample, _chrom, start, end in rows:
            hits = by_label.get(label, [])
            lengths = sorted(length for _target, length in hits)
            targets = {target for target, _length in hits}
            prof = profiles.get(label, {})
            mean = f"{statistics.fmean(lengths):.2f}" if lengths else "0"
            out.write(
                f"{label}\t{category}\t{end-start}\t{len(lengths)}\t{len(targets)}"
                f"\t{lengths[0] if lengths else 0}"
                f"\t{quantile(lengths, 0.05)}"
                f"\t{quantile(lengths, 0.50)}"
                f"\t{quantile(lengths, 0.95)}"
                f"\t{lengths[-1] if lengths else 0}"
                f"\t{mean}"
                f"\t{prof.get('seeds', '')}"
                f"\t{prof.get('located_hits', '')}"
                f"\t{prof.get('chains', '')}"
                f"\t{prof.get('after_filter', '')}"
                f"\t{prof.get('span_cap', '')}"
                f"\t{int(stdout_has_diag)}\n"
            )

    print(f"elapsed_seconds\t{elapsed:.2f}")
    print(f"panel_bed\t{panel_bed}")
    print(f"output_bed\t{output_bed}")
    print(f"stderr_log\t{stderr_log}")
    print(f"summary_tsv\t{summary_tsv}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
