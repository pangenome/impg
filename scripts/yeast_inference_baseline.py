#!/usr/bin/env python3
"""Inventory syng source paths and run a bounded explicit-query candidate pilot.

Only preliminary candidate intervals are reported: no homology, spanning,
source-spelling, chromosome, ploidy, or genotype validation is performed.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shlex
import shutil
import subprocess
import sys
import time


class BaselineError(ValueError):
    """Invalid input or an unsuccessful pilot; never accept partial results."""


CHECKS = {
    "homology": "pending", "spanning_support": "pending",
    "source_spelling": "pending", "chromosome_identity": "pending",
    "ploidy": "pending", "assembly_provenance": "pending",
}


def metadata():
    return {"schema_version": 1, "interval_semantics": "preliminary candidate intervals",
            "coordinates": "0-based half-open source coordinates",
            "validation_checks": dict(CHECKS)}


def uint(text, context, maximum=2**64 - 1):
    if not re.fullmatch(r"[0-9]+", text) or int(text) > maximum:
        raise BaselineError(f"{context}: invalid unsigned integer {text!r}")
    return int(text)


def load_names(path):
    """Read current 7-column and legacy 3/6-column syng tables without renaming."""
    records, ids, names = [], set(), set()
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            context = f"{path}:{line_no}"
            if len(fields) not in (3, 6, 7):
                raise BaselineError(f"{context}: expected 3, 6, or 7 names columns")
            source_id = uint(fields[0], context, 2**32 - 1)
            name = fields[1]
            length = uint(fields[2], context)
            if not name.strip() or any(ord(c) < 32 or ord(c) == 127 for c in name):
                raise BaselineError(f"{context}: empty name or control character")
            if length == 0:
                raise BaselineError(f"{context}: length must be positive")
            if source_id in ids or name in names:
                raise BaselineError(f"{context}: duplicate source ID or exact path name")
            ids.add(source_id)
            names.add(name)
            parts = name.split("#")
            # Conservative whole-name PanSN recognition; never strip descriptions.
            valid = len(parts) == 3 and all(parts) and not any(c.isspace() for c in name)
            pansn = dict(zip(("sample", "haplotype", "contig"), parts)) if valid else None
            records.append({"source_id": source_id, "name": name, "length": length,
                            "pansn": pansn, "names_columns": len(fields)})
    if not records:
        raise BaselineError(f"{path}: empty names table")
    return records


def fresh_directory(path):
    path = Path(path).absolute()
    if path.is_symlink() or (path.exists() and (not path.is_dir() or any(path.iterdir()))):
        raise BaselineError(f"refusing to overwrite nonempty or unsafe output directory: {path}")
    path.mkdir(parents=True, exist_ok=True)
    # Atomic reservation closes the check/mkdir race between cooperating runs.
    # Keep the marker even on failure: interrupted outputs are never reusable.
    try:
        with (path / ".baseline-reserved").open("x", encoding="utf-8") as marker:
            marker.write("Reserved by yeast_inference_baseline.py; use a fresh directory.\n")
    except FileExistsError as exc:
        raise BaselineError(f"output directory already reserved: {path}") from exc
    return path


def write_json(path, data):
    path = Path(path)
    temporary = path.with_name(path.name + ".tmp")
    try:
        temporary.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        temporary.replace(path)
    finally:
        temporary.unlink(missing_ok=True)


def write_tsv(path, columns, rows):
    with Path(path).open("w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")
        writer.writerow(columns)
        writer.writerows(rows)


def fingerprint(path, checksum=False):
    path = Path(path).resolve(strict=True)
    if not path.is_file():
        raise BaselineError(f"not a regular file: {path}")
    stat = path.stat()
    result = {"path": str(path), "size_bytes": stat.st_size, "mtime_ns": stat.st_mtime_ns}
    if checksum:
        digest = hashlib.sha256()
        with path.open("rb") as handle:
            for block in iter(lambda: handle.read(1024 * 1024), b""):
                digest.update(block)
        result["sha256"] = digest.hexdigest()
    return result


def inventory_data(records):
    parsed = [r for r in records if r["pansn"] is not None]
    samples = Counter(r["pansn"]["sample"] for r in parsed)
    haplotypes = {(r["pansn"]["sample"], r["pansn"]["haplotype"]) for r in parsed}
    return {**metadata(), "source_path_count": len(records),
            "source_bases": sum(r["length"] for r in records),
            "parsed_pansn_path_count": len(parsed),
            "parsed_sample_count": len(samples), "parsed_sample_path_counts": dict(sorted(samples.items())),
            "parsed_sample_haplotype_count": len(haplotypes),
            "unparsed_names": [r["name"] for r in records if r["pansn"] is None],
            "legacy_anchor_offset_names": [r["name"] for r in records if r["names_columns"] == 6],
            "source_ids_match_row_order": all(r["source_id"] == i for i, r in enumerate(records)),
            "ambiguities": ["All chromosome/component classifications and biological ploidy remain unknown.",
                            "Parsed sample/haplotype tokens are naming identities, not verified assembly counts.",
                            "Non-PanSN and whitespace-bearing names are retained exactly, with identity unparsed.",
                            "Legacy six-column names lack first-syncmer offsets; inventory is not index validation."],
            "sources": records}


def write_inventory(out, records, names_path):
    data = inventory_data(records)
    data["names_input"] = fingerprint(names_path, checksum=True)
    write_json(out / "inventory.json", data)
    write_tsv(out / "sequences.tsv",
              ["source_id", "name", "length", "sample", "haplotype", "contig", "identity_status"],
              ([r["source_id"], r["name"], r["length"],
                *((r["pansn"] or {}).get(k, "") for k in ("sample", "haplotype", "contig")),
                "parsed_pansn" if r["pansn"] else "unparsed"] for r in records))


def make_seeds(records, reference, start, window_size, max_windows):
    by_name = {r["name"]: r for r in records}
    if reference not in by_name:
        raise BaselineError(f"unknown exact reference name: {reference!r}")
    length = by_name[reference]["length"]
    if start < 0 or start >= length or window_size <= 0 or max_windows <= 0:
        raise BaselineError("start must be in bounds; window size and max windows must be positive")
    seeds = []
    for i in range(max_windows):
        left = start + i * window_size
        if left >= length:
            break
        seeds.append({"name": reference, "start": left, "end": min(left + window_size, length),
                      "label": f"seed_{i:06d}"})
    return seeds


def bed_rows(path, records, columns):
    by_name = {r["name"]: r for r in records}
    with Path(path).open(encoding="utf-8") as handle:
        for line_no, line in enumerate(handle, 1):
            if not line.strip():
                continue
            fields = line.rstrip("\r\n").split("\t")
            # A legal exact source name may itself begin with '#'.
            if line.startswith("#") and fields[0] not in by_name:
                continue
            context = f"{path}:{line_no}"
            if len(fields) != columns:
                raise BaselineError(f"{context}: expected BED{columns}")
            name, start_text, end_text, label = fields[:4]
            if name not in by_name:
                raise BaselineError(f"{context}: unknown exact source name {name!r}")
            start, end = uint(start_text, context), uint(end_text, context)
            if not 0 <= start < end <= by_name[name]["length"]:
                raise BaselineError(f"{context}: interval out of bounds or empty")
            if not label or any(c.isspace() or ord(c) < 32 for c in label):
                raise BaselineError(f"{context}: invalid seed label")
            row = {"name": name, "source_id": by_name[name]["source_id"],
                   "start": start, "end": end, "label": label, "line_number": line_no}
            if columns == 6:
                if fields[4] != ".":
                    uint(fields[4], context, 1000)
                if fields[5] not in ("+", "-"):
                    raise BaselineError(f"{context}: expected + or - strand")
                row.update(score=fields[4], strand=fields[5])
            yield row


def load_seeds(path, records):
    seeds = list(bed_rows(path, records, 4))
    if not seeds or len({s["label"] for s in seeds}) != len(seeds):
        raise BaselineError("seed BED4 must be nonempty with unique labels")
    return seeds


def union_bases(intervals):
    total, right = 0, -1
    for start, end in sorted(intervals):
        total += max(0, end - max(start, right))
        right = max(right, end)
    return total


def summarize_data(bed, records, seeds, max_candidate_rows):
    if max_candidate_rows <= 0:
        raise BaselineError("max candidate rows must be positive")
    by_label = {s["label"]: [] for s in seeds}
    candidates = []
    for row in bed_rows(bed, records, 6):
        if row["label"] not in by_label:
            raise BaselineError(f"{bed}:{row['line_number']}: unknown seed label {row['label']!r}")
        if len(candidates) >= max_candidate_rows:
            raise BaselineError(f"candidate row limit exceeded ({max_candidate_rows})")
        # Each emitted row is an occurrence, even if its interval repeats exactly.
        row["occurrence_id"] = f"candidate_{len(candidates):09d}"
        row["length"] = row["end"] - row["start"]
        candidates.append(row)
        by_label[row["label"]].append(row)
    per_seed = []
    for seed in seeds:
        hits = by_label[seed["label"]]
        lengths = [r["length"] for r in hits]
        exact = sum(r["name"] == seed["name"] and r["start"] == seed["start"]
                    and r["end"] == seed["end"] for r in hits)
        failures = ["no_candidates"] if not hits else []
        if hits and exact == len(hits):
            failures.append("only_exact_seed_intervals")
        per_seed.append({**seed, "candidate_count": len(hits),
                         "unique_path_count": len({r["name"] for r in hits}),
                         "unique_paths": sorted({r["name"] for r in hits}),
                         "strand_counts": {s: sum(r["strand"] == s for r in hits) for s in ("+", "-")},
                         "length_min": min(lengths, default=None), "length_max": max(lengths, default=None),
                         "length_sum": sum(lengths), "length_mean": sum(lengths) / len(lengths) if lengths else None,
                         "exact_seed_interval_count": exact, "failures": failures})
    seed_intervals = {r["name"]: [] for r in records}
    hit_intervals = {r["name"]: [] for r in records}
    for row in seeds:
        seed_intervals[row["name"]].append((row["start"], row["end"]))
    for row in candidates:
        hit_intervals[row["name"]].append((row["start"], row["end"]))
    coverage = []
    for source in records:
        seed_bp = union_bases(seed_intervals[source["name"]])
        hit_bp = union_bases(hit_intervals[source["name"]])
        coverage.append({"source_id": source["source_id"], "name": source["name"], "length": source["length"],
                         "seed_union_bases": seed_bp, "unseeded_bases": source["length"] - seed_bp,
                         "candidate_union_bases": hit_bp,
                         "candidate_uncovered_bases": source["length"] - hit_bp})
    return {**metadata(), "candidate_count": len(candidates), "per_seed": per_seed,
            "source_coverage": coverage,
            "coverage_caveat": "Coordinate union only, not homologous, callable, or spanning coverage.",
            "diagnostic_caveat": "No candidates or only exact seed intervals do not prove biological absence.",
            "seed_membership_policy": "Seed presence among hits is not required: the current refined CLI can omit the already-visited seed tuple. Seeds remain separate from emitted candidates; a future inference catalog must explicitly retain its seed candidate."}, candidates


def discard_summary(out):
    for name in ("summary.json", "summary.json.tmp", "candidates.tsv", "candidates.tsv.tmp"):
        (out / name).unlink(missing_ok=True)


def write_summary(out, summary, candidates):
    columns = ["occurrence_id", "label", "source_id", "name", "start", "end", "strand", "length", "score", "line_number"]
    try:
        write_tsv(out / "candidates.tsv.tmp", columns + ["eligibility_status", "homology_check", "spanning_check"],
                  ([r[c] for c in columns] + ["preliminary_candidate_interval", "pending", "pending"] for r in candidates))
        (out / "candidates.tsv.tmp").replace(out / "candidates.tsv")
        # Publish acceptance only after the complete candidate table exists.
        write_json(out / "summary.json", summary)
    except (OSError, UnicodeError):
        discard_summary(out)
        raise


def run_command(argv, stdout_path, stderr_path, timeout, manifest, stage, out):
    manifest["stage"] = stage
    write_json(out / "manifest.json", manifest)
    started = time.monotonic()
    try:
        with stdout_path.open("wb") as stdout, stderr_path.open("wb") as stderr:
            result = subprocess.run(argv, stdout=stdout, stderr=stderr, timeout=timeout, check=False)
        manifest[f"{stage}_returncode"] = result.returncode
        if result.returncode:
            raise BaselineError(f"{stage} failed with exit code {result.returncode}; see {stderr_path}")
    except subprocess.TimeoutExpired as exc:
        raise BaselineError(f"{stage} timed out after {timeout} seconds; see {stderr_path}") from exc
    finally:
        manifest[f"{stage}_elapsed_seconds"] = time.monotonic() - started


def run_pilot(args):
    if args.threads <= 0 or not 0 <= args.merge_distance <= 2**31 - 1:
        raise BaselineError("threads must be positive and merge distance in [0, 2147483647]")
    if not math.isfinite(args.timeout) or args.timeout <= 0 or args.max_candidate_rows <= 0:
        raise BaselineError("timeout must be finite and positive; max candidate rows must be positive")
    prefix = str(Path(args.index_prefix).resolve())
    sidecar = prefix if prefix.endswith(".syng") else prefix + ".syng"
    names_path = Path(sidecar + ".names")
    records = load_names(names_path)
    if not all(r["source_id"] == i for i, r in enumerate(records)):
        raise BaselineError("pilot requires source IDs to match row order (impg's names loader uses row order)")
    seeds = make_seeds(records, args.reference, args.start, args.window_size, args.max_windows)
    inputs = [fingerprint(prefix + ext) for ext in (".1gbwt", ".1khash")]
    inputs += [fingerprint(sidecar + ext) for ext in (".spos", ".pstep", ".meta")]
    inputs += [fingerprint(names_path, checksum=True), fingerprint(args.agc)]
    executable = shutil.which(args.impg)
    if executable is None:
        raise BaselineError(f"executable not found: {args.impg}")
    executable = str(Path(executable).resolve())
    out = fresh_directory(args.out_dir)
    seed_path = out / "seeds.bed"
    seed_path.write_text("".join(f"{s['name']}\t{s['start']}\t{s['end']}\t{s['label']}\n" for s in seeds), encoding="utf-8")
    argv = [executable, "query", "-a", prefix, "-b", str(seed_path),
            "--sequence-files", str(Path(args.agc).resolve()), "-o", "bed",
            "-t", str(args.threads), "-d", str(args.merge_distance), "--consider-strandness"]
    manifest = {**metadata(), "status": "running", "stage": "preparation", "argv": argv, "version_argv": [executable, "--version"],
                "cwd": str(Path.cwd()), "inputs": inputs, "seeds": seeds,
                "parameters": {k: v for k, v in vars(args).items() if k != "handler"},
                "harness": fingerprint(__file__, checksum=True), "executable": fingerprint(executable, checksum=True),
                "seed_bed": fingerprint(seed_path, checksum=True),
                "identity_caveat": "Large index/AGC files identified by path, size, mtime only; reuse the panel build checksum manifest.",
                "environment_caveat": "Inherits caller environment; not a hermetic execution.",
                "query_scope": "Explicit seeds against supplied full index, not partition discovery or panel completeness validation."}
    write_json(out / "argv.json", argv)
    (out / "command.sh").write_text(f"#!/bin/sh\ncd {shlex.quote(str(Path.cwd()))} && " + shlex.join(argv)
                                     + f" > {shlex.quote(str(out / 'candidates.bed'))} 2> {shlex.quote(str(out / 'query.stderr.log'))}\n", encoding="utf-8")
    try:
        write_inventory(out, records, names_path)
        run_command(manifest["version_argv"], out / "version.txt", out / "version.stderr.log",
                    args.timeout, manifest, "version", out)
        manifest["version"] = (out / "version.txt").read_text(encoding="utf-8").strip()
        if not manifest["version"]:
            raise BaselineError("version command produced no version text")
        run_command(argv, out / "candidates.bed", out / "query.stderr.log",
                    args.timeout, manifest, "query", out)
        manifest["stage"] = "output_validation"
        summary, candidates = summarize_data(out / "candidates.bed", records, seeds, args.max_candidate_rows)
        summary["command_status"] = "succeeded"
        write_summary(out, summary, candidates)
        manifest["status"] = "succeeded"
    except (BaselineError, OSError, UnicodeError) as exc:
        discard_summary(out)
        manifest.update(status="failed", error=str(exc))
        write_json(out / "failure.json", {**metadata(), "stage": manifest["stage"], "error": str(exc),
                                          "per_seed": [{**s, "candidate_count": None,
                                                        "failures": ["pilot_failed_results_not_accepted"]} for s in seeds]})
        raise
    finally:
        try:
            write_json(out / "manifest.json", manifest)
        except (OSError, UnicodeError):
            discard_summary(out)
            raise


def run_inventory(args):
    records = load_names(args.names)
    write_inventory(fresh_directory(args.out_dir), records, args.names)


def run_summary(args):
    records = load_names(args.names)
    seeds = load_seeds(args.seeds, records)
    out = fresh_directory(args.out_dir)
    try:
        summary, candidates = summarize_data(args.bed, records, seeds, args.max_candidate_rows)
        summary["command_status"] = "not_verified_external_bed"
        summary["inputs"] = [fingerprint(p, checksum=True) for p in (args.names, args.seeds, args.bed)]
        write_summary(out, summary, candidates)
    except (BaselineError, OSError, UnicodeError) as exc:
        write_json(out / "failure.json", {**metadata(), "error": str(exc),
                                          "per_seed": [{**s, "candidate_count": None,
                                                        "failures": ["invalid_output"]} for s in seeds]})
        raise


def parser():
    result = argparse.ArgumentParser(description=__doc__)
    commands = result.add_subparsers(dest="command", required=True)
    inv = commands.add_parser("inventory", help="Inventory exact names and naming ambiguity")
    inv.add_argument("--names", required=True)
    inv.set_defaults(handler=run_inventory)
    pilot = commands.add_parser("pilot", help="Generate BED4 seeds and run one batched impg query")
    pilot.add_argument("--index-prefix", required=True)
    pilot.add_argument("--agc", required=True)
    pilot.add_argument("--impg", default="impg")
    pilot.add_argument("--reference", required=True, help="Exact complete name from the names table")
    for name in ("start", "window-size", "max-windows", "threads", "merge-distance"):
        pilot.add_argument("--" + name, type=int, required=True)
    pilot.add_argument("--timeout", type=float, required=True, help="Seconds per subprocess, including version")
    pilot.set_defaults(handler=run_pilot)
    summary = commands.add_parser("summarize", help="Validate and summarize external labelled BED6 candidates")
    for name in ("names", "seeds", "bed"):
        summary.add_argument("--" + name, required=True)
    summary.set_defaults(handler=run_summary)
    for cmd in (pilot, summary):
        cmd.add_argument("--max-candidate-rows", type=int, default=1_000_000,
                         help="Reject more than this many output rows (default: 1000000)")
    for cmd in (inv, pilot, summary):
        cmd.add_argument("--out-dir", required=True, help="New or empty output directory; never overwrite")
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    try:
        args.handler(args)
    except (BaselineError, OSError, UnicodeError) as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    print(f"{args.command}: {Path(args.out_dir).absolute()}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
