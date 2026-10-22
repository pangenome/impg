#!/usr/bin/env python3
"""THE COMPLEXITY LEDGER's receipt miner: parses the call's own stderr
receipt lines (the serial tranche callC/callA logs, the bounded prod
logs, the local-frame CLI logs) into a per-locus TSV + a per-run phase
TSV, for the exponent fits. One row per locus:
  territory, form, locus, members, folds_pre, folds, records, units,
  matrix_pairs, scored_placed, no_pin, edge_discarded, merge_groups,
  admitted_pairs, equal_col_candidates, classes, scan_pairs, scan_kept,
  scan_skipped, scan_pruned_pct, stage_one_s, scan_wall_s, locus_wall_s.
Receipt-side only; reads logs, never runs anything."""
import sys, re, os, json, argparse

RE_MEMBERS = re.compile(r"locus (\d+): partition \d+, (\d+) members -> (\d+) folds \[([\d.]+)s\]")
RE_CLASSIFY = re.compile(r"locus (\d+): (\d+) records, (\d+) units \(touching")
RE_MATRIX = re.compile(r"locus (\d+): matrix (\d+) folds x (\d+) units = (\d+) pairs, (\d+) live MB")
RE_FACT = re.compile(r"locus (\d+): factorized .*?edge-discarded placements (\d+), no-pin \(unit,fold\) (\d+)")
RE_MERGE = re.compile(r"locus (\d+): alignment-connected folds: (\d+) merge groups \((\d+) admitted pairs of (\d+) equal-column candidates")
RE_BOUNDED = re.compile(r"locus (\d+): class scan bounded: (\d+) pairs in ([\d.]+)s \(stage one ([\d.]+)s; kept (\d+), skipped (\d+) \(([\d.]+)% pruned\)")
RE_DONE = re.compile(r"locus (\d+) done: (\d+) folds, (\d+) units, (\d+) classes,.*\[([\d.]+)s\]")
RE_PHASE = re.compile(r"phase summary: inputs ([\d.]+)s, quality \(FASTQ scan\) ([\d.]+)s, elsewhere ([\d.]+)s, derive cache ([\d.]+)s, binding ([\d.]+)s, variants ([\d.]+)s, rows/context assembly ([\d.]+)s, loci \(folds\+units\+scoring\+classes\+qual\+exactness\+receipts\) ([\d.]+)s; AGC fetches (\d+) total")
RE_INPUTS = re.compile(r"phase inputs: .*?\((\d+) census records, (\d+) windows\) \[([\d.]+)s\]")
RE_BIND = re.compile(r"phase binding: (\d+) records / (\d+) occurrences / (\d+) candidate keys .*? steps ([\d.]+)s .*? verify ([\d.]+)s, fetches (\d+) \(([\d.]+) MB\), lane loads (\d+) \(([\d.]+) MB\), phase ([\d.]+)s")
RE_ROWS = re.compile(r"phase rows: (\d+) distinct rows over (\d+) partitions .*? \[([\d.]+)s\]")
RE_COMPLETE = re.compile(r"complete: loci \[([^\]]*)\].*?total ([\d.]+)s")
RE_PRUNE = re.compile(r"far-bin prune: (\d+) records, (\d+) of (\d+) mult")

def parse_log(path):
    loci = {}
    run = {"log": path}
    with open(path, errors="replace") as fh:
        for line in fh:
            m = RE_MEMBERS.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["members"], d["folds_pre"], d["folds_build_s"] = int(m.group(2)), int(m.group(3)), float(m.group(4))
                continue
            m = RE_CLASSIFY.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["records"], d["units"] = int(m.group(2)), int(m.group(3))
                continue
            m = RE_MATRIX.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["matrix_pairs"], d["live_mb"] = int(m.group(4)), int(m.group(5))
                continue
            m = RE_FACT.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["edge_discarded"], d["no_pin"] = int(m.group(2)), int(m.group(3))
                continue
            m = RE_MERGE.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["merge_groups"], d["admitted_pairs"], d["equal_col_candidates"] = int(m.group(2)), int(m.group(3)), int(m.group(4))
                continue
            m = RE_BOUNDED.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["scan_pairs"], d["scan_wall_s"], d["stage_one_s"] = int(m.group(2)), float(m.group(3)), float(m.group(4))
                d["scan_kept"], d["scan_skipped"], d["scan_pruned_pct"] = int(m.group(5)), int(m.group(6)), float(m.group(7))
                continue
            m = RE_DONE.search(line)
            if m:
                d = loci.setdefault(int(m.group(1)), {})
                d["folds"], d["units"], d["classes"], d["locus_wall_s"] = int(m.group(2)), int(m.group(3)), int(m.group(4)), float(m.group(5))
                continue
            m = RE_PHASE.search(line)
            if m:
                run["inputs_s"], run["quality_s"], run["elsewhere_s"], run["derive_cache_s"] = m.group(1), m.group(2), m.group(3), m.group(4)
                run["binding_s"], run["variants_s"], run["rows_s"], run["loci_s"] = m.group(5), m.group(6), m.group(7), m.group(8)
                run["agc_fetches"] = m.group(9)
                continue
            m = RE_INPUTS.search(line)
            if m:
                run["census_records"], run["windows"], run["inputs_wall_s"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
                continue
            m = RE_BIND.search(line)
            if m:
                run["bind_records"], run["bind_occurrences"], run["bind_keys"] = int(m.group(1)), int(m.group(2)), int(m.group(3))
                run["bind_steps_s"], run["bind_verify_s"] = float(m.group(4)), float(m.group(5))
                run["bind_fetches"], run["bind_fetch_mb"], run["bind_lane_loads"], run["bind_lane_mb"], run["bind_phase_s"] = int(m.group(6)), float(m.group(7)), int(m.group(8)), float(m.group(9)), float(m.group(10))
                continue
            m = RE_ROWS.search(line)
            if m:
                run["distinct_rows"], run["row_partitions"], run["rows_wall_s"] = int(m.group(1)), int(m.group(2)), float(m.group(3))
                continue
            m = RE_COMPLETE.search(line)
            if m:
                run["complete_loci"], run["complete_total_s"] = m.group(1).replace(" ", ""), float(m.group(2))
                continue
            m = RE_PRUNE.search(line)
            if m:
                run["prune_records"], run["prune_mult"], run["prune_mult_total"] = int(m.group(1)), int(m.group(2)), int(m.group(3))
                continue
    run["loci"] = loci
    return run

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("logs", nargs="+")
    args = ap.parse_args()
    cols = ["territory", "form", "log", "locus", "members", "folds_pre", "folds", "records", "units",
            "matrix_pairs", "no_pin", "edge_discarded", "merge_groups", "admitted_pairs",
            "equal_col_candidates", "classes", "scan_pairs", "scan_kept", "scan_skipped",
            "scan_pruned_pct", "stage_one_s", "scan_wall_s", "folds_build_s", "locus_wall_s"]
    rcols = ["territory", "form", "log", "windows", "census_records", "inputs_wall_s", "quality_s",
             "elsewhere_s", "derive_cache_s", "binding_s", "variants_s", "rows_s", "loci_s",
             "agc_fetches", "bind_records", "bind_occurrences", "bind_keys", "bind_steps_s",
             "bind_verify_s", "bind_fetches", "bind_fetch_mb", "bind_lane_loads", "bind_lane_mb",
             "bind_phase_s", "distinct_rows", "rows_wall_s", "complete_total_s",
             "prune_records", "prune_mult", "prune_mult_total"]
    with open(args.out_prefix + "-loci.tsv", "w") as lf, open(args.out_prefix + "-runs.tsv", "w") as rf:
        lf.write("\t".join(cols) + "\n"); rf.write("\t".join(rcols) + "\n")
        for path in args.logs:
            name = os.path.basename(path)
            territory = name.replace("callC-", "").replace("callA-", "").replace("callB-", "").replace(".err", "")
            form = "C" if "callC" in name else ("A" if "callA" in name else ("B" if "callB" in name else "?"))
            run = parse_log(path)
            if "loci_s" not in run:
                continue  # incomplete run (mid-flight master)
            rf.write("\t".join([territory, form, path] + [str(run.get(c, "")) for c in rcols[3:]]) + "\n")
            for locus, d in sorted(run["loci"].items()):
                if "locus_wall_s" not in d:
                    continue
                lf.write("\t".join(str({"territory": territory, "form": form, "log": path, "locus": locus}[c] if c in ("territory", "form", "log", "locus") else d.get(c, "")) for c in cols) + "\n")
    print(f"mined {len(args.logs)} logs -> {args.out_prefix}-loci.tsv / -runs.tsv")

if __name__ == "__main__":
    main()
