#!/usr/bin/env python3
# THE DOSAGE-SURFACE GATE (per component): the thin-layer proof +
# the QC arithmetic + the agreement table. The dosage layer
# (`impg genome-infer emit-dosage`) consumes the chain's
# molecules.jsonl and emits the per-panel-segment copy counts; this
# gate proves the thin-layer property and renders the evaluation
# table from the test-mode artifact.
#
#   1. THE THIN-LAYER PROOF: the dosage is DERIVED from the
#      molecules — every segment, visit, copy count, per-locus
#      dosage view, and per-path summary is RE-DERIVED HERE in
#      independence from molecules.jsonl (+ the axis windows from
#      the partition-graphs, the chain's own derivation) and
#      compared field-identically with dosage.jsonl.
#      (The census-derived observed mass is the layer's own new
#      measurement — the same discipline as the chain gate, which
#      does not re-derive census votes; its ARITHMETIC is proven
#      below: the sums, the derived mu, and every expected/ratio
#      field re-computed bit-exactly from the emitted integers.)
#   2. THE QC ARITHMETIC: mu == attributed mass / emitted copy-bp
#      (the DERIVED per-covered-copy expectation, no constant);
#      expected_s == mu * dosage_s * length_s; ratio_s == observed_s
#      / expected_s (null when expected is 0); SUM(observed) ==
#      attributed mass; SUM(dosage * length) == copy-bp == the
#      summed member-row lengths (the copy-bp identity).
#   3. ZERO TRUTH KEYS in the default emission (no key naming
#      truth/rank1/perfect/assignment anywhere in dosage.jsonl).
#   4. THE AGREEMENT TABLE (assessment-side, the test-mode
#      artifact): the carried product columns match the product
#      run's own truth-QV artifact; the class counts sum to the
#      loci; the per-segment classes sum to the union universe;
#      the correct-dosage fraction and the rank-1 agreement are
#      re-tallied from the per-locus lines.
#
# Usage: dosage-gate.py <component> <molecules.jsonl> <dosage-dir>
#                        <partition-graphs-dir> <truth-qv.jsonl>
import json, sys

def fnv1a64(path):
    h = 0xcbf29ce484222325
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            for b in chunk:
                h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def merge_intervals(spans):
    merged = []
    for lo, hi in spans:
        if merged and lo <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], hi)
        else:
            merged.append([lo, hi])
    return [(lo, hi) for lo, hi in merged]

def elementary_segments(per_path):
    """The covered consecutive breakpoint pairs per path (the same
    rule the emission states: the gaps between disjoint rows are not
    material)."""
    segments = []
    for path in sorted(per_path):
        rows = per_path[path]
        bounds = sorted({b for row in rows for b in row})
        covered = merge_intervals(sorted(rows))
        for lo, hi in covered:
            inside = [b for b in bounds if lo <= b <= hi]
            for a, b in zip(inside, inside[1:]):
                segments.append((path, a, b))
    return segments

def rows_by_path(members):
    per_path = {}
    for member in members:
        if member["start"] < member["end"]:
            per_path.setdefault(member["path_name"], []).append(
                (member["start"], member["end"]))
    return per_path

def locus_dosage_view(fold_rows):
    per_path = {}
    for fold in fold_rows:
        for path, rows in fold.items():
            per_path.setdefault(path, []).extend(rows)
    segments = elementary_segments(per_path)
    histogram, covered_bp, flat2 = {}, 0, bool(segments)
    for path, start, end in segments:
        copies = sum(containing_rows(fold.get(path, []), start, end)
                     for fold in fold_rows)
        covered_bp += end - start
        histogram[copies] = histogram.get(copies, 0) + 1
        if copies != 2:
            flat2 = False
    return {"segment_count": len(segments), "covered_bp": covered_bp,
            "flat_class_of_2": flat2, "dosage_histogram": histogram}

def containing_rows(rows, start, end):
    return sum(1 for s, e in rows if s <= start and e >= end)

def main():
    component, molecules_path, dosage_dir, graphs_dir, truth_path = sys.argv[1:6]
    checks = 0
    failures = []
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            failures.append(label)

    m = json.loads(open(molecules_path).read().strip())
    check(m["component"].endswith("#" + component), "component identity")
    d = json.loads(open(f"{dosage_dir}/dosage.jsonl").read().strip())
    check(d["model"] == "diplotype-dosage-emission-v1", "model identity")

    # ---- 1. the thin-layer proof (bit-reproducible from molecules)
    fp = fnv1a64(molecules_path)
    check(d["molecules"]["fnv1a64"] == format(fp, "016x"),
          f"molecules fingerprint ({d['molecules']['fnv1a64']} vs {fp:016x})")
    check(d["molecules"]["bytes"] == len(open(molecules_path, "rb").read()),
          "molecules byte count")

    # the axis windows (the chain's own derivation)
    axis_rows = []
    import os
    for name in os.listdir(graphs_dir):
        if not name.endswith(".gfa.map.json"):
            continue
        with open(os.path.join(graphs_dir, name)) as f:
            amap = json.load(f)
        for axis in [x for x in amap["members"] if x["path_name"] == m["component"]]:
            axis_rows.append((amap["partition"], axis["start"], axis["end"]))
    axis_rows.sort(key=lambda x: (x[1], x[2], x[0]))
    window_axis = [(s, e) for _, s, e in axis_rows]

    # the derivation itself
    fold_rows = [
        [rows_by_path(l["members"]) for l in m["molecules"][0]["loci"]],
        [rows_by_path(l["members"]) for l in m["molecules"][1]["loci"]],
    ]
    universe_paths = {}
    for molecule in fold_rows:
        for locus_rows in molecule:
            for path, rows in locus_rows.items():
                universe_paths.setdefault(path, []).extend(rows)
    universe = elementary_segments(universe_paths)
    per_path = {path: sorted({b for s, e in rows for b in (s, e)})
                for path, rows in universe_paths.items()}
    covered_union = {path: merge_intervals(sorted(rows))
                     for path, rows in universe_paths.items()}
    global_index = {}
    slot_base, slot = {}, {}
    for index, (path, start, end) in enumerate(universe):
        if path not in slot:
            slot[path] = len(slot)
            slot_base[path] = index
        global_index[(path, start)] = index
    visits = [[] for _ in universe]
    slot_base, slot_starts = {}, {}
    for index, (path, start, end) in enumerate(universe):
        if path not in slot_base:
            slot_base[path] = index
            slot_starts[path] = []
        slot_starts[path].append(start)
    for molecule_number, molecule in enumerate(fold_rows, start=1):
        for locus_index, locus_rows in enumerate(molecule):
            entry = m["molecules"][molecule_number - 1]["loci"][locus_index]
            axis = window_axis[entry["locus"]]
            for path in sorted(locus_rows):
                rows = locus_rows[path]
                starts = slot_starts[path]
                ends = [universe[slot_base[path] + i][2]
                        for i in range(len(starts))]
                for row_start, row_end in rows:
                    lo = sum(1 for b in starts if b < row_start)
                    hi = sum(1 for b in starts if b < row_end)
                    for position in range(lo, hi):
                        index = slot_base[path] + position
                        if ends[position] <= row_end:
                            visits[index].append({
                                "molecule": molecule_number,
                                "locus": entry["locus"],
                                "partition": entry["partition"],
                                "fold_index": entry["fold_index"],
                                "axis_window": [axis[0], axis[1]],
                            })
    dosage = [len(v) for v in visits]
    check([[s["path_name"], s["start"], s["end"]] for s in d["segments"]] ==
          [[p, s, e] for p, s, e in universe],
          "the segment universe re-derived from molecules.jsonl")
    check([s["dosage"] for s in d["segments"]] == dosage,
          "the per-segment copy counts re-derived (repeat-visit multiplicity)")
    for index, segment in enumerate(d["segments"]):
        got = segment["visits"]
        check(got == visits[index],
              f"segment {segment['path_name']}:{segment['start']} visits")
        check(segment["visits_molecule1"] == sum(1 for v in got if v["molecule"] == 1)
              and segment["visits_molecule2"] == sum(1 for v in got if v["molecule"] == 2),
              f"segment {segment['path_name']}:{segment['start']} per-molecule counts")
        check(segment["axis_extent"] == [
            min(v["axis_window"][0] for v in got),
            max(v["axis_window"][1] for v in got),
        ] if got else [segment["start"], segment["end"]],
              f"segment {segment['path_name']}:{segment['start']} axis extent")
    # the copy-bp identity
    row_bp = sum(e - s for molecule in fold_rows for locus in molecule
                 for rows in locus.values() for s, e in rows)
    copy_bp = sum(c * (e - s) for c, (_, s, e) in zip(dosage, universe))
    check(copy_bp == row_bp == d["depth_qc"]["emitted_copy_bp"]
          == d["depth_qc"]["row_bp"], "the copy-bp identity")
    # the per-locus views
    for index, meta in enumerate(m["loci"]):
        emitted = d["loci"][index]
        check(emitted["locus"] == meta["locus"], "the locus order")
        view = locus_dosage_view([fold_rows[0][index], fold_rows[1][index]])
        emitted_histogram = {int(k): v
                             for k, v in emitted["dosage_histogram"].items()}
        check(emitted["segment_count"] == view["segment_count"]
              and emitted["covered_bp"] == view["covered_bp"]
              and emitted["flat_class_of_2"] == view["flat_class_of_2"]
              and emitted_histogram == view["dosage_histogram"],
              f"locus {meta['locus']}: the emitted dosage view re-derived")
        check(emitted["homozygous_pair"] ==
              (meta["fold_indices"][0] == meta["fold_indices"][1]),
              f"locus {meta['locus']}: the homozygous-pair flag")
    # the per-path summaries
    for entry in d["paths"]:
        rows = universe_paths[entry["path_name"]]
        covered = sum(hi - lo for lo, hi in merge_intervals(sorted(rows)))
        entry_bp = sum(c * (e - s) for c, (p, s, e) in
                       zip(dosage, universe) if p == entry["path_name"])
        check(entry["covered_bp"] == covered,
              f"path {entry['path_name']}: covered bp")
        check(entry["copy_bp"] == entry_bp,
              f"path {entry['path_name']}: copy bp")
        check(entry["zero_copy_bp"] ==
              max(0, entry["path_length"] - entry["covered_bp"]),
              f"path {entry['path_name']}: the zero-copy complement")

    # ---- 2. the QC arithmetic (from the emitted integers)
    qc = d["depth_qc"]
    attributed = sum(s["observed_mass"] for s in d["segments"])
    check(attributed == qc["attributed_mass"],
          "the attributed mass equals the summed observed mass")
    check(qc["unattributed_material_mass"] ==
          qc["material_mass"] - qc["attributed_mass"],
          "the unattributed material mass")
    mu = qc["per_covered_copy_bp_expectation"]
    check(mu == qc["attributed_mass"] / qc["emitted_copy_bp"],
          "mu is the derived per-covered-copy expectation")
    for segment in d["segments"]:
        expected = mu * (segment["dosage"] * segment["length"])
        check(segment["expected_mass"] == expected,
              f"segment {segment['path_name']}:{segment['start']} expected mass")
        ratio = (segment["observed_mass"] / expected) if expected > 0 else None
        check(segment["mass_ratio"] == ratio,
              f"segment {segment['path_name']}:{segment['start']} ratio stated not binned")

    # ---- 3. zero truth keys in the default emission
    def truth_keys(value):
        found = []
        if isinstance(value, dict):
            for key, item in value.items():
                if any(t in key for t in ("truth", "rank1", "perfect", "assignment")):
                    found.append(key)
                found.extend(truth_keys(item))
        elif isinstance(value, list):
            for item in value:
                found.extend(truth_keys(item))
        return found
    check(truth_keys(d) == [], "zero truth keys in dosage.jsonl")

    # ---- 4. the agreement table (the test-mode artifact)
    lines = [json.loads(l) for l in open(f"{dosage_dir}/dosage.jsonl.truth-qv.jsonl")]
    summary = lines[-1]
    per_locus = [l for l in lines if not l.get("summary")]
    product = {}
    for line in open(truth_path):
        r = json.loads(line)
        product[r["locus"]] = r
    check(len(per_locus) == len(m["loci"]), "the per-locus table covers the loci")
    for row in per_locus:
        entry = product[row["locus"]]
        check(row["truth_pair_expressible"] == entry["truth_pair_expressible"]
              and row["rank1"] == entry.get("rank1")
              and row["perfect"] == entry.get("perfect"),
              f"locus {row['locus']}: carried product columns unchanged")
    class_counts = {}
    for row in per_locus:
        class_counts[row["class"]] = class_counts.get(row["class"], 0) + 1
    check(summary["class_counts"] == class_counts, "the class counts sum to the loci")
    check(sum(class_counts.values()) == summary["loci"] == len(m["loci"]),
          "the class arithmetic")
    correct = sum(1 for row in per_locus
                  if row["truth_pair_expressible"] and row["dosage_equal"])
    check(summary["correct_dosage"] == correct, "the correct-dosage count")
    check(summary["correct_dosage_fraction"] ==
          correct / summary["truth_pair_expressible"],
          "the correct-dosage fraction")
    rank1_agree = sum(1 for row in per_locus
                      if row["rank1"] and row["dosage_equal"])
    check(summary["rank1_dosage_agree"] == rank1_agree, "the rank-1 agreement")
    check(summary["emitted_flat_class_of_2_loci"] ==
          sum(1 for row in per_locus if row["emitted"]["flat_class_of_2"]),
          "the emitted flat-class-of-2 count")
    check(summary["truth_flat_class_of_2_loci"] ==
          sum(1 for row in per_locus
              if row["truth"] and row["truth"]["flat_class_of_2"]),
          "the truth flat-class-of-2 count")
    check(sum(summary["segment_classes"].values()) == summary["union_segments"],
          "the per-segment classes sum to the union universe")

    # ---- the table
    print(f"# dosage-surface gate: {component}")
    print(f"# thin-layer proof: molecules.jsonl FNV {d['molecules']['fnv1a64']} "
          f"({d['molecules']['bytes']} B) consumed unchanged; every segment, "
          f"visit, copy count, and per-locus view re-derived from the molecules "
          f"bit-identically")
    print(f"# the emission: {len(d['segments'])} segments, copy-bp "
          f"{qc['emitted_copy_bp']} (the row-bp identity exact), "
          f"{len(d['paths'])} panel paths, "
          f"{sum(1 for l in d['loci'] if l['flat_class_of_2'])} flat-class-of-2 loci "
          f"of {len(d['loci'])}")
    print(f"# the depth QC: mu = {mu} per covered copy-bp (derived), attributed "
          f"{qc['attributed_mass']} of {qc['material_mass']} census material mass "
          f"({qc['unattributed_material_mass']} unattributed, the multi-matching "
          f"spread, stated not hidden)")
    print(f"# the agreement table (assessment-side, the test mode): loci "
          f"{summary['loci']}, expressible {summary['truth_pair_expressible']}, "
          f"rank-1 {summary['truth_rank1']}, correct dosage "
          f"{summary['correct_dosage']} ({summary['correct_dosage_fraction']}), "
          f"rank-1 agree {summary['rank1_dosage_agree']}")
    print(f"# classes: {json.dumps(summary['class_counts'])}")
    print(f"# segment classes over the union universe "
          f"({summary['union_segments']}): {json.dumps(summary['segment_classes'])}")
    print("#")
    print("# locus  class                              mechanism             equal  flat2  rank1")
    for row in per_locus:
        print(f"# {row['locus']!s:6}  {row['class']:34}  {row['mechanism']:20}  "
              f"{str(row['dosage_equal']):5}  {str(row['emitted']['flat_class_of_2']):5}  "
              f"{str(row['rank1']):5}")
    print(f"TOTAL CHECKS: {checks}, FAILURES: {len(failures)}")
    for f in failures:
        print(f"FAIL: {f}")
    print("ALL CHECKS PASS" if not failures else "GATE FAILED")
    return 0 if not failures else 1

if __name__ == "__main__":
    sys.exit(main())
