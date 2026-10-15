#!/usr/bin/env python3
# THE CHAIN-MOLECULES GATE (per component): the thin-layer proof +
# the switch table. The chain (`impg genome-infer chain-molecules`)
# consumes the CLI product's CLOSED calls.jsonl and emits the two
# molecules; this gate proves the thin-layer property and renders
# the evaluation table from the test-mode artifact.
#
#   1. THE THIN-LAYER PROOF: the per-locus calls are bit-identical
#      pre/post chain — the chain's emitted calls fingerprint (FNV,
#      computed at run time from the file it consumed) equals the
#      fingerprint of calls.jsonl NOW (the chain never wrote it), and
#      every molecule locus entry is EXACTLY the emitted called pair
#      (the fold indices, the member rows), reoriented only.
#   2. THE ORIENTATION SELF-CONSISTENCY: the emitted per-locus
#      orientation bits accumulate exactly the emitted per-boundary
#      orientations, and every boundary's orientation is the stated
#      lexicographic decision of its own emitted evidence.
#   3. THE SWITCH TABLE (assessment-side, the test mode's artifact):
#      per-boundary chain vs truth orientation, the honest brackets,
#      and the per-locus accuracy under best-case assignment CARRIED
#      UNCHANGED from the product run's own truth-QV artifact —
#      separate columns, never conflated (the owner's evaluation
#      ruling: per-partition accuracy AND switch errors).
#
# Usage: chain-molecules-gate.py <component> <calls.jsonl> <molecules-dir> <truth-qv.jsonl>
import json, sys

def fnv1a64(path):
    h = 0xcbf29ce484222325
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            for b in chunk:
                h = ((h ^ b) * 0x100000001b3) & 0xFFFFFFFFFFFFFFFF
    return h

def main():
    component, calls_path, chain_dir, truth_path = sys.argv[1:5]
    checks = 0
    failures = []
    def check(ok, label):
        nonlocal checks
        checks += 1
        if not ok:
            failures.append(label)

    calls = {}
    for line in open(calls_path):
        r = json.loads(line)
        calls[r["locus"]] = r

    m = json.loads(open(f"{chain_dir}/molecules.jsonl").read().strip())
    check(m["component"].endswith("#" + component), "component identity")
    check(m["model"] == "diplotype-chain-orientation-v1", "model identity")

    # ---- 1. the thin-layer proof
    fp = fnv1a64(calls_path)
    check(m["calls"]["fnv1a64"] == format(fp, "016x"),
          f"calls fingerprint pre/post chain ({m['calls']['fnv1a64']} vs {fp:016x})")
    check(m["calls"]["bytes"] == len(open(calls_path, "rb").read()),
          "calls byte count")
    loci_ids = sorted(calls)
    check([l["locus"] for l in m["loci"]] == loci_ids, "the called locus set")
    for entry in m["loci"]:
        call = calls[entry["locus"]]
        check(entry["fold_indices"] == call["diplotype"]["fold_indices"],
              f"locus {entry['locus']}: the emitted pair unchanged")
        state = entry["orientation_bit"]
        expected = call["diplotype"]["fold_indices"][state & 1]
        check(entry["molecule1_fold"] == expected,
              f"locus {entry['locus']}: molecule-1 fold follows the orientation bit")
        check(entry["qual"] == call.get("qual"),
              f"locus {entry['locus']}: QUAL carried unchanged")
        check(entry["tied_classes"] == (call["called_class_count"] > 1),
              f"locus {entry['locus']}: tied-class flag")
    # every molecule locus entry is exactly the emitted fold's members
    for mol in m["molecules"]:
        check(len(mol["loci"]) == len(loci_ids), f"molecule {mol['molecule']} locus count")
        for le in mol["loci"]:
            call = calls[le["locus"]]
            side = [i for i, p in enumerate(call["diplotype"]["paths"])
                    if p["fold_index"] == le["fold_index"]]
            # (A homozygous emitted pair carries the same fold twice;
            # the entry's members must equal that fold's members.)
            check(len(side) >= 1, f"locus {le['locus']}: fold index in the emitted pair")
            if not side:
                continue
            members = call["diplotype"]["paths"][side[0]]["members"]
            got = [(x["path_name"], x["start"], x["end"]) for x in le["members"]]
            want = [(x["path_name"], x["start"], x["end"]) for x in members]
            check(got == want, f"locus {le['locus']}: member rows bit-equal to the call")

    # ---- 2. the orientation self-consistency
    states = [l["orientation_bit"] for l in m["loci"]]
    check(states[0] == 0, "the chain state starts at 0 (the molecule labels' symmetry)")
    for i, b in enumerate(m["boundaries"]):
        want = states[i] ^ (1 if b["orientation"] == "swap" else 0)
        check(states[i + 1] == want,
              f"boundary {b['loci']}: orientation accumulates into the locus bits")
        vs, vx = b["crossing"]["votes_same"], b["crossing"]["votes_swap"]
        asame = (b["adjacency"]["same"]["touches"], b["adjacency"]["same"]["overlap_bp"])
        aswap = (b["adjacency"]["swap"]["touches"], b["adjacency"]["swap"]["overlap_bp"])
        if vs != vx:
            want_orient = "swap" if vx > vs else "same"
            want_by = "crossing"
        elif asame != aswap:
            want_orient = "swap" if aswap > asame else "same"
            want_by = "adjacency"
        else:
            want_orient, want_by = "same", "tie"
        check(b["orientation"] == want_orient and b["decided_by"] == want_by,
              f"boundary {b['loci']}: the stated lexicographic decision")
        check(b["unidentifiable"] == (want_by == "tie"),
              f"boundary {b['loci']}: the unidentifiable flag")

    # ---- 3. the switch table (the test-mode artifact)
    truth = {}
    for line in open(f"{chain_dir}/molecules.jsonl.truth-qv.jsonl"):
        r = json.loads(line)
        if r.get("summary"):
            summary = r
        else:
            truth[tuple(r["loci"])] = r
    product_truth = {}
    for line in open(truth_path):
        r = json.loads(line)
        product_truth[r["locus"]] = r
    expressible = sum(1 for r in product_truth.values() if r["truth_pair_expressible"])
    rank1 = sum(1 for r in product_truth.values() if r.get("rank1"))
    perfect = sum(1 for r in product_truth.values()
                  if r["truth_pair_expressible"] and r.get("perfect"))
    check(summary["boundaries"] == len(m["boundaries"]), "summary boundary count")
    check(summary["truth_pair_expressible"] == expressible and
          summary["truth_rank1"] == rank1 and summary["perfect"] == perfect,
          "the per-locus accuracy summary carried unchanged from the product artifact")
    check(summary["assessable_boundaries"] +
          sum(summary["bracketed_boundaries"].values()) == summary["boundaries"],
          "assessable + brackets = boundaries")
    counted = sum(1 for r in truth.values() if r["switch_error"])
    check(counted == summary["switches"], "the switch count arithmetic")
    for b in m["boundaries"]:
        row = truth[tuple(b["loci"])]
        check(row["chain_orientation"] == b["orientation"],
              f"boundary {b['loci']}: chain orientation matches the product record")

    # ---- the table
    print(f"# chain-molecules gate: {component}")
    print(f"# thin-layer proof: calls.jsonl FNV {m['calls']['fnv1a64']} "
          f"({m['calls']['bytes']} B) bit-identical pre/post chain; "
          f"every molecule locus = the emitted called pair, reoriented only")
    print(f"# per-locus product accuracy (carried unchanged): loci {summary['loci']}, "
          f"expressible {expressible}, truth rank-1 {rank1}, perfect {perfect}")
    print(f"# boundaries {summary['boundaries']}: assessable {summary['assessable_boundaries']}, "
          f"SWITCHES {summary['switches']}, brackets {json.dumps(summary['bracketed_boundaries'])}")
    decided = {}
    for b in m["boundaries"]:
        decided[b["decided_by"]] = decided.get(b["decided_by"], 0) + 1
    print(f"# decided by: {json.dumps(dict(sorted(decided.items())))}")
    print("#")
    print("# loci        chain  by          truth        switch  assessable")
    for b in m["boundaries"]:
        row = truth[tuple(b["loci"])]
        t = row["truth_orientation"]
        tshow = t if isinstance(t, str) else f"bracket:{t['bracket']}"
        print(f"# {b['loci']!s:12} {row['chain_orientation']:5s}  {row['chain_decided_by']:10s}  "
              f"{tshow:20s} {str(row['switch_error']):6s}  {row['assessable']}")
    print(f"TOTAL CHECKS: {checks}, FAILURES: {len(failures)}")
    for f in failures:
        print(f"FAIL: {f}")
    print("ALL CHECKS PASS" if not failures else "GATE FAILED")
    return 0 if not failures else 1

if __name__ == "__main__":
    sys.exit(main())
