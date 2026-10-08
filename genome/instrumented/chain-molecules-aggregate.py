#!/usr/bin/env python3
# THE CHAIN-MOLECULES AGGREGATE: the 17-component fleet table of the
# chain/phasing layer over the CLI product's closed calls. Per
# component: the per-locus product accuracy (carried UNCHANGED from
# the product run's own truth-QV artifact) AND the chain's own
# boundary evaluation (assessable boundaries, switches, brackets) —
# separate columns, NEVER conflated (the owner's evaluation ruling:
# per-partition accuracy under best-case assignment AND switch
# errors). The aggregate lands at chain-fleet-AGGREGATE.txt.
import json, os, sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
COMPONENTS = ["chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
              "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII",
              "chrXIII", "chrXIV", "chrXV", "chrXVI"]
TRUTH = {"chrMT": f"{D}/cli-product-chrMT-truthqv/calls.jsonl.truth-qv.jsonl"}

def product_counts(component):
    path = TRUTH.get(component, f"{D}/cli-product-{component}/calls.jsonl.truth-qv.jsonl")
    records = [json.loads(line) for line in open(path)]
    expressible = sum(1 for r in records if r["truth_pair_expressible"])
    rank1 = sum(1 for r in records if r.get("rank1"))
    perfect = sum(1 for r in records if r["truth_pair_expressible"] and r.get("perfect"))
    return len(records), expressible, rank1, perfect

def main():
    rows = []
    for component in COMPONENTS:
        chain_dir = f"{D}/chain-molecules-{component}"
        m = json.loads(open(f"{chain_dir}/molecules.jsonl").read().strip())
        summary = None
        for line in open(f"{chain_dir}/molecules.jsonl.truth-qv.jsonl"):
            r = json.loads(line)
            if r.get("summary"):
                summary = r
        decided = {}
        unidentifiable = 0
        orphan = 0
        for b in m["boundaries"]:
            decided[b["decided_by"]] = decided.get(b["decided_by"], 0) + 1
            unidentifiable += 1 if b["unidentifiable"] else 0
            orphan += b["crossing"]["orphan_occurrences"]
        loci, expressible, rank1, perfect = product_counts(component)
        assert summary["truth_pair_expressible"] == expressible
        assert summary["truth_rank1"] == rank1
        assert summary["perfect"] == perfect
        assert summary["boundaries"] == loci - 1
        rows.append({
            "component": component,
            "loci": loci, "expressible": expressible,
            "rank1": rank1, "perfect": perfect,
            "boundaries": summary["boundaries"],
            "assessable": summary["assessable_boundaries"],
            "switches": summary["switches"],
            "brackets": summary["bracketed_boundaries"],
            "decided": dict(sorted(decided.items())),
            "unidentifiable": unidentifiable,
            "orphan_crossing_occurrences": orphan,
        })
    out = sys.stdout
    print("# THE CHAIN/PHASING LAYER AGGREGATE (the 17-component fleet)", file=out)
    print("# The per-locus product accuracy columns are the COMMITTED numbers, "
          "carried unchanged;", file=out)
    print("# the switch columns are the chain layer's own evaluation over the "
          "same sample —", file=out)
    print("# separate columns, never conflated.", file=out)
    print("#", file=out)
    header = (f"{'component':8s} {'loci':>5s} {'expr':>5s} {'rank1':>6s} {'perf':>5s} "
              f"{'bnd':>4s} {'assess':>7s} {'SWITCH':>7s} {'brackets':>26s} "
              f"{'decided(cross/adj/tie)':>22s} {'unid':>5s} {'orphan':>7s}")
    print(header, file=out)
    tot = {k: 0 for k in ("loci", "expressible", "rank1", "perfect", "boundaries",
                          "assessable", "switches", "unidentifiable",
                          "orphan_crossing_occurrences")}
    for r in rows:
        brackets = ",".join(f"{k}:{v}" for k, v in sorted(r["brackets"].items()))
        dec = "/".join(str(r["decided"].get(k, 0)) for k in ("crossing", "adjacency", "tie"))
        print(f"{r['component']:8s} {r['loci']:5d} {r['expressible']:5d} {r['rank1']:6d} "
              f"{r['perfect']:5d} {r['boundaries']:4d} {r['assessable']:7d} "
              f"{r['switches']:7d} {brackets:26s} {dec:>22s} {r['unidentifiable']:5d} "
              f"{r['orphan_crossing_occurrences']:7d}", file=out)
        for k in tot:
            tot[k] += r[k]
    brackets_total = {}
    for r in rows:
        for k, v in r["brackets"].items():
            brackets_total[k] = brackets_total.get(k, 0) + v
    bt = ",".join(f"{k}:{v}" for k, v in sorted(brackets_total.items()))
    print("-" * len(header), file=out)
    print(f"{'TOTAL':8s} {tot['loci']:5d} {tot['expressible']:5d} {tot['rank1']:6d} "
          f"{tot['perfect']:5d} {tot['boundaries']:4d} {tot['assessable']:7d} "
          f"{tot['switches']:7d} {bt:26s} {'':22s} {tot['unidentifiable']:5d} "
          f"{tot['orphan_crossing_occurrences']:7d}", file=out)
    print("#", file=out)
    print(f"# per-locus product accuracy (UNCHANGED by the chain): "
          f"truth rank-1 {tot['rank1']}/{tot['expressible']} = "
          f"{100.0 * tot['rank1'] / tot['expressible']:.1f}%, perfect {tot['perfect']}",
          file=out)
    print(f"# the chain layer's own number: {tot['switches']} switches over "
          f"{tot['assessable']} assessable boundaries "
          f"({tot['boundaries'] - tot['assessable']} bracketed: {bt})",
          file=out)
    print(f"# unidentifiable boundaries (the flagged ties, carried by the stated rule): "
          f"{tot['unidentifiable']}", file=out)
    print(f"# orphan crossing occurrences (one-sided material the called pairs do not "
          f"express): {tot['orphan_crossing_occurrences']}", file=out)

if __name__ == "__main__":
    main()
