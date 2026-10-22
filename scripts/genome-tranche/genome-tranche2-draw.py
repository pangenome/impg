#!/usr/bin/env python3
"""THE GENOME-WIDE STRATIFIED SAMPLE - the TRANCHE-2 draw (the tranche-1
checkpoint's own value order, executed: S_EXP at depth, the S_BIND class
census, and the marathon roster).

THE DRAW RULES (the tranche-1 discipline, extended and stated; no
cherry-picking - the draw is a pure function of the committed ledgers,
the seed, and the tranche-1 manifest):
  * the strata machinery is the tranche-1 draw script VERBATIM (the same
    committed ledgers, the same disjoint priority-ordered strata, the
    same runnable-universe exclusions: poison 5, marathon roster 16,
    cut.exit terminal 16).
  * the POOLS EXCLUDE the tranche-1 drawn members (the 56 rows of the
    tranche-1 manifest, including its one prep-failed member T09025) -
    tranche 2 samples the strata BEYOND tranche 1, never re-draws a
    member that tranche 1 already owns.
  * S_EXP: 22 fresh members, random.Random(SEED=20261020), pool sorted
    by TID, random.sample(pool, 22) - the same stated-seed discipline.
  * S_BIND: the CLASS CENSUS, not a sample - every remaining member of
    the ledger-named 7479-class (the RV campaign's exit-1 bind-ambiguity
    fail ledger; 79 in the runnable universe, 8 already drawn by
    tranche 1 and owned by it, so 71 here). The tranche-2 question is
    the resolve-or-fail-closed verdict PER MEMBER of the whole class.
  * S_MAR: the 16-member marathon roster (the deep-lane memory tail,
    datasheet-runtime.md; a named census, NOT part of the runnable
    universe and NOT part of any stratified estimate - tranche 2 runs
    them through the same machinery to receipt the walls or the honest
    guard-kill class).
  * row order = stage priority: the S_EXP block first (the universe's
    76% stratum, zero tranche-1 receipts), then S_BIND, then S_MAR.

usage: genome-tranche2-draw.py <RV> <GRID_BED> <TRANCHE1_TSV> <OUT.tsv>
"""
import json
import glob
import os
import random
import sys

RV = sys.argv[1]
GRID = sys.argv[2]
T1 = sys.argv[3]
OUT = sys.argv[4]
SEED = 20261020

S_EXP_N = 22

# the named ledger classes (provenance: rollup-final.json + datasheet-runtime.md)
MARATHON = ["T00129", "T04052", "T01895", "T03855", "T04063", "T04013",
            "T03803", "T03915", "T03919", "T02722", "T04059", "T03411",
            "T04038", "T04066", "T02866", "T03641"]
POISON = {"T00197", "T01411", "T02562", "T00024", "T03679"}


def read_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if a and not line.startswith("#"):
                rows.append(a)
    return rows


def main():
    # ---- the bed4 grid + window order (tranche-1 draw VERBATIM)
    bed4 = {}
    by_chrom = {}
    for a in read_tsv(f"{RV}/bed4-genome.tsv"):
        chrom, start, end, tid = a[0], int(a[1]), int(a[2]), a[3]
        bed4[tid] = (chrom, start, end)
        by_chrom.setdefault(chrom, []).append((start, tid))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    order = {}
    for chrom, lst in by_chrom.items():
        for i, (_, tid) in enumerate(lst):
            order[tid] = (chrom, i)
    n_by_index = {v: k for k, v in order.items()}

    def neighbors(tid, radius=2):
        chrom, i = order[tid]
        for j in range(i - radius, i + radius + 1):
            if j != i and (chrom, j) in n_by_index:
                yield n_by_index[(chrom, j)]

    accepted = [l.strip() for l in open(f"{RV}/prep/accepted.txt") if l.strip()]

    expansion_rejected = set()
    for a in read_tsv(f"{RV}/prep/rejected-incremental.tsv"):
        if len(a) >= 2 and a[1] == "transitive-expansion":
            expansion_rejected.add(a[0])

    segdup = {}
    for a in read_tsv(GRID):
        if len(a) >= 5:
            segdup[a[0]] = float(a[4])

    split_rows = {}
    for T in accepted:
        try:
            split_rows[T] = json.load(open(f"{RV}/prep/{T}/split.info.json"))
        except Exception:
            split_rows[T] = {}

    cutstats = {}
    for p in glob.glob(f"{RV}/receipts/*/cut-stats.json"):
        T = os.path.basename(os.path.dirname(p))
        try:
            cutstats[T] = json.load(open(p))
        except Exception:
            pass
    for p in glob.glob(f"{RV}/map/*/*/cut-stats.json"):
        T = p.split("/")[-3]
        if T in accepted and T not in cutstats:
            try:
                cutstats[T] = json.load(open(p))
            except Exception:
                pass

    cut_exit = {T for T in accepted if os.path.exists(f"{RV}/prep/{T}/cut.exit")}
    cut_poison = {T for T in accepted if os.path.exists(f"{RV}/prep/{T}/cut.poison")}
    bind_fail = set()
    callok = set()
    for T in accepted:
        if os.path.exists(f"{RV}/prep/{T}/call.done"):
            callok.add(T)
        elif os.path.exists(f"{RV}/prep/{T}/call.exit.code"):
            code = open(f"{RV}/prep/{T}/call.exit.code").read().strip()
            if code == "1":
                bind_fail.add(T)

    MARATHON_SET = set(MARATHON)
    universe = [T for T in accepted
                if T not in POISON and T not in cut_poison
                and T not in MARATHON_SET and T not in cut_exit]

    ratios = {}
    for T, s in cutstats.items():
        rk = s.get("reads_kept", 0)
        if T in universe and rk and rk > 0:
            ratios[T] = s["census_mult"] / rk
    ranked = sorted(ratios.values())
    n = len(ranked)
    rep_boundary = ranked[int(n * 0.9)] if n else float("inf")

    strata = {"S_EXP": [], "S_SEG": [], "S_REP": [],
              "S_BIND": [], "S_ORD": []}
    for T in sorted(universe):
        nb = set(neighbors(T))
        if nb & expansion_rejected:
            strata["S_EXP"].append(T)
        elif segdup.get(T, 0.0) > 0.02 or (nb & MARATHON_SET):
            strata["S_SEG"].append(T)
        elif ratios.get(T, 0.0) >= rep_boundary:
            strata["S_REP"].append(T)
        elif T in bind_fail:
            strata["S_BIND"].append(T)
        else:
            strata["S_ORD"].append(T)

    # ---- tranche-1 ownership: its drawn members leave the pools
    t1_drawn = {a[0] for a in read_tsv(T1)}

    rng = random.Random(SEED)
    exp_pool = sorted(set(strata["S_EXP"]) - t1_drawn)
    exp_draw = rng.sample(exp_pool, min(S_EXP_N, len(exp_pool)))
    bind_pool = sorted(set(strata["S_BIND"]) - t1_drawn)   # the census: ALL
    mar_pool = sorted(MARATHON_SET & set(accepted))

    print(f"# THE TRANCHE-2 DRAW (seed {SEED}; rules in the module docstring)")
    print(f"# strata census (runnable universe {len(universe)} of {len(accepted)} accepted): "
          + ", ".join(f"{k} {len(strata[k])}" for k in
                      ("S_ORD", "S_SEG", "S_REP", "S_EXP", "S_BIND")))
    print(f"# tranche-1 owns {len(t1_drawn)} members (excluded from pools)")
    print(f"# S_EXP: {len(exp_pool)} available, {len(exp_draw)} drawn (sampled)")
    print(f"# S_BIND: {len(bind_pool)} of {len(strata['S_BIND'])} class members "
          f"(census; {len(strata['S_BIND']) - len(bind_pool)} owned by tranche 1)")
    print(f"# S_MAR: {len(mar_pool)} marathon roster members (census, "
          "outside the runnable universe)")

    hdr = ["tid", "stratum", "chrom", "start", "end", "segdup",
           "hg002_2_rows", "has_cut_receipt", "census_mult", "reads_kept",
           "census_mult_per_read", "rv_callok", "rv_bind_fail", "axis_path"]

    def row(T, tag):
        chrom, start, end = bed4[T]
        return [T, tag, chrom, str(start), str(end),
                str(segdup.get(T, 0.0)),
                str(split_rows.get(T, {}).get("hg002_2_rows")),
                str(T in cutstats),
                str(cutstats.get(T, {}).get("census_mult")),
                str(cutstats.get(T, {}).get("reads_kept")),
                str(round(ratios[T], 3) if T in ratios else "None"),
                str(T in callok), str(T in bind_fail),
                open(f"{RV}/prep/{T}/axis-path.txt").read().strip()
                if os.path.exists(f"{RV}/prep/{T}/axis-path.txt") else ""]

    rows = ([row(T, "S_EXP") for T in exp_draw]
            + [row(T, "S_BIND") for T in bind_pool]
            + [row(T, "S_MAR") for T in mar_pool])
    with open(OUT, "w") as f:
        f.write("#" + "\t".join(hdr) + "\n")
        for r in rows:
            f.write("\t".join(r) + "\n")
    print(f"# wrote {len(rows)} rows -> {OUT}")


if __name__ == "__main__":
    main()
