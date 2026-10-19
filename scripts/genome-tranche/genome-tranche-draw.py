#!/usr/bin/env python3
"""THE GENOME-WIDE STRATIFIED SAMPLE - the tranche-1 draw (the owner's
methodology ruling executed: "this is dangerous to be drawing lots of
broad conclusions off of six random loci" - the fix: measure the
product decisions at SAMPLED GENOME SCALE).

THE STRATA come from the COMMITTED LEDGERS ONLY (no new runs feed the
selection; every input is a receipt the campaign banked):

  * the split ledger          RV/prep/accepted.txt (5,844 = the product
                              universe) + bed4-genome.tsv (9,891 windows,
                              chrom/start/end/TID) + rejected-incremental.tsv
                              (transitive-expansion 3,252, the expansion
                              gate's rejected class)
  * the runtime receipts      RV/prep/<T>/split.info.json (the per-window
                              axis/truth-row counts) + receipts/*/cut-stats.json
                              + map/<T>/<T>/cut-stats.json (~2,078 committed
                              cuts inside the accepted universe) + the call
                              exit codes (prep/<T>/call.{exit.code,done})
  * the grid's segdup record  ../grid/territories.bed column 5 (the
                              per-territory segdup fraction; the committed
                              rollup's own class boundary 0.02)
  * the marathon roster        the 16 named deep-lane memory-tail members
                              (datasheet-runtime.md: T00129 T04052 T01895
                              T03855 T04063 T04013 T03803 T03915 T03919
                              T02722 T04059 T03411 T04038 T04066 T02866
                              T03641)
  * the poison census          the 5-member syng segdup-tail prep-terminal
                              class (cuts/poison.tsv: T00197 T01411 T02562
                              T00024 T03679)

THE DRAW RULES (no cherry-picking; the seed and the rules are the doc's):
  * the runnable universe U = accepted - the 5 poison members (the syng
    chain fails deterministically) - the 16 marathon roster members (the
    deep memory-tail class; tranche 2's named addition, carried in the
    wall estimate as the committed bound instead) - the 16 cut.exit
    terminal members (the chr2 AGC-offset class 12 + the rest; the RV cut
    fails deterministically).
  * the strata are DISJOINT, assigned in priority order (the specific
    classes first, the bulk last):
      S_EXP  EXPANSION-GATE-ADJACENT: a bed4 genomic neighbor (same chrom,
             window index +/-2) is a transitive-expansion rejection.
      S_SEG  SEGID-ADJACENT: the grid's own segdup class (fraction > 0.02,
             rv-rollup.py's committed boundary) OR a marathon roster member
             within window index +/-2 on the same chrom.
      S_REP  REPEAT-HEAVY: among U's committed cut receipts, the top decile
             of census_mult/reads_kept (the high-mult census signature; the
             boundary is COMPUTED from the receipts and stated, never
             tuned).
      S_BIND BIND-AMBIGUITY: the ledger-named 7479-class - the RV campaign
             call failed exit 1 (the two-origin-shift / no-unique-key bind
             ambiguity family).
      S_ORD  ORDINARY: everything else.
  * allocation: ORDINARY 22, SEGID 10, REPEAT 8, EXPANSION 8, BIND 8
    (~56 territories, >=8 per stratum, roughly proportional to the
    strata's own sizes with the minimum enforced).
  * the draw: python random.Random(SEED=20261017), each stratum's member
    list SORTED BY TID (the draw is a pure function of the committed
    ledgers + the seed), random.sample(list, n).

usage: genome-tranche-draw.py <RV> <GRID_BED> <OUT.tsv>
"""
import json
import glob
import os
import random
import sys

RV = sys.argv[1]
GRID = sys.argv[2]
OUT = sys.argv[3]
SEED = 20261017

ALLOCATION = [("S_ORD", 22), ("S_SEG", 10), ("S_REP", 8),
              ("S_EXP", 8), ("S_BIND", 8)]

# the named ledger classes (provenance: rollup-final.json + datasheet-runtime.md)
MARATHON = {"T00129", "T04052", "T01895", "T03855", "T04063", "T04013",
            "T03803", "T03915", "T03919", "T02722", "T04059", "T03411",
            "T04038", "T04066", "T02866", "T03641"}
POISON = {"T00197", "T01411", "T02562", "T00024", "T03679"}


def read_tsv(path):
    rows = []
    with open(path) as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if a:
                rows.append(a)
    return rows


def main():
    # ---- the bed4 grid: TID -> (chrom, start, end) + the chrom window order
    bed4 = {}
    by_chrom = {}
    for a in read_tsv(f"{RV}/bed4-genome.tsv"):
        chrom, start, end, tid = a[0], int(a[1]), int(a[2]), a[3]
        bed4[tid] = (chrom, start, end)
        by_chrom.setdefault(chrom, []).append((start, tid))
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    order = {}      # TID -> (chrom, index) in chrom order
    for chrom, lst in by_chrom.items():
        for i, (_, tid) in enumerate(lst):
            order[tid] = (chrom, i)
    n_by_index = {}  # (chrom, index) -> TID
    for tid, (chrom, i) in order.items():
        n_by_index[(chrom, i)] = tid

    def neighbors(tid, radius=2):
        chrom, i = order[tid]
        for j in range(i - radius, i + radius + 1):
            if j != i and (chrom, j) in n_by_index:
                yield n_by_index[(chrom, j)]

    # ---- the accepted universe
    accepted = [l.strip() for l in open(f"{RV}/prep/accepted.txt") if l.strip()]

    # ---- the rejection ledger (the expansion gate's class)
    expansion_rejected = set()
    for a in read_tsv(f"{RV}/prep/rejected-incremental.tsv"):
        if len(a) >= 2 and a[1] == "transitive-expansion":
            expansion_rejected.add(a[0])

    # ---- the grid's segdup fractions
    segdup = {}
    for a in read_tsv(GRID):
        if len(a) >= 5:
            segdup[a[0]] = float(a[4])

    # ---- the per-window split info (axis rows, hg002#2 rows)
    split_rows = {}
    for T in accepted:
        try:
            split_rows[T] = json.load(open(f"{RV}/prep/{T}/split.info.json"))
        except Exception:
            split_rows[T] = {}

    # ---- the committed cut receipts (the census signatures)
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

    # ---- the terminal / failure ledgers
    cut_exit = {T for T in accepted if os.path.exists(f"{RV}/prep/{T}/cut.exit")}
    cut_poison = {T for T in accepted if os.path.exists(f"{RV}/prep/{T}/cut.poison")}
    bind_fail = set()
    guard_fail = set()
    callok = set()
    for T in accepted:
        if os.path.exists(f"{RV}/prep/{T}/call.done"):
            callok.add(T)
        elif os.path.exists(f"{RV}/prep/{T}/call.exit.code"):
            code = open(f"{RV}/prep/{T}/call.exit.code").read().strip()
            if code == "1":
                bind_fail.add(T)
            elif code == "137":
                guard_fail.add(T)

    # ---- the runnable universe
    excluded = {
        "poison": POISON & set(accepted),
        "cut_poison": cut_poison,
        "marathon_roster": MARATHON & set(accepted),
        "cut_exit_terminal": cut_exit,
    }
    universe = [T for T in accepted
                if T not in POISON and T not in cut_poison
                and T not in MARATHON and T not in cut_exit]

    # ---- the repeat-heavy boundary: the top decile of census_mult/reads_kept
    ratios = {}
    for T, s in cutstats.items():
        rk = s.get("reads_kept", 0)
        if T in universe and rk and rk > 0:
            ratios[T] = s["census_mult"] / rk
    ranked = sorted(ratios.values())
    n = len(ranked)
    rep_boundary = ranked[int(n * 0.9)] if n else float("inf")

    # ---- the strata (disjoint, priority order)
    strata = {tag: [] for tag, _ in ALLOCATION}
    for T in sorted(universe):
        nb = set(neighbors(T))
        if nb & expansion_rejected:
            strata["S_EXP"].append(T)
        elif segdup.get(T, 0.0) > 0.02 or (nb & MARATHON):
            strata["S_SEG"].append(T)
        elif ratios.get(T, 0.0) >= rep_boundary:
            strata["S_REP"].append(T)
        elif T in bind_fail:
            strata["S_BIND"].append(T)
        else:
            strata["S_ORD"].append(T)

    # ---- the draw
    rng = random.Random(SEED)
    draw = []
    for tag, count in ALLOCATION:
        pool = sorted(strata[tag])
        picked = rng.sample(pool, min(count, len(pool)))
        for T in picked:
            chrom, start, end = bed4[T]
            draw.append({
                "tid": T, "stratum": tag, "chrom": chrom,
                "start": start, "end": end,
                "segdup": segdup.get(T, 0.0),
                "hg002_2_rows": split_rows.get(T, {}).get("hg002_2_rows"),
                "has_cut_receipt": T in cutstats,
                "census_mult": cutstats.get(T, {}).get("census_mult"),
                "reads_kept": cutstats.get(T, {}).get("reads_kept"),
                "census_mult_per_read": (round(ratios[T], 3)
                                         if T in ratios else None),
                "rv_callok": T in callok,
                "rv_bind_fail": T in bind_fail,
                "axis_path": open(f"{RV}/prep/{T}/axis-path.txt").read().strip()
                             if os.path.exists(f"{RV}/prep/{T}/axis-path.txt") else "",
            })

    # ---- the receipt: the strata census + the draw, printed and saved
    print(f"# THE TRANCHE-1 DRAW (seed {SEED}; the rules in the module docstring)")
    print(f"# accepted {len(accepted)}; excluded: " +
          ", ".join(f"{k} {len(v)}" for k, v in excluded.items() if v) +
          f"; runnable universe {len(universe)}")
    print(f"# repeat-heavy boundary: census_mult/reads_kept >= "
          f"{rep_boundary:.3f} (top decile of {n} committed cut receipts in U)")
    for tag, _ in ALLOCATION:
        print(f"# stratum {tag}: {len(strata[tag])} members in U, "
              f"{sum(1 for d in draw if d['stratum'] == tag)} drawn")
    hdr = ["tid", "stratum", "chrom", "start", "end", "segdup",
           "hg002_2_rows", "has_cut_receipt", "census_mult", "reads_kept",
           "census_mult_per_read", "rv_callok", "rv_bind_fail", "axis_path"]
    with open(OUT, "w") as f:
        f.write("#" + "\t".join(hdr) + "\n")
        for d in draw:
            f.write("\t".join(str(d[k]) for k in hdr) + "\n")
    print(f"# wrote {len(draw)} rows -> {OUT}")


if __name__ == "__main__":
    main()
