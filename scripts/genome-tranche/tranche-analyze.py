#!/usr/bin/env python3
"""THE GENOME TRANCHE - the analysis (the owner's decision-grade receipts):
  (a) THE WASH VERDICT AT SCALE - per stratum, per profile (A = restrict,
      B = global-native): rank-1 rate, the official-allwave QV distribution
      (median / >=QV30 / at-cap, the likegt max-qv convention, cap 50 - the
      same assessment as every banked receipt), the differing loci named and
      truth-judged; THE CONSERVED-ELSEWHERE RATE - the L2-17-class losses
      (B calls the truth pair, A does not) per stratum and per N loci.
  (b) THE CPU-WORK VARIANCE - per-territory serial (RAYON=1) call phases
      from the stage logs (the autopsy's serial-census method: inputs =
      the ex-load panel load, binding/variants/rows/loci/quality the true
      per-territory CPU), per-locus core-s by stratum, and the honest
      GENOME-WALL ESTIMATE WITH BOUNDS from the stratified mean/variance
      (normal-approximation 95% CI on the stratified total; the heavy
      tail receipted by median/p90/max beside it).
  (c) THE FAR-BIN MASS SHARE BY STRATUM - from the restricted (A) cut
      censuses: the records whose placements scatter into NON-ADJACENT
      windows (the prune's exact class; diameter > 1, the census's own
      pair-bins adjacency constant), their multiplicity share per
      territory and stratum - the prune's real genome-scale value.
Every conclusion is labeled with its sample size; nothing extrapolates
beyond its stratum's own receipts except the labeled stratified totals.
usage: tranche-analyze.py <WT> [out.txt]
"""
import json
import math
import os
import re
import statistics
import sys

WT = sys.argv[1]
OUT = sys.argv[2] if len(sys.argv) > 2 else f"{WT}/analysis.txt"

STRATA = ["S_ORD", "S_SEG", "S_REP", "S_EXP", "S_BIND"]
STRATUM_NAMES = {
    "S_ORD": "ORDINARY (the bulk class)",
    "S_SEG": "SEGID-ADJACENT (segdup-family + marathon neighbors)",
    "S_REP": "REPEAT-HEAVY (top-decile census_mult/reads_kept)",
    "S_EXP": "EXPANSION-GATE-ADJACENT",
    "S_BIND": "BIND-AMBIGUITY (the ledger-named 7479-class)",
}
# the strata's sizes in the runnable universe (the draw script's own census)
U_SIZE = {"S_ORD": 1064, "S_SEG": 217, "S_REP": 34, "S_EXP": 4413, "S_BIND": 79}
UNIVERSE = 5807
# the committed exclusions (named in the doc; carried separately)
EXCLUDED = {"marathon_roster": 16, "poison": 5, "cut_exit_terminal": 16}


def load_tranche():
    rows = []
    with open(f"{WT}/tranche.tsv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            rows.append({"tid": a[0], "stratum": a[1]})
    return rows


def load_truth(profile, T):
    """per-locus records: expressible/rank/qv/called/truth, the A/B's own
    convention (rulings-ab-analysis.py VERBATIM, transposed to territory)."""
    path = f"{WT}/call{profile}/{T}/calls.jsonl.truth-qv.jsonl"
    if not os.path.exists(path):
        return None
    out = {}
    for line in open(path):
        r = json.loads(line)
        off = r.get("official_allwave") or {}
        out[r["locus"]] = {
            "expressible": r.get("truth_pair_expressible"),
            "rank": r.get("truth_rank"),
            "rank1": r.get("rank1"),
            "qv": off.get("qv"),
            "called": sorted((r["called_folds"][0]["index"],
                              r["called_folds"][1]["index"])),
            "truth": sorted((r["truth_folds"][0]["index"],
                             r["truth_folds"][1]["index"]))
                     if r.get("truth_folds") else None,
        }
    return out


PHASE_RE = re.compile(
    r"phase summary: inputs ([0-9.]+)s, quality \(FASTQ scan\) ([0-9.]+)s, "
    r"elsewhere ([0-9.]+)s, derive cache ([0-9.]+)s, binding ([0-9.]+)s, "
    r"variants ([0-9.]+)s, rows/context assembly ([0-9.]+)s, "
    r"loci \(folds[^)]*\) ([0-9.]+)s")


def load_phases(profile, T):
    """the serial call's own phase summary (RAYON=1: walls are true CPU)."""
    path = f"{WT}/logs/call{profile}-{T}.err"
    if not os.path.exists(path):
        return None
    for line in open(path, errors="replace"):
        if "phase summary:" in line:
            m = PHASE_RE.search(line)
            if m:
                return {"inputs": float(m.group(1)),
                        "quality": float(m.group(2)),
                        "binding": float(m.group(5)),
                        "variants": float(m.group(6)),
                        "rows": float(m.group(7)),
                        "loci": float(m.group(8))}
    return None


LOCUS_DONE_RE = re.compile(r"locus (\d+) done:.*\[([0-9.]+)s\]")


def load_locus_walls(profile, T):
    path = f"{WT}/logs/call{profile}-{T}.err"
    if not os.path.exists(path):
        return []
    return [float(m.group(2)) for m in
            (LOCUS_DONE_RE.search(l) for l in open(path, errors="replace"))
            if m]


def load_fiterm(profile, T):
    """the truth-QV instrumentation's own allwave seconds, if present."""
    path = f"{WT}/call{profile}/{T}/instrument-receipt.jsonl"
    return None  # (the QV walls are inside the phase 'receipts' class)


def farbin_share(T, cutdir="cutA"):
    """the far-bin scattered mass share of a territory's census (the
    prune's own geometry: window diameter > 1 across all occurrences)."""
    path = f"{WT}/{cutdir}/{T}/cut-census.jsonl"
    if not os.path.exists(path):
        return None
    total = pruned = 0
    nrec = npruned = 0
    for line in open(path):
        r = json.loads(line)
        mult = r["multiplicity"]
        lo, hi = 1 << 30, -1
        for occ in r["occurrences"]:
            for w in occ["partitions"]:
                if w < lo:
                    lo = w
                if w > hi:
                    hi = w
        total += mult
        nrec += 1
        if hi - lo > 1:
            pruned += mult
            npruned += 1
    if not total:
        return None
    return {"mult_total": total, "mult_pruned": pruned,
            "share": pruned / total, "records": nrec, "records_pruned": npruned}


def load_cutstats(cutdir, T):
    path = f"{WT}/{cutdir}/{T}/cut-stats.json"
    if not os.path.exists(path):
        return None
    s = json.load(open(path))
    return {"census_mult": s.get("census_mult"),
            "reads_kept": s.get("reads_kept"),
            "wall_s": s.get("wall_seconds")}


def qv_stats(vals):
    if not vals:
        return None
    return {"n": len(vals), "median": statistics.median(vals),
            "ge30": sum(1 for q in vals if q >= 30.0),
            "at_cap": sum(1 for q in vals if q >= 50.0),
            "min": min(vals)}


def fmt_qv(s):
    if not s:
        return "no QV"
    return (f"QV median {s['median']:.1f}, >=QV30 {s['ge30']}/{s['n']}, "
            f"at-cap {s['at_cap']}/{s['n']}, min {s['min']:.1f}")


def stratified_total(samples_by_stratum, key):
    """the stratified total estimate + normal 95% CI half-width."""
    total = 0.0
    var = 0.0
    per = {}
    for s, vals in samples_by_stratum.items():
        if len(vals) == 0:
            per[s] = None
            continue
        mean = statistics.mean(vals)
        per[s] = {"n": len(vals), "mean": mean,
                  "median": statistics.median(vals),
                  "p90": sorted(vals)[min(len(vals) - 1,
                                          int(0.9 * len(vals)))],
                  "max": max(vals), "min": min(vals)}
        if len(vals) > 1:
            v = statistics.variance(vals)
            # finite-population-corrected stratum variance of the total
            N = U_SIZE[s]
            var += (N * N) * v / len(vals) * (1 - len(vals) / N)
        total += U_SIZE[s] * mean
    ci = 1.96 * math.sqrt(var) if var else 0.0
    return total, ci, per


def main():
    out = open(OUT, "w")

    def w(s=""):
        print(s)
        out.write(s + "\n")

    tranche = load_tranche()
    strata = {s: [r["tid"] for r in tranche if r["stratum"] == s] for s in STRATA}

    # ---------------- (a) the wash verdict --------------------------
    w("=" * 78)
    w("THE WASH VERDICT AT SCALE (per stratum, per profile; the official")
    w("allwave QV machinery, cap 50 - the same assessment as every banked")
    w("receipt; every rate carries its own denominator)")
    w("=" * 78)
    losses = []      # the L2-17-class losses (B right, A wrong)
    gains = []       # (A right, B wrong)
    for s in STRATA:
        members = strata[s]
        w(f"\n--- {s}: {STRATUM_NAMES[s]} (n={len(members)} territories drawn "
          f"of {U_SIZE[s]} in the runnable universe)")
        for prof in ("A", "B"):
            n_loci = n_ex = rank1 = 0
            qvs = []
            done = 0
            for T in members:
                recs = load_truth(prof, T)
                if recs is None:
                    continue
                done += 1
                n_loci += len(recs)
                for r in recs.values():
                    if r["expressible"]:
                        n_ex += 1
                        rank1 += bool(r["rank1"])
                    if r["qv"] is not None:
                        qvs.append(r["qv"])
            st = qv_stats(qvs)
            w(f"  profile {prof}: {done} territories called, {n_loci} loci, "
              f"truth-expressible {n_ex}, rank-1 {rank1}/{n_ex}"
              + (f" ({100 * rank1 / n_ex:.0f}%)" if n_ex else "")
              + f", {fmt_qv(st)}")
        # the differing loci + the conserved-elsewhere class
        n_diff = n_a_right = n_b_right = n_equal_qv = 0
        for T in members:
            A, B = load_truth("A", T), load_truth("B", T)
            if A is None or B is None:
                continue
            for locus in sorted(A):
                if locus not in B:
                    continue
                a, b = A[locus], B[locus]
                if a["called"] != b["called"]:
                    n_diff += 1
                    a_right = a["called"] == a["truth"]
                    b_right = b["called"] == b["truth"]
                    if b_right and not a_right:
                        losses.append((s, T, locus, a["qv"], b["qv"]))
                        w(f"    {T} locus {locus}: LOSS (B calls the truth "
                          f"pair, A does not): A qv {a['qv']}, B qv {b['qv']}")
                    elif a_right and not b_right:
                        gains.append((s, T, locus, a["qv"], b["qv"]))
                        w(f"    {T} locus {locus}: GAIN (A right, B wrong): "
                          f"A qv {a['qv']}, B qv {b['qv']}")
                    elif a_right and b_right:
                        n_a_right += 1  # both right, different material
                    else:
                        if a["qv"] is not None and b["qv"] is not None:
                            if abs(a["qv"] - b["qv"]) > 1e-9:
                                w(f"    {T} locus {locus}: neither matches "
                                  f"truth; A qv {a['qv']:.1f} vs B qv "
                                  f"{b['qv']:.1f}")
                            else:
                                n_equal_qv += 1
        w(f"  differing loci: {n_diff} (B-right-A-wrong {len([l for l in losses if l[0]==s])}, "
          f"A-right-B-wrong {len([g for g in gains if g[0]==s])}, "
          f"sequence-equivalent/equal-QV {n_equal_qv})")
    w(f"\nTHE CONSERVED-ELSEWHERE RATE (the L2-17-class: B calls the truth")
    w(f"pair, A does not): {len(losses)} losses across the tranche")
    # per stratum, per N loci
    for s in STRATA:
        members = strata[s]
        nloci = 0
        for T in members:
            A = load_truth("A", T)
            if A:
                nloci += len(A)
        n = len([l for l in losses if l[0] == s])
        if nloci:
            w(f"  {s}: {n} loss{'es' if n != 1 else ''} per {nloci} loci "
              f"({100 * n / nloci:.2f}% of loci)")
    n_g = len(gains)
    w(f"  (the reverse class - A right where B is wrong: {n_g})")

    # ---------------- (b) the CPU-work variance ----------------------
    w("")
    w("=" * 78)
    w("THE CPU-WORK VARIANCE (serial RAYON=1 truth calls: walls are true")
    w("CPU seconds; the ex-load convention: the inputs phase (the panel")
    w("load) is excluded, the fleet's fork form amortizes it)")
    w("=" * 78)
    by_stratum_exload = {s: [] for s in STRATA}
    by_stratum_loci = {s: [] for s in STRATA}
    per_locus = {s: [] for s in STRATA}
    for s in STRATA:
        for T in strata[s]:
            ph = load_phases("A", T)
            recs = load_truth("A", T)
            if ph is None or recs is None:
                continue
            ex = ph["binding"] + ph["variants"] + ph["rows"] + \
                ph["loci"] + ph["quality"]
            n = len(recs)
            by_stratum_exload[s].append(ex)
            if n:
                by_stratum_loci[s].append(ph["loci"] / n)
            per_locus[s].extend(load_locus_walls("A", T))
    w("\nper-territory ex-load call CPU-s, by stratum:")
    for s in STRATA:
        v = by_stratum_exload[s]
        if v:
            w(f"  {s}: n={len(v)} median {statistics.median(v):.1f} "
              f"mean {statistics.mean(v):.1f} p90 "
              f"{sorted(v)[int(0.9*len(v))-1 if len(v)>1 else 0]:.1f} "
              f"max {max(v):.1f}")
    w("\nper-locus core-s (the loci phase / loci), by stratum:")
    for s in STRATA:
        v = by_stratum_loci[s]
        if v:
            w(f"  {s}: n={len(v)} median {statistics.median(v):.2f} "
              f"mean {statistics.mean(v):.2f} max {max(v):.2f}")
    w("\nthe per-locus wall distribution (serial, profile A):")
    for s in STRATA:
        v = per_locus[s]
        if v:
            v.sort()
            top5 = sum(v[int(0.95 * len(v)):]) / max(1, len(v))
            w(f"  {s}: n={len(v)} median {v[len(v)//2]:.2f}s "
              f"p90 {v[int(0.9*len(v))]:.2f}s max {v[-1]:.2f}s "
              f"worst-5% share {100*top5/max(1e-9,sum(v)):.1f}%")
    # the genome-wall estimate with bounds
    total, ci, per = stratified_total(by_stratum_exload, "exload")
    w(f"\nTHE GENOME-WALL ESTIMATE (stratified, the 5,807-territory runnable")
    w(f"universe; every stratum's own receipts; 95% normal CI):")
    w(f"  ex-load call CPU total: {total/1e6:.2f}e6 core-s "
      f"(+/- {ci/1e6:.2f}e6, {100*ci/total:.0f}% relative)")
    w(f"  at 256 cores (the box's own floor arithmetic): "
      f"{total/256/3600:.1f}h (+/- {ci/256/3600:.1f}h)")
    w("  EXCLUDED, carried as the committed bound instead: the marathon")
    w(f"  roster (16 members, the RV-frame p90 2,912s-deepest 11,084s "
      f"serial lanes), the 5 poison, the 16 cut.exit - their mass is NOT")
    w("  in the stratified estimate; the tranche's S_SEG receipts bound it.")
    for s in STRATA:
        if per.get(s):
            p = per[s]
            w(f"  {s}: mean {p['mean']:.1f}s x {U_SIZE[s]} territories "
              f"(n={p['n']}, median {p['median']:.1f}, p90 {p['p90']:.1f}, "
              f"max {p['max']:.1f})")

    # the cut side
    w("\nthe cut side (per-territory cut walls, contended; the census mult")
    w("ratio B/A per stratum):")
    for s in STRATA:
        wa, wb, ratios = [], [], []
        for T in strata[s]:
            a, b = load_cutstats("cutA", T), load_cutstats("cutB", T)
            if a and b and a["census_mult"] and b["census_mult"]:
                ratios.append(b["census_mult"] / a["census_mult"])
                wa.append(a["wall_s"])
                wb.append(b["wall_s"])
        if ratios:
            w(f"  {s}: n={len(ratios)} census-mult B/A median "
              f"{statistics.median(ratios):.2f}x (min {min(ratios):.2f}x "
              f"max {max(ratios):.2f}x); cutA wall median "
              f"{statistics.median(wa):.0f}s, cutB wall median "
              f"{statistics.median(wb):.0f}s")

    # ---------------- (c) the far-bin mass share ---------------------
    w("")
    w("=" * 78)
    w("THE FAR-BIN MASS SHARE BY STRATUM (the prune's exact class: census")
    w("records whose placements scatter into NON-ADJACENT windows - the")
    w("owner's routing principle, the --far-bin-prune's own geometry,")
    w("measured on the restricted (A) censuses - the product census)")
    w("=" * 78)
    fb = {s: [] for s in STRATA}
    for s in STRATA:
        for T in strata[s]:
            r = farbin_share(T)
            if r:
                fb[s].append(r)
    for s in STRATA:
        v = [r["share"] for r in fb[s]]
        if v:
            v.sort()
            w(f"  {s}: n={len(v)} territories, share median "
              f"{statistics.median(v)*100:.2f}% p90 {v[int(0.9*len(v))]*100:.2f}% "
              f"max {v[-1]*100:.2f}% zero-share "
              f"{sum(1 for x in v if x == 0)}/{len(v)}")
    total_share, ci_share, per_share = stratified_total(
        {s: [r["share"] * r["mult_total"] for s, lst in fb.items() for r in lst]
         for s in []}, "x")  # placeholder: not used; shares are per-territory
    # the honest aggregate: the mult-weighted share over the tranche
    tm = sum(r["mult_total"] for lst in fb.values() for r in lst)
    pm = sum(r["mult_pruned"] for lst in fb.values() for r in lst)
    if tm:
        w(f"  TRANSCHE AGGREGATE (mult-weighted): {pm}/{tm} = "
          f"{100*pm/tm:.2f}% of the restricted census mass is far-bin "
          f"scattered (the mass the prune removes)")
    out.close()
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
