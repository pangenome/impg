#!/usr/bin/env python3
"""THE GENOME TRANCHE - the analysis (the owner's decision-grade receipts;
three census/profile forms per the owner's far-bin ruling, commit 2a3dd46):
  (a) THE WASH VERDICT AT SCALE - per stratum, per form:
      A = restricted (the old default), B = global-native (the record
      form), C = restricted + far-bin-prune (THE NEW PRODUCT DEFAULT).
      Rank-1 rates, the official-allwave QV distribution (the likegt
      max-qv convention, cap 50), the differing loci named and
      truth-judged; THE CONSERVED-ELSEWHERE RATE - the L2-17-class
      losses (the record form calls the truth pair, the default does
      not) per stratum, under BOTH defaults (B-vs-A and B-vs-C - does
      the wash conclusion hold under the new default?).
  (a') THE PRUNE'S OWN VERDICT AT SCALE (the owner's product question):
      C-vs-A truth wins/losses (does the five-locality four-truth-win
      hold at sampled genome scale? the gate bar: zero truth losses),
      the QV changes, and the per-territory loci-phase wall drop C/A -
      the far-bin mass share BY STRATUM is the value receipt, tied to
      form C's own prune receipts (the machinery's pruned-records/mass
      log lines) and cross-checked against the census-geometry share.
  (b) THE CPU-WORK VARIANCE - per-territory serial (RAYON=1) call
      phases from the stage logs (the autopsy's serial-census method),
      per-locus core-s by stratum, the GENOME-WALL ESTIMATE WITH BOUNDS
      from the stratified mean/variance - AT THE NEW DEFAULT (form C)
      beside the old default (form A).
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
PROFS = ("A", "B", "C")
PROF_NAMES = {"A": "A restricted (the old default)",
              "B": "B global-native (the record form)",
              "C": "C restricted+far-bin-prune (THE NEW DEFAULT)"}
# the strata's sizes in the runnable universe (the draw script's own census)
U_SIZE = {"S_ORD": 1064, "S_SEG": 217, "S_REP": 34, "S_EXP": 4413, "S_BIND": 79}
UNIVERSE = 5807
EXCLUDED = {"marathon_roster": 16, "poison": 5, "cut_exit_terminal": 16}
PHASES_KEY = {"A": "phA", "C": "phC"}


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
    """per-locus records: the A/B analysis's own convention
    (rulings-ab-analysis.py VERBATIM, transposed to territory)."""
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
LOCUS_DONE_RE = re.compile(r"locus (\d+) done:.*\[([0-9.]+)s\]")
PRUNE_RE = re.compile(
    r"far-bin prune: (\d+) records pruned, census mult (\d+) of (\d+) "
    r"\(([0-9.]+)%\)")


def load_phases(profile, T):
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


def load_locus_walls(profile, T):
    path = f"{WT}/logs/call{profile}-{T}.err"
    if not os.path.exists(path):
        return []
    return [float(m.group(2)) for m in
            (LOCUS_DONE_RE.search(l) for l in open(path, errors="replace"))
            if m]


def load_prune_receipt(T):
    """form C's own machinery receipt: the pruned records/mass line."""
    path = f"{WT}/logs/callC-{T}.err"
    if not os.path.exists(path):
        return None
    for line in open(path, errors="replace"):
        m = PRUNE_RE.search(line)
        if m:
            return {"records": int(m.group(1)), "mult": int(m.group(2)),
                    "total": int(m.group(3)), "pct": float(m.group(4))}
    return None


def farbin_share(T, cutdir="cutA"):
    """the census-geometry share (the prune's exact class: records whose
    placements scatter into NON-ADJACENT windows, diameter > 1)."""
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


def stratified_total(samples_by_stratum):
    """the stratified total estimate + finite-population-corrected
    normal 95% CI half-width + per-stratum receipts."""
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
            N = U_SIZE[s]
            var += (N * N) * v / len(vals) * (1 - len(vals) / N)
        total += U_SIZE[s] * mean
    ci = 1.96 * math.sqrt(var) if var else 0.0
    return total, ci, per


def pairwise(recs, X, Y, T, out, label):
    """the differing-loci verdict between two forms at one territory;
    returns (losses, gains, n_diff, equal_qv) where losses = Y calls the
    truth pair and X does not (the conserved-elsewhere class for that
    direction)."""
    losses, gains, n_diff, n_eq = [], [], 0, 0
    for locus in sorted(X):
        if locus not in Y:
            continue
        x, y = X[locus], Y[locus]
        if x["called"] != y["called"]:
            n_diff += 1
            x_right = x["called"] == x["truth"]
            y_right = y["called"] == y["truth"]
            if y_right and not x_right:
                losses.append((T, locus))
                out.append(f"    {T} locus {locus} [{label}]: Y={y['called']} "
                           f"calls the TRUTH pair; X={x['called']} does not "
                           f"(X qv {x['qv']}, Y qv {y['qv']})")
            elif x_right and not y_right:
                gains.append((T, locus))
            elif x_right and y_right:
                pass
            elif x["qv"] is not None and y["qv"] is not None:
                if abs(x["qv"] - y["qv"]) > 1e-9:
                    out.append(f"    {T} locus {locus} [{label}]: neither "
                               f"matches truth; X qv {x['qv']:.1f} vs "
                               f"Y qv {y['qv']:.1f}")
                else:
                    n_eq += 1
    return losses, gains, n_diff, n_eq


def main():
    out = open(OUT, "w")

    def w(s=""):
        print(s)
        out.write(s + "\n")

    tranche = load_tranche()
    strata = {s: [r["tid"] for r in tranche if r["stratum"] == s] for s in STRATA}

    # cache the per-territory per-form truth records
    truth = {p: {} for p in PROFS}
    for p in PROFS:
        for s in STRATA:
            for T in strata[s]:
                r = load_truth(p, T)
                if r is not None:
                    truth[p][T] = r

    # ---------------- (a) the wash verdict, three forms ------------------
    w("=" * 78)
    w("THE WASH VERDICT AT SCALE - THREE FORMS (the owner's far-bin ruling")
    w("applied: A restricted [the old default], B global-native [the record")
    w("form], C restricted+far-bin-prune [THE NEW PRODUCT DEFAULT]; the")
    w("official allwave QV machinery, cap 50 - the same assessment as every")
    w("banked receipt; every rate carries its own denominator)")
    w("=" * 78)
    for s in STRATA:
        members = strata[s]
        w(f"\n--- {s}: {STRATUM_NAMES[s]} (n={len(members)} territories "
          f"drawn of {U_SIZE[s]} in the runnable universe)")
        for p in PROFS:
            n_loci = n_ex = rank1 = 0
            qvs = []
            done = 0
            for T in members:
                recs = truth[p].get(T)
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
            w(f"  form {PROF_NAMES[p]}: {done} territories called, "
              f"{n_loci} loci, truth-expressible {n_ex}, "
              f"rank-1 {rank1}/{n_ex}"
              + (f" ({100 * rank1 / n_ex:.0f}%)" if n_ex else "")
              + f", {fmt_qv(st)}")

    # ---------------- the conserved-elsewhere rate, both directions ------
    w("")
    w("=" * 78)
    w("THE CONSERVED-ELSEWHERE RATE (the L2-17 class: the RECORD form calls")
    w("the truth pair, the DEFAULT does not) - per stratum, per N loci,")
    w("under BOTH defaults: B-vs-A (the old default) and B-vs-C (the new)")
    w("=" * 78)
    named = []
    for s in STRATA:
        members = strata[s]
        n_loci = 0
        bva = []   # B right, A wrong
        bvc = []   # B right, C wrong
        for T in members:
            A, B, C = (truth["A"].get(T), truth["B"].get(T),
                       truth["C"].get(T))
            if A is None:
                continue
            n_loci += len(A)
            if B is not None:
                l, _, _, _ = pairwise(A, A, B, T, named, "B-vs-A")
                bva.extend(l)
            if C is not None and B is not None:
                # B-vs-C: the wash cost under the NEW default = B right
                # where C is wrong (X=C, Y=B: losses = B-right-C-wrong)
                rl, _, _, _ = pairwise(C, C, B, T, named, "B-vs-C")
                bvc.extend(rl)
        if n_loci:
            w(f"  {s}: {n_loci} loci assessed; B-right-A-wrong "
              f"{len(bva)} ({100*len(bva)/n_loci:.2f}%); B-right-C-wrong "
              f"{len(bvc)} ({100*len(bvc)/n_loci:.2f}% of loci)")
    for line in named[:40]:
        w(line)
    if len(named) > 40:
        w(f"    ... ({len(named) - 40} more named differing loci)")

    # ---------------- (a') the prune's own verdict at scale -------------
    w("")
    w("=" * 78)
    w("THE PRUNE'S OWN VERDICT AT SCALE (the owner's product question: does")
    w("the five-locality four-truth-win hold at sampled genome scale? the")
    w("gate bar: ZERO truth losses C-vs-A)")
    w("=" * 78)
    c_wins, a_wins, qv_changed = [], [], []
    n_pairs = 0
    for s in STRATA:
        for T in strata[s]:
            A, C = truth["A"].get(T), truth["C"].get(T)
            if A is None or C is None:
                continue
            for locus in sorted(A):
                if locus not in C:
                    continue
                n_pairs += 1
                a, c = A[locus], C[locus]
                a_right = a["called"] == a["truth"]
                c_right = c["called"] == c["truth"]
                if c_right and not a_right:
                    c_wins.append((s, T, locus))
                    w(f"    {s} {T} locus {locus}: PRUNE TRUTH WIN "
                      f"(C calls the truth pair, A does not)")
                elif a_right and not c_right:
                    a_wins.append((s, T, locus))
                    w(f"    {s} {T} locus {locus}: PRUNE TRUTH LOSS "
                      f"(A calls the truth pair, C does not)")
                if (a["qv"] is not None and c["qv"] is not None
                        and abs(a["qv"] - c["qv"]) > 1e-9):
                    qv_changed.append((s, T, locus, a["qv"], c["qv"]))
    w(f"  C-vs-A across {n_pairs} paired loci: PRUNE WINS "
      f"{len(c_wins)}, LOSSES {len(a_wins)}"
      + (" - THE GATE BAR (zero truth losses) HOLDS" if not a_wins
         else " - THE GATE BAR IS VIOLATED (named above)"))
    up = sum(1 for *_, qa, qc in qv_changed if qc > qa)
    w(f"  QV-visible changes at {len(qv_changed)} loci "
      f"({up} higher under C, {len(qv_changed)-up} higher under A)")

    # ---------------- (b) the CPU-work variance -------------------------
    w("")
    w("=" * 78)
    w("THE CPU-WORK VARIANCE (serial RAYON=1 truth calls: walls are true")
    w("CPU seconds; the ex-load convention: the inputs phase (the panel")
    w("load) is excluded - the fleet's fork form amortizes it)")
    w("=" * 78)
    by_form = {p: {s: [] for s in STRATA} for p in ("A", "C")}
    per_locus = {p: {s: [] for s in STRATA} for p in ("A", "C")}
    for p in ("A", "C"):
        for s in STRATA:
            for T in strata[s]:
                ph = load_phases(p, T)
                recs = truth[p].get(T)
                if ph is None or recs is None:
                    continue
                ex = ph["binding"] + ph["variants"] + ph["rows"] + \
                    ph["loci"] + ph["quality"]
                by_form[p][s].append(ex)
                if recs:
                    per_locus[p][s].append(ph["loci"] / len(recs))
    for p in ("A", "C"):
        w(f"\nper-territory ex-load call CPU-s, form "
          f"{PROF_NAMES[p]}:")
        for s in STRATA:
            v = by_form[p][s]
            if v:
                w(f"  {s}: n={len(v)} median {statistics.median(v):.1f} "
                  f"mean {statistics.mean(v):.1f} "
                  f"p90 {sorted(v)[int(0.9*len(v))]:.1f} max {max(v):.1f}")
    w("\nper-locus core-s (the loci phase / loci), by stratum:")
    for p in ("A", "C"):
        for s in STRATA:
            v = per_locus[p][s]
            if v:
                w(f"  form {p} {s}: n={len(v)} median "
                  f"{statistics.median(v):.2f} mean "
                  f"{statistics.mean(v):.2f} max {max(v):.2f}")
    # the prune's loci-phase drop, per territory, tied to the far-bin share
    drops = []
    for s in STRATA:
        for T in strata[s]:
            pa, pc = load_phases("A", T), load_phases("C", T)
            fb = farbin_share(T)
            if pa and pc and fb:
                drops.append((s, T, pa["loci"], pc["loci"], fb["share"]))
    if drops:
        ratios = [c / a if a else 1.0 for _, _, a, c, _ in drops]
        w("\nTHE PRUNE'S LOCI-PHASE WALL DROP (serial A vs C, same "
          "census, same territory):")
        w(f"  n={len(drops)} median C/A {statistics.median(ratios):.2f}x "
          f"min {min(ratios):.2f}x max {max(ratios):.2f}x")
        hi = [d for d in drops if d[4] > 0.01]
        lo = [d for d in drops if d[4] <= 0.0001]
        if hi:
            rh = [d[3] / d[2] if d[2] else 1.0 for d in hi]
            w(f"  at far-bin share > 1% (n={len(hi)}): median C/A "
              f"{statistics.median(rh):.2f}x")
        if lo:
            rl = [d[3] / d[2] if d[2] else 1.0 for d in lo]
            w(f"  at far-bin share ~0 (n={len(lo)}): median C/A "
              f"{statistics.median(rl):.2f}x (the no-prune-to-do class)")
    # the per-locus wall distribution at the new default
    w("\nthe per-locus wall distribution (serial, form C):")
    for s in STRATA:
        v = []
        for T in strata[s]:
            v.extend(load_locus_walls("C", T))
        if v:
            v.sort()
            top5 = sum(v[int(0.95 * len(v)):])
            w(f"  {s}: n={len(v)} median {v[len(v)//2]:.2f}s "
              f"p90 {v[int(0.9*len(v))]:.2f}s max {v[-1]:.2f}s "
              f"worst-5% share {100*top5/max(1e-9, sum(v)):.1f}%")
    # the genome-wall estimate with bounds, at BOTH defaults
    for p in ("C", "A"):
        total, ci, per = stratified_total(by_form[p])
        w(f"\nTHE GENOME-WALL ESTIMATE (stratified, the {UNIVERSE}-territory "
          f"runnable universe, form {PROF_NAMES[p]};")
        w("  95% normal CI, finite-population-corrected):")
        w(f"  ex-load call CPU total: {total/1e6:.2f}e6 core-s "
          f"(+/- {ci/1e6:.2f}e6, {100*ci/max(1e-9,total):.0f}% relative)")
        w(f"  at 256 cores (the box's own floor arithmetic): "
          f"{total/256/3600:.1f}h (+/- {ci/256/3600:.1f}h)")
        w("  EXCLUDED, carried as the committed bound instead: the marathon")
        w("  roster (16 members, the RV-frame p90 2,912s-deepest 11,084s")
        w("  serial lanes), the 5 poison, the 16 cut.exit - their mass is")
        w("  NOT in the stratified estimate; S_SEG's receipts bound it.")
        for s in STRATA:
            if per.get(s):
                q = per[s]
                w(f"  {s}: mean {q['mean']:.1f}s x {U_SIZE[s]} territories "
                  f"(n={q['n']}, median {q['median']:.1f}, p90 "
                  f"{q['p90']:.1f}, max {q['max']:.1f})")
    # the cut side
    w("\nthe cut side (per-territory cut walls, contended; the census")
    w("mult ratio B/A per stratum):")
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

    # ---------------- (c) the far-bin mass share, tied to form C --------
    w("")
    w("=" * 78)
    w("THE FAR-BIN MASS SHARE BY STRATUM (the prune's own value receipt:")
    w("the census-geometry share - records whose placements scatter into")
    w("NON-ADJACENT windows, the pure diameter>1 rule - measured on the")
    w("restricted (A) censuses, cross-checked against form C's own")
    w("machinery receipts: the pruned records/mass log lines)")
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
              f"{statistics.median(v)*100:.2f}% p90 "
              f"{v[int(0.9*len(v))]*100:.2f}% max {v[-1]*100:.2f}% "
              f"zero-share {sum(1 for x in v if x == 0)}/{len(v)}")
    tm = sum(r["mult_total"] for lst in fb.values() for r in lst)
    pm = sum(r["mult_pruned"] for lst in fb.values() for r in lst)
    if tm:
        w(f"  TRANSCHE AGGREGATE (mult-weighted): {pm}/{tm} = "
          f"{100*pm/tm:.2f}% of the restricted census mass is far-bin "
          f"scattered (the mass the prune removes at every territory)")
    # the machinery's own form-C receipts, cross-checked
    agree = disagree = 0
    for s in STRATA:
        for T in strata[s]:
            r = farbin_share(T)
            m = load_prune_receipt(T)
            if r and m:
                if abs(r["share"] * 100 - m["pct"]) < 0.05:
                    agree += 1
                else:
                    disagree += 1
                    w(f"    CROSS-CHECK DISAGREE {s} {T}: census geometry "
                      f"{r['share']*100:.3f}% vs form C receipt "
                      f"{m['pct']:.3f}%")
    if agree or disagree:
        w(f"  the form-C machinery receipts vs the census geometry: "
          f"{agree} agree, {disagree} disagree")
    out.close()
    print(f"\nwrote {OUT}")


if __name__ == "__main__":
    main()
