#!/usr/bin/env python3
"""THE COMPLEXITY LEDGER's exponent fitter: log-log OLS regressions of the
per-phase work against the receipted per-locus/per-run shape counters,
with fit quality (R^2, n) named honestly for every fit. Inputs:
  --loci  mined-tranche-loci.tsv   (the serial tranche callC corpus, the
        current committed stack: 2,800+ per-locus rows, RAYON=1 walls =
        core-seconds, the exhaustive assessment scan)
  --runs  mined-tranche-runs.tsv   (the per-territory phase summaries)
  --autopsy  autopsy-census.txt    (the five-locality per-class serial
        timers, the SS10.6.14 census: per-locus folds_res/scoring/
        classes/qual/exact + anchor/placement/assertion, timer-inflated
        upper bounds scaled by the receipted 164.8/182.1 factor)
  --probe  class-bound-probe/*-final.jsonl (the eight deep-shape loci's
        fold/units/records/pinned anatomy)
Output: a text ledger of fits (slope(s), intercept, R^2, n) per phase,
the model formulas W = c * X^p * Y^q with MEASURED exponents, and the
median-vs-deep shares. Receipt-side only."""
import sys, os, re, json, argparse
import numpy as np

def ols(X, y, names):
    """X: (n,k) design in LOG space incl. intercept; y: log(t). Returns coef, r2, n."""
    X = np.asarray(X, float); y = np.asarray(y, float)
    keep = np.isfinite(X).all(axis=1) & np.isfinite(y) & (y > -30) & (X[:, 1] > 0)
    X, y = X[keep], y[keep]
    if len(y) < 8:
        return None
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    pred = X @ beta
    ss_res = ((y - pred) ** 2).sum()
    ss_tot = ((y - y.mean()) ** 2).sum()
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else float("nan")
    return beta, r2, len(y)

def fit_single(x, t, label, out):
    x = np.asarray(x, float); t = np.asarray(t, float)
    keep = (x > 0) & (t > 0) & np.isfinite(x) & np.isfinite(t)
    x, t = x[keep], t[keep]
    if len(t) < 8:
        out.append(f"  {label}: n={len(t)} TOO SMALL to fit")
        return None
    X = np.column_stack([np.ones_like(x), np.log(x)])
    res = ols(X, np.log(t), None)
    if res is None:
        out.append(f"  {label}: fit failed")
        return None
    (b0, p), r2, n = res
    c = np.exp(b0)
    out.append(f"  {label}: t = {c:.3g} * X^{p:.3f}   R^2={r2:.4f}  n={n}")
    return (c, p, r2, n)

def fit_biv(x1, x2, t, label, out):
    x1 = np.asarray(x1, float); x2 = np.asarray(x2, float); t = np.asarray(t, float)
    keep = (x1 > 0) & (x2 > 0) & (t > 0)
    x1, x2, t = x1[keep], x2[keep], t[keep]
    if len(t) < 8:
        return None
    X = np.column_stack([np.ones_like(x1), np.log(x1), np.log(x2)])
    res = ols(X, np.log(t), None)
    if res is None:
        return None
    (b0, p, q), r2, n = res
    out.append(f"  {label}: t = {np.exp(b0):.3g} * X1^{p:.3f} * X2^{q:.3f}   R^2={r2:.4f}  n={n}")
    return (np.exp(b0), p, q, r2, n)

def load_tsv(path):
    with open(path) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        rows = [dict(zip(header, line.rstrip("\n").split("\t"))) for line in fh if line.strip()]
    for r in rows:
        for k, v in list(r.items()):
            if k in ("territory", "form", "log"):
                continue
            try:
                r[k] = float(v) if v != "" else np.nan
            except ValueError:
                pass
    return rows

def col(rows, name):
    return np.array([r.get(name, np.nan) for r in rows], float)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--loci", required=True)
    ap.add_argument("--runs", required=True)
    ap.add_argument("--autopsy", default="")
    ap.add_argument("--probe", nargs="*", default=[])
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    out = []
    loci = load_tsv(args.loci)
    runs = load_tsv(args.runs)
    L = lambda name: col(loci, name)
    R = lambda name: col(runs, name)

    # derived counters
    folds = L("folds"); classes = L("classes"); units = L("units")
    members = L("members"); recs = L("records"); mpairs = L("matrix_pairs")
    nopin = L("no_pin"); wall = L("locus_wall_s")
    pinned = np.where(np.isfinite(mpairs) & np.isfinite(nopin), mpairs - nopin, np.nan)
    scan_model = classes * units          # the exhaustive unit-pair walk
    fold_pairs = folds * (folds + 1) / 2  # the structural quadratic (verify vs classes)
    resid = classes - fold_pairs

    out.append("== A. THE STRUCTURAL QUADRATIC RECEIPT (classes == folds(folds+1)/2?) ==")
    ok = np.isfinite(resid) & np.isfinite(classes)
    out.append(f"  classes == folds*(folds+1)/2 exactly at {int((np.abs(resid[ok])<1e-9).sum())}/{int(ok.sum())} loci"
               f" (max |resid| {np.nanmax(np.abs(resid[ok])) if ok.any() else float('nan'):.3g})")

    out.append("== B. PER-LOCUS SERIAL WHOLE-LOCUS WALL vs the candidate work terms (the tranche corpus, current stack, exhaustive scan) ==")
    fit_single(scan_model, wall, "whole-locus wall vs (classes x units) [the scan's unit-pair product]", out)
    fit_single(pinned, wall, "whole-locus wall vs pinned pairs [the pair loop]", out)
    fit_single(mpairs, wall, "whole-locus wall vs (folds x units) [the matrix]", out)
    fit_single(units * recs ** 0, wall, "whole-locus wall vs units", out)
    fit_single(recs, wall, "whole-locus wall vs records", out)
    fit_biv(scan_model, pinned, wall, "wall vs (classes*units) x (pinned)^q", out)
    fit_biv(fold_pairs, units, wall, "wall vs fold_pairs^p x units^q", out)

    out.append("== C. THE FOLD-ENUMERATION CLASS (members -> folds_pre -> folds; the merge pass) ==")
    fpre = L("folds_pre")
    fit_single(members, fpre, "folds_pre vs members", out)
    fit_single(members, folds, "folds vs members", out)
    fit_single(recs, folds, "folds vs records", out)
    eqc = L("equal_col_candidates"); mgroups = L("merge_groups"); adm = L("admitted_pairs")
    fit_single(fpre, eqc, "equal-column candidates vs folds_pre (the pairwise equal-column scan: quadratic?)", out)
    fit_single(fpre ** 2, eqc, "equal-column candidates vs folds_pre^2", out)
    fbs = L("folds_build_s")
    fit_single(members, fbs, "folds_build_s vs members", out)
    fit_single(fpre, fbs, "folds_build_s vs folds_pre", out)
    fit_single(recs, fbs, "folds_build_s vs records", out)
    fit_biv(members, recs, fbs, "folds_build_s vs members^p x records^q", out)
    mg_rate = np.nansum(adm) / max(np.nansum(eqc), 1)
    out.append(f"  merge-pass admission rate: admitted/candidates = {np.nansum(adm):.0f}/{np.nansum(eqc):.0f} = {mg_rate:.4f} over the corpus")
    out.append(f"  fold collapse: median folds_pre/folds = {np.nanmedian(fpre[np.isfinite(fpre)&(folds>0)]/folds[np.isfinite(fpre)&(folds>0)]):.4f};"
               f"  members->folds_pre median = {np.nanmedian(fpre[np.isfinite(fpre)&(members>0)]/members[np.isfinite(fpre)&(members>0)]):.4f}")

    out.append("== D. THE PER-CLASS AUTOPSY FITS (L1, serial, the SS10.6.14 census timers, scaled 0.905) ==")
    if args.autopsy and os.path.exists(args.autopsy):
        rows = []
        with open(args.autopsy) as fh:
            for line in fh:
                m = re.match(r"^\s*(\d+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+\|\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+([\d.]+)\s+\|\s+(\d+) x\s+([\d,]+) =\s+([\d,]+) \(\s*([\d,]+)\)\s+(\d+)", line)
                if m:
                    g = m.groups()
                    rows.append(dict(locus=int(g[0]), tot=float(g[1]), folds_res=float(g[2]), scoring=float(g[3]),
                                     classes=float(g[4]), qual=float(g[5]), exact=float(g[6]), anchor=float(g[7]),
                                     placement=float(g[8]), assertion=float(g[9]), marginal=float(g[10]),
                                     folds=int(g[11]), units=int(g[12].replace(",", "")),
                                     pairs=int(g[13].replace(",", "")), scored=int(g[14].replace(",", "")), recs=int(g[15])))
        A = lambda name: np.array([r[name] for r in rows], float)
        S = 164.8 / 182.1  # the receipted timer-overhead scale-out
        for name in ("tot", "folds_res", "scoring", "classes", "anchor", "placement", "assertion", "marginal"):
            A(name).__mul__(S)
        out.append(f"  (n={len(rows)} L1 loci; timers scaled by {S:.4f}; the assertion + per-variant anchor classes are the PRE-HOIST era)")
        cl = np.array([r["folds"] * (r["folds"] + 1) / 2 * r["units"] for r in rows], float)
        fit_single(cl, np.array([r["classes"] for r in rows]) * S, "class-scan timer vs (f(f+1)/2 x units) [the unit-pair walk]", out)
        sc = np.array([r["scored"] for r in rows], float)
        fit_single(sc, np.array([r["scoring"] for r in rows]) * S, "scoring timer vs scored/placed pairs", out)
        fit_single(np.array([r["pairs"] for r in rows], float), np.array([r["scoring"] for r in rows]) * S, "scoring timer vs all (fold,unit) pairs", out)
        an = np.array([r["folds"] * r["units"] for r in rows], float)
        fit_single(an, np.array([r["anchor"] for r in rows]) * S, "anchor timer vs (folds x units) [the pre-hoist per-pair class]", out)
        fit_single(sc, np.array([r["placement"] for r in rows]) * S, "placement timer vs scored pairs", out)
        frx = np.array([r["folds"] for r in rows], float)
        fit_single(np.array([r["recs"] for r in rows], float), np.array([r["folds_res"] for r in rows]) * S, "folds_res timer vs records", out)
        fit_biv(frx, np.array([r["recs"] for r in rows], float), np.array([r["folds_res"] for r in rows]) * S, "folds_res vs folds^p x records^q", out)
        # per-pair rates
        tot_scoring = np.sum([r["scoring"] for r in rows]) * S
        tot_scored = np.sum([r["scored"] for r in rows])
        tot_classes = np.sum([r["classes"] for r in rows]) * S
        tot_unitpairs = sum(r["folds"] * (r["folds"] + 1) / 2 * r["units"] for r in rows)
        out.append(f"  rates: scoring {tot_scoring/tot_scored*1e9:.2f} ns/pinned-pair; class scan {tot_classes/tot_unitpairs*1e9:.2f} ns/unit-pair")

    out.append("== E. THE PER-RUN FIXED CLASSES (the tranche phase summaries; the grain-lever's amortizable terms) ==")
    wins = R("windows"); inps = R("inputs_wall_s"); binds = R("bind_phase_s"); rows_s = R("rows_wall_s"); loci_s = R("loci_s")
    crecs = R("census_records")
    fit_single(wins, inps, "inputs wall vs windows [expect ~0: the panel load is per-call]", out)
    fit_single(crecs, binds, "binding phase vs census records", out)
    fit_single(R("bind_occurrences"), R("bind_verify_s"), "bind verify vs occurrences [the seek/probe class]", out)
    fit_single(R("distinct_rows"), rows_s, "rows/context vs distinct rows", out)
    fit_single(crecs, loci_s, "loci phase vs census records", out)
    fit_single(wins, loci_s, "loci phase vs windows", out)
    # medians and deep tails
    for name, arr in (("inputs_s", inps), ("binding_s", binds), ("rows_s", rows_s), ("loci_s", loci_s)):
        a = arr[np.isfinite(arr)]
        if len(a):
            out.append(f"  {name}: median {np.median(a):.1f}s  p90 {np.percentile(a,90):.1f}s  max {a.max():.1f}s  n={len(a)}")

    out.append("== F. THE MEDIAN-vs-DEEP LOCI (the tranche corpus shares) ==")
    for label, sel in (("median locus", np.nanpercentile(wall, 50)), ):
        pass
    sw = np.nanpercentile(wall, 50)
    dw = np.nanpercentile(wall, 95)
    tot = np.nansum(wall)
    ow = wall[np.isfinite(wall)]
    out.append(f"  per-locus serial wall: median {sw:.2f}s  p95 {dw:.2f}s  max {np.nanmax(ow):.1f}s  n={len(ow)}")
    sorted_w = np.sort(ow)[::-1]
    out.append(f"  worst-5% share of the corpus per-locus wall: {sorted_w[:max(1,len(sorted_w)//20)].sum()/tot*100:.1f}%")
    # scan-model share at median vs deep
    sm_ok = np.isfinite(scan_model) & np.isfinite(wall)
    out.append(f"  (classes x units) at the median locus: {np.nanmedian(scan_model):.3g}; at p95: {np.nanpercentile(scan_model,95):.3g}; max {np.nanmax(scan_model):.3g}")
    out.append(f"  folds: median {np.nanmedian(folds):.0f}  p95 {np.nanpercentile(folds,95):.0f}  max {np.nanmax(folds):.0f}")
    out.append(f"  units: median {np.nanmedian(units):.0f}  p95 {np.nanpercentile(units,95):.0f}  max {np.nanmax(units):.0f}")
    out.append(f"  records: median {np.nanmedian(recs):.0f}  p95 {np.nanpercentile(recs,95):.0f}  max {np.nanmax(recs):.0f}")

    if args.probe:
        out.append("== G. THE PROBE'S EIGHT DEEP SHAPES (the fold anatomy inputs) ==")
        for p in args.probe:
            for line in open(p):
                d = json.loads(line)
                if "folds" not in d or "locus" not in d:
                    continue
                out.append(f"  {d['member']} l{d['locus']}: folds {d['folds']} units {d['units']} pairs {d['pairs']}"
                           f" pinned {d['pinned']['pinned_share']:.2f} distinct_rows {d['memoization_exact_preserving']['distinct_rows']}"
                           f" scan_unit_pairs {d['counterfactual_bound']['scan_unit_pairs']:.3g}"
                           f" kept_prod {d['counterfactual_bound']['kept_production']}"
                           f" bound_walk_over_scan {d['counterfactual_bound']['bound_walk_over_scan']:.4f}")

    with open(args.out, "w") as fh:
        fh.write("\n".join(out) + "\n")
    print("\n".join(out))

if __name__ == "__main__":
    main()
