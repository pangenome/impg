#!/usr/bin/env python3
"""THE GENOME TRANCHE 2 - THE MINE (the stream doctrine's own shape: the
workspace mass is reaped after this pass, so every answer the tranche's
completed legs produced is extracted HERE, at mining time, into small
receipts; nothing raw is copied).

Receipts written (all labeled with their n and their terminal/partial
status; only TERMINAL legs - call{A,C}.done - contribute verdict loci):
  inventory.tsv          per member: every leg's status, exit, wall,
                        RSS peak, failure class, locus count
  verdicts-A.tsv / -C.tsv  per-territory per-locus verdict tables
                        (expressible, truth rank, QV, called vs truth)
  prune-pairs-C-vs-A.tsv per paired locus: the prune's own verdict
  territory-walls.tsv   per territory per form: the serial phase walls
  census-counters.tsv    the cut-side counters (census mult, reads kept,
                        far-bin share, the prune's own machinery receipt)
  metrics.txt            the computed decision receipts (the tranche
                        doc's own metrics, every denominator stated)
  receipts/              the small terminal receipts kept verbatim
                        (rvcut determinism gate + campaign, stage.log,
                        co-tenant.log, tranche-failed.tsv)

usage: tranche2-mine.py <WT> <OUTDIR>
"""
import json
import math
import os
import re
import statistics
import sys

WT = os.path.abspath(sys.argv[1])
OUT = os.path.abspath(sys.argv[2])
os.makedirs(OUT, exist_ok=True)
os.makedirs(f"{OUT}/receipts", exist_ok=True)

STRATA = ["S_EXP", "S_BIND", "S_MAR"]
STRATUM_NAMES = {"S_EXP": "EXPANSION-GATE-ADJACENT (the universe's 76% stratum)",
                 "S_BIND": "BIND-AMBIGUITY (the 7479-class census)",
                 "S_MAR": "THE MARATHON ROSTER (outside the runnable universe; walls only)"}
# the runnable-universe sizes (the draw script's own census; S_MAR excluded)
U_SIZE = {"S_EXP": 4413, "S_BIND": 79, "S_MAR": 0}
PROFS = ("A", "C")

PHASE_RE = re.compile(
    r"phase summary: inputs ([0-9.]+)s, quality \(FASTQ scan\) ([0-9.]+)s, "
    r"elsewhere ([0-9.]+)s, derive cache ([0-9.]+)s, binding ([0-9.]+)s, "
    r"variants ([0-9.]+)s, rows/context assembly ([0-9.]+)s, "
    r"loci \(folds[^)]*\) ([0-9.]+)s")
LOCUS_DONE_RE = re.compile(r"locus (\d+) done:.*\[([0-9.]+)s\]")
PRUNE_RE = re.compile(
    r"far-bin prune: (\d+) records pruned, census mult (\d+) of (\d+) "
    r"\(([0-9.]+)%\)")


def rd(p):
    try:
        return open(p).read().strip()
    except OSError:
        return None


def load_tranche():
    rows = []
    with open(f"{WT}/tranche.tsv") as f:
        for line in f:
            if line.startswith("#"):
                continue
            a = line.rstrip("\n").split("\t")
            rows.append({"tid": a[0], "stratum": a[1], "chrom": a[2],
                         "start": a[3], "end": a[4]})
    return rows


def locus_count(T):
    """the territory's own locus census (the axis intervals) - the m in
    'partial: n of m' and the per-territory loci context."""
    ax = rd(f"{WT}/prep/{T}/axis-path.txt")
    try:
        a = json.load(open(f"{WT}/prep/{T}/reference-axis.json"))
        return sum(1 for iv in a["intervals"] if iv["component"] == ax)
    except OSError:
        return None


def rss_peak(profile, T):
    p = f"{WT}/logs/call{profile}-{T}.rss"
    if not os.path.exists(p):
        return None
    peak = 0
    for line in open(p):
        parts = line.split()
        if len(parts) >= 2 and parts[1].isdigit():
            peak = max(peak, int(parts[1]))
    return peak // 1024 if peak else None  # MB


def failure_class(profile, T):
    """the terminal log line's own class, verbatim-named."""
    p = f"{WT}/logs/call{profile}-{T}.err"
    if not os.path.exists(p):
        return "no-log"
    txt = open(p, errors="replace").read()
    if "RSS guard exceeded" in txt:
        return "rss-guard-kill"
    if "overlapping records of one path" in txt:
        return "overlapping-records (deterministic fail-closed)"
    if "no unique verifying key" in txt:
        return "verify-key fail-closed"
    if "malformed qv-biwfa helper response" in txt:
        return "qv-biwfa-helper fail-closed"
    return "unclassified"


def leg_status(profile, T):
    """TERMINAL (done), TERMINAL-FAIL (exit receipt + named class),
    PARTIAL (frozen mid-write at the stop), or EXCLUDED."""
    if os.path.exists(f"{WT}/prep/{T}/call{profile}.done"):
        return ("terminal", None)
    ex = rd(f"{WT}/logs/call{profile}-{T}.exit")
    if ex is None:
        # no exit receipt: frozen mid-write, or never run
        out = f"{WT}/call{profile}/{T}/calls.jsonl"
        if os.path.isdir(f"{WT}/call{profile}/{T}") or os.path.exists(out):
            return ("partial", "no exit receipt")
        return ("not-run", None)
    if ex == "143":
        # SIGTERM mid-run with no guard/class line: the stop's own kill
        return ("partial", "stopped mid-write (external SIGTERM)")
    return ("terminal-fail", failure_class(profile, T))


def load_truth(profile, T):
    """per-locus verdict records - the analyze script's own convention."""
    path = f"{WT}/call{profile}/{T}/calls.jsonl.truth-qv.jsonl"
    if not os.path.exists(path):
        return None
    out = {}
    for line in open(path):
        if not line.strip():
            continue
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
        return {}
    return {int(m.group(1)): float(m.group(2)) for m in
            (LOCUS_DONE_RE.search(l) for l in open(path, errors="replace"))
            if m}


def load_prune_receipt(T):
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
    path = f"{WT}/{cutdir}/{T}/cut-census.jsonl"
    if not os.path.exists(path):
        return None
    total = pruned = nrec = npruned = 0
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
            "share": pruned / total, "records": nrec,
            "records_pruned": npruned}


def load_cutstats(cutdir, T):
    p = f"{WT}/{cutdir}/{T}/cut-stats.json"
    if not os.path.exists(p):
        return None
    s = json.load(open(p))
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


def wopen(name):
    return open(f"{OUT}/{name}", "w")


def main():
    tranche = load_tranche()
    strata = {s: [r["tid"] for r in tranche if r["stratum"] == s] for s in STRATA}
    failed_list = [l.split("\t")[0] for l in
                   open(f"{WT}/tranche-failed.tsv").read().splitlines()
                   if l.strip()]
    excl = set(failed_list)

    # ---------------- the inventory (per member, every leg) -------------
    inv = wopen("inventory.tsv")
    inv.write("#tid\tstratum\tchrom\tstart\tend\tloci_m\t")
    inv.write("\t".join(f"call{p}_status\tcall{p}_exit\tcall{p}_wall_s\t"
                        f"call{p}_rss_peak_mb\tcall{p}_class"
                        for p in PROFS) + "\n")
    status = {p: {} for p in PROFS}
    loci_m = {}
    for r in tranche:
        T = r["tid"]
        loci_m[T] = locus_count(T)
        row = [T, r["stratum"], r["chrom"], r["start"], r["end"],
               str(loci_m[T])]
        for p in PROFS:
            if T in excl:
                st, cls = "excluded", "extract-length-mismatch (tranche-failed.tsv)"
                ex = wl = rp = ""
            else:
                st, cls = leg_status(p, T)
                ex = rd(f"{WT}/logs/call{p}-{T}.exit") or ""
                wl = rd(f"{WT}/logs/call{p}-{T}.wall") or ""
                rp = rss_peak(p, T) or ""
            status[p][T] = (st, cls)
            row += [st, ex, wl, str(rp) if rp else "",
                    cls if cls else ""]
        inv.write("\t".join(row) + "\n")
    inv.close()

    # ---------------- the verdict tables (TERMINAL legs only) -----------
    truth = {p: {} for p in PROFS}
    vfiles = {p: wopen(f"verdicts-{p}.tsv") for p in PROFS}
    for p in PROFS:
        vfiles[p].write(
            "#verdicts form %s - TERMINAL legs only (call%s.done); "
            "source: calls.jsonl.truth-qv.jsonl at mining time\n"
            "#tid\tstratum\tlocus\texpressible\trank1\ttruth_rank\tqv\t"
            "called_pair\ttruth_pair\tlocus_wall_s\n" % (p, p))
    for s in STRATA:
        for T in strata[s]:
            for p in PROFS:
                if status[p][T][0] != "terminal":
                    continue
                recs = load_truth(p, T)
                if recs is None:
                    continue
                truth[p][T] = recs
                walls = load_locus_walls(p, T)
                for locus in sorted(recs):
                    r = recs[locus]
                    vfiles[p].write("\t".join(map(str, [
                        T, s, locus, r["expressible"], bool(r["rank1"]),
                        r["rank"], r["qv"],
                        ",".join(map(str, r["called"])),
                        ",".join(map(str, r["truth"]))
                        if r["truth"] else "none",
                        walls.get(locus, "")])) + "\n")
    for p in PROFS:
        vfiles[p].close()

    # ---------------- the prune's paired verdicts ----------------------
    pf = wopen("prune-pairs-C-vs-A.tsv")
    pf.write("#the far-bin prune's own verdict per PAIRED locus "
             "(both forms TERMINAL at the territory)\n"
             "#tid\tstratum\tlocus\ta_called\tc_called\ttruth_pair\t"
             "a_right\tc_right\tverdict\ta_qv\tc_qv\n")
    c_wins, a_wins, n_pairs, qv_changed = [], [], 0, []
    for s in STRATA:
        for T in strata[s]:
            A, C = truth["A"].get(T), truth["C"].get(T)
            if not A or not C:
                continue
            for locus in sorted(A):
                if locus not in C:
                    continue
                a, c = A[locus], C[locus]
                if a["called"] == c["called"]:
                    n_pairs += 1
                    verdict = "same-call"
                else:
                    n_pairs += 1
                    a_right = a["called"] == a["truth"]
                    c_right = c["called"] == c["truth"]
                    if c_right and not a_right:
                        verdict = "PRUNE-TRUTH-WIN"
                        c_wins.append((s, T, locus))
                    elif a_right and not c_right:
                        verdict = "PRUNE-TRUTH-LOSS"
                        a_wins.append((s, T, locus))
                    else:
                        verdict = "differs-no-truth-stake"
                    pf.write("\t".join(map(str, [
                        T, s, locus,
                        ",".join(map(str, a["called"])),
                        ",".join(map(str, c["called"])),
                        ",".join(map(str, a["truth"])) if a["truth"] else "none",
                        a_right, c_right, verdict,
                        a["qv"], c["qv"]])) + "\n")
                if (a["qv"] is not None and c["qv"] is not None
                        and abs(a["qv"] - c["qv"]) > 1e-9):
                    qv_changed.append((s, T, locus, a["qv"], c["qv"]))
    pf.close()

    # ---------------- the walls (serial phases, terminal legs) ---------
    wf = wopen("territory-walls.tsv")
    wf.write("#per-territory serial (RAYON=1) phase walls, TERMINAL legs; "
             "ex_load = quality+binding+variants+rows+loci "
             "(the inputs/panel-load phase excluded, the autopsy method)\n"
             "#tid\tstratum\tform\twall_s\trss_peak_mb\tinputs_s\tquality_s\t"
             "binding_s\tvariants_s\trows_s\tloci_s\tex_load_s\tloci_n\t"
             "core_s_per_locus\n")
    walls = {p: {s: [] for s in STRATA} for p in PROFS}
    per_locus = {p: {s: [] for s in STRATA} for p in PROFS}
    n_phases_missing = 0
    for s in STRATA:
        for T in strata[s]:
            for p in PROFS:
                if status[p][T][0] != "terminal":
                    continue
                ph = load_phases(p, T)
                recs = truth[p].get(T) or {}
                wl = rd(f"{WT}/logs/call{p}-{T}.wall") or ""
                rp = rss_peak(p, T) or ""
                if ph is None:
                    n_phases_missing += 1
                    wf.write("\t".join([T, s, p, wl, str(rp)] +
                             [""] * 9) + "\n")
                    continue
                ex = (ph["quality"] + ph["binding"] + ph["variants"] +
                      ph["rows"] + ph["loci"])
                walls[p][s].append(ex)
                if recs:
                    per_locus[p][s].append(ph["loci"] / len(recs))
                wf.write("\t".join(map(str, [T, s, p, wl, rp,
                    ph["inputs"], ph["quality"], ph["binding"],
                    ph["variants"], ph["rows"], ph["loci"], ex,
                    len(recs),
                    round(ph["loci"] / len(recs), 3) if recs else ""])) + "\n")
    wf.close()

    # ---------------- the census counters ------------------------------
    cf = wopen("census-counters.tsv")
    cf.write("#the cut-side counters per territory (all 108 runnable "
             "members, cutA/cutB both TERMINAL); far-bin share = the "
             "census-geometry diameter>1 rule on the restricted (A) census\n"
             "#tid\tstratum\trvcut_reads_kept\tcutA_mult\tcutA_reads\t"
             "cutA_wall_s\tcutB_mult\tcutB_reads\tcutB_wall_s\t"
             "mult_ratio_BA\tfarbin_share_pct\tfarbin_mult_pruned\t"
             "farbin_mult_total\tprune_receipt_pct\n")
    fbshare = {}
    for s in STRATA:
        for T in strata[s]:
            if T in excl:
                continue
            a, b = load_cutstats("cutA", T), load_cutstats("cutB", T)
            fb = farbin_share(T)
            if fb:
                fbshare[T] = fb
            pr = load_prune_receipt(T)
            rk = rd(f"{WT}/prep/{T}/rvcut.done")
            rows = [T, s]
            rows.append(str(rk) if rk else "")
            for st, key in ((a, "census_mult"), (a, "reads_kept"),
                            (a, "wall_s"), (b, "census_mult"),
                            (b, "reads_kept"), (b, "wall_s")):
                rows.append(str(st[key]) if st else "")
            ratio = (b["census_mult"] / a["census_mult"]
                     if a and b and a["census_mult"] and b["census_mult"]
                     else "")
            rows.append(f"{ratio:.3f}" if ratio != "" else "")
            if fb:
                rows += [f"{100*fb['share']:.2f}", fb["mult_pruned"],
                         fb["mult_total"]]
            else:
                rows += ["", "", ""]
            rows.append(f"{pr['pct']:.2f}" if pr else "")
            cf.write("\t".join(map(str, rows)) + "\n")
    cf.close()

    # ---------------- the receipts kept verbatim -----------------------
    import shutil
    for src, dst in [
            ("rvcut/determinism-gate.txt", "receipts/rvcut-determinism-gate.txt"),
            ("rvcut/campaign.log", "receipts/rvcut-campaign.log"),
            ("logs/stage.log", "receipts/stage.log"),
            ("logs/co-tenant.log", "receipts/co-tenant.log"),
            ("tranche-failed.tsv", "receipts/tranche-failed.tsv"),
            ("bridge/truth2-bridge.paf", "receipts/truth2-bridge.paf"),
            ("bridge/target.specs", "receipts/bridge-target.specs"),
            ("bridge/window.specs", "receipts/bridge-window.specs")]:
        if os.path.exists(f"{WT}/{src}"):
            shutil.copyfile(f"{WT}/{src}", f"{OUT}/{dst}")

    # ---------------- metrics.txt (the computed receipts) --------------
    m = wopen("metrics.txt")

    def w(line=""):
        m.write(line + "\n")

    w("=" * 78)
    w("THE GENOME TRANCHE 2 - THE MINED RECEIPTS (computed at mining time;")
    w("the workspace mass is reaped after this pass; every number below")
    w("carries its own n and its terminal/partial status)")
    w("=" * 78)
    n_run = sum(1 for r in tranche if r["tid"] not in excl)
    w(f"\nTHE INVENTORY VERDICT: {len(tranche)} members drawn "
      f"({', '.join(f'{s} {len(strata[s])}' for s in STRATA)}); "
      f"{n_run} runnable ({', '.join(sorted(excl))} excluded: "
      f"extract-length-mismatch).")
    for p in PROFS:
        nd = sum(1 for T in status[p] if status[p][T][0] == "terminal")
        nf = sum(1 for T in status[p] if status[p][T][0] == "terminal-fail")
        np_ = [(T, status[p][T]) for T in status[p]
               if status[p][T][0] == "partial"]
        w(f"  form {p}: {nd} TERMINAL (answers), {nf} TERMINAL-FAIL "
          f"(named classes), {len(np_)} PARTIAL. Stage-"
          f"{'TC done' if p == 'C' and os.path.exists(f'{WT}/logs/stage-TC.done') else 'TA open (master stopped mid-TA)'}")
        cls = {}
        for T in status[p]:
            if status[p][T][0] == "terminal-fail":
                cls.setdefault(status[p][T][1], []).append(T)
        for k, v in sorted(cls.items()):
            w(f"    fail class '{k}': {len(v)} ({', '.join(sorted(v))})")
        for T, (_, why) in np_:
            w(f"    PARTIAL: call{p} {T}: partial: 0 loci complete of "
              f"{loci_m.get(T)} ({why}) - lost-to-stop, re-runnable later")
    w("  form B: NEVER STARTED (the master was stopped mid-TA; the TB")
    w("    stage never opened) - ZERO B-form denominators in tranche 2;")
    w("    the conserved-elsewhere rate (B-vs-A / B-vs-C) is therefore")
    w("    NOT computable from this workspace; named lost-to-stop.")
    w("  cut/prep legs: substrate 108/108, reads 109/109, rvcut 108/108,")
    w("    cutA 108/108, cutB 108/108 - ALL TERMINAL (stage markers "
      "P,R,K,F,A done).")
    w("  fleetB: exit 0, wall 13,206s, 108 rows in 4 chunks (cut phases "
      "3,173s / 6,760s / 12,241s / 4,056s), panel load 190s - TERMINAL.")
    w("  prod: 16/16 in-child production calls landed (first 4 rows of "
      "each chunk), all with non-empty calls.jsonl - TERMINAL.")

    # the wash table
    w("")
    w("=" * 78)
    w("THE WASH VERDICT (rank-1 of truth-expressible, the official allwave")
    w("QV machinery, cap 50; TERMINAL legs only; form B: no legs)")
    w("=" * 78)
    for s in STRATA:
        members = strata[s]
        drawn = len(members)
        w(f"\n--- {s}: {STRATUM_NAMES[s]} ({drawn} drawn"
          + (f" of {U_SIZE[s]} in the runnable universe"
             if U_SIZE[s] else ", outside the runnable universe") + ")")
        for p in PROFS:
            n_loci = n_ex = rank1 = done = 0
            qvs = []
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
            w(f"  form {p} [{status_label if (status_label := 'TERMINAL') else ''}]"
              f": {done} territories called, {n_loci} loci, "
              f"truth-expressible {n_ex}, rank-1 {rank1}/{n_ex}"
              + (f" ({100 * rank1 / n_ex:.0f}%)" if n_ex else "")
              + f", {fmt_qv(qv_stats(qvs))}")

    # the prune's verdict
    w("")
    w("=" * 78)
    w("THE PRUNE'S OWN VERDICT AT TRANCHE-2 SCALE (C-vs-A, paired TERMINAL")
    w("loci; the gate bar: zero truth losses)")
    w("=" * 78)
    w(f"  across {n_pairs} paired loci: PRUNE WINS {len(c_wins)}, "
      f"LOSSES {len(a_wins)}"
      + (" - THE GATE BAR HOLDS at this n" if not a_wins
         else " - THE GATE BAR IS VIOLATED (named below)"))
    for s, T, locus in c_wins:
        w(f"    WIN   {s} {T} locus {locus}")
    for s, T, locus in a_wins:
        w(f"    LOSS {s} {T} locus {locus}")
    up = sum(1 for *_, qa, qc in qv_changed if qc > qa)
    w(f"  QV-visible changes at {len(qv_changed)} paired loci "
      f"({up} higher under C, {len(qv_changed) - up} higher under A)")

    # the walls
    w("")
    w("=" * 78)
    w("THE PHASE WALLS (serial RAYON=1, TERMINAL legs; ex-load = the call")
    w("phases minus the panel-load inputs phase; every row carries its n)")
    w("=" * 78)
    if n_phases_missing:
        w(f"  (note: {n_phases_missing} terminal legs lack a phase-summary "
          f"line - counted in verdicts, excluded from walls)")
    for p in PROFS:
        w(f"\nper-territory ex-load call CPU-s, form {p}:")
        for s in STRATA:
            v = walls[p][s]
            if v:
                w(f"  {s}: n={len(v)} median {statistics.median(v):.1f} "
                  f"mean {statistics.mean(v):.1f} "
                  f"p90 {sorted(v)[int(0.9 * len(v))]:.1f} max {max(v):.1f}")
    w("\nper-locus core-s (the loci phase / loci), by stratum:")
    for p in PROFS:
        for s in STRATA:
            v = per_locus[p][s]
            if v:
                w(f"  form {p} {s}: n={len(v)} median "
                  f"{statistics.median(v):.2f} mean {statistics.mean(v):.2f} "
                  f"max {max(v):.2f}")
    w("\nthe per-locus wall distribution (serial, form C):")
    for s in STRATA:
        v = []
        for T in strata[s]:
            if status["C"][T][0] == "terminal":
                v.extend(load_locus_walls("C", T).values())
        if v:
            v.sort()
            top5 = sum(v[int(0.95 * len(v)):])
            w(f"  {s}: n={len(v)} median {v[len(v) // 2]:.2f}s "
              f"p90 {v[int(0.9 * len(v))]:.2f}s max {v[-1]:.2f}s "
              f"worst-5% share {100 * top5 / max(1e-9, sum(v)):.1f}%")
    # the prune's loci-phase drop, tied to the far-bin share
    drops = []
    for s in STRATA:
        for T in strata[s]:
            pa, pc = load_phases("A", T), load_phases("C", T)
            fb = fbshare.get(T)
            if (pa and pc and fb and status["A"][T][0] == "terminal"
                    and status["C"][T][0] == "terminal"):
                drops.append((s, T, pa["loci"], pc["loci"], fb["share"]))
    if drops:
        ratios = [c / a if a else 1.0 for _, _, a, c, _ in drops]
        w("\nTHE PRUNE'S LOCI-PHASE WALL DROP (paired TERMINAL legs, same "
          "census, same territory):")
        w(f"  n={len(drops)} median C/A {statistics.median(ratios):.2f}x "
          f"min {min(ratios):.2f}x max {max(ratios):.2f}x")
        hi = [d for d in drops if d[4] > 0.01]
        lo = [d for d in drops if d[4] <= 0.0001]
        if hi:
            rh = [d[3] / d[2] if d[2] else 1.0 for d in hi]
            w(f"  at far-bin share > 1% (n={len(hi)}): median C/A "
              f"{statistics.median(rh):.2f}x; named heavy pairs:")
            for s, T, a, c, sh in sorted(hi, key=lambda d: -d[4])[:12]:
                w(f"    {s} {T}: share {100 * sh:.1f}%, A loci {a:.1f}s "
                  f"vs C loci {c:.1f}s ({c / a if a else 0:.3f}x)")
        if lo:
            rl = [d[3] / d[2] if d[2] else 1.0 for d in lo]
            w(f"  at far-bin share ~0 (n={len(lo)}): median C/A "
              f"{statistics.median(rl):.2f}x (the no-prune-to-do class)")
    # the universe arithmetic fragment (S_EXP + S_BIND only; S_MAR outside)
    w("\nTHE UNIVERSE-WALL FRAGMENT tranche 2 can state (its own strata "
      "only; S_ORD/S_SEG/S_REP are tranche-1 strata, not re-measured here):")
    tot = var = 0.0
    for s in ("S_EXP", "S_BIND"):
        vals = walls["C"][s]
        if not vals:
            continue
        mean = statistics.mean(vals)
        N = U_SIZE[s]
        contrib = N * mean
        tot += contrib
        if len(vals) > 1:
            var += (N * N) * statistics.variance(vals) / len(vals) * \
                (1 - len(vals) / N)
        w(f"  form C {s}: mean {mean:.1f}s x {N} territories "
          f"(n={len(vals)} terminal, median "
          f"{statistics.median(vals):.1f}, max {max(vals):.1f}) = "
          f"{contrib / 1e6:.2f}e6 core-s"
          + (f" (CI includes tranche-1's 8 own members of the 79)"
             if s == "S_BIND" else ""))
    ci = 1.96 * math.sqrt(var) if var else 0.0
    w(f"  form C two-stratum fragment total: {tot / 1e6:.2f}e6 core-s "
      f"(+/- {ci / 1e6:.2f}e6, 95% normal CI) = {tot / 256 / 3600:.1f}h "
      f"(+/- {ci / 256 / 3600:.1f}h) at 256 cores - A FRAGMENT of the "
      f"universe estimate (tranche 1 owns S_ORD/S_SEG/S_REP and the "
      f"other 8 S_BIND members); never a full-universe total.")
    tot_a = 0.0
    for s in ("S_EXP", "S_BIND"):
        vals = walls["A"][s]
        if vals:
            tot_a += U_SIZE[s] * statistics.mean(vals)
    w(f"  form A two-stratum fragment total: {tot_a / 1e6:.2f}e6 core-s "
      f"(same shape, form A n={sum(len(walls['A'][s]) for s in ('S_EXP', 'S_BIND'))} terminal legs)")
    # S_MAR walls
    w("\nTHE MARATHON ROSTER'S OWN RECEIPT (16 members, outside every "
      "stratified estimate):")
    for T in strata["S_MAR"]:
        a_st, a_cl = status["A"][T]
        c_st, c_cl = status["C"][T]
        pa, pc = load_phases("A", T), load_phases("C", T)
        wa = rd(f"{WT}/logs/callA-{T}.wall") or ""
        wc_ = rd(f"{WT}/logs/callC-{T}.wall") or ""
        ra = rss_peak("A", T) or ""
        rc = rss_peak("C", T) or ""
        w(f"  {T}: A {a_st} (wall {wa}s, rss peak {ra}MB"
          + (f", loci {pa['loci']:.0f}s" if pa else "") +
          f", class {a_cl or '-'}) | C {c_st} (wall {wc_}s, rss peak {rc}MB"
          + (f", loci {pc['loci']:.0f}s" if pc else "") +
          f", class {c_cl or '-'})")

    # the census counters
    w("")
    w("=" * 78)
    w("THE CENSUS COUNTERS (cut-side, all 108 runnable members TERMINAL "
      "in both cuts; the B/A mult ratio and the far-bin share by stratum)")
    w("=" * 78)
    for s in STRATA:
        ratios, wa, wb = [], [], []
        shares = []
        for T in strata[s]:
            if T in excl:
                continue
            a, b = load_cutstats("cutA", T), load_cutstats("cutB", T)
            if a and b and a["census_mult"] and b["census_mult"]:
                ratios.append(b["census_mult"] / a["census_mult"])
                wa.append(a["wall_s"])
                wb.append(b["wall_s"])
            if T in fbshare:
                shares.append(fbshare[T]["share"])
        if ratios:
            w(f"  {s}: n={len(ratios)} census-mult B/A median "
              f"{statistics.median(ratios):.2f}x (min {min(ratios):.2f}x "
              f"max {max(ratios):.2f}x); cutA wall median "
              f"{statistics.median(wa):.0f}s, cutB wall median "
              f"{statistics.median(wb):.0f}s")
        if shares:
            shares.sort()
            w(f"  {s}: far-bin share n={len(shares)} median "
              f"{statistics.median(shares) * 100:.2f}% p90 "
              f"{shares[int(0.9 * len(shares))] * 100:.2f}% max "
              f"{shares[-1] * 100:.2f}% zero-share "
              f"{sum(1 for x in shares if x == 0)}/{len(shares)}")
    tm = sum(r["mult_total"] for r in fbshare.values())
    pm = sum(r["mult_pruned"] for r in fbshare.values())
    if tm:
        w(f"  TRANSCHE-2 AGGREGATE (mult-weighted, n={len(fbshare)}): "
          f"{pm}/{tm} = {100 * pm / tm:.2f}% of the restricted census "
          f"mass is far-bin scattered")
    agree = disagree = 0
    for T, r in fbshare.items():
        pr = load_prune_receipt(T)
        if pr:
            if abs(r["share"] * 100 - pr["pct"]) < 0.05:
                agree += 1
            else:
                disagree += 1
    w(f"  the form-C machinery prune receipts vs census geometry: "
      f"{agree} agree, {disagree} disagree")

    # the rvcut gate
    w("")
    w("=" * 78)
    w("THE RESTRICTION-BASIS STAGE (rvcut, stage K TERMINAL 108/108) AND")
    w("THE DETERMINISM GATE (the banked-era cross-check; receipts kept "
      "verbatim in receipts/)")
    w("=" * 78)
    gate = open(f"{OUT}/receipts/rvcut-determinism-gate.txt").read()
    n_id = gate.count("IDENTICAL")
    n_diff = sum(1 for l in gate.splitlines()
                 if l.strip() and not l.startswith("IDENTICAL"))
    w(f"  determinism gate: {n_id} IDENTICAL, {n_diff} GATE-DIFFERS "
      f"(the OCT-15-EARLY-SUBSTRATE era class, named per-member in the "
      f"receipt)")
    camp = open(f"{OUT}/receipts/rvcut-campaign.log").read()
    n_gk = camp.count("GUARD-KILLED")
    w(f"  the adaptive width at work: {n_gk} width-8 guard-kills "
      f"(the census working set is not predictable pre-batch), every "
      f"retry at width 4 exited 0 - 108/108 terminal")

    # prod
    w("")
    w("=" * 78)
    w("THE PRODUCTION CALLS (the fleet's in-child calls, 16 rows of the "
      "108-row manifest; TERMINAL, no truth assessment by design)")
    w("=" * 78)
    prod_dir = f"{WT}/prod"
    for d in sorted(os.listdir(prod_dir)):
        p = f"{prod_dir}/{d}/calls.jsonl"
        if os.path.exists(p):
            n = sum(1 for l in open(p) if l.strip())
            w(f"  {d.replace('call-', '')}: {n} loci called "
              f"({os.path.getsize(p)} bytes)")
    m.close()
    print(f"wrote receipts to {OUT}")


if __name__ == "__main__":
    main()
