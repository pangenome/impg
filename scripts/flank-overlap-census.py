#!/usr/bin/env python3
"""THE FLANK-OVERLAP CENSUS (measure the object before designing the
lever; the owner's bet, stated then tested: many of the error/low-QV
windows are fragmentary-member windows).

For every window of the committed five-locality receipts (pilot +
L1-L4, both depths), classify each member row of the window's
partition map FULL-OVERLAP vs FRAGMENTARY at the two owner example
settings (--flank-overlap 10000 bp and --flank-overlap 10%), then
correlate fragmentary-member presence (in the window's candidate
universe, and in the CALLED pair itself) with the error/QV classes of
the committed baselines (rank-1, the low-QV tail, the rank-2+ misses,
the inexpressible loci). Assessment-side only; reads no committed
receipt for writing, writes only its own work directory.

Usage: flank-overlap-census.py <out-dir>
"""
import json
import math
import os
import sys

GEN = "/home/erikg/hprcv2/generic-sampling"
PILOT = "/home/erikg/hprcv2/hg002-locus-pilot"
PROJ = "/home/erikg/hprcv2/hg002-projection"
MASK = "/home/erikg/hprcv2/hg002-masking"

LOCALITIES = {
    "L1": {"maps": f"{GEN}/L1/locality-graphs",
           "component": "HG002#1#chr1:4632437-4883593"},
    "L2": {"maps": f"{GEN}/L2/locality-graphs",
           "component": "HG002#1#chr13:28850052-29100914"},
    "L3": {"maps": f"{GEN}/L3/locality-graphs",
           "component": "HG002#1#chr4:87197838-87447808"},
    "L4": {"maps": f"{GEN}/L4/locality-graphs",
           "component": "HG002#1#chr8:31206967-31457011"},
    "pilot": {"maps": f"{PILOT}/locality-graphs-hap1",
               "component": "HG002#1#chr1:148974789-149224773"},
}

# The committed baseline receipts per locality and depth.
def baselines(loc):
    if loc == "pilot":
        return {
            "deep": (f"{MASK}/pilot-baseline/cli/calls.jsonl",
                     f"{MASK}/pilot-baseline/cli/calls.jsonl.truth-qv.jsonl",
                     None),
            "35x": (None, None, None),
        }
    return {
        "deep": (f"{GEN}/{loc}/cli-product/calls.jsonl",
                 f"{GEN}/{loc}/cli-product/calls.jsonl.truth-qv.jsonl",
                 f"{GEN}/{loc}/cli-product-official-qv/calls.jsonl.truth-qv.jsonl"),
        "35x": (f"{PROJ}/cut-{loc}-cli/calls.jsonl",
                 f"{PROJ}/cut-{loc}-cli/calls.jsonl.truth-qv.jsonl",
                 f"/home/erikg/hprcv2/writeup-codification/chain-{loc}/calls.jsonl.truth-qv.jsonl"),
    }


def flank_bp(setting, window_len):
    if setting.endswith("%"):
        return math.ceil((float(setting[:-1]) / 100.0) * window_len)
    return int(setting)


def covers(member_start, member_end, axis_start, axis_end, setting):
    """The ruled PRESENCE semantics: the row's aligned hull must
    OVERLAP both flank intervals (the containment reading was measured
    and refuted - the hulls structurally stop short of the edges); the
    degenerate case (window shorter than N) reduces to simple presence
    in the window."""
    n = flank_bp(setting, axis_end - axis_start)
    if axis_end - axis_start < n:
        return member_start < axis_end and member_end > axis_start
    return (member_start < axis_start + n and member_end > axis_start
            and member_start < axis_end and member_end > axis_end - n)


def read_jsonl(path):
    if path is None or not os.path.exists(path):
        return []
    out = []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            out.append(json.loads(line))
    return out


def load_maps(maps_dir, component):
    """Per window id: (partition, axis (start,end), member rows); the
    window id is the map's own window<N> index (the locus id)."""
    windows = {}
    names = [n for n in os.listdir(maps_dir) if n.endswith(".gfa.map.json")]
    def order(name):
        return int(name.replace("window", "").replace(".gfa.map.json", ""))
    for name in sorted(names, key=order):
        m = json.load(open(os.path.join(maps_dir, name)))
        window_id = order(name)
        for row in m["members"]:
            if row["path_name"] == component:
                windows[window_id] = (m["partition"], (row["start"], row["end"]),
                                      m["members"])
                break
    return windows


def main():
    out_dir = sys.argv[1]
    os.makedirs(out_dir, exist_ok=True)
    settings = ["10000", "10%"]
    per_window_rows = []
    for loc, spec in sorted(LOCALITIES.items()):
        component = spec["component"]
        windows = load_maps(spec["maps"], component)
        depth_classes = {}
        for depth, (calls_path, qv_path, official_path) in baselines(loc).items():
            if qv_path is None:
                continue
            qv = {r["locus"]: r for r in read_jsonl(qv_path)}
            off = {}
            if official_path:
                for r in read_jsonl(official_path):
                    o = r.get("official_allwave")
                    if isinstance(o, dict):
                        off[r["locus"]] = o
            calls = {r["locus"]: r for r in read_jsonl(calls_path)} if calls_path else {}
            depth_classes[depth] = (qv, off, calls)
        for window_id in sorted(windows):
            partition, (axis_start, axis_end), members = windows[window_id]
            frag = {s: 0 for s in settings}
            frag_called = {s: None for s in settings}
            truth_row_frag = {s: None for s in settings}
            for row in members:
                if row["path_name"] == component:
                    continue  # the axis row is never excluded (the frame)
                for s in settings:
                    if not covers(row["start"], row["end"], axis_start, axis_end, s):
                        frag[s] += 1
                        if row["path_name"].split("#")[0:2] == ["HG002", "2"]:
                            truth_row_frag[s] = True
                    elif row["path_name"].split("#")[0:2] == ["HG002", "2"] \
                            and truth_row_frag[s] is None:
                        truth_row_frag[s] = False
            for depth, (qv, off, calls) in depth_classes.items():
                call = calls.get(window_id)
                if call is None:
                    continue
                fold_rows = []
                for path in call.get("diplotype", {}).get("paths", []):
                    fold_rows.extend(path.get("members", []))
                for s in settings:
                    frag_called[s] = int(all(
                        not covers(r["start"], r["end"], axis_start, axis_end, s)
                        for r in fold_rows)) if fold_rows else None
                rec = {"locality": loc, "depth": depth, "locus": window_id,
                       "partition": partition,
                       "window": [axis_start, axis_end],
                       "window_len": axis_end - axis_start,
                       "members": len(members),
                       "frag_10kb": frag["10000"], "frag_10pct": frag["10%"],
                       "truth_row_frag_10kb": truth_row_frag["10000"],
                       "truth_row_frag_10pct": truth_row_frag["10%"],
                       "called_frag_10kb": frag_called["10000"],
                       "called_frag_10pct": frag_called["10%"]}
                q = qv.get(window_id)
                if q is not None:
                    rec.update({"rank1": q.get("rank1"),
                                "truth_rank": q.get("truth_rank"),
                                "biwfa_qv": q.get("qv"),
                                "perfect": q.get("perfect"),
                                "expressible": q.get("truth_pair_expressible")})
                o = off.get(window_id)
                if o is not None:
                    rec["official_qv"] = o.get("qv")
                per_window_rows.append(rec)
    with open(os.path.join(out_dir, "census-windows.jsonl"), "w") as fh:
        for r in per_window_rows:
            fh.write(json.dumps(r) + "\n")

    # ---------------- the aggregate correlation tables
    def agg(rows, key):
        n = len(rows)
        return (n, sum(1 for r in rows if r[key] and r[key] > 0))

    lines = []
    for depth in ("deep", "35x"):
        rows = [r for r in per_window_rows if r["depth"] == depth]
        if not rows:
            continue
        lines.append(f"== THE CENSUS, {depth} (per-window fragmentary-member prevalence by call class) ==")
        classes = [
            ("rank-1 (truth_rank 1)", lambda r: r.get("rank1") is True),
            ("rank-1 + biWFA QV60 PERFECT", lambda r: r.get("rank1") is True and r.get("perfect") is True),
            ("miss (truth_rank > 1)", lambda r: r.get("truth_rank") not in (None, 1)),
            ("low-QV (official < 30)", lambda r: r.get("official_qv") is not None and r["official_qv"] < 30),
            ("low-QV (official < 20)", lambda r: r.get("official_qv") is not None and r["official_qv"] < 20),
            ("inexpressible", lambda r: r.get("expressible") is False),
            ("TRUTH ROW fragmentary @10%", lambda r: r.get("truth_row_frag_10pct") is True),
            ("TRUTH ROW full @10%", lambda r: r.get("truth_row_frag_10pct") is False),
        ]
        lines.append(f"{'class':34s} {'n':>4s} {'win>=1frag@10%':>14s} {'frac':>6s} "
                     f"{'win>=1frag@10kb':>15s} {'frac':>6s} {'called-frag@10%':>16s}")
        for name, pred in classes:
            sub = [r for r in rows if pred(r)]
            n, f10pct = agg(sub, "frag_10pct")
            _, f10kb = agg(sub, "frag_10kb")
            called_frag = sum(1 for r in sub if r.get("called_frag_10pct"))
            have_called = sum(1 for r in sub if r.get("called_frag_10pct") is not None)
            lines.append(f"{name:34s} {n:>4d} {f10pct:>14d} "
                         f"{(f10pct/n if n else 0):>6.1%} {f10kb:>15d} "
                         f"{(f10kb/n if n else 0):>6.1%} "
                         f"{(f'{called_frag}/{have_called}' if have_called else '-'):>16s}")
        # the overall universe prevalence
        tot = sum(r["members"] - 1 for r in rows)
        f10pct = sum(r["frag_10pct"] for r in rows)
        f10kb = sum(r["frag_10kb"] for r in rows)
        lines.append(f"universe total member rows (axis excluded): {tot}; "
                     f"fragmentary @10%: {f10pct} ({f10pct/tot:.1%}); "
                     f"@10kb: {f10kb} ({f10kb/tot:.1%})")
        lines.append("")
    # per-locality detail
    for loc in sorted(LOCALITIES):
        for depth in ("deep", "35x"):
            rows = [r for r in per_window_rows
                    if r["locality"] == loc and r["depth"] == depth]
            if not rows:
                continue
            lines.append(f"-- {loc} {depth}: {len(rows)} windows, "
                         f"windows with >=1 fragmentary member @10%: "
                         f"{sum(1 for r in rows if r['frag_10pct'])}, @10kb: "
                         f"{sum(1 for r in rows if r['frag_10kb'])}; "
                         f"called pairs fragmentary @10%: "
                         f"{sum(1 for r in rows if r.get('called_frag_10pct'))}")
            for r in rows:
                if ((r.get("truth_rank") not in (None, 1)) or
                        (r.get("official_qv") is not None and r["official_qv"] < 30) or
                        r.get("expressible") is False or
                        r.get("truth_row_frag_10pct") is True):
                    lines.append(
                        f"   locus {r['locus']:>2d} (win {r['window'][0]}-{r['window'][1]}, "
                        f"len {r['window_len']}): members {r['members']}, "
                        f"frag @10% {r['frag_10pct']}, @10kb {r['frag_10kb']}, "
                        f"truth-row-frag @10% {r.get('truth_row_frag_10pct')}, "
                        f"called-frag @10% {r.get('called_frag_10pct')}, "
                        f"rank {r.get('truth_rank')}, biWFA QV {r.get('biwfa_qv')}, "
                        f"official QV {r.get('official_qv')}, "
                        f"expressible {r.get('expressible')}")
        lines.append("")
    table = "\n".join(lines)
    with open(os.path.join(out_dir, "census-tables.txt"), "w") as fh:
        fh.write(table)
    print(table)


if __name__ == "__main__":
    main()
