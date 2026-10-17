#!/usr/bin/env python3
"""THE WHOLE-GENOME RUNTIME VERIFICATION - THE ROLLUP (the measured runtime
datasheet's generator; every number from the receipts/markers on disk, never
re-run): the stage walls, the per-territory cost distribution (the ordinary
66% vs the segdup tails, joined against the grid's own segdup fraction), the
failure-class census, the lane-packing reality (RSS peaks), the disk
high-water (the heartbeat's minimum df), and the end-to-end wall vs the
claims (projection 20.79min committed; per-territory cut ~387s serial; call
~23s; the wave ~19-21min at 230-256 lanes EXTRAPOLATED; serial-sum
~40-42min).
usage: rv-rollup.py <RV> [--json out.json]

CLOSE-OUT ADAPTATION (2026-10-17T12:43Z, the STOPPED-BY-OWNER partial-drain close -
receipted here per the close order): the rollup does NOT demand all-terminal state;
it emits the partial-drain close honestly from whatever markers exist (2,189 of 5,844
call-terminal at the owner's stop; the never-called 3,655 counted as exactly that;
batch-127's wrapper-side markers absent with its work complete - counted by its four
map/<T>/<T>/cut-stats.json receipts, never fabricated). The one measurement fix at
the close: the two split-rejection ledgers (rejected.tsv + rejected-incremental.tsv)
record the SAME windows (the incremental pass and the final reconciliation agree) -
deduped by territory so each rejection counts ONCE (3,252 transitive-expansion +
795 ambiguous-axis-bridge = 4,047). Verified idempotent at the close: consecutive
runs byte-identical (rollup-final.json / rollup-final.txt, heartbeat_lines 3,353
final after the heartbeat daemon was stopped at the close).
"""
import json, os, sys, glob, statistics as st

RV = sys.argv[1]
OUT = {}
def j(path, default=None):
    try:
        return json.load(open(path))
    except Exception:
        return default

# ---- the grid context
segdup = {}
excluded_bp = 0
for line in open(f"{RV}/../grid/territories.bed"):
    a = line.rstrip("\n").split("\t")
    if len(a) >= 5:
        segdup[a[0]] = float(a[4])
grid = j(f"{RV}/../grid/grid-stats.txt") or {}

# ---- stage walls
q_walls = [int(open(f).read()) for f in glob.glob(f"{RV}/query/shard-*.wall")]
OUT["query"] = {"shards": len(q_walls), "wall_s_sum": sum(q_walls),
                "wall_s_max": max(q_walls) if q_walls else None,
                "wall_s_at_width": max(q_walls) if q_walls else None}

prep_steps = {}
for step in ["extract", "syng", "beds", "export", "walkreq", "walks", "maps"]:
    vals = [int(open(f).read()) for f in glob.glob(f"{RV}/prep/*/{step}.wall")]
    prep_steps[step] = {"n": len(vals), "sum_s": sum(vals),
                        "median_s": st.median(vals) if vals else None,
                        "p90_s": sorted(vals)[int(0.9*len(vals))] if vals else None,
                        "max_s": max(vals) if vals else None}
OUT["prep"] = prep_steps

batch_walls = [int(open(f).read()) for f in glob.glob(f"{RV}/cuts/batch-*/batch.wall")]
batch_logs = glob.glob(f"{RV}/cuts/batch-*/batch.log")
cut_stats = {}
for f in glob.glob(f"{RV}/map/*/*/cut-stats.json"):
    d = j(f)
    if d:
        cut_stats[os.path.basename(os.path.dirname(os.path.dirname(f)))] = d
# also the reaped receipts
for f in glob.glob(f"{RV}/receipts/*/cut-stats.json"):
    d = j(f)
    if d:
        cut_stats[os.path.basename(os.path.dirname(f))] = d
OUT["cut"] = {
    "batches": len(batch_walls), "batch_wall_s_sum": sum(batch_walls),
    "batch_wall_s_median": st.median(batch_walls) if batch_walls else None,
    "territories_cut": len(cut_stats),
    "cut_wall_s_sum": sum(v.get("wall_seconds", 0) for v in cut_stats.values()),
    "scan_s_sum": sum(v.get("scan_seconds", 0) for v in cut_stats.values()),
    "load_s_sum": sum(v.get("load_seconds", 0) for v in cut_stats.values()),
    "reads_kept_sum": sum(v.get("reads_kept", 0) for v in cut_stats.values()),
}

# ---- the wave: per-territory costs, the cost distribution
calls = []
def n_int(path):
    try:
        return int(open(path).read().strip() or 0)
    except Exception:
        return None

for f in glob.glob(f"{RV}/receipts/*/call.wall"):
    T = os.path.basename(os.path.dirname(f))
    r = {"territory": T, "call_s": int(open(f).read()),
         "segdup": segdup.get(T, 0.0)}
    p = os.path.dirname(f)
    r["reads_s"] = n_int(f"{p}/reads.wall")
    r["gpartition_s"] = n_int(f"{p}/gpartition.wall")
    r["peak_rss_kb"] = n_int(f"{p}/call.peak_rss_kb")
    nloci = 0
    try:
        nloci = sum(1 for _ in open(f"{p}/calls.jsonl"))
    except Exception:
        pass
    r["loci"] = nloci
    calls.append(r)
for f in glob.glob(f"{RV}/prep/*/call.wall"):
    T = os.path.basename(os.path.dirname(f))
    if any(c["territory"] == T for c in calls):
        continue
    p = os.path.dirname(f)
    r = {"territory": T, "call_s": int(open(f).read()), "segdup": segdup.get(T, 0.0),
         "reads_s": None, "gpartition_s": None, "peak_rss_kb": None, "loci": None}
    for m, k in [("reads.wall", "reads_s"), ("gpartition.wall", "gpartition_s"),
                 ("call.peak_rss_kb", "peak_rss_kb")]:
        if os.path.exists(f"{p}/{m}"):
            r[k] = n_int(f"{p}/{m}")
    calls.append(r)

def dist(vals):
    if not vals:
        return None
    s = sorted(vals)
    return {"n": len(s), "median": st.median(s), "p10": s[int(0.1*len(s))],
            "p90": s[int(0.9*len(s))], "max": s[-1], "sum": sum(s)}

wave_calls = [c for c in calls if c["call_s"] is not None]
zero = [c for c in wave_calls if c["segdup"] == 0.0]
tail = [c for c in wave_calls if c["segdup"] > 0.02]
OUT["wave"] = {
    "territories_called": len(wave_calls),
    "call_wall_s": dist([c["call_s"] for c in wave_calls]),
    "call_wall_s_zero_segdup": dist([c["call_s"] for c in zero]),
    "call_wall_s_segdup_tail": dist([c["call_s"] for c in tail]),
    "reads_fetch_s": dist([c["reads_s"] for c in wave_calls if c["reads_s"] is not None]),
    "gpartition_s": dist([c["gpartition_s"] for c in wave_calls if c["gpartition_s"] is not None]),
    "peak_rss_kb": dist([c["peak_rss_kb"] for c in wave_calls if c["peak_rss_kb"]]),
    "loci_sum": sum(c["loci"] or 0 for c in wave_calls),
}

# ---- the failure-class census
rej = {}
rej_seen = set()  # the two ledgers record the same windows (the incremental
# pass and the final reconciliation agree); count each territory ONCE
exp_rows = exp_paths = 0
for rf in (f"{RV}/prep/rejected.tsv", f"{RV}/prep/rejected-incremental.tsv"):
    if not os.path.exists(rf):
        continue
    for line in open(rf):
        a = line.rstrip("\n").split("\t")
        if len(a) < 2 or a[0] in rej_seen:
            continue
        rej_seen.add(a[0])
        cls = a[1].split(" rows ")[0]
        rej[cls] = rej.get(cls, 0) + 1
        if cls == "transitive-expansion" and " rows " in a[1]:
            try:
                exp_rows += int(a[1].split(" rows ")[1].split(" paths ")[0])
                exp_paths += int(a[1].split(" paths ")[1])
            except Exception:
                pass
accepted = [l.strip() for l in open(f"{RV}/prep/accepted.txt")] if os.path.exists(f"{RV}/prep/accepted.txt") else []
# the per-territory split.rejected bp census lives in the per-territory split.rejected
# (only for territories the deduped ledgers above did not already count)
for f in glob.glob(f"{RV}/prep/*/split.rejected"):
    T = os.path.basename(os.path.dirname(f))
    if T in rej_seen:
        continue
    a = open(f, errors="replace").read().split()
    if len(a) >= 5 and a[0] == "transitive-expansion" and a[1] == "rows" and a[3] == "paths":
        try:
            exp_rows += int(a[2]); exp_paths += int(a[4])
        except Exception:
            pass
# the guard-kill class (the batch-width packaging events, never territory strikes)
guard_kills = {}
batch_exits = {}
for f in glob.glob(f"{RV}/cuts/batch-*/batch.exit"):
    st = open(f).read().strip()
    batch_exits[st] = batch_exits.get(st, 0) + 1
for line in open(f"{RV}/cuts/batch-failures.log") if os.path.exists(f"{RV}/cuts/batch-failures.log") else []:
    if "GUARD-KILLED at width" in line:
        w = line.split("at width ")[1].split(" ")[0]
        guard_kills[w] = guard_kills.get(w, 0) + 1
widths = {}
for f in glob.glob(f"{RV}/prep/*/cut.width"):
    w = open(f).read().strip()
    widths[w] = widths.get(w, 0) + 1
OUT["failures"] = {
    "split_rejected": rej,
    "split_accepted": len(accepted),
    "transitive_expansion_rows_sum": exp_rows,
    "transitive_expansion_paths_sum": exp_paths,
    "cut_poison": [line.rstrip("\n").split("\t")[0] for line in
                   open(f"{RV}/cuts/poison.tsv")] if os.path.exists(f"{RV}/cuts/poison.tsv") else [],
    "no_evidence": len(glob.glob(f"{RV}/prep/*/call.skipped-no-evidence")),
    "call_exit": len(glob.glob(f"{RV}/prep/*/call.exit")),
    "cut_exit": len(glob.glob(f"{RV}/prep/*/cut.exit")),
    "batch_exit_census": batch_exits,
    "guard_kill_widths": guard_kills,
    "adaptive_width_pending": widths,
}

# ---- the disk high-water + lane-packing from the heartbeat
hb = []
for line in open(f"{RV}/wave-logs/heartbeat.log") if os.path.exists(f"{RV}/wave-logs/heartbeat.log") else []:
    a = line.split()
    try:
        hb.append({"t": a[0], "df_gb": int(a[2][:-1]), "memavail_gib": int(a[4][:-3])})
    except Exception:
        pass
OUT["disk"] = {"min_df_gb": min((h["df_gb"] for h in hb), default=None),
               "min_memavail_gib": min((h["memavail_gib"] for h in hb), default=None),
               "heartbeat_lines": len(hb)}

print(json.dumps(OUT, indent=1))
if "--json" in sys.argv:
    open(sys.argv[sys.argv.index("--json")+1], "w").write(json.dumps(OUT, indent=1))
