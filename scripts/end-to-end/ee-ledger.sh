#!/usr/bin/env bash
# THE END-TO-END RECEIPT - THE LEDGER MINER (receipts at write time,
# the stream law's own shape): scans the marker set and emits the
# per-leg TSV + the stage sums (walls, CORE-SECONDS from the .cpu
# receipts, the fail-closed census with classes, the guard census).
# usage: ee-ledger.sh   (stdout = the TSV; sums on stderr)
source /home/erikg/impg-genome-inference/scripts/end-to-end/ee-env.sh

cs() { awk '{print ($1+$2)}' "$1" 2>/dev/null || echo 0; }
wl() { cat "$1" 2>/dev/null || echo 0; }

# the per-leg TSV (one python pass; the .cpu receipts carry floats)
python3 - "$EE" <<'PYT'
import glob, json, os, re, sys
ee = sys.argv[1]
steps = ("extract","syng","beds","export","walkreq","walks","maps")
def num(path):
    # the .cpu receipts may carry the command's stderr before the
    # TIMEFORMAT line (the time-capture subshell shares stderr); the
    # LAST two-float line is the core-s receipt
    try:
        for line in reversed(open(path).read().splitlines()):
            a=line.split()
            if len(a)==2:
                try: return float(a[0])+float(a[1])
                except ValueError: continue
        return 0.0
    except Exception:
        return 0.0
def wall(path):
    try: return float(open(path).read().strip())
    except Exception: return 0.0
print("territory\tsubstrate_wall_s\tsubstrate_core_s\tpanelprep_wall_s\tpanelprep_core_s\treads_wall_s\tmap_wall_s\tmap_guard_events\tcall_done\tcalls_bytes\tloci\tchunk")
for t in open(f"{ee}/chunks/universe.order"):
    t=t.strip()
    if not t: continue
    p, l = f"{ee}/prep/{t}", f"{ee}/ledger/{t}"
    sw = sum(wall(f"{p}/{s}.wall") for s in steps)
    sc = sum(num(f"{p}/{s}.cpu") for s in steps)
    g = 0
    try:
        g = sum(1 for line in open(f"{l}/map.guard.txt") if "GUARD kill" in line)
    except Exception: pass
    cb = loci = 0; chunk = "-"
    try:
        for line in open(f"{l}/leg.json"):
            j = json.loads(line)
            cb += j.get("calls_bytes",0); loci += j.get("loci",0); chunk = j.get("chunk",chunk)
    except Exception: pass
    print(f"{t}\t{sw:.0f}\t{sc:.0f}\t{wall(f'{l}/panelprep.wall'):.0f}\t{num(f'{l}/panelprep.cpu'):.0f}\t"
          f"{wall(f'{l}/reads.wall'):.0f}\t{wall(f'{l}/map.wall'):.0f}\t{g}\t"
          f"{1 if os.path.exists(f'{l}/call.done') else 0}\t{cb}\t{loci}\t{chunk}")
PYT

# ---- the sums (stderr) -----------------------------------------------
python3 - "$EE" "$RV" 1>&2 <<'PY'
import glob, json, os, re, sys
ee, rv = sys.argv[1], sys.argv[2]
order = [t.strip() for t in open(f"{ee}/chunks/universe.order") if t.strip()]
def cpu(path):
    # the .cpu receipts may carry the command's stderr before the
    # TIMEFORMAT line (the time-capture subshell shares stderr); the
    # LAST two-float line is the core-s receipt
    try:
        for line in reversed(open(path).read().splitlines()):
            a = line.split()
            if len(a) == 2:
                try: return float(a[0]) + float(a[1])
                except ValueError: continue
        return 0.0
    except Exception:
        return 0.0
def wall(path):
    try: return float(open(path).read().strip())
    except Exception: return 0.0
steps = ("extract","syng","beds","export","walkreq","walks","maps")
sums = {s: [0.0,0.0] for s in steps}
sub_w = sub_c = pp_w = pp_c = m_w = m_c = r_w = 0.0
calls_done = calls_bytes = loci_done = 0
fail = {"substrate": [], "panelprep": [], "map": [], "call": []}
guards = []
for t in order:
    p, l = f"{ee}/prep/{t}", f"{ee}/ledger/{t}"
    for s in steps:
        sums[s][0] += wall(f"{p}/{s}.wall"); sums[s][1] += cpu(f"{p}/{s}.cpu")
    sub_w += sum(wall(f"{p}/{s}.wall") for s in steps)
    sub_c += sum(cpu(f"{p}/{s}.cpu") for s in steps)
    pp_w += wall(f"{l}/panelprep.wall"); pp_c += cpu(f"{l}/panelprep.cpu")
    m_w += wall(f"{l}/map.wall"); m_c += cpu(f"{l}/map.cpu")
    r_w += wall(f"{l}/reads.wall")
    if os.path.exists(f"{p}/substrate.exit") and not os.path.exists(f"{p}/substrate.done"): fail["substrate"].append(t)
    if os.path.exists(f"{p}/substrate.skipped"): fail["substrate"].append(t + "(skipped)")
    if os.path.exists(f"{l}/map.exit"): fail["map"].append(t)
    if os.path.exists(f"{l}/call.exit"): fail["call"].append(t)
    if os.path.exists(f"{l}/call.done"):
        calls_done += 1
        try:
            for line in open(f"{l}/leg.json"):
                j = json.loads(line)
                calls_bytes += j.get("calls_bytes", 0); loci_done += j.get("loci", 0)
        except Exception: pass
for f in glob.glob(f"{ee}/receipts/call-chunk-*.txt") + glob.glob(f"{ee}/ledger/*/map.guard.txt"):
    try:
        for line in open(f):
            if "GUARD kill" in line: guards.append(line.strip())
    except Exception: pass
n = len(order)
print(f"SUMS over {n} territories ({calls_done} called, {loci_done} loci, {calls_bytes} calls bytes):")
for s in steps:
    print(f"  substrate/{s}: wall {sums[s][0]:.0f}s core-s {sums[s][1]:.0f}")
print(f"  substrate total: wall {sub_w:.0f}s core-s {sub_c:.0f}")
print(f"  panel-prep: wall {pp_w:.0f}s core-s {pp_c:.0f}")
print(f"  reads fetch: wall {r_w:.0f}s (io-bound, core-s in reads.cpu)")
print(f"  map cut: wall {m_w:.0f}s core-s {m_c:.0f}")
print(f"  FAIL-CLOSED census: substrate {len(set(fail['substrate']))} {sorted(set(fail['substrate']))[:12]}")
print(f"                 map {len(set(fail['map']))} {sorted(set(fail['map']))[:12]}")
print(f"                 call {len(set(fail['call']))} {sorted(set(fail['call']))[:12]}")
print(f"  GUARD events: {len(guards)}")
for g in guards[:10]: print(f"    {g}")
# the call core-s: sum the chunk cpu receipts over completed chunks
cc = 0.0; cw = 0.0; nch = 0
for f in glob.glob(f"{ee}/ledger/chunk-*/call.cpu"):
    cc += cpu(f); nch += 1
for f in glob.glob(f"{ee}/ledger/chunk-*/call.wall"):
    cw += wall(f)
print(f"  call chunks: {nch} chunks wall {cw:.0f}s core-s {cc:.0f} (the chunk walls include the parent panel load + the child cut+call)")
PY
