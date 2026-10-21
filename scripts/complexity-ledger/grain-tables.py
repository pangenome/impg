#!/usr/bin/env python3
"""THE COMPLEXITY LEDGER's grain-curve tabulator: reads the grain
re-measure work dirs (L1-w*/, T03226-w*/) and emits the cost-vs-grain
table (per-stage walls, the serial CLI callC phase decomposition, the
per-window fold/unit shapes from the call's own stderr receipts, and
the truth sidecar's rank-1 counts) - the lever bought or refused with
accuracy receipts, mechanically. Receipt-side only."""
import sys, os, re, json, argparse

def read_wall(d, name):
    p = os.path.join(d, name + ".wall")
    if os.path.exists(p):
        try:
            return int(open(p).read().strip())
        except ValueError:
            return None
    return None

def phase_summary(d, call="s9-cli-deep"):
    p = os.path.join(d, call + ".err")
    if call == "callC":
        p = os.path.join(d, "callC.err")
    if not os.path.exists(p):
        return None
    m = re.search(r"phase summary: inputs ([\d.]+)s.*?binding ([\d.]+)s.*?rows/context assembly ([\d.]+)s, loci \(.*?\) ([\d.]+)s", open(p, errors="replace").read())
    if not m:
        return None
    return dict(inputs=float(m.group(1)), binding=float(m.group(2)),
                rows=float(m.group(3)), loci=float(m.group(4)))

def locus_shapes(d, call="s9-cli-deep"):
    p = os.path.join(d, call + ".err")
    if call == "callC":
        p = os.path.join(d, "callC.err")
    if not os.path.exists(p):
        return None
    txt = open(p, errors="replace").read()
    folds = [int(m) for m in re.findall(r"locus \d+ done: (\d+) folds", txt)]
    units = [int(m) for m in re.findall(r"locus \d+ done: \d+ folds, (\d+) units", txt)]
    classes = [int(m) for m in re.findall(r"locus \d+ done: \d+ folds, \d+ units, (\d+) classes", txt)]
    walls = [float(m) for m in re.findall(r"locus \d+ done: .*?\[([\d.]+)s\]", txt)]
    if not folds:
        return None
    return dict(n=len(folds), folds_max=max(folds), folds_med=sorted(folds)[len(folds)//2],
                units_max=max(units), classes_max=max(classes),
                loci_wall_sum=sum(walls) if walls else None,
                loci_wall_med=sorted(walls)[len(walls)//2] if walls else None)

def truth_stats(d, call="s9-cli-deep"):
    base = os.path.join(d, call if call != "callC" else "callC")
    side = os.path.join(base, "calls.jsonl.truth-qv.jsonl")
    if not os.path.exists(side):
        return None
    n = r1 = expr = 0
    for line in open(side):
        try:
            d0 = json.loads(line)
        except json.JSONDecodeError:
            continue
        n += 1
        rank = d0.get("truth_rank")
        if rank is not None:
            expr += 1
            if rank == 1:
                r1 += 1
    return dict(loci=n, rank1=r1, expr=expr)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("dirs", nargs="+")
    args = ap.parse_args()
    print("| dir | substrate walls (s1/s4/s7/s8 or prep/cutA) | census s | CLI serial wall s | phase: inputs/bind/rows/loci | windows | folds max/med | units max | rank1/loci |")
    print("|---|---|---|---|---|---|---|---|---|")
    for d in args.dirs:
        name = os.path.basename(d)
        if not os.path.isdir(d):
            print(f"| {name} | MISSING | | | | | | | |")
            continue
        if name.startswith("L1"):
            sub = "/".join(str(x) for x in [read_wall(d, "s1-partition"), read_wall(d, "s4-pggb")] if x)
            call = "s9-cli-deep"
        else:
            prep = [read_wall(d, "prep/beds"), read_wall(d, "prep/export"), read_wall(d, "prep/maps")]
            cut = read_wall(d, "cutA")
            sub = "prep " + "/".join(str(x) for x in prep if x is not None) + f" cutA {cut}"
            call = "callC"
        cens = read_wall(d, "s7-census") or (read_wall(d, "s8-derive-cache") or 0)
        cliw = read_wall(d, call)
        ph = phase_summary(d, call) or {}
        sh = locus_shapes(d, call) or {}
        tr = truth_stats(d, call) or {}
        print(f"| {name} | {sub} | {cens} | {cliw} | "
              f"{ph.get('inputs','')}/{ph.get('binding','')}/{ph.get('rows','')}/{ph.get('loci','')} | "
              f"{sh.get('n','')} | {sh.get('folds_max','')}/{sh.get('folds_med','')} | {sh.get('units_max','')} | "
              f"{tr.get('rank1','')}/{tr.get('loci','')} |")

if __name__ == "__main__":
    main()
