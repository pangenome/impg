#!/usr/bin/env python3
"""THE FLANK-OVERLAP DELTA TABLES: score a --flank-overlap run against
its committed baseline (rank-1, PERFECT, the biWFA + official QV
distributions, the low-QV tail's fate, the changed loci named, the
called representative strains where calls are given). Assessment-side
only; reads committed receipts read-only.

usage: flank-overlap-delta.py <label> <baseline-truth-qv.jsonl> <new-truth-qv.jsonl> \
        [--baseline-official <jsonl>] [--baseline-calls <calls.jsonl>] \
        [--new-calls <calls.jsonl>]
"""
import json
import sys


def read_jsonl(path):
    out = []
    opener = open
    if path.endswith(".gz"):
        import gzip
        opener = gzip.open
    with opener(path, "rt") as fh:
        for line in fh:
            if line.strip():
                out.append(json.loads(line))
    return out


def official_of(rec):
    o = rec.get("official_allwave")
    return o.get("qv") if isinstance(o, dict) else None


def dist(qvs):
    qvs = sorted(q for q in qvs if q is not None)
    n = len(qvs)
    if not n:
        return "n=0"
    med = qvs[n // 2] if n % 2 else (qvs[n // 2 - 1] + qvs[n // 2]) / 2
    cap50 = sum(1 for q in qvs if q >= 49.999)
    return (f"n={n} median {med:6.2f} QV50 {cap50/n:6.1%} "
            f">=40 {sum(1 for q in qvs if q >= 40)/n:6.1%} "
            f">=30 {sum(1 for q in qvs if q >= 30)/n:6.1%} "
            f"<30 {sum(1 for q in qvs if q < 30)/n:6.1%} "
            f"<20 {sum(1 for q in qvs if q < 20)/n:6.1%} min {qvs[0]:6.2f}")


def strains_of_call(call):
    out = []
    for path in (call or {}).get("diplotype", {}).get("paths", []):
        names = sorted({m["path_name"].split("#")[0] for m in path.get("members", [])})
        out.append(",".join(names))
    return out


def main():
    label = sys.argv[1]
    base_path, new_path = sys.argv[2], sys.argv[3]
    opts = sys.argv[4:]
    base_official = new_official = None
    base_calls = new_calls = None
    it = iter(opts)
    for o in it:
        v = next(it)
        if o == "--baseline-official":
            base_official = v
        elif o == "--new-calls":
            new_calls = v
        elif o == "--baseline-calls":
            base_calls = v
    base = {r["locus"]: r for r in read_jsonl(base_path)}
    new = {r["locus"]: r for r in read_jsonl(new_path)}
    bo = {}
    if base_official:
        for r in read_jsonl(base_official):
            q = official_of(r)
            if q is not None:
                bo[r["locus"]] = q
    lines = [f"== {label} =="]
    common = sorted(set(base) & set(new))
    lines.append(f"windows: baseline {len(base)}, new {len(new)}, common {len(common)}")
    b_rank1 = sum(1 for l in common if base[l].get("rank1"))
    n_rank1 = sum(1 for l in common if new[l].get("rank1"))
    b_perf = sum(1 for l in common if base[l].get("perfect"))
    n_perf = sum(1 for l in common if new[l].get("perfect"))
    lines.append(f"rank-1: {b_rank1} -> {n_rank1}   PERFECT: {b_perf} -> {n_perf}")
    bq = [base[l]["qv"] for l in common]
    nq = [new[l]["qv"] for l in common]
    lines.append(f"biWFA QV  base: {dist(bq)}")
    lines.append(f"biWFA QV  new : {dist(nq)}")
    no = {l: official_of(new[l]) for l in common if official_of(new[l]) is not None}
    if bo:
        lines.append(f"official  base: {dist([bo[l] for l in common if l in bo])}")
    if no:
        lines.append(f"official  new : {dist([no[l] for l in common if l in no])}")
    # the low-QV tail's fate (official where available, else biWFA)
    def qv_b(l):
        return bo.get(l) if bo else None
    tail = [l for l in common if (bo[l] if bo else base[l]["qv"]) < 30]
    if tail:
        fates = []
        for l in tail:
            bqv = bo[l] if bo else base[l]["qv"]
            nqv = (no.get(l) if no else None)
            nqv = new[l]["qv"] if nqv is None else nqv
            def f1(q):
                return "None" if q is None else f"{q:.1f}"
            fates.append(f"L{l}: {f1(bqv)}->{f1(nqv)}"
                         f"{' RANK1' if new[l].get('rank1') else ''}")
        lines.append("low-QV tail (<30) fate: " + "; ".join(fates))
    # the changed loci
    changed = [l for l in common
               if base[l].get("truth_rank") != new[l].get("truth_rank")
               or (base[l]["qv"] is None) != (new[l]["qv"] is None)
               or (base[l]["qv"] is not None
                   and abs(base[l]["qv"] - new[l]["qv"]) > 1e-9)]
    ncalls, bcalls = None, None
    if new_calls:
        ncalls = {r["locus"]: r for r in read_jsonl(new_calls)}
    if base_calls:
        bcalls = {r["locus"]: r for r in read_jsonl(base_calls)}
    lines.append(f"changed loci: {len(changed)} of {len(common)}")
    for l in changed:
        def fmt(q):
            return "None" if q is None else f"{q:.2f}"
        line = (f"  L{l}: rank {base[l].get('truth_rank')} -> {new[l].get('truth_rank')}, "
                f"biWFA QV {fmt(base[l]['qv'])} -> {fmt(new[l]['qv'])}")
        if no or bo:
            b_o = bo.get(l)
            n_o = no.get(l)
            if b_o is not None or n_o is not None:
                line += f", official {b_o if b_o is not None else '-'} -> {n_o if n_o is not None else '-'}"
        if ncalls and bcalls and l in ncalls and l in bcalls:
            bs, ns = strains_of_call(bcalls[l]), strains_of_call(ncalls[l])
            if bs != ns:
                line += f", reps {bs} -> {ns}"
        lines.append(line)
    print("\n".join(lines))


if __name__ == "__main__":
    main()
