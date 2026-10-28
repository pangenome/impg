#!/usr/bin/env python3
"""EXHAUSTIVE exactness census: for every locus whose (fold,unit) cell
count is small enough, run BOTH accumulations at EVERY unordered class
pair (i<=j, incl. the hom pairs):

  committed : total = 0.0 ; for u ascending: total += count_u*mix(a_u,b_u)
  reduced   : for the non-constant vector-groups in first-occurrence
              order: total += mass_g*mix(a_rep,b_rep) ; then total + S
              (S = the constant-vector scalar, one per locus)

and compare bits. Also compare the argmax pair of each form.
"""
import json
import math
import multiprocessing as mp
import re
import sys
import numpy as np

ING = sys.argv[2] if len(sys.argv) > 2 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl"
REC = sys.argv[3] if len(sys.argv) > 3 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl"
OUT = sys.argv[4] if len(sys.argv) > 4 else "census.json"
LINES = [l for l in open(ING, "rb").readlines()]
INF = float("inf")


def mix(a, b):
    if a == -INF and b == -INF:
        return -INF
    m = a if a >= b else b
    return m + math.log(0.5 * math.exp(a - m) + 0.5 * math.exp(b - m))


def bits(x):
    return np.float64(x).view(np.uint64).item()


def key(x):
    b = bits(x)
    return b | 0x8000000000000000 if (b >> 63) == 0 else (~b) & 0xFFFFFFFFFFFFFFFF


def ulp(a, b):
    return abs(key(a) - key(b))


REMITTED = None


def emitted():
    global REMITTED
    if REMITTED is None:
        REMITTED = {}
        for l in open(REC, "rb").readlines():
            r = json.loads(l)
            REMITTED[r["locus"]] = r["class_log_likelihoods"]
    return REMITTED


def work(arg):
    lineno, max_pairs = arg
    line = LINES[lineno]
    i = line.index(b'"ll_matrix":')
    j = line.index(b',"locus":')
    sub = line[i + len(b'"ll_matrix":'):j]
    folds = sub.count(b'[') - 1
    flat = np.fromstring(
        sub.replace(b'[', b' ').replace(b']', b' ').replace(b',', b' '), sep=' ')
    u = line.index(b'"units":[')
    counts = np.array([int(x) for x in re.findall(rb'"count":(\d+)', line[u:])],
                      dtype=np.int64)
    locus = int(re.search(rb',"locus":(\d+)', line).group(1))
    pairs = folds * (folds + 1) // 2
    if counts.size == 0 or pairs * counts.size > max_pairs:
        return {'locus': locus, 'line': lineno, 'skipped': True,
                'folds': int(folds), 'units': int(counts.size)}
    assert flat.size == folds * counts.size, (locus, flat.size)
    M = flat.reshape(folds, counts.size)

    Mt = np.ascontiguousarray(M.T)
    colview = Mt.view(np.dtype((np.void, folds * 8))).ravel()
    uniq, first_idx, inv, gsize = np.unique(
        colview, return_index=True, return_inverse=True, return_counts=True)
    reps_bits = Mt.view(np.uint64)[first_idx]
    is_const = reps_bits.min(axis=1) == reps_bits.max(axis=1)
    gmass = np.bincount(inv, weights=counts.astype(np.float64),
                        minlength=uniq.size)
    order = np.argsort(first_idx)
    const_g = order[is_const[order]]
    nonconst_idx = order[~is_const[order]].tolist()
    const_units = np.nonzero(np.isin(inv, const_g))[0] if const_g.size else \
        np.array([], dtype=np.int64)
    const_val = 0.0
    S = 0.0
    for g in const_g.tolist():
        v = np.array([reps_bits[g, 0]], dtype=np.uint64).view(np.float64)[0]
        const_val = v
        for uu in np.nonzero(inv == g)[0].tolist():
            S += float(counts[uu]) * float(v)
    nonconst_units = np.nonzero(~np.isin(inv, const_g))[0].tolist() \
        if const_g.size else list(range(counts.size))
    rep_cols = [int(first_idx[g]) for g in nonconst_idx]
    rep_mass = [float(gmass[g]) for g in nonconst_idx]
    cf = [float(c) for c in counts]
    nonconst_cols = nonconst_units

    em = emitted()[locus]
    rep_pairs = rep_mism = 0
    rep_maxulp = 0
    n_pairs = checks = mism = 0
    max_ulp = 0
    max_delta = 0.0
    best_c = (-INF, None)
    best_r = (-INF, None)
    for fi in range(folds):
        a = M[fi]
        for fj in range(fi, folds):
            b = M[fj]
            tc = 0.0
            for k in range(counts.size):
                tc += cf[k] * mix(a[k], b[k])
            tr = 0.0
            for g in range(len(rep_cols)):
                c = rep_cols[g]
                tr += rep_mass[g] * mix(a[c], b[c])
            tr += S
            n_pairs += 1
            flat_k = fj * (fj + 1) // 2 + fi
            ev = em[flat_k]
            if ev is not None:
                rep_pairs += 1
                if bits(ev) != bits(tc):
                    rep_mism += 1
                    rep_maxulp = max(rep_maxulp, ulp(ev, tc))
            if bits(tc) != bits(tr):
                mism += 1
                d = ulp(tc, tr)
                max_ulp = max(max_ulp, d)
                max_delta = max(max_delta, abs(tc - tr))
            if tc > best_c[0]:
                best_c = (tc, (fi, fj))
            if tr > best_r[0]:
                best_r = (tr, (fi, fj))
    return {'locus': locus, 'line': lineno, 'skipped': False,
            'folds': int(folds), 'units': int(counts.size),
            'distinct': int(uniq.size), 'const_units': int(const_units.size),
            'nonconst_units': len(nonconst_cols), 'nonconst_groups': len(rep_cols),
            'pairs': n_pairs, 'bit_mismatch': mism,
            'repl_pairs': rep_pairs, 'repl_mismatch': rep_mism,
            'repl_maxulp': rep_maxulp,
            'max_ulp': max_ulp, 'max_abs_delta': max_delta,
            'argmax_same': best_c[1] == best_r[1],
            'best_committed': repr(best_c[0]), 'best_reduced': repr(best_r[0]),
            'best_pair': list(best_c[1]) if best_c[1] else None,
            'S': repr(S), 'const_val': repr(const_val)}


if __name__ == '__main__':
    max_pairs = int(sys.argv[1]) if len(sys.argv) > 1 else 60_000_000
    args = [(n, max_pairs) for n in range(len(LINES))]
    with mp.Pool(min(len(args), 24)) as pool:
        out = pool.map(work, args)
    out.sort(key=lambda r: r['locus'])
    json.dump(out, open(OUT, 'w'), indent=1)
    print(f"{'loc':>3} {'F':>4} {'U':>6} {'N':>5} {'cU':>6} {'pairs':>7} "
          f"{'mism':>7} {'maxulp':>6} {'max|d|':>9} {'argmax':>6} "
          f"{'repl':>7} {'replmis':>7} {'replulp':>7}")
    tb = tt = 0
    for r in out:
        if r['skipped']:
            print(f"{r['locus']:>3} {r['folds']:>4} {r['units']:>6} "
                  f"{'-- skipped (too many pairs)':>40}")
            continue
        tb += r['bit_mismatch']
        tt += r['pairs']
        print(f"{r['locus']:>3} {r['folds']:>4} {r['units']:>6} {r['distinct']:>5} "
              f"{r['const_units']:>6} {r['pairs']:>7} {r['bit_mismatch']:>7} "
              f"{r['max_ulp']:>6} {r['max_abs_delta']:>9.3e} "
              f"{str(r['argmax_same']):>6} {r['repl_pairs']:>7} "
              f"{r['repl_mismatch']:>7} {r['repl_maxulp']:>7}")
    print(f"\nEXHAUSTIVE over the censused loci: {tb} bit-mismatched pairs "
          f"of {tt} ({100.0*tb/max(1,tt):.2f}%)")
    rp = sum(r.get('repl_pairs', 0) for r in out)
    rm = sum(r.get('repl_mismatch', 0) for r in out)
    print(f"REPLICATION GATE vs the emitted class_log_likelihoods: "
          f"{rp - rm}/{rp} pairs bit-identical ({rm} mismatches)")
