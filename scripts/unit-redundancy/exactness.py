#!/usr/bin/env python3
"""Bit-exactness experiment on the DEEP locus (the ingredients line with
the most (fold, unit) cells).

Ground truth: instrument-receipt.jsonl's class_log_likelihoods for the
same locus = THE committed unit-order accumulation
`total += count_u * mix_logsumexp(ll_a[u], ll_b[u])` at every kept
(bounded-scan) pair.

Measured here, for the kept pairs:
  * REPLICATION gate: re-run the committed loop in pure python (libm
    exp/ln) and check it reproduces the emitted bits (proves the python
    mix == rust mix and that unit_counts were read correctly);
  * V2 dedup-only, V3 const-drop+dedup (the full reduction), and (on a
    subsample) V1 const-drop-only: compare bits against the emitted
    committed value; report mismatches and ulp distances.
"""
import json
import math
import re
import sys
import numpy as np

ING = sys.argv[1] if len(sys.argv) > 1 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl"
REC = sys.argv[2] if len(sys.argv) > 2 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl"
OUT = sys.argv[3] if len(sys.argv) > 3 else "exactness.json"
DEEP_LINE = 8
SAMPLE_COMMITTED = 400
SAMPLE_V1 = 4000

INF = float("inf")


def mix(a, b):
    if a == -INF and b == -INF:
        return -INF
    m = a if a >= b else b
    return m + math.log(0.5 * math.exp(a - m) + 0.5 * math.exp(b - m))


def bits(x):
    return np.float64(x).view(np.uint64).item()


def _key(x):
    b = bits(x)
    return b | 0x8000000000000000 if (b >> 63) == 0 else (~b) & 0xFFFFFFFFFFFFFFFF


def ulp(a, b):
    return abs(_key(a) - _key(b))


# ---------------------------------------------------------------- load
line = open(ING, "rb").readlines()[DEEP_LINE]
i = line.index(b'"ll_matrix":')
j = line.index(b',"locus":')
sub = line[i + len(b'"ll_matrix":'):j]
folds = sub.count(b'[') - 1
flat = np.fromstring(
    sub.replace(b'[', b' ').replace(b']', b' ').replace(b',', b' '), sep=' ')
u = line.index(b'"units":[')
counts_i = np.array([int(x) for x in
                     re.findall(rb'"count":(\d+)', line[u:])], dtype=np.int64)
assert flat.size == folds * counts_i.size, (flat.size, folds, counts_i.size)
M = flat.reshape(folds, counts_i.size)
print(f"locus {DEEP_LINE}: folds={folds} units={counts_i.size} "
      f"cells={folds*counts_i.size}", flush=True)

rec = json.loads(open(REC, "rb").readlines()[DEEP_LINE])
E = rec["scoring"]["elsewhere_log_prob"]
cll = rec["class_log_likelihoods"]
cfp = rec["class_fold_pairs"]
assert len(cll) == len(cfp) == folds * (folds + 1) // 2
kept = [k for k, v in enumerate(cll) if v is not None]
print(f"E={E!r} kept_pairs={len(kept)} of {len(cll)}", flush=True)

# ---------------------------------------------------------------- groups
Mt = np.ascontiguousarray(M.T)
colview = Mt.view(np.dtype((np.void, folds * 8))).ravel()
uniq, first_idx, inv, gsize = np.unique(
    colview, return_index=True, return_inverse=True, return_counts=True)
bitsM = Mt.view(np.uint64)
reps = bitsM[first_idx]
is_const = reps.min(axis=1) == reps.max(axis=1)
gmass = np.bincount(inv, weights=counts_i.astype(np.float64),
                    minlength=uniq.size)
order = np.argsort(first_idx)                    # group order = first occurrence
const_g = order[is_const[order]]
nonconst_g = order[~is_const[order]]
assert const_g.size == 1, const_g
assert int(gsize[const_g].sum()) == int((inv == const_g[0]).sum())
const_unit_mask = np.zeros(counts_i.size, dtype=bool)
const_unit_mask[inv == const_g[0]] = True
const_units_idx = np.nonzero(const_unit_mask)[0]
nonconst_units_idx = np.nonzero(~const_unit_mask)[0]
print(f"distinct={uniq.size} const_groups={const_g.size} "
      f"const_units={const_units_idx.size} nonconst_units={nonconst_units_idx.size}",
      flush=True)

const_val = np.array([reps[const_g[0], 0]], dtype=np.uint64).view(np.float64)[0]
assert const_val == E, (const_val, E)
# the scalar: sum of count_u * v over the constant units, ascending unit order
S = 0.0
for uu in const_units_idx.tolist():
    S += float(counts_i[uu]) * float(const_val)
print(f"const value bits==E bits: {bits(const_val)==bits(E)}  "
      f"mix(v,v)==v bits: {bits(mix(const_val, const_val))==bits(const_val)}  "
      f"S={S!r}", flush=True)
# mix(v,v)==v over every distinct group value
bad = 0
for g in range(uniq.size):
    v = np.array([reps[g, 0]], dtype=np.uint64).view(np.float64)[0]
    if bits(mix(v, v)) != bits(v):
        bad += 1
print(f"mix(v,v) bit-equal to v for all {uniq.size} distinct values: {bad==0} "
      f"(violations {bad})", flush=True)

# ---------------------------------------------------------------- loops
nf = folds
res = {'V1': [], 'V2': [], 'V3': [], 'REPL': []}


def committed(i, j):
    a, b = M[i], M[j]
    t = 0.0
    c = counts_i
    for k in range(c.size):
        t += float(c[k]) * mix(a[k], b[k])
    return t


def v1_const_drop(i, j):
    a, b = M[i], M[j]
    t = 0.0
    for k in nonconst_units_idx.tolist():
        t += float(counts_i[k]) * mix(a[k], b[k])
    return t + S


def v2_dedup(i, j):
    a, b = M[i], M[j]
    t = 0.0
    for g in order.tolist():
        t += float(gmass[g]) * mix(a[first_idx[g]], b[first_idx[g]])
    return t


def v3_both(i, j):
    a, b = M[i], M[j]
    t = 0.0
    for g in nonconst_g.tolist():
        t += float(gmass[g]) * mix(a[first_idx[g]], b[first_idx[g]])
    return t + S


rng = np.random.default_rng(20251027)
sample = sorted(rng.choice(len(kept), size=min(SAMPLE_COMMITTED, len(kept)),
                          replace=False).tolist())
sample_v1 = sorted(rng.choice(len(kept), size=min(SAMPLE_V1, len(kept)),
                             replace=False).tolist())

for n, k in enumerate(kept):
    i, j = cfp[k]
    got = v3_both(i, j)
    res['V3'].append((k, i, j, cll[k], got))
    if n % 4 == 0:
        res['V2'].append((k, i, j, cll[k], v2_dedup(i, j)))
    if n in sample_v1:
        res['V1'].append((k, i, j, cll[k], v1_const_drop(i, j)))
    if n in sample:
        res['REPL'].append((k, i, j, cll[k], committed(i, j)))
    if n % 10000 == 0:
        print(f"  pair {n}/{len(kept)}", flush=True)

print("\n=== results (committed value = the emitted class_log_likelihoods) ===")
summary = {}
for name in ('REPL', 'V1', 'V2', 'V3'):
    rows = res[name]
    mism = [(k, i, j, c, g) for (k, i, j, c, g) in rows if bits(c) != bits(g)]
    ulps = [ulp(c, g) for (k, i, j, c, g) in rows]
    deltas = [abs(c - g) for (k, i, j, c, g) in rows]
    summary[name] = {
        'pairs': len(rows),
        'bit_mismatches': len(mism),
        'max_ulp': max(ulps) if ulps else 0,
        'mean_ulp': float(np.mean(ulps)) if ulps else 0.0,
        'max_abs_delta': max(deltas) if deltas else 0.0,
        'examples': [[m[0], m[1], m[2], repr(m[3]), repr(m[4]), int(ulp(m[3], m[4]))]
                     for m in mism[:5]],
    }
    print(f"{name:5s} pairs={len(rows):6d} bit_mismatch={len(mism):6d} "
          f"max_ulp={max(ulps) if ulps else 0:6d} "
          f"max|delta|={max(deltas) if deltas else 0.0:.3e}")
    for e in summary[name]['examples']:
        print('    ', e)

# argmax stability over the kept pairs, V3 vs emitted
emit = [(k, cll[k]) for k in kept]
red = [(v[0], v[4]) for v in res['V3']]
assert len(red) == len(res['V3'])
be = max(emit, key=lambda kv: kv[1])
br = max(red, key=lambda kv: kv[1])
srt = sorted((v for _, v in emit), reverse=True)
print(f"\nMAP kept-pair: emitted flat={be[0]} {be[1]!r} | reduced flat={br[0]} "
      f"{br[1]!r} | top1-top2 gap = {srt[0]-srt[1]!r}")
top5e = sorted(emit, key=lambda kv: -kv[1])[:5]
redmap = dict(red)
top5r = sorted(red, key=lambda kv: -kv[1])[:5]
print("top5 emitted  :", [(k, repr(v)) for k, v in top5e])
print("top5 reduced  :", [(k, repr(v)) for k, v in top5r])
print("reduced @ top5 emitted flats:", [(k, repr(redmap[k])) for k, _ in top5e])
print("argmax unchanged (MAP flat):", be[0] == br[0])
gap_ulp = ulp(srt[0], srt[1]) if len(srt) > 1 else 0
print(f"top1-top2 committed gap: {srt[0]-srt[1]:.6e} abs = {gap_ulp} ulp; "
      f"max observed deviation = {summary['V3']['max_ulp']} ulp "
      f"(margin ratio {gap_ulp/max(1,summary['V3']['max_ulp']):.3g}x)")
json.dump({'locus_line': DEEP_LINE, 'folds': int(folds),
           'units': int(counts_i.size), 'distinct': int(uniq.size),
           'const_units': int(const_units_idx.size), 'const_value': repr(float(const_val)),
           'S': repr(S), 'kept_pairs': len(kept), 'summary': summary,
           'map': [int(be[0]), int(br[0]), repr(be[1]), repr(br[1]),
                   repr(srt[0] - srt[1])]},
          open(OUT, 'w'), indent=1)
print("wrote", OUT)
