#!/usr/bin/env python3
"""The pinned-UNION inner length of the EXISTING sparse-skip index.

`ClassScanFastIndex::class_pair_sum` (realign_score.rs:5927) walks the
union of the two rows' non-`E` (pinned) units, not all `n_units`. This
script measures, at one locus, how long that walk actually is, so the
unit-dedup's 842-term form can be compared against what the tree already
does rather than only against the committed dense walk.

usage: union_probe.py [line_index] [samples] [ingredients] [receipt]
"""
import json
import sys
import numpy as np

LINE = int(sys.argv[1]) if len(sys.argv) > 1 else 8
SAMPLES = int(sys.argv[2]) if len(sys.argv) > 2 else 3000
ING = sys.argv[3] if len(sys.argv) > 3 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl"
REC = sys.argv[4] if len(sys.argv) > 4 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl"

line = open(ING, "rb").readlines()[LINE]
i = line.index(b'"ll_matrix":')
j = line.index(b',"locus":')
sub = line[i + len(b'"ll_matrix":'):j]
folds = sub.count(b'[') - 1
M = np.fromstring(sub.replace(b'[', b' ').replace(b']', b' ')
                  .replace(b',', b' '), sep=' ').reshape(folds, -1)
rec = json.loads(open(REC, "rb").readlines()[LINE])
E = np.float64(rec["scoring"]["elsewhere_log_prob"])
ebits = E.view(np.uint64)
pinned = [np.nonzero(M[f].view(np.uint64) != ebits)[0] for f in range(folds)]
n_units = M.shape[1]
print(f"locus {rec['locus']}: folds={folds} units={n_units} E={float(E)!r}")
print(f"mean pinned per row: {np.mean([len(p) for p in pinned]):.1f}")

rng = np.random.default_rng(7)
un = np.empty(SAMPLES, dtype=np.int64)
for s in range(SAMPLES):
    a, b = rng.integers(folds), rng.integers(folds)
    un[s] = len(np.union1d(pinned[a], pinned[b]))
print(f"pinned-union over {SAMPLES} random pairs: mean {un.mean():.1f} "
      f"median {np.median(un):.0f} p95 {np.percentile(un, 95):.0f} "
      f"max {un.max()}")
print(f"sparse-skip inner shrink vs the committed dense walk: "
      f"{n_units / un.mean():.2f}x")
