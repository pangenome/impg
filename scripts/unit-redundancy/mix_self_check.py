#!/usr/bin/env python3
"""mix(v,v) == v bit-for-bit over EVERY distinct ll value in the file.

The reduction "drop the constant-vector units" is only sound because the
pair term of a constant unit is exactly the constant: for any finite v,
mix_logsumexp(v,v) = v + ln(0.5*exp(0) + 0.5*exp(0)) = v + ln(1) = v + 0.0
= v, bit for bit (exp(0)=1 and ln(1)=+0.0 are exact in IEEE-754). This
script re-derives that by evaluation over the file's actual values.
"""
import json
import math
import re
import sys
import numpy as np

PATH = sys.argv[1] if len(sys.argv) > 1 else \
    "/tmp/probe-fast/T00048/instrument-receipt.jsonl.ingredients.jsonl"
INF = float("inf")


def mix(a, b):
    if a == -INF and b == -INF:
        return -INF
    m = a if a >= b else b
    return m + math.log(0.5 * math.exp(a - m) + 0.5 * math.exp(b - m))


total = viol = 0
with open(PATH, "rb") as f:
    for line in f:
        i = line.index(b'"ll_matrix":')
        j = line.index(b',"locus":')
        sub = line[i + len(b'"ll_matrix":'):j]
        flat = np.fromstring(sub.replace(b'[', b' ').replace(b']', b' ')
                             .replace(b',', b' '), sep=' ')
        if flat.size == 0:
            continue
        for v in np.unique(flat).tolist():
            total += 1
            if np.float64(mix(v, v)).view(np.uint64) != \
                    np.float64(v).view(np.uint64):
                viol += 1
                print(f"VIOLATION v={v!r}", file=sys.stderr)
        del flat, line
print(f"mix(v,v) bit-identical to v: {total - viol}/{total} distinct "
      f"values (violations {viol})")
assert viol == 0, "mix(v,v) != v for some value; the constant drop is unsound"
