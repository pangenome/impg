#!/usr/bin/env python3
"""Measure whether a territory's GLOBAL node ids cluster in the panel's
node-id space — the knife edge the projection index lives on.

An inverted token->record index wins only if a territory's node ids fall
in FEW index buckets: then a query decodes a small slice instead of all
236M postings. Measured on synthetic shapes the difference is 900x
(contiguous band, 1 of 1,039 buckets) vs 0.8x (scattered, all 1,039).

Input: the `IMPG_DUMP_TERRITORY_NODES` dump (one sorted global node id
per line). Output: bucket coverage at several shifts, and the predicted
postings decoded per query, which is the number that decides the design.
"""
import sys
from collections import Counter

TOTAL_POSTINGS = 236_000_000  # 2.6 postings/record x 90.9M records (measured)
NODE_SPACE = 1_088_855_468    # the panel's node-id span (measured)


def load(path):
    with open(path) as handle:
        return [int(line) for line in handle if line.strip()]


def main(path):
    ids = load(path)
    n = len(ids)
    span = ids[-1] - ids[0]
    print(f"territory node set: {n:,} ids")
    print(f"  id range        {ids[0]:,} .. {ids[-1]:,}   (span {span:,})")
    print(f"  panel node span {NODE_SPACE:,}  -> this territory covers {100*n/NODE_SPACE:.3f}% of it")

    # contiguity: how many maximal runs of consecutive-ish ids?
    runs, prev = 1, ids[0]
    gap = max(1, span // 1000)
    for value in ids[1:]:
        if value - prev > gap:
            runs += 1
        prev = value
    print(f"  runs (gap > {gap:,})   {runs:,}  -> mean run {n/max(1,runs):,.0f} ids")

    print("\n  bucket coverage (a query decodes every touched bucket):")
    print("  shift   buckets_total   touched   %touched   predicted postings decoded   vs 236M scan")
    for shift in (20, 18, 16, 14, 12):
        mask = (1 << (64 - shift)) - 1
        touched = {((v & 0xFFFFFFFFFFFFFFFF) >> shift) & mask for v in ids}
        total_buckets = (NODE_SPACE >> shift) + 1
        frac = len(touched) / total_buckets
        postings = frac * TOTAL_POSTINGS
        print(
            f"  2^{shift:<4}  {total_buckets:>10,}   {len(touched):>7,}   {100*frac:6.1f}%   "
            f"{postings:>24,.0f}   {TOTAL_POSTINGS/max(1,postings):6.1f}x"
        )

    # the decision
    mask = (1 << (64 - 16)) - 1
    touched = {((v & 0xFFFFFFFFFFFFFFFF) >> 16) & mask for v in ids}
    frac = len(touched) / ((NODE_SPACE >> 16) + 1)
    if frac < 0.05:
        verdict = "CLUSTERED - the inverted index wins big; build it"
    elif frac < 0.30:
        verdict = "PARTIAL - index helps; pick the shift by this table"
    else:
        verdict = "SCATTERED - the token index does NOT beat the scan; the win must come from a different key"
    print(f"\n  VERDICT (2^16): {verdict}")


if __name__ == "__main__":
    main(sys.argv[1])
