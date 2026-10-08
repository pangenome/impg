#!/usr/bin/env python3
"""A-vs-B telescoped-extend re-baseline study (supervisor decision point 3).

Compares a DEV IMPG_DEV_TELESCOPED_EXTEND run (A: per-pair constants +
O(overlap) corrections, regrouped f64 additions) against the production
fold-path run (B: bit-identical to the v19-lineage reference) field by
field and verifies:
  (i)   identical INTEGER accounting/count signatures for every selection;
  (ii)  all loss (f64) fields within 1e-9 absolute;
  (iii) identical selected routes, enumerating any differences with the
        competing margins (f64 pseudo-tie flips);
  (iv)  unchanged tie sets (count-equivalence-class integer fields).
Promotion of A to production is an owner decision on this report.
"""
import json
import re
import sys

LOSS_ABS_TOL = 1e-9

# Machine/environment fields (wall-clock stage timers, probes, RSS, sync
# timing) are not part of the study criteria (i)-(iv) and necessarily differ
# between two separate process runs; excluded from the comparison.
TIMING_RE = re.compile(r"seconds|elapsed|peak_rss|rss|timestamp|sync_|memo_")


def load(path):
    with open(path) as handle:
        return json.load(handle)


def classify(path):
    if TIMING_RE.search(path) or "/stage_probe" in path:
        return "timing"
    if any(
        tag in path
        for tag in ("loss", "score", "margin", "incumbent")
    ):
        return "loss"
    if any(
        tag in path
        for tag in (
            "work",
            "occupancy",
            "state_bytes",
            "dropped",
            "pruned",
            "queries",
            "resets",
            "windows",
            "classes",
            "collisions",
            "bytes",
            "count",
            "entries",
            "hashes",
            "backpointers",
            "cuts",
            "beam_width",
            "locus",
            "sources",
            "ranges",
            "status",
        )
    ):
        return "int-or-count"
    return "other"


def walk(x, y, path=""):
    if isinstance(x, dict):
        if not isinstance(y, dict):
            yield (path, "type", x, y)
            return
        for key in x:
            yield from walk(x[key], y.get(key), path + "/" + key)
        for key in y:
            if key not in x:
                yield (path + "/" + key, "only-in-B", None, y[key])
    elif isinstance(x, list):
        if not isinstance(y, list) or len(x) != len(y):
            yield (path, "len", x, y)
            return
        for index, (u, v) in enumerate(zip(x, y)):
            yield from walk(u, v, f"{path}[{index}]")
    else:
        if x != y:
            yield (path, "value", x, y)


def main():
    a_path, b_path = sys.argv[1], sys.argv[2]
    a, b = load(a_path), load(b_path)
    diffs = list(walk(a, b))
    int_bad, loss_ok, loss_bad, route_bad, other, timing = [], [], [], [], [], []
    for path, kind, x, y in diffs:
        cls = classify(path)
        if cls == "timing":
            timing.append((path, kind, x, y))
        elif cls == "loss":
            if isinstance(x, (int, float)) and isinstance(y, (int, float)):
                if abs(float(x) - float(y)) <= LOSS_ABS_TOL:
                    loss_ok.append((path, x, y))
                else:
                    loss_bad.append((path, x, y))
            else:
                other.append((path, kind, x, y))
        elif cls == "int-or-count":
            int_bad.append((path, x, y))
        elif "selected" in path or "choices" in path or "alleles" in path:
            route_bad.append((path, x, y))
        else:
            other.append((path, kind, x, y))
    print(f"total field diffs: {len(diffs)} (excluded timing/probe/RSS fields: {len(timing)})")
    print(f"(i) integer/count signature diffs: {len(int_bad)}")
    for item in int_bad[:20]:
        print("   ", item)
    print(f"(ii) loss diffs within {LOSS_ABS_TOL}: {len(loss_ok)}; VIOLATIONS: {len(loss_bad)}")
    for item in loss_bad[:20]:
        print("   ", item, "delta", float(item[1]) - float(item[2]))
    print(f"(iii) selected-route diffs: {len(route_bad)}")
    for item in route_bad[:20]:
        print("   ", item)
    print(f"other diffs: {len(other)}")
    for item in other[:20]:
        print("   ", item)
    verdict = not int_bad and not loss_bad and not route_bad and not other
    print("VERDICT:", "A matches B within the study criteria" if verdict else "CRITERIA VIOLATIONS PRESENT")
    return 0 if verdict else 1


if __name__ == "__main__":
    sys.exit(main())
