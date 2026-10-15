#!/usr/bin/env python3
"""territory-normalization-predict.py -- the receipt-side prediction of
STAGE 2, the extent-normalization design (the owner's approved lever,
measured at 307 loci by the causal autopsy).

THE DERIVED RULE (documented in full in the report of record and in
the docs section; this script is its MEASUREMENT side, receipt-side
only, no instrument changes):

  The locus's normalized comparison domain is the locus's TERRITORY
  EXTENT -- the window the locus is defined on, in the routing's own
  committed territory coordinates (the census's per-occurrence
  touched-window sets). A record is IN-TERRITORY evidence of locus L
  iff at least one of its occurrences has touched-window set EXACTLY
  [L]: the occurrence's anchored material lies wholly within L's
  territory (the census's own anchor-span territory-touch rule,
  candidate-independent, committed). OUT-OF-TERRITORY units (records
  that touch L but have no occurrence wholly inside L's territory)
  enter the likelihood at the generated-elsewhere branch E(read) for
  EVERY candidate identically -- the marginalization's candidate-
  independent branch -- so they contribute ZERO differential and
  cannot vote for a longer spanner.

Because the normalized convention scores every in-territory unit
exactly as the committed convention does and pushes only
out-of-territory units to the (candidate-independent) E branch, the
prediction is EXACT up to a per-locus constant: each class LL shifts
by sum(count * E over out-of-territory units) -- the same constant for
every class -- so every ranking, truth rank, log-gap (up to the
constant, which cancels in the gap too) and winner is predicted here
bit-faithfully from the committed receipts + the committed census.

Usage: territory-normalization-predict.py COMPONENT [COMPONENT ...]
  (--out FILE to write the per-locus JSONL receipt; default stdout
   for the table only)
"""
import json
import math
import os
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"

import numpy as np


def census_witnesses(component):
    """Stream the committed multi-matching census; per record id, the
    set of windows L for which some occurrence's touched-window set is
    EXACTLY [L] (the in-territory witness)."""
    path = f"{D}/cosine-diagnostic-scratch/{component}/cosine-multi-census.jsonl"
    witnesses = {}
    record = 0
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            r = d.get("record", record)
            single = set()
            for occ in d.get("occurrences", []):
                parts = occ.get("partitions", [])
                if len(parts) == 1:
                    single.add(parts[0])
            witnesses[r] = single
            record += 1
    return witnesses


def predict_component(component, out=None):
    witnesses = census_witnesses(component)
    receipt_path = f"{D}/realign-exhaustive-{component}.jsonl"
    ingredients_path = f"{D}/realign-exhaustive-{component}.jsonl.ingredients.jsonl"
    rows = []
    converted = 0
    regressed = 0
    with open(receipt_path) as freceipt, open(ingredients_path) as fingredients:
        for rline, iline in zip(freceipt, fingredients):
            d = json.loads(rline)
            ing = json.loads(iline)
            locus = d["locus"]
            assert ing["locus"] == locus, "receipt/ingredients locus mismatch"
            units = ing["units"]
            folds = ing["folds"]
            matrix = np.array(ing["ll_matrix"], dtype=float)
            n_folds = len(folds)
            n_units = len(units)
            assert matrix.shape == (n_folds, n_units)
            counts = np.array([u["count"] for u in units], dtype=float)
            # The in-territory mask (the derived rule's witness).
            in_mask = np.array(
                [locus in witnesses.get(u["record"], ()) for u in units], dtype=bool
            )
            out_units = int((~in_mask).sum())
            out_mass = float(counts[~in_mask].sum()) if out_units else 0.0
            in_counts = counts[in_mask]
            # The normalized class table: in-territory units only (the
            # out-of-territory units contribute the class-independent
            # constant count*E and cancel).
            class_lls = np.empty(n_folds * (n_folds + 1) // 2)
            idx = 0
            sub = matrix[:, in_mask]
            for j in range(n_folds):
                mixed = np.logaddexp(sub[j], sub[: j + 1]) - math.log(2.0)
                for i in range(j + 1):
                    class_lls[idx] = float(np.dot(in_counts, mixed[i]))
                    idx += 1
            if d["truth_folds"] is None:
                rows.append(
                    {
                        "component": component,
                        "locus": locus,
                        "expressible": False,
                        "units": n_units,
                        "units_out_of_territory": out_units,
                        "mass_out_of_territory": out_mass,
                    }
                )
                continue
            ti, tj = d["truth_folds"]
            truth_flat = tj * (tj + 1) // 2 + ti
            truth_ll = float(class_lls[truth_flat])
            rank = 1 + int((class_lls > truth_ll).sum())
            best = float(class_lls.max())
            winner_flat = int(np.argmax(class_lls))
            wi = winner_flat
            wj = 0
            # invert the flat index
            while wi > wj:
                wi -= wj + 1
                wj += 1
            gap = best - truth_ll
            before_rank = d["truth_rank"]
            before_gap = d["log_gap"]
            verdict = "unchanged"
            if before_rank == 1 and rank == 1:
                verdict = "HOLD"
            elif before_rank == 1:
                verdict = "REGRESS"
                regressed += 1
            elif rank == 1:
                verdict = "CONVERT"
                converted += 1
            rows.append(
                {
                    "component": component,
                    "locus": locus,
                    "expressible": True,
                    "units": n_units,
                    "units_out_of_territory": out_units,
                    "mass_out_of_territory": out_mass,
                    "before": {
                        "truth_rank": before_rank,
                        "log_gap": before_gap,
                        "winner": d["best_fold_indices"],
                    },
                    "after": {
                        "truth_rank": rank,
                        "log_gap": gap,
                        "winner": [wi, wj],
                    },
                    "verdict": verdict,
                }
            )
    return rows, converted, regressed


def main():
    args = [a for a in sys.argv[1:] if not a.startswith("--")]
    out = None
    if "--out" in sys.argv:
        out = sys.argv[sys.argv.index("--out") + 1]
    total_converted = 0
    total_regressed = 0
    all_rows = []
    for component in args:
        rows, converted, regressed = predict_component(component)
        total_converted += converted
        total_regressed += regressed
        all_rows.extend(rows)
        print(f"== {component}: predicted CONVERT {converted}, REGRESS {regressed}")
        for row in rows:
            if not row["expressible"]:
                print(
                    f"  L{row['locus']}: inexpressible "
                    f"(out-of-territory units {row['units_out_of_territory']}/{row['units']})"
                )
                continue
            b = row["before"]
            a = row["after"]
            print(
                f"  L{row['locus']}: rank {b['truth_rank']} -> {a['truth_rank']}  "
                f"gap {b['log_gap']:.1f} -> {a['log_gap']:.1f}  "
                f"out-territory {row['units_out_of_territory']}/{row['units']} units "
                f"({row['mass_out_of_territory']:.0f} mass)  {row['verdict']}"
            )
    print(f"== TOTAL: CONVERT {total_converted}, REGRESS {total_regressed}")
    if out:
        with open(out, "w") as f:
            for row in all_rows:
                f.write(json.dumps(row) + "\n")
        print(f"receipt rows -> {out}")


if __name__ == "__main__":
    main()
