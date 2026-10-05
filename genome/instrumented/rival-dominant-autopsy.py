#!/usr/bin/env python3
"""STAGE 4 — THE RIVAL-DOMINANT AUTOPSY under the normalized rule
(the owner's approved four, stage 4; receipt-side, assessment-only).

THE COHORT, re-derived under the territory convention (not inherited
from the pre-normalization classing): every expressible locus whose
truth is NOT called and whose BOTH-PLACED EVIDENCE GAP is positive —
the rival out-places the truth on the material BOTH pairs place (the
shared evidence). The unplaced-mass term may still favor the winner
or not; the defining property of this cohort is that NO mass
re-attribution repair can convert the locus (the prior autopsy's
"no mass-re-attribution lever" class).

PER LOCUS, WHY the rival out-places the truth on shared material —
named from the receipts (the territory scoring receipt joined with
the territory sequence-QV receipt; the mechanism features are all
derived, no thresholds-as-constants: every cut below is either an
exact structural equality or the stated convention):

  - twin_indistinguishable: identity >= 99.9% (the shared-evidence
    rival signature — the called and truth diplotype spell
    effectively the same sequence; the panel genuinely cannot
    distinguish them at this locus).
  - interior_substitution: gap_columns > mismatches (the rival's
    shared-evidence win is carried by interior gap material — the
    diverged-copy class: the reads match the rival's different-
    length interior better than the truth's).
  - substitution_better_match: mismatches >= gap_columns (the win
    is carried by base-substitution columns — the rival copy's
    bases genuinely match the reads better).
  - zygosity_dosage: the winner is homozygous while the truth pair
    is heterozygous, or the reverse (the mixture geometry itself:
    a pair of two identical copies explains a read population
    differently from two diverged homologs at equal per-unit
    scores).
  - unplaced_assist: the residual (one-side-placed) term is still
    winner-favored — the shared-evidence deficit is compounded by
    an unplaced asymmetry the normalization did not fully remove.

The primary mechanism is the first named in that order (the
strongest honest verdict first: indistinguishability, then the
material class, then geometry, then the compounding term).

THE FIX MAPPING, honest: no geometry (mass-re-attribution) lever
applies to this cohort BY DEFINITION (the shared evidence already
prefers the rival). The named levers that remain are model-side
and are stated with the measured constituency each would address;
the no-fix verdict is stated where the evidence genuinely prefers
the rival (what that means for the product: at those loci the
truth is not the maximum-likelihood answer under any comparison-
domain repair — the panel's actual copy ambiguity bounds the
call).

Usage: rival-dominant-autopsy.py [--components C1,C2,...]
Writes realign-rival-autopsy-<C>.jsonl per component and the
aggregate realign-rival-autopsy-tables.txt at the validation dir.
"""
import argparse
import json
import statistics
import sys
from pathlib import Path

D = Path("/home/erikg/yeast/genome-balanced-diploid-validation-20260930")
ALL_COMPONENTS = [
    "chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
    "chrVIII", "chrIX", "chrXI", "chrXIII", "chrXIV", "chrXVI",
    "chrVII", "chrX", "chrXII", "chrXV",
]


def say(msg):
    print(msg, flush=True)


def median(values):
    return statistics.median(values) if values else None


def load_qv(component):
    qv = {}
    path = D / f"realign-sequence-qv-territory-{component}.jsonl"
    if not path.exists():
        return qv
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            qv[d["locus"]] = d
    return qv


def load_prior(component):
    """The pre-normalization autopsy's per-locus records (the
    S1/S2 signatures for the prior-cause comparison; the census
    class is call-independent)."""
    prior = {}
    path = D / f"realign-nonrank1-autopsy-{component}.jsonl"
    if not path.exists():
        return prior
    with open(path) as f:
        for line in f:
            d = json.loads(line)
            prior[d["locus"]] = d
    return prior


def fold_len(receipt, index):
    return receipt["fold_identities"][index]["length"]


def autopsy_locus(component, rec, qv_rec, prior_rec):
    lg = rec["log_gap_decomposition"]
    both_placed = lg["both_placed_evidence_gap"]
    log_gap = rec["log_gap"]
    residual = log_gap - both_placed
    winner = rec["best_fold_indices"]
    truth = rec["truth_folds"]
    winner_homozygous = winner[0] == winner[1]
    truth_homozygous = truth[0] == truth[1]
    winner_len_sum = fold_len(rec, winner[0]) + fold_len(rec, winner[1])
    truth_len_sum = fold_len(rec, truth[0]) + fold_len(rec, truth[1])
    identity = qv_rec.get("identity") if qv_rec else None
    mismatches = qv_rec.get("mismatches") if qv_rec else None
    gap_columns = qv_rec.get("gap_columns") if qv_rec else None
    twin = identity is not None and identity >= 0.999
    interior_substitution = (
        mismatches is not None and gap_columns is not None and gap_columns > mismatches
    )
    substitution_better = (
        mismatches is not None and gap_columns is not None and mismatches >= gap_columns
    )
    zygosity_dosage = winner_homozygous != truth_homozygous
    unplaced_assist = residual > 0
    if twin:
        mechanism = "twin_indistinguishable"
    elif interior_substitution:
        mechanism = "interior_substitution"
    elif substitution_better:
        mechanism = "substitution_better_match"
    elif zygosity_dosage:
        mechanism = "zygosity_dosage"
    else:
        mechanism = "shared_evidence_only"
    prior_rival = None
    if prior_rec is not None:
        # The prior autopsy's rival-dominant signature: S1 (unplaced
        # advantage) with NOT S2 (shared evidence NOT truth-favored).
        prior_rival = bool(prior_rec.get("s1_unplaced_mass_advantage")) and not bool(
            prior_rec.get("s2_both_placed_truth_favored")
        )
    return {
        "component": component,
        "locus": rec["locus"],
        "partition": rec["partition"],
        "census_class": (qv_rec or {}).get("class") or (prior_rec or {}).get("class"),
        "truth_rank": rec["truth_rank"],
        "log_gap": log_gap,
        "both_placed_evidence_gap": both_placed,
        "residual_unplaced_gap": residual,
        "unplaced_assist": unplaced_assist,
        "winner_folds": winner,
        "winner_homozygous": winner_homozygous,
        "winner_len_sum": winner_len_sum,
        "truth_folds": truth,
        "truth_homozygous": truth_homozygous,
        "truth_len_sum": truth_len_sum,
        "called_class_count": rec["called_class_count"],
        "truth_in_called_set": rec.get("truth_in_called_set"),
        "winner_unplaced_mass": lg["winner_unplaced_mass"],
        "truth_unplaced_mass": lg["truth_unplaced_mass"],
        "winner_unplaced_units": lg["winner_unplaced_units"],
        "truth_unplaced_units": lg["truth_unplaced_units"],
        "qv": (qv_rec or {}).get("qv"),
        "identity": identity,
        "mismatches": mismatches,
        "gap_columns": gap_columns,
        "total_edits": (qv_rec or {}).get("total_edits"),
        "perfect": (qv_rec or {}).get("perfect"),
        "twin_indistinguishable": twin,
        "interior_substitution": interior_substitution,
        "substitution_better_match": substitution_better,
        "zygosity_dosage": zygosity_dosage,
        "mechanism": mechanism,
        "prior_rival_dominant": prior_rival,
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--components", default=None)
    ap.add_argument("--out-suffix", default="")
    args = ap.parse_args()
    components = args.components.split(",") if args.components else ALL_COMPONENTS

    records = []
    per_component = {}
    for component in components:
        receipt_path = D / f"realign-exhaustive-territory-{component}.jsonl"
        if not receipt_path.exists():
            say(f"== {component}: territory receipt MISSING, skipped")
            continue
        qv = load_qv(component)
        prior = load_prior(component)
        kept = 0
        nonrank1 = 0
        out_path = D / f"realign-rival-autopsy-{component}{args.out_suffix}.jsonl"
        with open(receipt_path) as f, open(out_path, "w") as out:
            for line in f:
                rec = json.loads(line)
                if not rec.get("truth_pair_expressible"):
                    continue
                if rec["truth_rank"] == 1:
                    continue
                nonrank1 += 1
                lg = rec["log_gap_decomposition"]
                if lg["both_placed_evidence_gap"] <= 0:
                    continue
                kept += 1
                row = autopsy_locus(component, rec, qv.get(rec["locus"]), prior.get(rec["locus"]))
                records.append(row)
                out.write(json.dumps(row) + "\n")
        per_component[component] = (kept, nonrank1)
        say(f"== {component}: {kept} rival-dominant of {nonrank1} non-rank-1 expressible")

    tables_path = D / f"realign-rival-autopsy-tables{args.out_suffix}.txt"
    with open(tables_path, "w") as g:
        g.write("STAGE 4 — THE RIVAL-DOMINANT AUTOPSY UNDER THE NORMALIZED RULE\n")
        g.write("(cohort: expressible, truth not rank-1, both_placed_evidence_gap > 0)\n\n")
        g.write("PER COMPONENT (cohort / non-rank-1 expressible):\n")
        for component, (kept, nonrank1) in per_component.items():
            g.write(f"  {component}: {kept} / {nonrank1}\n")
        total = len(records)
        prior_rival_count = sum(1 for r in records if r["prior_rival_dominant"])
        g.write(f"\nCOHORT: {total} loci; carrying the pre-normalization rival-dominant "
                f"signature: {prior_rival_count}\n\n")
        g.write("MECHANISM CLASSES (primary, stated precedence: twin-indistinguishable, "
                "interior-substitution, substitution-better-match, zygosity-dosage, "
                "shared-evidence-only):\n")
        by_mech = {}
        for r in records:
            by_mech.setdefault(r["mechanism"], []).append(r)
        for mech in sorted(by_mech, key=lambda m: -len(by_mech[m])):
            rows = by_mech[mech]
            g.write(
                f"  {mech}: {len(rows)} "
                f"(median both_placed gap {median([r['both_placed_evidence_gap'] for r in rows]):.1f}, "
                f"median QV {median([r['qv'] for r in rows if r['qv'] is not None])}, "
                f"median identity {median([r['identity'] for r in rows if r['identity'] is not None])})\n"
            )
            g.write("    by census class: " + json.dumps({
                c: sum(1 for r in rows if r["census_class"] == c)
                for c in sorted({r["census_class"] for r in rows})
            }) + "\n")
        g.write("\nFEATURE CROSS-TAB (not exclusive):\n")
        for feature in ["twin_indistinguishable", "interior_substitution",
                        "substitution_better_match", "zygosity_dosage", "unplaced_assist"]:
            g.write(f"  {feature}: {sum(1 for r in records if r[feature])}\n")
        homo_w = sum(1 for r in records if r["winner_homozygous"])
        g.write(f"  winner homozygous: {homo_w} / {total}\n")
        g.write("\nPER-LOCUS RECORDS (component locus partition census_class mechanism "
                "both_placed residual qv identity winner truth):\n")
        for r in sorted(records, key=lambda r: -r["both_placed_evidence_gap"]):
            g.write(
                f"  {r['component']} L{r['locus']} p{r['partition']} "
                f"{r['census_class']} {r['mechanism']} "
                f"bp={r['both_placed_evidence_gap']:.1f} res={r['residual_unplaced_gap']:.1f} "
                f"qv={r['qv']} id={r['identity']} "
                f"winner={r['winner_folds']} ({'hom' if r['winner_homozygous'] else 'het'}, "
                f"{r['winner_len_sum']}bp) truth={r['truth_folds']} "
                f"({'hom' if r['truth_homozygous'] else 'het'}, {r['truth_len_sum']}bp)\n"
            )
    say(f"TOTAL COHORT {total}; tables at {tables_path}")


if __name__ == "__main__":
    main()
