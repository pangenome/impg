#!/usr/bin/env python3
"""Independent receipt checker for the read-matched observed-convention stage
(owner-approved option b, the likegt fix).

Validates, from the receipts on disk and nothing else:
  phase 1 — the run markers (the gate-off identity run; the read-matched
            chrMT/chrI runs: exit 0, walls, RSS, stage heartbeats);
  phase 2 — GATE-OFF IDENTITY: the identity receipt (new binary, gate unset)
            reproduces the committed remedy receipt on every semantic field
            (only the documented HashMap-order observed-norm ULPs differ);
            sidecars and the admission diagnostic byte-identical;
  phase 3 — THE SMEAR MEASUREMENT FIRST: per locus, observed node+edge mass
            before (remedy) vs after (read-matched) — asserted ZERO change at
            every locus; the emission census from the run logs (gapped
            occurrences 0, span-bp == matched-bp: no record's chained span
            contains unmatched material); per-key equality via the
            byte-identical ingredients/records sidecars; and the SEATS: the
            off-truth-route winner-extra mass per locus, re-attributed as
            GENUINELY MATCHED sample evidence (not span-containment credit);
  phase 4 — the independent likelihood re-derivation from the readmatched
            ingredients (the closed form already validated bit-close against
            the committed receipts), reproducing winner and truth rank;
  phase 5 — the before/after truth-rank/log-gap/QUAL table per locus and the
            NO-REGRESSION GUARD: chrMT L3/L5 and chrI L4/5/6/L12 rank 1, the
            chrI L14 bit-tie with the truth class;
  phase 6 — the residual classification per still-open locus with the
            corrected attribution (matched off-truth-route evidence vs
            cross-window frame mass vs door policy);
  phase 7 — walls and RSS peaks under the 64GiB guard.

Assessment-side only. No product file is read or written.
"""

import json
import math
import os
import re
import sys

DATA = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
FAILURES = []
GUARD_KIB = 64 * 1024 * 1024
NO_REGRESSION_RANK1 = [("chrMT", 3), ("chrMT", 5), ("chrI", 4), ("chrI", 5), ("chrI", 6), ("chrI", 12)]
TIE_LOCUS = ("chrI", 14)


def fail(message):
    FAILURES.append(message)
    print(f"FAIL: {message}")


def check(condition, message):
    if condition:
        print(f"ok: {message}")
    else:
        fail(message)


def load_census(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            out[row["locus"]] = row
    return out


def load_diagnostic(path):
    header = None
    loci = {}
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            if row.get("kind") == "admission-diagnostic-header":
                header = row
            else:
                loci[row["locus"]] = row
    return header, loci


def load_ingredients(path):
    out = {}
    with open(path) as handle:
        for line in handle:
            row = json.loads(line)
            out[row["locus"]] = row
    return out


def rank_locus(ingredients):
    """The full per-locus ranking recomputed from the ingredients with the
    exact closed form (validated bit-close against the committed receipts).
    Returns (ranking[(row_i,row_j)] -> logL sorted desc, mass, universe)."""
    rows = ingredients["rows"]
    un = set()
    ue = set()
    for row in rows:
        for entry in row["nodes"]:
            un.add(entry[0])
        for entry in row["edges"]:
            ue.add(entry[0])
    counts = {}
    for record in ingredients["records"]:
        share = record["share"]
        for key in record["nodes"]:
            if key in un:
                counts[key] = counts.get(key, 0.0) + share
        for key in record["edges"]:
            if key in ue:
                counts[key] = counts.get(key, 0.0) + share
    mass = sum(counts.values())
    universe = len(un) + len(ue)
    if mass <= 0:
        return None, 0.0, universe
    gamma = sum(math.lgamma(1.0 + count) for count in counts.values())
    beta = mass / universe
    ranking = []
    for i in range(len(rows)):
        for j in range(i, len(rows)):
            merged = {}
            for index in (i, j):
                for kind in ("nodes", "edges"):
                    for key, _mult, exposure in rows[index][kind]:
                        merged.setdefault(key, [0.0])[0] += exposure
            explainable = True
            a_term = 0.0
            covered = 0.0
            keys = 0
            for key, (exposure,) in merged.items():
                keys += 1
                count = counts.get(key, 0.0)
                if count > 0:
                    if exposure > 0:
                        a_term += count * math.log(exposure)
                        covered += count
                    else:
                        explainable = False
                        break
            if not explainable:
                continue
            total_exposure = sum(value[0] for value in merged.values())
            if covered > 0 and total_exposure <= 0:
                continue
            explained = mass if total_exposure > 0 else 0.0
            constant = explained + gamma + universe * beta - mass * math.log(beta)
            mass_term = covered * math.log(total_exposure * beta / mass) if covered > 0 else 0.0
            nll = constant + mass_term - a_term - keys * beta
            ranking.append(((i, j), -nll))
    ranking.sort(key=lambda entry: -entry[1])
    return ranking, mass, universe


def parse_census(err_path):
    with open(err_path) as handle:
        text = handle.read()
    match = re.search(
        r"\[read-matched\] records (\d+) occurrences (\d+) "
        r"gapped-occurrences (\d+) span-bp (\d+) matched-bp (\d+)",
        text,
    )
    if not match:
        return None
    values = tuple(int(value) for value in match.groups())
    return {
        "records": values[0],
        "occurrences": values[1],
        "gapped": values[2],
        "span_bp": values[3],
        "matched_bp": values[4],
    }


def semantic_diffs(old, new):
    """Field-level comparison ignoring the documented HashMap-order
    observed-norm ULPs; returns (other-diffs, norm-ulp-diffs)."""
    other = []
    ulp = []
    for locus, row in new.items():
        assert locus in old
        for key, value in row.items():
            if key in ("observed_node_norm", "observed_edge_norm"):
                if old[locus].get(key) != value:
                    ulp.append((locus, key))
            elif old[locus].get(key) != value:
                other.append((locus, key))
    return other, ulp


def window_of(axis_windows, bp):
    for index, (start, end) in enumerate(axis_windows):
        if start <= bp < end:
            return index
    return None


def main():
    components = ("chrMT", "chrI")
    before = {c: load_census(f"{DATA}/cosine-graph-likelihood-remedy-{c}.jsonl") for c in components}
    after = {
        c: load_census(f"{DATA}/cosine-graph-likelihood-readmatched-{c}.jsonl") for c in components
    }
    identity = load_census(f"{DATA}/cosine-graph-likelihood-identity-chrMT.jsonl")
    diag = {
        c: load_diagnostic(f"{DATA}/cosine-admission-diagnostic-readmatched-{c}.jsonl") for c in components
    }
    ingredients = {
        c: load_ingredients(f"{DATA}/cosine-graph-likelihood-readmatched-{c}.jsonl.ingredients.jsonl")
        for c in components
    }
    ingredients_before = {
        c: load_ingredients(f"{DATA}/cosine-graph-likelihood-remedy-{c}.jsonl.ingredients.jsonl")
        for c in components
    }

    # ---------------- phase 1: run markers.
    print("== phase 1: run markers")
    runs = [("run-cosine-identity-pilot", ("chrMT",)), ("run-cosine-readmatched-pilot", components)]
    for prefix, comps in runs:
        for component in comps:
            exit_path = f"{DATA}/{prefix}-{component}.exit"
            check(os.path.exists(exit_path), f"{prefix}-{component}.exit exists")
            with open(exit_path) as handle:
                check(handle.read().strip() == "0", f"{prefix}-{component} exit 0")
            for sidecar in ("done", "wall", "rss", "stages"):
                check(
                    os.path.exists(f"{DATA}/{prefix}-{component}.{sidecar}"),
                    f"{prefix}-{component}.{sidecar}",
                )

    # ---------------- phase 2: gate-off identity.
    print("\n== phase 2: gate-off identity (new binary, gate unset)")
    other, ulp = semantic_diffs(before["chrMT"], identity)
    check(not other, f"identity receipt matches the committed remedy receipt semantically ({other})")
    check(ulp, "the only identity differences are the documented observed-norm ULPs")
    for sidecar in ("ties", "competitors", "records", "ingredients"):
        same = os.path.exists(f"{DATA}/cosine-graph-likelihood-identity-chrMT.jsonl.{sidecar}.jsonl") and os.path.exists(
            f"{DATA}/cosine-graph-likelihood-remedy-chrMT.jsonl.{sidecar}.jsonl"
        )
        if same:
            with open(f"{DATA}/cosine-graph-likelihood-identity-chrMT.jsonl.{sidecar}.jsonl", "rb") as a:
                with open(f"{DATA}/cosine-graph-likelihood-remedy-chrMT.jsonl.{sidecar}.jsonl", "rb") as b:
                    same = a.read() == b.read()
        check(same, f"identity {sidecar} sidecar byte-identical with the remedy receipt")
    with open(f"{DATA}/cosine-admission-diagnostic-identity-chrMT.jsonl", "rb") as a:
        with open(f"{DATA}/cosine-admission-diagnostic-remedy-chrMT.jsonl", "rb") as b:
            check(a.read() == b.read(), "identity admission diagnostic byte-identical")

    # ---------------- phase 3: THE SMEAR MEASUREMENT FIRST.
    print("\n== phase 3: the smear measurement (observed mass, before -> after)")
    for component in components:
        total_before = total_after = 0.0
        for locus in sorted(after[component]):
            old_row, new_row = before[component][locus], after[component][locus]
            old_mass = old_row["observed_node_mass"] + old_row["observed_edge_mass"]
            new_mass = new_row["observed_node_mass"] + new_row["observed_edge_mass"]
            total_before += old_mass
            total_after += new_mass
            check(
                old_mass == new_mass,
                f"{component} L{locus} observed mass unchanged under the convention "
                f"({old_mass:.1f} -> {new_mass:.1f})",
            )
        print(
            f"   {component}: total observed mass {total_before:.1f} -> {total_after:.1f} "
            f"(removed {100.0 * (total_before - total_after) / total_before:.1f}%)"
        )
        census = parse_census(f"{DATA}/run-cosine-readmatched-pilot-{component}.err")
        check(census is not None, f"{component} run log carries the read-matched emission census")
        if census:
            check(
                census["gapped"] == 0,
                f"{component} census: gapped occurrences 0 (of {census['occurrences']})",
            )
            check(
                census["span_bp"] == census["matched_bp"],
                f"{component} census: span-bp == matched-bp ({census['span_bp']} == {census['matched_bp']}) "
                "- no record's chained span contains unmatched material",
            )
            print(
                f"   {component} census: records {census['records']}, occurrences {census['occurrences']}, "
                f"gapped {census['gapped']}, span-bp {census['span_bp']} == matched-bp {census['matched_bp']}"
            )
        for sidecar in ("ties", "competitors", "records", "ingredients"):
            with open(f"{DATA}/cosine-graph-likelihood-readmatched-{component}.jsonl.{sidecar}.jsonl", "rb") as a:
                with open(f"{DATA}/cosine-graph-likelihood-remedy-{component}.jsonl.{sidecar}.jsonl", "rb") as b:
                    check(
                        a.read() == b.read(),
                        f"{component} {sidecar} sidecar byte-identical (per-record covered key sets unchanged)",
                    )
        with open(f"{DATA}/cosine-admission-diagnostic-readmatched-{component}.jsonl", "rb") as a:
            with open(f"{DATA}/cosine-admission-diagnostic-remedy-{component}.jsonl", "rb") as b:
                check(a.read() == b.read(), f"{component} read-matched admission diagnostic byte-identical")

    # The seats: per still-open locus, the winner's extra-key mass that sits
    # OFF the truth's routes (the remedy stage's "smear" seats) - reported with
    # the corrected attribution: GENUINELY MATCHED sample evidence.
    print("\n== phase 3c: the smear seats, re-attributed (winner-extra mass off the truth's routes)")
    seat_totals = {}
    for component in components:
        header, loci = diag[component]
        axis_windows = [(entry["start"], entry["end"]) for entry in header["axis"]]
        for locus in sorted(after[component]):
            census = after[component][locus]
            if not census.get("truth_pair_expressible") or census["truth_rank"] == 1:
                seat_totals[(component, locus)] = None
                continue
            ingredient = ingredients[component][locus]
            truth_rows = set(census["truth_graph_rows"])
            winner_rows = census["qual_called_classes"][0]["row_indices"]

            def usage(index):
                row = ingredient["rows"][index]
                return (
                    set(entry[0] for entry in row["nodes"]),
                    set(entry[0] for entry in row["edges"]),
                )

            truth_nodes, truth_edges = set(), set()
            for index in truth_rows:
                nodes, edges = usage(index)
                truth_nodes |= nodes
                truth_edges |= edges
            winner_nodes, winner_edges = set(), set()
            for index in winner_rows:
                nodes, edges = usage(index)
                winner_nodes |= nodes
                winner_edges |= edges
            node_mass = {entry[0]: entry[1] for entry in loci[locus]["universe_nodes"]}
            edge_mass = {entry[0]: entry[1] for entry in loci[locus]["universe_edges"]}
            node_pos = {entry[0]: entry[2] for entry in loci[locus]["universe_nodes"]}
            edge_pos = {entry[0]: entry[2] for entry in loci[locus]["universe_edges"]}
            off_truth = 0.0
            cross_window = {}
            local = 0.0
            for key in list(winner_nodes - truth_nodes) + list(winner_edges - truth_edges):
                mass = node_mass.get(key, 0.0) + edge_mass.get(key, 0.0)
                positions = node_pos.get(key, []) or edge_pos.get(key, []) or []
                if not positions:
                    off_truth += mass
                    continue
                windows = set(window_of(axis_windows, position[1]) for position in positions)
                if windows == {locus}:
                    local += mass
                else:
                    for window in windows:
                        cross_window[window] = cross_window.get(window, 0.0) + mass / len(windows)
            locus_mass = census["observed_node_mass"] + census["observed_edge_mass"]
            seat_totals[(component, locus)] = (off_truth, cross_window, local, locus_mass)
            print(
                f"   {component} L{locus}: winner-extra off-truth-route mass {off_truth:.1f} "
                f"({100.0 * off_truth / locus_mass:.1f}% of locus mass) - MATCHED sample evidence "
                "(the convention removes none of it; no unmatched span credit exists)"
            )

    # ---------------- phase 4: independent re-derivation.
    print("\n== phase 4: independent likelihood re-derivation (read-matched ingredients)")
    for component in components:
        for locus, census in sorted(after[component].items()):
            ranking, _mass, _universe = rank_locus(ingredients[component][locus])
            best = ranking[0][1] if ranking else None
            if census.get("best_log_likelihood") is not None and best is not None:
                scale = abs(census["best_log_likelihood"])
                delta = abs(census["best_log_likelihood"] - best)
                check(delta <= 1e-6 * scale, f"{component} L{locus} winner logL re-derived")
            if census.get("truth_pair_expressible"):
                truth_rows = census["truth_graph_rows"]
                truth_logl = None
                for pair, logl in ranking:
                    if set(pair) == set(truth_rows):
                        truth_logl = logl
                        break
                check(truth_logl is not None, f"{component} L{locus} truth class re-derived")
                if truth_logl is not None:
                    scale = abs(truth_logl)
                    strictly_above = sum(1 for _pair, logl in ranking if logl > truth_logl + 1e-9 * scale)
                    check(
                        strictly_above + 1 == census["truth_rank"],
                        f"{component} L{locus} truth rank {census['truth_rank']} re-derived",
                    )

    # ---------------- phase 5: before/after table + no-regression guard.
    print("\n== phase 5: before/after table (truth rank / log-gap / QUAL; identical by measurement)")
    for component in components:
        print(f"   {component}:")
        print("   locus  rank before->after   gap before->after    QUAL before->after")
        for locus in sorted(after[component]):
            old_row, new_row = before[component][locus], after[component][locus]

            def summary(row):
                if not row.get("truth_pair_expressible"):
                    return "non-expr", "non-expr", str(row.get("qual"))
                gap = row["best_log_likelihood"] - row["truth_log_likelihood"]
                return row["truth_rank"], f"{gap:.1f}", str(row.get("qual"))

            print(f"   L{locus:>2}  {summary(old_row)[0]} -> {summary(new_row)[0]}   "
                  f"{summary(old_row)[1]} -> {summary(new_row)[1]}   "
                  f"{summary(old_row)[2]} -> {summary(new_row)[2]}")
    print("\n== no-regression guard (the current wins)")
    for component, locus in NO_REGRESSION_RANK1:
        old_row, new_row = before[component][locus], after[component][locus]
        stays = (
            new_row.get("truth_pair_expressible")
            and new_row["truth_rank"] == 1
            and old_row.get("truth_rank") == 1
        )
        check(stays, f"{component} L{locus} rank-1 win preserved under the read-matched convention")
    component, locus = TIE_LOCUS
    tie_row = after[component][locus]
    check(
        tie_row.get("ties_with_truth_streamed", 0) >= 1,
        f"{component} L{locus} bit-tie with the truth class preserved",
    )

    # ---------------- phase 6: residual classification, corrected attribution.
    print("\n== phase 6: residual classification (still-open loci, corrected attribution)")
    for (component, locus), seats in seat_totals.items():
        if seats is None:
            continue
        off_truth, cross_window, local, locus_mass = seats
        census = after[component][locus]
        gap = census["best_log_likelihood"] - census["truth_log_likelihood"]
        parts = [
            f"MATCHED off-truth-route evidence {off_truth:.1f} "
            f"({100.0 * off_truth / locus_mass:.1f}% of locus mass)"
        ]
        if local > 0:
            parts.append(f"door-policy(local truth mass un-spelled) {local:.1f}")
        for window, mass_value in sorted(cross_window.items(), key=lambda kv: -kv[1]):
            parts.append(f"frame(window {window}) {mass_value:.1f}")
        print(f"   {component} L{locus}: rank {census['truth_rank']}, gap {gap:.1f} — " + "; ".join(parts))

    # ---------------- phase 7: walls and RSS.
    print("\n== phase 7: walls and RSS (the 64GiB guard)")
    for prefix, comps in runs:
        for component in comps:
            with open(f"{DATA}/{prefix}-{component}.wall") as handle:
                wall = int(handle.read().strip())
            peak = 0
            with open(f"{DATA}/{prefix}-{component}.rss") as handle:
                for line in handle:
                    fields = line.split()
                    if len(fields) >= 3:
                        peak = max(peak, int(fields[2]))
            check(peak < GUARD_KIB, f"{prefix}-{component}: peak RSS {peak} kB under the 64GiB guard")
            print(f"   {prefix}-{component}: wall {wall}s, external RSS peak {peak} kB")

    print()
    if FAILURES:
        print(f"CHECK FAILURES: {len(FAILURES)}")
        return 1
    print("ALL PHASES PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
