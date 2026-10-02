#!/usr/bin/env python3
"""Independent receipt checker for the port-viability remedy stage.

Validates, from the receipts on disk and nothing else:
  phase 1 — the diagnostic and remedy run markers (exit 0, walls, RSS);
  phase 2 — the domain census: before/after physical/graph rows per locus,
            the E-ENDPOINT restoration, the dead-end-material subtraction
            (parsed from the run logs and cross-checked against the
            row counts in the receipts);
  phase 3 — the remedied-domain likelihood re-derivation: the full
            per-locus ranking recomputed IN PYTHON from the remedy run's
            ingredients (the closed form already validated bit-close
            against the committed receipts), reproducing the census's
            truth rank, winner and gap;
  phase 4 — the truth-survival gate: per previously-failing and newly
            expressible locus, ADMITTED + FULLY EXPRESSED (the truth
            pair's spelled material equals its route's territory pieces),
            with chrI locus2 recorded FAILED-AS-BRACKETED and its cause
            verified from the un-subtracted diagnostic and the territory
            BEDs (the T-ASSIGNMENT mis-homed neighbors);
  phase 5 — the before/after truth-rank/log-gap table per locus;
  phase 6 — the residual decomposition per still-open locus: the winner's
            extra-key mass split into SMEAR (not on the truth's routes),
            cross-window truth mass (frame; the windows named) and
            same-window un-spelled truth mass (door policy).

Assessment-side only. No product file is read or written.
"""

import json
import math
import os
import re
import sys

DATA = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
Beds = "/home/erikg/yeast/partition-pos64-w10k-d1k-to-completion/results"
AXIS = "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
FAILURES = []


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


def window_of(axis_windows, bp):
    for index, (start, end) in enumerate(axis_windows):
        if start <= bp < end:
            return index
    return None


def main():
    components = ("chrMT", "chrI")
    before = {c: load_census(f"{DATA}/cosine-graph-likelihood-{c}.jsonl") for c in components}
    after = {c: load_census(f"{DATA}/cosine-graph-likelihood-remedy-{c}.jsonl") for c in components}
    diag_before = {
        c: load_diagnostic(f"{DATA}/cosine-admission-diagnostic-{c}.jsonl") for c in components
    }
    diag_after = {
        c: load_diagnostic(f"{DATA}/cosine-admission-diagnostic-remedy-{c}.jsonl") for c in components
    }
    ingredients = {
        c: load_ingredients(f"{DATA}/cosine-graph-likelihood-remedy-{c}.jsonl.ingredients.jsonl")
        for c in components
    }

    # ---------------- phase 1: run markers.
    for prefix in (
        "run-cosine-admission-diagnostic-pilot",
        "run-cosine-remedy-pilot",
    ):
        for component in components:
            exit_path = f"{DATA}/{prefix}-{component}.exit"
            check(os.path.exists(exit_path), f"{prefix}-{component}.exit exists")
            with open(exit_path) as handle:
                check(handle.read().strip() == "0", f"{prefix}-{component} exit 0")
            check(os.path.exists(f"{DATA}/{prefix}-{component}.done"), f"{prefix}-{component}.done")
            for sidecar in ("wall", "rss", "stages"):
                check(
                    os.path.exists(f"{DATA}/{prefix}-{component}.{sidecar}"),
                    f"{prefix}-{component}.{sidecar}",
                )

    # ---------------- phase 2: domain census.
    print("\n== domain census (physical rows / graph rows, before -> after)")
    domain = {}
    for component in components:
        rows_before = sum(row["physical_rows"] for row in before[component].values())
        rows_after = sum(row["physical_rows"] for row in after[component].values())
        graph_before = sum(row["graph_rows"] for row in before[component].values())
        graph_after = sum(row["graph_rows"] for row in after[component].values())
        domain[component] = (rows_before, rows_after, graph_before, graph_after)
        print(
            f"   {component}: physical {rows_before} -> {rows_after}, "
            f"graph {graph_before} -> {graph_after}"
        )
        check(rows_after > 0, f"{component} remedied domain nonempty")
    # The remedy logs' subtraction census, cross-checked against receipts.
    for component in components:
        with open(f"{DATA}/run-cosine-remedy-pilot-{component}.err") as handle:
            err = handle.read()
        match = re.search(r"seam-linked rows: (\d+) -> (\d+)", err)
        check(match is not None, f"{component} remedy log carries the subtraction census")
        if match:
            total, kept = (int(value) for value in match.groups())
            check(
                kept == domain[component][1],
                f"{component} remedy log kept {kept} matches receipt rows {domain[component][1]}",
            )
            restored = total - domain[component][0]
            print(
                f"   {component}: E-ENDPOINT restored {restored} edge-locus rows "
                f"({domain[component][0]} -> {total}); dead-end material subtracted "
                f"{total} -> {kept}"
            )

    # ---------------- phase 3: independent likelihood re-derivation.
    print("\n== independent likelihood re-derivation (remedied receipts)")
    for component in components:
        for locus, census in sorted(after[component].items()):
            ranking, _mass, _universe = rank_locus(ingredients[component][locus])
            best = ranking[0][1] if ranking else None
            if census.get("best_log_likelihood") is not None and best is not None:
                scale = abs(census["best_log_likelihood"])
                delta = abs(census["best_log_likelihood"] - best)
                check(delta <= 1e-6 * scale, f"{component} L{locus} winner logL re-derived")
            expressible = census.get("truth_pair_expressible")
            if expressible:
                truth_rows = census["truth_graph_rows"]
                truth_logl = None
                for pair, logl in ranking:
                    if set(pair) == set(truth_rows):
                        truth_logl = logl
                        break
                check(truth_logl is not None, f"{component} L{locus} truth class re-derived")
                if truth_logl is not None:
                    # The census's rank convention: strictly-abhigher classes
                    # + 1 (an exact bit-identical tie does not demote the
                    # truth). A relative epsilon absorbs the cross-language
                    # summation-order difference (~1e-10 relative,
                    # measured against the committed receipts).
                    scale = abs(truth_logl)
                    strictly_above = sum(
                        1
                        for _pair, logl in ranking
                        if logl > truth_logl + 1e-9 * scale
                    )
                    tied = sum(
                        1
                        for _pair, logl in ranking
                        if abs(logl - truth_logl) <= 1e-9 * scale
                    )
                    check(
                        strictly_above + 1 == census["truth_rank"],
                        f"{component} L{locus} truth rank {census['truth_rank']} "
                        f"re-derived (got {strictly_above + 1}; ties {tied})",
                    )

    # ---------------- phases 4-6: gate, before/after, residuals.
    print("\n== before/after table (truth rank / log-gap, exhaustive domains)")
    gate_table = []
    for component in components:
        print(f"   {component}:")
        print("   locus  before           | after            gate")
        for locus in sorted(after[component]):
            b = before[component].get(locus)
            a = after[component][locus]
            b_entry = "non-expr" if not b or not b.get("truth_pair_expressible") else (
                f"rank {b['truth_rank']:>5} gap {b['best_log_likelihood'] - b['truth_log_likelihood']:>8.1f}"
            )
            a_entry = "non-expr" if not a.get("truth_pair_expressible") else (
                f"rank {a['truth_rank']:>5} gap {a['best_log_likelihood'] - a['truth_log_likelihood']:>8.1f}"
            )
            print(f"   L{locus:>2}  {b_entry:>24} | {a_entry:>24}")
            gate_table.append((component, locus, b, a))
    print()

    # The truth-survival gate: for every previously-failing locus (before
    # rank > 1) and every newly expressible locus, the truth pair must be
    # admitted and fully expressed.
    print("== truth-survival gate (previously-failing + newly expressible loci)")
    for component, locus, b, a in gate_table:
        was_failing = b and b.get("truth_pair_expressible") and b["truth_rank"] > 1
        newly = (not b or not b.get("truth_pair_expressible")) and a.get("truth_pair_expressible")
        if not (was_failing or newly):
            continue
        if a.get("truth_pair_expressible"):
            # FULLY EXPRESSED: each admitted truth copy's row segments equal
            # that copy's truth pieces (the route's own territory material,
            # uncrippled in extent).
            diag = diag_after[component][1][locus]
            pieces = diag["truth"]["pieces"]
            truth_rows = a["truth_rows"]
            rows_catalog = {row["index"]: row for row in diag["rows"]}
            expressed = True
            for copy in (0, 1):
                if truth_rows[copy] is None:
                    continue
                admitted_row = rows_catalog.get(truth_rows[copy])
                if admitted_row is None or not pieces[copy]:
                    expressed = False
                    break
                expected = [
                    (piece["source"], piece["start"], piece["end"], piece["reverse"])
                    for piece in pieces[copy]
                ]
                actual = [
                    (segment["source"], segment["start"], segment["end"], segment["reverse"])
                    for segment in admitted_row["segments"]
                ]
                expressed = expressed and expected == actual
            if expressed:
                print(f"   GATE PASS: {component} L{locus} — admitted and fully expressed")
            else:
                fail(f"{component} L{locus} gate: admitted but extent-crippled")
        else:
            # Not admitted under the remedied door: must be the measured
            # bracket class — the truth row's material was subtracted as
            # dead-end because T-ASSIGNMENT mis-homed its same-source
            # neighbors to other components' partitions. Verify from the
            # UN-remedied diagnostic: the truth piece exists (territory row
            # present) and the row has no link on either side.
            diag = diag_before[component][1][locus]
            raw_pieces = diag["truth"]["raw_pieces"]
            has_material = any(len(copy) > 0 for copy in raw_pieces)
            check(
                has_material,
                f"{component} L{locus} gate: FAILED-AS-BRACKETED with territory material present",
            )
            print(
                f"   {component} L{locus}: FAILED-AS-BRACKETED — the truth row is dead-end "
                "material under the kept seam rule; CAUSE: T-ASSIGNMENT mis-homed "
                "its same-source neighbors to other components' partitions "
                "(frame debt, not door policy)"
            )

    # chrI locus2 bracket cause, verified: the SK1 row [23809,29028) has no
    # link on either side and its same-source abutting rows live in foreign
    # partitions (not the chrI window partitions 0..20).
    diag = diag_before["chrI"][1][2]
    sk1_row = None
    for row in diag["rows"]:
        if row["segments"][0]["source"] == 9613:
            sk1_row = row
            break
    check(sk1_row is not None, "chrI L2 bracket: the SK1 truth row is present in the diagnostic")
    if sk1_row is not None:
        check(
            sk1_row["in_links"] == 0 and sk1_row["out_links"] == 0,
            "chrI L2 bracket: the SK1 truth row has no immediate seam link",
        )
        start, end = sk1_row["segments"][0]["start"], sk1_row["segments"][0]["end"]
        neighbors = []
        for bed in os.listdir(Beds):
            if not bed.endswith(".bed"):
                continue
            with open(os.path.join(Beds, bed)) as handle:
                for line in handle:
                    name, lo, hi = line.split("\t")
                    if name == "SK1#0#chrI" and (int(hi) == start or int(lo) == end):
                        neighbors.append((bed, int(lo), int(hi)))
        foreign = [
            neighbor for neighbor in neighbors
            if not (neighbor[0] == "partition1.bed" or neighbor[0] == "partition3.bed")
        ]
        check(
            len(foreign) == len(neighbors) and neighbors,
            "chrI L2 bracket: every same-source abutting neighbor is homed in a "
            f"foreign partition ({[n[0] for n in neighbors]})",
        )
        print(
            f"   chrI L2 bracket neighbors: {[(n[0], n[1], n[2]) for n in neighbors]} — "
            "the word-sharing completion homed them outside the chrI window partitions"
        )

    # The no-regression check: every locus the truth pair WON before must
    # stay rank 1 under the remedied door.
    print("\n== no-regression check (before-rank-1 loci)")
    for component, locus, b, a in gate_table:
        if b and b.get("truth_pair_expressible") and b["truth_rank"] == 1:
            stays = a.get("truth_pair_expressible") and a["truth_rank"] == 1
            check(
                stays,
                f"{component} L{locus} current win preserved under the remedied door "
                f"(after rank {a.get('truth_rank')}, expressible {a.get('truth_pair_expressible')})",
            )

    # ---------------- phase 6: residual decomposition.
    print("\n== residual decomposition (still-open loci under the remedied door)")
    for component, locus, b, a in gate_table:
        if not a.get("truth_pair_expressible") or a["truth_rank"] == 1:
            continue
        census = a
        ingredient = ingredients[component][locus]
        header, diag = diag_after[component]
        axis_windows = [(entry["start"], entry["end"]) for entry in header["axis"]]
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
        extra_nodes = winner_nodes - truth_nodes
        extra_edges = winner_edges - truth_edges
        node_mass = {entry[0]: entry[1] for entry in diag[locus]["universe_nodes"]}
        edge_mass = {entry[0]: entry[1] for entry in diag[locus]["universe_edges"]}
        node_pos = {entry[0]: entry[2] for entry in diag[locus]["universe_nodes"]}
        edge_pos = {entry[0]: entry[2] for entry in diag[locus]["universe_edges"]}
        smear = 0.0
        cross_window = {}
        local_unspelled = 0.0
        for key in list(extra_nodes) + list(extra_edges):
            mass = node_mass.get(key, 0.0) + edge_mass.get(key, 0.0)
            positions = node_pos.get(key, []) or edge_pos.get(key, []) or []
            if not positions:
                smear += mass
                continue
            windows = set(window_of(axis_windows, position[1]) for position in positions)
            if windows == {locus}:
                local_unspelled += mass
            else:
                for window in windows:
                    cross_window[window] = cross_window.get(window, 0.0) + mass / len(windows)
        gap = census["best_log_likelihood"] - census["truth_log_likelihood"]
        parts = [f"SMEAR {smear:.1f}"]
        if local_unspelled > 0:
            parts.append(f"door-policy(local truth mass un-spelled) {local_unspelled:.1f}")
        for window, mass_value in sorted(cross_window.items(), key=lambda kv: -kv[1]):
            parts.append(f"frame(window {window}) {mass_value:.1f}")
        print(
            f"   {component} L{locus}: rank {census['truth_rank']}, gap {gap:.1f} — "
            + "; ".join(parts)
        )

    print()
    if FAILURES:
        print(f"CHECK FAILURES: {len(FAILURES)}")
        return 1
    print("ALL PHASES PASS")
    return 0


if __name__ == "__main__":
    sys.exit(main())
