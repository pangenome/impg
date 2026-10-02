#!/usr/bin/env python3
"""Independently reconcile the graph-space condensation pilot receipts.

Checks, per component (chrMT, chrI), fail-closed and without recomputing
scores: receipt markers; locus alignment with the balanced Stage-1 census and
the row-space pilot; pure-route physical domain identity; truth survival masks;
exhaustive pair-count bounds for BOTH arms (nodes-only and nodes+edges);
every streamed competitor's arm flags, scores and truth references; the
path-continuity audit (the chain-junction channel's measured emptiness); and
the arm-delta summary (the edge layer's contribution as separate numbers).

QUAL phase (material-class product semantics, owner ruling 2026-10-01): the
new receipts must reproduce every stage-1 measurement bit-identically, and
per locus the checker re-derives the class share (called class value over the
total similarity of every distinct material class), the Phred QUAL
-10*log10(1 - share), the called-set size/signatures/physical-pair
multiplicity, truth-in-called-set, the k-way-tie derived bound
Q <= -10*log10(1 - 1/k) (equality iff the tie holds all the domain's mass),
and reconciles the tie-evidence stream against the stage-1 tied-pairs
counter, including bit-exactness and node/edge usage distances versus the
truth class.
"""
import json
import math
from pathlib import Path

D = Path('/home/erikg/yeast/genome-balanced-diploid-validation-20260930')
# Stage-1 census fields the re-run must reproduce exactly: integers,
# booleans, masks, row-index vectors. (identity_kinds compared as sorted
# multisets below: it iterates a Rust HashMap whose order is not stable.)
STAGE1_EXACT = [
    'physical_rows', 'graph_rows', 'multi_segment_rows', 'disjoint_material_rows',
    'record_count', 'node_universe_count', 'edge_universe_count',
    'observed_node_mass', 'observed_edge_mass', 'eligible_pairs_nodes',
    'eligible_pairs_combined', 'competitor_entries', 'truth_piece_presence',
    'truth_pair_expressible', 'truth_rows', 'truth_graph_rows',
    'nodes_truth_rank', 'nodes_truth_tied_pairs', 'nodes_higher_competitors',
    'nodes_best_row_indices', 'combined_truth_rank',
    'combined_truth_tied_pairs', 'combined_higher_competitors',
    'combined_best_row_indices', 'rank_delta_combined_vs_nodes',
]
# Floating-point fields accumulated over HashMap iteration (the observed
# norm-of-mass sums) or derived from them (every cosine): the summation
# ORDER varies across processes, so cross-run values can differ by a few
# f64 ULPs (measured max 4e-16 relative) without any rank, tie or count
# changing. Reconciled, not equated: this is a receipt tolerance, not a
# product constant.
STAGE1_FLOAT_ULP = [
    'observed_node_norm', 'observed_edge_norm', 'nodes_truth_cosine',
    'combined_truth_cosine', 'nodes_best_cosine', 'combined_best_cosine',
]

for component, expected in [('chrMT', 8), ('chrI', 17)]:
    prefix = D / f'run-cosine-graph-pilot-{component}'
    assert Path(f'{prefix}.done').exists()
    assert Path(f'{prefix}.exit').read_text().strip() == '0'
    old = json.loads((D / f'run-balanced-{component}.log').read_text())['spine']['stage1_sweep']['loci']
    row_space = [json.loads(line) for line in (D / f'cosine-exhaustive-{component}.jsonl').open()]
    rows = [json.loads(line) for line in (D / f'cosine-graph-exhaustive-{component}.jsonl').open()]
    assert len(rows) == len(row_space) == len(old)
    assert all(row['locus'] == prior['locus'] == sweep['locus']
               for row, prior, sweep in zip(rows, row_space, old))
    # The SAME pure-route domain: identical physical row counts and truth masks.
    assert all(row['physical_rows'] == prior['physical_rows'] for row, prior in zip(rows, row_space))
    assert all(row['truth_pair_expressible'] == prior['truth_pair_expressible']
               for row, prior in zip(rows, row_space))
    assert all(row['truth_piece_presence'] == prior['truth_piece_presence']
               for row, prior in zip(rows, row_space))
    assert sum(row['truth_pair_expressible'] for row in rows) == expected
    # Fail-closed exhaustive enumeration bounds, both arms; graph coalescing
    # never enlarges the domain.
    for row in rows:
        pairs = row['graph_rows'] * (row['graph_rows'] + 1) // 2
        assert row['graph_rows'] <= row['physical_rows']
        assert row['eligible_pairs_nodes'] <= pairs
        assert row['eligible_pairs_combined'] <= pairs
        # The pure-route domain's path-continuity audit: no candidate row can
        # spell a junction adjacency that no single record spans.
        assert row['disjoint_material_rows'] == 0, row['locus']
        assert all(count > 0 for _, count in row['identity_kinds'])
    by_locus = {row['locus']: row for row in rows}
    entries = 0
    per_arm = {'nodes': 0, 'combined': 0}
    empty_usage_competitors = 0
    with (D / f'cosine-graph-exhaustive-{component}.jsonl.competitors.jsonl').open() as stream:
        for line in stream:
            competitor = json.loads(line)
            truth = by_locus[competitor['locus']]
            assert truth['truth_pair_expressible']
            for arm in ('nodes', 'combined'):
                score = competitor[f'{arm}_cosine']
                truth_score = competitor[f'truth_{arm}_cosine']
                assert truth_score == truth[f'{arm}_truth_cosine']
                better = score > truth_score + 1e-12
                assert better == (arm in competitor['arms']), (arm, competitor['locus'])
                per_arm[arm] += int(better)
            assert len(competitor['identities']) == len(competitor['row_indices']) == 2
            # An empty graph row (no fully-contained syncmer window in its
            # material) contributes no expected mass; a pair with it scores
            # as the other row alone. Counted, not excluded.
            if any(count == 0 for count in competitor['node_counts']):
                empty_usage_competitors += 1
            entries += 1
    assert entries == sum(row['competitor_entries'] for row in rows)
    assert per_arm['nodes'] == sum(row['nodes_higher_competitors'] for row in rows
                                   if row['nodes_higher_competitors'] is not None)
    assert per_arm['combined'] == sum(row['combined_higher_competitors'] for row in rows
                                      if row['combined_higher_competitors'] is not None)
    survivors = [row for row in rows if row['truth_pair_expressible']]
    for arm in ('nodes', 'combined'):
        rank1 = sum(row[f'{arm}_truth_rank'] == 1 for row in survivors)
        rank2 = sum(row[f'{arm}_truth_rank'] <= 2 for row in survivors)
        ranks = {}
        for row in survivors:
            ranks[row[f'{arm}_truth_rank']] = ranks.get(row[f'{arm}_truth_rank'], 0) + 1
        print(f'{component} {arm} arm: rank1 {rank1}/{expected}, rank<=2 {rank2}/{expected}, '
              f'higher-scoring pairs {per_arm[arm]}')
        print(f'  truth ranks by locus: {dict(sorted(ranks.items()))}')
    changed = [row['locus'] for row in survivors
               if row['rank_delta_combined_vs_nodes'] != 0]
    edge_mass = sum(row['observed_edge_mass'] for row in rows)
    node_mass = sum(row['observed_node_mass'] for row in rows)
    print(f'{component}: edge-layer observed mass {edge_mass:.1f} vs node mass {node_mass:.1f}; '
          f'{len(changed)} truth rank(s) changed by the edge layer: {changed}; '
          f'{empty_usage_competitors} competitor pair(s) involve an empty-usage row')

    # --- QUAL receipts phase (material-class product semantics) ---
    qprefix = D / f'run-cosine-graph-qual-pilot-{component}'
    assert Path(f'{qprefix}.done').exists()
    assert Path(f'{qprefix}.exit').read_text().strip() == '0'
    qrows = [json.loads(line) for line in (D / f'cosine-graph-qual-{component}.jsonl').open()]
    assert len(qrows) == len(rows)
    for prior, row in zip(rows, qrows):
        assert prior['locus'] == row['locus']
        for field in STAGE1_EXACT:
            assert prior[field] == row[field], (component, field, row['locus'])
        for field in STAGE1_FLOAT_ULP:
            old, new = prior[field], row[field]
            if old is None or new is None:
                assert old is None and new is None, (component, field, row['locus'])
            else:
                assert abs(new - old) <= 1e-14 * max(1.0, abs(old)), (component, field, row['locus'])
        # identity_kinds iterates a Rust HashMap: order is not stable across
        # processes, so compare as sorted multisets.
        assert sorted(map(list, prior['identity_kinds'])) == sorted(map(list, row['identity_kinds']))
    tie_stream = {}
    with (D / f'cosine-graph-qual-{component}.jsonl.ties.jsonl').open() as stream:
        for line in stream:
            entry = json.loads(line)
            tie_stream.setdefault(entry['locus'], []).append(entry)
    ties_total = bitexact_total = 0
    for row in qrows:
        eligible = row['eligible_pairs_combined']
        classes = row['qual_called_classes']
        k = row['qual_called_class_count']
        if eligible == 0:
            assert k == 0 and classes == []
            assert row['qual_similarity_total'] is None
            assert row['qual_share'] is None and row['qual'] is None
        else:
            assert k >= 1 and len(classes) == k
            best = row['qual_best_similarity']
            total = row['qual_similarity_total']
            assert best == row['combined_best_cosine']
            assert total >= best
            share = row['qual_share']
            assert share == best / total
            if share == 1.0:
                assert row['qual'] is None and row['qual_unbounded'] is True
            else:
                expected_qual = -10.0 * math.log10(1.0 - share)
                assert row['qual'] is not None
                assert abs(row['qual'] - expected_qual) <= 1e-9 * max(1.0, abs(expected_qual))
                assert row['qual_unbounded'] is False
            physical = 0
            for position, cls in enumerate(classes):
                (a, b) = cls['row_indices']
                assert cls['combined_cosine'] == best
                (ma, mb) = cls['row_member_counts']
                pairs = ma * (ma + 1) // 2 if a == b else ma * mb
                assert cls['physical_pair_members'] == pairs
                physical += pairs
                if position == 0:
                    assert cls['node_distance_to_first_called_class'] == [0, 0]
                    assert cls['edge_distance_to_first_called_class'] == [0, 0]
                    assert cls['node_differing_observed_mass_to_first_called'] == 0.0
                    assert cls['edge_differing_observed_mass_to_first_called'] == 0.0
                else:
                    # Fail-closed audit mirroring the truth-window stream: a
                    # called-set tie between distinct classes is expected
                    # only with zero observed mass on the differing segments.
                    assert cls['combined_cosine'] == best
                    assert cls['node_differing_observed_mass_to_first_called'] == 0.0
                    assert cls['edge_differing_observed_mass_to_first_called'] == 0.0
                if row['truth_pair_expressible']:
                    assert cls['is_truth_class'] == (
                        sorted(cls['row_indices']) == sorted(row['truth_graph_rows']))
            assert row['qual_called_physical_pairs'] == physical
            truth_score = row['combined_truth_cosine']
            if truth_score is None:
                expected_in = False
            else:
                expected_in = any(cls['is_truth_class'] for cls in classes)
                # The called set is every class at the bit-identical maximum,
                # so truth-in-called-set is exactly truth score == best.
                assert expected_in == (truth_score == best)
            if row['truth_pair_expressible']:
                assert row['qual_truth_in_called_set'] == expected_in
            else:
                assert row['qual_truth_in_called_set'] is None
            # Derived bound: for a k-way class tie, share <= 1/k (the k tied
            # classes alone contribute k*best to the total), hence
            # Q <= -10*log10(1 - 1/k), equality iff the tie holds all mass.
            if k >= 2:
                bound = -10.0 * math.log10(1.0 - 1.0 / k)
                assert row['qual'] is not None
                assert row['qual'] <= bound + 1e-9, (component, row['locus'])
        # Tie-evidence stream reconciliation (the stage-1 1e-12 window around
        # the truth class, including the truth class itself).
        entries_here = tie_stream.get(row['locus'], [])
        assert row['combined_ties_with_truth_streamed'] == len(entries_here)
        ties_total += len(entries_here)
        expected_count = row['combined_truth_tied_pairs']
        assert len(entries_here) == (expected_count or 0)
        if expected_count:
            truth_score = row['combined_truth_cosine']
            truth_rows = sorted(row['truth_graph_rows'])
            truth_entries = 0
            for entry in entries_here:
                assert abs(entry['combined_cosine'] - truth_score) <= 1e-12
                assert entry['truth_combined_cosine'] == truth_score
                assert entry['combined_bit_exact'] == (entry['combined_ulp_delta'] == 0)
                if entry['combined_bit_exact']:
                    bitexact_total += 1
                    assert entry['combined_cosine'] == truth_score
                    # Fail-closed audit: a bit-exact score tie between
                    # DISTINCT classes is only expected when every differing
                    # segment carries zero observed mass; a tie with observed
                    # mass on the differing segments would need exact f64
                    # cancellation and must be flagged, not accepted.
                    if not entry['is_truth_class']:
                        assert entry['node_differing_observed_mass_vs_truth'] == 0.0
                        assert entry['edge_differing_observed_mass_vs_truth'] == 0.0
                if entry['is_truth_class']:
                    truth_entries += 1
                    assert sorted(entry['row_indices']) == truth_rows
                    assert entry['node_distance_vs_truth_class'] == [0, 0]
                    assert entry['edge_distance_vs_truth_class'] == [0, 0]
                    assert entry['combined_ulp_delta'] == 0
            assert truth_entries == 1
    qsurvivors = [row for row in qrows if row['truth_pair_expressible']]
    print(f'{component} QUAL table (material classes, combined arm), truth-expressible loci:')
    for row in qsurvivors:
        qual = 'inf' if row['qual_unbounded'] else ('%.2f' % row['qual'] if row['qual'] is not None else 'null')
        print(f"  locus {row['locus']}: called classes {row['qual_called_class_count']}, "
              f"physical pairs in called set {row['qual_called_physical_pairs']}, "
              f"truth-in-called-set {'YES' if row['qual_truth_in_called_set'] else 'no'}, "
              f"share {row['qual_share']:.6f}, QUAL {qual}")
    finite_quals = sorted(row['qual'] for row in qsurvivors if row['qual'] is not None)
    unbounded = sum(1 for row in qsurvivors if row['qual_unbounded'])
    in_set = sum(1 for row in qsurvivors if row['qual_truth_in_called_set'])
    tied_max = sum(1 for row in qsurvivors if row['qual_called_class_count'] >= 2)
    if finite_quals:
        median = finite_quals[len(finite_quals) // 2] if len(finite_quals) % 2 else \
            0.5 * (finite_quals[len(finite_quals) // 2 - 1] + finite_quals[len(finite_quals) // 2])
        summary = (f'min {finite_quals[0]:.2f}, median {median:.2f}, max {finite_quals[-1]:.2f}')
    else:
        summary = 'none finite'
    print(f'{component} QUAL summary: truth-in-called-set {in_set}/{expected}, '
          f'loci with >=2 called classes {tied_max}/{expected}, '
          f'unbounded QUAL {unbounded}/{expected}; finite QUAL distribution: {summary}')
    print(f'{component} tie evidence: {ties_total} streamed entries within the 1e-12 window '
          f'of truth, {bitexact_total} bit-exact')
print('all checks passed')
