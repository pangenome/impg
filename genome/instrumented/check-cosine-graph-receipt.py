#!/usr/bin/env python3
"""Independently reconcile the graph-space condensation pilot receipts.

Checks, per component (chrMT, chrI), fail-closed and without recomputing
scores: receipt markers; locus alignment with the balanced Stage-1 census and
the row-space pilot; pure-route physical domain identity; truth survival masks;
exhaustive pair-count bounds for BOTH arms (nodes-only and nodes+edges);
every streamed competitor's arm flags, scores and truth references; the
path-continuity audit (the chain-junction channel's measured emptiness); and
the arm-delta summary (the edge layer's contribution as separate numbers).
"""
import json
from pathlib import Path

D = Path('/home/erikg/yeast/genome-balanced-diploid-validation-20260930')
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
