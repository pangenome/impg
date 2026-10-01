#!/usr/bin/env python3
"""Independently reconcile the two fail-closed full-domain pilot receipts."""
import json
from pathlib import Path

D = Path('/home/erikg/yeast/genome-balanced-diploid-validation-20260930')
for component, expected in [('chrMT', 8), ('chrI', 17)]:
    prefix = D / f'run-cosine-exhaustive-pilot-{component}'
    assert Path(f'{prefix}.done').exists() and Path(f'{prefix}.exit').read_text().strip() == '0'
    old = json.loads((D / f'run-balanced-{component}.log').read_text())['spine']['stage1_sweep']['loci']
    rows = [json.loads(line) for line in (D / f'cosine-exhaustive-{component}.jsonl').open()]
    assert len(rows) == len(old)
    assert all(row['locus'] == previous['locus'] for row, previous in zip(rows, old))
    assert all(row['truth_pair_expressible'] == (previous['truth_pair_in_domain'] == [True, True])
               for row, previous in zip(rows, old))
    assert sum(row['truth_pair_expressible'] for row in rows) == expected
    assert all(row['physical_rows'] == previous['alleles'] for row, previous in zip(rows, old))
    assert all(row['eligible_pairs'] <= row['material_vectors'] * (row['material_vectors'] + 1) // 2
               for row in rows)
    by_locus = {row['locus']: row for row in rows}
    count = 0
    with (D / f'cosine-exhaustive-{component}.jsonl.competitors.jsonl').open() as stream:
        for line in stream:
            competitor = json.loads(line)
            truth = by_locus[competitor['locus']]
            assert truth['truth_pair_expressible']
            assert competitor['cosine'] > competitor['truth_cosine'] + 1e-12
            assert competitor['truth_cosine'] == truth['truth_cosine']
            assert len(competitor['identities']) == len(competitor['material']) == 2
            count += 1
    assert count == sum(row['higher_cosine_competitors'] for row in rows)
    rank1 = sum(row['truth_rank'] == 1 for row in rows if row['truth_pair_expressible'])
    print(f'{component}: {expected}/{expected} truth survivors, {rank1}/{expected} rank1, '
          f'{count} named higher-scoring physical pairs')
