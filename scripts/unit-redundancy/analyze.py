#!/usr/bin/env python3
"""Offline unit-redundancy analysis of one ingredients JSONL file.

Per row (one locus): quantify the two exact reductions over the
ll_matrix's UNIT dimension (columns):

  1. constant-vector units (the column's bits all equal across folds)
  2. identical-vector units (bit-identical full columns, exact grouping)

Nothing is rounded: grouping is by the exact f64 bit pattern of the whole
column (np.void view compare) — the same predicate the codebase's
ll_bits_equal uses.
"""
import json
import re
import sys
import numpy as np

PATH = sys.argv[1]
OUT = sys.argv[2]


def load_row(line):
    i = line.index(b'"ll_matrix":')
    j = line.index(b',"locus":')
    sub = line[i + len(b'"ll_matrix":'):j]
    folds = sub.count(b'[') - 1
    flat = np.fromstring(
        sub.replace(b'[', b' ').replace(b']', b' ').replace(b',', b' '),
        sep=' ')
    u = line.index(b'"units":[')
    counts = np.array([int(x) for x in
                       re.findall(rb'"count":(\d+)', line[u:])],
                      dtype=np.int64)
    locus = int(re.search(rb',"locus":(\d+)', line).group(1))
    partition = int(re.search(rb',"partition":(\d+)', line).group(1))
    if counts.size == 0:
        return locus, partition, folds, counts, np.zeros((folds, 0))
    assert flat.size == folds * counts.size, (
        f"locus {locus}: {flat.size} != {folds}*{counts.size}")
    return locus, partition, folds, counts, flat.reshape(folds, counts.size)


def degenerate(locus, partition, folds, counts):
    return {
        'line': -1, 'locus': locus, 'partition': partition,
        'folds': int(folds), 'units': int(counts.size),
        'cells': 0, 'unit_count_mass': int(counts.sum()),
        'distinct_rows': folds, 'x_row_memo': 1.0,
        'inner_terms_rowmemo': 0, 'inner_terms_rowmemo_units': 0,
        'x_combined': 0.0, 'nonconst_units': 0, 'x_nonconst_dedup': 0.0,
        'const_value_mass_all': {},
        'distinct_vectors': 0, 'const_groups': 0, 'const_units': 0,
        'const_units_pct': 0.0, 'const_values': 0,
        'const_value_hist': {}, 'const_value_mass': {},
        'nonconst_groups': 0, 'mult_hist_all': {},
        'mult_hist_nonconst': {}, 'max_multiplicity': 0,
        'pairs_incl_self': folds * (folds + 1) // 2,
        'inner_terms_full': 0, 'inner_terms_reduced': 0,
        'x_unit_term': 0.0, 'x_identical_only': 0.0, 'empty': True,
    }


def group_columns(M):
    """Exact bit grouping of M's columns. -> (n_distinct, first_idx,
    group_size, inverse, is_const_group, const_value)"""
    folds = M.shape[0]
    Mt = np.ascontiguousarray(M.T)
    colview = Mt.view(np.dtype((np.void, folds * 8))).ravel()
    uniq, first_idx, inv, group_size = np.unique(
        colview, return_index=True, return_inverse=True, return_counts=True)
    bits = Mt.view(np.uint64)
    reps = bits[first_idx]                       # (n_distinct, folds)
    is_const = reps.min(axis=1) == reps.max(axis=1)
    const_value = np.zeros(reps.shape[0], dtype=np.float64)
    const_value.view(np.uint64)[:] = np.where(
        is_const, reps[:, 0], np.uint64(0))
    return (uniq.size, first_idx, group_size, inv, is_const, const_value, Mt)


rows = []
with open(PATH, "rb") as f:
    for lineno, line in enumerate(f):
        locus, partition, folds, counts, M = load_row(line)
        if counts.size == 0:
            r = degenerate(locus, partition, folds, counts)
            r['line'] = lineno
            rows.append(r)
            print(f"line {lineno:2d} locus {locus:2d}: F={folds:4d} "
                  f"U=0 (empty matrix, no unit term)", flush=True)
            continue
        (n_distinct, first_idx, group_size, inv, is_const,
         const_value, Mt) = group_columns(M)

        const_group_idx = np.nonzero(is_const)[0]
        nonconst_group_idx = np.nonzero(~is_const)[0]
        const_units = int(group_size[const_group_idx].sum())

        # count-mass per group (units sharing a bit-identical column)
        mass = np.bincount(inv, weights=counts.astype(np.float64),
                           minlength=n_distinct)

        # distinct constant VALUES: each is its own vector group already,
        # but report the value's own multiplicity (units at that value)
        val_hist = {}
        val_mass = {}
        for g in const_group_idx.tolist():
            key = repr(float(const_value[g]))
            val_hist[key] = val_hist.get(key, 0) + int(group_size[g])
            val_mass[key] = val_mass.get(key, 0.0) + float(mass[g])

        def hist_of(idx):
            h = {}
            for g in idx.tolist():
                s = int(group_size[g])
                h[s] = h.get(s, 0) + 1
            return dict(sorted(h.items()))

        # --- fold-ROW dedup (the already-implemented row memo) ----------
        rowview = np.ascontiguousarray(M).view(
            np.dtype((np.void, counts.size * 8))).ravel()
        distinct_rows = int(np.unique(rowview).size)

        term_full = folds * (folds + 1) // 2
        term_reduced = n_distinct - const_group_idx.size + 1

        rows.append({
            'line': lineno,
            'locus': locus,
            'partition': partition,
            'folds': int(folds),
            'units': int(counts.size),
            'cells': int(folds * counts.size),
            'unit_count_mass': int(counts.sum()),
            'distinct_vectors': int(n_distinct),
            'const_groups': int(const_group_idx.size),
            'const_units': const_units,
            'const_units_pct': 100.0 * const_units / counts.size,
            'const_values': len(val_hist),
            'const_value_hist': dict(sorted(val_hist.items(),
                                            key=lambda kv: -kv[1])[:12]),
            'const_value_mass': dict(sorted(val_mass.items(),
                                            key=lambda kv: -kv[1])[:12]),
            'nonconst_groups': int(nonconst_group_idx.size),
            'mult_hist_all': hist_of(np.arange(n_distinct)),
            'mult_hist_nonconst': hist_of(nonconst_group_idx),
            'max_multiplicity': int(group_size.max()),
            'pairs_incl_self': int(term_full),
            'inner_terms_full': int(term_full * counts.size),
            'inner_terms_reduced': int(term_full * term_reduced),
            'x_unit_term': float(counts.size / term_reduced),
            'x_identical_only': float(counts.size / n_distinct),
            'distinct_rows': distinct_rows,
            'x_row_memo': float(term_full /
                                (distinct_rows * (distinct_rows + 1) // 2)),
            'inner_terms_rowmemo': int(
                (distinct_rows * (distinct_rows + 1) // 2) * counts.size),
            'inner_terms_rowmemo_units': int(
                (distinct_rows * (distinct_rows + 1) // 2) * term_reduced),
            'x_combined': float(term_full * counts.size /
                                ((distinct_rows * (distinct_rows + 1) // 2)
                                 * term_reduced)),
            'nonconst_units': int(counts.size - const_units),
            'x_nonconst_dedup': float((counts.size - const_units) /
                                      max(1, nonconst_group_idx.size)),
            'const_value_mass_all': {repr(float(const_value[g])):
                                     [int(group_size[g]), float(mass[g])]
                                     for g in const_group_idx.tolist()},
        })
        print(f"line {lineno:2d} locus {locus:2d}: F={folds:4d} "
              f"U={counts.size:6d} distinct={n_distinct:6d} "
              f"const_units={const_units:6d} "
              f"X={counts.size/term_reduced:8.3f}", flush=True)
        del line, M, Mt

with open(OUT, "w") as f:
    json.dump(rows, f, indent=1)
print("wrote", OUT)
