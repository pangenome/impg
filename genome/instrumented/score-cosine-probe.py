#!/usr/bin/env python3
"""Assessment-only cosine of two saved-data read-placement diagnostic arms.

Rank is strictly among the materialized truth, selected, eight best unrestricted
and eight best viable old-objective class pairs, not all diplotypes.
"""
import json
import math
import sys
from pathlib import Path


def cosine(dot, observed_norm, expected_norm):
    if observed_norm <= 0 or expected_norm <= 0:
        return None
    return dot / math.sqrt(observed_norm * expected_norm)


def score_locus(row):
    alleles = {entry['allele']: entry for entry in row['alleles']}
    observed = {tuple(feature): mass for feature, mass in row['observed_feature_mass']}
    observed_norm_b = sum(value * value for value in observed.values())
    bins = row['coverage_bins']
    observed_norm_a = sum((hi - lo) * value * value for _, lo, hi, value in bins)
    scores = []
    for candidate in row['candidate_pairs']:
        first, second = (alleles[index] for index in candidate['alleles'])
        intervals = first['material_intervals'] + second['material_intervals']
        dot_a = expected_norm_a = 0.0
        depth = row['depth_per_copy']
        for path, lo, hi, observed_coverage in bins:
            # All candidate material endpoints are bin boundaries. Two copies
            # on one path add at the same base, never overwrite each other.
            dosage = sum(1 for p, start, end in intervals
                         if p == path and start <= lo and hi <= end)
            expected = depth * dosage
            width = hi - lo
            dot_a += width * observed_coverage * expected
            expected_norm_a += width * expected * expected
        profile = {}
        for entry in (first, second):
            for feature, count in entry['feature_profile']:
                key = tuple(feature)
                profile[key] = profile.get(key, 0.0) + depth * count
        dot_b = sum(expected * observed.get(feature, 0.0)
                    for feature, expected in profile.items())
        expected_norm_b = sum(expected * expected for expected in profile.values())
        scores.append({'label': candidate['label'],
                       'alleles': candidate['alleles'],
                       'old_objective_loss': candidate['old_objective_loss'],
                       'coverage_cosine': cosine(dot_a, observed_norm_a, expected_norm_a),
                       'feature_cosine': cosine(dot_b, observed_norm_b, expected_norm_b)})
    # Identical physical allele pairs from several old-objective class rows
    # count once; keep all their labels in the diagnostic for provenance.
    unique = {}
    for item in scores:
        unique.setdefault(tuple(sorted(item['alleles'])), item)
    truth = unique[tuple(sorted(scores[0]['alleles']))]
    result = {'locus': row['locus'], 'record_count': row['record_count'],
              'old_selected_truth_pair': len(scores) > 1 and
                  tuple(sorted(scores[0]['alleles'])) == tuple(sorted(scores[1]['alleles'])),
              'truth_alleles_port_viable': row.get('truth_alleles_port_viable'),
              'truth_classes_port_viable': row.get('truth_classes_port_viable'),
              'unique_candidates': len(unique), 'scores': scores}
    for arm in ('coverage', 'feature'):
        key = f'{arm}_cosine'
        value = truth[key]
        if value is None:
            result[arm + '_rank'] = None
            result[arm + '_tied_best'] = None
            continue
        better = sum(item[key] is not None and item[key] > value + 1e-12
                     for item in unique.values())
        best = max(item[key] for item in unique.values() if item[key] is not None)
        result[arm + '_rank'] = 1 + better
        selected = scores[1][key] if len(scores) > 1 else None
        result[arm + '_truth_minus_selected'] = value - selected if selected is not None else None
        result[arm + '_tied_best'] = better == 0 and sum(
            item[key] is not None and abs(item[key] - best) <= 1e-12
            for item in unique.values()) > 1
    return result


def main(paths):
    rows = []
    total = 0
    for path in paths:
        component = Path(path).stem.rsplit('-', 1)[-1]
        for line in Path(path).open():
            row = json.loads(line)
            if not row['truth_pair_expressible']:
                continue
            total += 1
            scored = score_locus(row)
            scored['component'] = component
            rows.append(scored)
    summary = {'expressible_loci': total,
               'comparison_set': 'truth + selected + eight best unrestricted + eight best viable old-loss class pairs',
               'old_selected_truth_pairs': sum(row['old_selected_truth_pair'] for row in rows),
               'truth_alleles_port_viable': sum(row['truth_alleles_port_viable'] is True for row in rows),
               'truth_classes_port_viable': sum(row['truth_classes_port_viable'] is True for row in rows),
               'loci': rows}
    for arm in ('coverage', 'feature'):
        ranks = [item[arm + '_rank'] for item in rows]
        wrong = [row for row in rows if not row['old_selected_truth_pair']]
        summary[arm] = {
            'wrong_selected_loci': len(wrong),
            'truth_beats_selected_on_wrong': sum(row.get(arm + '_truth_minus_selected') is not None
                                                 and row[arm + '_truth_minus_selected'] > 1e-12 for row in wrong),
            'truth_ties_selected_on_wrong': sum(row.get(arm + '_truth_minus_selected') is not None
                                               and abs(row[arm + '_truth_minus_selected']) <= 1e-12 for row in wrong),
            'rank_1_on_wrong_selected': sum(row[arm + '_rank'] == 1 for row in wrong),
            'rank_at_most_2_on_wrong_selected': sum(row[arm + '_rank'] is not None and row[arm + '_rank'] <= 2 for row in wrong),
            'scorable': sum(rank is not None for rank in ranks),
            'rank_1_including_ties': sum(rank == 1 for rank in ranks),
            'rank_1_unique': sum(rank == 1 and not item[arm + '_tied_best']
                                 for item, rank in zip(rows, ranks)),
            'rank_at_most_2': sum(rank is not None and rank <= 2 for rank in ranks),
            'rank_gt_2': sum(rank is not None and rank > 2 for rank in ranks),
            'unscorable': sum(rank is None for rank in ranks),
        }
    print(json.dumps(summary, indent=2))


if __name__ == '__main__':
    main(sys.argv[1:])
