#!/usr/bin/env python3
"""Coverage-vector COSIGT probe arithmetic, not a production call test."""
import importlib.util
import unittest
from pathlib import Path

spec = importlib.util.spec_from_file_location('probe', Path(__file__).with_name('score-cosine-probe.py'))
probe = importlib.util.module_from_spec(spec)
spec.loader.exec_module(probe)


class CosineProbeTest(unittest.TestCase):
    def test_row_factored_pair_cosine_matches_direct_sum(self):
        observed = [15.0, 10.0, 3.0]
        first = [15.0, 15.0, 0.0]
        second = [0.0, 15.0, 15.0]
        dot = lambda a, b: sum(x * y for x, y in zip(a, b))
        pair = [a + b for a, b in zip(first, second)]
        factored = probe.cosine(dot(observed, first) + dot(observed, second),
                                dot(observed, observed),
                                dot(first, first) + dot(second, second) + 2 * dot(first, second))
        direct = probe.cosine(dot(observed, pair), dot(observed, observed), dot(pair, pair))
        self.assertAlmostEqual(factored, direct)

    def test_placement_and_feature_arms_penalize_wrong_material(self):
        row = {
            'locus': 4, 'record_count': 3, 'depth_per_copy': 15,
            'coverage_bins': [[0, 0, 10, 15], [1, 0, 10, 15], [2, 0, 10, 0]],
            'observed_feature_mass': [[[1], 15], [[2], 15], [[3], 0]],
            'alleles': [
                {'allele': i, 'material_intervals': [[i, 0, 10]],
                 'feature_profile': [[[i + 1], 1]]} for i in range(3)],
            'candidate_pairs': [
                {'label': label, 'alleles': pair, 'old_objective_loss': None}
                for label, pair in [('truth', [0, 1]), ('selected', [0, 2]),
                                    ('old_objective_top_1', [0, 2])]],
        }
        result = probe.score_locus(row)
        self.assertEqual((result['coverage_rank'], result['feature_rank']), (1, 1))
        self.assertEqual(result['unique_candidates'], 2)
        self.assertAlmostEqual(result['scores'][0]['coverage_cosine'], 1.0)
        self.assertAlmostEqual(result['scores'][0]['feature_cosine'], 1.0)
        self.assertLess(result['scores'][1]['coverage_cosine'], 1.0)
        self.assertLess(result['scores'][1]['feature_cosine'], 1.0)

    def test_dosage_twice_on_one_path_and_zero_observations(self):
        row = {'locus': 0, 'record_count': 0, 'depth_per_copy': 15,
               'coverage_bins': [[0, 0, 10, 0]], 'observed_feature_mass': [],
               'alleles': [{'allele': 0, 'material_intervals': [[0, 0, 10]],
                            'feature_profile': [[[1], 1]]}],
               'candidate_pairs': [{'label': 'truth', 'alleles': [0, 0],
                                    'old_objective_loss': None}]}
        result = probe.score_locus(row)
        self.assertIsNone(result['coverage_rank'])
        self.assertIsNone(result['feature_rank'])


if __name__ == '__main__':
    unittest.main()
