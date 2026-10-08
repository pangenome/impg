#!/usr/bin/env python3
"""Source-named regression gates for global two-slot diploid configurations."""
import importlib.util
import random
import unittest
from pathlib import Path

import mappy

spec = importlib.util.spec_from_file_location(
    'score_diploid_distance', Path(__file__).with_name('score-diploid-distance.py'))
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)


def truth_pair():
    rng = random.Random(20260930)
    first = ''.join(rng.choices('ACGT', k=6000))
    second = list(first)
    for index in range(0, len(first), 13):
        second[index] = 'ACGT'[('ACGT'.index(first[index]) + 1) % 4]
    return first, ''.join(second)


class FourWayDiploidDistanceTest(unittest.TestCase):
    def setUp(self):
        self.truths = truth_pair()
        self.aligners = [mappy.Aligner(seq=seq, preset='asm5') for seq in self.truths]
        self.assertTrue(all(self.aligners))

    def test_truth_self_and_swapped_slots_have_zero_sequence_distance(self):
        for inferred in [
            [(0, 0, self.truths[0]), (0, 1, self.truths[1])],
            [(0, 0, self.truths[1]), (0, 1, self.truths[0])],
        ]:
            result = score.best_configuration(inferred, self.truths, self.aligners)
            self.assertEqual(result['error_bp'], 0)
            self.assertEqual(result['unaligned_query_bp'], 0)
            self.assertEqual(result['correct_truth_bp'], {'S288C': 6000, 'SK1': 6000})
            self.assertEqual(result['chosen_configuration'], 0 if inferred[0][2] == self.truths[0] else 1)
            self.assertGreater(result['configurations'][1 - result['chosen_configuration']]['error_bp'], 0)

    def test_missing_slot_is_not_free_accuracy(self):
        result = score.best_configuration([(0, 0, self.truths[0])], self.truths, self.aligners)
        self.assertEqual(result['distance'], 0.5)
        self.assertEqual(result['chosen_configuration'], 0)
        self.assertEqual(result['uncovered_truth_bp']['SK1'], 6000)
        self.assertIsNone(result['assigned_truth_concordance']['SK1'])
        self.assertIsNone(result['configurations'][1]['assigned_truth_concordance']['S288C'])
        self.assertEqual(result['events']['two_slot_phase_switch'], [])

    def test_single_mosaic_slot_cannot_cover_both_truth_molecules(self):
        inferred = [(0, 0, self.truths[0][:3000]),
                    (1, 0, self.truths[1][3000:])]
        result = score.best_configuration(inferred, self.truths, self.aligners)
        assigned = result['slot_to_truth'][0]
        other = result['slot_to_truth'][1]
        self.assertEqual(result['covered_truth_bp'][other], 0)
        self.assertEqual(result['uncovered_truth_bp'][other], len(self.truths[0]))
        self.assertGreater(result['other_truth_block_count'], 0)
        self.assertGreater(result['other_truth_query_bp'], 0)
        self.assertEqual(len(result['events']['single_slot_ancestry_or_crossover']), 1)
        self.assertEqual(result['events']['two_slot_phase_switch'], [])
        self.assertGreater(result['error_bp'], len(self.truths[0]))

    def test_reverse_fragment_competes_in_both_global_configurations(self):
        inferred = [(0, 0, score.rc(self.truths[1][1000:5000]))]
        result = score.best_configuration(inferred, self.truths, self.aligners)
        self.assertEqual(result['slot_to_truth'], ['SK1', 'S288C'])
        self.assertEqual(result['unaligned_query_bp'], 0)
        self.assertTrue(any(block['strand'] == -1 and block['truth'] == 'SK1'
                            and block['edits'] == 0 for block in result['assignments']))


if __name__ == '__main__':
    unittest.main()
