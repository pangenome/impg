#!/usr/bin/env python3
"""The degenerate zero-unit locus convention of the locality domain's
E audit and class table (the chrVII L56 class): the bounding property
is vacuous on an empty matrix, and the scorer's empty-product
convention gives every class LL exactly 0.0, the winner the single
class. Assessment-side; these test the checker's own mirrors of the
committed scorer conventions."""
import importlib.util
import sys
import unittest
from pathlib import Path

import numpy as np

_argv = sys.argv
sys.argv = [_argv[0]]  # the checker parses its committed CLI at import
try:
    spec = importlib.util.spec_from_file_location(
        'realign_scoring_checker',
        Path(__file__).with_name('check-realign-scoring.py'))
    checker = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(checker)
finally:
    sys.argv = _argv

# a representative elsewhere branch (the chrVII gate's own E)
E = -17.870035


class DegenerateZeroUnitLocusTest(unittest.TestCase):
    def test_bounding_property_is_vacuous_on_the_empty_matrix(self):
        # chrVII L56: the 15bp contig-end window whose universe holds
        # only the axis fold - one fold, zero units, an EMPTY matrix.
        # matrix.min() would raise on it; the guard must not.
        matrix = np.zeros((1, 0))
        self.assertTrue(checker.e_bounding_property_holds(matrix, E))

    def test_bounding_property_still_binds_nonempty_matrices(self):
        # no behavior change elsewhere: entries at exactly E pass (the
        # no-pin convention), an entry below E still fails.
        at_e = np.full((3, 5), E)
        self.assertTrue(checker.e_bounding_property_holds(at_e, E))
        below = at_e.copy()
        below[1, 2] = E - 1.0
        self.assertFalse(checker.e_bounding_property_holds(below, E))

    def test_zero_unit_locus_class_lls_are_exactly_zero(self):
        # the empty product: every class LL exactly 0.0, every class
        # ties, the winner the single class by flat order.
        for n_folds in (1, 3):
            matrix = np.zeros((n_folds, 0))
            counts = np.zeros((0,))
            table = checker.class_ll_table(matrix, counts)
            self.assertEqual(
                table.shape, (n_folds * (n_folds + 1) // 2,)
            )
            self.assertTrue(all(ll == 0.0 for ll in table))
            self.assertEqual(int(np.argmax(table)), 0)

    def test_class_table_matches_the_direct_computation_when_stocked(self):
        # the extraction is a no-op: the phase-5 arithmetic equals the
        # direct mix_logsumexp sum over unordered fold pairs.
        matrix = np.array(
            [
                [-20.0, -17.5, -30.0],
                [-25.0, E, E],
                [-19.0, -22.0, -40.0],
                [-18.0, -21.0, -26.0],
            ]
        )
        counts = np.array([1.0, 2.0, 1.0])
        table = checker.class_ll_table(matrix, counts)
        expected = []
        for j in range(4):
            for i in range(j + 1):
                total = 0.0
                for u in range(3):
                    total += counts[u] * checker.mix_logsumexp(
                        matrix[i][u], matrix[j][u]
                    )
                expected.append(total)
        self.assertEqual(len(table), len(expected))
        for got, want in zip(table, expected):
            self.assertAlmostEqual(got, want, places=9)

    def test_held_rank1_assertion_not_applicable_at_inexpressible_truth(self):
        # the same degenerate locus's second seam: a before-rank-1
        # locus whose truth pair became INEXPRESSIBLE under the new
        # domains (chrVII L56, the over-join separation class) - the
        # 8N winner-equality assertion is NOT APPLICABLE (None), the
        # regression is already named by the gate itself.
        def fold_member_key(d, index):
            return tuple(
                sorted(
                    (m["path_name"], m["start"], m["end"])
                    for m in d["fold_identities"][index]["members"]
                )
            )

        winner = [{"path_name": "SK1#0#chrVII", "start": 0, "end": 15}]
        held = {
            "truth_folds": None,
            "best_fold_indices": [0],
            "fold_identities": [{"members": winner}],
        }
        self.assertIsNone(
            checker.held_rank1_winner_assertion(held, fold_member_key)
        )
        # no behavior change elsewhere: winner == truth passes,
        # winner != truth fails
        truth = [{"path_name": "SK1#0#chrVII", "start": 0, "end": 15}]
        ok = {
            "truth_folds": [0, 0],
            "best_fold_indices": [0, 0],
            "fold_identities": [{"members": winner}, {"members": truth}],
        }
        self.assertTrue(checker.held_rank1_winner_assertion(ok, fold_member_key))
        rival = [{"path_name": "CENPK#0#chrVII", "start": 1, "end": 16}]
        bad = {
            "truth_folds": [1, 1],
            "best_fold_indices": [0, 0],
            "fold_identities": [{"members": winner}, {"members": rival}],
        }
        self.assertFalse(checker.held_rank1_winner_assertion(bad, fold_member_key))


if __name__ == "__main__":
    unittest.main()
