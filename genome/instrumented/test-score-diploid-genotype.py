#!/usr/bin/env python3
"""Regression tests for the attested per-locus injective diploid yardstick."""
import importlib.util
import unittest
from pathlib import Path
from types import SimpleNamespace

spec = importlib.util.spec_from_file_location(
    'score_diploid_genotype', Path(__file__).with_name('score-diploid-genotype.py'))
score = importlib.util.module_from_spec(spec)
spec.loader.exec_module(score)


class StubMapper:
    def __init__(self, windows):
        self.windows = windows

    def map(self, query):
        target = self.windows.get(query)
        if target is None:
            return []
        return [SimpleNamespace(q_st=0, q_en=len(query), r_st=target[0], r_en=target[1],
                                is_primary=True, mapq=60,
                                strand=target[2] if len(target) == 3 else 1)]


class InjectiveGenotypeDistanceTest(unittest.TestCase):
    def test_two_calls_cannot_both_claim_the_same_truth_allele(self):
        row = score.match_genotype(('AAAA', 'AAAA'), ('AAAA', 'CCCC'))
        self.assertEqual(row['error_bp'], 4)
        self.assertEqual(row['truth_bp'], 8)
        self.assertEqual(row['direct_error_bp'], row['swapped_error_bp'])
        self.assertFalse(row['assignment_identifiable'])
        self.assertEqual(score.match_genotype(('CCCC', 'AAAA'), ('AAAA', 'CCCC'))['error_bp'], 0)

    def test_missing_call_pays_its_assigned_truths_full_length(self):
        row = score.match_genotype(('AAAA', None), ('AAAA', 'CCCCCC'))
        self.assertEqual(row['error_bp'], 6)
        self.assertEqual(row['direct_error_bp'], 6)
        self.assertEqual(row['swapped_error_bp'], 10)
        self.assertEqual(row['called_slot_error_bp'], [0, 6])

    def test_attested_loci_switch_assignment_and_unattested_locus_is_bracketed(self):
        first = 'A' * 100
        second = 'C' * 100
        truths = (first + first + first, second + second + second)
        loci = [{'locus': index, 'axis_interval': [index * 100, (index + 1) * 100],
                 'diplotype': {'dosage': [{'dosage': 2}]}}
                for index in range(3)]
        source = {0: first, 1: second}
        route = [[{'segments': [{'source': index, 'start': 0, 'end': 100, 'reverse': False}]}, None]
                 for index in (0, 1, 0)]
        # First two orthologs are attested. The third is NOT a zero-length
        # truth allele even though the call there resembles S288C perfectly.
        mapper = StubMapper({first: (0, 100)})
        # Identical S288C windows cannot distinguish the third by sequence;
        # use a distinct window to test the correspondence bracket.
        truths = (first + first + 'G' * 100, truths[1])
        result = score.score_loci(loci, route, truths, source, mapper)
        self.assertEqual(result['assessable_loci'], 2)
        self.assertEqual(result['bracketed_loci'], 1)
        self.assertEqual(result['bracketed_S288C_axis_bp'], 100)
        self.assertEqual(result['assignment_track'], ['direct', 'swapped', None])
        self.assertEqual(result['assignment_transitions'], 1)
        self.assertEqual(result['informative_adjacent_boundaries'], 1)
        self.assertEqual(result['error_bp'], 200)
        self.assertEqual(result['truth_bp'], 400)
        self.assertIsNone(result['rows'][2]['truth_allele_bp']['SK1'])
        self.assertNotIn('error_bp', result['rows'][2])

    def test_reverse_strand_ortholog_is_oriented_before_genotype_matching(self):
        first = 'A' * 100
        # Genomic SK1 C*100 reverses into the G*100 local allele.
        mapper = StubMapper({first: (0, 100, -1)})
        loci = [{'locus': 0, 'axis_interval': [0, 100],
                 'diplotype': {'dosage': [{'dosage': 2}]}}]
        route = [[{'segments': [{'source': 1, 'start': 0, 'end': 100}]}, None]]
        result = score.score_loci(loci, route, (first, 'C' * 100),
                                  {1: 'G' * 100}, mapper)
        row = result['rows'][0]
        self.assertEqual(row['SK1_truth_strand'], -1)
        self.assertEqual(row['assignment'], 'swapped')
        self.assertEqual(row['called_slot_error_bp'], [0, 100])
        self.assertEqual(result['error_bp'], 100)


if __name__ == '__main__':
    unittest.main()
