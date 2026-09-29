#!/usr/bin/env python3
"""Assessment-only diploid score: full unique homolog mappings, no proportional coordinates.

Usage: score-diploid.py chrMT smoke2 (writes JSON on stdout). An absent
second inferred slot is candidate-only, not a scored diploid error.
"""
import json
import subprocess
import sys
from pathlib import Path

import edlib
import mappy

BASE = Path('/home/erikg/yeast/genome-diploid-validation-20260929')
COMP, TAG = sys.argv[1:3]
truth = {name.rsplit(':', 1)[0]: seq for name, seq, *_ in
         mappy.fastx_read(str(BASE / 'private-truth/truth-source-paths.fa'))}
source = {}
for fasta in (BASE / 'sources-all.fa', Path('/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/genome-genotyping-v1-item1/sources-core.fa')):
    source.update({int(name.split(':', 1)[0]): seq for name, seq, *_ in mappy.fastx_read(str(fasta))})
names = {}
for line in Path('/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names').read_text().splitlines():
    fields = line.split('\t')
    if fields[0].isdigit():
        names[int(fields[0])] = fields[1]
        if fields[1] in truth:
            source[int(fields[0])] = truth[fields[1]]
assert all(len(seq) > 0 for seq in source.values())
s288c = truth['S288C#0#' + COMP]
sk1 = truth['SK1#0#' + COMP]
sample_truth = json.loads((BASE / 'private-truth/truth-diploid.json').read_text())
doses = {slot: next(row['dose'] for row in sample_truth[key] if row['component'] == COMP)
         for slot, key in (('S288C', 'slot1'), ('SK1', 'slot2'))}
truth_dosage = {slot: len(doses) * dose / sum(doses.values()) for slot, dose in doses.items()}
run = json.loads((BASE / f'run-{TAG}-{COMP}.log').read_text())
rescore = run['spine']['phasing']['rescore']
assert rescore['references'][0]['label'] == 'truth'
assert rescore['references'][0]['m1_oracle_rescore'] is not None
for name, sequence in (('S288C', s288c), ('SK1', sk1)):
    recorded = json.loads((BASE / 'routes' / f'{name}-{COMP}.json').read_text())
    assert len(recorded['segments']) == 1
    segment = recorded['segments'][0]
    assert not segment['reverse'] and segment['start'] == 0 and segment['end'] == len(sequence)
    assert source[segment['source']] == sequence
loci = run['spine']['genotype_calls']['loci']
route = rescore['selected_route']
assert len(loci) == len(route)
mapper = mappy.Aligner(seq=sk1, preset='asm5')
assert mapper


def rc(seq):
    return seq.translate(str.maketrans('ACGTacgt', 'TGCAtgca'))[::-1]


missing = sorted({part['source'] for entry in route
                  for slot in (entry if isinstance(entry, list) else [entry])
                  if isinstance(slot, dict) for part in slot.get('segments', [])} - source.keys())
if missing:
    extra = BASE / f'sources-score-{COMP}.fa'
    subprocess.run(['/home/erikg/impg-genome-inference/target/release/examples/genotyping_spell_dump',
                    '/home/erikg/yeast/yeast235.agc', str(extra),
                    *(f'{i}:{names[i]}' for i in missing)], check=True, stdout=subprocess.DEVNULL)
    source.update({int(name.split(':', 1)[0]): seq for name, seq, *_ in mappy.fastx_read(str(extra))})
assert not (set(missing) - source.keys())


def spell(slot):
    if not isinstance(slot, dict):
        return None
    parts = []
    for part in slot.get('segments', []):
        seq = source[part['source']][part['start']:part['end']]
        parts.append(rc(seq) if part.get('reverse') else seq)
    return ''.join(parts) if parts else None


def accuracy(query, target):
    # Global edit distance: do not let semi-global placement hide unmatched ends.
    return 1 - edlib.align(query, target, mode='NW', task='distance')['editDistance'] / max(len(query), len(target))


def full_unique_homolog(window):
    hits = [h for h in mapper.map(window) if h.q_st == 0 and h.q_en == len(window)
            and h.is_primary and h.mapq > 0]
    return hits[0] if len(hits) == 1 else None


rows = []
for locus, inferred_route in zip(loci, route):
    a, b = locus['axis_interval']
    t1 = s288c[a:b]
    hit = full_unique_homolog(t1)
    slots = inferred_route if isinstance(inferred_route, list) else [inferred_route]
    inferred = [s for s in (spell(slot) for slot in slots) if s is not None]
    row = {'locus': locus['locus'], 'axis_interval': [a, b],
           'truth_dosage': truth_dosage,
           'selected_copy_count': len(inferred),
           'called_dosage_classes': [x['dosage'] for x in locus['diplotype']['dosage']],
           'exact_homolog': hit is not None}
    if hit is not None:
        t2 = sk1[hit.r_st:hit.r_en]
        if hit.strand < 0:
            t2 = rc(t2)
        # Gate: self-score the TWO actual truth sequences independently and
        # through the two-slot best-matching calculation, not the called route.
        self_scores = [[accuracy(t, target) for target in (t1, t2)] for t in (t1, t2)]
        row['truth_self_accuracy'] = max(self_scores[0][0] + self_scores[1][1],
                                          self_scores[0][1] + self_scores[1][0]) / 2
        row['truth_self_switches'] = 0
        assert self_scores[0][0] == 1.0 and self_scores[1][1] == 1.0
        assert row['truth_self_accuracy'] == 1.0
        scores = [[accuracy(s, target) for target in (t1, t2)] for s in inferred]
        if len(scores) >= 2:
            direct = scores[0][0] + scores[1][1]
            swapped = scores[0][1] + scores[1][0]
            row['pair_accuracy'] = max(direct, swapped) / 2
            row['assignment'] = 'direct' if direct >= swapped else 'swapped'
        elif len(scores) == 1:
            row['candidate_only_accuracy'] = max(scores[0])
            row['candidate_truth_match'] = ('S288C' if scores[0][0] > scores[0][1]
                                            else 'SK1' if scores[0][1] > scores[0][0]
                                            else 'ambiguous')
            row['candidate_truth_accuracies'] = dict(zip(('S288C', 'SK1'), scores[0]))
            row['candidate_assigned_dosage'] = ({name: (len(doses) if name == row['candidate_truth_match'] else 0)
                                                 for name in doses}
                                                if row['candidate_truth_match'] != 'ambiguous' else None)
    rows.append(row)

exact = [x for x in rows if x['exact_homolog']]
# Stage-1's existential predicate accepts any row expressing the route's
# local pieces, not only the genotype table's exact-piece representative.
truth_domain = [locus['truth_pair_in_domain'] for locus in run['spine']['stage1_sweep']['loci']]
assert len(truth_domain) == len(loci)
assert all(len(flags) == len(doses) for flags in truth_domain)
paired = [x for x in exact if 'pair_accuracy' in x]
candidates = [x for x in exact if 'candidate_only_accuracy' in x]
# Switching is meaningful only across adjacent, fully typed, two-slot windows.
switches = [0, 0]
comparisons = 0
for left, right in zip(rows, rows[1:]):
    if right['locus'] == left['locus'] + 1 and 'assignment' in left and 'assignment' in right:
        comparisons += 1
        if left['assignment'] != right['assignment']:
            switches[0] += 1
            switches[1] += 1
summary = {'component': COMP, 'tag': TAG, 'exit': (BASE / f'run-{TAG}-{COMP}.exit').read_text().strip(),
           'wall_seconds': int((BASE / f'run-{TAG}-{COMP}.wall').read_text()),
           'selected_ploidy': rescore['selected_ploidy'], 'windows': len(rows),
           'truth_in_domain_by_slot': {'S288C': sum(flags[0] for flags in truth_domain),
                                       'SK1': sum(flags[1] for flags in truth_domain)},
           'truth_pair_in_domain_windows': sum(all(flags) for flags in truth_domain),
           'truth_pair_vs_selected_m1_nats': (rescore['references'][0]['m1_oracle_rescore']
                                              - rescore['selected_chain_self_check']['external_m1_oracle']),
           'exact_homolog_windows': len(exact), 'bracketed_homolog_windows': len(rows) - len(exact),
           'truth_self_test': {'accuracy': 1.0 if exact else None, 'switches': 0 if exact else None,
                               'windows': len(exact)},
           'two_slot_windows': len(paired), 'candidate_only_windows': len(candidates),
           'pair_accuracy_attested_only': sum(x['pair_accuracy'] for x in paired) / len(paired) if paired else None,
           'candidate_accuracy_attested_only': sum(x['candidate_only_accuracy'] for x in candidates) / len(candidates) if candidates else None,
           'candidate_truth_match_counts': {slot: sum(x['candidate_truth_match'] == slot for x in candidates)
                                            for slot in ('S288C', 'SK1', 'ambiguous')},
           'switches_per_haplotype': switches if comparisons else None,
           'switch_boundaries': comparisons,
           'truth_dosage': truth_dosage,
           'selected_single_class_dosage': float(len(doses)) if rescore['selected_ploidy'] == 'haploid' else None,
           'rows': rows}
assert summary['exit'] == '0'
print(json.dumps(summary, indent=2))
