#!/usr/bin/env python3
"""Assessment-only diploid distance under the best global two-slot assignment.

Usage: score-diploid-distance.py chrMT smoke2. Both slot-to-truth permutations
are scored separately; other-truth crossover blocks are diagnostics, not
coverage credit for a second molecule. Production outputs are untouched.
"""
import bisect
import json
import re
import subprocess
import sys
from pathlib import Path

import edlib
import mappy

BASE = Path('/home/erikg/yeast/genome-diploid-validation-20260929')
CORE = Path('/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/genome-genotyping-v1-item1/sources-core.fa')
NAMES = Path('/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng.names')
TRUTH_NAMES = ('S288C', 'SK1')
TRANSLATE = str.maketrans('ACGTacgt', 'TGCAtgca')


def rc(seq):
    return seq.translate(TRANSLATE)[::-1]


def cigar_ops(cigar):
    return [(int(n), op) for n, op in re.findall(r'(\d+)([=XID])', cigar)]


def candidate(query, truth_index, truth, hit):
    q_start, q_end = hit.q_st, hit.q_en
    t_start, t_end = hit.r_st, hit.r_en
    target = truth[t_start:t_end]
    fragment = query[q_start:q_end]
    if hit.strand < 0:
        fragment = rc(fragment)
    if fragment == target:
        ops, edits = [(len(fragment), '=')], 0
    else:
        result = edlib.align(fragment, target, mode='NW', task='path')
        if result['editDistance'] < 0 or not result['cigar']:
            raise ValueError('local block failed global edit alignment')
        ops, edits = cigar_ops(result['cigar']), result['editDistance']
    reference_consumed = sum(n for n, op in ops if op in '=XD')
    query_consumed = sum(n for n, op in ops if op in '=XI')
    assert reference_consumed == t_end - t_start and query_consumed == q_end - q_start
    correct = sum(n for n, op in ops if op == '=')
    inserted = sum(n for n, op in ops if op == 'I')
    # Local scheduling gain before cross-window truth-position deduplication:
    # a mapped query block avoids an unaligned-query charge and each matched
    # truth base avoids an omission. The final union is scored separately.
    gain = correct + (q_end - q_start) - inserted
    return {'q_start': q_start, 'q_end': q_end, 't_start': t_start, 't_end': t_end,
            'truth_index': truth_index, 'strand': hit.strand, 'ops': ops,
            'edits': edits, 'correct': correct, 'inserted': inserted, 'gain': gain}


def select_blocks(query, aligners, truths, allowed_truth_indices=(0, 1)):
    # Assembly mappers may omit a perfectly identical full chromosome over
    # repetitive rDNA (chrXII), even when the two literal strings are equal.
    # Exact identity is a valid alignment certificate without seed mapping.
    for truth_index in allowed_truth_indices:
        truth = truths[truth_index]
        if query == truth:
            length = len(query)
            return [{'q_start': 0, 'q_end': length, 't_start': 0, 't_end': length,
                     'truth_index': truth_index, 'strand': 1,
                     'ops': [(length, '=')], 'edits': 0, 'correct': length,
                     'inserted': 0, 'gain': 2 * length}], 1
    options = []
    for truth_index in allowed_truth_indices:
        aligner, truth = aligners[truth_index], truths[truth_index]
        for hit in aligner.map(query):
            if hit.q_en > hit.q_st and hit.r_en > hit.r_st:
                options.append(candidate(query, truth_index, truth, hit))
    options.sort(key=lambda c: (c['q_end'], c['q_start'], c['truth_index'],
                                c['t_start'], c['strand']))
    # Weighted interval scheduling over the reported mapper blocks. The
    # distance call restricts this to ONE assigned truth per inferred slot;
    # the separate diagnostic call may identify other-truth ancestry.
    ends = [c['q_end'] for c in options]
    best = [0] * (len(options) + 1)
    pred = [0] * len(options)
    take = [False] * len(options)
    for i, c in enumerate(options):
        pred[i] = bisect.bisect_right(ends, c['q_start'], 0, i)
        chosen = best[pred[i]] + c['gain']
        if chosen > best[i]:
            best[i + 1] = chosen
            take[i] = True
        else:
            best[i + 1] = best[i]
    selected = []
    i = len(options)
    while i:
        if take[i - 1]:
            selected.append(options[i - 1])
            i = pred[i - 1]
        else:
            i -= 1
    selected.sort(key=lambda c: c['q_start'])
    assert sum(c['q_end'] - c['q_start'] for c in selected) <= len(query)
    return selected, len(options)


def score_blocks(inferred, truths, aligners, slot_to_truth=(0, 1), ancestry=False):
    # Each inferred entry is (locus, slot, spelled sequence). A whole slot
    # has exactly ONE assigned truth for this configuration, including when
    # it is empty. Physical truth bases receive credit once across windows.
    correct_positions = [bytearray(len(seq)) for seq in truths]
    covered_positions = [bytearray(len(seq)) for seq in truths]
    inserted = unaligned_query = duplicated_truth = candidate_hits = 0
    assignments = []
    other_truth_blocks = []
    ancestry_assignments = []
    for locus, slot, sequence in inferred:
        assert slot in (0, 1)
        assigned = slot_to_truth[slot]
        blocks, count = select_blocks(sequence, aligners, truths, (assigned,))
        candidate_hits += count
        if ancestry:
            diagnostic, _ = select_blocks(sequence, aligners, truths)
            for block in diagnostic:
                ancestry_assignments.append((locus, slot, block))
                if block['truth_index'] != assigned:
                    other_truth_blocks.append({
                        'locus': locus, 'slot': slot,
                        'assigned_truth': TRUTH_NAMES[assigned],
                        'other_truth': TRUTH_NAMES[block['truth_index']],
                        'query_interval': [block['q_start'], block['q_end']],
                        'other_truth_interval': [block['t_start'], block['t_end']],
                        'strand': block['strand'], 'correct_bp': block['correct'],
                        'edits': block['edits'],
                        'charged_only_by_assigned_truth_comparison': True})
        unaligned_query += len(sequence) - sum(c['q_end'] - c['q_start'] for c in blocks)
        for block in blocks:
            target = block['truth_index']
            pos = block['t_start']
            for n, op in block['ops']:
                if op == '=':
                    for index in range(pos, pos + n):
                        if correct_positions[target][index]:
                            duplicated_truth += 1
                        correct_positions[target][index] = 1
                        covered_positions[target][index] = 1
                    pos += n
                elif op == 'X' or op == 'D':
                    covered_positions[target][pos:pos + n] = b'\x01' * n
                    pos += n
                elif op == 'I':
                    inserted += n
            assert pos == block['t_end']
            assignments.append({'locus': locus, 'slot': slot, 'truth': TRUTH_NAMES[target],
                                'query_interval': [block['q_start'], block['q_end']],
                                'truth_interval': [block['t_start'], block['t_end']],
                                'strand': block['strand'], 'correct_bp': block['correct'],
                                'edits': block['edits']})
    correct = [sum(positions) for positions in correct_positions]
    covered = [sum(positions) for positions in covered_positions]
    total_truth = sum(map(len, truths))
    # All unaccounted ASSIGNED truth bases cost one: mismatches, deletions
    # and the empty second molecule alike. Query-only bases are additional
    # errors. Other-truth diagnostic blocks never contribute coverage.
    # Selected routes are window-indexed, not physical assemblies. Overlap
    # of their truth bases is procedural and receives no second error charge;
    # the haploid scorer also unions truth positions once per material.
    errors = total_truth - sum(correct) + inserted + unaligned_query
    events = {'single_slot_ancestry_or_crossover': [], 'two_slot_phase_switch': []}
    if ancestry:
        for slot in (0, 1):
            ordered = sorted(((locus, block) for locus, s, block in ancestry_assignments
                              if s == slot), key=lambda entry: (entry[0], entry[1]['q_start']))
            for (prev_locus, prev), (curr_locus, curr) in zip(ordered, ordered[1:]):
                if prev['truth_index'] != curr['truth_index']:
                    # A truth-label transition in an isolated slot is ancestry,
                    # never a measured switch between two inferred molecules.
                    events['single_slot_ancestry_or_crossover'].append({
                        'slot': slot, 'before_locus': prev_locus, 'after_locus': curr_locus,
                        'from': TRUTH_NAMES[prev['truth_index']],
                        'to': TRUTH_NAMES[curr['truth_index']]})
    return {'distance': errors / total_truth, 'error_bp': errors,
            'truth_bp': total_truth, 'slot_to_truth': [TRUTH_NAMES[i] for i in slot_to_truth],
            'correct_truth_bp': dict(zip(TRUTH_NAMES, correct)),
            'covered_truth_bp': dict(zip(TRUTH_NAMES, covered)),
            'uncovered_truth_bp': dict(zip(TRUTH_NAMES,
                                          [len(seq) - bp for seq, bp in zip(truths, covered)])),
            'inserted_query_bp': inserted, 'unaligned_query_bp': unaligned_query,
            'duplicated_truth_bp': duplicated_truth, 'candidate_alignment_blocks': candidate_hits,
            'selected_alignment_blocks': len(assignments), 'events': events,
            'other_truth_blocks': other_truth_blocks,
            'other_truth_block_count': len(other_truth_blocks),
            'other_truth_query_bp': sum(block['query_interval'][1] - block['query_interval'][0]
                                        for block in other_truth_blocks),
            'assignments': assignments}


def best_configuration(inferred, truths, aligners):
    configurations = [score_blocks(inferred, truths, aligners, mapping)
                      for mapping in ((0, 1), (1, 0))]
    chosen = min(range(2), key=lambda index: (configurations[index]['error_bp'], index))
    scored = score_blocks(inferred, truths, aligners,
                          (0, 1) if chosen == 0 else (1, 0), ancestry=True)
    assert scored['error_bp'] == configurations[chosen]['error_bp']
    scored['chosen_configuration'] = chosen
    configuration_rows = []
    for result in configurations:
        row = {key: result[key] for key in
               ('slot_to_truth', 'distance', 'error_bp', 'truth_bp',
                'correct_truth_bp', 'covered_truth_bp', 'uncovered_truth_bp',
                'inserted_query_bp', 'unaligned_query_bp')}
        row['assigned_truth_concordance'] = {
            name: result['correct_truth_bp'][name] / result['covered_truth_bp'][name]
            if result['covered_truth_bp'][name] else None for name in TRUTH_NAMES}
        configuration_rows.append(row)
    scored['configurations'] = configuration_rows
    scored['assigned_truth_concordance'] = configuration_rows[chosen]['assigned_truth_concordance']
    scored['other_truth_charge_treatment'] = (
        'other-truth alignments identify ancestry only; they never cover the other truth '
        'or credit a second molecule; differences and unaligned query against the '
        'globally assigned truth remain charged by the fixed-configuration distance')
    return scored


def main(component, tag):
    truth_fasta = {name.rsplit(':', 1)[0]: seq for name, seq, *_ in
                   mappy.fastx_read(str(BASE / 'private-truth/truth-source-paths.fa'))}
    truths = tuple(truth_fasta[f'{name}#0#{component}'] for name in TRUTH_NAMES)
    assert all(truths)
    source = {}
    for fasta in (BASE / 'sources-all.fa', CORE):
        source.update({int(name.split(':', 1)[0]): seq for name, seq, *_ in
                       mappy.fastx_read(str(fasta))})
    run = json.loads((BASE / f'run-{tag}-{component}.log').read_text())
    route = run['spine']['phasing']['rescore']['selected_route']
    missing = sorted({piece['source'] for locus in route for slot in locus
                      if isinstance(slot, dict) for piece in slot.get('segments', [])}
                     - source.keys())
    if missing:
        names = {int(fields[0]): fields[1] for line in NAMES.read_text().splitlines()
                 if (fields := line.split('\t'))[0].isdigit()}
        extra = BASE / f'sources-score-{component}.fa'
        subprocess.run(['/home/erikg/impg-genome-inference/target/release/examples/genotyping_spell_dump',
                        '/home/erikg/yeast/yeast235.agc', str(extra),
                        *(f'{index}:{names[index]}' for index in missing)],
                       check=True, stdout=subprocess.DEVNULL)
        source.update({int(name.split(':', 1)[0]): seq for name, seq, *_ in
                       mappy.fastx_read(str(extra))})
    assert all(piece['source'] in source for locus in route for slot in locus
               if isinstance(slot, dict) for piece in slot.get('segments', []))
    inferred = []
    for locus, slots in enumerate(route):
        for slot_index, slot in enumerate(slots):
            if isinstance(slot, dict):
                pieces = [rc(segment) if part.get('reverse') else segment
                          for part in slot['segments']
                          if (segment := source[part['source']][part['start']:part['end']])]
                if pieces:
                    inferred.append((locus, slot_index, ''.join(pieces)))
    aligners = [mappy.Aligner(seq=seq, preset='asm5', best_n=20) for seq in truths]
    assert all(aligners)
    # Truth self/permutation gates use the same complete-haplotype mapper and
    # edit parser; skipped query spans, target omissions or inversions fail.
    self_score = best_configuration([(0, 0, truths[0]), (0, 1, truths[1])], truths, aligners)
    swap_score = best_configuration([(0, 0, truths[1]), (0, 1, truths[0])], truths, aligners)
    assert self_score['error_bp'] == swap_score['error_bp'] == 0
    assert self_score['chosen_configuration'] == 0 and swap_score['chosen_configuration'] == 1
    assert self_score['unaligned_query_bp'] == swap_score['unaligned_query_bp'] == 0
    scored = best_configuration(inferred, truths, aligners)
    scored.update({'component': component, 'tag': tag, 'run_exit': (BASE / f'run-{tag}-{component}.exit').read_text().strip(),
                   'original_wall_seconds': int((BASE / f'run-{tag}-{component}.wall').read_text()),
                   'selected_ploidy': run['spine']['phasing']['rescore']['selected_ploidy'],
                   'inferred_slot_count_by_locus': [sum(isinstance(slot, dict) for slot in slots)
                                                    for slots in route],
                   'truth_self_error_bp': self_score['error_bp'],
                   'truth_swapped_error_bp': swap_score['error_bp']})
    assert scored['run_exit'] == '0'
    print(json.dumps(scored, indent=2))


if __name__ == '__main__':
    main(*sys.argv[1:3])
