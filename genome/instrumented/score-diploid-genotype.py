#!/usr/bin/env python3
"""Assessment-only, injective per-locus matching of two called and truth alleles.

Score only loci with a full, unique SK1 ortholog of the S288C axis window;
missing correspondence is bracketed, not silently treated as an empty allele.
Usage: score-diploid-genotype.py chrMT smoke2
"""
import json
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


def rc(sequence):
    return sequence.translate(TRANSLATE)[::-1]


def local_distance(call, truth):
    if call is None:
        return len(truth)
    distance = edlib.align(call, truth, mode='NW', task='distance')['editDistance']
    if distance < 0:
        raise ValueError('local global edit alignment failed')
    return distance


def match_genotype(calls, truths):
    assert len(calls) == len(truths) == 2 and all(truths)
    costs = [[local_distance(call, truth) for truth in truths] for call in calls]
    direct = costs[0][0] + costs[1][1]
    swapped = costs[0][1] + costs[1][0]
    # Equal totals do not identify ancestry; the deterministic direct
    # representative is retained without inventing an informative switch.
    assignment = 'direct' if direct <= swapped else 'swapped'
    indices = [0, 1] if assignment == 'direct' else [1, 0]
    error = min(direct, swapped)
    return {'error_bp': error, 'truth_bp': sum(map(len, truths)),
            'direct_error_bp': direct, 'swapped_error_bp': swapped,
            'assignment': assignment, 'assignment_identifiable': direct != swapped,
            'called_slot_to_truth': [TRUTH_NAMES[index] for index in indices],
            'called_slot_error_bp': [costs[slot][index] for slot, index in enumerate(indices)],
            'truth_allele_bp': dict(zip(TRUTH_NAMES, map(len, truths)))}


def full_unique_homolog(window, mapper):
    hits = [hit for hit in mapper.map(window)
            if hit.q_st == 0 and hit.q_en == len(window)
            and hit.is_primary and hit.mapq > 0]
    return hits[0] if len(hits) == 1 else None


def spell(slot, source):
    if not isinstance(slot, dict):
        return None
    pieces = []
    for segment in slot.get('segments', []):
        sequence = source[segment['source']][segment['start']:segment['end']]
        pieces.append(rc(sequence) if segment.get('reverse') else sequence)
    return ''.join(pieces) or None


def load_source(route, truths, component):
    source = {}
    for fasta in (BASE / 'sources-all.fa', CORE):
        source.update({int(name.split(':', 1)[0]): sequence
                       for name, sequence, *_ in mappy.fastx_read(str(fasta))})
    names = {int(fields[0]): fields[1] for line in NAMES.read_text().splitlines()
             if (fields := line.split('\t'))[0].isdigit()}
    truth_fasta = {f'{name}#0#{component}': sequence
                   for name, sequence in zip(TRUTH_NAMES, truths)}
    for index, name in names.items():
        if name in truth_fasta:
            source[index] = truth_fasta[name]
    missing = sorted({part['source'] for slots in route for slot in slots
                      if isinstance(slot, dict) for part in slot.get('segments', [])}
                     - source.keys())
    if missing:
        extra = BASE / f'sources-score-{component}.fa'
        subprocess.run(['/home/erikg/impg-genome-inference/target/release/examples/genotyping_spell_dump',
                        '/home/erikg/yeast/yeast235.agc', str(extra),
                        *(f'{index}:{names[index]}' for index in missing)],
                       check=True, stdout=subprocess.DEVNULL)
        source.update({int(name.split(':', 1)[0]): sequence
                       for name, sequence, *_ in mappy.fastx_read(str(extra))})
    assert not (set(missing) - source.keys())
    return source


def score_loci(loci, route, truths, source, mapper, legacy_exact=None):
    assert len(loci) == len(route)
    rows = []
    attested_error = attested_truth_bp = 0
    bracketed_axis_bp = transitions = informative_boundaries = 0
    previous = None
    for locus, (meta, slots) in enumerate(zip(loci, route)):
        assert len(slots) == 2 and meta['locus'] == locus
        start, end = meta['axis_interval']
        first = truths[0][start:end]
        assert first and end > start
        homolog = full_unique_homolog(first, mapper)
        if legacy_exact is not None:
            assert (homolog is not None) == legacy_exact[locus], locus
        calls = [spell(slot, source) for slot in slots]
        row = {'locus': locus, 'axis_interval': [start, end],
               'called_allele_bp': [len(call) if call is not None else 0 for call in calls],
               'called_copy_count': sum(call is not None for call in calls),
               'called_dosage_classes': [entry['dosage'] for entry in meta['diplotype']['dosage']]}
        if homolog is None:
            row.update({'assessable': False, 'bracket_reason': 'no_full_unique_SK1_ortholog',
                        'assignment': None, 'truth_allele_bp': {'S288C': len(first), 'SK1': None}})
            bracketed_axis_bp += len(first)
            previous = None
        else:
            second = truths[1][homolog.r_st:homolog.r_en]
            if homolog.strand < 0:
                second = rc(second)
            assert second
            match = match_genotype(calls, (first, second))
            assert match_genotype((first, second), (first, second))['error_bp'] == 0
            assert match_genotype((second, first), (first, second))['error_bp'] == 0
            match.update({'assessable': True, 'SK1_truth_interval': [homolog.r_st, homolog.r_en],
                          'SK1_truth_strand': homolog.strand})
            row.update(match)
            attested_error += match['error_bp']
            attested_truth_bp += match['truth_bp']
            if match['assignment_identifiable']:
                if previous is not None and previous['locus'] + 1 == locus:
                    informative_boundaries += 1
                    if previous['assignment'] != match['assignment']:
                        transitions += 1
                previous = {'locus': locus, 'assignment': match['assignment']}
            else:
                previous = None
        rows.append(row)
    return {'error_bp': attested_error, 'truth_bp': attested_truth_bp,
            'distance': attested_error / attested_truth_bp,
            'assessable_loci': sum(row['assessable'] for row in rows),
            'bracketed_loci': sum(not row['assessable'] for row in rows),
            'bracketed_S288C_axis_bp': bracketed_axis_bp,
            'assignment_transitions': transitions,
            'informative_adjacent_boundaries': informative_boundaries,
            'assignment_ties': sum(row['assessable'] and not row['assignment_identifiable']
                                   for row in rows),
            'assignment_track': [row['assignment'] if row['assessable'] else None for row in rows],
            'rows': rows}


def main(component, tag):
    truth_fasta = {name.rsplit(':', 1)[0]: sequence for name, sequence, *_ in
                   mappy.fastx_read(str(BASE / 'private-truth/truth-source-paths.fa'))}
    truths = tuple(truth_fasta[f'{name}#0#{component}'] for name in TRUTH_NAMES)
    run = json.loads((BASE / f'run-{tag}-{component}.log').read_text())
    spine = run['spine']
    route = spine['phasing']['rescore']['selected_route']
    loci = spine['genotype_calls']['loci']
    source = load_source(route, truths, component)
    mapper = mappy.Aligner(seq=truths[1], preset='asm5')
    assert mapper
    old = json.loads((BASE / f'score-{tag}-{component}.json').read_text())
    assert len(old['rows']) == len(loci)
    scored = score_loci(loci, route, truths, source, mapper,
                        [row['exact_homolog'] for row in old['rows']])
    assert scored['assessable_loci'] == old['exact_homolog_windows']
    assert scored['bracketed_loci'] == old['bracketed_homolog_windows']
    assert spine['phasing']['rescore']['references'][0]['label'] == 'truth'
    scored.update({'component': component, 'tag': tag,
                   'product_run_exit': (BASE / f'run-{tag}-{component}.exit').read_text().strip(),
                   'product_wall_seconds': int((BASE / f'run-{tag}-{component}.wall').read_text()),
                   'selected_ploidy': spine['phasing']['rescore']['selected_ploidy'],
                   'truth_dosage': old['truth_dosage'],
                   'selected_single_class_dosage': old['selected_single_class_dosage'],
                   'yardstick': 'injective_per_locus_genotype_on_attested_SK1_homologs',
                   'bracket_policy': 'unattested SK1 orthology has no assigned local allele; not scored'})
    assert scored['product_run_exit'] == '0'
    print(json.dumps(scored, indent=2))


if __name__ == '__main__':
    main(*sys.argv[1:3])
