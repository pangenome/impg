#!/usr/bin/env python3
"""Generate the balanced 15x/15x diploid sample with the original read model.

Pure S288C and SK1 routes, 150 bp, uniform start, 50% reverse complement,
error-free all-I FASTQ and deterministic seed. Reuses the existing yeast panel;
no new syng or truth-derived selection feature is constructed.
"""
import gzip
import hashlib
import json
import math
import random
import subprocess
from pathlib import Path

PANEL = Path('/home/erikg/yeast/syng-k63-s8-seed7-acgt-only-pos64/yeast235.syng')
ARCHIVE = Path('/home/erikg/yeast/yeast235.agc')
EXTRACTOR = Path('/home/erikg/impg/target/experiments/yeast-partition-renderings/extract-bed')
OUT = Path('/home/erikg/yeast/genome-balanced-diploid-validation-20260930')
READ_LENGTH = 150
SEED = 20260929
HOMOLOG_DEPTH = 15.0


def main():
    (OUT / 'private-truth').mkdir(parents=True, exist_ok=True)
    metadata = {}
    for line in Path(str(PANEL) + '.names').read_text().splitlines():
        fields = line.split('\t')
        metadata[fields[1]] = int(fields[2])
    first = sorted(name for name in metadata if name.startswith('S288C#0#chr'))
    second = sorted(name for name in metadata if name.startswith('SK1#0#chr'))
    assert len(first) == len(second) == 17
    assert [name.split('#')[-1] for name in first] == [name.split('#')[-1] for name in second]
    names = sorted(first + second)
    bed = OUT / 'private-truth/truth-source-paths.bed'
    bed.write_text(''.join(f'{name}\t0\t{metadata[name]}\n' for name in names))
    fasta = OUT / 'private-truth/truth-source-paths.fa'
    with (OUT / 'extract.log').open('wb') as log:
        subprocess.run([str(EXTRACTOR), str(ARCHIVE), str(bed), str(fasta)],
                       stdout=log, stderr=subprocess.STDOUT, check=True)
    sequences = {}
    current = None
    for line in fasta.read_bytes().splitlines():
        if line.startswith(b'>'):
            current = line[1:].decode().rsplit(':', 1)[0]
            sequences[current] = bytearray()
        else:
            sequences[current].extend(line)
    assert set(sequences) == set(names)

    rng = random.Random(SEED)
    complement = bytes.maketrans(b'ACGTacgt', b'TGCAtgca')
    per_path = []
    total = 0
    reads = OUT / 'reads.fastq.gz'
    with reads.open('xb') as raw:
        with gzip.GzipFile(filename='', mode='wb', fileobj=raw, mtime=0,
                           compresslevel=6) as fq:
            for name in names:
                sequence = bytes(sequences[name])
                length = len(sequence)
                count = math.ceil(HOMOLOG_DEPTH * length / READ_LENGTH)
                for _ in range(count):
                    start = rng.randrange(length - READ_LENGTH + 1)
                    read = sequence[start:start + READ_LENGTH]
                    if rng.getrandbits(1):
                        read = read.translate(complement)[::-1]
                    fq.write(b'@sim-' + str(total).encode() + b'\n' + read
                             + b'\n+\n' + b'I' * READ_LENGTH + b'\n')
                    total += 1
                per_path.append({'path': name, 'slot': 1 if name in first else 2,
                                 'dose': HOMOLOG_DEPTH, 'length': length,
                                 'reads': count, 'nominal_depth': count * READ_LENGTH / length,
                                 'sequence_sha256': hashlib.sha256(sequence).hexdigest()})
    truth = {'read_length': READ_LENGTH, 'seed': SEED,
             'slot1': [{'component': name.split('#')[-1], 'path': name,
                        'dose': HOMOLOG_DEPTH} for name in first],
             'slot2': [{'component': name.split('#')[-1], 'path': name,
                        'dose': HOMOLOG_DEPTH} for name in second],
             'per_path': per_path, 'total_reads': total}
    (OUT / 'private-truth/truth-diploid.json').write_text(json.dumps(truth, indent=1))
    print(f'wrote {total} reads across 34 pure routes '
          f'({HOMOLOG_DEPTH}x S288C + {HOMOLOG_DEPTH}x SK1)')


if __name__ == '__main__':
    main()
