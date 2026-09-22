# Full C4 Whole-Region SweepGA Seed Plus 1kb POA Crush

Task: `c4-full-sweepga-seed-poa1kb`
Date: 2026-06-05

## Summary

Ran the requested POA-only small-bubble crush on the full C4 whole-region
SweepGA/seqwish seed graph:

`/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa`

Output directory:

`/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z`

The initial `--max-iterations 5` run materially converged before the limit:
rounds 1-3 selected and resolved candidates, and round 4 found no eligible
selected candidates. No continuation pass was needed.

## Artifacts

- Polished GFA:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa`
- Ygs-sorted polished GFA:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.gfa`
- Seed graph-report TSV:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/seed.graph-report.tsv`
- Polished graph-report TSV:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/poa1kb.graph-report.tsv`
- Crush timing/logs:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/crush.stderr.log`
  and `crush.stdout.log`
- Sort timing/logs:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/gfasort.stderr.log`
  and `gfasort.stdout.log`
- Render timing/logs:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/gfalook.stderr.log`
  and `gfalook.stdout.log`
- Path validation logs:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/compare_gfa_paths.stdout.log`
  and `compare_gfa_paths.stderr.log`
- Uploaded PNG:
  `http://hypervolu.me/~erik/impg/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.mean-depth.png`

The public PNG URL returned HTTP 200 on 2026-06-05 10:51 UTC.

## Software

- `impg`: `/home/erikg/impg/target/release/impg`, version `0.4.1`
- `gfasort`: `/home/erikg/.cargo/bin/gfasort`
- `gfalook`: `/home/erikg/.cargo/bin/gfalook`
- `compare_gfa_paths`:
  `/home/erikg/impg/target/release/examples/compare_gfa_paths`

`gfasort` and `gfalook` do not expose a `--version` flag; exact command lines
and logs are captured below and in the output directory.

## Exact Commands

```bash
/usr/bin/time -v /home/erikg/impg/target/release/impg crush -g /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa -o /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa --method poa --max-iterations 5 --max-traversal-len 1k --max-median-traversal-len 1k --max-total-sequence 1m --max-traversals 10k -t 32 -v 1
```

```bash
/usr/bin/time -v /home/erikg/impg/target/release/impg graph-report -g /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa -o /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/seed.graph-report.tsv --format tsv --povu -t 32 -v 1
```

```bash
/usr/bin/time -v /home/erikg/impg/target/release/impg graph-report -g /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa -o /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/poa1kb.graph-report.tsv --format tsv --povu -t 32 -v 1
```

```bash
/usr/bin/time -v /home/erikg/.cargo/bin/gfasort -i /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa -o /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.gfa -p Ygs -t 32 -v 1
```

```bash
/usr/bin/time -v /home/erikg/.cargo/bin/gfalook -i /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.gfa -o /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.mean-depth.png -x 2400 -y 1200 -m -t 32 -v 1
```

```bash
/usr/bin/time -v /home/erikg/impg/target/release/examples/compare_gfa_paths /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa
```

```bash
/usr/bin/time -v scp -o BatchMode=yes -o ConnectTimeout=20 /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.mean-depth.png erik@hypervolu.me:www/impg/
```

## Timing

| Step | Wall time | Max RSS |
| --- | ---: | ---: |
| POA crush | 0:34.64 | 1,898,296 KB |
| seed graph-report | 0:12.23 | 1,935,784 KB |
| polished graph-report | 0:11.52 | 1,986,336 KB |
| gfasort `Ygs` | 0:18.46 | 265,372 KB |
| gfalook `-m` render | 0:02.17 | 208,896 KB |
| path comparison | 0:00.58 | 356,352 KB |
| PNG upload | 0:01.44 | 8,192 KB |

## Crush Result

The crush log reports:

- Initial POVU on input: 3,080 sites, 2,575 polymorphic candidates, 1,742 roots.
- Round 1: 748 selected, 748 resolved.
- Round 2: 597 selected, 597 resolved.
- Round 3: 151 selected, 151 resolved.
- Round 4: no eligible candidates from 3,161 POVU sites.
- Final summary: 1,496 resolved, 0 bailed, 1,496 candidates seen across 3
  resolving rounds.

Compression quality thresholds were disabled as requested; no replacement was
rejected by ratio and no alternate aligner was attempted.

## Path Validation

`compare_gfa_paths` passed:

```text
expected_paths	465
observed_paths	465
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

## Graph-Report Delta

| Metric | Seed | POA 1kb |
| --- | ---: | ---: |
| segments | 8,186 | 9,079 |
| links | 11,303 | 12,339 |
| paths | 465 | 465 |
| path steps | 4,632,725 | 4,843,152 |
| total segment bp | 232,093 | 240,589 |
| path jump p99 | 6,687 | 3,682 |
| path white-space bp p99 | 189,754 | 100,506 |
| path white-space bp max | 217,619 | 240,054 |
| duplicate sequence fraction | 0.704618 | 0.750743 |
| local repeat context nodes | 572 | 657 |
| local repeat context occurrences | 997 | 1,106 |
| POVU sites | 3,080 | 3,161 |
| POVU leaf sites | 3,027 | 2,990 |

## Interpretation

The 1kb POA pass did exactly the local polishing it was allowed to do: it
exhausted the eligible direct small-bubble frontier, resolved 1,496 candidates,
and preserved every path spelling. The p99 path jump and p99 path white-space
metrics improved substantially, which is consistent with cleaning many local
small bubbles.

This was not a global C4 clean-up. The output remains in `REVIEW` status,
POVU sites increased slightly, duplicate sequence fraction and local repeat
contexts increased, and the maximum path white-space bridge worsened. Visually,
the uploaded whole-region PNG should be read as a polished local-bubble graph
that still has unresolved long-range tangles/tails and high-depth repeat-driven
structure. That remaining structure is outside the direct 1kb POA recipe and
would need a larger-scale method to address.

No source code was modified for this experiment.
