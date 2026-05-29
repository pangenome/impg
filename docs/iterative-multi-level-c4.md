# Iterative multi-level C4 crush-to-size

Task: `iterative-multi-level`

Output directory:

```text
/home/erikg/impg/data/c4_iterative_multi_level_20260528T235244Z
```

Best PNG:

```text
https://hypervolu.me/~erik/impg/c4-iterative-multi-level-final1.png
```

Upload confirmation:

```text
-rw-r--r-- 1 erik erik 965K May 29 00:39 www/impg/c4-iterative-multi-level-final1.png
```

## Implementation

`method=iterative-multi-level` adds an explicit window-crush search mode to
`impg crush`. Each round re-runs POVU/flubble decomposition on the current graph,
generates contained regions from multiple views, builds replacements through the
existing transparent condenser path, validates path spellings, and accepts only
candidate sets that improve the configured graph-size objective.

Candidate sources implemented:

| Source | Purpose |
| --- | --- |
| `top-level` | Existing top-level flubble regions. |
| `sibling-run` | Same-parent/same-level neighboring site runs. |
| `parent-descendants` | Parent regions plus contained descendants up to target size. |
| `sliding-window` | Reference/path-order neighboring flubble windows. |
| `level-window` | Multi-level windows across adjacent flubble levels. |
| `stringy-neighborhood` | Local neighborhoods around trivial-stringy/under-aligned signatures. |

CLI controls:

```text
--method iterative-multi-level
--window-mode sibling|sliding|combined
--window-target-bp N
--max-window-sites N
--candidate-limit N
--min-objective-delta N
```

The objective used for acceptance is:

```text
segment_bp * 1024 + segments * 16 + path_steps
```

The resolver logs generated candidates, candidates evaluated after the cap,
accepted candidates, source counts, local objective deltas, global objective
deltas, segments, segment_bp, trivial_stringy, build/rewrite wall time, and
path-preserving rewrite validation.

Two real-C4 hardening fixes landed during this task:

- Multi-level generated windows now materialize traversal sequences before
  invoking POASTA/SPOA/SweepGA. Without this, real C4 windows reached aligners as
  empty traversal records.
- Candidate condenser panics are caught and counted as failed candidates for this
  aggressive search mode. Accepted replacements still must pass exact path
  validation.

## Commands

Input for the final bounded experiment was the prior accepted iterative-multi C4
graph used by the neighbor-merge baseline:

```text
/home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.nosort.gfa
```

Path spellings were validated against the latest no-crush C4 reference:

```text
/home/erikg/impg/data/c4_low_min_match_20260528T163541Z/ref/run.nosort.gfa
```

Final bounded variant command shape:

```bash
target/release/impg crush \
  --gfa /home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.nosort.gfa \
  --output "$out/<variant>/run.nosort.gfa" \
  --method iterative-multi-level \
  --window-mode sibling \
  --window-target-bp 12000 \
  --max-window-sites 4 \
  --candidate-limit 16 \
  --min-objective-delta 1 \
  --max-iterations 1 \
  --auto-spoa-max-traversal-len 0 \
  --auto-poasta-max-traversal-len 0 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

Only `--window-mode` changed across `sibling`, `sliding`, and `combined`.
The one-round cap was chosen because per-candidate whole-graph objective trials
are expensive on real C4. An exploratory max-3 sibling run completed round 1 and
was in round 2 when stopped; the final artifact records the bounded cap-1 runs
that completed cleanly.

The best output was rendered with:

```bash
gfasort -i "$out/sibling_final1_iterative/run.nosort.gfa" \
  -o "$out/sibling_final1_iterative/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/sibling_final1_iterative/run.Ygs.gfa" \
  -o "$out/sibling_final1_iterative/c4-iterative-multi-level-final1.png" \
  -m -x 2200 -y 1200
scp "$out/sibling_final1_iterative/c4-iterative-multi-level-final1.png" \
  erik@hypervolu.me:www/impg/c4-iterative-multi-level-final1.png
ssh erik@hypervolu.me 'ls -lh www/impg/c4-iterative-multi-level-final1.png'
```

## Final C4 Results

All final runs preserved 465/465 path spellings against the no-crush reference
when keyed by the PanSN prefix before the coordinate suffix. Missing paths,
extra paths, and spelling mismatches were all zero.

| Variant | Mode | POVU sites | Polymorphic sites | Generated | Evaluated | Accepted | Global delta | Segments | Segment bp | Trivial-stringy | Wall |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `sibling_final1` | sibling | 4,590 | 4,320 | 14,222 | 16 | 1 | 5,462,331 | 16,407 | 391,391 | 101 | 2:06.93 |
| `sliding_final1` | sliding | 4,590 | 4,320 | 9,945 | 16 | 1 | 2,251,580 | 16,703 | 394,476 | 112 | 1:23.01 |
| `combined_final1` | combined | 4,590 | 4,320 | 15,798 | 16 | 1 | 5,462,331 | 16,407 | 391,391 | 101 | 2:17.31 |

Round-by-round accepted summaries:

| Variant | Round | Accepted source | Local segment delta | Local segment-bp delta | Before | After |
| --- | ---: | --- | ---: | ---: | --- | --- |
| `sibling_final1` | 1 | `sibling-run=1` | 923 | 6,890 | 16,547 S / 396,438 bp / 103 stringy | 16,407 S / 391,391 bp / 101 stringy |
| `sliding_final1` | 1 | `sliding-window=1` | 810 | 6,852 | 16,547 S / 396,438 bp / 103 stringy | 16,703 S / 394,476 bp / 112 stringy |
| `combined_final1` | 1 | `sibling-run=1` | 923 | 6,890 | 16,547 S / 396,438 bp / 103 stringy | 16,407 S / 391,391 bp / 101 stringy |

The best strict output is `sibling_final1`/`combined_final1`. It improves the
iterative-multi input by 140 segments, 5,047 segment-bp, and 2 trivial-stringy
neighborhoods while preserving every path spelling.

## Relaxed Probe

Before the global objective gate was tightened, a relaxed sliding run forced
acceptance with `--min-objective-delta=-1000000`:

| Output | Accepted | Segments | Segment bp | Trivial-stringy | Path mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| Raw relaxed sliding | 4 | 16,853 | 397,639 | 119 | 0 |
| Relaxed sliding + `gfaffix` | 4 | 16,849 | 397,607 | not rerun | 0 |

This confirmed the objective gate is necessary: locally compact windows can
still grow the whole graph, and gfaffix did not turn that relaxed rewrite into a
size win.

## Comparisons

| Graph | Segments | Links | Paths | Segment bp | Trivial-stringy | Notes |
| --- | ---: | ---: | ---: | ---: | ---: | --- |
| No-crush reference | 18,048 | 20,943 | 465 | 389,354 | 506 | Latest reference in `c4_low_min_match.../ref`. |
| Iterative-multi baseline | 16,547 | 23,610 | 465 | 396,438 | 103 | Prior accepted input for this search. |
| Current best neighbor-merge | 18,761 | 25,763 | 465 | 461,241 | 72 | Better stringy, worse size. |
| Low-min + POASTA PNG | 21,540 before smoothing | n/a | 465 | 797,812 before smoothing | 491 before smoothing | Upstream low-min default/off top-flubble run; smoothing preserved paths. |
| `iterative-multi-level` best | 16,407 | 23,200 | 465 | 391,391 | 101 | This task, sibling/combined final1. |
| PGGB control | 13,288 | 16,240 | 465 | 234,524 | 12 | Prior comparison target. |

## Interpretation

Multi-level/window iteration does find a real path-preserving size win that the
neighbor-merge baseline missed, but it does not fix C4 under-alignment or the
large remaining PGGB gap. The best one-round result is only a modest
condensation: segment count drops below the iterative-multi baseline, and
segment-bp moves close to the no-crush reference, but trivial-stringy count only
moves from 103 to 101.

The central finding is that region discovery is not the only bottleneck. Many
multi-level windows look strongly compressive inside their local replacement
graph, but whole-graph rewrite validation shows they can still increase global
segments or segment-bp. A default strategy should therefore keep the explicit
global objective gate enabled and use `combined` or `sibling` windows with a
small candidate cap as an optional aggressive pass after the existing
iterative-multi stage. It should not replace the current default pipeline until
post-rewrite compaction or a cheaper whole-graph objective search closes more of
the PGGB gap.

## Validation

- `cargo test --all -- --test-threads=1` passed.
- `cargo install --path .` completed and replaced the local `impg` and
  `gfaffix` binaries from this worktree.
- `cargo test iterative_multi_level --lib -- --nocapture` passed.
- `cargo test test_gfa_output_format_accepts_iterative_multi_level_crush_params --bin impg -- --nocapture` passed.
- C4 fragment regressions passed earlier in this task:
  `cargo test --test test_crush_integration c4_fragment -- --test-threads=1 --nocapture`
  and
  `cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --test-threads=1 --nocapture`.
- Real C4 final variants completed for sibling, sliding, and combined
  one-round bounded search.
- Final variants preserve 465/465 paths and have zero spelling mismatches
  against the no-crush reference.
- Best PNG was generated, uploaded to hypervolu.me, and confirmed with `ssh ls`.
