# Expand multi-bubble C4 SweepGA windows

Task: `expand-multi-bubble`

Output directory:

```text
/home/erikg/impg/data/c4_expand_multi_bubble_20260529T020548Z
```

Best PNG:

```text
https://hypervolu.me/~erik/impg/c4-expand-multi-bubble-best.png
```

Upload confirmation:

```text
-rw-r--r-- 1 erik erik 971K May 29 02:51 www/impg/c4-expand-multi-bubble-best.png
```

## Lineage

This task is a follow-up to `iterative-multi-level` and used
`c4-iterative-multi-level-final1.png` plus its metrics
`16,407 S / 391,391 bp / 101 trivial-stringy` as the visual and numeric
baseline. The implementation keeps the prior `method=iterative-multi-level`
primitive, but changes its multi-bubble proposal and selection behavior.

| Modification | Lineage | Rationale |
| --- | --- | --- |
| Default max window sites increased from 4 to 8. | `iterative-multi-level` sibling/sliding windows. | The remaining C4 under-alignment appears to span runs of adjacent flubbles, so the generator should be able to propose wider runs. |
| Default candidate cap increased from 96 to 192. | Prior bounded candidate queue. | Wider windows generate more sibling/sliding combinations; the cap needs room for these without making unbounded C4 runs the default. |
| Multi-site windows use SweepGA first. | Existing per-candidate `auto` method selection. | Initial multi-bubble runs should be aligned as a collection before trusting lower-level bubble boundaries. |
| Added explicit total-sequence cap for iterative windows. | Existing `max_total_sequence` candidate bound. | Aggressive C4 windows otherwise include multi-megabase traversal collections that can dominate wall time. |
| Replaced per-candidate whole-graph scoring in the hot path with local objective sorting plus non-overlapping batch selection. | Prior one-at-a-time whole-graph objective trial. | Keeps the explicit graph-size objective while avoiding the stall seen when every candidate required a full rewrite trial. |
| Added bounded whole-graph fallback. | Prior strict global objective gate. | If the cheap batch fails globally, try up to 24 locally promising non-overlapping alternatives so one bad wide window does not block a good later candidate. |
| Expanded logs with window sites, root-span bp, region bp, cheap-selection counters, fallback counters, objective deltas, graph metrics, and wall time. | Prior iterative logs. | Real C4 experiments need enough telemetry to distinguish search, build, and global-objective failure modes. |
| Added CLI aliases `multi-bubble-window-mode`, `wide-window-sites`, `multi-bubble-window-sites`, and `max-multi-bubble-window-sites`. | Existing `window-mode` and `max-window-sites`. | Makes the new multi-bubble interpretation addressable without changing the stable old option names. |

## Commands

Real C4 input:

```text
/home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.nosort.gfa
```

Path-spelling reference:

```text
/home/erikg/impg/data/c4_low_min_match_20260528T163541Z/ref/run.nosort.gfa
```

Wide capped variants used this command shape:

```bash
target/release/impg crush \
  --gfa /home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.nosort.gfa \
  --output "$out/<variant>/run.nosort.gfa" \
  --method iterative-multi-level \
  --window-mode sibling|sliding|combined \
  --window-target-bp 30000 \
  --max-window-sites 8 \
  --candidate-limit 20 \
  --max-total-sequence 2500000 \
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

The best bounded prior-shape/fallback probe used:

```bash
target/release/impg crush \
  --gfa /home/erikg/impg/data/c4_crush_iterative_multi_20260527T050056Z/run.nosort.gfa \
  --output "$out/sibling_prior_shape_fb_r3/run.nosort.gfa" \
  --method iterative-multi-level \
  --window-mode sibling \
  --window-target-bp 12000 \
  --max-window-sites 4 \
  --candidate-limit 16 \
  --max-total-sequence 2500000 \
  --min-objective-delta 1 \
  --max-iterations 3 \
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

Final bounded POASTA smoothing probe:

```bash
target/release/impg crush \
  --gfa "$out/sibling_prior_shape_fb_r3/run.nosort.gfa" \
  --output "$out/sibling_prior_shape_fb_r3_poasta_smooth/run.nosort.gfa" \
  --method poasta \
  --max-iterations 1 \
  --max-span 1000 \
  --max-traversal-len 1000 \
  --max-median-traversal-len 500 \
  --max-total-sequence 250000 \
  --max-traversals 512 \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

The best output was rendered and uploaded with:

```bash
gfasort -i "$out/sibling_prior_shape_fb_r3/run.nosort.gfa" \
  -o "$out/sibling_prior_shape_fb_r3/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/sibling_prior_shape_fb_r3/run.Ygs.gfa" \
  -o "$out/sibling_prior_shape_fb_r3/c4-expand-multi-bubble-best.png" \
  -m -x 2200 -y 1200
scp "$out/sibling_prior_shape_fb_r3/c4-expand-multi-bubble-best.png" \
  erik@hypervolu.me:www/impg/c4-expand-multi-bubble-best.png
ssh erik@hypervolu.me 'ls -lh www/impg/c4-expand-multi-bubble-best.png'
```

## Real C4 Results

All completed outputs preserve 465/465 paths against the no-crush C4 reference.
Missing paths, extra paths, and spelling mismatches were all zero. The
`trivial-stringy` values below are the project graph-quality counters logged by
`impg::resolution`, not a raw count of low-coverage segments.

| Graph | Segments | Links | Paths | Segment bp | Trivial-stringy | Matched | Missing | Extra | Mismatches |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Iterative-multi input | 16,547 | 23,610 | 465 | 396,438 | 103 | 465 | 0 | 0 | 0 |
| `iterative-multi-level-final1` | 16,407 | 23,200 | 465 | 391,391 | 101 | 465 | 0 | 0 | 0 |
| Neighbor-merge baseline | 18,761 | 25,763 | 465 | 461,241 | 72 | 465 | 0 | 0 | 0 |
| `sibling_wide_cap_fb_r1` | 16,547 | 23,610 | 465 | 396,438 | 103 | 465 | 0 | 0 | 0 |
| `sliding_wide_cap_fb_r1` | 16,547 | 23,610 | 465 | 396,438 | 103 | 465 | 0 | 0 | 0 |
| `combined_wide_cap_fb_r1` | 16,547 | 23,610 | 465 | 396,438 | 103 | 465 | 0 | 0 | 0 |
| `sibling_prior_shape_fb_r1` | 16,227 | 22,882 | 465 | 392,671 | 100 | 465 | 0 | 0 | 0 |
| `combined_prior_shape_fb_r1` | 16,227 | 22,882 | 465 | 392,671 | 100 | 465 | 0 | 0 | 0 |
| `sibling_prior_shape_fb_r3` | 16,227 | 22,882 | 465 | 392,671 | 100 | 465 | 0 | 0 | 0 |
| `sibling_prior_shape_fb_r3_poasta_smooth` | 17,692 | 25,041 | 465 | 450,460 | 268 | 465 | 0 | 0 | 0 |

Main run telemetry:

| Variant | Window mode | Target bp | Max sites | Iterations requested | Generated | Evaluated | Accepted | Fallback trials | Global delta | Region bp distribution | Before | After | Wall |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- | ---: |
| `sibling_wide_r1` | sibling | 30,000 | 8 | 1 | 29,148 | 32 | stopped | n/a | n/a | n=32, min=1,446,010, p50=3,497,896, p90=3,533,099, max=3,587,846 | 16,547 S / 396,438 bp / 103 stringy | stopped by agent after the wide first batch proved too expensive | about 7m |
| `sibling_wide_cap_fb_r1` | sibling | 30,000 | 8 | 1 | 29,009 | 20 | 0 | 20 | -6,553,381 | n=20, min=1,087,428, p50=1,989,792, p90=2,110,529, max=2,155,412 | 16,547 S / 396,438 bp / 103 stringy | unchanged | 3:00.37 |
| `sliding_wide_cap_fb_r1` | sliding | 30,000 | 8 | 1 | 19,678 | 20 | 0 | 20 | -5,072,016 | n=20, min=728,567, p50=1,681,040, p90=2,384,540, max=2,388,821 | 16,547 S / 396,438 bp / 103 stringy | unchanged | 2:12.80 |
| `combined_wide_cap_fb_r1` | combined | 30,000 | 8 | 1 | 33,137 | 20 | 0 | 20 | -6,553,381 | n=20, min=1,087,428, p50=1,989,792, p90=2,110,529, max=2,155,412 | 16,547 S / 396,438 bp / 103 stringy | unchanged | 3:33.16 |
| `sibling_prior_shape_fb_r1` | sibling | 12,000 | 4 | 1 | 14,222 | 16 | 2 | 3 | 4,556,217 | n=16, min=945,324, p50=2,008,953, p90=3,481,625, max=4,820,734 | 16,547 S / 396,438 bp / 103 stringy | 16,227 S / 392,671 bp / 100 stringy | 1:59.01 |
| `combined_prior_shape_fb_r1` | combined | 12,000 | 4 | 1 | 15,798 | 16 | 2 | 3 | 4,556,217 | n=16, min=945,324, p50=2,008,953, p90=3,481,625, max=4,820,734 | 16,547 S / 396,438 bp / 103 stringy | 16,227 S / 392,671 bp / 100 stringy | 2:08.95 |
| `sibling_prior_shape_fb_r3` round 1 | sibling | 12,000 | 4 | 3 | 14,222 | 16 | 2 | 3 | 4,556,217 | n=16, min=945,324, p50=2,008,953, p90=3,481,625, max=4,820,734 | 16,547 S / 396,438 bp / 103 stringy | 16,227 S / 392,671 bp / 100 stringy | 3:27.26 total |
| `sibling_prior_shape_fb_r3` round 2 | sibling | 12,000 | 4 | 3 | 14,625 | 16 | 0 | 13 | -30,975,533 | n=16, min=217,621, p50=2,891,207, p90=4,791,461, max=4,820,734 | 16,227 S / 392,671 bp / 100 stringy | unchanged | included above |
| `sibling_prior_shape_fb_r3_poasta_smooth` | final smoothing | n/a | n/a | 1 | 3,699 | 1,408 | 1,408 | n/a | n/a | total max=2,186,334 | 16,227 S / 392,671 bp / 100 stringy | 17,692 S / 450,460 bp / 268 stringy | 1:29.49 |

The accepted prior-shape/fallback batch had local objective summary:

```text
n=2, source-sites n=2, min=4, p50=4, p90=4, max=4, total=8,
score-delta sum=8,036,505, segment-delta sum=1,710,
segment-bp-delta sum=7,144, input_segments=2,599,
output_segments=889, input_bp=18,110, output_bp=10,966
```

## Interpretation

The wider 30 kb / 8-site multi-bubble SweepGA windows did not solve the
remaining C4 under-alignment. They generated many large regions and produced
locally compressive replacements, but the explicit whole-graph objective
rejected every completed wide variant. The visible best PNG still has broad
under-aligned vertical bands and large gap columns, so the fundamental clustering
problem remains.

The useful improvement came from the bounded global fallback, not from wider
windows. It recovers two non-overlapping sibling-run SweepGA candidates from the
prior 12 kb / 4-site shape and improves the previous best segment count and
trivial-stringy count:

| Comparison | Segments | Segment bp | Trivial-stringy | Notes |
| --- | ---: | ---: | ---: | --- |
| `iterative-multi-level-final1` | 16,407 | 391,391 | 101 | Prior best visual baseline. |
| `expand-multi-bubble` best | 16,227 | 392,671 | 100 | 180 fewer segments and 1 fewer stringy counter than final1, but 1,280 more segment bp. |
| Neighbor-merge | 18,761 | 461,241 | 72 | Still better on stringy count, much worse on size. |

Final POASTA smoothing after the wide SweepGA rounds is not a good default for
this C4 fragment. It preserved paths, but grew the graph to
`17,692 S / 450,460 bp / 268 stringy`.

Suggested default strategy:

- Keep the explicit graph-size gate enabled.
- Keep multi-site windows SweepGA-first.
- Keep the cheap local objective plus non-overlapping batch selection, with the
  bounded whole-graph fallback as the safety valve.
- Do not make the 30 kb / 8-site C4-wide window shape the default. It is useful
  as an experiment mode, but the completed real C4 variants were all globally
  rejected.
- Use the 12 kb / 4-site sibling or combined shape as the current aggressive C4
  post-pass if a slight segment-count improvement is worth the small segment-bp
  regression relative to `iterative-multi-level-final1`.

## Validation

- `cargo build --release` passed.
- Focused iterative multi-level unit tests passed with
  `cargo test iterative_multi_level --lib -- --nocapture`.
- CLI alias parser tests passed with
  `cargo test test_gfa_output_format_accepts_iterative_multi_level_crush_params --bin impg -- --nocapture`
  and
  `cargo test test_gfa_output_format_accepts_wide_multi_bubble_window_aliases --bin impg -- --nocapture`.
- C4 fragment regressions passed with
  `cargo test --test test_crush_integration c4_fragment -- --test-threads=1 --nocapture`
  and
  `cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --test-threads=1 --nocapture`.
- Real C4 variants completed for wider sibling windows, wider sliding windows,
  combined sibling plus sliding windows, repeated best variant, and final
  bounded POASTA smoothing.
- All completed real C4 outputs preserve 465/465 paths with zero spelling
  mismatches.
- Best PNG was generated, uploaded to hypervolu.me, and confirmed with `ssh ls`.
- `cargo test --all -- --test-threads=1` and `cargo install --path .` were run
  after the final patch before commit.
