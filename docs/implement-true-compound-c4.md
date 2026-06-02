# Implement true compound-region C4 crush

Date: 2026-06-02

Status: failed the C4 quality criterion after implementing and validating compound-region selection. The selected C4 region is now a true compound run, but the direct SweepGA/wfmash/seqwish induction for that extracted sequence set expands the local graph before lacing.

## Implementation

- Added `compound-run` iterative multi-level candidates in `src/resolution.rs`.
- Compound runs merge contiguous, overlapping, and nested POVU sites by root-path order with a small gap allowance and a span cap derived from `--window-target-bp`.
- Compound-run candidates are generated before top-level, parent-descendant, sibling, and sliding windows in `largest`, `parent`, and `combined` modes.
- Compound runs are routed to `ResolutionMethod::Sweepga` rather than the small direct POA/POASTA tier.
- Candidate logging now reports `source_max_span` so the resolver can prefer a compound run that contains the largest unresolved source flubble instead of a generic all-path window.
- Replacement lacing now reuses surviving original segment IDs when a replacement segment has an identical sequence to a segment that remains elsewhere in the graph. This prevents avoidable duplicate S-lines when the SYNG input has repeat-glued segments outside the replaced interval.

## Synthetic tests

The focused resolution suite now includes these targeted compound-region tests:

- `compound_run_merges_adjacent_bubbles_before_single_flubbles`
- `compound_run_condenses_correlated_adjacent_bubbles_more_than_separate_lacing`
- `compound_run_selects_largest_nested_region`
- `compound_replacement_lacing_preserves_exact_path_spellings`
- `compound_run_keeps_short_internal_anchor_inside_region`
- `replacement_lacing_reuses_surviving_original_segments_by_sequence`
- `compound_run_priority_prefers_largest_source_site_span`

These cover adjacent underaligned bubbles, nested/largest-region selection, exact path spelling preservation after whole-region replacement, repeated short anchors inside the region, and repeat-glued sequence reuse during lacing.

## Validation commands

```bash
CPATH=/home/erikg/wfmash/build/vendored_htslib/include \
CPLUS_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
C_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib \
LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:${LD_LIBRARY_PATH:-} \
CMAKE_PREFIX_PATH=/home/erikg/wfmash/build/vendored_htslib:/home/erikg/micromamba \
cargo test --lib resolution::
```

Result: 87 passed, 0 failed.

```bash
CPATH=/home/erikg/wfmash/build/vendored_htslib/include \
CPLUS_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
C_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib \
LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:${LD_LIBRARY_PATH:-} \
CMAKE_PREFIX_PATH=/home/erikg/wfmash/build/vendored_htslib:/home/erikg/micromamba \
cargo build --release --bin impg
```

Result: passed with existing dead-code warnings.

`cargo test --bin impg` was not required because CLI syntax was not changed.

## C4 command

Input:

```text
/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa
```

Final ranked run:

```bash
OUT=/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z
TMPDIR="$OUT/tmp" IMPG_CRUSH_DEBUG_DIR="$OUT/debug" \
LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:${LD_LIBRARY_PATH:-} \
target/release/impg crush \
  --gfa /home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa \
  --output "$OUT/compound_largest_wfmash.gfa" \
  --method iterative-multi-level \
  --window-mode largest \
  --window-target-bp 30k \
  --max-window-sites 8 \
  --candidate-limit 999 \
  --max-iterations 8 \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --max-transclosure-cells 0 \
  --max-poasta-cells 0 \
  --sweepga-aligner wfmash \
  --sweepga-no-filter true \
  --min-match-length off \
  --polish-method poasta \
  --polish-rounds until-done \
  -t 32 -v 1
```

Result: exit 0, 8 resolved, 0 bailed, 8 candidates seen across 8 rounds. Wall time was 1:15.41 and max RSS was 1,858,124 KB.

## Artifacts

- Final GFA: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/compound_largest_wfmash.gfa`
- Final graph-report TSV: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/compound_largest_wfmash.graph-report.povu.tsv`
- Final gfalook `-m` PNG: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/compound_largest_wfmash.mean-depth.png`
- Input graph-report TSV: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/input_syng_control.graph-report.povu.tsv`
- PGGB graph-report TSV: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/pggb_control.graph-report.povu.tsv`
- Crush stderr log: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/crush.stderr.log`
- Debug directory: `/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/debug`

## Metric comparison

All rows use `target/release/impg graph-report --format tsv --povu`.

| Graph | Segments | Links | Segment bp | bp-weighted cov | Singleton bp | Link jump p99 | Link jump max | Path jump p99 | Path jump max | WS p99 | WS max | WS bridges >= threshold | POVU sites | POVU leaf sites |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Input SYNG control | 18970 | 21863 | 425742 | 250.596544 | 81420 | 16082 | 18831 | 5580 | 18831 | 130878 | 422205 | 168197 | 2479 | 2217 |
| This compound-run output | 19080 | 22539 | 537822 | 198.373205 | 151503 | 16161 | 18943 | 7407 | 18943 | 180678 | 526163 | 198526 | 2477 | 2208 |
| PGGB control | 13288 | 16240 | 234524 | 454.919215 | 2890 | 5 | 3488 | 5 | 3488 | 16 | 73181 | 22851 | 2870 | 2833 |

The compound-run output is not materially close to PGGB. It is worse than the input on segment count, links, segment bp, bp-weighted coverage, singleton bp, path jump p99, path white-space p99/max, and long white-space bridge count.

Visual assessment of the `gfalook -m` PNG: the output still has a large sparse white region and dense long cross-links. It does not resemble the compact PGGB C4 rendering.

## Debugging answers

Why do residual C4 bubbles remain split into A/B-looking clusters instead of being aligned into one condensed local graph?

The resolver now selects the C4 compound run as one region, but the replacement inducer itself does not condense the extracted spellings. Round 1 selected:

```text
source=compound-run sites=106 source_max_span=32834 root_span=36515bp traversals=45 max=30148 median=3777 total=697376 unique_steps=1070
```

The replacement built from those exact extracted sequences had only 62 segments, but it grew the local sequence mass and singleton mass:

```text
input_segments=1070, output_segments=62
input_bp=28755, output_bp=64122
input_singleton_bp=577, output_singleton_bp=33998
```

So the A/B-looking residual structure is not primarily a lacing failure. The all-vs-all induction is aligning a subset of shared material while leaving large singleton-heavy haplotype blocks, and exact spelling-preserving lacing keeps those blocks.

Are the unresolved structures single POVU flubbles, adjacent flubbles that need to be merged, or artifacts of repeat-glue in the input SYNG graph?

They are adjacent/nested POVU flubbles around the C4 residual cluster. Admission-only parent-mode diagnostics showed C4 compound candidates with many merged sites, including the selected `sites=106` run, while individual top-level C4 sites remained visible. The input also has repeat-glued/reused sequence context, which is why the lacing change to reuse surviving original segment sequences is useful, but that is not enough to fix the C4 quality because the induced replacement graph is already poor before global lacing.

Does a direct compound-region SweepGA/seqwish induction on the same extracted sequences produce a condensed graph before lacing?

No, not with the currently working backend/settings. The first C4 compound replacement expanded from 28,755 bp to 64,122 bp and from 577 singleton bp to 33,998 singleton bp before it was laced into the global graph. That points at induction parameters/input preparation/backend behavior.

FastGA was also tested on a large compound selection and failed in `FAtoGDB` with a buffer overflow. Wfmash completes, but with the current wrapper settings it under-condenses the C4 compound sequence set. The next concrete fix should focus on reproducing PGGB-like all-vs-all induction for the extracted compound FASTA, not on additional suppression or candidate skipping.
