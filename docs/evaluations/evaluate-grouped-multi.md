# Grouped Multi-Bubble C4 Seed Evaluation

Task: `evaluate-grouped-multi`
Date: 2026-06-04

## Inputs And Artifacts

Input graph:

```text
/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_scaffold0.initial.gfa
```

Experiment output directory:

```text
/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z
```

Rendered Ygs PNGs were generated for the baseline and all eight crush outputs
under:

```text
/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/renders
```

Representative uploaded PNGs:

| View | URL |
| --- | --- |
| Baseline seed | https://hypervolu.me/~erik/impg/c4-grouped-multi-seed-baseline.png |
| Compact grouped top-1 | https://hypervolu.me/~erik/impg/c4-grouped-multi-seed-compact-top1.png |
| Wide grouped top-1 | https://hypervolu.me/~erik/impg/c4-grouped-multi-seed-wide-top1.png |

The remote upload confirmation is in
`/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/upload_confirmation.txt`.

## Experiment Matrix

I ran `impg crush --method iterative-multi-level` across four
`--window-mode` values:

```text
sibling
sliding
combined
outward
```

Each mode was tested with two window settings:

| Setting | `--window-target-bp` | `--max-window-sites` |
| --- | ---: | ---: |
| compact | 12000 | 4 |
| wide | 30000 | 8 |

For this scheduled evaluation I bounded each run to the top selected candidate
with `--candidate-limit 1 --max-iterations 1`, because an initial
`--candidate-limit 20` probe in
`/home/erikg/impg/data/c4_grouped_multi_bubble_seed_20260604T120526Z`
started with full 465-path complete-homologous windows and was not practical
for the whole matrix. The completed matrix kept exact path preservation as the
hard gate and did not change code.

## Crush Results

All eight completed crush runs exited 0. Every run resolved one candidate and
bailed zero candidates. In every case, the applied candidate source was
`complete-homologous-window`, not the mode-specific sibling/sliding/outward
window type. The window modes still changed candidate generation counts, but
the top candidate chosen by the current ordering was the same within each
target/site setting.

| Variant | Runtime | Max RSS KB | Generated | Selected source | Resolved | Bailed |
| --- | ---: | ---: | ---: | --- | ---: | ---: |
| sibling compact | 1:10.94 | 3602776 | 9097 | complete-homologous-window | 1 | 0 |
| sliding compact | 1:13.77 | 3562740 | 6946 | complete-homologous-window | 1 | 0 |
| combined compact | 1:18.18 | 3580460 | 9377 | complete-homologous-window | 1 | 0 |
| outward compact | 1:08.52 | 3448776 | 5002 | complete-homologous-window | 1 | 0 |
| sibling wide | 2:07.37 | 6390832 | 18569 | complete-homologous-window | 1 | 0 |
| sliding wide | 1:59.87 | 6042412 | 13823 | complete-homologous-window | 1 | 0 |
| combined wide | 2:12.78 | 5801700 | 19699 | complete-homologous-window | 1 | 0 |
| outward wide | 1:58.41 | 5502176 | 10710 | complete-homologous-window | 1 | 0 |

Within each target/site setting, the final graph metrics were identical across
all four window modes.

## Graph Metrics

Metrics below are from `impg graph-report` on the rendered/final GFAs. The
compact row applies to all compact mode variants; the wide row applies to all
wide mode variants.

| Graph | Segments | Links | Paths | Path steps | Segment bp | bp-weighted cov | Singleton bp | High-cov bp | White p99 | White max | White bridges |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| baseline | 7411 | 10072 | 465 | 3459347 | 234828 | 454.330293 | 2922 | 197414 | 66 | 33624 | 8790 |
| compact grouped | 7340 | 10015 | 465 | 3765901 | 234321 | 455.313327 | 2961 | 196930 | 105109 | 233944 | 208954 |
| wide grouped | 7624 | 10419 | 465 | 4018627 | 233456 | 457.000351 | 2915 | 196113 | 101857 | 233079 | 297324 |

Other reported metrics:

| Graph | Segment whitespace bp total | Duplicate sequence frac | Local repeat context nodes |
| --- | ---: | ---: | ---: |
| baseline | 17110105 | 0.669950 | 583 |
| compact grouped | 17096508 | 0.671390 | 576 |
| wide grouped | 17070754 | 0.680220 | 587 |

## Path Preservation

Exact path preservation passed for all outputs:

| Check | Result |
| --- | --- |
| Expected paths | 465 |
| Observed paths | 465 |
| Missing paths | 0 |
| Extra paths | 0 |
| Path spelling mismatches | 0 |

`validation/render_validation_status.tsv` records exit status 0 for
`graph_report`, `gfasort`, `gfalook`, and `path_compare` for the baseline plus
all eight grouped outputs.

## Verdict

Grouping did not improve the residual C4 underalignment in this evaluation.
The exact path-preservation gate passed, and the compact setting slightly
reduced segment count and segment bp, but both compact and wide grouped outputs
substantially worsened path white-space metrics relative to the seed:

```text
white-space p99: 66 baseline -> 105109 compact / 101857 wide
white-space bridges: 8790 baseline -> 208954 compact / 297324 wide
```

The wide setting reduced total segment bp more than the compact setting, but it
increased segment count, path steps, duplicate sequence fraction, and long
white-space bridges. The mode comparison is therefore inconclusive for the
mode-specific grouping behavior: sibling, sliding, combined, and outward
generated different candidate pools, but the current candidate ordering selected
the same complete-homologous candidate for all modes at a given
window-target-bp/max-window-sites setting.

Validation note: I attempted the required build checks before committing this
docs-only report. `cargo build` failed in the vendored `wfmash-rs` dependency
because `htslib/faidx.h` is not available in this environment. `cargo test
--no-run` failed before tests executed because this checkout is also missing
the `vendor/syng/*.c` sources, and the same `wfmash-rs` header failure was hit
again. Since no code was changed, final exact-path validation used the
standalone path spelling comparator recorded in the experiment output.
