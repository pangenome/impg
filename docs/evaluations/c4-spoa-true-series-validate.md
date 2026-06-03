# C4 true SPOA threshold series validation

Task: `c4-spoa-true-series-validate`
Date: 2026-06-03

## Inputs

- True SPOA series: `/home/erikg/impg/data/c4_one_many_initial_spoa_true_series_20260603T193935Z`
- Seed report: `/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/reports/one_many_scaffold0.initial.graph-report.tsv`
- Previous all-candidate 5k/10k run: `/home/erikg/impg/data/c4_one_many_initial_spoa_large_20260603T180335Z`

## Validation results

All required true-series thresholds are present in `series.tsv`: `100`, `250`, `500`, `1k`, `2k`, `5k`, `10k`.

For every threshold:

- GFA exists under `graphs/`.
- Graph report TSV exists under `reports/`.
- Sorted GFA exists under `sorted/`.
- PNG exists under `renders/`.
- Extracted path list has exactly 465 names.
- Path-name diff against the seed path list is 0 bytes.
- Uploaded PNG URL returns HTTP 200 with nonzero content.

HTTP checks:

| Threshold | HTTP | Bytes |
| --- | ---: | ---: |
| 100 | 200 | 1037901 |
| 250 | 200 | 1149411 |
| 500 | 200 | 1164473 |
| 1k | 200 | 1156942 |
| 2k | 200 | 1162555 |
| 5k | 200 | 1153131 |
| 10k | 200 | 1170938 |

## Metric comparison

Key graph-report fields:

| Output | Segments | Links | Paths | Total bp | Bp-weighted cov | Singleton bp | Path jump p99 | White p99 | Bridges >= threshold | White frac | Duplicate sequence frac |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| seed initial | 7411 | 10072 | 465 | 234828 | 454.330293 | 2922 | 5 | 66 | 8790 | 0.156693 | 0.669950 |
| true SPOA 100 | 8119 | 10904 | 465 | 239264 | 445.906923 | 3895 | 3454 | 111591 | 212876 | 0.171538 | 0.714127 |
| true SPOA 250 | 8936 | 11747 | 465 | 243592 | 437.984310 | 4916 | 3721 | 109862 | 213772 | 0.186229 | 0.744181 |
| true SPOA 500 | 8317 | 11132 | 465 | 243630 | 437.915996 | 4745 | 3428 | 111829 | 214187 | 0.186321 | 0.725983 |
| true SPOA 1k | 8317 | 11132 | 465 | 243528 | 438.099414 | 4604 | 3428 | 111829 | 214187 | 0.185980 | 0.725863 |
| true SPOA 2k | 8317 | 11132 | 465 | 243528 | 438.099414 | 4604 | 3428 | 111829 | 214187 | 0.185980 | 0.725863 |
| true SPOA 5k | 8317 | 11132 | 465 | 243528 | 438.099414 | 4604 | 3428 | 111829 | 214187 | 0.185980 | 0.725863 |
| true SPOA 10k | 8387 | 11227 | 465 | 250010 | 426.740826 | 4609 | 3445 | 112207 | 214719 | 0.207074 | 0.726243 |
| all-candidate 5k | 8387 | 11227 | 465 | 250010 | 426.740826 | 4609 | 3445 | 112207 | 214719 | 0.207074 | 0.726243 |
| all-candidate 10k | 8387 | 11227 | 465 | 250010 | 426.740826 | 4609 | 3445 | 112207 | 214719 | 0.207074 | 0.726243 |

Relative to the seed, every true SPOA threshold preserves path names but worsens the quality signals that motivated the check:

- Path-jump p99 increases from `5` to `3428`-`3721`.
- White-space p99 increases from `66` bp to `109862`-`112207` bp.
- Bridges above threshold increase from `8790` to `212876`-`214719`.
- Singleton bp increases from `2922` to `3895`-`4916`.
- Duplicate-sequence fraction increases from `0.669950` to `0.714127`-`0.744181`.

The previous all-candidate 5k and 10k outputs are metric-identical to the true-series 10k output. The true-series 5k output is smaller than that prior all-candidate 5k/10k result, but it remains much worse than the seed on path-jump and white-space metrics.

## Code patch validation

Base patch commit: `cb94d10dfdf4d21ce5022a2c86196b518f73df48` (`Respect direct SPOA crush length budgets`).

This validation also updated one stale unit-test expectation:
`resolution::tests::all_candidates_processed_regardless_of_budget` still asserted that explicit direct POA ignored traversal budgets. It now asserts that direct POA skips the over-budget candidate while still processing an eligible candidate in the same graph.

Validated with:

```bash
env C_INCLUDE_PATH=/home/erikg/htslib-local/include CPLUS_INCLUDE_PATH=/home/erikg/htslib-local/include CPATH=/home/erikg/htslib-local/include LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib LD_LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib CFLAGS=-I/home/erikg/htslib-local/include CXXFLAGS=-I/home/erikg/htslib-local/include LDFLAGS=-L/home/erikg/htslib-local/lib CMAKE_PREFIX_PATH=/home/erikg/htslib-local CMAKE_INCLUDE_PATH=/home/erikg/htslib-local/include CMAKE_LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib cargo test direct_poa --lib
```

Result: passed. Cargo ran 6 filtered tests:

- `resolution::tests::direct_poa_respects_median_traversal_budget`
- `resolution::tests::direct_poa_respects_max_traversal_budget`
- `resolution::tests::direct_poa_processes_candidate_within_median_traversal_budget`
- `resolution::tests::direct_poa_processes_candidate_within_max_traversal_budget`
- `resolution::tests::direct_poa_skips_over_budget_candidate_but_processes_eligible_candidate`
- `resolution::tests::direct_poa_and_poasta_accept_topology_changing_replacements`

The fresh worktree needed `git submodule update --init --recursive vendor/syng vendor/gfaffix` before the test could build. A plain test run also failed until local htslib and pinned jemalloc paths were provided to satisfy the vendored `wfmash-rs` build.

Additional validation with the same environment:

- `cargo build`: passed.
- `cargo test --lib`: passed, 374 passed, 0 failed.

## Recommendation

Do not promote any direct true SPOA threshold as an improvement over the seed. All seven thresholds preserve the 465 path names and upload/render correctly, but none improves on the seed; each substantially worsens path-jump p99, white-space p99, bridge count, singleton bp, and duplicate-sequence fraction.

If a SPOA-only diagnostic output is still needed, use the true-series 100 or 1k/2k/5k outputs depending on the metric being inspected, but they should remain diagnostic only. The seed remains the recommended graph.
