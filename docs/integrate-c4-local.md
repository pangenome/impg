# Integrate C4 Local Replacement Fixes

Task: `integrate-c4-local`

Branch state: the active PR branch already contained tree-equivalent versions of the two requested commits:

- `eb3d937` via local branch commit `3393d5d`
- `131f10e` via local branch commit `e8060bb`

This integration tightened the naming contract so local replacement FASTA/PAF/GFA-visible names remain semantic source names. Local subsets now adjust only the trailing source range where possible. Synthetic `__impg_bubble_*` names are no longer appended to visible local replacement headers.

## Validation

Focused tests:

- `cargo test --lib candidate_sequence_headers_preserve_source_path_names_when_available`
- `cargo test --lib local_replacement_visible_names_adjust_only_ranges`
- `cargo test --lib candidate_from_root_interval_carries_graph_path_names`
- `cargo test --lib local_wfmash`
- `cargo test --lib compound_seqwish_tail_keeps_raw_many_to_many_condensation_path_correct`
- `cargo test --lib wfmash_raw_mapping_multiplicity_matches_local_filter_mode`
- `cargo test --bin impg`
- `cargo test`

All passed. The full test suite includes the long C4 slice path-preservation integration tests.

Installed binary:

- `/home/erikg/.cargo/bin/impg`

The install required seeding the `wfmash-rs` release build-script output with the compatible local `d47b7e3` wfmash binary because the vendored C++ build still cannot find `htslib/faidx.h` in this environment.

## Exact Region

Command output:

- GFA: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/exact_region/exact_region.seqwish.gfa`
- Debug directory: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/exact_region/debug`

Debug artifacts include:

- `combined.fa`
- `path_spellings.tsv`
- `raw.paf`
- `filtered.paf`
- `final.paf`
- `raw.coverage.tsv`
- `filtered.coverage.tsv`
- `final.coverage.tsv`
- `raw.metrics.tsv`
- `filtered.metrics.tsv`
- `final.metrics.tsv`

Exact-region metrics:

| metric | value |
|---|---:|
| final PAF records | 880 |
| covered sequences | 45 / 45 |
| zero-coverage sequences | 0 |
| min coverage fraction | 1.000000 |
| mean coverage fraction | 1.000000 |
| replacement segments | 123 |
| replacement segment bp | 30,186 |
| singleton bp | 7 |

This matches the fixed behavior and does not regress to the old 64,206 bp / 34,049 singleton bp replacement.

## Full C4

Input:

- `/home/erikg/impg/data/explore_c4_parent_20260601T212246Z/masked_syng_control/run.nosort.gfa`

Outputs:

- Final GFA: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfa`
- Graph report TSV: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.graph-report.tsv`
- `gfalook -m` PNG: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfalook-m.png`
- Crush stderr log: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/crush.stderr.log`
- Debug directory: `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/debug`

Full-run debug artifacts include:

- `debug/graph_build_0000/combined.fa`
- `debug/graph_build_0000/path_spellings.tsv`
- `debug/graph_build_0000/raw.paf`
- `debug/graph_build_0000/filtered.paf`
- `debug/graph_build_0000/final.paf`
- `debug/graph_build_0000/raw.coverage.tsv`
- `debug/graph_build_0000/filtered.coverage.tsv`
- `debug/graph_build_0000/final.coverage.tsv`
- `debug/replacement_0000_sweepga_seqwish/seqwish.gfa`
- `debug/replacement_0000_sweepga_seqwish/unchopped.gfa`

The round-1 local wfmash replacement used local multiplicity and identity as intended:

- `-n117`
- identity `70`
- 117 local sequences

The full crush resolved 8 candidates over 8 configured rounds with no bails. Round 1 condensed the previously bad large local site from 44,848 input bp / 4,972 singleton bp to 33,278 output bp / 376 singleton bp. Later small residual rounds did not recreate the old singleton expansion at that site.

Comparison using `impg graph-report` fields:

| run | segments | links | segment bp | bp-weighted coverage | singleton bp | long whitespace bridges |
|---|---:|---:|---:|---:|---:|---:|
| PGGB baseline | 13,288 | 16,240 | 234,524 | 454.919 | 2,890 | 11,876 |
| previous failed largest-iter SYNG crush | 19,051 | 22,137 | 464,409 | 229.732 | 79,905 | 163,276 |
| integrated full C4 crush | 19,026 | 22,085 | 462,940 | 230.461 | 78,435 | 163,336 |
| starting masked SYNG control | 18,970 | 21,863 | 425,742 | 250.597 | 81,420 | 168,197 |

The integrated run improves modestly over the previous failed largest-iter SYNG crush: 25 fewer segments, 1,469 fewer segment bp, and 1,470 fewer singleton bp. It also improves singleton bp versus the starting masked SYNG control by 2,985 bp. It is still not PGGB-like: segment bp remains about 2x PGGB, singleton bp remains much larger than PGGB, and long whitespace bridges remain high.

## Why This Was Missed

Earlier sequence-preservation tests verified that local replacement graph paths spelled the right sequences, but they did not require visible source-name semantics to survive into local FASTA/PAF/GFA names. That meant alignment and grouping filters received valid sequence records with bad metadata. For C4, losing PanSN/source identity can change grouping and filtering behavior even when path sequence preservation still passes.
