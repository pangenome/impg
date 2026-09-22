# Settle Local Replacement Filtering and Poasta Residual Diagnosis

Task: `settle-local-replacement`

Date: 2026-06-02

## Summary

The exact fixed C4 compound-region failure was primarily raw local wfmash coverage, not a second hidden loss in SweepGA scaffold filtering. On the reproduced 45-sequence C4 region, the canonical raw local PAF has full coverage and the shared SweepGA PAF filter leaves it unchanged: 880 raw records, 690 undirected pairs, 45/45 sequences covered, and the same 123-segment / 30,186 bp seqwish graph.

There was, however, a separate local replacement naming contract bug. The actual debug `combined.fa` from the reproduced fixed region used synthetic names such as `__impg_bubble_path8_0`, and `path_spellings.tsv` only recorded synthetic ID plus length. That drops original PanSN/source path identity before wfmash/SweepGA sees the local records. For the exact C4 region below, this did not explain the measured coverage result because the filter retained all canonical records. As a general local replacement contract, it is wrong: PanSN/haplotype grouping is no longer recoverable from local FASTA/PAF names.

This task fixes that contract in `src/resolution.rs`: local replacement records now keep source/query-compatible visible path names when available. For local subsets, only the range portion of the name is adjusted to the selected interval; no synthetic local suffix is appended. For example, a selected sub-interval of `HG001#1#chr6:100-103` can be emitted as:

```text
HG001#1#chr6:101-102
```

That keeps PanSN/source identity recoverable by SweepGA-style name parsing and keeps local FASTA/PAF/GFA path names semantically compatible with `impg query -o fasta` / `impg query -o gfa` output names.

Poasta is not the main blocker. The residual windows being handed to Poasta late in the C4 after-fix run are heterogeneous repeat/CNV residuals, not clean coherent homologous intervals. The extracted final residual `>272218192>272218467` has 333 traversals, median traversal length 366 bp, max 6,734 bp, and total traversal sequence 987,793 bp. A standalone one-round Poasta pass on that saved subgraph did compress it from 446 segments / 15,819 bp to 102 segments / 6,763 bp while preserving paths, so the writer/global/two-piece affine machinery is mechanically functional. The near-identity late full-run Poasta rounds are better explained by stage/scale/selection: after prior replacement, the frontier repeatedly selects already-fragmented residual child windows with little remaining coherent shared structure.

## Reproduction

Exact fixed C4 compound region from the previous `diagnose-and-fix` task:

```text
/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_after_fix/debug/combined.fa
```

New diagnostic root for this task:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/
```

Variant artifact root:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/local_variants/
```

Poasta residual artifact root:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/
```

## Local Variants

The measured metrics are in:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/local_variants/variant_metrics.tsv
```

| Variant | Raw records | Raw pairs | Raw covered | Final records | Final pairs | Final covered | GFA bp | Singleton bp | Shared bp | Path preserved |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|
| `current_fixed_no_filter` | 880 | 690 | 45/45 | 880 | 690 | 45/45 | 30,186 | 7 | 30,179 | yes |
| `wfmash_internal_filter` | 380 | 190 | 20/45 | 380 | 190 | 20/45 | 128,103 | 95,666 | 32,437 | yes |
| `raw_then_sweepga_filter` | 880 | 690 | 45/45 | 880 | 690 | 45/45 | 30,186 | 7 | 30,179 | yes |
| `sweepga_direct_coupled` | 880 | 690 | 45/45 | 880 | 690 | 45/45 | 30,186 | 7 | 30,179 | yes |

Important details:

- `current_fixed_no_filter` uses the preserved canonical raw local PAF from the prior exact after-fix debug directory and disables the seqwish-tail filter. It is the known fixed local wfmash path.
- `raw_then_sweepga_filter` feeds the same raw PAF through impg's shared SweepGA `filter_generated_paf` / `apply_paf_filter` path with `1:1` mapping, `1:1` scaffold filtering, `scaffold_mass=0`, `min_map_length=1`, and `min_identity=0`. It retains all 880 records.
- `wfmash_internal_filter` used the system `/home/erikg/bin/wfmash` v0.24.1 path. In that binary, `-f` means `--no-filter`, not "filter". Both the `-f` raw/no-filter attempt and the default wfmash-filtered output emitted only 380 records and covered 20/45 sequences. This is not an acceptable replacement contract for the fixed C4 block.
- `sweepga_direct_coupled` proves that an impg direct graph build can reach the same result with the vendored SweepGA/wfmash stack on this exact block. It is not yet the fully expressive contract, because the direct graph CLI couples raw wfmash `-n` to the requested filter mode; it cannot express "raw high multiplicity, then 1:1 SweepGA filter" as cleanly as the replacement path.

## Naming Contract

The user-suspected naming issue was real.

Observed pre-fix diagnostic files:

```text
/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_after_fix/debug/combined.fa
/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_after_fix/debug/path_spellings.tsv
```

The `combined.fa` headers were synthetic:

```text
>__impg_bubble_path8_0
>__impg_bubble_path126_9
```

The sidecar only had:

```text
sequence    length
__impg_bubble_path8_0    30148
...
```

Impact:

- Original source path names were not passed to local wfmash, SweepGA, seqwish, or Poasta replacement inputs.
- PanSN grouping information was therefore unavailable from the local FASTA/PAF names.
- The shared SweepGA PAF filter path uses SweepGA's `filter_config_from_align_cfg`; that filter extracts group prefixes from PAF record names when `#` is present. Synthetic `__impg...` IDs have no `#`, so each local traversal becomes its own group.
- On the exact 45-sequence C4 block, this did not change the observed final PAF because raw-then-SweepGA filtering retained all records. The bug is still important for local compound-region replacement generally, especially wherever SweepGA/wfmash PanSN-aware grouping or haplotype-level grouping is semantically relevant.

Fix:

- `PathRange` now carries `source_path_name: Option<String>`.
- All production candidate constructors that have access to `graph.paths[path_idx].name` populate that field.
- `candidate_named_sequences()` now emits the source/query-compatible visible name, adjusting only the range specification for the selected local interval when possible.
- The in-memory BiWFA path was moved to the same header helper so replacement methods stay consistent.
- Tests added:
  - `resolution::tests::candidate_sequence_headers_preserve_source_path_names_when_available`
  - `resolution::tests::candidate_from_root_interval_carries_graph_path_names`
  - `resolution::tests::local_replacement_visible_names_adjust_only_ranges`

## Contract Answers

1. **Was the previous failure also grouping/scaffold-chain filtering, or purely raw local mapping coverage?**

   For the exact fixed C4 compound region, it was raw local mapping coverage. With the corrected local raw PAF, SweepGA filtering/chaining did not remove records and did not change coverage or the induced graph. The newly verified naming bug means the previous local inputs were also not carrying PanSN/source identity, but the reproduction does not show it as the cause of this exact coverage failure.

2. **Should local compound-region replacement be raw/local wfmash all-vs-all, wfmash `-f`/filter output, SweepGA filter/chaining, then seqwish?**

   The right shape is: raw local all-vs-all wfmash with local identity defaults and high-enough raw multiplicity, then upstream SweepGA plane-sweep/scaffold filtering/chaining, then seqwish induction. It should not rely on "wfmash `-f`" as a filtering stage: in the tested wfmash binary, `-f` is `--no-filter`.

3. **Can/should we use SweepGA directly as the filtering/chaining layer after raw wfmash?**

   Yes. The filtering/chaining layer should be SweepGA's own `PafFilter` path exposed through `filter_generated_paf` / `apply_paf_filter` or through `sweepga_align` when that API is semantically configured for the replacement. impg should not duplicate SweepGA plane-sweep/scaffold logic. The remaining API/CLI caution is that direct `impg graph --aligner wfmash --num-mappings 1:1` currently couples raw wfmash multiplicity to the filter request; the replacement path is better because it can run raw high multiplicity first and filter later.

4. **What exactly is the current Poasta problem?**

   The dominant problem is wrong stage/scale/region selection, with already-fragmented heterogeneous input late in the run. The residual candidates mix very short and long traversal classes in a repeat/CNV context, then later rounds reselect small child residuals after previous Poasta output. The evidence does not point to Poasta's global/two-piece affine parameters or exact graph writer as the main failure: the standalone residual probe compressed substantially and the resolution tests cover global two-piece alignment and exact path normalization.

## Poasta Residual Probe

Full after-fix graph report:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/c4_full_after_fix.graph-report.md
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/c4_full_after_fix.flubble-paths.gfa
```

Extracted residual site:

```text
/home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/site_272218192_272218467/
```

Saved files:

```text
input.fa
path_spellings.tsv
input.gfa
input.summary.tsv
poasta.output.gfa
poasta.stderr.log
poasta.stdout.log
```

Input summary:

| Boundary start | Boundary end | Traversals | Segments | Links | Min len | Median len | Max len | Total len |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| `272218192+` | `272218467+` | 333 | 446 | 479 | 366 | 366 | 6,734 | 987,793 |

Standalone one-round Poasta result from `poasta.stderr.log`:

| Metric | Before | After |
|---|---:|---:|
| Segments | 446 | 102 |
| Segment bp | 15,819 | 6,763 |
| Links | 479 | 138 |
| Path steps | 39,696 | 12,356 |
| Singleton bp | 1,202 | 13 |
| Path preservation | accepted | accepted |

Interpretation:

- This residual site is not a small clean homologous interval. It is a broad root-like repeat/CNV window with a large short-vs-long traversal split.
- Poasta can still compress a broad saved residual when run as one selected root, but the full after-fix run's late near-identity rounds show that repeatedly descending into residual children after prior replacement does not create enough new coherent material to condense.
- The blocker is therefore not a graph-quality gate issue. The only hard rejection remains path sequence preservation, and the relevant fix direction is better residual selection/staging rather than rejecting weak-looking graphs after the fact.

## Validation

Commands run:

```text
rustfmt src/resolution.rs
cargo test --lib candidate_
cargo test --lib resolution::
cargo install --path .
```

`cargo test --lib resolution::` passed: 85 tests passed, 0 failed. The worktree needed temporary local symlinks from `vendor/syng` and `vendor/gfaffix` to the main checkout vendor directories because the WG worktree had empty vendor placeholders; those symlinks were used only for validation and are not part of the committed change.

`cargo install --path .` completed and replaced `/home/erikg/.cargo/bin/impg` and `/home/erikg/.cargo/bin/gfaffix` from this worktree.

`cargo test --bin impg` was not run because no CLI syntax changed.
