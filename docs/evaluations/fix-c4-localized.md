# Fix C4 Localized Polish Degradation

Task: `fix-c4-localized`

## Diagnosis

The C4 degradation in `run_fuller_c4_20260609T160000Z` was caused at the
localized-polish chunk selection boundary, before any graph-quality diagnostic
metric was consulted. `src/localized_polish.rs::select_chunks` selected dirty
chunks independently from a pre-iteration report, then
`resolve_gfa_dirty_chunks` resolved those chunks sequentially against a graph
that was being rewritten after each accepted replacement. That made adjacent
dirty sites stale after the first replacement and let later chunks reprocess
already-polished interiors.

Concrete old-run degradation boundaries:

| area | old chunks / ranges | observed effect |
|---|---|---|
| left frontier, steps 45-94 | `chunk_0029..chunk_0031` on `GRCh38#0#chr6:31949370-32068272`, path `1252..6317`, orders `44..2749`, nodes `45,47,67,69,77,79,92,93` | Applied as three separate SPOA replacements. Segments rose `3872 -> 3897`, segment bp `98149 -> 98843`, compression fell `541.685804 -> 537.882501` with no path-step improvement. |
| left CHM13 frontier | `chunk_0000` path `4292..6316`, orders `91..92`, nodes `92,93`; `chunk_0001` path `5080..7095`, order `100`, node `101` | Neighboring flanked windows were selected independently instead of as one anchored compound window. |
| larger early span, steps 188-586 | `chunk_0039` on `GRCh38#0#chr6:31949370-32068272`, path `27003..31446`, orders `538..2754`, nodes including `539,541,550,552,561,563,567,569,576,577,584,586` | One SPOA replacement increased segments `3948 -> 3995`, segment bp `101073 -> 103517`, links `5335 -> 5384`, while compression fell `526.015058 -> 513.596028`. |
| late/right span | `chunk_0072`, path `95407..97551`, orders `1858..2661`, nodes `1859,1861,1867,1869,1985,2662` | Iteration 1 accepted another exact-path-preserving replacement that increased segments `4120 -> 4158`, segment bp `108234 -> 110122`, and reduced compression `491.212743 -> 482.791086`. |
| late/right candidate band around 3243-3748 | examples include `chunk_0191` orders `2037..3759`, `chunk_0280` orders `86..3661`, `chunk_0739` orders `444..3847`, `chunk_0812` orders `1487..3751` | These were present as broad crossing candidate windows in the strong dirty report. The old selector would eventually feed the same stale split-window/reprocessing pattern into these bands as budgets increased. |
| duplicate interiors | `chunk_0087` duplicates `chunk_0029` nodes on `HG00097#1`; `chunk_0097` duplicates `chunk_0039`; `chunk_0120` and `chunk_0127` later re-hit already rewritten interiors | Several were reported `applied` with before/after graph metrics unchanged, proving the loop spent budget reprocessing provenance-equivalent interiors. |

The resolver was still enforcing the correct hard gate: exact path names and
spellings. The bad graph shape came from accepting exact-preserving but stale,
partial replacements. Graph metrics remain diagnostic only.

## Implementation

The fix is in `src/localized_polish.rs`.

- Before budget selection, same-path dirty chunks with overlapping or nearby
  flanked path spans are merged into one `DirtyChunk` compound window. The
  compound chunk preserves the union of path/core spans, site ids, kinds, nodes,
  severity, flank, and diagnostic ids.
- Selection now tracks dirty-window provenance. A chunk is deferred if its core
  interval or node set is already covered by an applied chunk from an earlier
  iteration or by a chunk selected earlier in the same iteration. This prevents
  the HG00097-style duplicate interior reprocessing seen in the strong C4 run.
- Replacement acceptance is unchanged. Exact path corruption remains the hard
  rejection gate; graph-quality metrics are not used to accept or reject a
  replacement.

## Regression Tests

Added synthetic unit regressions in `src/localized_polish.rs`:

- `localized_polish_compounds_overlapping_dirty_windows_before_resolution`
- `localized_polish_defers_duplicate_dirty_interiors_by_provenance`

The first test was run before the implementation and failed with the old
selector returning two windows instead of one:

```text
assertion failed: overlapping flanked windows from one scar should be resolved as a single anchored compound window
left: 2
right: 1
```

After the fix, both tests pass.

## C4 Rerun

Post-fix artifacts are outside the repository:

`/home/erikg/impg/data/fix_c4_localized_20260610T004245Z`

Both C4 reruns used this worktree's `target/release/impg` and the same C4
query profile as the anchored syng rerun:

- syng index: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`
- AGC: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`
- BED: `/home/erikg/impg/data/fix_c4_localized_20260610T004245Z/c4.bed`
- common query settings: `-d 100k`, `--num-mappings 1:many`,
  `--scaffold-filter 1:many`, `--scaffold-jump 0`, `--min-match-len 311`,
  `-t 16 -v 1`

Localized configs:

| run | config |
|---|---|
| one-chunk baseline | `iterations=1,max-chunks-per-iteration=1,max-total-chunks=1,max-total-bp=200k,max-runtime-secs=300,flank=1k,method=poa,polish-rounds=0` |
| expanded | `iterations=2,max-chunks-per-iteration=32,max-total-chunks=64,max-total-bp=2m,max-runtime-secs=1800,flank=1k,method=poa,polish-rounds=0` |

Runtime and localized-polish summary:

| run | exit | wall | RSS KB | collected paths | initial candidates | selected | selected bp | applied | final status | exact paths |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---|---|
| one-chunk baseline | 0 | `5:57.80` | 48,384,436 | 6,719 | 6,719 | 1 | 1,492 | 0 | `stalled-no-change` | pass |
| expanded 2x32/64 | 0 | `6:14.88` | 48,395,384 | 6,659 | 6,659 | 32 | 71,982 | 0 | `stalled-no-change` | pass |

Within-run graph diagnostics did not degrade because no localized replacement
was applied:

| run | segments before/after | path steps before/after | singleton bp before/after | white-space bp before/after |
|---|---:|---:|---:|---:|
| one-chunk baseline | `4445 / 4445` | `859842 / 859842` | `4508997 / 4508997` | `6990356186 / 6990356186` |
| expanded 2x32/64 | `4419 / 4419` | `871731 / 871731` | `4496684 / 4496684` | `7089234729 / 7089234729` |

PNG/render outputs, useful for visual comparison:

- baseline:
  `/home/erikg/impg/data/fix_c4_localized_20260610T004245Z/baseline/renders/C4_GRCh38_53kb.png`
- expanded:
  `/home/erikg/impg/data/fix_c4_localized_20260610T004245Z/expanded/renders/C4_GRCh38_53kb.png`

Result: expanded localized polish is non-degrading after the selection fix, but
the current anchored syng C4 graph still does not improve because the selected
frontier chunks all lack both shared flank anchors. This is a remaining
anchor-availability limitation, not the old exact-preserving scar introduced by
stale split-window reprocessing.

## Validation

Commands run:

```bash
source ./env.sh >/tmp/impg_env.log
CARGO_BUILD_JOBS=8 cargo test localized_polish_ --lib
CARGO_BUILD_JOBS=8 cargo test --test test_localized_polish --test test_localized_resolver
CARGO_BUILD_JOBS=8 cargo build --release
```

Full `cargo build`, full `cargo test`, `git diff --check`, and artifact/blob
scans were run after this report was written and are recorded in the task log.
