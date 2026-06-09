# Investigate HG01960 Left-Span Anchor Support

Task: `investigate-hg01960-left`
Date: 2026-06-09

## Question

`implement-anchored-syng-2` fixed anchored syng interval emission for C4
long-tail paths, but `HG01960#2#JBHIHN010000011.1` still did not return to
the older one-chunk baseline left coordinate. This pass checks whether the
missing `32045898-32078636` span is lost during transitive collection, chain
filtering, boundary refinement, or whether it lacks syng anchor support under
the current C4 query.

The validated `implement-anchored-syng-2` C4 run was:

```text
/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z
```

The recorded command did not pass `--transitive`, so the syng-local path used
`syng_max_depth = 1`. That rules out cross-hop transitive frontier collection
as the cause for this C4 row.

## Current Emitted Intervals

Reproduced from:

```text
/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z/debug/c4_tuned_pipeline/syng_local_intervals.tsv
```

`HG01960#2#JBHIHN010000011.1` emitted intervals:

| Source | Intervals | Total bp |
| --- | --- | ---: |
| One-chunk baseline | `32045898-32158486` | 112,588 |
| Moderate clipped case | `32078636-32158486` | 79,850 |
| `implement-anchored-syng-2` emitted | `32079320-32081611`; `32081713-32084801`; `32084879-32087009`; `32087117-32091393`; `32091842-32099042`; `32099505-32116754`; `32118213-32125412`; `32125875-32135203` | 52,761 |

The first retained emitted row has:

```text
pre_merge_interval=32078636-32081833
anchor_count=48
query_anchor_min=31982860
query_anchor_max=32017146
target_anchor_min=32079440
target_anchor_max=32081428
emitted_fetch_interval=32079320-32081611
```

So the final C4 output is not restored to the baseline left coordinate. It
begins near the clipped coordinate because the retained chain anchors start at
`32079440`.

## Raw Syng Collection Probe

To separate raw syng support from SweepGA chain filtering and BiWFA boundary
refinement, I reran the same query as FASTA with `--syng-raw`, `-d 100k`, the
same HPRC syng index, the same AGC sequence source, and a debug directory:

```text
/home/erikg/impg/target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --syng-raw \
  --debug-dir /tmp/hg01960_syng_raw_probe \
  -o fasta \
  -O /tmp/hg01960_syng_raw_probe/out
```

This bypasses transitive depth, chain filtering, and boundary realignment.
For `HG01960#2#JBHIHN010000011.1`, raw current-settings collection produced:

| Interval | Anchor count | Query anchors | Target anchors | Emitted interval | Reason |
| --- | ---: | --- | --- | --- | --- |
| `32077322-32103548` | 1,435 | `31982236-32035323` | `32077442-32103365` | `32077322-32103548` | `anchor_supported` |
| `32103913-32106095` | 20 | `32007452-32009331` | `32104033-32105912` | `32103913-32106095` | `merged_chain_continuous` |
| `32106271-32113998` | 300 | `31982089-32017234` | `32106391-32113815` | `32106271-32113998` | `split_chain_discontinuous` |
| `32114048-32116874` | 155 | `31991216-32026477` | `32114168-32116691` | `32114048-32116874` | `split_chain_discontinuous` |
| `32117332-32127347` | 500 | `31994500-32035323` | `32117452-32127164` | `32117332-32127347` | `split_chain_discontinuous` |
| `32127831-32136156` | 125 | `32005001-32012905` | `32127951-32135973` | `32127831-32136156` | `merged_chain_continuous` |

Under the current default syng query settings, the earliest raw target anchor
on this path is `32077442`. That is already about 31.5 kb to the right of the
one-chunk baseline start `32045898`.

## Seed-Filter Sensitivity Probe

I also repeated the raw probe with seed filtering disabled:

```text
/home/erikg/impg/target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31982056-32035418' \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --syng-raw \
  --syng-seed-drop-top-fraction 0 \
  --syng-seed-walk-anchors 1 \
  --debug-dir /tmp/hg01960_syng_raw_nofilter_probe \
  -o fasta \
  -O /tmp/hg01960_syng_raw_nofilter_probe/out
```

The debug TSV was written before the huge FASTA stream completed; I stopped
the scratch FASTA write after extracting the TSV rows. With the seed filter
disabled, `HG01960#2#JBHIHN010000011.1` still started in the same right-shifted
neighborhood:

| Interval | Anchor count | Query anchors | Target anchors | Emitted interval | Reason |
| --- | ---: | --- | --- | --- | --- |
| `32077003-32136360` | 3,055 | `31982089-32035323` | `32077123-32136177` | `32077003-32136360` | `anchor_supported` |

This proves the broad missing baseline-left block is not caused by the default
top-frequency seed drop or five-syncmer walk seeding. Even with single-syncmer
seeds and no top-frequency drop, the earliest target anchor is `32077123`, not
near `32045898`.

## Stage Finding

The missing one-chunk baseline span is lost before chain filtering and before
boundary refinement:

| Stage | Earliest HG01960 support/interval |
| --- | --- |
| Raw syng collection, current settings | interval `32077322-32103548`, first target anchor `32077442` |
| Raw syng collection, seed filter disabled | interval `32077003-32136360`, first target anchor `32077123` |
| Boundary-realigned/chained C4 output | first emitted interval `32079320-32081611`, first retained target anchor `32079440` |

There is a smaller post-collection effect: chain/scaffold filtering removes the
earliest raw anchors between about `32077123` and `32079440`, and boundary
refinement then emits from `32079320`. That accounts for roughly 2.3 kb near
the clipped edge, but it does not explain or recover the larger
`32045898-32077003` baseline-left span.

The larger left span is unsupported by the current syng index/query evidence
for `GRCh38#0#chr6:31982056-32035418`. Restoring `HG01960#2#JBHIHN010000011.1`
to `32045898` would require reintroducing unanchored coordinate carryover, the
behavior that caused the C4 long-tail artifacts fixed by anchored emission.

## Outcome

No code change was made. The current anchored emission behavior is conservative
for this path: it does not restore `HG01960#2#JBHIHN010000011.1` to the
one-chunk baseline left coordinate because no syng anchors support that left
span under either current seed settings or a no-seed-filter raw query.

Because no code changed, no new synthetic regression test was added.

## Validation Notes

- Reproduced the current `implement-anchored-syng-2` emitted intervals from
  `syng_local_intervals.tsv`.
- Identified the loss point as raw syng collection / absent direct syng support
  for the broad baseline-left span, with only a small additional chain-filter
  trim near `32077k-32079k`.
- The initial fresh-target `cargo test effective_min_chain_anchors_scales_with_query_span --lib`
  attempt failed while compiling `wfmash-rs` because the local worktree lacks
  `htslib/faidx.h`; final validation was rerun using the shared project target
  cache.
- `CARGO_TARGET_DIR=/home/erikg/impg/target cargo build` passed after
  initializing the vendored `syng` and `gfaffix` submodules in this worktree.
- `CARGO_TARGET_DIR=/home/erikg/impg/target cargo test` passed with no
  regressions.
