# Partition-observation replay: full-panel controls

## Status

The partition-major implementation is executable, independently reviewed, and
passes 590 parent release tests. Commit `9e93575985a88debf9733b3eccf2cd1dad145331`
was normally installed with `cargo install`; executable SHA256:
`a12d10785cbe80edc26e4f150e9fa18f7859b8e6594a55aae19c6610300b0a42`.

The new error-free controls improve SK1's **alignment-derived QV**, but reduce
whole-truth coverage for both samples. **This is not an unqualified accuracy
improvement, complete genome recovery, or approval for default adoption.**
The old mosaics and successful whole-template comparator remain frozen.

See [implementation and model](partition-observations.md),
[causal diagnosis](observation-registry-repair.md), and
[original sequence controls](whole-genome-sequence-results.md).

## Unchanged inputs and policy

The same panel, AGC, sample MEM-BWTs, S288C reference axis, and complete truth
FASTAs were used. Physical ownership remains exactly once: 19,421 partitions,
383,849 original intervals, 9,901 paths, and 3,336,976,856 source bp. Neither
raw observation orientation nor duplicate descriptions create copies.

Depth 10, background 0.1, maximum mean deviance 10, switch penalty 10, and
`gap-policy=split` were fixed. Both inference outputs and both reconstructions
were frozen before truth evaluation. The original pinned wfmash binary and
one-to-one alignment-selection/evaluation policy were unchanged. All input and
frozen-output checksum checks passed. Runs used CPUs 252–255, nice +10 and at
most four configured threads; commands were serialized without a timeout.

## First all-source, partition986-feature gate

All 2,484 original partition986 features, including zeros, were compiled across
**every panel partition/source**. Compilation took 767.03 seconds, with
1,479,096 KiB peak RSS. All 170 previously diagnosed missing SK1 physical
witnesses were recovered. The crossing witness and both overlapping-MEM
count-two witnesses reproduced their frozen native totals.

SK1's best local bundle changed from `CBM#2` to `SK1#0`; its new source-versus-
second-source score gap was 204.565867. S288C retained the same eight-way tie
including S288C. Old/new objective values are not directly subtracted.

Importantly, **all 166 formerly CBM-exclusive positive features are now shared
between partitions and excluded from local scoring**. The result is not simply
rescoring those observations as local CBM/SK1 co-support. Overall, 1,112 features
were retained and 1,372 excluded: 122 of the original 993 SK1 occurrences remain
in local scoring, with 871 explicitly recorded in excluded factors. Observation
support, exact exposure, and conservative validity changed together.

## Full compilation and scope parity

The full compiler processed all **11,896,777 features**, producing
**232,757,273 distinct physical incidences**. Three duplicate physical records
were merged without summing their occurrence profiles. All 66,051 incidences
from the partition986-feature run are byte-identical in full scope, and the
ownership artifact is identical. All input fingerprints remained unchanged.

- Compilation: **4,239.45 seconds** (1h10m39s); **2,138,072 KiB** peak RSS.
- Exact profile events: 2,174,498,237.
- Core context fetches: 383,849; 3,449,161,256 context bp including read-only halos.
- Reused profiles: 232,198,656; bounded fallback fetches: 557,888.
- Merged incidence stream: 40,486,117,243 bytes.
- Retained output after compilation: 82,213,012,305 bytes / 401,381 files.
  This is retained disk usage, **not measured peak scratch usage**.

Both callers used 8,882,931 local factors and excluded 3,013,846. Exclusion reason
counts overlap: shared-between-groups 2,817,636; context-nonlocal 462,310;
boundary-crossing 97,719; source-terminal-context 69,443. Exclusion is not deletion
or evidence of biological absence. Calling plus axis threading took 322.06 seconds
for S288C and 324.42 seconds for SK1, with respective peak RSS 795,560 and
794,028 KiB. Full-feature compilation is measured here; the small feature-scope
runtime was not extrapolated as a whole-panel performance claim.

## Sequence results: QV and coverage together

| Control | Model | Emitted bp | Assessed error columns | Alignment QV | Whole-truth coverage | Query coverage |
|---|---|---:|---:|---:|---:|---:|
| S288C | Frozen original | 11,482,108 | 0 | zero observed errors | 93.7823% | 99.9965% |
| S288C | Partition replay | 11,289,316 | 0 | zero observed errors | **92.2086%** | 99.9976% |
| SK1 | Frozen original | 11,068,241 | 3,013 | 35.6337 | 90.7458% | 99.5963% |
| SK1 | Partition replay | 10,926,109 | 768 | **41.5311** | **89.9411%** | 99.9988% |

The new SK1 errors are 4 mismatch + 382 insertion + 382 deletion columns, all
within one selected 60,026 bp `DBVPG6044#0#chrVII` block. This is an alignment
column accounting result, not 768 independent variants or proof that the
alignment is optimal. The block covers axis partitions459–464. It remains a
concrete downstream audit target, not a donor-specific repair instruction.

The frozen-call audit narrows this further: SK1 and DBVPG6044 have exactly tied
local scores in all six partitions, and both sets of source intervals are
contiguous. SK1 has known `+` orientation in partitions459–463, but
`conflicting-signed-contexts` / unknown orientation in partition464; DBVPG6044
retains known `+` and is the unique selected DP state throughout the block.
Thus the retained local likelihood does not favor the erroneous donor; the
orientation/continuity constraints distinguish these alternatives.

An exact downstream string comparison finds equal flanks surrounding a **415 bp
reverse-complement difference** (block offsets 55,475–55,890; SK1 truth
861,167–861,582). This explains why feature-polarity conflicts need to be
distinguished from whole-member orientation. It does not authorize overriding
that state with truth, choosing a majority orientation, or changing the frozen
alignment-derived error accounting. Direct core-shard evidence gives SK1 469
same-polarity and eight reverse-polarity features versus the reference; all eight
reverse witnesses overlap the internal variant or its boundary context (up to
21 bp beyond the minimal differing window). DBVPG6044 has 477 same-polarity
features and no reverse witnesses. This is not a majority-vote orientation rule.

A stronger local observability check enumerates every 150 bp DNA window in the
two 60,026 bp spans, canonicalizes reverse-complement orbits, and preserves start
multiplicity. The resulting **canonical DNA multisets are exactly equal**:
59,877 starts and 59,817 distinct canonical reads each, with zero difference in
either direction. Consequently, the internal variant is not distinguishable by
these local strand-invariant 150 bp observations alone. This is not a claim of
whole-genome indistinguishability; source/haplotype linkage assumptions remain
separate. More permissive thresholds cannot create missing local information.

S288C has 11,289,047 exact assessed columns and 953,895 unaligned truth bp.
SK1 has 10,926,361 assessed columns, 10,925,593 matches and 1,221,944 unaligned
truth bp. Full truth denominators remain 12,242,942 and 12,147,923 bp,
respectively, including 17,357 non-ACGT bp in each. Unaligned query bp are 269
and 130. Zero observed errors has null numeric QV, not a claim of a perfect genome.

S288C emits 147 fragments with 165 unresolved intervals; SK1 emits 170 fragments
with 198 unresolved intervals. No source-gap imputation or cross-block source
reuse was reported. Zero split-query structural diagnostics does not establish
chromosome-scale correctness. More exclusions and fragmentation can trade
coverage for conditional accuracy; the coverage decline must not be hidden.

Compared with the frozen original, 29 S288C intervals newly become unresolved
and seven become emitted; the corresponding SK1 counts are 27 and nine. Newly
unresolved no-positive-local-evidence intervals account for 21 / 181,397 axis bp
in S288C and 17 / 168,082 axis bp in SK1. Other new unresolved cases include
nonidentical alternatives and orientation uncertainty. These axis bp are not
additive missing-truth bp. The 111 repeated-axis appearances covering 553,177
axis bp remain excluded in both controls.

## Evidence and remaining gates

Local immutable run roots under `/home/erikg/yeast/`:

- `partition986-observation-repair-20260911T224838Z/`: input checksums, incidence
  audit, frozen calls/comparison, and downstream exclusion breakdown.
- `genome-partition-observations-20260911T231459Z/`: full compiler manifest,
  resources, complete shards/registry, and streamed scope-parity validation.
- `genome-partition-replay-controls-20260912T003055Z/`: frozen inference and
  reconstructions, native PAF/base-replay evaluations, whole-denominator summaries,
  and downstream coverage/residual-error audits.

The parent test, timeout-recovery, independent-review and normal-install records
are retained under `target/experiments/genome-mem-bwt-pipeline/candidate-copy-repair/`.
PR #243 remains draft and unmerged; its head was checked at `9e93575` with no
reported status checks, which is not a CI pass.

Next decisions require explaining the remaining error block and newly unresolved
regions without tuning to truth or relaxing conservative validity. General shared-
factor/junction inference, repeated/off-axis copy placement, strict held-outs,
recombinants, noisy reads, and coupled ploidy remain separate gates. The successful
observation compiler does not certify partition membership as homology or dosage.
