# Joint diploid mosaic search: observable dosage recovery

Status: **reviewed and parent-reproduced finite prototype**, 2026-09-15. This is the first automatic test in which both selected primary haplotypes must be non-native mosaics. It is not production diploid calling, exhaustive support or genome-scale validation.

## What problem this solves

A sample MEM-BWT records substring **occurrence counts**, not read placements. For two complete haplotype assignments `G1` and `G2`, the search predicts their integer feature counts, adds them, and evaluates one coupled loss:

\[
q_f=q_f(G_1)+q_f(G_2),\qquad
s_f=10q_f/150,\qquad
R=\sum_f[s_f-C_f\log(1+s_f/0.1)].
\]

The observations and background occur once. Adding two independently computed haploid losses would be a different, incorrect model. Both copies remain complete physical panel routes with independently checked topology, orientation, joins and within-copy source-span capacity. Reusing the same panel source across the two copies is legal and represents dosage.

Earlier search changed one copy at a time. It could find useful single-copy mosaics, but it missed a state requiring compatible evidence to be assembled across both copies. The new operator combines already generated complete physical alternatives before applying the unchanged joint scorer.

## Deterministic mosaic-required test

The transparent development fixture contains two 1,536 bp primary panel alleles, A and B, plus identical 400 bp passive molecules. A and B differ in four 24 bp tracts:

- `[256,280)`
- `[576,600)`
- `[896,920)`
- `[1216,1240)`

The two truth primaries carry A/B masks 10 and 12:

| tract | copy 1 | copy 2 | B dosage |
|---|---:|---:|---:|
| 1 | A | A | 0 |
| 2 | B | A | 1 |
| 3 | A | B | 1 |
| 4 | B | B | 2 |

Both copies therefore contain A and B tracts and are non-native mosaics. Native A/A, A/B and B/B have constant dosage 0, 1 or 2 and cannot explain the sample. The sample has 224 anonymous L150 reads / 33,600 bases, using start stride 15 plus the final start and alternating reverse complements. Inference receives only panel, routes, sample and ordinary fixed options—never the masks, tract coordinates, finite domain or truth.

Assessment alone constructs the 16 legal A-endpoint tract masks plus native B, giving 17 physical haplotypes and 153 unordered pairs. That restricted scan checks the count-vector optimum and all three native exclusions. It is not a certificate over every possible panel route.

## Baseline failure

The unchanged paired search stopped at its bounded `frontier_no_progress` condition:

| measurement | unchanged search |
|---|---:|
| best loss | -1645.1430224277594 |
| restricted finite minimum | -1679.8676690136972 |
| recovered dosage | `[0,1,1,1]` |
| complete paired objectives | 56 |
| public haploid profiles | 112 |

This is a real candidate-generation failure: a better feasible count-exact state existed in the declared physical domain, but search did not construct it.

## Joint physical projection

For each returned complete single-copy candidate, the search examines the existing exploratory alternatives from the same immutable paired baseline:

1. Collect traversal-coordinate breakpoints from their complete routes.
2. Slice generated physical source segments at adjacent breakpoint windows. Reverse segments are sliced in traversal order while retaining their source coordinates and orientation.
3. Project a segment into the other copy at the same route slot and traversal interval.
4. Coalesce only exactly adjacent pieces from the same source and orientation; splitting storage cannot create a new physical copy.
5. Reject invalid orientation seams and publicly validate the projected complete assignment.
6. Evaluate the complete two-copy pair with the unchanged joint scorer.
7. Feed the best completed pair back into the existing baseline/refinement machinery.

The triggering single-copy change need not improve by itself. This permits a jointly useful two-copy state without introducing truth positions, a copying HMM, local read placements, artificial linkage, a family prior or a closest-reference rule.

Each returned candidate scores at most 16 joint projections. Deduplication retains only normalized, publicly valid, score-eligible pairs after successful scoring. Runtime invariants bound the temporary set to 17 full pair keys: the current pair plus 16 projections. Fixed state accounting reserves all 17 keys (626,688 bytes at the configured maximum geometry). No configured scientific or search cap changed.

## Recovery result

| measurement | joint projection |
|---|---:|
| best loss | **-1679.8676690136972** |
| restricted finite minimum | **-1679.8676690136972** |
| retained best pairs | 2 |
| recovered dosage | **`[0,1,1,2]`** |
| retained pairs with exact truth count vector | 2 / 2 |
| retained pairs with two non-native primaries | 2 / 2 |
| complete paired objectives | 262 |
| public haploid profiles | 524 |
| public validations | 2,642 |
| logical peak state | 121,881,016 / 134,217,728 bytes |
| work | 3,993,623 / 4,000,000 |

Native losses are -1272.2600081607065 (A/A), -1625.7156340324018 (A/B) and -1451.4166984878595 (B/B). Their complete count vectors differ from the truth vector and all are worse than the recovered state.

Truth/query denominators are 3,872/3,872 total bp, 3,072/3,072 primary bp and 192/192 tract bp. These are exact fixed-coordinate/whole-route checks, not alignment QV.

## Why phase is not called

The two dosage-one tracts are 320 bp apart, beyond L150 linkage. The selected masks `[8,14]` and the arbitrary truth masks `[10,12]` have the same complete occurrence-count vector. Both retained results recover the observable dosage and likelihood optimum, but neither matches the arbitrary whole-molecule phase.

That is the correct disposition: sequence/dosage supported by the observation operator is recovered, while unsupported long-range phase remains ambiguous. It does not establish that phase is generally unidentifiable; closer informative variants still require correct phase recovery.

## Validation and limitations

Fresh independent review found no issue after one correction: invalid projections had initially been inserted into an unbounded deduplication map. Keys are now retained only after public validation and successful scoring, with explicit 17-key accounting.

Parent reproduction reported:

- focused mosaic dosage: 1 passed;
- paired search: 14 passed, 3 ignored;
- B1: 7 passed, 1 ignored;
- B2: 4 passed, 1 ignored;
- coupling: 1 passed, 1 ignored;
- full portable: **726 passed, 0 failed, 27 ignored**;
- historical exhaustive CIS/TRANS test: unchanged expected failure.

The first parent wrapper exited after all tests because its report requested `meter.work` rather than the actual `meter.work_used`. The failure is preserved; the corrected summary was generated from the unchanged complete outputs and matching pre/post source hashes.

Remaining limitations:

- the operator nearly exhausts the 4M work cap on this small fixture;
- support remains incomplete and the global bound is null;
- the successful fixture is forward-oriented, although existing reverse-orientation regressions pass;
- the 153-pair result is only a restricted finite reference;
- the implementation remains an experimental example, not a production `impg` command;
- no genome-scale, noisy-read, CNV/ploidy, calibrated-confidence or sequence-emission claim follows.

## Code and reproduction

Core implementation:

- `examples/panel_route_diploid_search/mod.rs`
- `examples/panel_route_diploid_search/score.rs`

Fixture and assessment:

- `tests/test_panel_route_diploid_search.rs`
- `tests/panel_route_diploid_search/helpers.rs`

Focused reproduction:

```sh
cargo test --offline --locked --release -j4 \
  --test test_panel_route_diploid_search \
  paired_automatic_mosaic_required_dosage -- \
  --exact --nocapture --test-threads=1
```

Evidence:

- worker red/green: `candidate-copy-repair/core-diploid-dosage-v1/`
- parent reproduction: `candidate-copy-repair/core-diploid-dosage-parent-v1/`
- implementation review: managed output `3d6db6b4-6123-4ec2-ae45-5dc09d747129`
- boundedness correction/review: managed output `77f649fb-f279-44c6-bb56-3285b5cb8a03`

## Next engineering step

Do not extend this global projection operator into the genome-scale architecture. Its role was to prove joint two-copy count arithmetic and a real coupled candidate-generation failure. Follow the [partition-local diplotype inference and physical chaining guide](partition-local-diplotype-inference.md): exhaustively genotype complete spanning allele pairs within bounded syng partitions, use boundary MEMs for physical phase transitions, reconstruct two complete assignments, and retain the unchanged exact whole-genome joint score as authoritative.
