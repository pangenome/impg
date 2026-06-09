# Implement Anchored Syng Interval Emission

Task: `implement-anchored-syng-2`
Date: 2026-06-09

## Code Changes

- Preserved syng anchor provenance through the GFA/FASTA local-seed interval path by carrying `HomologousIntervalWithAnchors` through transitive query collection.
- Added anchored interval emission in `src/syng.rs`. Emission trims fetch spans to anchor-supported target coordinates plus configured padding/query extension, and only merges nearby same-path intervals when target and query anchor gaps are collinear.
- Kept boundary-refined anchored outputs split, but restored the legacy coordinate-merged per-hop frontier for transitive traversal. This lets subsequent hops discover the same coordinate span as the old query path while still emitting only anchor-qualified fetch intervals.
- Routed `gfa:syng-local` and syng FASTA output through anchored emission. Other syng query outputs keep the existing coordinate-only merge behavior.
- Added `syng_local_intervals.tsv` debug output for local-seed GFA/FASTA collection with:
  `path_key`, `pre_merge_interval`, `post_merge_interval`, `strand`, `anchor_count`,
  `query_anchor_min`, `query_anchor_max`, `target_anchor_min`, `target_anchor_max`,
  `emitted_fetch_interval`, and `reason`.

## Synthetic Tests

New tests in `src/syng.rs` cover:

- Trimming an unanchored terminal suffix.
- Trimming an unanchored left prefix.
- Splitting two target blocks separated by a 32,738 bp target-only gap.
- Splitting two target blocks separated by a 65,476 bp target-only gap.
- Merging when query and target gaps are both collinear.
- Keeping anchored output candidates split while using coordinate-merged spans for transitive BFS traversal.

Targeted validation:

```text
cargo test test_anchored_emission --lib
5 passed
cargo test anchored_frontier_merges_for_traversal_without_collapsing_output_hits --lib
1 passed
```

## C4 Rerun

Output tree:

```text
/home/erikg/impg/data/implement_anchored_syng_2_frontier_20260609T224925Z
```

Command profile: the same C4 one-chunk localized syng-local profile used by the diagnosis pass, with `--debug-dir` enabled for both the pipeline and localized polish stages.

Run result:

- Exit status: 0
- Wall clock: 6:51.06
- Syng-local interval collection: 6,793 sequences, 26,166,770 bp
- Local seed graph build: `path_validation=pass`
- Graph report: 4,244 segments, 2,474 links, 6,793 paths, 890,777 path steps
- Graph report status: `REVIEW` for diagnostic topology metrics only; exact path spelling remained valid
- Debug TSV: `/home/erikg/impg/data/implement_anchored_syng_2_frontier_20260609T224925Z/debug/c4_tuned_pipeline/syng_local_intervals.tsv`
- Path spellings: `/home/erikg/impg/data/implement_anchored_syng_2_frontier_20260609T224925Z/debug/c4_tuned_pipeline/graph_build_0000/path_spellings.tsv`
- Render: `/home/erikg/impg/data/implement_anchored_syng_2_frontier_20260609T224925Z/renders_tuned/C4_GRCh38_53kb.png`
- Uploaded PNG: http://hypervolu.me/~erik/impg/implement-anchored-syng-2-C4_GRCh38_53kb-frontier-20260609T224925Z.png
- Upload verification: HTTP 200, `Content-Type: image/png`, `Content-Length: 1324813`

Generated GFA/GFA-derived debug/render/log artifacts are under `/home/erikg/impg/data/...` and were not committed.

## Before/After Selected Intervals

The before intervals are from `docs/evaluations/diagnose-c4-query-2.md`.
The after intervals are the emitted fetch intervals from the new debug TSV and match `path_spellings.tsv`.

| Path | Before selected interval | After emitted intervals | Result |
| --- | --- | --- | --- |
| `NA18952#1#JBHIJS010000065.1` | `31855636-32040012`, 184,376 bp | 31 intervals, 122,156 bp: `31887445-31888937`; `31889006-31891078`; `31891268-31892206`; `31892415-31894486`; `31894564-31896540`; `31896802-31899871`; `31900646-31905294`; `31905448-31907091`; `31907193-31908874`; `31909190-31921675`; `31921744-31924944`; `31925153-31927224`; `31927302-31929278`; `31929540-31932609`; `31933384-31938032`; `31938186-31939829`; `31939931-31941612`; `31941928-31954413`; `31954482-31959962`; `31960040-31962016`; `31962278-31965347`; `31966122-31970770`; `31970924-31971540`; `31971787-31972567`; `31972669-31974350`; `31974666-31992700`; `31992778-31994754`; `31995016-32004278`; `32004525-32005305`; `32005407-32007088`; `32007404-32016729` | Corrected: the single 184 kb long-tail fetch is split into anchor-bounded intervals; terminal tail is not absorbed into one selected span. |
| `NA18948#2#CM101588.1` | `31979894-32164269`, 184,375 bp | 30 intervals, 121,960 bp: `32011703-32016464`; `32016673-32018744`; `32018822-32020798`; `32021060-32024129`; `32024904-32029552`; `32029706-32033132`; `32033448-32048074`; `32048264-32049202`; `32049411-32051482`; `32051560-32053536`; `32053798-32056867`; `32057642-32063060`; `32063307-32064087`; `32064189-32065870`; `32066186-32080812`; `32081002-32081940`; `32082149-32084220`; `32084298-32086274`; `32086536-32089605`; `32090380-32095028`; `32095182-32096825`; `32096927-32098608`; `32098924-32111409`; `32111478-32113550`; `32113740-32116958`; `32117036-32119012`; `32119274-32128536`; `32128783-32129563`; `32129665-32131346`; `32131662-32140986` | Corrected: the single 184 kb long-tail fetch is split into anchor-bounded intervals. |
| `HG02392#1#CM099882.1` | `31884401-32056195`, 171,794 bp | 28 intervals, 106,988 bp: `31916948-31917754`; `31917823-31919895`; `31920216-31921023`; `31921232-31923304`; `31923382-31925512`; `31925620-31928569`; `31928689-31934112`; `31934266-31934882`; `31935129-31937545`; `31938008-31950493`; `31950562-31954940`; `31955060-31956267`; `31956716-31958184`; `31958227-31963915`; `31964378-31976863`; `31976932-31977937`; `31978056-31978714`; `31979325-31982413`; `31982491-31984621`; `31984729-31987678`; `31987798-31989005`; `31989454-31996654`; `31997117-32007719`; `32009084-32009762`; `32009888-32015778`; `32015922-32021485`; `32021587-32026764`; `32029073-32032912` | Corrected: the single 171 kb long-tail fetch is split/trimmed by anchor support. |
| `HG01960#2#JBHIHN010000011.1` | One-chunk baseline: `32045898-32158486`, 112,588 bp. Moderate clipped case: `32078636-32158486`, 79,850 bp | 13 emitted intervals, 51,962 bp: `32079320-32081611`; `32081713-32082520`; `32082729-32084801`; `32084879-32087009`; `32087117-32090066`; `32090186-32091393`; `32091842-32099042`; `32099505-32111990`; `32112059-32116754`; `32118213-32121979`; `32122133-32122749`; `32122996-32125412`; `32125875-32135203`. Debug TSV pre-merge now includes the baseline-left coordinate in the first row: `32045898-32081833`, with target anchors only at `32079440-32081428`. | Partially corrected at provenance/traversal level: the baseline-left coordinate is now preserved in `pre_merge_interval`, proving transitive collection can see it. The emitted fetch remains clipped by design because the left span is unanchored under the retained chain, so the visual left clip is not fully corrected. |

## Interpretation

The high-confidence 3-prime long-tail candidates are corrected under the anchored emission rule: they no longer appear as one coordinate-merged fetch interval across repeated C4 blocks. Rows with `split_chain_discontinuous` in the debug TSV mark same-path blocks that are within `-d` but whose query/target anchor gaps are not collinear.

The HG01960-style left clip is only partially corrected. The frontier patch restores the baseline-left coordinate in `pre_merge_interval`, but the final emitted fetch remains close to the clipped start. The debug TSV gives the immediate reason: the restored left span lacks retained anchor support in the syng chain reaching `HG01960#2#JBHIHN010000011.1`, so the anchor-provenance merge/trim rule does not reintroduce it into the fetched sequence. Follow-up `investigate-hg01960-left` is queued to decide whether that anchorless left span should be recovered by a separate boundary-continuity rule or is genuinely unsupported.

## Validation

- `cargo test test_anchored_emission --lib` passed.
- `cargo test query_scaled_chain_limits_use_merge_distance_before_refinement --lib` passed.
- `cargo build` and full `cargo test` passed for the main anchored-emission commit.
- After the frontier patch, targeted tests and `cargo build --release` passed using `CARGO_TARGET_DIR=/home/erikg/impg/target` because the WG cleanup pass had removed this worktree's local `target/`.
- Final `cargo build`, full `cargo test`, `cargo install --path .`, `git diff --check`, and staged artifact/blob scan were rerun before the follow-up commit; see the WG task log for final command results.
