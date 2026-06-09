# Implement Anchored Syng Interval Emission

Task: `implement-anchored-syng-2`
Date: 2026-06-09

## Code Changes

- Preserved syng anchor provenance through the GFA/FASTA local-seed interval path by carrying `HomologousIntervalWithAnchors` through transitive query collection.
- Added anchored interval emission in `src/syng.rs`. Emission trims fetch spans to anchor-supported target coordinates plus configured padding/query extension, and only merges nearby same-path intervals when target and query anchor gaps are collinear.
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

Targeted validation:

```text
cargo test test_anchored_emission --lib
5 passed
```

## C4 Rerun

Output tree:

```text
/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z
```

Command profile: the same C4 one-chunk localized syng-local profile used by the diagnosis pass, with `--debug-dir` enabled for both the pipeline and localized polish stages.

Run result:

- Exit status: 0
- Wall clock: 6:18.09
- Syng-local interval collection: 6,694 sequences, 26,184,511 bp
- Local seed graph build: `path_validation=pass`
- Graph report: 4,488 segments, 2,440 links, 6,694 paths, 866,888 path steps
- Graph report status: `REVIEW` for diagnostic topology metrics only; exact path spelling remained valid
- Debug TSV: `/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z/debug/c4_tuned_pipeline/syng_local_intervals.tsv`
- Path spellings: `/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z/debug/c4_tuned_pipeline/graph_build_0000/path_spellings.tsv`
- Render: `/home/erikg/impg/data/implement_anchored_syng_2_refined_20260609T222143Z/renders_tuned/C4_GRCh38_53kb.png`
- Uploaded PNG: http://hypervolu.me/~erik/impg/implement-anchored-syng-2-C4_GRCh38_53kb-20260609T222143Z.png
- Upload verification: HTTP 200, `Content-Type: image/png`, `Content-Length: 1294973`

Generated GFA/GFA-derived debug/render/log artifacts are under `/home/erikg/impg/data/...` and were not committed.

## Before/After Selected Intervals

The before intervals are from `docs/evaluations/diagnose-c4-query-2.md`.
The after intervals are the emitted fetch intervals from the new debug TSV and match `path_spellings.tsv`.

| Path | Before selected interval | After emitted intervals | Result |
| --- | --- | --- | --- |
| `NA18952#1#JBHIJS010000065.1` | `31855636-32040012`, 184,376 bp | 26 intervals, 122,595 bp: `31887445-31892206`; `31892415-31894486`; `31894564-31896540`; `31896802-31899871`; `31900646-31906064`; `31906311-31908874`; `31909190-31924944`; `31925153-31927224`; `31927302-31929278`; `31929540-31932609`; `31933384-31938802`; `31939049-31941612`; `31941928-31959962`; `31960040-31962016`; `31962278-31965347`; `31966122-31970770`; `31970924-31971540`; `31971787-31972567`; `31972669-31974350`; `31974666-31987151`; `31987220-31992700`; `31992778-31994754`; `31995016-32003508`; `32003662-32005305`; `32005407-32007088`; `32007404-32016729` | Corrected: the single 184 kb long-tail fetch is split into anchor-bounded intervals; terminal tail is not absorbed into one selected span. |
| `NA18948#2#CM101588.1` | `31979894-32164269`, 184,375 bp | 33 intervals, 121,205 bp: `32011703-32015336`; `32015526-32016464`; `32016673-32018744`; `32018822-32020798`; `32021060-32024129`; `32024904-32029552`; `32029706-32030322`; `32030569-32033132`; `32033448-32048074`; `32048264-32049202`; `32049411-32051482`; `32051560-32053536`; `32053798-32056867`; `32057642-32063060`; `32063307-32064087`; `32064189-32065870`; `32066186-32078671`; `32078740-32081940`; `32082149-32084220`; `32084298-32086274`; `32086536-32089605`; `32090380-32095028`; `32095182-32095798`; `32096045-32098608`; `32098924-32113550`; `32113740-32114678`; `32114887-32116958`; `32117036-32119012`; `32119274-32127766`; `32127920-32128536`; `32128783-32129563`; `32129665-32131346`; `32131662-32140986` | Corrected: the single 184 kb long-tail fetch is split into anchor-bounded intervals. |
| `HG02392#1#CM099882.1` | `31884401-32056195`, 171,794 bp | 30 intervals, 106,118 bp: `31916948-31917754`; `31917823-31919895`; `31920216-31923304`; `31923382-31925512`; `31925620-31929896`; `31930345-31934112`; `31934266-31934882`; `31935129-31937545`; `31938008-31951567`; `31951686-31954940`; `31955060-31956267`; `31956716-31958184`; `31958227-31963915`; `31964378-31976863`; `31976932-31977937`; `31978056-31978714`; `31979325-31982413`; `31982491-31984621`; `31984729-31987678`; `31987798-31989005`; `31989454-31993221`; `31993375-31993991`; `31994238-31996654`; `31997117-32007719`; `32009084-32009762`; `32009888-32015778`; `32015922-32019688`; `32019842-32020458`; `32020705-32026764`; `32029073-32032912` | Corrected: the single 171 kb long-tail fetch is split/trimmed by anchor support. |
| `HG01960#2#JBHIHN010000011.1` | One-chunk baseline: `32045898-32158486`, 112,588 bp. Moderate clipped case: `32078636-32158486`, 79,850 bp | 8 intervals, 52,761 bp: `32079320-32081611`; `32081713-32084801`; `32084879-32087009`; `32087117-32091393`; `32091842-32099042`; `32099505-32116754`; `32118213-32125412`; `32125875-32135203` | Not restored to the one-chunk baseline left coordinate. The emitted intervals begin near the clipped coordinate because the retained syng anchors for this path start at `32079440`; the new anchored rule does not reintroduce the unanchored left span. |

## Interpretation

The high-confidence 3-prime long-tail candidates are corrected under the anchored emission rule: they no longer appear as one coordinate-merged fetch interval across repeated C4 blocks. Rows with `split_chain_discontinuous` in the debug TSV mark same-path blocks that are within `-d` but whose query/target anchor gaps are not collinear.

The HG01960-style left clip is not corrected by this patch. The final selected intervals remain close to the clipped start rather than the previous one-chunk baseline start. The debug TSV gives the immediate reason: the left span lacks retained anchor support in the syng chain reaching `HG01960#2#JBHIHN010000011.1`, so an anchor-provenance merge/trim rule cannot safely restore it without falling back to the coordinate-only behavior that caused the long-tail artifacts.

## Validation

- `cargo test test_anchored_emission --lib` passed.
- `cargo test query_scaled_chain_limits_use_merge_distance_before_refinement --lib` passed during implementation.
- `cargo build` passed during implementation.
- Full `cargo test`, final `cargo build`, `cargo install --path .`, `git diff --check`, and the staged artifact/blob scan were run after this report was written; see the WG task log for final command results.
