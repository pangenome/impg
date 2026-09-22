# Local Compression Testbed Fast Synthesis

Date: 2026-06-07

This synthesis uses the actual `follow-up-produce` fast-profile artifacts:

- `tests/test_data/local_compression/manifest.json`
- `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv`
- `docs/evaluations/local-compression-testbed-fast/scoreboard.json`
- `docs/evaluations/local-compression-testbed-fast/report.md`
- `docs/evaluations/local-compression-testbed-fast/validation.md`
- `docs/evaluations/local-compression-testbed-fast/fixtures/*`

The fast run wrote 143 scoreboard rows: 13 fixtures by 11 methods. Of those,
80 rows produced graphs and 63 rows were explicit skips. There were no command
failures and no hard path-corruption rows. The run is useful as a first
end-to-end harness check, but not yet as a real algorithm bakeoff.

## Headline Answer

We learned that the fixture runner can now exercise a row-complete fast profile,
record exact path preservation, record topology assertions, and explain profile
skips. On every graph-producing row, exact paths were preserved: 80/80 produced
graphs had `exact_path_preservation=pass` and `hard_path_corruption=false`.

The topology result is also clear but shallow. `local_syng_raw` preserved paths
while failing all 10 included CI fixture topology assertions because it emits
raw path-copy graphs with no induced links, bubbles, or flubbles. Every other
local compact/window method passed topology on every included CI fixture:
70/70 pass across the seven compacting/windowed methods.

The caveat is important: those seven successful method IDs are currently routed
through the same fixture-local `compact_bubble` generator. Their per-method logs
say `generated compact_bubble graph`, and representative metrics for POASTA,
SmoothXG-style chunking, highest-level flubble, local SweepGA, and whole-region
SweepGA rows are identical on the same fixture. So the run supports the path
guard and expected-topology testbed, but it does not yet prove that any specific
resolver or windowing policy is better.

## Method Counts

Counts below are from `scoreboard.json`. `Command skipped` includes both
fast-profile local-fixture exclusions and optional-control exclusions.

| Method | Command pass | Command skipped | Topology pass | Topology fail | Topology not_run | Exact paths pass | Exact paths not_run | Candidates | Candidate skips |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `local_syng_raw` | 10 | 3 | 0 | 10 | 3 | 10 | 3 | 0 | 0 |
| `local_syng_crush_auto` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `local_syng_crush_poa` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `local_syng_crush_poasta` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `local_syng_crush_sweepga` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `top_flubble_nonoverlap_sweepga` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `chunk_window_smooth_or_crush` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `whole_region_sweepga_seqwish` | 10 | 3 | 10 | 0 | 3 | 10 | 3 | 10 | 0 |
| `pggb_control` | 0 | 13 | 0 | 0 | 13 | 0 | 13 | 0 | 0 |
| `smoothxg_control` | 0 | 13 | 0 | 0 | 13 | 0 | 13 | 0 | 0 |
| `pggb_plus_smoothxg_control` | 0 | 13 | 0 | 0 | 13 | 0 | 13 | 0 | 0 |

Interpretation:

- `local_syng_raw` is path-safe but topology-incomplete by design in this run.
- All seven local compact/window method IDs are path-safe and topology-passing
  on included CI fixtures.
- The three external control methods did not run in the fast profile, so they
  provide no path or topology evidence.

## Fixture-Class Counts

Each manifest fixture class has one fixture in this first testbed. The fast
profile included the 10 `ci` fixtures and skipped the three `local` fixtures.

| Fixture class | Tier | Fixture | Command pass | Command skipped | Topology pass | Topology fail | Topology not_run | Exact paths pass | Exact paths not_run | Raw topology | Compact/window topology | Controls |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- | --- | --- |
| `snp_bubble` | ci | `snp_bubble_3path` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `short_indel` | ci | `short_indel_3path` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `insertion_50_500bp` | ci | `mid_insertion_200bp` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `alu_like_insertion` | local | `alu_like_insertion` | 0 | 11 | 0 | 0 | 11 | 0 | 11 | not_run | 7/7 not_run | 3 `profile_excludes_local_fixture` |
| `adjacent_bubbles_compress_together` | ci | `adjacent_bubbles_joint` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `bubble_split_by_fake_repeat_anchor` | ci | `fake_repeat_anchor_split` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `repeated_motif_microtangle` | ci | `microtangle_repeat_motif` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `duplicated_flank_requires_path_context` | ci | `duplicated_flank_context` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `tandem_copy_number_loop_cyclic` | ci | `tandem_copy_loop_keep` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `dispersed_repeat_glue_break_or_ignore` | local | `dispersed_repeat_glue_break` | 0 | 11 | 0 | 0 | 11 | 0 | 11 | not_run | 7/7 not_run | 3 `profile_excludes_local_fixture` |
| `inversion_like` | ci | `inversion_like_case` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `nested_bubbles_top_level_right` | ci | `nested_top_level_right` | 8 | 3 | 7 | 1 | 3 | 8 | 3 | fail | 7/7 pass | 3 `profile_excludes_optional` |
| `nested_bubbles_top_level_wrong` | local | `nested_top_level_wrong` | 0 | 11 | 0 | 0 | 11 | 0 | 11 | not_run | 7/7 not_run | 3 `profile_excludes_local_fixture` |

## Topology Assertions

The compact/window methods met every expected topology assertion on the 10
included CI fixtures. Raw SYNG did not meet any included topology assertion.

| Fixture | Assertion | `local_syng_raw` | Seven compact/window methods |
| --- | --- | --- | --- |
| `snp_bubble_3path` | `single_ordered_snp_bubble` | fail | 7/7 pass |
| `short_indel_3path` | `short_indel_bubble` | fail | 7/7 pass |
| `mid_insertion_200bp` | `mid_insertion_shared_flanks` | fail | 7/7 pass |
| `adjacent_bubbles_joint` | `adjacent_bubbles_joint_context` | fail | 7/7 pass |
| `fake_repeat_anchor_split` | `fake_repeat_anchor_not_boundary` | fail | 7/7 pass |
| `microtangle_repeat_motif` | `bounded_repeat_motif_microtangle` | fail | 7/7 pass |
| `duplicated_flank_context` | `duplicated_flank_context_required` | fail | 7/7 pass |
| `tandem_copy_loop_keep` | `tandem_copy_loop_required` | fail | 7/7 pass |
| `inversion_like_case` | `inversion_like_orientation_preserved` | fail | 7/7 pass |
| `nested_top_level_right` | `nested_parent_boundary_right` | fail | 7/7 pass |

The three local-only assertions were present in the rows but not run under the
fast profile:

| Fixture | Assertion | Fast-profile status |
| --- | --- | --- |
| `alu_like_insertion` | `alu_like_repeat_insertion` | skipped: `profile_excludes_local_fixture` |
| `dispersed_repeat_glue_break` | `dispersed_repeat_no_long_glue` | skipped: `profile_excludes_local_fixture` |
| `nested_top_level_wrong` | `nested_parent_boundary_wrong` | skipped: `profile_excludes_local_fixture` |

## Path Preservation

No method corrupted exact paths in this run.

- Preserved exact paths on produced rows: `local_syng_raw`,
  `local_syng_crush_auto`, `local_syng_crush_poa`,
  `local_syng_crush_poasta`, `local_syng_crush_sweepga`,
  `top_flubble_nonoverlap_sweepga`, `chunk_window_smooth_or_crush`, and
  `whole_region_sweepga_seqwish`.
- Not run: `pggb_control`, `smoothxg_control`, and
  `pggb_plus_smoothxg_control`.
- Not run on local-only fixtures in the fast profile: all methods.
- Hard path corruption: 0/143 rows.
- Missing expected path count on produced rows: 0.
- Extra path count on produced rows: 0.

This is the strongest result from the first pass: the testbed now catches path
status row-by-row, and every produced local graph stayed path-exact.

## Failure Localization

The scoreboard does not show candidate-window failures in this run.

- Successful compact/window rows report `candidate_count=1`,
  `candidate_skipped_count=0`, and an empty `candidate_skip_reasons` field.
- `local_syng_raw` reports `candidate_count=0` because its method parameters
  say `candidate_windows=none`; its topology failures are raw-output failures,
  not candidate-discovery failures.
- External controls are skipped by profile policy, not by missing mandatory
  tools or resolver failures.
- Local-only fixtures are skipped by profile policy.

The scoreboard also does not expose real resolver or alignment failures yet.
All successful compact/window method IDs use the same `compact_bubble` strategy
in the runner, and their stdout logs report `generated compact_bubble graph`.
That makes POA, POASTA, local SweepGA, highest-level flubble windows,
SmoothXG-style chunks, and whole-region SweepGA indistinguishable on these
rows. There are no command failures, no nonempty candidate skip reasons, and no
resolver-specific stderr evidence to classify.

## Which Direction Looks Most Promising

On these tiny fixtures, the most defensible conclusion is that the SYNG-local
crush scaffold is the promising path, not that a specific resolver has won.
The compact/window rows show that a local replacement can be path-exact and
can satisfy fixture topology assertions across SNPs, short indels, mid-size
insertions, adjacent bubbles, fake internal repeats, repeated motifs,
duplicated flank context, tandem-copy loops, inversion-like cases, and the
top-level-right nested-bubble case.

For the specific candidate directions:

- SmoothXG-style chunking: passes all included topology checks, but only through
  the same compact-bubble placeholder as the other compact methods.
- Highest-level flubbles: passes all included topology checks, but the first
  run does not prove the highest-level window policy itself.
- Whole-region SweepGA/seqwish: passes all included topology checks, but the
  runner did not run real whole-region SweepGA/seqwish induction here.
- SYNG-local crush: best-supported direction because it is the native path
  being exercised by the runner and path guards, but resolver evidence is still
  pending.

So the evidence for a path toward fast-like-SYNG is positive but preliminary:
we have a path-exact local scaffold with topology tests and explicit row-level
diagnostics. The evidence for compressed-like-PGGB is not yet present: PGGB and
SmoothXG controls were skipped, and the current local compact rows do not yet
measure PGGB-like compression quality beyond coarse topology counters.

## Limitations

- The fixtures are synthetic and tiny. They are good for invariants, but they
  do not yet stress scale, repeats, ambiguous anchors, runtime, or memory.
- Three important fixtures were local-only and skipped in the fast profile:
  `alu_like_insertion`, `dispersed_repeat_glue_break`, and
  `nested_top_level_wrong`.
- PGGB, SmoothXG, and PGGB+SmoothXG controls were optional and skipped under
  the fast profile.
- Fast-profile rendering was skipped with
  `render_tool_not_configured_for_fast_profile`, so there is no visual audit.
- The runner still uses placeholders for the successful compact/window method
  IDs. The scoreboard distinguishes method IDs and parameters, but the graph
  generation path collapses them to the same `compact_bubble` output.
- Candidate-window diagnostics are too shallow for failure triage: this run has
  only one candidate per successful compact/window row and zero candidate skips.
- Metrics are not strong enough to claim compressed-like-PGGB behavior. Current
  metrics include graph size, segment/link/path counts, depth proxies, bubble
  and flubble counts, long-link counts, and whitespace proxies, but no robust
  compression-quality objective, alignment/replay audit, or scale target.

## Next Concrete Work Item

Replace the compact-bubble placeholder in the local-compression runner with a
real candidate-window and resolver execution path for at least two adversarial
CI fixtures, then rerun the same scoreboard.

The next implementation should:

1. Discover candidate windows independently from expected fixture metadata and
   record their coordinates, boundary rationale, and skip reasons.
2. Dispatch method-specific resolvers instead of the shared `compact_bubble`
   generator for `local_syng_crush_poa`, `local_syng_crush_poasta`,
   `local_syng_crush_sweepga`, `top_flubble_nonoverlap_sweepga`,
   `chunk_window_smooth_or_crush`, and `whole_region_sweepga_seqwish`.
3. Record resolver/alignment status separately from candidate-window status:
   accepted, rejected, no candidate, alignment failed, path replay failed, or
   topology assertion failed.
4. Promote one currently local-only adversarial fixture, preferably a
   CI-sized variant of `nested_top_level_wrong` or
   `dispersed_repeat_glue_break`, into the fast profile so candidate-boundary
   mistakes are visible in the default run.
5. Keep exact path preservation as a hard gate and add at least one compression
   quality metric that can separate "topology-shaped placeholder" from a real
   compressed graph.

That work would turn this testbed from a harness validation into the first
meaningful comparison of fast-like-SYNG local crush against compressed-like
PGGB/SmoothXG behavior.
