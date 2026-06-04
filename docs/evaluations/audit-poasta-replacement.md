# Audit: POASTA Replacement Cycles In C4 Crush

Task: `audit-poasta-replacement`

Branch/worktree run: `wg/agent-478/audit-poasta-replacement`

Diagnostic run directory:
`data/audit_poasta_replacement_20260604T182011Z`

Input seed:
`/home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/sorted/combined_compact_12000bp_4sites_top1_r1.Ygs.gfa`

## Instrumentation

The scoped debug path is enabled by `IMPG_CRUSH_DEBUG_DIR`.

- `src/resolution.rs:779` stores replacement method/source/objective metadata on `ReplacementPlan`; multi-level candidates attach it at `src/resolution.rs:3082`.
- `src/resolution.rs:9027` writes raw POASTA `graph_to_gfa` output to `replacement_NNNN_poasta/raw.graph_to_gfa.gfa`.
- `src/resolution.rs:9034` writes the exact normalized graph after `poasta_gfa_to_exact_graph` and path ordering to `replacement_NNNN_poasta/exact.normalized.gfa`.
- `src/resolution.rs:7713` still performs exact path-sequence validation before writing the final debug frontier.
- `src/resolution.rs:10231` writes `applied_frontier_0000/final_laced.gfa`, `applied-candidates.tsv`, and one `replacement.gfa` plus `path-lengths.tsv` per applied candidate.
- `scripts/audit_poasta_replacement_cycles.py:237` computes directed SCC/self-loop summaries; `scripts/audit_poasta_replacement_cycles.py:285` compares exact path sequences; `scripts/audit_poasta_replacement_cycles.py:339` analyzes the final laced graph, applied replacements, and raw/exact POASTA GFAs; `scripts/audit_poasta_replacement_cycles.py:512` summarizes final-cycle node origins by inferred replacement ID ranges.

The cap64 rerun used the same key crush parameters as the prior diagnostic run:

```bash
IMPG_CRUSH_DEBUG_DIR=data/audit_poasta_replacement_20260604T182011Z/debug \
target/release/impg crush \
  --gfa /home/erikg/impg/data/c4_grouped_multi_bubble_seed_top1_20260604T121002Z/sorted/combined_compact_12000bp_4sites_top1_r1.Ygs.gfa \
  --output data/audit_poasta_replacement_20260604T182011Z/graphs/multilevel_cap64_spoa2k_poasta15k_r1.gfa \
  --method iterative-multi-level \
  --window-mode combined \
  --window-target-bp 12000 \
  --max-window-sites 4 \
  --candidate-limit 64 \
  --max-iterations 1 \
  --max-span 0 \
  --max-traversal-len 100k \
  --max-median-traversal-len 50k \
  --max-total-sequence 500m \
  --max-traversals 100k \
  --auto-spoa-max-traversal-len 2k \
  --auto-poasta-max-traversal-len 15k \
  --poa-scoring 1,4,6,2,26,1 \
  --min-match-length off \
  --max-pair-alignments 0 \
  --max-replacement-paf-bytes 0 \
  --sweepga-no-filter true \
  --polish-rounds 0 \
  --threads 32 \
  --verbose 1
```

The run completed successfully in 4:17.81 wall-clock seconds according to
`data/audit_poasta_replacement_20260604T182011Z/logs/multilevel_cap64_spoa2k_poasta15k_r1.stderr.log`.

## Artifacts

- `data/audit_poasta_replacement_20260604T182011Z/debug/applied_frontier_0000/applied-candidates.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/debug/applied_frontier_0000/final_laced.gfa`
- `data/audit_poasta_replacement_20260604T182011Z/analysis/gfa-cycle-summary.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/analysis/path-preservation.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/analysis/candidate-path-preservation.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/analysis/poasta-build-matches.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/analysis/final-cycle-origins.tsv`
- `data/audit_poasta_replacement_20260604T182011Z/renderings/Poasta_crush02_path-labels.svg`
- `data/audit_poasta_replacement_20260604T182011Z/renderings/Poasta_crush04_path-labels.svg`

The two renderings were produced with `gfalook` without `--hide-path-names`.
The SVGs contain path-name text labels, e.g. `CHM13#0#chr6...` in
`Poasta_crush02_path-labels.svg`.

## Exact Path Preservation

No quality metric was used as an acceptance gate. The only hard gate was exact
path-sequence preservation.

Final graph preservation, from `analysis/path-preservation.tsv`:

| Observed GFA | Ref paths | Observed paths | Common | Missing | Extra | Mismatched |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `graphs/multilevel_cap64_spoa2k_poasta15k_r1.gfa` | 465 | 465 | 465 | 0 | 0 | 0 |
| `debug/applied_frontier_0000/final_laced.gfa` | 465 | 465 | 465 | 0 | 0 | 0 |

Focused replacement preservation, from
`analysis/candidate-path-preservation.tsv`:

| Replacement | Paths checked | Failed |
| --- | ---: | ---: |
| `Poasta_crush02_chw_span14980bp_med14980bp_cov465of465_sites152_steps2555-2864` | 465 | 0 |
| `Poasta_crush04_top_span6460bp_med6461bp_cov47of465_sites1_steps3320-3523` | 47 | 0 |
| `Poasta_crush22_chw_span14944bp_med14949bp_cov465of465_sites154_steps7377-7833` | 465 | 0 |
| `Poasta_crush23_chw_span14934bp_med14934bp_cov465of465_sites163_steps6059-6650` | 465 | 0 |

## Focused Cycle Counts

Counts from `analysis/gfa-cycle-summary.tsv`.

| GFA | Segments | Links | Paths | Path steps | Directed SCCs >1 | Directed self-loop edges |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Poasta_crush02.../replacement.gfa` | 447 | 603 | 465 | 140835 | 0 | 0 |
| `Poasta_crush04.../replacement.gfa` | 72 | 96 | 47 | 1900 | 0 | 0 |
| `Poasta_crush22.../replacement.gfa` | 497 | 694 | 465 | 153406 | 0 | 0 |
| `Poasta_crush23.../replacement.gfa` | 607 | 818 | 465 | 189277 | 0 | 0 |

The focused POASTA replacements are directed acyclic under the oriented-link
analysis. There are also no segment-level self-loop links for those four rows.

## POASTA Raw/Exact Match

The focused applied replacements matched these raw/exact POASTA debug dirs by
path-sequence digest, from `analysis/poasta-build-matches.tsv`:

| Applied replacement | Build dir | Raw SCCs | Exact SCCs | Applied SCCs | Raw self-loops | Exact self-loops | Applied self-loops |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| `Poasta_crush02...` | `replacement_0002_poasta` | 0 | 0 | 0 | 0 | 0 | 0 |
| `Poasta_crush04...` | `replacement_0012_poasta` | 0 | 0 | 0 | 0 | 0 | 0 |
| `Poasta_crush22...` | `replacement_0020_poasta` | 0 | 0 | 0 | 0 | 0 | 0 |
| `Poasta_crush23...` | `replacement_0017_poasta` | 0 | 0 | 0 | 0 | 0 | 0 |

Across the whole debug set:

| Kind | GFAs checked | Cyclic GFAs |
| --- | ---: | ---: |
| Applied replacements | 23 | 2 |
| Raw POASTA `graph_to_gfa` GFAs | 21 | 0 |
| Exact normalized POASTA GFAs | 21 | 0 |

The only cyclic applied replacements are `Sweepga_crush01` and
`Sweepga_crush03`. They account for the two cyclic applied-replacement rows in
`analysis/gfa-cycle-summary.tsv`.

## Final Laced Graph

The final laced graph is cyclic:

| Graph | Segments | Links | Paths | Path steps | Directed SCCs >1 | Directed self-loop edges | Largest SCC | Cyclic oriented nodes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Input seed | 7340 | 10015 | 465 | 3765901 | 72 | 152 | 2020 | 5596 |
| Final laced/debug | 7405 | 10122 | 465 | 3619762 | 60 | 136 | 2090 | 4924 |
| Final output | 7405 | 10122 | 465 | 3619762 | 60 | 136 | 2090 | 4924 |

The seed graph was already cyclic. The final output exactly matches the debug
`final_laced.gfa` counts, so these cycles exist before any external sort,
gfaffix, or post-processing.

`analysis/final-cycle-origins.tsv` infers replacement ID ranges from the input
max segment ID and applied-candidate order. It shows the final SCC/self-loop
origins:

| Cycle type | Origin class | Count | Max size | Notes |
| --- | --- | ---: | ---: | --- |
| Directed self-loop | `input_seed` | 104 | 1 | Inherited seed loops remain. |
| Directed self-loop | `Sweepga` | 24 | 1 | From `Sweepga_crush01`. |
| Directed self-loop | `Poasta` | 8 | 1 | All in `Poasta_crush22` after final lacing. |
| Directed SCC | `input_seed` | 40 | 222 | Seed-only cycles remain. |
| Directed SCC | `Sweepga` | 10 | 9 | From `Sweepga_crush01`. |
| Directed SCC | `Sweepga,input_seed` | 2 | 2090 | Large SCC crosses `Sweepga_crush03` and seed nodes. |
| Directed SCC | `Poasta` | 2 | 4 | `Poasta_crush22` nodes after final lacing. |
| Directed SCC | `Poasta` | 4 | 8 | Crosses `Poasta_crush22` and `Poasta_crush23` nodes after final lacing. |
| Directed SCC | `Poasta,Sweepga` | 2 | 31 | Crosses `Poasta_crush02` and `Sweepga_crush01` after final lacing. |

This means final cycles can include POASTA replacement node IDs, but those
cycles are not present in the POASTA replacement GFAs themselves.

## Conclusion

For this cap64 C4 run, POASTA itself did not emit directed cycles or self-loops
in the replacement GFAs. The evidence separates the stages:

- POASTA `graph_to_gfa`: all 21 raw POASTA debug GFAs are acyclic.
- W-walk clipping/slice interning and exact normalization: all 21 exact
  normalized POASTA GFAs are acyclic.
- Applied POASTA replacements: the focused replacements and every applied
  POASTA replacement are acyclic and preserve candidate path sequences exactly.
- Final lacing: the final graph is cyclic and some final SCCs/self-loops include
  POASTA replacement node IDs only after lacing.
- Seed graph: the input seed is already cyclic, and many final cycles are
  inherited seed-only cycles.
- External sort/gfaffix/unchop: not implicated by this run. The debug
  `final_laced.gfa` and direct final output have identical counts and cycle
  summaries before any external post-processing.

The apparent POASTA loops are therefore final-lacing/seed/neighbor-replacement
artifacts, not POASTA `graph_to_gfa` or POASTA normalization artifacts. The most
specific stage where POASTA node IDs become cyclic is the final
`render_rewritten_graph`/path-link synthesis step that laces all accepted
replacements back into full paths while preserving exact path sequences.
