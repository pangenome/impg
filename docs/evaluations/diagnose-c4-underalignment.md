# Diagnose C4 Underalignment After POA/abPOA Polish

Task: `diagnose-c4-underalignment`
Date: 2026-06-06

## Verdict

The stage that makes the graph look underaligned relative to the SweepGA seed is
**abPOA motif-local replacement graph construction before lacing**. The accepted
abPOA replacements expand small local windows into near base-level graphs, and
lacing then preserves the paths through that expanded representation.

This is not primarily a Ygs rendering artifact, not primarily missing gfaffix,
not primarily lacing, and not primarily repeated no-op candidate application.
Those other issues are real but explain different symptoms:

- Ygs sorting/rendering changes layout whitespace, not topology or repeat-run
  counts.
- Missing explicit normalization leaves seed-origin self-loop repeat runs
  visible, but POA does not create that problem and normalization removes it
  path-preservingly.
- Lacing is preserving paths; the bad object being laced is the expanded abPOA
  replacement graph.
- Repeated no-op candidates waste late iterations, especially in motif mode, but
  they are secondary to abPOA's first-round segment/path-step explosion.

Non-Ygs fallback PNGs were not used as evidence. For abPOA, the Ygs render
manifests are evidence that Ygs sort failed (`exit 143`) before a PNG could be
made, so the diagnosis uses graph-report topology plus crush logs rather than
fallback images.

## Ygs Metrics

The one-many seed, POA 1kb, and self-loop-normalized POA are all measured on
Ygs-sorted GFAs for the layout-sensitive whitespace columns.

| graph | segments | segment bp | path steps | node coverage p10/median/p90 | bp-weighted coverage | whitespace p99/max bp | direct self-loops | adjacent same-node steps | path preservation |
| --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- |
| SweepGA seed, one-many minmatch1, Ygs | 8,186 | 232,093 | 4,632,725 | 1 / 464 / 931 | 459.684152 | 93 / 37,037 | 152 | 617,670 | baseline |
| POA 1kb, Ygs | 9,079 | 240,589 | 4,843,152 | 1 / 464 / 930 | 443.451172 | 220 / 41,992 | 144 | 604,675 | 465/465 exact |
| POA 1kb self-loop-normalized, Ygs | 9,416 | 243,782 | 4,238,477 | 1 / 464 / 929 | 437.642951 | 98 / 45,061 | 0 | 0 | 465/465 exact |

Sources: `data/c4_self_loop_normalization_20260605T130000Z/seed.Ygs.graph-report.tsv:2`,
`data/c4_self_loop_normalization_20260605T130000Z/poa1kb.Ygs.graph-report.tsv:2`,
`data/c4_self_loop_normalization_20260605T130000Z/normalized.Ygs.graph-report.tsv:2`,
`data/c4_full_poa1kb_20260605T104518Z/compare_gfa_paths.stdout.log:1`,
and `data/c4_self_loop_normalization_20260605T130000Z/compare_gfa_paths.seed-vs-normalized.stdout.log:1`.

This chain rules out POA 1kb and Ygs sorting as the cause of a new
underaligned-looking graph. POA slightly increases segments and segment bp, but
it slightly reduces direct self-loops and adjacent same-node steps. The explicit
self-loop normalizer then removes direct self-loop edges and adjacent same-node
path steps while preserving all 465 paths.

## Motif POA vs abPOA

| graph | Ygs render status | segments | segment bp | path steps | node coverage p10/median/p90 | bp-weighted coverage | whitespace p99/max bp | direct self-loops | adjacent same-node steps | path preservation |
| --- | --- | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | --- |
| k311 SweepGA seed, Ygs | Ygs graph-report available | 7,151 | 250,580 | 3,129,613 | 1 / 462 / 930 | 425.770109 | 234 / 36,303 | 60 | 263,737 | baseline |
| k311 POA motif25k, Ygs | ok | 10,201 | 298,971 | 3,460,617 | 1 / 452 / 466 | 356.855595 | 366 / 45,117 | 24 | 40,952 | 465/465 exact |
| k311 abPOA motif10k from seed | Ygs sort failed, no Ygs PNG | 178,142 | 280,038 | 80,156,435 | 1 / 454 / 466 | 380.982131 | Ygs unavailable | 24 | 40,209 | 465/465 exact |
| k311 abPOA motif10k from POA2kb input | Ygs sort failed, no Ygs PNG | 177,970 | 283,273 | 79,717,063 | 1 / 452 / 466 | 376.631285 | Ygs unavailable | 27 | 41,956 | 465/465 exact |

Sources: `data/c4_k311_poa_threshold_20260605T140018Z/reports/c4.k311.seed.Ygs.graph-report.tsv:2`,
`data/c4_aggressive_motif_matrix_20260605T210000Z/reports/c4.k311.poa.motif25k.cand512.iter8.Ygs.graph-report.tsv:2`,
`data/c4_sweepga_seed_abpoa10k_20260605T223500Z/reports/c4.k311.seed.abpoa.motif10k.cand512.iter8.graph-report.tsv:2`,
`data/c4_aggressive_motif_matrix_20260605T210000Z/reports/c4.k311.abpoa.motif10k.cand512.iter8.graph-report.tsv:2`,
`data/c4_aggressive_motif_matrix_20260605T210000Z/reports/rendered-best.tsv:2`,
`data/c4_aggressive_motif_matrix_20260605T210000Z/reports/rendered-best.tsv:3`,
`data/c4_sweepga_seed_abpoa10k_20260605T223500Z/reports/rendered.tsv:2`,
`data/c4_aggressive_motif_matrix_20260605T210000Z/reports/c4.k311.poa.motif25k.cand512.iter8.compare_gfa_paths.stdout.txt:1`,
`data/c4_sweepga_seed_abpoa10k_20260605T223500Z/reports/c4.k311.seed.abpoa.motif10k.cand512.iter8.compare_gfa_paths.stdout.txt:1`,
and `data/c4_aggressive_motif_matrix_20260605T210000Z/reports/c4.k311.abpoa.motif10k.cand512.iter8.compare_gfa_paths.stdout.txt:1`.

The contrast is decisive: POA motif-local changes the graph moderately and
reduces repeat-run artifacts; abPOA motif-local expands segment count by about
25x and path steps by about 25.6x while preserving path spellings. That shape is
what makes the output appear underaligned.

## Stage Attribution

Replacement graph construction is the first bad stage. In the abPOA seed run,
round 1 built a candidate from 5 input segments / 305 bp into 330 output
segments / 330 bp, and another from 4 input segments / 238 bp into 252 output
segments / 252 bp. Both were then applied in the same round. See
`data/c4_sweepga_seed_abpoa10k_20260605T223500Z/logs/c4.k311.seed.abpoa.motif10k.cand512.iter8.crush.stderr.log:38`,
`:121`, `:237`, and `:245`.

Lacing is not the primary fault. The compare outputs report 465 expected paths,
465 observed paths, and zero spelling mismatches for POA, normalized POA, and
abPOA. That means lacing is integrating the replacements without corrupting the
path spellings; it is not deciding whether the replacement representation is
too fragmented.

Missing gfaffix/normalization is not the primary underalignment cause. Direct
self-loops are already present in the SweepGA seed, POA reduces them slightly in
the 1kb chain, and explicit normalization removes them exactly. The normalizer
is therefore the right fix for residual self-loop repeat artifacts, but it does
not explain abPOA's 178k-segment graph.

Repeated no-op candidate application is present but secondary. Late POA motif
rounds apply many sparse-offshoot candidates with `objective_score_delta=0` and
zero segment deltas, for example
`data/c4_aggressive_motif_matrix_20260605T210000Z/logs/c4.k311.poa.motif25k.cand512.iter8.crush.stderr.log:556`
through `:571`. The same pattern appears in late abPOA rounds. This wastes
iterations and should be fixed, but it does not account for the initial abPOA
explosion already visible in round 1.

## Recommended Next Change

Harden candidate acceptance before lacing:

1. Normalize/compact every replacement graph locally before scoring.
2. Reject any replacement whose local graph-report deltas expand segments or
   path steps beyond the source window unless it also has a positive objective
   improvement large enough to justify the expansion.
3. Make `objective_score_delta <= 0` candidates non-applicable rather than
   diagnostic-only, so repeated no-op sparse-offshoot candidates are skipped
   before lacing.

For the residual loop artifact specifically, keep the explicit self-loop
normalizer as a final path-preserving postprocess after POA/Ygs output. Do not
expect gfaffix or Ygs sorting to solve that topology issue.

## Validation

- Used Ygs-sorted graph-report metrics for seed, POA 1kb, normalized POA, k311
  seed, and POA motif25k wherever Ygs sort completed.
- Treated abPOA Ygs `gfasort-failed` rows as failed render evidence and did not
  use fallback Yg/non-Ygs PNGs.
- Checked graph-report requirements: segment count, total segment bp, path
  steps, node coverage distribution, whitespace distribution, direct
  self-loops, adjacent same-node path steps, and path preservation.
- Identified the exact failing stage as abPOA replacement graph construction,
  with lacing, missing normalization/gfaffix, and no-op application separated as
  secondary/non-primary causes.
