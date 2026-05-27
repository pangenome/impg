# crush: smoothxg between rounds

User insight (task `crush-smoothxg-between`):

> If mistakes get made at first pass they can stick as we dig into the
> hierarchy.

Before this change, the `crush` pipeline ran N rounds of per-bubble crush and
then a single full-graph `smoothxg` pass at the very end (or no smoothing at
all). When round 1 produced stringy structure inside a replacement subgraph,
round-2 POVU decomposed that stringiness into bogus child bubbles, and the
error propagated downward through the descent tree.

## What changed

`resolve_graph_bubbles` (`src/resolution.rs`) now optionally runs a `smoothxg`
pass on the current working graph **between every round** of per-bubble crush.
Each round's output is cleaned before the next round's POVU decomposition
runs, so stringy artifacts from one round can't seed spurious children in the
next.

### Code layout

- `ResolutionConfig::smooth_between_rounds: Option<SmoothBetweenRoundsConfig>`
  in `src/resolution.rs` controls whether the per-round pass runs and how it
  is parameterised. `None` (the default) preserves prior behaviour exactly.
- Inside the round loop in `resolve_graph_bubbles`, immediately after
  `apply_replacement_frontier` produces `next_graph` and assigns it to
  `graph`, the working graph is rendered to GFA, passed through
  `crate::smooth::smooth_gfa` (the same function the post-crush `:smooth`
  pipeline stage uses), and re-parsed. The resulting `Graph` becomes the
  round-(N+1) discovery input.
- The fresh-id allocator is bumped past the maximum segment id present in the
  smoothed graph so the strict-monotonic id invariant survives the rewrite.

### Why path bp-signatures stay stable

The explicit bubble tree's `bp_signature` is bp-coordinate keyed on path
names; it does not depend on segment ids or graph topology. `smoothxg` is
sequence-preserving, so path bytes (and therefore path bp positions) round
trip unchanged. `strip_full_range_path_names` is called after smoothing so a
bare path name like `chr1` survives the `name → name:0-N → name` round
trip; a query-suffixed name like `chr1:100-200` survives because the
emit-side `:start-end` reconstruction reproduces the same `start-end`.

The library unit test
`resolution::tests::smooth_between_rounds_preserves_path_sequences`
guards this contract on the small insertion-bubble fixture.

## CLI

The per-round pass is opt-in via the `smooth-between` parameter on the
`:crush` pipeline stage:

| Form | Effect |
|---|---|
| `smooth-between=true` | enable with default `target_poa_lengths=[700]` |
| `smooth-between=700` | single-pass at 700 bp |
| `smooth-between=300/700` | two passes (300 bp, then 700 bp) |
| `smooth-between=false` | explicit disable (same as omitting) |

`smooth-between-max-node-length=N` and
`smooth-between-padding-fraction=F` tune the auxiliary smoothxg knobs.

### Aliases

- `smooth-between`, `smooth-between-rounds`, `smoothxg-between` — all equivalent.

### Tuned for speed

Defaults are deliberately lighter than the post-crush `:smooth` stage,
which uses pggb's `[700, 1100]` two-pass invocation: per-round smoothing
runs once per crush round, so a single short pass keeps the per-round
budget below the 2-min target.

## Validation

Real C4 GRCh38 (`GRCh38#0#chr6:31891045-32123783`, 465 paths). Final metrics
and per-round wall times appear in
`data/c4_crush_smoothxg_between_rounds_<timestamp>/`:

- `command.txt` — invocation
- `time.txt` — wall clock and peak RSS
- `stderr.log` — per-round log including new `smoothxg-between` lines
- `run.nosort.gfa` — final graph

PNG: `data/c4_crush_smoothxg_between_rounds_<timestamp>/c4-crush-smoothxg-between-rounds.png`

### Comparison

| Run | Segs | Segment bp | Trivial-stringy | Wall |
|---|---|---|---|---|
| crush+:smooth (`crush-smoothxg-on`, baseline) | 27 k | 410 kb | 191 | 60 min |
| crush+:smooth-between+:smooth (this) | _filled by validation_ | _filled_ | _filled_ | _filled_ |
| PGGB | 13.3 k | 234 kb | 12 | 13 min |

The bottom row is the long-term target; this change is a step toward it by
cleaning the working graph between rounds instead of only at the end.
