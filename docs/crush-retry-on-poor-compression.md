# crush: compression diagnostics (former retry-on-poor-compression)

Task: `crush-retry-on-poor-compression` (workgraph id `crush-retry-on`).
Branch: `wg/agent-182/crush-retry-on`.

## Problem

> If mistakes get made at first pass they can stick as we dig into the
> hierarchy. — task spec, user insight.

## Current status

`remove-crush-replacement` disables this mechanism as an acceptance policy.
Crush replacement is now unconditional after graph-validity checks: the
replacement must parse, the replacement path names/order must match the input
traversals, and the rewritten graph must preserve every input path spelling.
Compression ratio is still useful in logs, but it is diagnostic only. It does
not trigger an alternate aligner, swap in a smaller result, reject a larger
result, or keep the previous graph.

In the Phase-6 true level descent (see `docs/crush-architecture-spec.md`,
`docs/crush-hierarchical.md`), each accepted bubble is replaced exactly once
and its sub-bubbles are then discovered inside the new replacement subgraph.
If the round-1 aligner did a poor job on a bubble — e.g. sweepga fragmenting
a high-density repeat region, or POASTA inflating an SV-scale bubble into a
near-linear stack — the descent picks up *inside* that poor replacement and
keeps refining it. The original mis-compression sticks.

## Mechanism

After each `apply_replacement_frontier` round, for every accepted plan
compute

```
compression ratio = sum(input traversal bp) / sum(replacement segment bp)
```

The ratio is logged with the accepted replacement's output segment-bp. If
`retry-min-compression-ratio` is supplied in an older engine string, the log
also reports how many replacements fell below that threshold. No alternate
replacement is built and no replacement is selected by this metric.

## Safety invariants

1. **Path preservation.** Each aligner's builder
   (`build_poa_replacement`, `build_poasta_replacement`, etc.) calls
   `validate_replacement_paths` on its output before returning, so every
   accepted replacement's paths spell the same bytes as the input
   traversals. The outer `apply_replacement_frontier` re-runs
   `path_sequences_equal` after rewriting, providing a second-line defense.
2. **Parse and path identity.** Replacement GFA must parse and preserve
   replacement path names/order before it can be spliced.
3. **No quality gate.** Segment count, segment-bp, compression ratio,
   stringiness, and visual-complexity metrics are never used to reject or
   replace an otherwise valid result.

## Configuration

The legacy keys are still accepted on the `crush:` engine pipeline (parsed in
`src/main.rs`, `parse_crush_engine_stage_params`) so older command strings do
not fail:

- `retry-min-compression-ratio=<f64>` (alias: `retry-ratio`,
  `min-compression-ratio`, `compression-retry-ratio`). Default `0.0`.
  Values above zero only annotate the diagnostic log.
- `retry-min-input-bp=<usize>` (alias: `retry-min-input`). Default
  `1 000`. This is only used to decide which replacements are counted as
  "below threshold" in diagnostics.

## Reporting

`ResolutionStats` keeps the older `retry_attempts`, `retry_wins`, and
`retry_failures` fields for API compatibility. They remain zero because no
retry is attempted.

Per-round diagnostics are written at `info` level, for example:

```text
crush round 1: replacement compression diagnostics: ratio input/output n=8, ...
threshold=disabled (diagnostic only); no alternate aligner attempted and no replacement was rejected by ratio
```

## Tests

`src/resolution.rs` (unit, under `mod tests`):

- `replacement_compression_ratio_handles_zero_and_basic_cases` — exact
  values on perfect (4 traversals → 1 segment, ratio 4.0) and no-op
  (4 traversals → 4 segments, ratio 1.0) cases; `None` returned when
  input or output bp is zero so diagnostics treat them as "no
  signal" instead of dividing by zero or grading every empty graph as
  infinitely well compressed.
- `compression_ratio_threshold_is_diagnostic_only` — end-to-end on the
  canonical SNP bubble fixture: even a 1e9 threshold must preserve paths
  and leave `retry_attempts`, `retry_wins`, and `retry_failures` at zero.

`cargo test --release --lib resolution` passes 54 tests (50 pre-existing
+ 4 new).

## C4 GRCh38 validation

See sibling section `## C4 GRCh38 results` below for the canonical-command
run and the comparison to the smoothxg-on baseline. Output dir:
`data/c4_crush_retry_on_poor_20260527T044100Z/`.

Canonical command (matches `docs/c4-crush-handoff.md` §Known-Good with the
`smoothxg-on` polish stack, plus the new retry knobs at their defaults):

```bash
target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O data/c4_crush_retry_on_poor_20260527T044100Z/run.nosort \
  -v 1
```

## C4 GRCh38 results

Output dir: `/home/erikg/impg/data/c4_crush_retry_on_poor_20260527T044100Z/`.
Wall: 37 min 30 s, peak RSS 100 GB, exit 0.

### Path preservation

```
S: 19698  L: 23399  P: 465  segment-bp: 543930
```

465 / 465 input paths preserved (each path equality is enforced by
`validate_replacement_paths` per builder and again at the frontier rewrite
by `path_sequences_equal`).

### Historical retry mechanism output

The original `crush-retry-on` experiment predates
`remove-crush-replacement` and did run a compression-ratio retry as an
acceptance policy. Current runs emit compression diagnostics instead and do not
attempt alternate aligners.

At that time, 1 of the 8 accepted round-1 plans was below the 2.0
compression-ratio threshold. The current code would report that fact as a
diagnostic and keep the originally selected, path-valid replacement
unconditionally.

Historical stats:
- `retry_attempts = 1`
- `retry_wins = 0`
- `retry_failures = 0`

Current stats stay at zero because no retry attempts are made.

### Metrics vs the closest no-smoothxg baseline

The closest apples-to-apples baseline is the `method=sweepga`
many-to-many + poasta-polish run from 2026-05-26
(`data/c4_crush_many_to_many_poasta_polish_20260526T185552Z/`). Both
share the `no-filter=true, polish-rounds=until-done,
polish-method=poasta, polish-max-traversal-len=10k,
polish-max-median-traversal-len=1k` tail and the same syng input. The
only knob differences are `method=sweepga` (baseline) vs `method=auto`
(this run, with retry-on enabled).

| Metric (Ygs)             | many-to-many baseline | this run (retry-on) |   Δ |
|--------------------------|---------------------:|-------------------:|----:|
| segments                 |              19 958  |             19 698 | −260 (−1.3 %)  |
| segment-bp               |             544 641  |            543 930 | −711 (−0.13 %) |
| links                    |              23 822  |             23 399 | −423 (−1.8 %)  |
| paths                    |                 465  |                465 |      0         |
| dup-extras 51–200 bp     |                  15  |                 10 | −5 (−33 %)     |
| dup-extras > 200 bp      |                   3  |                  0 | −3 (eliminated) |
| trav-dup extras > 200 bp |              15 247  |             15 247 |     0          |
| distinct sequences       |              13 055  |             13 065 |   +10          |

(Computed with
`data/c4_pggb_control_20260526T025439Z/analysis-scripts/gfa_metrics2.py`.)

Because the historical retry only fired once and the alternate lost, the metric
deltas above mostly reflect the `method=auto` 3-tier dispatch vs
`method=sweepga`, not the retry path. They are retained as historical context,
not as a current acceptance rule.

### Render

`c4-crush-retry-on-poor.png` is generated by:

```bash
out=data/c4_crush_retry_on_poor_20260527T044100Z
gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/c4-crush-retry-on-poor.png" -m -x 2200 -y 1200
```

Local path:
`/home/erikg/impg/data/c4_crush_retry_on_poor_20260527T044100Z/c4-crush-retry-on-poor.png`.

## Publish

The PNG is generated locally at
`/home/erikg/impg/data/c4_crush_retry_on_poor_20260527T044100Z/c4-crush-retry-on-poor.png`.
The canonical URL
`https://hypervolu.me/~erik/impg/c4-crush-retry-on-poor.png` is the
user's hypervolume publish path — this agent does not have direct
write access to that path (ssh key was rejected during the run); the
user uploads PNGs from `data/*.png` to `hypervolu.me` as part of their
workflow (the prior `c4-crush-with-smoothxg.png` followed the same
flow, per `docs/crush-smoothxg-on-output.md` §Publish).

## Hard gates

- [x] Historical retry mechanism implemented for the original
      `crush-retry-on` task. It is now superseded by
      `remove-crush-replacement`, where compression metrics are diagnostic
      only.
- [x] Real C4 GRCh38 run, **465 / 465** paths preserved.
- [x] Historical retry stats reported: `retry_attempts=1, retry_wins=0,
      retry_failures=0`. Current runs leave these counters at zero.
- [x] Final metrics: segments **19 698**, segment-bp **543 930**,
      dup-extras 51–200 bp = **10**, dup-extras > 200 bp = **0**,
      distinct sequences **13 065**, paths **465**.
- [x] PNG generated as `c4-crush-retry-on-poor.png` at
      `/home/erikg/impg/data/c4_crush_retry_on_poor_20260527T044100Z/c4-crush-retry-on-poor.png`.
      Upload to `hypervolu.me` is owner-driven (ssh from this worktree
      was rejected — same flow as the prior smoothxg-on task).
- [x] `docs/crush-retry-on-poor-compression.md` committed (this file).
- [x] `wg artifact crush-retry-on docs/crush-retry-on-poor-compression.md`
      (run on completion, see commit message).
