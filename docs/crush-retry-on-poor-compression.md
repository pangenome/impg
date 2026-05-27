# crush: retry-on-poor-compression

Task: `crush-retry-on-poor-compression` (workgraph id `crush-retry-on`).
Branch: `wg/agent-182/crush-retry-on`.

## Problem

> If mistakes get made at first pass they can stick as we dig into the
> hierarchy. — task spec, user insight.

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

If the ratio is below `retry_min_compression_ratio` (default `2.0`) **and**
the bubble's input size is above `retry_min_input_bp` (default `1 000` bp,
to keep noisy 1-bp bubbles out of the retry path), rebuild the same bubble
with an *alternate* aligner. Keep whichever replacement is smaller in
segment-bp.

`alternate_replacement_method` in `src/resolution.rs` enforces a fixed
mapping so the cascade is deterministic and auditable:

| Original method   | Alternate |
|-------------------|-----------|
| `Sweepga`         | `Poasta`  |
| `Wfmash`          | `Poasta`  |
| `Poasta`          | `Allwave` |
| `Allwave`         | `Poasta`  |
| `Poa` (sPOA)      | `Poasta`  |
| `StarBiwfa`       | `Poasta`  |
| `Auto` / `Hierarchical` | (none — these are routing modes that always resolve to a concrete aligner before this point) |

## Safety invariants

1. **Bounded wall.** Each bubble gets at most ONE alternate attempt
   (original + 1 retry). There is no recursive cascade and no per-bubble
   loop.
2. **Path preservation.** Each aligner's builder
   (`build_poa_replacement`, `build_poasta_replacement`, etc.) calls
   `validate_replacement_paths` on its output before returning, so the
   retry can only swap in a replacement whose paths spell the same bytes
   as the input traversals. The outer `apply_replacement_frontier` re-runs
   `path_sequences_equal` after rewriting, providing a second-line
   defense.
3. **Never lose ground.** If the alternate's build errors, returns an
   empty replacement, or merely fails to beat the original on segment-bp,
   the original plan is kept. The retry path can only make a bubble's
   replacement smaller; it can never make it larger or invalid.

## Configuration

Exposed on the `crush:` engine pipeline (parsed in `src/main.rs`,
`parse_crush_engine_stage_params`):

- `retry-min-compression-ratio=<f64>` (alias: `retry-ratio`,
  `min-compression-ratio`, `compression-retry-ratio`). Default `2.0`.
  Set to `0.0` to disable the mechanism entirely.
- `retry-min-input-bp=<usize>` (alias: `retry-min-input`). Default
  `1 000`. Bubbles whose total input traversal bp is below this floor
  are not retried regardless of ratio — small bubbles where neither
  aligner can meaningfully compress the input would waste wall on
  retries that cannot improve the result.

## Reporting

`ResolutionStats` carries three new counters:

- `retry_attempts` — number of bubbles below threshold whose alternate
  was actually built (not just considered).
- `retry_wins` — alternate compressed strictly better than the original
  and was swapped in.
- `retry_failures` — alternate's build returned an error, an empty
  graph, or otherwise could not be evaluated; original was kept.

Per-bubble decisions are also written to the log at `info` level with the
short signature of the bubble, the original→alternate method pair, the
old/new segment-bp, and old/new ratios. A round summary line at the end
of each retry pass reports `attempts/wins/failures` for that round and
the running totals.

## Tests

`src/resolution.rs` (unit, under `mod tests`):

- `alternate_replacement_method_picks_a_different_aligner` — enforces
  the mapping table above and the "alternate must differ from original"
  invariant.
- `replacement_compression_ratio_handles_zero_and_basic_cases` — exact
  values on perfect (4 traversals → 1 segment, ratio 4.0) and no-op
  (4 traversals → 4 segments, ratio 1.0) cases; `None` returned when
  input or output bp is zero so the retry path treats them as "no
  signal" instead of dividing by zero or grading every empty graph as
  infinitely well compressed.
- `retry_preserves_path_sequences_on_simple_bubble` — end-to-end on the
  canonical SNP bubble fixture: with a 1e9 threshold every bubble
  becomes a retry candidate; output paths must still spell the input
  byte-for-byte.
- `retry_skips_bubbles_below_min_input_bp` — verifies the input-bp
  floor suppresses retries on tiny fixtures even when the ratio
  threshold is aggressive.

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

### Retry mechanism in action

From `stderr.log`:

```
crush round 1: retry-on-poor-compression: 1/8 bubble(s) below ratio=2.00 (min-input-bp=1000), attempting alternate aligner
crush round 1 retry: bubble sig=CHM13#0#chr6:31744284-31976975:143224-143389|...
  Poa→Poasta alternate LOST (orig=25230bp ratio=1.30, alt=25230bp ratio=1.30); keeping original
crush round 1 retry summary: attempts=1 wins=0 failures=0 (running totals: attempts=1 wins=0 failures=0)
```

Of the 8 accepted round-1 plans, 1 was below the 2.0 compression-ratio
threshold (input traversal bp / output segment bp = 1.30). The original
aligner was `Poa` (sPOA, picked by the median-based 3-tier dispatch for
a small-median bubble); the retry rebuilt it with `Poasta` and produced
the same output bp (25 230). Because the alternate did not strictly
beat the original, the original plan was kept (`AlternateLost`).
Rounds 2 and 3 saw no candidates below the threshold.

Total stats:
- `retry_attempts = 1`
- `retry_wins = 0`
- `retry_failures = 0`

The retry mechanism is conservative on this C4 fixture: the worst
round-1 compression (Poa, ratio 1.30) was already at the Poa/Poasta
crossover where neither aligner has an obvious advantage on a 165 bp
bubble with 47 traversals. With deeper / different inputs the
mechanism can fire more often; the same code path will swap in
whichever aligner produces the smaller segment-bp.

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

Because the retry only fired once and the alternate lost, the metric
deltas above mostly reflect the `method=auto` 3-tier dispatch vs
`method=sweepga`, not the retry path. The important property is that
enabling retry did not regress anything (paths preserved, segment-bp
went down not up, dup-extras either improved or stayed flat).

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

- [x] Retry mechanism implemented; `cargo test --release --lib resolution`
      passes 54 tests (50 pre-existing + 4 new). `cargo test --all`
      passes everything except the pre-existing
      `nested_bubble_level_descent_actually_descends` red test
      (documented as RED on HEAD before this branch and unrelated to
      the retry code path).
- [x] Real C4 GRCh38 run, **465 / 465** paths preserved.
- [x] Retry stats reported: `retry_attempts=1, retry_wins=0,
      retry_failures=0` (logged at info level; also accessible via
      `ResolutionStats`).
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
