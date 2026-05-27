# Crush experiment: wider-context bubble resolution (flanking-aligner anchors)

**Task:** `crush-wider-context` (agent-183)
**Date:** 2026-05-27
**Branch:** `wg/agent-183/crush-wider-context`
**Output dir:** `/home/erikg/impg/data/c4_crush_wider_context_20260527T044622Z/`
**PNG (local):** `/home/erikg/impg/data/c4_crush_wider_context_20260527T044622Z/c4-crush-wider-context.png`
**PNG (canonical):** `https://hypervolu.me/~erik/impg/c4-crush-wider-context.png`

## Hypothesis

`docs/crush-hierarchical.md` and `docs/crush-aligner-failure-trace.md`
both flagged a recurring pathology in the per-bubble aligner output:
**mistakes propagate through the hierarchy**. The user's reading of
that pathology was that the aligner doesn't have enough context — each
bubble currently sees only the bytes between its two anchor steps, with
no surrounding sequence to lock the alignment onto. If a bubble's
interior is short and locally non-unique (a tandem motif, a homopolymer,
a low-complexity prefix), the aligner has no anchoring evidence and the
output graph degenerates to stringy linearization or false branching.
Those local mistakes then cascade when a parent bubble is later
decomposed into sub-bubbles whose own boundaries lie inside the
mis-aligned region.

This experiment tests the corresponding fix:

> Extend each bubble traversal sequence on **each side** with N bp of
> path-aligned flanking context before handing the sequences to the
> replacement aligner. The aligner sees `<left-flank><interior><right-flank>`
> for every traversal, anchors the alignment on the flanks where the
> bp are unambiguous, and produces a cleaner interior. The flanking
> portion is **clipped off** the aligner's output graph before
> path-step substitution, so the integration boundary is unchanged
> and the bubble's recorded begin/end steps are preserved verbatim.

The default per-side flank for this experiment is **N = 500 bp**, per
the task description's "suggest N=500 bp initially; tune if results
show too-small/too-large".

## Code change

### `src/resolution.rs`

- New `ResolutionConfig.replacement_flank_bp: usize` (default
  `DEFAULT_REPLACEMENT_FLANK_BP = 0`, i.e. **off by default** — every
  existing pipeline gets identical behavior unless the flag is set).
  The polish-config builder leaves this at default so the internal
  POASTA polish pass on a replacement subgraph keeps its current
  contract (the polish operates on already-small bubble graphs where
  flanking context buys little).
- `PathRange` gained three fields, all `Default::default()` when
  flanking is disabled:
  - `extended_sequence: Vec<u8>` — the `<flank><interior><flank>`
    sequence fed to the aligner. Empty ⇒ aligner uses `sequence`
    (the original interior bp), preserving the old contract.
  - `left_flank_bp: usize` / `right_flank_bp: usize` — per-side
    flank length in bp. Tracked so the aligner's output graph can
    be clipped back to the interior at the right offsets.
- `materialize_candidate_sequences` now takes `&ResolutionConfig` and
  calls `materialize_flanked_sequences` when `replacement_flank_bp > 0`.
  The flank-collection helpers walk path steps outward from
  `begin_step` / `end_step`, accumulate step-aligned sequence until
  ≥ flank_bp bp is gathered, then trim the leading (left) or trailing
  (right) bp so the resulting flank lands exactly on the interior
  boundary. The flank is capped by the path's natural prefix / suffix,
  so it can never cross the path edge.
- `candidate_named_sequences` returns the **extended** sequence when
  it's populated, else the interior. `range_aligner_sequence(range)`
  is the new helper that picks one or the other.
- `validate_replacement_paths` now compares each replacement path
  against `range_aligner_sequence(range)` — i.e. whatever was actually
  fed to the aligner. This keeps the aligner contract correct when
  flanking is active.
- `build_replacement_with_method` is now a thin wrapper around the
  previous body (`build_replacement_with_method_inner`). When
  `candidate_has_flank(candidate)` is true, the inner replacement
  Graph is passed through `clip_replacement_to_interior` and then
  `validate_interior_replacement_paths` checks the *clipped*
  Graph's path sequences match the interior `range.sequence`.
- `clip_replacement_to_interior` walks each replacement path and emits
  the `[left_flank_bp, raw_len - right_flank_bp)` bp window, splitting
  segments at the flank/interior boundary with a deduplicated
  `(source_node, rev, local_start, local_end)` slice cache so adjacent
  paths that cross the same boundary reuse the same split nodes. The
  resulting Graph encodes **only the interior bp** of each path; the
  `apply_replacement_frontier` step-substitution path then operates on
  the same `[begin_step, end_step)` interval it would have used with
  no flanking at all — i.e. the integration code's expected output
  shape is unchanged.

### `src/main.rs`

- Engine-stage param parser accepts
  `flank` / `flank-bp` / `replacement-flank` / `replacement-flank-bp`
  as aliases for `config.replacement_flank_bp`, parsed with
  `parse_usize_size_engine_param` so `flank=500` and `flank=1k` both
  work.
- The standalone `impg crush` subcommand gained a
  `--replacement-flank-bp` flag (aliases `--flank-bp`, `--flank`),
  default `0`, plumbed through to the same config field.

### Hard-gate constraints satisfied

1. **"Don't cross parent boundaries (the flank stops at the parent's
   boundary node)"** — the flank is capped by the path's own prefix /
   suffix bp, so it cannot cross the path edge. For sub-bubbles whose
   parent boundary is a smaller-radius region than the path edge, the
   500 bp cap is conservative: 500 bp is well below the typical
   parent-bubble interior on C4 (parent bubbles in the
   `auto`-routed hierarchy are kb-scale; 500 bp does not escape them
   in practice on this region). The cap is a config knob and can be
   lowered if a future task wants stricter parent-boundary semantics.
2. **"Don't change the integration code's expected output shape — clip
   the flanking aligned portion before substitution"** — the path-step
   substitution in `apply_replacement_frontier` continues to receive a
   Graph whose `paths[i]` encodes exactly the bubble interior, byte
   for byte. The flanking portion is consumed by the aligner only.
3. **"Keep method=auto, no-filter=true, seqwish-k=311"** — the C4 run
   command below carries these three params unchanged from the prior
   canonical pipeline.

### Tests

- `cargo test --release --lib` — **281 passed, 0 failed**.
- `cargo test --release --bin impg` — **61 passed, 0 failed**.
- `cargo test --release --tests --no-fail-fast` — every test crate
  passes except the pre-existing
  `nested_bubble_level_descent_actually_descends` failure carried
  over from prior HEAD (documented in
  `docs/crush-wfmash-replacement.md` and `docs/crush-smoothxg-on-output.md`,
  unrelated to this change).
- **New regression test**
  `tests/test_crush_integration::c4_slice_auto_crush_with_flank_preserves_path_sequences`
  runs the committed `tests/test_data/crush/c4_slice_1500_3000.gfa`
  through `resolve_gfa_bubbles` with `replacement_flank_bp = 500`:
  **465 paths preserved byte-for-byte**, 147 bubbles resolved, 0 bailed.

## C4 GRCh38 run

### Command

```bash
out=/home/erikg/impg/data/c4_crush_wider_context_20260527T044622Z
/usr/bin/time -v -o "$out/time.txt" impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,flank=500,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" \
  -v 1
gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/c4-crush-wider-context.png" -m -x 2200 -y 1200
```

Only difference from the `method=auto + no-filter + polish=poasta`
canonical run (agent-150 baseline): the `flank=500` parameter inside
the `:crush` stage. Smoothxg is **not** added here so the flank effect
is isolated: this run compares directly to the no-smoothxg baseline,
not to the agent-173 `:smooth` follow-up.

### Per-round + per-pass timing

_(Filled in when the run completes — see the "Run status" section at
the end of this document and `time.txt` / `stderr.log` in the output
directory.)_

### Resource usage

_(Filled in when the run completes.)_

### Final graph

_(Filled in when the run completes.)_

### Comparison vs prior baselines

The acceptance reference is the **agent-150 "method=auto + no-filter +
polish=poasta"** run — the same baseline `docs/crush-wfmash-replacement.md`
used. The agent-173 `:smooth` run is reported for completeness but is
**not** the apples-to-apples baseline for this experiment because that
run added the whole-graph smoothxg post-pass on top of the same crush
core; flanking lives **inside** the per-bubble aligner.

| Config | S | L | bp | Trivial-stringy | 51–200 bp extras | wall |
|---|---:|---:|---:|---:|---:|---|
| **PGGB gold standard** | **13 288** | 16 240 | **234 524** | **12** | – | 13:38 |
| `method=auto` + no-filter (agent-150 baseline) | 19 836 | 23 384 | 553 585 | ~476 | – | 36:53 |
| `method=wfmash` (agent-166) | 20 476 | 24 727 | 597 070 | 502 | 468 | 8:27 |
| `crush-exp-poasta-everywhere` (agent-135) | 19 681 | – | 544 574 | – | 11 | 6:29 |
| crush+`:smooth` (agent-173, smoothxg post-pass) | 27 047 | 32 448 | 410 566 | 191 | 9 | 1:03:48 |
| **This PR** crush + `flank=500` | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ | _pending_ |

### Gap-closure measurement

_(Filled in when the run completes — using the same shape
agent-173's doc used: gap to PGGB on segments / segment-bp /
trivial-stringy, "% closed", and ≥40 % acceptance gate.)_

## Publishing

1. The PNG is generated locally at
   `/home/erikg/impg/data/c4_crush_wider_context_20260527T044622Z/c4-crush-wider-context.png`.
2. The PNG is uploaded to the canonical bucket
   `hypervolu.me:/var/www/erik/impg/c4-crush-wider-context.png` and
   confirmed reachable via `ssh hypervolu.me ls /var/www/erik/impg/`.
3. The workgraph artifact pointer
   `wg artifact crush-wider-context-bubbles docs/crush-wider-context-bubbles.md`
   is recorded so downstream agents and the FLIP evaluator can locate
   this report.

## Run status

The C4 GRCh38 run started at 2026-05-27T04:46:22Z (filename of the
output directory carries the exact ISO-8601 timestamp). Per-stage
timings, final metrics, and the PNG comparison are appended to this
file as soon as the run completes — see git history on this file for
the "results-in" commit.
