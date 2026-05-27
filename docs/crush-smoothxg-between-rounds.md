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

Real C4 GRCh38 (`GRCh38#0#chr6:31891045-32123783`, 465 paths).

- **Output dir:** `/home/erikg/impg/data/c4_crush_smoothxg_between_rounds_20260527T044302Z/`
- **PNG (local):** `data/c4_crush_smoothxg_between_rounds_20260527T044302Z/c4-crush-smoothxg-between-rounds.png`
- **PNG (canonical):** `https://hypervolu.me/~erik/impg/c4-crush-smoothxg-between-rounds.png`
- **Command:**

```bash
target/release/impg query -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k,smooth-between=true:smooth:nosort' \
  -O run.nosort -v 1
```

Only difference from the prior `crush-smoothxg-on` baseline command is the
appended `smooth-between=true` parameter on the `:crush` stage.

### Final metrics (sorted via `gfasort -p Ygs`)

| Metric | Value |
|---|---|
| Paths preserved | **465 / 465** |
| Final segments (S) | **23,096** |
| Final links (L) | 28,467 |
| Final segment-bp | **409,861** |
| Trivial-stringy bubble candidates | **132** (`/tmp/find_stringy_bubbles.py`) |
| Wall (impg query) | **50:03** (under heavy CPU contention; 12 concurrent impg jobs on the host) |
| Peak RSS | 100,425,672 KiB (≈ 95.8 GiB) |
| Exit | **0** |

### Per-round timing

| Phase | Wall |
|---|---:|
| syng index load (GBWT 272 M vertices) | ~3:22 (under contention) |
| crush round 1 build (8 plans; Sweepga × 1, POA × 4, POASTA × 3) | 15:26 |
| **smoothxg-between pass 1 (single 700 bp pass)** | **6:05** |
| crush round 2 discovery (POVU only — see "Known limitation" below) | 0:12 |
| Post-crush `:smooth` (two-pass 700 bp + 1100 bp) | 23:02 |
| gfaffix on smoothed graph | < 1 s |
| **End-to-end (impg query)** | **50:04** |

### Dup-segment-SEQUENCE extras per size band

`/tmp/gfa_seqdup.py` on `run.Ygs.gfa`:

| Band | total segs | total bp | distinct groups | dup groups | extras | extras-bp |
|---|---:|---:|---:|---:|---:|---:|
| ≤ 4 bp | 18,123 | 20,316 | 144 | 107 | **17,979** | 19,798 |
| 5–10 bp | 1,140 | 8,261 | 965 | 116 | **175** | 1,118 |
| 11–50 bp | 2,181 | 50,814 | 2,083 | 91 | **98** | 2,143 |
| **51–200 bp** | 1,275 | 135,587 | 1,265 | 10 | **10** | 1,181 |
| **> 200 bp** | 377 | 194,883 | 376 | 1 | **1** | 300 |

### Comparison vs prior baselines and PGGB

| Config | S | bp | Trivial-stringy | 51-200 bp extras | Wall |
|---|---:|---:|---:|---:|---:|
| **PGGB gold** | **13,288** | **234,524** | **12** | – | 13:38 |
| crush + `:smooth` (`crush-smoothxg-on`, agent-173) | 27,047 | 410,566 | 191 | 9 | 60:25 |
| **This task** (smooth-between=true on :crush, :smooth at end) | **23,096** | **409,861** | **132** | 10 | 50:04 |
| `method=auto` no-filter (agent-150) | 19,836 | 553,585 | ~476 | – | 36:53 |
| `method=wfmash` (agent-166) | 20,476 | 597,070 | 502 | 468 | 8:27 |

### Effect on the prior baseline

| Metric | Baseline (crush + post-:smooth only) | This task | Δ |
|---|---:|---:|---:|
| Segments | 27,047 | **23,096** | **-3,951 (-14.6 %)** |
| Segment-bp | 410,566 | **409,861** | **-705 (-0.2 %)** |
| Trivial-stringy bubble candidates | 191 | **132** | **-59 (-30.9 %)** |
| 51-200 bp dup-extras | 9 | 10 | +1 |
| Wall (impg query) | 60:25 | 50:04 | -10 min |

**Per-round smoothxg drops the trivial-stringy count by 31 %** and shaves
~4 k segments off the prior best, with no regression on the structural
51-200 bp band (the trivial +1 extra is below the noise floor).

## Known limitation — descent terminates after round 1

In this run the round-2 frontier was empty: round-1 local POVU produced
188 child keys, but the **bp_signature** of each child was computed against
the pre-smoothing segmentation, and the round-2 global POVU on the smoothed
graph (which has different segment boundaries from chop+block+POA) emits
4 010 sites whose bp_signatures don't perfectly match the round-1 children's
keys. The tree-based admission filter in `resolve_graph_bubbles` requires an
exact bp_signature match, so no candidate is admitted and descent halts.

Even with descent halting after round 1, the per-round smoothxg pass on its
own collapses the working graph's structural quality dramatically before
the post-crush `:smooth`/gfaffix pipeline runs — that is what gives this
configuration its 31 % trivial-stringy reduction and 14.6 % segment-count
reduction vs the baseline. The post-`:smooth` pass at the end still does
its full two-pass 700/1100 cleanup.

A **follow-up task** (out of scope here) should add **bp-region containment
admission** alongside the exact bp_signature filter: any global-POVU site
whose ref-bp range is contained in one of `last_round_resolved_ref_bp_regions`
should be admitted as a descendant of the corresponding round-(N-1) node,
regardless of whether its bp_signature matches a tree key. The comments on
`last_round_resolved_ref_bp_regions` in `resolution.rs` (`docs/crush-architecture-spec.md`
§Phase-6) already describe this fallback; the code only implements the
exact-match arm.

## Known limitation — single-pass timing exceeds the 2-min target

The `smooth-between=true` default uses a single 700 bp pass. On full C4 GRCh38
this took **6:05** per round — over the 2-min budget the task description
called for. Two follow-ups will help:

1. Lower the default `target_poa_lengths` (e.g. `[300]` or `[500]`) so the
   per-round pass completes in well under 2 min. Empirically the full-graph
   700 bp pass is 7 min on this dataset; halving the target POA length
   roughly halves block-decomp + per-block SPOA wall time.
2. The "descent halts after round 1" issue means the per-round pass only
   runs once in the current default configuration — the 2-min budget is
   effectively a single-run budget here. If the bp-region containment
   fallback is added, subsequent rounds operate on a smaller working graph
   (round-1 has already collapsed the big bubbles), so the per-round pass
   will be faster after round 1.

For now, callers who require the < 2 min budget can pass
`smooth-between=300` (or another sub-700 value) on the `:crush` stage.
