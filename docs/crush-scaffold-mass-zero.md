# Crush experiment: `scaffold-mass=0` (disable the sweepga scaffold-mass filter for bubble-local replacement induction)

**Task:** `crush-scaffold-mass` (workgraph id `crush-scaffold-mass`)
**Date:** 2026-05-26
**Branch:** `wg/agent-160/crush-scaffold-mass`
**Output dir:** `/home/erikg/impg/data/c4_crush_scaffold_mass_zero_20260526T165458Z/`
**PNG:** `https://hypervolu.me/~erik/impg/c4-crush-scaffold-mass-zero.png`

## Hypothesis

User tip from a student: the SweepGA `scaffold_mass` default is set to
`10_000` (10 kb) for whole-genome alignment, and inside
`filter_config_from_align_cfg` it is then adaptively clamped to
`avg_seq_len × 3/5` for short bubble-local inputs (≈ 419 bp for the
plan-7-class bimodal inputs traced in
`docs/crush-fragment-source-trace.md`). The deep-diag noted that the
adaptive clamp drops 99.6 %–100 % of the raw PAF for those short
inputs.

`no-filter=true` already short-circuits the post-alignment plane-sweep
/ scaffold filter (`filter_generated_paf` returns its input unchanged
at `src/commands/graph.rs:677-685`), so the scaffold-mass filter is
already inactive on the canonical run. The student's hypothesis is
that the scaffold-mass filter should be disabled at the
**resolution-config layer** — so that even if a user toggles
`no-filter=false`, the bubble-local pipeline doesn't silently lose
99 %+ of its PAF to a whole-genome default that was never tuned for
~700 bp medians.

## Inventory: every `scaffold_mass` reference

### impg side (this repo)

| # | File:line | Kind | Current value / formula |
|---|---|---|---|
| 1 | `src/commands/mod.rs:153` | `FilterParams.scaffold_mass: u64` field (CLI-level filter params shared across `impg align`/`impg graph`) | — |
| 2 | `src/commands/mod.rs:163-198` | `build_filter_config`: derives `min_scaffold_length = round_nice(params.scaffold_mass.min(avg_seq_len * 3 / 5))` when `avg_seq_len > 0`, otherwise `params.scaffold_mass` | adaptive shrink for short inputs |
| 3 | `src/commands/graph.rs:78` | `GraphBuildConfig.scaffold_mass: u64` field (whole-graph builder config) | — |
| 4 | `src/commands/graph.rs:122` | `GraphBuildConfig::default()` initialiser | `10_000` (10 kb) |
| 5 | `src/commands/graph.rs:693` | `filter_generated_paf` constructs a `SweepgaAlignConfig` from `GraphBuildConfig` (only reached when `config.no_filter == false`) | passes `config.scaffold_mass` through to sweepga |
| 6 | `src/commands/align.rs:43,69,702,739,773,856` | `AlignConfig.scaffold_mass: u64` field + default(s) `10_000` + propagation into sweepga's `SweepgaAlignConfig` | `10_000` (10 kb) |
| 7 | `src/main.rs:8113,9038,9087` | CLI plumbing: `impg align` and `impg graph` forward `aln.sw.scaffold_mass` from the clap-parsed `--scaffold-mass` option into `AlignConfig`/`GraphBuildConfig` | user-controllable via `--scaffold-mass <bp>` |
| 8 | `src/resolution.rs:134` (**this PR**) | `ResolutionConfig.replacement_scaffold_mass: u64` field | **default `0`** (disabled — see `DEFAULT_REPLACEMENT_SCAFFOLD_MASS` at `src/resolution.rs:264`) |
| 9 | `src/resolution.rs:2480` (**this PR**) | `seqwish_replacement_config` wires `replacement_scaffold_mass` into `GraphBuildConfig.scaffold_mass` for the crush dispatch | passes `0` by default |
| 10 | `src/main.rs:2708-2710` (**this PR**) | Engine-string parser entry `replacement-scaffold-mass` / `scaffold-mass` inside the `crush` stage | user-overridable via `-o gfa:…:crush,scaffold-mass=N,…` |
| 11 | `src/syng.rs:4427` | `chain_anchors_with_sweepga_scaffold_mass(raw_hits, syncmer_len, chain_budget, min_anchors * syncmer_len)` — **syng-map** anchor chaining (`impg map` path; not used by crush) | passes `min_anchors * syncmer_len` |
| 12 | `src/lib.rs:348` | Same syng-map call, second use site | `effective_min * syncmer_len` |
| 13 | `src/syng_transitive.rs:104,113,1012,1655,1674` | Public `chain_anchors_with_sweepga_scaffold_mass` helper that constructs a `FilterConfig` with `min_scaffold_length = min_scaffold_length.max(syncmer_len)` (so 0 floored to syncmer length) | used only by `impg map` / `impg query` syng-anchor refinement |

The syng-anchor refinement path (rows 11-13) is **not on the crush
replacement-induction path**. The crush path goes
`build_sweepga_seqwish_replacement` (`src/resolution.rs:2752`) → raw
PAF → `build_gfa_from_paf_and_sequences`
(`src/syng_graph.rs:697`) → `filter_generated_paf`
(`src/commands/graph.rs:671`) → optional sweepga
`filter_config_from_align_cfg`.

### sweepga side (vendored, `ddd31d39b6a68fc972025b048076032341b66835`)

| # | File:line | Kind | Current value / formula |
|---|---|---|---|
| 1 | `src/library_api.rs:130` | `SweepgaAlignConfig.scaffold_mass: u64` field | — |
| 2 | `src/library_api.rs:164` | `SweepgaAlignConfig::default()` | `10_000` (10 kb) |
| 3 | `src/library_api.rs:223-258` | `filter_config_from_align_cfg`: clamps via `clamp_scaffold_params(jump, mass, Some(avg_seq_len), true)` then sets `FilterConfig.min_scaffold_length = scaffold_mass` | adaptive shrink |
| 4 | `src/library_api.rs:321-322` | `sweepga_align`: **returns raw PAF unchanged when `no_filter == true`** | scaffold-mass filter never runs |
| 5 | `src/pansn.rs:194-225` | `clamp_scaffold_params`: `mass = round_nice(min(user_mass, avg.saturating_mul(3) / 5))` when adaptive. **`round_nice(0) = 0`** (so `user_mass = 0` propagates through the clamp unchanged) | adaptive shrink |
| 6 | `src/paf_filter.rs:39` | `FilterConfig.min_scaffold_length: u64 // -S/--scaffold-mass` field | — |
| 7 | `src/paf_filter.rs:452` | Scaffold-chain length filter: `chain.total_length >= self.config.min_scaffold_length` — at `min_scaffold_length = 0` every chain passes | gate of the scaffold-mass filter |
| 8 | `src/paf_filter.rs:461` | Plane-sweep guarded by `if self.config.min_scaffold_length > 0 || self.config.min_scaffold_identity > 0.0` — at `min_scaffold_length = 0` the scaffold-mass branch is skipped entirely | gate of the scaffold-mass filter |

## Code change (this PR)

Two surfaces, neither of which affects the existing `impg align` /
`impg graph` defaults:

```rust
// src/resolution.rs (additions)
pub struct ResolutionConfig {
    …
    /// Minimum scaffold chain length (`min_scaffold_length` in the sweepga PAF
    /// filter) for pairwise replacement graph induction. Set to 0 to disable
    /// the scaffold-mass filter entirely. The crush default is 0.
    pub replacement_scaffold_mass: u64,
}

pub const DEFAULT_REPLACEMENT_SCAFFOLD_MASS: u64 = 0;

impl Default for ResolutionConfig {
    fn default() -> Self {
        Self {
            …
            replacement_scaffold_mass: DEFAULT_REPLACEMENT_SCAFFOLD_MASS,
            …
        }
    }
}

fn seqwish_replacement_config(
    config: &ResolutionConfig,
) -> crate::commands::graph::GraphBuildConfig {
    …
    crate::commands::graph::GraphBuildConfig {
        …
        scaffold_filter: config.replacement_scaffold_filter.clone(),
        scaffold_mass: config.replacement_scaffold_mass,   // ← new wiring
        …
    }
}
```

```rust
// src/main.rs (engine-string parser, parse_crush_stage)
"replacement-scaffold-mass" | "scaffold-mass" => {
    config.replacement_scaffold_mass =
        parse_usize_size_engine_param(raw, &param.key, &param.value)? as u64;
}
```

### Test coverage

Two new tests in `resolution::tests`:

- `replacement_scaffold_mass_defaults_to_zero_and_propagates` — asserts
  `ResolutionConfig::default().replacement_scaffold_mass == 0` and
  that `seqwish_replacement_config` propagates that to
  `GraphBuildConfig.scaffold_mass == 0`.
- `replacement_scaffold_mass_user_override_propagates` — asserts that
  setting `replacement_scaffold_mass = 500` carries through to
  `GraphBuildConfig.scaffold_mass == 500`.

`cargo test --release --lib` → **281 passed; 0 failed**.

## Run

### Command (verbatim)

```bash
out=/home/erikg/impg/data/c4_crush_scaffold_mass_zero_20260526T165458Z
mkdir -p "$out"
source ./env.sh
/usr/bin/time -v ./target/release/impg \
  query -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,scaffold-mass=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 \
  > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"
```

Differences from `crush-sweepga-everywhere-unfiltered.md`'s canonical
command: a single new engine-string parameter, `scaffold-mass=0`. The
binary additionally encodes that as the new **default** of
`ResolutionConfig::replacement_scaffold_mass`, so even an
identical-looking command line on this binary would carry the new
default; the explicit `scaffold-mass=0` is for documentation /
auditability.

### Results

| Field | This run (`scaffold-mass=0`) | Baseline (`crush-sweepga-everywhere-unfiltered`, agent-150) |
|---|---:|---:|
| Elapsed wall | **5:24.00** | 5:51.75 |
| User time | 2 770.28 s | 2 765.91 s |
| System time | 254.59 s | 268.20 s |
| Percent of CPU | 933 % | 862 % |
| Exit status | **0** | 0 |
| Maximum resident set size | 55 365 256 KiB (52.80 GiB) | 55 402 268 KiB (52.83 GiB) |
| Paths preserved | **465 / 465** | 465 / 465 |
| Per-plan aligner distribution | **8/8 Sweepga** (round 1; round 2 found 0 candidates and the loop terminated) | 8/8 Sweepga (identical) |
| Final segments (S) | **20 040** | 20 040 |
| Final links (L) | 23 904 | 23 904 |
| Final total bp | **544 479** | 544 479 |
| Distinct canonical sequences | 13 152 | 13 152 |
| 51–200 bp dup extras (sequence-dedup) | **14** | 14 |
| > 200 bp dup extras | **2** | 2 |
| 11–50 bp dup extras | 1 534 | 1 534 |
| 5–10 bp dup extras | 405 | 405 |
| ≤ 4 bp dup extras | 4 933 | 4 933 |
| seq-extras (segments − distinct) | 6 888 | 6 888 |

(`gfa_metrics2.py` from `data/c4_pggb_control_20260526T025439Z/analysis-scripts`; both
`run.nosort.gfa` files diff to the byte across the S/L/P-line ordering used by
the script.)

The "51–200 bp dup extras = 14 (vs the 15 reported in agent-150's
synthesis)" is a counting-script difference, not a regression: agent-150's
synthesis-table figure came from a separate ad-hoc band counter; the
canonical `gfa_metrics2.py` (used here and re-run against
agent-150's GFA to confirm) returns the same 14 for both runs.

### PAF retention

Per-plan PAF retention is **100 %** at every filter stage for both
runs (this PR and the agent-150 baseline). The proof is the same as
in `crush-sweepga-everywhere-unfiltered.md` §"PAF retention proof":
when `no-filter=true`, `filter_generated_paf` returns its input
unchanged at `src/commands/graph.rs:677-685`, so the scaffold-mass
filter is never invoked, so a `scaffold_mass=0` flag cannot remove
any additional lines. The "≥ no-filter baseline" gate in the task
brief is satisfied by equality.

## Finding

The task brief anticipated this outcome:

> If scaffold_mass=0 produces NO improvement: that's a finding (the
> filter was already disabled by no-filter=true). Document it.

**Confirmed.** For the canonical crush-everywhere run
(`no-filter=true`), `filter_generated_paf` returns early without
invoking `filter_config_from_align_cfg` at all
(`src/commands/graph.rs:677-685`). The scaffold-mass field is dead
code in that branch. The "this run" vs "agent-150 baseline" columns
in the Results table are **byte-identical on every metric the task
brief asks for**: paths (465/465), segments (20 040), segment-bp
(544 479), dup extras in the 51–200 bp band (14), > 200 bp dup
extras (2). The 27-second wall-clock difference (5:24 vs 5:51) is
FastGA scheduling jitter — both runs used 32 rayon threads and
otherwise identical engine settings.

This change is still **load-bearing for any future caller that turns
`no-filter=false`**: at that point the resolution-config default of
`replacement_scaffold_mass = 0` skips the scaffold-mass filter
entirely, instead of clamping to ≈ 419 bp for plan-7-class inputs
and dropping 99.6 %–100 % of the raw PAF (the failure documented in
`docs/crush-fragment-source-trace.md` §"Parameter values at the
failure point").

## Files changed

- `src/resolution.rs` — `ResolutionConfig.replacement_scaffold_mass`
  field, `DEFAULT_REPLACEMENT_SCAFFOLD_MASS = 0` constant,
  `Default::default()` initialiser, `seqwish_replacement_config`
  wiring, two new tests.
- `src/main.rs` — `parse_crush_stage` engine-string parser entry for
  `replacement-scaffold-mass` / `scaffold-mass`.
- `docs/crush-scaffold-mass-zero.md` — this document.
