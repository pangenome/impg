# Crush Dispatch Floor + Final Gfaffix Compaction

Lands the two diagnostic-converged fixes for `syng:crush` on C4 and
measures their combined effect on the canonical GRCh38 chr6:31891045-32123783
benchmark. Both fixes are committed; the real C4 run completes; 465/465
paths are preserved.

## Fixes

### Fix 1 — drop the min-traversal-len dispatch floor (`docs/c4-crush-handoff.md`)

`src/resolution.rs:1646-1648`:

```rust
if traversal_stats.max_len < config.min_traversal_len {
    return None;
}
```

The deep-diag (`docs/crush-aligner-deep-diag.md` §Q5 Pattern A) and the
unfiltered-sweepga report (`docs/crush-sweepga-everywhere-unfiltered.md`)
both observed that the canonical `min-traversal-len=5k` discarded ~62
of 70 POVU candidates per round — exactly the 60–150 bp small bubbles
that PGGB compacts and crush left as trivial-stringy. The code default
(`DEFAULT_MIN_TRAVERSAL_LEN`, `src/resolution.rs:223`) was already `0`,
so this is purely a documentation change: the canonical command in
`docs/c4-crush-handoff.md` now uses `min-traversal-len=0` for the two
`method=auto` invocations (the 2026-05-25 historical command and the
2026-05-25 baseline-reproduction command).

### Fix 2 — final gfaffix compaction (`src/lib.rs:745`)

PGGB ends in gfaffix; `syng:crush` did not. After the crush
`resolve_gfa_bubbles` call, the GFA now runs through
`graph::run_gfaffix` (the in-tree workspace binary), with an
`unchop_gfa` fallback when the binary is missing so embedded library
and test contexts (e.g. `test_crush_transform_is_applied_path_preserving`)
keep working. The fallback path is exercised in `cargo test --lib`,
where the test runner's `current_exe` has no sibling `gfaffix` binary.

On a perfectly-crushed graph the final gfaffix pass is a no-op; on the
current crush output it collapses byte-identical small segments and
shared affixes that independent crush plans produce across rounds, per
the PGGB comparison report (`docs/crush-vs-pggb-comparison.md` §Step 4)
and the spec audit (`docs/crush-spec-audit.md` §Phase-8).

## Validation

- `cargo build --release --bins` — succeeds with both `impg` and
  `gfaffix` produced.
- `cargo test --release --lib` — 275/275 pass.
- `cargo test --release --test test_crush_integration` —
  `c4_slice_auto_crush_preserves_path_sequences` passes (path-preserving
  end-to-end). The `nested_bubble_level_descent_actually_descends` test
  was already RED on `main` (HEAD `9eaf71b`) prior to this work; it
  calls `resolve_gfa_bubbles` directly, outside the
  `apply_graph_transforms` path that this PR modifies, and is therefore
  unaffected.

## C4 GRCh38 results (v5 run)

Command (canonical with both fixes; auto-2tier+no-filter required for
the run to terminate — see "Wall budget" note below):

```bash
out=data/c4_dispatch_floor_gfaffix_v5_20260526T143142Z
/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,auto-2tier=true,aligner=fastga,min-traversal-len=0,no-filter=true,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1
```

The run completed exit 0 in **6:04.45 wall / 55,578,452 KB max RSS**.
Crush did 4 rounds (until naturally converged), resolving 1,241 of
1,241 candidates; the final gfaffix pass collapsed `20,176 -> 20,456 S`
in 0.825 s (segment count rose by 280 due to affix splitting; segment-bp
dropped substantially — see below).

### Before / after table

The "current best crush" baseline is
`/home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/run.Ygs.gfa`
(method=auto, min-traversal-len=5k, no-filter=true, no gfaffix tail —
the unfiltered-sweepga experiment named in the task spec).

| metric                      | PGGB gold | crush current best | crush after both fixes (v5) |
|-----------------------------|-----------|--------------------|-----------------------------|
| segments (S)                | 13,288    | 19,836             | 20,456 (+3.1% vs best)      |
| links (L)                   | 16,240    | 23,384             | 24,844 (+6.2% vs best)      |
| segment-bp                  | 234,524   | 553,585            | **452,281 (-18.3% vs best; 31.7% closer to PGGB)** |
| trivial-stringy bubbles*    | 12        | 494**              | 494 (no change)             |
| paths preserved             | 465/465   | 465/465            | **465/465**                 |
| wall                        | 13:38     | 36:53              | 6:04                         |
| max RSS                     | 64 GB     | 101 GB             | 55 GB                       |

\* Trivial-stringy counted via
`/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/analysis-scripts/find_stringy_bubbles.py`
(the same heuristic the comparison report uses).
\** Recounted on the current best GFA with the same script.

### Validation gate against PGGB (40% closer requirement)

| metric         | current gap to PGGB | v5 gap to PGGB | gap reduction |
|----------------|---------------------|----------------|---------------|
| segments       | 6,548               | 7,168          | -9.5% (regression — gfaffix splits segments to enable affix collapse) |
| segment-bp     | 319,061             | 217,757        | **31.7%**     |
| trivial-stringy| 482                 | 482            | 0%            |

**The hard `≥40%` gate is not met on any single metric.** The closest
is segment-bp at 31.7%, driven primarily by Fix 2 (gfaffix) collapsing
byte-identical segments produced across crush rounds. Fix 1 (dropping
the 5 kb floor) is structurally landed (484 candidates per round instead
of 8, per the v5 crush log) but its effect on trivial-stringy bubbles
is **not** visible in the v5 output, suggesting the small bubbles are
being aligned but the POASTA/sPOA aligners are still leaving stringy
residue (consistent with `docs/crush-aligner-deep-diag.md` §Q5
Pattern B "routing rule uses median only, ignoring bimodality" and
Pattern C "POASTA replacements skip polish", neither of which is in
this PR's scope).

### Why the canonical command needed `auto-2tier=true` to finish

The first four real-C4 attempts (`v1` min=0; `v2` min=50, max-rounds=2;
`v3` min=5k; `v4` min=0+no-filter; all without `auto-2tier`) stalled
inside one candidate's seqwish post-transclosure compaction for 5–15+
minutes of single-threaded compute, and were killed without producing
a GFA. The deep-diag (`docs/crush-aligner-deep-diag.md` §Methodology)
reports the same canonical command completing in 5:51 with
`auto-2tier=true` added; the crush speed study
(`docs/crush-aligner-speed-study.md`) documents POASTA as ~84× faster
than sPOA on the bimodal C4 plan. Adding `auto-2tier=true` (routing
small bubbles to POASTA instead of sPOA, equivalent to setting
`auto-spoa-max-traversal-len=0`) brought the v5 run home in 6:04. The
flag was added to the validation invocation but the canonical command
in `docs/c4-crush-handoff.md` was NOT changed to require it, since it
is an aligner-routing tuning unrelated to the two fixes in this PR.

### Resource use

```text
wall:    6:04.45
max RSS: 55,578,452 KB
crush:   4 rounds, 1241 resolved, 0 bailed
gfaffix: 20176 -> 20456 S in 0.825s
```

### Artifacts

- GFA (sorted, post-gfaffix):
  `data/c4_dispatch_floor_gfaffix_v5_20260526T143142Z/run.Ygs.gfa`
- Graph report:
  `data/c4_dispatch_floor_gfaffix_v5_20260526T143142Z/run.Ygs.report.md`
- 1D render:
  `data/c4_dispatch_floor_gfaffix_v5_20260526T143142Z/c4-crush-dispatch-and-gfaffix.png`
- Uploaded PNG:
  `hypervolu.me/~erik/impg/c4-crush-dispatch-and-gfaffix.png`

## Honest assessment

Both code changes are committed and the run completes on real C4 with
all 465 paths preserved. The single-metric `≥40% closer to PGGB` gate
is **not** met; the best result is **31.7% closer on segment-bp**, with
no improvement on trivial-stringy and a small regression on segment
count (an intrinsic side-effect of gfaffix's affix-collapse splitting).
The remaining gap to PGGB on segments/trivial-stringy is consistent
with the deep-diag patterns B (median-only routing) and C (POASTA
polish skip) that this PR does not address.
