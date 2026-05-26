# Crush Dispatch Floor + Final Gfaffix Compaction

Lands the two diagnostic-converged fixes for `syng:crush` on C4:

1. **Drop the min-traversal-len dispatch floor.** The canonical command
   pinned `min-traversal-len=5k`, which discarded ~62 of 70 POVU
   candidates per round in the deep-diag
   (`docs/crush-aligner-deep-diag.md` §Q5) and the unfiltered-sweepga
   report (`docs/crush-sweepga-everywhere-unfiltered.md`) — exactly the
   60–150 bp small bubbles that PGGB compacts and crush left as
   trivial-stringy. The canonical command in `docs/c4-crush-handoff.md`
   now uses `min-traversal-len=0`, sending small bubbles through the
   routing dispatch (sPOA / POASTA / sweepga) per the `method=auto`
   rules. The code default (`DEFAULT_MIN_TRAVERSAL_LEN`,
   `src/resolution.rs:223`) was already `0`, so this is a docs/canonical
   change.

2. **Add a final gfaffix compaction after crush.** PGGB ends in gfaffix;
   `syng:crush` did not (`src/lib.rs:745`). The PGGB comparison report
   (`docs/crush-vs-pggb-comparison.md` §Step 4) and the spec audit
   (`docs/crush-spec-audit.md` §Phase-8) both point at this. Even with
   dense alignments, PGGB-style compaction collapses byte-identical
   segments that independent crush plans produce across rounds. The new
   final pass calls `graph::run_gfaffix` (the in-tree workspace binary)
   immediately after `resolve_gfa_bubbles`; if the binary is unavailable
   (e.g. embedded library use or tests), it falls back to in-process
   `graph::unchop_gfa` so callers that never had access to gfaffix keep
   working.

## Code changes

### `src/lib.rs:745` — `apply_graph_transforms`

After the crush block, before the optional sort pipeline, the
post-resolution GFA is now passed through `run_gfaffix` (with an
`unchop_gfa` fallback when the gfaffix binary is missing). The pass
should be a no-op on a perfectly-crushed graph; on the current crush
output it collapses byte-identical small segments produced across
rounds.

### `docs/c4-crush-handoff.md`

Both `method=auto` canonical command blocks (the recorded historical
command and the "Run current HEAD semantics first" baseline) updated:
`min-traversal-len=5k` → `min-traversal-len=0`. A migration note dated
2026-05-26 records the rationale and links the diagnostic chain.

## Validation

- `cargo build --release --bins` — succeeds with both `impg` and
  `gfaffix` produced.
- `cargo test --release --lib` — 275/275 pass.
- `cargo test --release --test test_crush_integration` —
  `c4_slice_auto_crush_preserves_path_sequences` passes (path-preserving
  end-to-end); the `nested_bubble_level_descent_actually_descends` test
  was already RED on HEAD `9eaf71b` prior to this work (calls
  `resolve_gfa_bubbles` directly, unaffected by `apply_graph_transforms`
  changes).
- Real C4 GRCh38 query run with both fixes — see "C4 GRCh38 results"
  below.

## C4 GRCh38 results

Command (canonical, with both fixes active):

```bash
out=data/c4_dispatch_floor_gfaffix_20260526T135042Z
/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,aligner=fastga,min-traversal-len=0,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/dispatch_floor_gfaffix.nosort" \
  -v 1
```

(See "Resource use" below for wall + RSS; metrics filled in after the
run completes.)

### Before / after table

| metric                       | PGGB gold | crush current best | crush after both fixes |
|------------------------------|-----------|--------------------|------------------------|
| segments                     | 13,288    | 19,836             | _filled at completion_ |
| segment-bp                   | 234,524   | 553,585            | _filled at completion_ |
| trivial-stringy bubbles      | 12        | 476–500            | _filled at completion_ |
| 51–200 bp dup-extras         | _n/a_     | 11                 | _filled at completion_ |
| paths preserved              | 465/465   | 465/465            | _filled at completion_ |

### Resource use

```text
wall: _filled at completion_
max RSS: _filled at completion_
```

### Rendering

PNG: `c4-crush-dispatch-and-gfaffix.png` (uploaded after gfasort + gfalook).
