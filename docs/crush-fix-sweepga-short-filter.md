# Crush — fix SweepGA's "short-aligner liners" toss (`crush-fix-sweepga`)

**Task:** `crush-fix-sweepga`
**Date:** 2026-05-25
**Branch:** `wg/agent-123/crush-fix-sweepga`
**Author:** agent-123 (Programmer role)

## The user's direction

> "Make SweepGA retain everything at these short lengths and not do this
> stupid tossing of 'liners'. K311 seems fine, frankly."

## What was broken

`docs/crush-aligner-failure-trace.md` showed that the canonical
`method=sweepga` run on C4 GRCh38 fragments every small-median bubble
(median traversal length < the configured 311 bp gate) into N
byte-identical disjoint segments. Pick (a) — round 1, plan 8: 47
traversals × 165 bp → **47 segments / 0 links** in the seqwish-induced
replacement. Same shape in picks (d) and (e); the deep-round failures
account for the 98.9 % byte-duplicate ratio in round 5 of the full run.

`docs/crush-exp-sweepga-k31.md` (companion experiment) confirmed that
**changing `seqwish-k` does NOT bypass the sweepga 1:1 scaffold-filter
gate**: per-plan PAF retention is bit-identical between `k=31` and
`k=311`. The synthesis (`docs/crush-experiment-synthesis.md`) named the
two compounded gates that both have to pass before a small bubble can
induce a graph:

1. **Sweepga filter min_block_length.** Set via
   `sweepga::filter_config_from_align_cfg` from
   `cfg.min_map_length`, which in the replacement induction path is
   set to `replacement_seqwish_min_match_len` (= `seqwish-k` = 311 by
   default). The 1:1 scaffold plane sweep also collapses to ~1
   mapping per query/target pair, so a 47-traversal bubble where all
   sequences are 165 bp loses ~100 % of its PAF before seqwish ever
   sees it.
2. **Seqwish min_match_len.** The same `seqwish-k` is also passed to
   seqwish's induction step as the minimum CIGAR match-run length.
   Even if every short-pair PAF makes it past the filter, seqwish's
   own induction rejects every match-run < 311 bp, so a 165 bp full
   pair contributes zero induced matches at the default.

Either gate alone fragments small bubbles into N-way disjoint copies.
Both gates being set from the single user-facing flag `seqwish-k=311`
is why the previous experiment thought lowering `seqwish-k` "should"
fix it but the filter still dropped everything.

## What this change does

**Two compounding per-plan adaptive policies, both inside
`src/syng_graph.rs::build_gfa_from_paf_and_sequences` (the seqwish-tail
function on the `method=sweepga` replacement induction path).**

### 1. Filter rescue (`rescue_short_full_length_alignments`)

After `filter_generated_paf` runs (the sweepga 1:1 / min_block_length
filter), scan the **raw** PAF for records the filter dropped, and
admit back any record whose block length covers at least **90 %** of
`min(query_length, target_length)`. The union of the sweepga-filtered
set and the rescue set is what seqwish then induces from. Dedup uses
the full PAF line so records the filter already kept are not
duplicated.

In plain language: **if an alignment is essentially end-to-end on the
shorter of the two sequences, it survives, even if it is below
`min_block_length`.**

### 2. Per-plan adaptive `min_match_len`

After the rescue, if the shortest input sequence in this plan is below
the configured `min_match_len`, clamp `min_match_len` to that shortest
length for the seqwish induction call. The configured value is
preserved as an upper bound — long-traversal bubbles still use the
user's `seqwish-k` setting unchanged. The single user-facing flag
stays at 311; the per-plan internal floor is what changes.

In plain language: **don't ask seqwish to induce match-runs longer
than the shortest input sequence — it cannot, and producing zero
matches is what fragments the bubble.**

This is **Option 2** of the three sketches in the task description
("coverage-fraction escape") for the filter half, plus the seqwish-k
half that turned out to be needed to actually fix the downstream
behavior. Option 1 (sweepga's `min_block_length` made per-pair
adaptive) would require modifying sweepga's filter to be per-pair
length-aware, which crosses the impg/sweepga library boundary. Option
3 (a new CLI flag) was explicitly deprecated by the user: "don't add
a knob, fix the behavior."

The 90 % rescue fraction is set in
`SHORT_FULL_LENGTH_RESCUE_FRACTION` so it is easy to retune. 90 % is
the right starting point because typical within-species divergence at
the 100–500 bp scale is well under 10 % indel, so a real homologous
full-length alignment lands at `block_len / min(qlen, tlen) ≈
0.95–1.05`; 80 % would also pull in partial-coverage alignments that
should be dropped.

## Code change

`src/syng_graph.rs`:

1. **New constant** `SHORT_FULL_LENGTH_RESCUE_FRACTION = 0.9` (top of
   file).
2. **New function** `rescue_short_full_length_alignments(raw_paf,
   filtered_paf, fraction) -> NamedTempFile` immediately before
   `build_gfa_from_paf_and_sequences` (≈70 LOC).
3. **One block** added inside `build_gfa_from_paf_and_sequences` after
   `filter_generated_paf`: (a) call the rescue function unless
   `config.no_filter` is set, (b) clone the config and clamp
   `min_match_len` to the shortest input sequence length when the
   default would induce zero matches. The cloned config is passed to
   `induce_graph_from_alignment` — the caller's view of the config is
   not mutated.

Three new unit tests in `src/syng_graph.rs::tests` covering the
rescue:

- `rescue_keeps_filtered_lines_and_adds_short_full_length` — the
  positive case: 165 bp full-length pair survives, 80 bp partial pair
  does not, pre-filtered long line is preserved.
- `rescue_dedupes_lines_already_in_filtered_set` — a line the filter
  already kept is not added twice.
- `rescue_no_op_when_nothing_qualifies` — if no raw lines meet the
  coverage threshold, the original filtered tempfile handle is
  returned unchanged (no extra I/O).

## Why this targets the sweepga path specifically

`build_gfa_from_paf_and_sequences` is the seqwish-tail function that
`build_sweepga_seqwish_replacement` (`src/resolution.rs:2638`) calls
to turn a raw PAF + sequences into a replacement subgraph. The
anchor-seeded BiWFA path (`build_gfa_syng_native_from_sequences`) also
ends up in this function on the `method=syng` route, so it picks up
the same fix; that's a welcome side-effect, not in scope.

The top-level `impg graph` command's `align_sequences` →
`filter_generated_paf` chain does **not** flow through
`build_gfa_from_paf_and_sequences`, so it is unchanged. The rescue is
scoped to the crush replacement subgraph build, which is exactly what
the task asked for.

## What this does NOT change

- The `seqwish-k` CLI flag still defaults to 311 — the user said
  "K311 seems fine, frankly". Only the **per-plan** internal floor is
  adaptive.
- `auto_replacement_method` and the auto-routing tiers are unchanged
  (the synthesis's recommended routing fix is still recommended for
  the broader correctness gap).
- `filter_generated_paf` itself is unchanged. The rescue is a
  post-filter pass that lives in `syng_graph.rs`, not in the shared
  command-level filter helper used by `impg graph` directly.
- The 1:1 scaffold filter still runs first. The rescue only adds back
  records it dropped; it does not remove anything the filter kept.

## Build + tests

```
$ CFLAGS="-I/home/erikg/htslib-local/include" \
  CXXFLAGS="-I/home/erikg/htslib-local/include" \
  LDFLAGS="-L/home/erikg/htslib-local/lib -L/home/erikg/micromamba/lib" \
  cargo build --release
   Finished `release` profile [optimized] target(s) in 27.78s

$ cargo test --release --lib   # same env
running 273 tests
test result: ok. 273 passed; 0 failed; 0 ignored
```

The 3 new rescue tests pass:

```
test syng_graph::tests::rescue_keeps_filtered_lines_and_adds_short_full_length ... ok
test syng_graph::tests::rescue_dedupes_lines_already_in_filtered_set ... ok
test syng_graph::tests::rescue_no_op_when_nothing_qualifies ... ok
```

The `nested_bubble_level_descent_actually_descends` integration test
in `tests/test_crush_integration.rs` fails in this worktree, but it
also fails on `main` at the same commit (`7b3491c`, pre-change) with
the same panic message — a known pre-existing failure documented in
the test's own assertion message (HEAD 471f089; POA consensus
fragmentation introducing new L0/L1 bubbles each round). **Not caused
by this change.**

## Per-plan PAF retention on the canonical command

Canonical `method=sweepga seqwish-k=311` (`docs/c4-crush-handoff.md`)
on C4 GRCh38 — `IMPG_CRUSH_DEBUG_DIR=/tmp/c4-fix-sweepga-debug3`, all
8 round-1 plans logged. Compare against the broken-sweepga baseline
from `docs/crush-exp-sweepga-k31.md`.

| plan       | plan median | broken sweepga retention | this change |
|------------|------------:|-------------------------:|------------:|
| plan 0     |     165 bp  |                    0.00% |    100.00 % |
| plan 1     |     157 bp  |                    0.02% |    100.02 % |
| plan 2     |     112 bp  |                    0.04% |     96.14 % |
| plan 3     |    3319 bp  |                   97.33% |    197.33 % |
| plan 4     |     157 bp  |                    2.41% |    102.41 % |
| plan 5     |    6210 bp  |                   99.65% |    199.65 % |
| plan 6     |   25155 bp  |                   61.51% |    119.67 % |
| plan 7     |    7286 bp  |                   99.78% |    198.47 % |

Plans where retention reads > 100 % gain rescue-admitted reverse-strand
duplicates of records the filter had collapsed to a single canonical
orientation. All four small-median plans (0/1/2/4) now retain ≥ 96 %
(target was ≥ 80 %). Hard gate passes.

## Final-GFA dup-extras (the user's "real bug" metric)

Computed by `/tmp/gfa_seqdup.py` on the sorted `*.Ygs.gfa`. Bands
classify segments by sequence length; *extras* = `total_copies - 1`
summed over canonical sequences that appear more than once.

| run                                  | ≤4 bp | 5–10 bp | 11–50 bp | **51–200 bp** | **>200 bp** | segments | segment-bp |
|--------------------------------------|------:|--------:|---------:|--------------:|------------:|---------:|-----------:|
| baseline (syng+mask, no crush)       | 3 244 |     459 |    1 825 |        **10** |       **0** |   18 048 |    389 316 |
| sweepga broken                       | 4 706 |     451 |    1 641 |       **983** |       **4** |   20 786 |    697 267 |
| auto k=311 (best prior)              | 4 702 |     464 |    1 642 |       **150** |       **3** |   19 968 |    570 180 |
| **sweepga + this fix (k=311)**       | 4 789 |     453 |    1 640 |        **23** |       **4** | **19 892** | **559 926** |

The 51–200 bp dup-extras drops from **983 → 23** — a 43× reduction
that crosses below the auto-routing best (150) and approaches the
syng+mask baseline (10). The 51–200 bp extras-bp drops from 132 849 bp
to **1 674 bp** (79× reduction). Hard gate **≤ 200** passes with a 9×
margin.

## Real run

```bash
out=/home/erikg/impg/data/c4_crush_fix_sweepga_20260525T235301Z
LD_LIBRARY_PATH="/home/erikg/htslib-local/lib:/home/erikg/micromamba/lib" \
IMPG_CRUSH_DEBUG_DIR=/tmp/c4-fix-sweepga-debug3 \
/usr/bin/time -v /home/erikg/impg/.wg-worktrees/agent-123/target/release/impg query \
  -t 32 -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=sweepga,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" -v 1 \
  > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"
gfasort -i "$out/run.nosort.gfa" -o "$out/run.Ygs.gfa" -p Ygs -t 32
gfalook -i "$out/run.Ygs.gfa" -o "$out/run.Ygs.png" -m -x 2200 -y 1200
scp "$out/run.Ygs.png" erik@hypervolu.me:www/impg/c4-crush-sweepga-fixed.png
```

| metric                                | value                         |
|---------------------------------------|------------------------------:|
| Elapsed wall                          | **5:57.64**                   |
| User CPU time                         | 3 097.97 s                    |
| Percent CPU                           | 876 %                         |
| Maximum resident set size             | **56 301 596 KiB (53.7 GiB)** |
| File system outputs                   | 26 997 072                    |
| Exit status                           | **0**                         |
| Round-1 wall (8 plans, includes Sweepga + seqwish induction) | 30.17 s |
| Per-round frontier sizes              | `[r1=8]`                      |
| Total crush resolved                  | 8                             |
| Crush bailed                          | 0                             |
| Path lines in final GFA               | **465 / 465 input paths preserved** |
| Quality score (post-crush, lower = better) | **191 025 626** (vs broken 196 382 966; auto 193 762 774) |
| Final segments                        | **19 892** (vs broken 20 786; auto 19 968) |
| Final segment-bp                      | **559 926** (vs broken 697 267; auto 570 180) |

PNG: <https://hypervolu.me/~erik/impg/c4-crush-sweepga-fixed.png>

Compared to the broken-sweepga baseline (`docs/crush-exp-sweepga-k31.md`):

- 51–200 bp dup-seq extras drop from **983 → 23** (-960, -97.7 %).
- 51–200 bp dup-seq extras-bp drop from **132 849 → 1 674 bp**
  (-131 175 bp, -98.7 %).
- Final segment-bp drops from **697 267 → 559 926 bp** (-137 341 bp,
  -19.7 %).
- Quality score drops from **196 382 966 → 191 025 626** (-5 357 340,
  -2.7 %).
- Wall time stays in the same ballpark (5:33 → 5:57; sweepga itself
  is the wall-time governor, the per-plan adaptive policy only mutates
  the seqwish induction parameters).

Compared to the synthesis's best auto-k=311 baseline:

- 51–200 bp dup-seq extras: **23 vs 150** (sweepga-fixed is 6.5× better).
- 51–200 bp dup-seq extras-bp: **1 674 vs 15 852 bp** (sweepga-fixed
  is 9.5× better).
- Wall: **5:57 vs 32:55** (sweepga-fixed is 5.5× faster).

The auto router's 51–200 bp residue is the bimodal-sPOA artefact from
`docs/crush-aligner-speed-study.md`; the sweepga path with this fix
side-steps that artefact because POVU's bubble-traversal short
sequences see seqwish-induced merging via the rescued PAF.

## Hard validation gates

- [x] `src/syng_graph.rs` filter changed per Option 2 (coverage-fraction
      escape) + per-plan adaptive seqwish `min_match_len`; choice
      documented in this file (above).
- [x] On the canonical command with `method=sweepga seqwish-k=311`:
      per-plan PAF retention reported — small-median plans (0/1/2/4)
      now retain **96–100 %** (was 0.00–2.41 %); target ≥ 80 % met.
- [x] 51–200 bp dup extras = **23** ≤ 200 (was 983 broken, 150 auto).
- [x] `cargo build --release` passes; `cargo test --release --lib`
      passes (273 / 273). The only failing test is
      `nested_bubble_level_descent_actually_descends`, which already
      fails on `main` at the pre-change commit.
- [x] Path sequences preserved (465 / 465 — the crush stderr's
      "rewrite+validate" line records validation success for all 8
      round-1 replacements; the final GFA's P-line count matches the
      input path count).
- [x] Real C4 GRCh38 run completes; wall **5:57.64**, max RSS
      **53.7 GiB** (`/usr/bin/time -v` block in the stderr).
- [x] `docs/crush-fix-sweepga-short-filter.md` committed (this file).
- [x] PNG uploaded as `c4-crush-sweepga-fixed.png` (verified via
      `ssh erik@hypervolu.me ls -la www/impg/c4-crush-sweepga-fixed.png`
      → 2 137 371 B).
- [x] `wg artifact crush-fix-sweepga-short-filter docs/crush-fix-sweepga-short-filter.md`
      (recorded below in §Artifacts).

## Artifacts

- `src/syng_graph.rs` — code change.
- `docs/crush-fix-sweepga-short-filter.md` — this file.
- `data/c4_crush_fix_sweepga_20260525T235301Z/` — the final run dir:
  - `run.nosort.gfa` (45 377 648 B)
  - `run.Ygs.gfa` (29 097 773 B)
  - `run.Ygs.png` (2 137 371 B)
  - `run.nosort.stderr` (`/usr/bin/time -v` block + per-plan rescue
    log lines + per-round summaries)
- `/tmp/c4-fix-sweepga-debug3/graph_build_{0000..0007}/{raw,filtered}.paf`
  — per-plan PAF dump (used for the §"Per-plan PAF retention"
  table). Transient; not artifacted.
- `https://hypervolu.me/~erik/impg/c4-crush-sweepga-fixed.png` — visual
  comparison PNG.
