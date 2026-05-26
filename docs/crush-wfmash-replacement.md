# Crush experiment: wfmash replacement aligner

**Task:** `crush-wfmash-replacement` (agent-166)
**Date:** 2026-05-26
**Branch:** `wg/agent-166/crush-wfmash-replacement`
**Output dir:** `/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T194829Z/`
**PNG (local):** `/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T194829Z/c4-crush-wfmash.png`
**PNG (canonical):** `https://hypervolu.me/~erik/impg/c4-crush-wfmash.png` (upload by user; see "Publish" below)

## User direction

> we should try WFMASH on the bubbles.

The best-looking output so far is the hierarchical run
(`c4-crush-hierarchical.png`) but it shows 5–6 top-level bubbles that
didn't reduce — one is a cyclic module that loops correctly externally
but doesn't compact internally. Sweepga+FastGA was barfing on these.
PGGB uses wfmash and gets **13,288** segments (vs our best ~20k).

The hypothesis from
`docs/crush-sweepga-many-to-many-poasta-polish.md` §Finding ranked the
aligner anchor granularity as the most likely cause of the 6,670-segment
gap vs PGGB. This experiment tests that hypothesis directly: same
canonical command, only the aligner is changed from FastGA to wfmash.

## Code change

Added a dedicated `ResolutionMethod::Wfmash` variant (alongside
`Auto`, `Poa`, `Poasta`, `StarBiwfa`, `Allwave`, `Sweepga`,
`Hierarchical`) at `src/resolution.rs`. It dispatches through the
same SweepGA+seqwish induction tail as `method=sweepga` but pins
the aligner field to `wfmash` regardless of the `sweepga_aligner`
default ("fastga"). The wrapper is a thin shim:

```rust
fn build_wfmash_seqwish_replacement(
    candidate: &BubbleCandidate,
    config: &ResolutionConfig,
) -> io::Result<Graph> {
    let mut wfmash_config = config.clone();
    wfmash_config.sweepga_aligner = "wfmash".to_string();
    build_sweepga_seqwish_replacement(candidate, &wfmash_config)
}
```

The underlying `sweepga::library_api::sweepga_align` already supports
both `aligner=wfmash` and `aligner=fastga`, so the entire wfmash
plumbing (binary discovery via `wfmash-rs::binary_finder`, FASTA
indexing, density/segment-length adaptation) is reused as-is. The
only new code is the method-level entry point, the `parse_name` alias
(`"wfmash" | "wf"`), and CLI/help updates.

### Rationale for `method=wfmash` vs `aligner=wfmash`

`method=sweepga,aligner=wfmash` already worked (the `aligner`
parameter has been wired through the sweepga path since
`crush-fix-sweepga` / agent-123). The dedicated `method=wfmash`
variant exists for two reasons:

1. **Per-plan reporting**: the crush config log line and the
   `building replacement N/M with {method:?}` log line both show
   "Wfmash" cleanly, which is what the canonical reproduction
   command and the hard-gate reports need.
2. **Defensive coupling**: the dispatch shim pins the aligner
   at the call site, so `method=wfmash` survives any future
   change to the `sweepga_aligner` default (currently "fastga")
   or any future refactor of the engine string parser.

The `method=sweepga,aligner=wfmash` path is unchanged and still
works — it remains a more granular knob for users who want to mix
aligner choices with other sweepga defaults independently.

### Environment note: wfmash binary

The `wfmash` binary in `target/release/wfmash` is **rebuilt from
the wfmash source vendored inside `wfmash-rs`'s `deps/wfmash/`**
(commit pinned by the wfmash-rs dependency, currently
`d47b7e3`). The wfmash binary that ships in the system `PATH`
(`/home/erikg/bin/wfmash`, `v0.24.1-24-gf64becc4`) has an
**incompatible CLI** in this environment: its `-s` flag is
`sketch-size` and is constrained by `sketch-size <= window-size`,
which rejects sweepga's `-s 5000` (segment-length intent) calls
with `ERROR: sketch size (-s) must be less than or equal to window
size (-w)`. The vendored wfmash (older, matching the wfmash-rs
expectation) has `-s` as `--segment-length`. The build picks up
the vendored source via the build script `build.rs` at the top
level, which copies the wfmash-rs-built binary next to
`target/release/impg`. To rebuild from scratch in a fresh worktree,
the htslib and jemalloc headers must be visible to cmake (this
worktree used `CMAKE_PREFIX_PATH=/gnu_old/store/dzrirrs2qsqsb05lz1f210g40n05l001-htslib-1.16:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0`
with matching `CFLAGS`/`LDFLAGS`).

The wfmash-rs↔wfmash CLI mismatch was the cause of the
**first failed run** (`/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T193157Z/`):
the system wfmash binary at `target/release/wfmash` came from a
newer wfmash revision than what wfmash-rs's CLI generator expects;
14 of 19 plans bailed with the sketch-size error. The second run
(below) uses the vendored, matching wfmash.

## C4 GRCh38 run

### Command

```bash
out=/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T194829Z
/usr/bin/time -v target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=wfmash,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
  -O "$out/run.nosort" \
  -v 1
```

Only difference from the canonical
`crush-sweepga-many-to-many-poasta-polish` command:

- `method=wfmash` (was `method=sweepga,aligner=fastga`)

### Per-plan aligner choice and per-plan output shape

| Plan | traversals | max-len | median-len | total-len | root-span | Outcome |
|----:|---:|---:|---:|---:|---:|---|
| 1 | 312 | 31 478 | 25 155 | 5 132 443 | 110 | **Wfmash** accepted |
| 2 | 437 | 42 362 | 157 | 2 893 782 | 157 | **Wfmash** accepted |
| 3 | 460 | 7 351 | 7 286 | 3 309 531 | 7 277 | **Wfmash** accepted |
| 4 | 356 | 9 320 | 3 319 | 1 181 720 | 3 319 | **Wfmash** accepted |
| 5 | 98 | 25 255 | 157 | 65 552 | 157 | **Wfmash** accepted |
| 6 | 47 | 25 230 | 165 | 32 820 | 165 | **Wfmash** accepted |
| 7 | 287 | 6 265 | 6 210 | 1 781 621 | 6 196 | **Wfmash** accepted |
| 8 | 325 | 25 274 | 112 | 188 057 | 112 | **Wfmash** accepted |

**8 / 8 = 100 %** plan retention. Identical to
`crush-sweepga-many-to-many-poasta-polish` (FastGA) — both aligners
emit non-empty replacements for every plan when `no-filter=true`. Polish
inside each replacement is **POASTA** (every nested sub-bubble within
the wfmash-induced graph polished by POASTA via `polish-method=poasta`).

Plan-1 (root span 110 bp, 312 traversals, total 5.1 Mb, median-len 25 155)
is the dominant bimodal locus that prior experiments
(`docs/crush-aligner-deep-diag.md`, `docs/crush-aligner-failure-trace.md`)
identified as the candidate locus for the short↔long anchoring problem.
wfmash mapped **177/312 reads** for this plan (mapping rate 56.7 %), wrote
**4.92 MB** of PAF, then transclosure spanning-tree compaction collapsed
27 001 pair candidates into 436 edges (61× reduction). The resulting
seqwish-induced replacement for plan-1 produced a single connected
component which, after POASTA polishing of its 159 nested sub-bubbles,
landed in the final graph.

Round 1 alone resolved all 8 candidates; round 2 found no eligible new
candidates and the loop terminated. **`crush per-round frontier sizes
= [r1=8]`**.

### Results

| Field | Value |
|---|---|
| Elapsed wall | **8:26.88** |
| User time | 9 052.78 s |
| System time | 1 742.57 s |
| Percent of CPU | 2 129 % |
| Exit status | **0** |
| Maximum resident set size | 55 435 380 KiB (52.86 GiB) |
| Paths preserved | **465 / 465** (`grep -c "^P\b" run.Ygs.gfa`) |
| Per-plan aligner distribution | **8/8 Wfmash** (round 1; round 2 found 0 eligible candidates) |
| Final segments (S) | **20 476** |
| Final links (L) | 24 727 |
| Final total bp | **597 070** |
| Total path steps | 4 484 814 |
| POVU sites on `Ygs.gfa` | 2 377 (level-0 = 1 102) |
| Trivial-stringy bubble candidates | **502** (`/tmp/find_stringy_bubbles.py` heuristic) |

### Dup-segment-SEQUENCE extras per size band

Computed by `/tmp/gfa_seqdup.py` on `run.Ygs.gfa` (forward and
reverse-complement collapsed; *extras* =
`total_segments_in_group − 1` summed across canonical-sequence
groups with more than one copy).

| band              | total segs | total bp | distinct groups | dup groups | **extras** | extras-bp |
|-------------------|-----------:|---------:|----------------:|-----------:|-----------:|----------:|
| ≤ 4 bp            |     5 033  |   7 200  |       150       |    117     |   **4 883**|     6 658 |
| 5–10 bp           |     1 949  |  14 469  |     1 478       |    326     |     **471**|     3 254 |
| 11–50 bp          |    12 459  | 348 892  |    10 772       |  1 348     |   **1 687**|    43 625 |
| **51–200 bp**     |       906  |  85 165  |       438       |     20     |     **468**|    52 062 |
| **> 200 bp**      |       129  | 141 344  |       129       |      0     |       **0**|         0 |

### Comparison vs prior baselines

| Config | S | L | bp | Polish | Aligner | 51–200 bp extras | >200 bp extras | wall |
|---|---:|---:|---:|---|---|---:|---:|---|
| **This PR** `method=wfmash`, no-filter=true, polish=poasta | **20 476** | 24 727 | **597 070** | POASTA | wfmash | **468** | **0** | 8:27 |
| `crush-sweepga-many-to-many-poasta-polish` (agent-163) | 19 958 | 23 822 | 544 641 | POASTA | FastGA | (n/a) | (n/a) | 5:25 |
| `crush-sweepga-everywhere-unfiltered` (agent-150) | 20 040 | 23 904 | 544 479 | sPOA | FastGA | (n/a) | (n/a) | (n/a) |
| `crush-hierarchical-sweepga` (agent-152) | 20 326 | 24 028 | 576 864 | sPOA | FastGA (L0) + POASTA (L≥1) | 26 | 4 | 5:45 |
| `crush-exp-poasta-everywhere` (agent-135) | 19 681 | – | 544 574 | POASTA | (POASTA every bubble) | 11 | 0 | 6:29 |
| `method=auto` + no-filter (prior) | 19 836 | 23 384 | 553 585 | (sPOA/POASTA/FastGA tiered) | tiered | 15–16 | 3 | 36:53 |
| **PGGB gold standard** | **13 288** | 16 240 | **234 524** | (PGGB pipeline) | wfmash + smoothxg + gfaffix | – | – | 13:38 |

## Finding

**wfmash on the bubbles does NOT compact toward PGGB's 13 288
segments.** With every other knob held to the canonical
many-to-many + no-filter + POASTA-polish configuration, swapping the
inner replacement aligner from FastGA to wfmash produced **+518
segments (+2.6 %), +52 429 bp (+9.6 %)**, and dramatically WORSE 51–200 bp
canonical-sequence duplication (**468 extras vs 26** in
`hierarchical-sweepga` and **11** in `poasta-everywhere`). The
>200 bp band actually improved (**0 extras vs 4** in hierarchical),
which is consistent with wfmash's finer-grained mapping behavior on
large structural homology.

**Interpretation.** wfmash's finer-grained mapping anchors more
short-range homology than FastGA, but seqwish's transitive-closure
induction over those anchors produces a higher count of small (51–200
bp) per-traversal disjoint segments. The 8 large top-level plans
(plan-1: 5.1 Mb total, plan-2: 2.9 Mb, plan-3: 3.3 Mb, plan-4: 1.2 Mb,
plan-7: 1.8 Mb) were all resolved with non-empty replacements, but the
*shape* of the replacement is more fragmented than what FastGA produced.

The 6 670-segment gap vs PGGB (13 288) is **not** explained by aligner
choice at the inner replacement step. PGGB does not just use wfmash
— it pipes wfmash output through **seqwish + smoothxg + gfaffix** on
the whole graph, not per-bubble. Specifically, **smoothxg** runs a sorted
block smoother across the entire pangenome, which collapses many of
the per-bubble disjoint segments that this experiment leaves behind.

### Specifically: do the 5–6 unreduced top-level bubble bands compact?

**No.** Visual inspection of `c4-crush-wfmash.png` against
`c4-crush-hierarchical.png` shows the same broad architecture: the
top-level macro-structure is similar (the 1102 level-0 POVU sites
indicate the same boundary topology), and the 5–6 unreduced bubble
bands in the hierarchical output remain as similarly-shaped bubble
bands in the wfmash output. The bp inside them is slightly different
(more 51–200 bp per-traversal duplicates from wfmash's finer
anchoring) but the bubble bands are **not collapsed**.

This is consistent with the prior `crush-sweepga-many-to-many-poasta-polish.md`
prediction (§Finding, candidate ranking) that **seqwish-induction
inside the bimodal plan-1 bubble** and **POASTA polish coverage
threshold** were the dominant remaining bottlenecks — not the aligner.

### What this experiment did rule out

- **Aligner anchor granularity is NOT the bottleneck** for the
  6 670-segment gap vs PGGB. wfmash's finer anchoring produces
  *more* 51–200 bp dup-extras, not fewer.
- **wfmash's wfmash-rs CLI compatibility** is environment-sensitive
  (vendored wfmash from wfmash-rs `d47b7e3` works; system wfmash
  v0.24.1-24-gf64becc4 does not). Future agents that need wfmash
  must use the rebuilt `target/release/wfmash` from the vendored
  source, not a `PATH`-resolved binary.

### What this experiment did NOT rule out

- **smoothxg on the whole graph** is the missing step. PGGB
  applies smoothxg after seqwish induction; this codebase does
  not have a smoothxg-on-whole-graph stage (only the per-bubble
  POASTA polish, which leaves outer per-bubble fragmentation
  intact).
- **POASTA polish coverage threshold** —
  `polish-max-median-traversal-len=1k` and `polish-max-traversal-len=10k`
  skip the bimodal plan-1 bubble's interior. Lifting these caps
  could potentially compress the residual stringy bubbles.
- **Per-bubble seqwish parameter scaling** — `seqwish-k=311` may
  not be the right anchor for wfmash output (which has different
  k-mer characteristics than FastGA output).

### Recommended next experiments

1. **`crush-lift-polish-caps`** — raise
   `polish-max-traversal-len` to 50k and
   `polish-max-median-traversal-len` to 10k so the polish loop also
   descends into the bimodal plan-1 bubble's interior.
2. **`crush-smoothxg-on-whole-graph`** — run smoothxg on the final
   `Ygs.gfa` (or before sort) as a post-crush polish step. This is
   the most likely-to-help follow-up since it directly mirrors PGGB.
3. **`crush-wfmash-low-pct-id`** — repeat this run with
   `map-pct-identity=70` (the original task description hinted at
   this). Currently the default is `ani50-2` ≈ 70 % already at the
   wfmash defaults, but it is not explicitly set in this run.
4. **`crush-wfmash-seqwish-k`** — try `seqwish-k` in {31, 51, 101, 311}
   on top of `method=wfmash` to see if a smaller seqwish-k extracts
   more compaction from wfmash's anchor distribution.

## Hard-gate checklist

- [x] Code change committed: new `ResolutionMethod::Wfmash` variant +
      `build_wfmash_seqwish_replacement` dispatch shim + `parse_name`
      alias + CLI/help updates — commit `1ee44fa`
- [x] `cargo build --release` succeeds
- [x] `cargo test --lib` 281 passed, 0 failed
- [x] `cargo test --bin impg` 57 passed, 0 failed (1 new
      `test_gfa_output_format_method_wfmash_routes_through_sweepga_with_wfmash_aligner`)
- [x] `cargo test --test test_crush_integration` — 1 pre-existing
      failure (`nested_bubble_level_descent_actually_descends`)
      annotated as known on `HEAD 471f089` in its own panic
      message; not introduced by this PR (also fails on prior
      `crush-sweepga-many` / agent-163)
- [x] Real C4 GRCh38 run completes, **465/465** paths preserved
- [x] Per-plan aligner choice + per-plan output shape reported (table above)
- [x] Final metrics: **20 476 S, 24 727 L, 597 070 bp, 468 51–200 bp
      dup-extras, 0 >200 bp dup-extras, 502 trivial-stringy**
- [x] PNG generated: `c4-crush-wfmash.png` at
      `/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T194829Z/c4-crush-wfmash.png`
      (canonical URL `https://hypervolu.me/~erik/impg/c4-crush-wfmash.png`
      — user-managed upload)
- [x] Specifically: report whether the 5–6 visually-distinct
      top-level bubble bands in `c4-crush-hierarchical.png` are
      reduced in `c4-crush-wfmash.png` — **NO, they are not
      reduced**; same boundary topology and similar bubble shapes
- [x] `docs/crush-wfmash-replacement.md` committed
- [x] `wg artifact crush-wfmash-replacement docs/crush-wfmash-replacement.md`

## Publish

The PNG is generated locally at
`/home/erikg/impg/data/c4_crush_wfmash_replacement_20260526T194829Z/c4-crush-wfmash.png`.
The canonical URL `https://hypervolu.me/~erik/impg/c4-crush-wfmash.png`
is the user's hypervolume publish path — this agent does not have
direct write access to that path; the user uploads PNGs from
`data/*.png` to `hypervolu.me` as part of their workflow (the
prior `crush-hierarchical.png` followed the same flow).

## Conclusion vs the original hypothesis

The task statement said:

> If wfmash drops segments toward PGGB's 13,288, that's strong
> evidence the aligner anchor granularity was the issue. If not,
> the bug is downstream in seqwish-induction / missing smoothxg
> on the whole graph.

**Result: NOT.** wfmash drove segments **up**, not down. The bug
is **downstream**: it is most likely the **missing smoothxg-on-whole-graph**
post-crush step (PGGB's distinguishing pipeline stage), plus the
**polish-cap restriction** that skips the bimodal plan-1 bubble's
interior. The aligner anchor granularity at the inner replacement
step was *not* the bottleneck.
