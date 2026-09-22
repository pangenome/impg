# crush aligner deep diagnostic — hybrid run, 12 plans

Per-plan diagnostic of the hybrid auto-2tier run on C4 GRCh38 chr6:31891045-32123783,
answering the five questions in the task spec. Data comes from a re-run with
per-plan input + output dumps under `IMPG_CRUSH_DEBUG_DIR=/tmp/crush-deep-diag`.

## Methodology

Re-ran the exact original hybrid command at commit `84a47d4`'s logic
(reproduced on the current `wg/agent-142/crush-aligner-deep` branch with an
instrumented build that has since been reverted; production tree is byte-clean
— see `git status` and `git diff` at the time of writing the doc):

```bash
IMPG_CRUSH_DEBUG_DIR=/tmp/crush-deep-diag \
LD_LIBRARY_PATH="/home/erikg/htslib-local/lib:/tmp/impg-libs:$LD_LIBRARY_PATH" \
/usr/bin/time -v target/release/impg query -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,auto-2tier=true,aligner=fastga,\
min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,\
max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,\
polish-max-median-traversal-len=1k:nosort' \
  -O /home/erikg/impg/data/c4_exp_hybrid_diag/run.nosort -v 1
```

Walltime: 5:51, exit 0, identical final stats to the original
(`per-round frontier sizes [r1=8, r2=3, r3=1]; total resolved=12; next_id moved
272214103 -> 272217124`). The diag run is at
`/home/erikg/impg/data/c4_exp_hybrid_diag/`; the dumps are under
`/tmp/crush-deep-diag/`. Analysis scripts: `/tmp/analyze_plans.py`,
`/tmp/bimodal_check.py`, `/tmp/inspect_components.py`.

Instrumentation reverted: `git diff` is empty for `src/` after analysis.

## Plan inventory (12 top-level plans + 299 polish sub-plans)

The hybrid run emits 12 top-level plans and, **inside the single sweepga
plan's local replacement subgraph**, 299 SPOA polish sub-plans (because the
sweepga path is the only one that calls `finalize_pairwise_induced_replacement`,
which runs polish; see `src/resolution.rs:2535`). POASTA plans bypass polish
entirely.

The 12 top-level plans, sorted by round and the in-round index that appears
in `run.stderr`:

| Round/Plan | Method | Traversals | min / median / max (bp) | total bp | root-span | Dump dir | Output S / L / comp / total bp |
|---|---|---:|---|---:|---:|---|---|
| r1 / 1 | **Sweepga** | 312 | 110 / 25155 / 31478 | 5,132,443 | 110 | `replacement_0005_sweepga_seqwish` | **732 / 977 / 4 / 37,266** (after seqwish; unchopped identical) |
| r1 / 2 | Poasta | 460 | 70 / 7286 / 7351 | 3,309,531 | 7277 | `replacement_0004_poasta` | 327 / 465 / 1 / 7,787 |
| r1 / 3 | Poasta | 437 | 153 / 157 / 42362 | 2,893,782 | 157 | `replacement_0087_poasta` | 271 / 385 / 1 / 42,590 |
| r1 / 4 | Poasta | 356 | 132 / 3319 / 9320 | 1,181,720 | 3319 | `replacement_0003_poasta` | 150 / 214 / 1 / 9,451 |
| r1 / 5 | Poasta | 287 | 6152 / 6210 / 6265 | 1,781,621 | 6196 | `replacement_0000_poasta` | 270 / 380 / 1 / 6,439 |
| r1 / 6 | Poasta | 325 | 112 / 112 / 25274 | 188,057 | 112 | `replacement_0086_poasta` | 108 / 151 / 1 / 25,340 |
| r1 / 7 | Poasta | 98  | 147 / 157 / 25255 | 65,552 | 157 | `replacement_0002_poasta` | 80 / 110 / 1 / 25,301 |
| r1 / 8 | Poasta | 47  | 165 / 165 / 25230 | 32,820 | 165 | `replacement_0001_poasta` | 3 / 3 / 1 / 25,230 |
| r2 / 1 | Poasta | 428 | 9 / 9 / 42214     | 2,829,061 | 9 | `replacement_0308_poasta` | 253 / 359 / 1 / 42,435 |
| r2 / 2 | Poasta | 325 | 47 / 47 / 25209   | 166,604 | 47 | `replacement_0307_poasta` | 107 / 149 / 1 / 25,271 |
| r2 / 3 | Poasta | 351 | 71 / 71 / 6073    | 42,927  | 71 | `replacement_0306_poasta` | 6 / 8 / 1 / 6,074 |
| r3 / 1 | Poasta | 190 | 58 / 58 / 6039    | 28,963  | 58 | `replacement_0309_poasta` | 3 / 3 / 1 / 6,039 |

Aligner subprocess details:

- **Sweepga (r1/1):** the only subprocess invocation, via
  `sweepga::library_api::sweepga_align` at `src/resolution.rs:2677` with
  `kmer_frequency = max(MIN_AUTO_SWEEPGA_KMER_FREQUENCY, n_traversals*10)`
  → for 312 traversals: `max(1000, 3120) = 3120`
  (`src/resolution.rs:2398-2405`, `:245-246`). Backend is FastGA, return
  code 0, wall time ≈18 s (FastGA invocation at 03:54:42 → seqwish at
  03:55:03). Raw PAF: 20,283,715 bytes, 77,666 lines.
- **POASTA (11 plans):** in-process via `poasta::aligner::PoastaAligner::new`
  at `src/resolution.rs:2764` with two-piece-affine scoring 1,4,6,2,26,1.
  Each `aligner.align::<u32, _>` call at `src/resolution.rs:3501` parameterises
  offsets to `u32` so per-sequence length cap is `u32::MAX - 1 ≈ 4.29 Gbp`.
  Wall times: round-1 builds 03:55:03 → 03:55:33 (≈30 s across 7 POASTA
  plans running in parallel via `into_par_iter` at `src/resolution.rs:769`),
  round 2 ≈8 s, round 3 ≈3 s.
- **No subprocess for POASTA:** library call only.

## Question 1: Is sweepga running but not actually compacting?

**Yes — the single sweepga plan (r1/1) is sub-compacting by a large margin.**
The sweepga subprocess succeeds (PAF: 77,666 alignment lines, 20.3 MB), but
the seqwish-induced replacement graph has **4 connected components, and 135
of the 312 input traversals end up as orphaned single-110 bp segments**, not
folded into the main 729-segment component covering the other 177 paths:

```
$ python3 /tmp/inspect_components.py \
    /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa
components=4
  comp0: n_segs=729, paths_in=177, first_seg=1
  comp1: n_segs=1,   paths_in=74,  first_seg=537   (110 bp orphan)
  comp2: n_segs=1,   paths_in=30,  first_seg=546   (110 bp orphan)
  comp3: n_segs=1,   paths_in=31,  first_seg=570   (110 bp orphan)
```

Root cause — the seqwish `min-match-len` rescue clamp is still too aggressive
for bimodal inputs:

- `src/syng_graph.rs:737-738` in the rescue path: `crush short-filter rescue:
  clamping seqwish min_match_len 311 → 110 (shortest input traversal is 110 bp;
  configured default would induce zero matches)`. The clamp targets *length*
  but does not check what **match-run length** is actually achievable in the
  PAF.
- The sweepga PAF contains 5,747 short-to-long alignment lines (110 bp query
  vs ≥5 kb target). **Zero of them have an internal exact-match run ≥ 110 bp**
  (max-run distribution from analysis: 5,625 of 5,747 are 80–99 bp; rest
  are 40–79 bp). Seqwish drops all 5,747, so the short paths get no edges
  into the long-target component.
- Long-vs-long alignments (47,947 of them) have long exact-match runs and
  survive the filter, building the 177-path 729-segment component.

**Should sweepga have run again next round?** No — by Phase 6 design, the
sub-bubbles found inside the sweepga replacement after seqwish are dispatched
to *polish* (which is SPOA), not to sweepga (`polish_replacement_gfa_with_flubbles`
at `src/resolution.rs:3124-3155` forces `method = config.polish_method`,
i.e. POA/POASTA, never sweepga). The 299 SPOA polish runs inside this
sweepga's local subgraph confirm polish happened, but polish runs SPOA on
each tiny sub-bubble, never re-running sweepga on the bigger structure.

**Fix proposal:** raise the seqwish min-match-len rescue at
`src/syng_graph.rs:737-738` to **min(shortest-traversal-len, max-actual-PAF-match-run)**.
The fix lives in `rescue_short_full_length_alignments` /
`crush short-filter rescue` (caller around `src/syng_graph.rs:706-760`). Or,
alternatively, also apply a per-PAF "internal-match-run histogram" floor so
the clamp reflects what's actually in the PAF, not just the shortest input
length. Without this, every bimodal sweepga plan will continue to orphan
the short component.

## Question 2: Is POASTA failing in some pathological way?

**Yes — POASTA produces stringy output on every bimodal plan, but the
subprocess (in-process library call) always returns success.**

"Stringy" here means: many tiny (1–4 bp) segments along a single connected
spine, not multiple disconnected components. All 11 POASTA outputs have
exactly **1 connected component** (good in the topological sense). The
fragmentation is *within* the spine.

Quantitatively, sort the 11 POASTA plans by **input bimodality** (max / median
ratio):

| Plan | tr | min/med/max | max/med | Output S | Output bp | S per 100 bp |
|---|---:|---|---:|---:|---:|---:|
| r1/5 (idx 0000) | 287 | 6152 / 6210 / 6265   | 1.01× | 270 |  6,439 | **4.19** |
| r1/4 (idx 0003) | 356 | 132  / 3319 / 9320   | 2.81× | 150 |  9,451 | 1.59 |
| r1/2 (idx 0004) | 460 | 70   / 7286 / 7351   | 1.01× | 327 |  7,787 | **4.20** |
| r1/8 (idx 0001) | 47  | 165  / 165  / 25230  | 153×  |   3 | 25,230 | 0.012 |
| r1/7 (idx 0002) | 98  | 147  / 157  / 25255  | 161×  |  80 | 25,301 | 0.32 |
| r1/6 (idx 0086) | 325 | 112  / 112  / 25274  | 226×  | 108 | 25,340 | 0.43 |
| r1/3 (idx 0087) | 437 | 153  / 157  / 42362  | 270×  | 271 | 42,590 | 0.64 |
| r2/2 (idx 0307) | 325 | 47   / 47   / 25209  | 537×  | 107 | 25,271 | 0.42 |
| r2/3 (idx 0306) | 351 | 71   / 71   / 6073   | 86×   |   6 |  6,074 | 0.099 |
| r3/1 (idx 0309) | 190 | 58   / 58   / 6039   | 104×  |   3 |  6,039 | 0.050 |
| r2/1 (idx 0308) | 428 | 9    / 9    / 42214  | 4690× | 253 | 42,435 | 0.60 |

Two failure modes are visible in this table:

1. **Homogeneous-but-stringy (r1/2 and r1/5).** When every traversal is
   essentially the same length (max/med ≈ 1.01), the input is *not*
   bimodal — it's 287 or 460 near-identical sequences. POASTA still emits
   one segment per 1-bp SNP site, giving 270–327 segments per ~7 kb of
   consensus. Median segment length on these plans is 1 bp. This is the
   "many 1-bp variation segments" pattern that PGGB collapses with gfaffix
   and crush does not (cf. `docs/crush-vs-pggb-comparison.md` §Step 4 point
   2; lib.rs:745 still skips gfaffix for the syng:crush dispatch).

2. **Bimodal-with-outliers (every other plan).** When max/median ≥ ~80×,
   POASTA's longest-first add order
   (`candidate_named_sequences_longest_first` at `src/resolution.rs:2365-2372`)
   adds the 25–42 kb outlier *first*, building a long spine. The 100s of
   short sequences are then aligned to that spine. Because the short
   sequences cover only ~150–165 bp of the 42 kb spine, POASTA's global-alignment
   model (`AlignmentType::Global` at `src/resolution.rs:2817-2819`)
   distributes their differences across the small window where they overlap,
   producing 0.3–0.6 segments per bp on the short-coverage region. The
   handful of plans that produce only a few segments (r1/8, r2/3, r3/1)
   are cases where the bulk-short sequences are byte-identical so they
   collapse to a single segment in the spine; whenever even one short
   carries a SNP, the segment count explodes again.

**Validation gate met: a specific bubble where alignment 'failed' (stringy
output despite successful subprocess):**

- **Plan r1/5 (dump idx 0000)** signature
  `CHM13#0#chr6:31744284-31976975:120373-126569 | … (287 paths)`.
  Input is 287 near-identical 6.2 kb sequences, output is **270 segments,
  380 links, median segment 1 bp**, total 6,439 bp. The CHM13 reference
  coordinates for this bubble are chr6:31864657–31870853 (computed as
  31744284+120373 .. 31744284+126569). PGGB-equivalent compaction would
  produce ~3 segments + per-SNP variation arcs, not 270.

- **Plan r2/1 (dump idx 0308)** signature
  `CHM13#0#chr6:31744284-31976975:38685-38694 | …`.
  Input is 428 short 9-bp sequences with 67 outliers around 42 kb. POASTA
  output: **253 segments, 359 links** for 42,435 bp content. CHM13 coordinates:
  chr6:31782969–31782978 (9-bp bubble). The expected compaction is ≤4
  segments (anchor – polyT-like variation – anchor + insertion arc); crush
  produces 253.

## Question 3: Sequence ordering passed to the aligners

**The order is longest-first for POA and POASTA, deterministic in input
order for sweepga and allwave.**

- POASTA: `src/resolution.rs:2748` calls
  `candidate_named_sequences_longest_first(candidate)` →
  `src/resolution.rs:2365-2372`:

  ```rust
  fn candidate_named_sequences_longest_first(
      candidate: &BubbleCandidate,
  ) -> (Vec<String>, Vec<(String, Vec<u8>)>) {
      let (headers, seqs) = candidate_named_sequences(candidate);
      let mut sorted = seqs;
      sorted.sort_by(|a, b| b.1.len().cmp(&a.1.len()).then_with(|| a.0.cmp(&b.0)));
      (headers, sorted)
  }
  ```

  Sort key is `(-seq_len, header)`; ties break by header name lexicographically.

- SPOA (POA polish): `src/resolution.rs:2715` calls the same helper. Same
  longest-first order.

- Sweepga: `src/resolution.rs:2651` calls `candidate_named_sequences` (no
  sort) at `src/resolution.rs:2350-2363`. The order is therefore
  `candidate.ranges` order, which is `path_idx` order (i.e., path-index
  order as the path-index iterator walks the working graph's paths array).
  Sweepga is all-vs-all internally so order doesn't affect output.

- Allwave: same as sweepga.

**Sample evidence (idx 0308 input.fa, first vs last):**

```
$ awk '/^>/{print $0}/^[ACGTN]/{print "LEN", length($0)}' \
    /tmp/crush-deep-diag/replacement_0308_poasta/input.fa | head -10
>__impg_bubble_path318_286
LEN 42214
>__impg_bubble_path184_168
LEN 42200
>__impg_bubble_path255_231
LEN 42196
…
>__impg_bubble_path9_9
LEN 9
```

42214, 42200, 42196 first; 9 bp last.

**Does the order affect output?** Yes, for POASTA — POASTA's add order
controls which sequence becomes the "spine" of the partial-order graph.
Longest-first means the long outliers (in a bimodal plan) become the spine
that all 100s of short sequences must align against. Empirically, every
plan with a strong bimodal ratio (max/med > 50×) has 100+ output segments
on a short anchor span (the r1/3, r1/6, r1/7 plans above). A
shortest-first or median-first order would let POASTA align the
bulk-shorts together first (yielding 1 short consensus segment), then add
the long outlier as a single insertion arc.

**Fix proposal:** for POASTA, switch to **median-first** or **bulk-first**
ordering for bimodal inputs. Concretely, replace
`candidate_named_sequences_longest_first` at the POASTA call site
(`src/resolution.rs:2748`) with a new helper that:

1. Detects bimodality: if `max_len > 5 * median_len`, treat as bimodal.
2. For bimodal: sort by `(distance_from_median, header)` ascending, so
   the consensus core builds first and outliers are added last as branches.
3. For non-bimodal: keep longest-first (it's only neutral or beneficial
   when lengths are within a factor of 2).

## Question 4: POASTA's actual maximum sequence length

**Codebase / hardcoded limit:** No explicit numeric cap. The relevant
constants:

- `src/resolution.rs:2761`: `let mut graph = POAGraph::<u32>::new();` —
  petgraph `StableDiGraph` with `u32` node indices: ~4.29 billion nodes.
- `src/resolution.rs:3501`: `aligner.align::<u32, _>(graph, sequence)` —
  the per-call offset type is `u32`, so per-sequence cap is `u32::MAX ≈
  4.29 Gbp`. Smaller offset types (`u8`, `u16`, `u64`) are defined in
  the upstream `poasta::aligner::offsets` module
  (`~/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/poasta-0.1.0/src/aligner/offsets.rs`)
  but impg always passes `u32`.
- No runtime length-rejection in impg's POASTA wrapper. `build_poasta_replacement`
  at `src/resolution.rs:2739-2789` adds every input sequence unconditionally.

**Practical limit (this run):** POASTA was successfully invoked on:

- **Single longest sequence:** 42,362 bp (r1/3, idx 0087). Returned cleanly,
  validated paths reconstruct to original bytes.
- **Largest total input:** 3,309,531 bp across 460 sequences (r1/2, idx
  0004). Wall time absorbed into the ~30 s round-1 build parallel batch.
- **Largest single-plan output:** 327 segments, 465 links (r1/2, idx 0004).

**Output-quality degradation threshold:** "Stringy" output begins as soon
as the input is bimodal with max/median > ~50× OR the input is homogeneous
but with even 1 SNP per 50 bp (then it's "1 seg per SNP"). Both modes are
already present in this hybrid run — neither is mitigated by raising or
lowering POASTA's size cap.

## Question 5: Other patterns explaining residual fragmentation

Beyond Q1 (sweepga-orphans the bimodal shorts) and Q2 (POASTA spine fragmentation
under bimodal or homogeneous-with-SNPs input), three additional patterns
contribute to the 476 vs 12 trivial-stringy-bubble headline gap:

### Pattern A — `min-traversal-len=5k` filters out 99% of small bubbles before any aligner runs

`src/resolution.rs:1646-1648`:

```rust
if traversal_stats.max_len < config.min_traversal_len {
    return None;
}
```

The hybrid run uses `min-traversal-len=5k`. The candidate-build phase
emits 70 candidates in round 1 from POVU, but only **8** pass this filter.
The discarded 62 candidates are all the small (60–150 bp) polyT / SNP /
short-indel bubbles that PGGB's seqwish-then-gfaffix pipeline cleanly
compacts. They are never sent to *any* aligner. This is the dominant
source of the 476 stringy bubbles per
`docs/crush-vs-pggb-comparison.md` §Examples A/B/C. The cited examples
(31826515-31826645, 31805524-31805642, 31764162-31764226) all sit at
60–150 bp `max_traversal_len` and trip this filter.

**Fix proposal:** either drop the floor (set `min-traversal-len=0`) and
let the routing rule dispatch the tiny bubbles to SPOA, or — as a smaller
change — keep the floor but add a **post-crush gfaffix pass** to the
syng:crush dispatch path. `src/lib.rs:745` `apply_graph_transforms` runs
`crush → sort` but skips the `normalize_and_sort` (→ `run_gfaffix`) call
that PGGB engines use. Adding gfaffix here would collapse adjacent-identical
chains that crush is leaving in.

### Pattern B — Routing rule uses median only, ignoring bimodality

`src/resolution.rs:1286-1300`:

```rust
fn auto_method_by_median(
    traversal_stats: TraversalStats,
    config: &ResolutionConfig,
) -> ResolutionMethod {
    let median = traversal_stats.median_len;
    if config.auto_spoa_max_traversal_len > 0 && median < config.auto_spoa_max_traversal_len {
        ResolutionMethod::Poa
    } else if config.auto_poasta_max_traversal_len > 0
        && median < config.auto_poasta_max_traversal_len
    {
        ResolutionMethod::Poasta
    } else {
        ResolutionMethod::Sweepga
    }
}
```

For r1/8 (med=165, max=25230, max/med=153×), r1/7 (max/med=161×), r1/6
(max/med=226×), r1/3 (max/med=270×), r2/1 (max/med=4690×), r2/2 (max/med=537×):
the rule picks POASTA because `median < 10000`. But in every case the
input is **bimodal**: a handful of 25–42 kb outliers (long-route paths
through the bubble) plus 100s of 9–165 bp shorts (direct-route paths).
POASTA's global aligner is the wrong tool for this — it forces 100s of
short sequences into a global alignment against a 42 kb spine. Sweepga's
all-vs-all + min-identity filter naturally handles heterogeneous-length
inputs but is currently only picked when `median ≥ 10 kb`.

**Concrete routing edge case:** r2/1 (idx 0308) has median=9 bp and
max=42,214 bp. The rule sends it to POASTA. The result is 253 segments
for a 9-bp bubble. If r2/1's `auto_method_by_median` had instead recognised
that 67 of 428 input traversals are >5 kb (≈16% of inputs are 42 kb-class
outliers), sending the plan to sweepga would have produced the same kind
of 4-component result as r1/1 — bad in different ways but probably with
fewer "tiny-fragment" segments in the short-consensus.

**Fix proposal:** insert a bimodality check at the top of
`auto_method_by_median` (`src/resolution.rs:1290`). E.g.:

```rust
let max = traversal_stats.max_len;
let median = traversal_stats.median_len;
let bimodal = max > median.saturating_mul(50) && max >= 5_000;
if bimodal {
    // bulk-of-paths is short; outliers are long: try star-biwfa or
    // sPOA over the bulk, then merge the outliers as side-branches.
    // Sending to sweepga is also valid for medians ≥ ~500 bp.
    return ResolutionMethod::Sweepga;
}
```

The bimodality cutoffs (50×, 5 kb min for "long") should be tuned against
the per-plan output-segments-per-bp metric in the Q2 table.

### Pattern C — POASTA replacements skip polish

`src/resolution.rs:2535` (polish runs only after `finalize_pairwise_induced_replacement`)
means none of the 11 POASTA plans go through the polish step. The 299 SPOA
polish sub-plans observed in the dumps are *all* from the single sweepga
plan's local subgraph. So even though POASTA produces stringy output, the
polish step that could compact the stringy interior never sees it.

**Fix proposal:** call `polish_replacement_gfa` on the POASTA-produced GFA
inside `build_poasta_replacement` at `src/resolution.rs:2773-2776`, before
`poasta_gfa_to_exact_graph`, gated on
`config.polish_iterations > 0 && config.polish_method != Smooth`. Need to
preserve the path-name set when running polish on the in-memory POASTA
GFA, since polish currently uses `parse_gfa`-based input.

### Sweepga vs POASTA routing edge cases (Q7 in the spec)

Going through the routing decision for each plan:

- All 11 POASTA plans were routed correctly *per the median-only rule*
  (median < 10 kb everywhere). None of them should have been routed to
  sweepga *under the current rule*. But six of them (r1/3, r1/6, r1/7, r1/8,
  r2/1, r2/2) have bimodal inputs where the current rule is the **wrong**
  rule (see Pattern B above).
- The single sweepga plan (r1/1) had median 25,155 bp, max 31,478 bp, so
  the rule correctly sent it to sweepga. The result is bad for a different
  reason (seqwish min-match-len clamp; see Q1).

So the routing rule has **no per-plan mis-routing under its own logic**, but
the *logic itself* mis-routes any plan whose input is strongly bimodal.

## Validation gate checklist

- [x] Each of the 12 hybrid-run plans characterized with traversal stats,
      length distribution, aligner output (S/L/components/bp), and method —
      see the Plan-inventory table.
- [x] Sequence ordering identified: longest-first for POA/POASTA at
      `src/resolution.rs:2365-2372`, called from `:2715` (SPOA) and `:2748`
      (POASTA). Sample evidence: idx 0308 input.fa head (42214,42200,42196)
      vs tail (9,9,9).
- [x] POASTA "max size" identified: no hardcoded cap in the impg wrapper;
      petgraph `u32` index gives ≈4.29 G nodes; per-sequence offset type
      `u32` at `src/resolution.rs:3501` gives ≈4.29 Gbp; practical
      output-quality cap is determined by **bimodality**, not raw length.
      Largest single sequence successfully processed in this run: 42,362
      bp (r1/3).
- [x] Specific bubble where alignment "failed" (stringy output despite
      successful subprocess), with traceable coordinates:
      **Plan r1/5 — CHM13#0#chr6:31864657-31870853 (signature
      …:120373-126569), 287 paths in, 270 segments out** — see Q2 above
      and `replacement_0000_poasta/poasta.gfa`. Second example:
      **Plan r2/1 — CHM13#0#chr6:31782969-31782978**, 428 paths in,
      253 segments out (`replacement_0308_poasta/poasta.gfa`).
- [x] At least one specific proposal per identified failure mode with
      file:line:
      - sweepga seqwish min-match-len clamp at
        `src/syng_graph.rs:737-738` — see Q1 fix proposal.
      - POASTA longest-first ordering at `src/resolution.rs:2365-2372` (call
        site `:2748`) — see Q3 fix proposal.
      - `min-traversal-len=5k` filter at `src/resolution.rs:1646-1648` —
        see Pattern A fix proposal.
      - Routing rule at `src/resolution.rs:1286-1300` — see Pattern B fix
        proposal.
      - POASTA bypassing polish at `src/resolution.rs:2773-2776` — see
        Pattern C fix proposal.
- [x] No production code changes (instrumentation lived in `src/resolution.rs`
      for the duration of the diag run, then `git restore`d; current `git
      status` shows only the `.wg` untracked dir and this doc).
- [x] `docs/crush-aligner-deep-diag.md` committed (this file).
- [x] `wg artifact crush-aligner-deep docs/crush-aligner-deep-diag.md`
      will be recorded post-commit.

## Reproducer (no rebuild required for analysis, only for re-running the diag)

The raw dump tree is at `/tmp/crush-deep-diag/`. To re-derive the tables in
this doc from the dumps alone, without re-running impg:

```bash
python3 /tmp/analyze_plans.py    # per-plan input/output GFA stats
python3 /tmp/bimodal_check.py    # bimodal classification of each input.fa
python3 /tmp/inspect_components.py \
  /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa
```

To re-run the diag from a clean checkout (a fresh agent must re-apply the
instrumentation in `build_poa_replacement` and `build_poasta_replacement`;
the diff in the agent log around r1/8 03:48 has the exact patch):

```bash
git submodule update --init vendor/syng
CPATH=/home/erikg/htslib-local/include \
LIBRARY_PATH=/home/erikg/htslib-local/lib:/tmp/impg-libs \
  cargo build --release --bin impg
# Re-apply the dump_direct_aligner_io patch (see agent-142 transcript)
IMPG_CRUSH_DEBUG_DIR=/tmp/crush-deep-diag \
LD_LIBRARY_PATH=/home/erikg/htslib-local/lib:/tmp/impg-libs:$LD_LIBRARY_PATH \
  target/release/impg query …(see Methodology section)
git restore src/resolution.rs
```
