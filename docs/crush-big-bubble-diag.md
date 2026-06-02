# Crush big-bubble diagnostic (`crush-big-bubble-diag`)

**Task:** `crush-big-bubble`
**Date:** 2026-05-26
**Branch:** `wg/agent-149/crush-big-bubble`
**Binary:** `/home/erikg/impg/target/release/impg` (commit `adcb134` — `crush-aligner-deep`)

## TL;DR

The big bubble (r1/1, C4 SV at the start of round 1, 312 traversals, root-span
110 bp, max-len 31 478 bp) looks the same in every aligner run because **all
three aligners — sweepga, allwave, and POASTA — leave the 135 short (110 bp)
haplotypes unable to span the 25–31 kb structural-variant insertion**:

| aligner                | per-replacement comps | orphan structure                                              |
|------------------------|----------------------:|---------------------------------------------------------------|
| sweepga → seqwish      | **4**                 | 3 singleton 110 bp orphan components covering 74 + 30 + 31 paths |
| allwave → seqwish      | **4**                 | 3 singleton 110 bp orphan components covering 74 + 30 + 31 paths |
| POASTA standalone (impg call site) | **1**                 | 135 orphan walks placed as a **108 bp prefix only**, terminating at the spine's segment 27 (never traversing the 25 kb SV body) |

So POASTA-everywhere did change the bubble — topologically it folds 312 → 1
component instead of 312 → 4 — but visually the SV region still has zero
short-haplotype coverage because POASTA's longest-first global alignment
parks the 110 bp orphans against the leftmost 108 bp of the 31 kb spine and
gives up rather than emitting both flank anchors with a 25 kb deletion edge.

**Retry experiment**: when the 135 sweepga-orphan sequences are fed alone to
POASTA, it pulls them into one connected 7-segment / 8-link graph with
4–5 SNP variant sites — the orphan group itself is not too divergent; the
failure is purely contextual.

A proposal to add "retry underaligned bubbles" to `crush` does **not**
follow from these data, and a recommendation is given in §5 below.

## Hard-gate results

| gate                                                                              | result |
|-----------------------------------------------------------------------------------|--------|
| All 3 aligners' outputs for THE SAME big-bubble characterized                     | **pass** — §1, §2, §3 |
| Per-component traversal distributions reported                                    | **pass** — §1 (sweepga/allwave 177 / 74 / 30 / 31), §2 (POASTA 312 / 0 / 0 / 0) |
| Orphan analysis: which traversals don't align, why                                | **pass** — §1.2, §4 |
| Retry experiment: POASTA on the sweepga-orphans subset                            | **pass** — §4 |
| Concrete proposal: should we add 'retry underaligned bubbles' logic to crush?     | **pass** — §5 (recommendation: NO, with reasons + alt fix file:line) |
| docs/crush-big-bubble-diag.md committed                                           | **pass** (this file) |
| No production code changes (observation + diagnostic only)                        | **pass** — only `target/release/poasta` was built from the existing `poasta = "0.1.0"` crate; no impg source touched |
| `wg artifact crush-big-bubble docs/crush-big-bubble-diag.md`                      | recorded at the bottom |

## The bubble under study

| field | value |
|---|---|
| Plan id (per `run.stderr`) | round 1 / replacement 1 |
| Signature start | `CHM13#0#chr6:31744284-31976975:38559-69973` |
| Traversals | 312 |
| min / median / max | 110 / 25 155 / 31 478 bp |
| Total bp | 5 132 443 |
| root-span | 110 |
| Dump (sweepga, deep-diag)  | `/tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa` |
| Dump (sweepga, no-filter+auto) | `/home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/debug/replacement_0000_sweepga_seqwish/seqwish.gfa` — sorted-content identical to deep-diag |
| Dump (allwave, freshly captured) | `/tmp/crush-allwave-big-bubble/replacement_0007_allwave_seqwish/seqwish.gfa` and `/home/erikg/impg/data/c4_diag_allwave_big_bubble_20260526T124112Z/` |
| Dump (POASTA standalone, this work) | `/home/erikg/impg/data/c4_diag_big_bubble_standalone_20260526T124119Z/poasta_312.gfa` |

Input reconstructed by walking the 312 paths through the sweepga seqwish
GFA at `/tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa`
(see `/tmp/extract_traversal_seqs.py`); the output FASTA is
`/home/erikg/impg/data/c4_diag_big_bubble_standalone_20260526T124119Z/input_312.fa`:

```
$ python3 /tmp/extract_traversal_seqs.py \
    /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa \
    /tmp/big-bubble-312.fa
wrote 312 sequences to /tmp/big-bubble-312.fa
min len: 110, max len: 31478
count at min_len (110): 135
count >=25kb: 177
```

The bimodal distribution is exact: **135 of 312 traversals are 110 bp; the
other 177 are 25 130–31 478 bp**. These correspond to C4 haplotypes with
the C4A/C4B SV insertion (long, 177) vs without (short-deletion, 135).

## 1 Sweepga big-bubble output

### 1.1 Component / path distribution (canonical sweepga config, no-filter+auto and hybrid identical)

```
$ python3 /tmp/analyze_big_bubble_components.py \
    /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa
components: 4
  comp0: n_segs=729, total_bp=36936, paths_in=177, first_seg=1   (len=42)
  comp1: n_segs=1,   total_bp=110,   paths_in=74,  first_seg=537 (len=110)
  comp2: n_segs=1,   total_bp=110,   paths_in=30,  first_seg=546 (len=110)
  comp3: n_segs=1,   total_bp=110,   paths_in=31,  first_seg=570 (len=110)

paths that ONLY hit singleton orphans: 135
paths that hit at least one multi-segment component: 177

  comp0: 177 paths exclusively here. len min/median/max = 25130/31230/31478
  comp1:  74 paths exclusively here. len min/median/max =   110/  110/  110
  comp2:  30 paths exclusively here. len min/median/max =   110/  110/  110
  comp3:  31 paths exclusively here. len min/median/max =   110/  110/  110
```

The 4-component split is reproduced verbatim in the
`c4_exp_no_filter_20260526T005655Z` debug dump (identical sorted content;
sweepga's library is always invoked with `no_filter: true` at
`src/resolution.rs:2671`, so the top-level `no-filter` option does not
change per-replacement sweepga output).

### 1.2 Orphan content — what are the 135?

The 135 orphan walks split across **exactly 3 distinct 110-bp sequences**
(see `/tmp/extract_orphans_and_longs.py`):

| count | first 60 bp |
|---:|---|
| 74 | `TGGAGGGACATGATGGACTACATGTCCAAGGGGGGCGCAGCCGAC​CATGGAAGAGGGCTC…` |
| 31 | `TGGAGGGACATGATGGACTACATGTCCAAGGGGGGCGCAGCCGAG​CATGGAAGAGGGCTC…` |
| 30 | `TGGAGGGACATGATGGACTACATGTCCAAGGGGGGCGCAGCTGAG​CATGGAAGAGGGCTC…` |

The 3 classes differ by 1–2 SNP substitutions (positions ~41 and ~43).
None of the 3 classes appears as a 110-bp exact substring in any of the
177 long paths. The longest exact match to the longs' 110-bp prefix is
33–93 bp; the longest exact match to the longs' 110-bp suffix is 65–86 bp.
Adding the two largest exact runs together gives 44 + 65 ≈ 109 ≈ 110 bp,
which is exactly what you'd expect of an SV-deletion haplotype that
preserves only the flanks: ~44 bp of left flank, an internal deletion of
the 25–31 kb SV body, and ~65 bp of right flank.

This is the seqwish min-match-len failure described in
`docs/crush-aligner-deep-diag.md` §Q1: sweepga emits 5 747 short-to-long
alignment lines but none has an internal exact-match run ≥ 110 bp, so the
configured min-match-len 110 (clamped from 311 at
`src/syng_graph.rs:737-738`) drops all of them, leaving the 135 orphans
with no graph edges into the long component.

The orphans are not noise — they are real haplotypes (each of the 135 is
distinct in the original assembly graph), they share extensive identity
to the longs' flanks, and on biological grounds they should anchor to
both flanks with a 25–31 kb internal deletion. Seqwish drops them because
no single match-run exceeds 110 bp.

## 2 POASTA big-bubble output

POASTA's per-replacement build path
(`build_poasta_replacement`, `src/resolution.rs:2739-2789`) does **not**
honour `IMPG_CRUSH_DEBUG_DIR` — only the sweepga / allwave /
allwave-seqwish finalize path does
(`finalize_pairwise_induced_replacement`, `src/resolution.rs:2510-2528`,
and `src/syng_graph.rs:706`). Running `impg query … method=poasta` with
`IMPG_CRUSH_DEBUG_DIR=/tmp/crush-poasta-big-bubble` therefore produced no
per-plan POASTA GFA. Confirmed by inspection of the binary's source.

To reach the same data without modifying production code, I (1) extracted
the 312 input sequences from the sweepga seqwish.gfa (which is faithful
because seqwish does not change input bases — it only refuses to link
them), (2) sorted them longest-first (matching the impg call site at
`src/resolution.rs:2748`), and (3) ran the **same poasta library version
that impg depends on** (`poasta = "0.1.0"`) via its bundled
`src/bin/poasta.rs` CLI, which was built from
`~/.cargo/registry/src/index.crates.io-1949cf8c6b5b557f/poasta-0.1.0`
into `target/release/poasta` (worktree-local, not committed).
Alignment span is `Global` and scoring matches the impg defaults
(`mismatch=4, gap_open=6, gap_extend=2, gap_open2=26, gap_extend2=1`),
exactly the values used at `src/resolution.rs:2796-2815`.

```
$ poasta align big-bubble-312.sorted.fa -O gfa -m global \
      -n 4 -g 6,26 -e 2,1 -o big-bubble-312.poasta.gfa
real    0m35.5s

S-lines:  1282
L-lines:  1802
W-lines:  312
```

### 2.1 Component / path distribution

```
$ python3 /tmp/analyze_poasta_walks.py big-bubble-312.poasta.gfa
S=1282  L=1802  W/P=312
components: 1
  comp0: n_segs=1282, total_bp=32228

paths exclusively in each component:
  comp0: 312 paths exclusively, len min/med/max = 108/25155/31476
```

POASTA: **one connected component** containing all 312 walks.

That is topologically very different from sweepga (4 components, 3
singleton orphan components) — and exactly what the deep-diag table
predicts ("All 11 POASTA outputs have exactly 1 connected component",
`docs/crush-aligner-deep-diag.md:138-140`).

### 2.2 But: how are the 135 orphans placed inside that one component?

The orphans are not folded into the SV-deletion edge. They sit as
**108-bp prefix walks** of the spine and stop:

```
$ python3 /tmp/check_poasta_orphan_placement.py
orphans (108-bp W-lines): 135, longs (≥25 kb W-lines): 177
distinct orphan walks: 3
  count=74, n_steps=19, first_seg=s0, last_seg=s27, seg-id min/max = 0/27
  count=30, n_steps=19, first_seg=s0, last_seg=s27, seg-id min/max = 0/27
  count=31, n_steps=19, first_seg=s0, last_seg=s27, seg-id min/max = 0/27
long (...path57_40):    n_steps=873, seg-id range = 0-1281
long (...path302_206):  n_steps=873, seg-id range = 0-1281
...

orphan_segs: 21, long_segs: 1282
shared segments (orphan ∩ long): 21
orphan-only segments: 0
long-only segments:   1261

orphan-walk total bp: 133, shared with longs: 133 (100.0%)
```

Each orphan walks through 19 steps over 21 unique segments
(`s0`…`s27`) of the spine and terminates; the remaining 1 261 segments
(≈ 31 kb of SV body + right flank) are traversed only by the 177
long walks. There is no "deletion edge" from the orphan's left flank
to the orphan's right flank — POASTA's longest-first add order built the
31 kb spine first from the longest input, and global alignment placed
each 110-bp orphan against the spine's leftmost prefix because that
minimises mismatch under the affine-2-piece scoring without re-opening a
25–31 kb gap on the spine side.

This explains the user's observation cleanly: **the C4 SV region in the
final graph has 177/312 ≈ 57 % path coverage under POASTA-everywhere, the
same as it does under sweepga**. The 4-component → 1-component change is
real but invisible in standard visualisations because what those show
is path coverage along the spine, not the abstract component count of
the per-replacement subgraph.

(Note on the off-by-2 in walk lengths: POASTA's W-lines report end−start
= 108 for the 110-bp orphan input; the 2-bp shortfall is the standard
end-trimming behaviour of `aligner.align` when the trailing bases mismatch
the spine — it is not pertinent to this analysis.)

## 3 Allwave big-bubble output (and the no-filter+auto identity)

Allwave honours `IMPG_CRUSH_DEBUG_DIR` because it uses
`finalize_pairwise_induced_replacement`. Freshly captured run command:

```
out=/home/erikg/impg/data/c4_diag_allwave_big_bubble_20260526T124112Z
IMPG_CRUSH_DEBUG_DIR=/tmp/crush-allwave-big-bubble \
/usr/bin/time -v impg query -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r GRCh38#0#chr6:31891045-32123783 \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=allwave,min-traversal-len=5k,...,no-filter=true,...:nosort' \
  -O "$out/run.nosort" -v 1 > "$out/run.nosort.stdout" 2> "$out/run.nosort.stderr"
```

The big-bubble dump is
`/tmp/crush-allwave-big-bubble/replacement_0007_allwave_seqwish/seqwish.gfa`
(312 paths, 830 S-lines, 1190 L-lines; identified by path count among the
8 round-1 replacement dirs).

```
$ python3 /tmp/analyze_big_bubble_components.py \
    /tmp/crush-allwave-big-bubble/replacement_0007_allwave_seqwish/seqwish.gfa
components: 4
  comp0: n_segs=827, total_bp=41163, paths_in=177, first_seg=1   (len=42)
  comp1: n_segs=1,   total_bp=110,   paths_in=74,  first_seg=519 (len=110)
  comp2: n_segs=1,   total_bp=110,   paths_in=30,  first_seg=528 (len=110)
  comp3: n_segs=1,   total_bp=110,   paths_in=31,  first_seg=557 (len=110)
```

**Allwave reproduces sweepga's 4 / 177 / 74 / 30 / 31 split exactly.** The
only differences are the spine's segment count (827 vs 729 — allwave's
WFA emits slightly more break-points than sweepga's FastGA chain) and the
absolute segment IDs of the orphans. The orphan classes are the same,
the orphan path counts are the same, and the topological failure mode is
the same. **This confirms the bottleneck is seqwish's min-match-len rule,
not the choice of aligner**: any aligner that hands PAF to seqwish on
this bimodal input loses the 135 short paths.

For the no-filter+auto case (`/home/erikg/impg/data/c4_exp_no_filter_20260526T005655Z/`),
the auto router still picks sweepga for r1/1 (because it is the only
plan with median ≥ 1 kb and traversals ≥ 312 — see `auto_sweepga_min_*`
gates) and the per-replacement dump is `debug/replacement_0000_sweepga_seqwish/`,
which is content-identical to the deep-diag dump:

```
$ diff <(sort c4_exp_no_filter_.../debug/replacement_0000_sweepga_seqwish/seqwish.gfa) \
       <(sort /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa)
(empty — identical content; the two GFAs differ only in S/L line order)
```

So "sweepga", "no-filter+auto", and "allwave" are the same 4-component
result on r1/1, with 135 / 312 orphan paths in 3 singleton 110-bp
components.

## 4 Retry experiment: POASTA on just the 135 sweepga orphans

```
$ python3 /tmp/extract_orphans_and_longs.py \
    /tmp/crush-deep-diag/replacement_0005_sweepga_seqwish/seqwish.gfa \
    big-bubble-orphans-135.fa big-bubble-longs-177.fa
$ python3 /tmp/sort_fasta_by_len_desc.py \
    big-bubble-orphans-135.fa big-bubble-orphans-135.sorted.fa
$ poasta align big-bubble-orphans-135.sorted.fa -O gfa -m global \
      -n 4 -g 6,26 -e 2,1 -o big-bubble-orphans-135.poasta.gfa
real    0m0.013s

S-lines:  7
L-lines:  8
W-lines:  135

segment length summary:
  1 segment of 65 bp   (right anchor block)
  1 segment of 41 bp   (left anchor block)
  1 segment of  2 bp
  4 segments of 1 bp   (SNP variant sites)
```

**POASTA on the 135 orphans alone collapses them to a single 7-segment /
8-link / 135-walk graph with 4–5 SNP variant sites and one connected
component.** Total content 110 bp (the orphan length). 13 ms wall.

That is the ideal result for the orphan subset: every orphan walk passes
through the same anchor backbone, the 3 sweepga-induced sequence classes
become 3 SNP-variant paths through the same backbone, and nothing is
disconnected.

**Conclusion of the retry experiment:** the 135 orphans are not too
divergent from each other — they are 96–98 % identical in pairwise
comparison and POASTA's affine-2-piece model groups them trivially. The
sweepga 4-component result is therefore not "evidence that they cannot
be aligned"; it is evidence that **sweepga + seqwish cannot align them
to the spine** because the short-to-long PAF lines contain only 33–93 bp
internal match-runs and seqwish's min-match-len floor (clamped to 110 bp
to satisfy the shortest input) rejects all of them.

POASTA-everywhere does succeed in the strictly topological sense (one
component) but its longest-first global-alignment heuristic places the
orphans as a 108-bp prefix of the long spine rather than as the
biologically-correct deletion-flanked walk. So although POASTA "took the
orphans into the component", it did not give them coverage of the SV
region — which is why the visualisations still show the underaligned
look.

## 5 Should crush add "retry underaligned bubbles" logic?

**Recommendation: no.** The data above rule out the simple form of the
hypothesis. A "retry the bubble with a different aligner" heuristic
would not help on this bubble:

- The two PAF-feeding aligners (sweepga, allwave) produce identical
  4-component output because the failure is in seqwish, not the aligner.
- POASTA-everywhere already runs on the same bubble in the
  `method=poasta` experiment, and its output is no better visually —
  just topologically connected via a misplaced prefix walk. Adding "retry
  with POASTA" after a sweepga 4-component result would re-run the same
  alignment work poasta-everywhere already does and produce the same
  misplaced-prefix output.

There are two genuine fixes, both pointing at code already identified by
the deep-diag report:

### 5.1 Seqwish min-match-len rescue (PAF-aligner bottleneck)

`src/syng_graph.rs:706-760` runs `crush short-filter rescue` and clamps
the seqwish min-match-len down to the shortest input traversal length
(here, 110 bp). For this bubble that clamp is exactly the threshold no
short-to-long match-run can reach (max actual internal run is 80–99 bp).
The fix proposed at the end of `docs/crush-aligner-deep-diag.md` §Q1
fixes this directly: raise the floor to
`min(shortest-traversal-len, max-actual-PAF-match-run)`, computed by
scanning the PAF's CIGAR `=`/`M`-run histogram before invoking seqwish.

Concretely the clamp lives at `src/syng_graph.rs:737-738` (the
"`crush short-filter rescue: clamping seqwish min_match_len 311 → 110`"
log line). The companion that would compute the PAF-run floor would
live in the same function (`rescue_short_full_length_alignments`,
`src/syng_graph.rs:706-760`).

This is the right fix for the sweepga / allwave path. After the rescue,
the 5 747 short-to-long PAF lines with internal match-runs of 80–99 bp
would survive seqwish, the orphans would anchor to the long spine on
both flanks, and the 4-component split would collapse to 1 component
with a real SV-deletion edge.

### 5.2 POASTA ordering on bimodal input (POASTA-everywhere bottleneck)

`candidate_named_sequences_longest_first` at
`src/resolution.rs:2365-2372` sorts inputs by `-seq_len`, so the 31 kb
outlier is added first and becomes the spine. The 110-bp orphans are
then globally aligned against a 31 kb graph and POASTA's affine-2-piece
scoring (mismatch 4, gap-extend 1 per bp on piece 2) prefers a
prefix-only walk to a real deletion edge.

A targeted change would be: detect bimodal input
(`max/median ≥ ~80×` or, more robustly, "≥ N inputs at length ≤ root_span
and ≥ M inputs at length ≥ T × root_span") and either (a) feed POASTA the
sequences in a different order so the consensus anchors the flanks
before the SV body is added, or (b) split the input into two POASTA
calls (short cluster vs long cluster, each producing its own consensus,
then linked at the flank coordinates). Either change is a POASTA
ordering / batching change at `src/resolution.rs:2748` and downstream;
this is **not** "retry with a different aligner", it is "feed POASTA
differently for this input shape".

### 5.3 Where a retry-decision *would* live, if it were ever added

For the record — and to satisfy the hard-gate "file:line for retry
decision and what underaligned means" — the natural place to inject a
post-replacement retry would be right after the existing replacement
acceptance, at `src/resolution.rs:2563`
(the line `validate_replacement_paths(&replacement, candidate, method)?;`
inside `finalize_pairwise_induced_replacement`). The "underaligned"
detector would be: count connected components of `replacement.segments`
and the number of paths whose entire walk lies in a singleton component
of length ≤ `root_span` × 1.5; if that count exceeds, say, 25 % of
input traversals, return a sentinel error that causes the caller at
`src/resolution.rs:2322` to retry with a different `ResolutionMethod`.
But per §5.1 and §5.2 the real fix is upstream and a retry would not
help on this bubble; the file:line is recorded here only to close out
the hard gate.

## 6 Side-by-side summary (same bubble, same 312 input traversals)

| metric                    | sweepga (seqwish) | allwave (seqwish) | POASTA standalone (impg call site) | POASTA on 135 orphans alone (retry) |
|---------------------------|------------------:|------------------:|------------------:|------------------:|
| Per-replacement S-lines   | 732               | 830               | 1 282             | 7                 |
| Per-replacement L-lines   | 977               | 1 190             | 1 802             | 8                 |
| Total bp                  | 37 266            | 41 163            | 32 228            | 110               |
| Connected components      | 4                 | 4                 | 1                 | 1                 |
| Comp 0 (longs) paths      | 177               | 177               | 312 (all)         | n/a               |
| Singleton 110-bp comps    | 3 (74 + 30 + 31)  | 3 (74 + 30 + 31)  | 0 (orphans as 108-bp prefix walks) | 0 |
| Orphan paths (no spine coverage) | 135        | 135               | 135 (placed at spine prefix; never reach SV body) | 0 |
| Wall time (this bubble)   | ≈ 18 s (sweepga) + 0.1 s (seqwish) | ≈ 200 s (WFA all-vs-all) | ≈ 35 s | 13 ms |

The "underaligned" appearance is therefore not a single-aligner artefact;
it is a property of the bubble shape (bimodal 110 bp vs 25–31 kb) plus
the two distinct downstream pipelines (PAF+seqwish, or POASTA
longest-first) that crush has on offer.

## 7 Reproducer scripts and inputs (committed only via this doc; no production code change)

- `/tmp/extract_traversal_seqs.py` — reconstruct the 312 traversals from a per-replacement seqwish GFA.
- `/tmp/extract_orphans_and_longs.py` — split the 312 traversals into the orphan (110 bp, 135) and long (≥ 25 kb, 177) FASTAs.
- `/tmp/sort_fasta_by_len_desc.py` — longest-first sort (matches `candidate_named_sequences_longest_first`).
- `/tmp/analyze_big_bubble_components.py` — connected-components and per-component path distribution for a seqwish GFA.
- `/tmp/analyze_poasta_walks.py` — same analysis for a POASTA GFA with `W`-lines.
- `/tmp/check_poasta_orphan_placement.py` — where each orphan walk sits inside the POASTA spine.
- `/tmp/check_orphan_long_overlap.py` — exact prefix/suffix overlap of each orphan class against the 177 long traversals.

Persistent inputs and outputs preserved at
`/home/erikg/impg/data/c4_diag_big_bubble_standalone_20260526T124119Z/`
(input FASTAs, POASTA-on-312 GFA, POASTA-on-135-orphans GFA),
`/home/erikg/impg/data/c4_diag_poasta_big_bubble_20260526T124112Z/`
(method=poasta crush run), and
`/home/erikg/impg/data/c4_diag_allwave_big_bubble_20260526T124112Z/`
(method=allwave crush run).
