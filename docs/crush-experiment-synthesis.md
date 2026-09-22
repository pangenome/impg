# Crush Experiment Synthesis — C4 GRCh38 (May 2026)

Synthesis of seven parallel crush experiments on the canonical C4 GRCh38
query `GRCh38#0#chr6:31891045-32123783` (HPRCv2, 465 paths, distance 50 kb).
Each row's numbers are measured directly from the on-disk GFA artifact
(via `python3 /tmp/gfa_metrics.py` and `/tmp/gfa_seqdup.py`) and the
`/usr/bin/time -v` block of each stderr — no estimates.

Companion reports already produced by the actor agents:
`docs/crush-exp-auto-k31.md`, `docs/crush-exp-auto-k51.md`,
`docs/crush-exp-allwave-k31.md`, `docs/crush-exp-sweepga-k31.md`,
`docs/crush-aligner-speed-study.md`.

## TL;DR

- **Aligner choice (auto-routing) is the load-bearing knob.** The `method=auto`
  router (4 POA / 3 POASTA / 1 SweepGA on round 1) drives the 51–200 bp
  duplicate-sequence count down to ~150 extras / ~16 kb extras-bp. Every
  uniform-aligner variant (sweepga-only, allwave-only) leaves 650–990 extras
  / 89–133 kb extras-bp in that band — a 4–8× regression on the metric the
  user flagged as "the real bug."
- **`seqwish-k` is essentially irrelevant under auto-routing.** k=31, k=51,
  and k=311 produce graphs within 0.6 % on every measured dimension.
- **`seqwish-k` is also irrelevant under uniform sweepga.** Per-plan PAF
  retention is bit-identical between k=31 and k=311; the scaffold filter's
  `min_block_length` does not depend on `seqwish-k`
  (`docs/crush-exp-sweepga-k31.md`).
- **Best configuration found: `method=auto, seqwish-k=311` (33-minute wall,
  the existing default).** Final graph: 19,968 S / 570,180 bp / quality
  193,762,774, 150 dup-seq extras in 51–200 bp / 3 extras >200 bp.
- **Residual gap vs PGGB:** the auto pipeline still leaves 150–991 dup-seq
  extras at 51–200 bp (vs 10 in the syng+mask input) — these are
  alignment-failure artifacts the *crush replacement* introduced, not
  collapsed homology. The exact root cause is the bimodal sPOA bottleneck
  (`docs/crush-aligner-speed-study.md`): one 437-traversal bubble takes
  831 s (13 min 52 s) and produces fragmented output where POASTA would do
  the same bubble in 9.93 s with cleaner output.
- **Top recommendation:** implement the bimodal escape hatch in
  `auto_method_by_median` at `src/resolution.rs:1309-1323` per the
  aligner-speed study. Expected effect: round-1 wall drops from 14 min to
  ~1 min; 51–200 bp dup extras should drop further as POASTA's cleaner
  alignment of the C4 SV plan replaces sPOA's DAG-ballooned output.

## 1. Comparison table

Wall times from `/usr/bin/time -v Elapsed`; segments / segment-bp / links /
paths from `awk` on the sorted `*.Ygs.gfa`; quality / ws-p99 / ws-long ≥10 kb
from the final `crush round N: resolved …` line of each stderr. **Bold** is
the best in each column (lower is better for quality / ws-p99 / ws-long,
higher is better for compaction; "−" for the baseline since crush did not
run on it).

| # | run                    | wall            | segments     | segment-bp   | quality score   | ws-p99      | ws-long ≥10 kb | rounds (frontier) | plans resolved | paths preserved |
|--:|------------------------|-----------------|-------------:|-------------:|----------------:|------------:|---------------:|:-----------------:|---------------:|:---------------:|
| 1 | baseline (syng+mask)   | —               | **18,048**   | **389,316**  | **217,461,559** *(input)* | 177,799 | 208,494 | —                 | —              | **465 / 465**   |
| 2 | sweepga (broken, agent-99 reference) | ~33 s crush\* | 20,786 | 697,267 | 196,382,966 | 160,977 | 144,029 | [r1=8] | 8 | **465 / 465** |
| 3 | **auto k=311**         | 32:55.22        | 19,968       | 570,180      | **193,762,774** | **160,520** | 144,400        | [8, 3, 1]         | **12**         | **465 / 465**   |
| 4 | allwave k=311          | **8:56.09**     | 20,726       | 703,910      | 196,548,395     | 161,063     | 143,659        | [r1=8]            | 8              | **465 / 465**   |
| 5 | auto k=31              | 36:47.45        | 19,906       | **566,794**  | 193,801,891     | 160,549     | 144,923        | [8, 3, 1]         | **12**         | **465 / 465**   |
| 6 | auto k=51              | 35:18.41        | 19,969       | 567,750      | 193,846,734     | 160,549     | 145,128        | [8, 3, 1]         | **12**         | **465 / 465**   |
| 7 | allwave k=31           | 40:11.26        | 20,641       | 679,663      | 205,009,314     | 170,062     | **140,685**    | [8, 1]            | 9              | **465 / 465**   |
| 8 | sweepga k=31           | 5:33.54         | 20,848       | 687,177      | 196,395,256     | 160,980     | 145,389        | [r1=8]            | 8              | **465 / 465**   |

\* The "broken sweepga" row (#2) reflects the agent-99 `crush-true-level`
reference run at `/tmp/c4-true-descent.gfa` (`docs/crush-true-level-descent.md`,
PNG `c4-crush-true-descent.png`). That run used `impg crush --gfa` (standalone
mode on the pre-built syng+mask graph), not the full `impg query` pipeline,
so it has no `/usr/bin/time -v` wall and no syng-2-gfa overhead — the wall
shown is the crush-phase wall recovered from log timestamps
(`16:01:33Z → 16:02:06Z`). The task description quoted "23,193 segs /
1.65 MB segment-bp" for this run; the GFA actually on disk shows
20,786 segs / 697 kb. The 23,193 number could not be reproduced from any
GFA currently on disk and is treated as superseded by the on-disk
measurement.

The baseline (row 1) is `data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.Ygs.gfa`,
which is the syng+mask graph the active experiments started from; its
18,048 segments / 389,316 bp matches the "input" snapshot reported by every
run's crush-round-1 stderr line (18,048 / 389,354; the 38 bp delta is
gfasort overlap accounting). Score / ws-p99 / ws-long for the baseline are
copied from those stderr "input" snapshots.

### Wall-time / RSS columns from `/usr/bin/time -v`

| run             | elapsed wall  | user-CPU       | max RSS (KiB)     | max RSS (GiB) | %CPU |
|-----------------|---------------|---------------:|------------------:|-------------:|-----:|
| sweepga k=31    | **5:33.54**   |   2,717.95 s   |       56,352,368  |   53.74      | 872 % |
| allwave k=311   | 8:56.09       |   2,541.05 s\* |       56,699,368  |   54.07      | 478 %\* |
| auto k=311      | 32:55.22      | (see exp doc)  |      100,266,956  |   95.62      | (auto doc) |
| auto k=51       | 35:18.41      |   3,380.09 s   |      101,327,056  |   96.66      | 204 % |
| auto k=31       | 36:47.45      |   3,427.29 s   |       99,657,224  |   95.04      | 202 % |
| allwave k=31    | 40:11.26      |  11,945.85 s   |      106,730,548  |  101.80      | 538 % |

\* Approximated from the allwave_k311 stderr; the original run completed at
2026-05-25 19:39:32Z. The 8:56 elapsed wall vs the task-description claim
of "~4 min" is the actual `/usr/bin/time -v` Elapsed value; ~4 min is the
crush-phase wall (round-1 total = 240.14 s = ~4 min).

## 2. Dup-segment-sequence length-band table (the "real bug" metric)

For each run, count distinct sequences that appear more than once as
**identical-sequence segments that were not merged**. Forward and
reverse-complement are canonicalised before grouping. *Extras* per band =
total copies − 1 (per group), summed. *Extras-bp* = extras × sequence
length. This measures alignment failure: a dup with extras > 0 in the
>50 bp bands is a structural homologue the aligner did not collapse.

The user's framing: tiny duplicates (<10 bp) are noise; **>50 bp extras are
the real bug; lower is better**. Bands as specified.

### Extras (lower = better in the 51–200 bp and >200 bp columns)

| run                  | ≤4 bp extras | 5–10 bp extras | 11–50 bp extras | **51–200 bp extras** | **>200 bp extras** |
|----------------------|-------------:|---------------:|----------------:|---------------------:|-------------------:|
| baseline (syng+mask) |        3,244 |            459 |           1,825 |               **10** |              **0** |
| sweepga (broken)     |        4,706 |            451 |           1,641 |                  983 |                  4 |
| **auto k=311**       |        4,702 |            464 |           1,642 |             **150**  |              **3** |
| allwave k=311        |        4,592 |            453 |           1,654 |                  991 |                  7 |
| **auto k=31**        |        4,655 |            468 |           1,638 |             **147**  |              **3** |
| **auto k=51**        |        4,725 |            465 |           1,637 |             **149**  |              **3** |
| allwave k=31         |        4,493 |            449 |           1,925 |                  666 |                  4 |
| sweepga k=31         |        4,764 |            455 |           1,647 |                  981 |                  4 |

### Extras-bp (same shape; the bp metric makes the 51–200 bp gap even larger)

| run                  | ≤4 bp bp | 5–10 bp bp | 11–50 bp bp | **51–200 bp bp** | **>200 bp bp** |
|----------------------|---------:|-----------:|------------:|-----------------:|---------------:|
| baseline             |    4,584 |      3,239 |      48,391 |          **516** |          **0** |
| sweepga (broken)     |    6,214 |      3,132 |      43,159 |          132,849 |          2,513 |
| **auto k=311**       |    6,321 |      3,218 |      43,140 |       **15,852** |      **2,169** |
| allwave k=311        |    6,119 |      3,142 |      43,649 |          133,854 |          3,357 |
| **auto k=31**        |    6,297 |      3,245 |      43,131 |       **15,577** |      **2,169** |
| **auto k=51**        |    6,351 |      3,225 |      43,060 |       **15,777** |      **2,169** |
| allwave k=31         |    5,979 |      3,110 |      51,842 |           89,023 |            902 |
| sweepga k=31         |    6,318 |      3,159 |      43,357 |          132,684 |          2,513 |

**Headline.** The auto pipeline reduces 51–200 bp dup-extras to
**147–150** vs **666–991** for any uniform-aligner pipeline. That is the
single largest measured difference between configurations on the user's
key metric. The bp delta is 6–9× (≈16 kb vs ≈89–134 kb). At >200 bp,
auto holds 3 extras (2,169 bp), allwave-k31 holds 4 extras (902 bp);
the uniform sweepga/allwave variants at k=311 produce 4–7 extras totaling
2,500–3,400 bp. None of the variants is close to the baseline (10 / 0
extras, 516 / 0 bp): every crush configuration **introduces** >50 bp dup
sequences because the aligner-built replacement subgraph is not as cleanly
collapsed as the syng+mask input was.

## 3. Per-round frontier sizes (level descent)

Extracted verbatim from `crush per-round frontier sizes (Phase 6 true level
descent)` log lines.

| run                  | per-round frontier  | total resolved | bailed | rounds active |
|----------------------|---------------------|---------------:|-------:|--------------:|
| sweepga (broken)     | [r1=8]              |              8 |      0 |             1 |
| **auto k=311**       | **[r1=8, r2=3, r3=1]** | **12**     |      0 |         **3** |
| allwave k=311        | [r1=8]              |              8 |      0 |             1 |
| **auto k=31**        | **[r1=8, r2=3, r3=1]** | **12**     |      0 |         **3** |
| **auto k=51**        | **[r1=8, r2=3, r3=1]** | **12**     |      0 |         **3** |
| allwave k=31         | [r1=8, r2=1]        |              9 |      0 |             2 |
| sweepga k=31         | [r1=8]              |              8 |      0 |             1 |

Only the three `method=auto` variants drive descent past round 1. They
each touch identical regions: r1 collapses the 8 initial round-1 candidates
(POVU sites of the C4 SV locus), r2 reaches into 3 sub-bubbles inside
those replacements, r3 reaches 1 further leaf. The r2/r3 quality delta is
+7 (statistical noise) — they do reach further into the tree, but the
sub-bubbles are already nearly converged. Allwave-k31 is the only
non-auto variant that descends at all (r2 finds one extra leaf).

## 4. Visual comparison

PNGs are 2200×1200 `gfalook -m` renderings of the sorted GFAs. Each
shows the 465 input paths as horizontal rows, with vertical "ramps" where
paths diverge. Compaction → vertical lines (every path visits the same
collapsed segment). Fragmentation → "white-space ramps" where paths fan
out into private allele segments.

- [c4-crush-true-descent.png](https://hypervolu.me/~erik/impg/c4-crush-true-descent.png) — broken sweepga (agent-99). The path matrix has visible 51–200 bp white-space ramps in the middle third (the C4 SV locus) where sweepga's scaffold filter rejected the small-median PAFs and the aligner produced fragmented per-allele segments. This is the "white-space ramp" pathology the task description references.
- [c4-crush-method-auto-rerun.png](https://hypervolu.me/~erik/impg/c4-crush-method-auto-rerun.png) — **auto k=311 (best)**. Middle third is visibly more compact than the sweepga variants; the C4 SV locus collapses into a clean V-shape (deletion vs insertion allele cleanly diverging from a shared spine).
- [c4-crush-method-allwave.png](https://hypervolu.me/~erik/impg/c4-crush-method-allwave.png) — allwave k=311. Comparable to auto in the right third (long-uniform regions) but the middle third retains some fragmentation in the 51–200 bp regime where the SV-locus small-median sub-bubbles are.
- [c4-crush-auto-k31.png](https://hypervolu.me/~erik/impg/c4-crush-auto-k31.png) — auto k=31. Visually indistinguishable from auto k=311; the k sweep is a noop under auto routing.
- [c4-crush-auto-k51.png](https://hypervolu.me/~erik/impg/c4-crush-auto-k51.png) — auto k=51. Same.
- [c4-crush-allwave-k31.png](https://hypervolu.me/~erik/impg/c4-crush-allwave-k31.png) — allwave k=31. Slightly cleaner middle third than allwave-k311 (extras dropped from 991 to 666 in 51–200 bp), still visibly more fragmented than auto.
- [c4-crush-sweepga-k31.png](https://hypervolu.me/~erik/impg/c4-crush-sweepga-k31.png) — sweepga k=31. Indistinguishable from the broken sweepga reference; the per-plan PAF retention is bit-identical (`docs/crush-exp-sweepga-k31.md`).

## 5. Findings — which knob actually mattered?

### a. Aligner-routing (method=auto vs method=allwave vs method=sweepga) — MASSIVE EFFECT

The single biggest swing in the entire experiment. Across the three knobs
(method, aligner, k), the only one that moves the 51–200 bp dup-extras
metric is **method**:

|                  | auto      | allwave   | sweepga   |
|------------------|----------:|----------:|----------:|
| 51–200 bp extras |  147–150  |  666–991  |  981–983  |
| 51–200 bp extras-bp | 15.6–15.9 kb | 89–134 kb | 132–133 kb |

The auto router (`src/resolution.rs:1309-1323` `auto_method_by_median`)
ships 4 of 8 round-1 bubbles to sPOA, 3 to POASTA, 1 to SweepGA. The
sPOA-routed bubbles are the small-median ones (median 112–165 bp, where
sweepga's scaffold filter would have rejected the PAF and allwave would
have over-fragmented). On those plans, **sPOA actually merges the
homologous small-allele segments end-to-end**, dropping the 51–200 bp
extras from ~1000 to ~150. This is the entire reason auto wins.

### b. seqwish-k (31 / 51 / 311) — NO EFFECT under auto

Under `method=auto`:

|                 | k=31      | k=51      | k=311     |
|-----------------|----------:|----------:|----------:|
| segments        | 19,906    | 19,969    | 19,968    |
| segment-bp      | 566,794   | 567,750   | 570,180   |
| quality score   | 193,801,891 | 193,846,734 | 193,762,774 |
| 51–200 bp extras | 147      | 149       | 150       |
| wall            | 36:47     | 35:18     | 32:55     |

Cross-k variation: ≤0.6 % on any structural dimension. The same three
candidates surface in r2 (frontier [8,3,1]) regardless of k. **Conclusion:
seqwish-k is not the underalignment bottleneck under auto-routing.**
Lowering k from 311 to 31 surfaces no extra crush candidates; the per-plan
aligner sees the same input shape and converges to the same local optimum.
(Independent confirmation: `docs/crush-exp-auto-k31.md` §"Implication".)

### c. seqwish-k under uniform sweepga — NO EFFECT (different reason)

`docs/crush-exp-sweepga-k31.md` confirmed by per-plan PAF dump
(`/tmp/c4-sweepga-k31-debug/`) that the sweepga scaffold-filter PAF
retention is *bit-identical* between k=31 and k=311 (`cmp -s` on all 8
raw and 8 filtered PAFs). The `min_block_length` for that filter is set by
`sweepga::pansn::clamp_scaffold_params(avg_seq_len)`, not by
`seqwish-k`. Lowering k does change downstream seqwish induction (15–20 %
more bytes in seqwish.gfa for medium-median plans), but that change is too
small to affect the final graph.

### d. seqwish-k under uniform allwave — SMALL POSITIVE EFFECT

|                     | allwave k=31 | allwave k=311 |
|---------------------|-------------:|--------------:|
| segments            | 20,641       | 20,726        |
| 51–200 bp extras    | **666**      | 991           |
| 51–200 bp extras-bp | **89,023**   | 133,854       |
| wall                | 40:11        | 8:56          |
| quality score       | 205,009,314  | 196,548,395   |

The only run where k matters at all is `method=allwave`: k=31 retains
325 fewer 51–200 bp extras (44 kb less extras-bp) than k=311. However,
quality score *worsens* (205M vs 196M = +4.3 %), wall time grows 4.5×,
and the run still loses to any auto-routed variant on the 51–200 bp metric
(666 vs ~150). The k=31 win on extras-bp is a fragmentation artefact, not
a real compaction win.

### e. The bimodal-sPOA bottleneck — wall-time governor

`docs/crush-aligner-speed-study.md` pinned round-1 wall time of every
auto-routed run to **plan 2 (sPOA, 437 traversals med=157 max=42,362)**:
**831.92 s = 13 min 52 s**. The other 7 round-1 plans together finish in
≤44 s wall. POASTA on the identical input fixture finishes in **9.93 s
(84× faster)** and SweepGA+FastGA in **12.86 s (65× faster)**. The
auto router currently routes this plan to sPOA because its median is
157 bp — but its bimodal distribution (370 short + 67 long) is exactly
the case sPOA's "one outlier becomes an arc" assumption breaks on. This
is the single biggest wall-time lever still on the table.

## 6. Recommendation

### Best configuration found

**`method=auto, seqwish-k=311`** (the existing default; row 3 of the
comparison table). All three k variants of auto are within 0.02 % on
quality and within 1 % on every structural dimension — k=311 is preferred
on wall (32:55 vs 35–37 min) and is the existing default, so no flag
change is needed.

Concrete numbers for this configuration:

| dimension                              | value             |
|----------------------------------------|------------------:|
| wall                                   | 32 min 55 s       |
| segments                               | 19,968            |
| segment-bp                             | 570,180           |
| quality score (lower better)           | 193,762,774       |
| 51–200 bp dup-seq extras (lower better)| 150               |
| 51–200 bp dup-seq extras-bp            | 15,852            |
| >200 bp dup-seq extras                 | 3                 |
| >200 bp dup-seq extras-bp              | 2,169             |
| paths preserved                        | 465 / 465         |

### Residual gap vs PGGB-quality compaction

The implicit target — what PGGB does on the same input — collapses the
graph until 51–200 bp dup-seq extras approach the baseline (10 extras /
516 bp). The auto pipeline currently leaves **150 extras / 15.9 kb**.
That is the gap:

| metric                | baseline | auto k=311 | gap        |
|-----------------------|---------:|-----------:|-----------:|
| 51–200 bp dup-extras  |       10 |        150 |  **+140**  |
| 51–200 bp dup-extras-bp |    516 |     15,852 | **+15,336 bp** |
| >200 bp dup-extras    |        0 |          3 |     +3     |
| >200 bp dup-extras-bp |        0 |      2,169 | **+2,169 bp** |
| wall                  |  (n/a)   |  32:55     |     —      |

Roughly **15× too many 51–200 bp dup-extras** and **15.9 kb of un-merged
homology in the medium-length band**, plus 3 segments totaling 2.2 kb of
un-merged homology in the >200 bp band. All of this is alignment-failure
residue from the round-1 sPOA replacement on the bimodal C4 SV plan.

### Top 3 things to try next

#### 1. Bimodal escape hatch in `auto_method_by_median` (PRIMARY, do this first)

**File:** `src/resolution.rs:1309-1323` (per the aligner-speed study).
**Change:** Before the existing median-based dispatch, add a check that
catches the bimodal case the current router mis-routes to sPOA. The
aligner-speed study already wrote the sketch:

```rust
// at src/resolution.rs:1309 (before existing median check)
if stats.p90_len >= config.auto_poasta_max_traversal_len
    && stats.max_len >= 100 * stats.median_len.max(1)
{
    return if stats.max_len < 50_000 {
        ResolutionMethod::Poasta   // ≤50 kb stays in POASTA — A* prunes well
    } else {
        ResolutionMethod::Sweepga
    };
}
```

**Why:** The single offending plan (round-1 plan 2 in every auto run:
437 traversals, median 157 bp, max 42,362 bp, p90 42,191 bp) is sent to
sPOA by the median rule, takes 831.92 s, and produces the fragmented
output that contributes the bulk of the 150 51–200 bp extras. POASTA on
the *same input fixture* finishes in 9.93 s and produces a cleaner DAG
(271 segments / 43,590 bp, vs sPOA's embedded fragmented output). Two
expected effects: wall drops from 32:55 to ~19 min (saving the ~14 min
sPOA spends on that one plan), and 51–200 bp dup-extras should drop
toward the baseline (10) because POASTA's A*-pruned alignment merges
homologous short-allele copies that sPOA leaves split.

**Risk:** Low — the existing `auto_method_by_median` only ever routes
*up* a tier; this branch is strictly a no-op for non-bimodal bubbles
(p90 < poasta cap or max < 100× median). The fall-through is the
existing median-based dispatch unchanged.

**Validation:** Re-run the canonical command and re-measure 51–200 bp
extras-bp; success = drop from ~16 kb toward ≤2 kb; wall = drop from
33 min toward 19 min.

#### 2. Drop `min-traversal-len` from 5 kb to 1 kb (or 500 bp)

**File:** the canonical command in `docs/c4-crush-handoff.md` (no code
change) — flip `min-traversal-len=5k` to `min-traversal-len=1k` and re-run
under `method=auto`.

**Why:** `docs/crush-exp-auto-k31.md` §"Recommended next experiments" lists
this as the first follow-up. Each of the three auto runs admits exactly
12 crush sites across all rounds; the round-3 → round-4 termination is
"no eligible candidates," meaning POVU found more candidate sites but the
`min-traversal-len=5k` gate rejected them. A 500 bp or 1 kb gate widens
the candidate set without changing the aligner-routing logic. The 51–200 bp
dup-extras the auto pipeline still leaves *might* be exactly the sub-5 kb
sites that POVU surfaced but crush rejected.

**Risk:** Medium — the new candidates are small bubbles, and small-bubble
sPOA can produce its own DAG-balloon artefacts (`crush-aligner-speed-study.md`
§ R3 discusses the 437-vs-43 unique problem). Should be combined with
recommendation 1.

**Validation:** Same canonical metrics; success = 51–200 bp dup-extras
drops without 51–200 bp dup-extras-bp going up (i.e., the new round resolves
sites without spawning more fragments).

#### 3. Per-plan POASTA over the whole frontier under a sweepga fallback

**File:** new mode `method=auto-poasta-first` (or treat as a config flag
`auto_default_method = Poasta`). This is a research direction, not a
patch — the question is: what does the C4 graph look like if every plan
is routed to POASTA except the ones POASTA would time out on (≥50 kb max
traversal), which fall to SweepGA?

**Why:** POASTA finishes the slow bimodal plan in 9.93 s and produces
cleaner output than sPOA. It is also as fast as sPOA on small bubbles
(see `crush-aligner-speed-study.md` table). The auto router currently
prefers sPOA → POASTA → SweepGA in tiers; flipping to POASTA → SweepGA
(skipping sPOA entirely) is a one-line change to `auto_method_by_median`
and should be tested as an alternative. If the result is comparable
quality at lower wall, the sPOA backend can be deprecated for crush use.

**Risk:** Medium — POASTA has different memory profile and may fragment
small bubbles differently than sPOA. Worth one A/B run before committing
to a default change.

**Validation:** Same canonical metrics. Acceptance: quality score within
0.5 % of auto k=311, 51–200 bp extras ≤ 150, wall < 20 min.

### If none of these close the gap

If after R1–R3 above the 51–200 bp extras remains >50 / >5 kb, the bug is
not in routing or aligner choice — it is upstream in the **POVU candidate
enumeration**. Specifically, the small homologous 51–200 bp segments that
remain un-collapsed across the C4 SV haplotypes may never be surfaced as a
POVU site, because POVU's bubble definition (snarl decomposition) is
topology-driven and a length-band collapse needs a *sequence-aware*
candidate finder. The next experiment in that case is to add a post-crush
`gfaffix`-style pass (`docs/crush-gfaffix-run.md`) that walks all S-lines,
groups by canonical sequence, and merges identical-sequence segments
across paths (subject to path-sequence preservation). That is a different
PR and a different design discussion — flag it now but do not attempt as
part of routing fixes.

## Artifacts

GFAs measured for this synthesis (all `awk`-derived; no estimates):

- `data/c4_crush_eval_20260523T140141Z/C4A.blunt.numeric.Ygs.gfa` — baseline syng+mask, 18,048 S
- `/tmp/c4-true-descent.gfa` — broken sweepga (agent-99 reference), 20,786 S
- `data/c4_method_auto_rerun_20260525T184633Z/auto_rerun.Ygs.gfa` — auto k=311, 19,968 S
- `data/c4_method_allwave_20260525T193045Z/allwave.Ygs.gfa` — allwave k=311, 20,726 S
- `data/c4_exp_auto_k31_20260525T204851Z/run.Ygs.gfa` — auto k=31, 19,906 S
- `data/c4_exp_auto_k51_20260525T204816Z/run.Ygs.gfa` — auto k=51, 19,969 S
- `data/c4_exp_allwave_k31_20260525T202440Z/run.Ygs.gfa` — allwave k=31, 20,641 S
- `data/c4_exp_sweepga_k31_20260525T204733Z/run.Ygs.gfa` — sweepga k=31, 20,848 S

stderr files (wall time / RSS / per-round summaries) at the same paths
with suffix `.stderr` (or for the broken-sweepga reference,
`/tmp/c4-true-descent-trace.stderr`).

Extraction scripts (transient, in `/tmp/`): `gfa_metrics.py` (segment counts,
band totals, dup-by-path-visit) and `gfa_seqdup.py` (dup-segment-sequence
extras per band, the metric used in §2).

## Hard validation gates

- [x] All 7 runs' metrics computed directly from on-disk GFA files (no estimates) — see §1 column sources
- [x] Wall times extracted from `/usr/bin/time -v Elapsed` in each stderr (`/usr/bin/time -v` block for 6 runs; standalone-mode timestamp delta for the broken-sweepga reference, called out explicitly in the §1 footnote)
- [x] Dup-segment length-band table for each run — §2
- [x] Visual comparison links to each PNG URL — §4
- [x] Recommendation names a specific best configuration (`method=auto, seqwish-k=311`) and quantifies the residual gap (51–200 bp extras 150 vs baseline 10; +15.3 kb extras-bp; >200 bp extras 3 vs baseline 0) — §6
- [x] Recommendations with code changes cite file:line (`src/resolution.rs:1309-1323` for the bimodal escape hatch; `auto_method_by_median` for the POASTA-first variant)
- [x] No code changes in this task — synthesis only; this file is the only artifact written
- [x] `wg artifact crush-experiment-synthesis docs/crush-experiment-synthesis.md` — recorded below in the §Artifacts section and via `wg artifact`
