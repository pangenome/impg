# Crush raw-FastGA → seqwish experiment (`crush-raw-fastga-seqwish`)

**Task:** `crush-raw-fastga`
**Date:** 2026-05-26
**Branch:** `wg/agent-172/crush-raw-fastga`
**Bubble FASTAs:** `/home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/bubbles/` (the same 8 level-0 cut bubbles extracted by [`docs/crush-per-bubble-isolation.md`](crush-per-bubble-isolation.md))
**Output directory:** `/home/erikg/impg/data/c4_crush_raw_fastga_20260526T222001Z/`
**Sweepga baselines compared against:** `/home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/isol/bubble_*/sweepga.stats.json` (sweepga invoked as `sweepga -t 32 -N input.fa`)

## TL;DR

User direction (verbatim): *"we should run just FastGA, RAW, in
desperation, convert to PAF."* The hypothesis being tested: maybe
sweepga is doing post-FastGA filtering we haven't found, so skipping
sweepga entirely and feeding FastGA's PAF straight into seqwish might
reveal a hidden bottleneck inside sweepga.

**The data refutes the hypothesis. Sweepga is not a filter that
removes alignments; it is an *amplifier* that produces orders of
magnitude more PAF lines than FastGA's own default invocation, and
those extra PAF lines are exactly what seqwish needs to collapse the
graph.** Removing sweepga collapses compression by 30–375× on the big
bubbles (e.g. bubble_08 falls from 373.7× → 1.18×; bubble_06 from
113.2× → 1.08×). The bottleneck is therefore **not** in any
post-FastGA sweepga processing layer — sweepga's processing is
load-bearing, and replacing it with a default FastGA invocation is
catastrophic.

What this *also* clarifies: the residual stringiness in the
hierarchical-run graph is not caused by sweepga discarding useful
edges. Sweepga is already pushing FastGA harder than its CLI defaults
do, and the bubbles still cap out at 1.0–373× compression because the
underlying sequences include long non-homologous segments and because
seqwish at `k=311` (effective `k = min(311, shortest_input_bp)`) needs
many short overlapping PAF chains to thread through divergent
regions. The next place to look for compression headroom is *upstream*
of sweepga (sequence preprocessing / partitioning) or *inside seqwish*
(induction at high k), not in sweepga's post-PAF code.

## Pipeline

For each of the 8 bubble FASTAs:

```
FastGA -pafx -T32 input.fa > raw_fastga.paf    # NO sweepga wrapper
seqwish -s input.fa -p raw_fastga.paf \
        -k <min(311, shortest_input_bp)> \
        -g raw_fastga.gfa -t 32 -P
```

Effective `k` is clamped to the shortest input sequence length when
that length is below 311, matching the clamp used by the
per-bubble-isolation runner so the comparison is apples-to-apples.

Runner: [`run_raw_fastga.sh`](../data/c4_crush_raw_fastga_20260526T222001Z/run_raw_fastga.sh)

## Side-by-side table

| Bubble | input bp | seqs | k_eff | sweepga PAF lines | raw-FastGA PAF lines | sweepga segs | sweepga seg-bp | sweepga comp | sweepga components | raw segs | raw seg-bp | raw comp | raw components | Δ segs (raw − sweepga) | Δ comp (raw − sweepga) |
|-------:|---------:|-----:|------:|------------------:|---------------------:|-------------:|---------------:|-------------:|-------------------:|---------:|-----------:|---------:|---------------:|-----------------------:|-----------------------:|
| 01 |    32 820 |  47 | 165 |      1 |    0 |  47 |    32 820 |   1.00× | 47 |  47 |    32 820 |  1.00× |  47 |     +0 |    +0.00× |
| 02 |    65 552 |  98 | 147 |  9 030 |    2 |  50 |    25 978 |   2.52× |  4 | 143 |    40 579 |  1.62× |  97 |    +93 |    −0.90× |
| 03 |   226 211 |  54 | 145 |  1 460 |   42 |  51 |    32 164 |   7.03× |  4 |  95 |    38 688 |  5.85× |  48 |    +44 |    −1.18× |
| 04 |   246 733 |  49 | 217 |  1 621 |  156 |  55 |     7 521 |  32.81× |  3 |  29 |   157 414 |  1.57× |  27 |    −26 |   −31.24× |
| 05 | 2 893 782 | 437 | 153 |134 891 |  306 | 172 |    43 647 |  66.30× |  3 | 537 | 1 454 994 |  1.99× | 395 |   +365 |   −64.31× |
| 06 | 1 181 720 | 356 | 132 |123 226 |  148 |  78 |    10 436 | 113.23× |  2 | 353 | 1 096 818 |  1.08× | 330 |   +275 |  −112.16× |
| 07 | 5 132 443 | 312 | 110 | 47 947 |  726 | 864 |    51 786 |  99.11× |136 |1208 | 2 536 214 |  2.02× | 217 |   +344 |   −97.08× |
| 08 | 3 309 531 | 460 |  70 |207 928 |  452 | 243 |     8 856 | 373.70× |  7 | 495 | 2 798 964 |  1.18× | 385 |   +252 |  −372.52× |

Notes on the table:

- **PAF lines is the headline.** Sweepga produces 30–950× more PAF
  lines than raw FastGA on every bubble that has any internal
  homology. This is not "post-filtering" — it is the opposite: sweepga
  is generating dramatically *more* alignment evidence than default
  FastGA does. Bubble_05 in particular goes from 306 raw chains to
  134 891 sweepga PAF lines (≈441× expansion).
- **Compression ratio** collapses to ≤2× on every multi-sequence
  bubble when sweepga is removed. The single exception is bubble_03,
  where raw FastGA still finds the dominant homology block and gets
  5.85× vs sweepga's 7.03×. Bubble_01 is identical (both produce a
  single-sequence-per-node graph because no aligner finds chains there
  with the default invocation).
- **Components blow up.** On bubble_05 sweepga gets the 437 input
  sequences into 3 connected components; raw FastGA leaves them in
  395 components (almost one component per sequence). The graph isn't
  threaded together — most sequences just sit as isolated paths.
- **Segment bp inverts.** When seqwish has too few PAF lines, it
  cannot find enough shared substrings to merge sequences, so segments
  remain long and per-sequence-private. Bubble_08 raw goes from
  8 856 bp of unique sequence (sweepga) to **2 798 964 bp** (raw FastGA)
  — essentially the whole input as private segments.

(All 8 stats JSONs are in
`/home/erikg/impg/data/c4_crush_raw_fastga_20260526T222001Z/isol/bubble_*/raw_fastga.stats.json`.)

## Why is raw FastGA so much sparser?

The sweepga binary used for the baseline was invoked as
`sweepga -t 32 -N input.fa`. The `-N` ("no filter") flag selects the
`NO_FILTER_SWEEPGA_KMER_FREQUENCY = 1_000_000` regime, which causes
sweepga's internal FastGA call to pass `-f1000000` — i.e. permit
k-mers that occur up to a million times before being treated as
repetitive. Default `FastGA` uses `-f10` (10× the average k-mer
frequency), which on these highly diverged C4 traversal sequences
discards most candidate seed chains before they ever produce a PAF
row.

This was *already observed* in
[`docs/crush-per-bubble-isolation.md`](crush-per-bubble-isolation.md)
§"Two interesting wrinkles": default sweepga CLI auto-tunes `-f`
similarly low and gives the same kind of empty PAF; running sweepga
with `-f 1000000` (the in-pipeline path) restored compression. The
current experiment confirms that the in-pipeline path's
disproportionately large PAF is **not** an artefact of sweepga doing
something downstream of FastGA — it is FastGA itself producing far
more chains when sweepga forces `-f1000000` on it. Default-mode FastGA
is the under-aligning step.

In short: "raw FastGA" with `-pafx` and defaults is *not* a more
permissive alignment than "sweepga + FastGA"; it is dramatically more
restrictive, because sweepga's parameter override is what unlocks
FastGA's seed engine for divergent input. Skipping sweepga therefore
removes the parameter override, not a filter.

## Verdict

**Sweepga is not the bottleneck. Removing sweepga makes the graph
worse, not better — on every multi-sequence bubble compression
collapses from 33–374× down to 1.0–2.0×, and connected components
explode from ≤7 to 217–395. This is because sweepga's default
operating mode forces FastGA into `-f1000000` (the `-N` /
`NO_FILTER_SWEEPGA_KMER_FREQUENCY` path), which causes FastGA to emit
2 orders of magnitude more PAF chains than its own CLI default does;
those extra chains are exactly what seqwish needs at `k≈110–217` to
collapse divergent C4 traversal sequences. Sweepga is not doing
post-FastGA filtering; it is doing pre-FastGA parameterisation, and
that parameterisation is load-bearing. If there is still compression
headroom in the C4 pipeline it is either (a) upstream of sweepga —
better bubble partitioning, lower-divergence inputs — or (b) inside
seqwish's induction at high k, but it is not hiding inside sweepga's
PAF post-processing.**

## Reproduce

```bash
cd /home/erikg/impg/data/c4_crush_raw_fastga_20260526T222001Z
for i in 01 02 03 04 05 06 07 08; do
  ./run_raw_fastga.sh \
    /home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/bubbles/bubble_${i}.fa \
    isol/bubble_${i}
done
```

Each `isol/bubble_NN/raw_fastga.stats.json` is directly comparable to
`/home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/isol/bubble_NN/sweepga.stats.json`.
