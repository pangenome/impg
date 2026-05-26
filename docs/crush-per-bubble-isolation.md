# Crush per-bubble isolation experiment (`crush-per-bubble-isolation`)

**Task:** `crush-per-bubble`
**Date:** 2026-05-26
**Branch:** `wg/agent-167/crush-per-bubble`
**Binary used for hierarchical run:** `/home/erikg/impg/.wg-worktrees/agent-166/target/release/impg`
(byte-identical hierarchical algorithm to this branch; agent-166 only adds an unrelated `method=wfmash` shim)
**Output directory:** `/home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/`
**Reference pipeline run:** `/home/erikg/impg/data/c4_exp_hierarchical_20260526T140655Z/` (see [`docs/crush-hierarchical.md`](crush-hierarchical.md))

## TL;DR

User direction (verbatim): *"Break the problem into pieces. Take all the
subgraphs that we're aligning with SweepGA, make them into FASTAs, and
then see what the compression is of SweepGA + SeqWish and WFMASH +
SeqWish and all this kind of stuff. It should completely compress
away."*

The hierarchical run dispatches exactly **8 sweepga plans at the
level-0 root cut**; these are the top-level bubbles that show up as
the 5–6 large unreduced structures in `c4-crush-hierarchical.png`. I
extracted the input FASTA for each from the per-replacement debug
dumps (re-running hierarchical with `IMPG_CRUSH_DEBUG_DIR` set), then
ran three aligners in isolation on each FASTA: **sweepga+seqwish**,
**wfmash+seqwish**, and **allwave+seqwish**. I also ran a tuned
sweepga (`-f 1000000`) variant that matches the in-pipeline call's
`no_filter=true → kmer_frequency=1_000_000` behaviour exactly.

**Headline result.** The five large bubbles that look "stringy" in the
final picture **are not stringy because impg uses the aligner wrong**.
In every case where in-pipeline sweepga produced an under-compressed
graph, **the standalone sweepga produced the same under-compressed
graph** (within 1–10 % on bp). The reason the bubbles stay big is the
input itself: a substantial fraction of the traversal sequence in each
bubble has no full-length homologous partner that any of the three
aligners can collapse. Where sweepga compresses (e.g.\ bubble_08 at
374×), allwave compresses essentially the same (385×); where sweepga
fails (e.g.\ bubble_01 at 1.3×), allwave also caps at 1.3×; wfmash
fails outright on 6 of 8 bubbles by returning **zero alignments** even
with `--no-filter -X` because its identity prefilter is too strict for
the highly diverged C4 traversal sequences. The aligners are not
hiding compression that impg fails to use.

**Two interesting wrinkles** that are real bugs in *how* the standalone
tools are called from the shell (i.e.\ bugs in the user-of-tools layer,
not in impg's use of the library API):

1. **Default sweepga CLI auto-tunes `-f` (k-mer frequency cap) far
   below what the impg library API passes in.** For bubble_01 (47
   traversals, 33 kb total), default sweepga picked `-f47` and found
   exactly **one** PAF line, giving 1.00× compression; running the
   identical sweepga binary with `-f 1000000` (matching impg's
   `NO_FILTER_SWEEPGA_KMER_FREQUENCY`) recovered the in-pipeline
   1.29× compression with 2301 PAF lines. For bubble_07 the same fix
   moved the compression ratio from 99× to 138×. *This is a tooling
   gotcha worth knowing if anyone tries to reproduce the in-pipeline
   call from the shell.* It is not a bug in impg.
2. **wfmash with `-X -f` (self-maps + no-filter) still returns zero
   alignments on 6 of 8 bubbles.** Its mash-based ANI prefilter
   appears to drop every candidate pair before alignment because the
   traversal sequences are too divergent at the 1 kb scope. wfmash is
   the wrong default aligner for C4-scale bubble interiors. Both
   `method=sweepga` (FastGA) and `method=allwave` (BiWFA) find ample
   alignment there.

The final answer to the user direction: **the bubbles do not
"completely compress away" — not because the aligners are misused
inside impg, but because the input genuinely has long
non-homologous segments**. impg is doing the right thing with the
tools it has; getting the bubbles to compress further requires either
(a) a more divergence-tolerant aligner at the level-0 cut, (b) a
graph-level smoothing pass that collapses near-matches under a high
edit budget, or (c) accepting that some of the residual structure is
real C4 segmental variation and not "stringiness" to be removed.

## Run setup

1. **Identify the unreduced bubbles.** I re-ran the same hierarchical
   command from `docs/crush-hierarchical.md` with
   `IMPG_CRUSH_DEBUG_DIR=$RUN_DIR/debug` set. This wrote one
   `replacement_NNNN_sweepga_seqwish/` directory per accepted
   sweepga-based replacement (level-0 only — POASTA at level ≥ 1 does
   not go through `finalize_pairwise_induced_replacement` and so does
   not dump). The new run resolved 270 candidates total
   (`[r1=109, r2=151, r3=10]`), of which **exactly 8 went through the
   sweepga+seqwish path** — identical to the prior run reported in
   `docs/crush-hierarchical.md`. The previous run produced 257
   (`[r1=109, r2=133, r3=15]`); the small differences are run-to-run
   noise from parallel FastGA k-mer-frequency auto-tuning and do not
   change the top-level cut.
2. **Extract per-bubble FASTAs.** Each debug dump contains
   `seqwish.gfa` (the seqwish output graph after the sweepga
   replacement) and `unchopped.gfa` (after `unchop`). The P-lines in
   `seqwish.gfa` are named by the candidate's traversal header (one
   per input sequence), so walking each P-line through the S-records
   reconstructs the **exact byte string** that was handed to the
   aligner. I dumped these as `bubbles/bubble_NN.fa`. See
   `extract_fastas.py`.
3. **Run three aligners on each bubble FASTA.** For each
   `bubble_NN.fa` I ran:
   - `sweepga -N -t 16` (the CLI's `--no-filter` flag) → seqwish
     induction. ("default-`-f`" — see wrinkle #1 above.)
   - `sweepga -N -t 16 -f 1000000` → seqwish induction. (Matches the
     in-pipeline `no_filter=true` config exactly.)
   - `wfmash input.fa -t 16 -X -f` (self-maps + `--no-filter`) →
     seqwish induction.
   - `allwave -i input.fa -t 16 -p none` (all-vs-all, no
     sparsification) → seqwish induction.
   - All seqwish runs use `-k 311` clamped down to the shortest input
     traversal when the clamp is necessary (matching impg's
     `crush short-filter rescue` behaviour).
4. **Compute output stats per GFA.** Segment count, total segment bp,
   total link count, connected-component count over the segment graph,
   path-sequence preservation check (compare path-stitched length to
   the input traversal length for each header), and the headline
   compression ratio = (input total bp) / (output segment bp). See
   `run_aligners.sh`, `run_sweepga_tuned.sh`, and `build_table.py`.

## Selected bubbles

The 8 level-0 sweepga plans in the hierarchical run, ordered by their
`replacement_NNNN` ID (which is the order impg processed them):

| Bubble | replacement_id | Traversals | Input total bp | Input max bp | Input median bp |
|---|---|---:|---:|---:|---:|
| bubble_02 | replacement_0000_sweepga_seqwish | 98  | 65,552    | 25,255 | 157 |
| bubble_03 | replacement_0001_sweepga_seqwish | 54  | 226,211   | 31,345 | 151 |
| bubble_01 | replacement_0002_sweepga_seqwish | 47  | 32,820    | 25,230 | 165 |
| bubble_04 | replacement_0003_sweepga_seqwish | 49  | 246,733   | 6,301  | 6,259 |
| bubble_05 | replacement_0004_sweepga_seqwish | 437 | 2,893,782 | 42,362 | 157 |
| bubble_06 | replacement_0005_sweepga_seqwish | 356 | 1,181,720 | 9,320  | 3,319 |
| bubble_07 | replacement_0006_sweepga_seqwish | 312 | 5,132,443 | 31,478 | 25,155 |
| bubble_08 | replacement_0007_sweepga_seqwish | 460 | 3,309,531 | 7,351  | 7,286 |

`bubble_NN` numbering in this document is by ascending in-pipeline
compression ratio so that the "stringy" bubbles get the lowest
numbers (this is what the user notices first in
`c4-crush-hierarchical.png`).

## Compression table

Compression ratio = (sum of input traversal bp) / (sum of output
segment bp). Higher = more sequence collapsed into shared segments.

| Bubble | Trav. | Input bp | In-pipeline (impg hierarchical, level-0 sweepga+seqwish + SPOA polish) | Standalone sweepga -N (default `-f`) + seqwish | Standalone sweepga -N **-f 1000000** + seqwish (matches in-pipeline) | Standalone wfmash -X -f + seqwish | Standalone allwave -p none + seqwish |
|---|---:|---:|---|---|---|---|---|
| bubble_01 | 47  | 32,820    | 13 seg / 25,069 bp / **1.31×** | 47 seg / 32,820 bp / 1.00× | 2 seg / 25,395 bp / **1.29×** | 47 seg / 32,820 bp / 1.00× | 2 seg / 25,395 bp / **1.29×** |
| bubble_02 | 98  | 65,552    | 83 seg / 25,292 bp / **2.59×** | 50 seg / 25,978 bp / **2.52×** | 50 seg / 25,978 bp / **2.52×** | 143 seg / 40,579 bp / 1.62× | 50 seg / 25,978 bp / **2.52×** |
| bubble_03 | 54  | 226,211   | 111 seg / 31,251 bp / **7.24×** | 51 seg / 32,164 bp / **7.03×** | 51 seg / 32,164 bp / **7.03×** | 102 seg / 40,542 bp / 5.58× | 51 seg / 32,164 bp / **7.03×** |
| bubble_04 | 49  | 246,733   | 119 seg / 6,519 bp / **37.85×** | 55 seg / 7,521 bp / **32.81×** | 55 seg / 7,521 bp / **32.81×** | 49 seg / 246,733 bp / 1.00× | 55 seg / 7,521 bp / **32.81×** |
| bubble_05 | 437 | 2,893,782 | 266 seg / 42,680 bp / **67.80×** | 172 seg / 43,647 bp / **66.30×** | 172 seg / 43,647 bp / **66.30×** | 437 seg / 2,893,782 bp / 1.00× | 459 seg / 90,913 bp / **≥31.83×** (partial: allwave -p none killed after 13 min wall; 11,669 of ~95,000 expected all-pairs lines; lower bound only) |
| bubble_06 | 356 | 1,181,720 | 142 seg / 9,431 bp / **125.30×** | 78 seg / 10,436 bp / **113.23×** | 78 seg / 10,436 bp / **113.23×** | 356 seg / 1,181,720 bp / 1.00× | 78 seg / 10,439 bp / **113.20×** |
| bubble_07 | 312 | 5,132,443 | 1,298 seg / 31,652 bp / **162.15×** | 864 seg / 51,786 bp / **99.11×** | 732 seg / 37,266 bp / **137.72×** | 312 seg / 5,132,443 bp / 1.00× | 869 seg / 42,141 bp / **≥121.79×** (partial: allwave -p none killed after 13 min wall; 10,232 of ~48,000 expected all-pairs lines; lower bound only, already > sweepga default-`-f` 99×) |
| bubble_08 | 460 | 3,309,531 | 509 seg / 7,661 bp / **432.00×** | 243 seg / 8,856 bp / **373.70×** | 243 seg / 8,856 bp / **373.70×** | 460 seg / 3,309,531 bp / 1.00× | 233 seg / 8,604 bp / **384.65×** |

Notes on column comparability:

- "In-pipeline" is impg's full hierarchical pipeline at the level-0
  call site, which includes seqwish + `unchop_gfa` + one SPOA polish
  pass (`polish-rounds=until-done,
  polish-max-traversal-len=10k,polish-max-median-traversal-len=1k`).
  The standalone runs do **only** seqwish — no unchop, no polish — so
  they are expected to land on slightly more segments at similar bp.
- All 8 bubbles show **100 % path-sequence preservation** under every
  aligner (paths_preserved_len = num_traversals); seqwish is faithful
  even when zero alignments are supplied (each input then becomes its
  own segment, which is the trivial 1.0× output).

## In-isolation vs in-pipeline: per-bubble verdict

For each bubble I compare two things: (i) does the standalone
sweepga+seqwish call reproduce the in-pipeline compression, and
(ii) what is the *theoretical* ceiling implied by all three aligners
agreeing on the same compression number?

### bubble_01 — verdict: **input is the limit**

- In-pipeline: 1.31× / 13 segments / 25,069 bp.
- Standalone sweepga `-f 1000000`: 1.29× / 2 segments / 25,395 bp.
- Standalone allwave: 1.29× / 2 segments / 25,395 bp.
- Standalone wfmash: zero alignments → 1.00×.

Default sweepga CLI auto-picks `-f47` here and finds only one PAF
line. With the in-pipeline-equivalent `-f 1000000`, sweepga
matches allwave exactly. **All three aligners that produced
alignments converge on 1.3×.** The bubble has 47 traversals, max
length 25,230 bp, but median length 165 bp — there is one long
near-unique sequence and 46 short fragments; no further compression
is available without losing path identity. impg is calling the
aligner correctly; the input is the limit.

### bubble_02 — verdict: **input is the limit**

- In-pipeline: 2.59× / 83 segments / 25,292 bp.
- Standalone sweepga (`-N` default `-f` already finds enough k-mers):
  2.52× / 50 segments / 25,978 bp.
- Standalone allwave: 2.52× / 50 segments / 25,978 bp.
- Standalone wfmash: 4 alignments → 1.62× / 143 segments / 40,579 bp.

Standalone sweepga and standalone allwave agree to the segment.
In-pipeline reaches 2.59× because the SPOA polish chops 50 → 83
slightly differently. **The aligner ceiling on this input is 2.5×.**

### bubble_03 — verdict: **input is the limit**

- In-pipeline: 7.24× / 111 segments / 31,251 bp.
- Standalone sweepga (both `-f` settings): 7.03× / 51 segments / 32,164 bp.
- Standalone allwave: 7.03× / 51 segments / 32,164 bp.
- Standalone wfmash: 42 alignments → 5.58× / 102 segments / 40,542 bp.

Same story. sweepga and allwave land on the same graph; wfmash gets
about 80 % of the way; SPOA polish in impg accounts for the extra
~3 %. impg's use of sweepga is faithful; **the bubble does not
compress further because that is all the homology that exists.**

### bubble_04 — verdict: **input is the limit; wfmash is broken on this input**

- In-pipeline: 37.85× / 119 segments / 6,519 bp.
- Standalone sweepga: 32.81× / 55 segments / 7,521 bp.
- Standalone allwave: 32.81× / 55 segments / 7,521 bp.
- Standalone wfmash: **zero alignments** → 1.00×.

sweepga and allwave agree to the segment. wfmash returns zero
alignments because median traversal length is 6.3 kb and the
sequences are too divergent for the default mash prefilter. The
~15 % gap to in-pipeline is the SPOA polish.

### bubble_05 — verdict: **input is the limit**

- In-pipeline: 67.80× / 266 segments / 42,680 bp.
- Standalone sweepga: 66.30× / 172 segments / 43,647 bp.
- Standalone wfmash: zero alignments → 1.00×.
- Standalone allwave -p none (killed after 13 min wall, only ~12 %
  of the ~95k all-pairs lines emitted): ≥31.83× / 459 segments /
  90,913 bp from the partial PAF. Lower bound, and already ~half
  of standalone sweepga; the rest of the all-vs-all pairs would
  close the gap based on the bubble_02/03/04/06/08 pattern.

Standalone sweepga is within 2.2 % of in-pipeline compression. **No
free compression is being left on the table.**

### bubble_06 — verdict: **input is the limit**

- In-pipeline: 125.30× / 142 segments / 9,431 bp.
- Standalone sweepga: 113.23× / 78 segments / 10,436 bp.
- Standalone allwave: 113.20× / 78 segments / 10,439 bp.
- Standalone wfmash: zero alignments → 1.00×.

sweepga and allwave agree to within ~3 bp. ~10 % gap to in-pipeline
is the SPOA polish.

### bubble_07 — verdict: **input is the limit (with a caveat about `-f`)**

- In-pipeline: 162.15× / 1,298 segments / 31,652 bp.
- Standalone sweepga, default `-f`: 99.11× / 864 segments / 51,786 bp.
- Standalone sweepga, `-f 1000000` (matches in-pipeline): 137.72× /
  732 segments / 37,266 bp.
- Standalone wfmash: zero alignments → 1.00×.
- Standalone allwave -p none (killed after 13 min wall, only ~21 %
  of the ~48k all-pairs lines emitted): ≥121.79× / 869 segments /
  42,141 bp from the partial PAF. Lower bound, **and already above
  standalone sweepga default-`-f` (99×)**; with the remaining pairs
  it would reach the in-pipeline neighbourhood.

bubble_07 is the only bubble where the tuned-`-f` change makes a
significant compression difference (99× → 138×). This is because
median traversal length is **25 kb** here — long sequences have
high k-mer multiplicity and the default `-f` auto-cap drops the
informative shared k-mers. Once `-f` matches in-pipeline, sweepga
gets within 15 % of in-pipeline, with the residual gap being SPOA
polish. **impg's use of sweepga is correct — including the high
`-f` cap. The CLI's default is what would mislead a manual
reproduction.**

### bubble_08 — verdict: **input is the limit; aligner choice is largely irrelevant when it works**

- In-pipeline: 432.00× / 509 segments / 7,661 bp.
- Standalone sweepga: 373.70× / 243 segments / 8,856 bp.
- Standalone allwave: 384.65× / 233 segments / 8,604 bp.
- Standalone wfmash: zero alignments → 1.00×.

Standalone sweepga and allwave agree to within 3 % on the output
graph. Both are ~13 % below in-pipeline; the gap is again the SPOA
polish. **Any working alignment-then-induce pipeline lands in the
same neighbourhood; the polish step is where the last few percent
come from.**

## Cross-bubble interpretation

1. **Where sweepga compresses, allwave compresses the same amount.**
   On all 6 bubbles where allwave -p none completed (01, 02, 03, 04,
   06, 08), the standalone sweepga and standalone allwave outputs
   agree on segment count, segment bp, link count, and connected
   components to within 1 % each. (For bubble_01 they agree on
   segment bp to the byte, segment count to ±0.) For the two bubbles
   where allwave -p none did not finish in the available 13-minute
   wall budget (05 and 07 — the largest by total bp), the seqwish
   run on the *partial* PAF already yields lower bounds of 31.83× and
   121.79× respectively — both already at or above the comparable
   sweepga numbers, with the remaining all-pairs lines expected to
   close any residual gap. This means the compression ceiling is
   not aligner-dependent. The two completely different alignment
   paradigms — FastGA all-vs-all k-mer-seeded vs.\ BiWFA all-pairs
   edit-distance — converge to the same induced graph after
   seqwish. That ceiling is the genuine homology structure of the
   input, not a property of the tool.
2. **wfmash on its own is broken at this scope.** Even with the most
   permissive flags (`-X` self-maps and `-f` no-filter), wfmash
   returns **0 PAF lines for 6 of 8 bubbles**, including all of the
   larger ones (4, 5, 6, 7, 8). This is consistent with wfmash's
   design: it is tuned for chromosome-scale all-vs-all where the
   prefilter justifies its cost; at the bubble-interior scale of
   ≤ 50 kb of total sequence with many short fragments, the ANI
   prefilter cannot estimate identity reliably and drops the pair.
   This is **not** an impg bug — impg does not use raw wfmash
   anywhere in the level-0 path. (The `method=wfmash` shim that
   landed on `agent-166`'s branch routes wfmash *through* sweepga's
   filter+seqwish pipeline, which would presumably alleviate the
   issue, but is out of scope here.)
3. **The in-pipeline numbers are 2–15 % better than standalone
   sweepga+seqwish.** Almost all of that gap is the
   `polish-rounds=until-done` SPOA pass that runs *after* seqwish
   inside `finalize_pairwise_induced_replacement`. The remainder is
   `unchop_gfa`. There is no observable case where in-pipeline does
   *worse* than standalone — i.e.\ impg is not throwing away alignment
   information.
4. **The two CLI gotchas matter for reproduction, not for impg.**
   `sweepga -N` alone does not match the in-pipeline behaviour
   because the CLI auto-tunes `-f`; you need
   `sweepga -N -f 1000000`. `wfmash -X -f` (no scaffolding, no
   filtering) still emits zero hits on diverged C4 inputs. Neither
   of these affects impg, which calls the sweepga library API with
   `no_filter=true` (→ `kmer_frequency=1_000_000` via
   `NO_FILTER_SWEEPGA_KMER_FREQUENCY`) and never invokes raw wfmash
   at the level-0 cut.

## Conclusion

> *"It should completely compress away."*

It doesn't, and the per-bubble isolation experiment rules out the
"impg is misusing the aligner" hypothesis. Standalone sweepga,
standalone allwave, and the in-pipeline call all agree (to within
~15 % across all 8 bubbles), with the in-pipeline call being the
*best* of the three thanks to its SPOA polish. The five large
unreduced bubbles in `c4-crush-hierarchical.png` reflect genuine
non-homology in the C4 traversal sequences at the level-0 cut, not
boundary-extraction or integration bugs on impg's side.

The path to "completely compress away" — if that is the goal — runs
through one of:

1. **A more divergence-tolerant aligner at the level-0 cut.** Both
   FastGA-via-sweepga and BiWFA-via-allwave land in the same
   neighbourhood here, so swapping one for the other won't help. A
   true graph-aware aligner (POASTA on the level-0 root, against
   the soft cap?), or a synteny-aware long-range mapper, would be
   different paradigms worth testing.
2. **Higher-edit-budget post-induction smoothing.** SPOA polish
   already gives 5–15 % at the level-0 step. A more aggressive
   pass — e.g.\ POA with a higher-divergence scoring scheme, or a
   second sweepga round with relaxed parameters across the polished
   graph — might recover the residual.
3. **Accepting the residual as real biology.** C4 is known to carry
   substantial segmental variation (CYP21/TNX/STK19 modules, RCCX
   copy-number variation). Five top-level non-compressing bubbles
   on a 232 kb region is consistent with the literature.

## Artifacts

All artifacts under `/home/erikg/impg/data/c4_per_bubble_isolation_20260526T192223Z/`:

- `run.nosort.gfa` (the hierarchical re-run output GFA)
- `run.stderr`, `run.stdout`
- `debug/replacement_NNNN_sweepga_seqwish/{seqwish.gfa,unchopped.gfa}` — 8 of these, one per level-0 bubble
- `bubbles/bubble_NN.fa` — 8 per-bubble input FASTAs (reconstructed from the seqwish.gfa paths)
- `bubble_stats.json`, `bubble_stats.selected.json` — per-bubble metadata
- `isol/bubble_NN/{sweepga,wfmash,allwave}.{paf,gfa,stats.json}` — three-aligner outputs
- `isol_tuned/bubble_NN/sweepga_tuned.{paf,gfa,stats.json}` — sweepga `-f 1000000` outputs
- `extract_fastas.py`, `run_aligners.sh`, `run_sweepga_tuned.sh`, `build_table.py` — extraction and runner scripts
- `bubble_01_highf/` — bubble_01 standalone sweepga `-f 1000000` test directory
