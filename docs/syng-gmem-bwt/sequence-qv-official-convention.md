# THE OFFICIAL SEQUENCE-QV CONVENTION (2026-10-10) — the owner's ruling, the likegt max-qv form, and the restatement of every committed QV table

The owner's ruling: the sequence QV is **quantity 1** — called-vs-truth sequence
accuracy — and its OFFICIAL convention is likegt's own `max-qv` form, surveyed
from likegt's implementation (`/home/erikg/likegt` `src/commands/max_qv.rs`, the
README's allwave aligner), not invented here. The committed `--truth-qv` test
mode's biWFA End2End injective-assignment computation (the yardstick of record
through the whole campaign) remains comparable BESIDE it; every committed QV
table is restated under the official convention below, with the per-window
deltas and their mechanisms named.

## 1. The convention, verbatim

* **Alignment**: allwave aligns the TRUTH haplotypes (the queries) against the
  CALLED haplotypes (the targets), all-vs-all, with likegt's own invocation —
  `allwave -i <fasta> -t 4 -p tree:5:0:0` (the max-qv command's default
  sparsification and thread count; allwave's OWN default dual-piece gap-affine
  scores 0,5,8,2,24,1, never overridden by likegt). Only truth-as-query lines
  enter the QV.
* **The PAF derivation** (likegt's `parse_paf_line`): `query_span =
  query_end − query_start`; `identity = matches / query_span` — the per-base
  error is measured in TRUTH sequence space **over the aligned span** (no
  forced End2End overhang); lines from one query to one target pool their
  (matches, span) sums first (likegt's per-target aggregation).
* **Per-truth-hap BEST** (likegt's most-similar-haplotype selection): for each
  truth haplotype the best alignment among the called homologs is taken —
  **no injective constraint**; both truth haps may take the same called
  homolog when it is the best for each.
* **QV** = −10·log10(1 − identity), clamped to [0, 50] — a perfect alignment
  carries **QV 50.0 exactly** (likegt's own cap; NOT the committed biWFA
  receipts' QV-60 convention for zero edits — the two caps differ BY
  CONVENTION, the zero-edit property does not). The diplotype sequence QV is
  the **MEAN** of the two per-truth-hap QVs (likegt's `avg_qv`).

## 2. The implementation (assessment-side, env-gated, additive)

`src/genome_inference/realign_score.rs` gains the official convention BESIDE
the committed computation, env-gated by `IMPG_TRUTH_QV_ALLWAVE=1`
(`IMPG_TRUTH_QV_ALLWAVE_BIN` overrides the allwave binary, default from PATH):
each truth-QV test record gains an additive `official_allwave` object — the
per-hap identities, QVs, best-called selections, the pooled per-pair
alignments, and the diplotype QV — while the committed biWFA fields stay
byte-identical. The env-unset default path is unchanged (proven by the
cut35-L3 control: the committed invocation with the gate unset reproduces
`calls.jsonl` and `calls.jsonl.truth-qv.jsonl` byte-identically). Three new
unit tests gate the PAF derivation, the per-hap-best pooling/selection, and the
convention end-to-end through the real allwave binary on a known diplotype
(`cargo test --release --lib realign_score`: 28 passed, 0 failed).

The yeast fleet's restatement instrument is
`genome/instrumented/realign-sequence-qv-allwave.py`: it imports the committed
`realign-sequence-qv.py` of record and reuses its panel fetch, fold criterion
and GFA/constituent containment gates VERBATIM — only the alignment and the
QV derivation are the official convention's — and joins the committed biWFA
per-locus receipts read-only for the per-window delta.

## 3. The restated tables (the gate: 30/30 CLI families PASS)

The CLI restatement re-ran every committed `--truth-qv` receipt family through
the SAME invocation VERBATIM plus the env gate, into separate `-official-qv`
out-dirs. The corrected comparability gate over all 30 families: **biWFA
fields 0 diffs everywhere; the product `calls.jsonl` byte-identical everywhere;
`official_allwave` absent only at the committed truth-pair-INEXPRESSIBLE loci**
(5 exempt: the pilot's 20kb locus 6 and 50kb segdup-block locus, L2's 20kb
locus, at their depths — the inexpressible class has no truth folds to align).
Every rank-1 call at every family and depth is **official-PERFECT** — both
truth haps' best alignments at identity 1.0, QV 50 — restating the biWFA
zero-edit QV-60 property one-for-one.

| family | n | rank-1 (unchanged) | official-perfect rank-1 | biWFA QV med → official QV med | misses' official QV med (min..max) |
|---|---:|---:|---:|---|---|
| pilot-deep | 22 | 4 | 4 | 26.67 → 35.65 | 34.11 (0.83..45.01) |
| L1-deep | 25 | 20 | 20 | 60.00 → 50.00 | 45.00 (35.28..50.00) |
| L2-deep | 25 | 22 | 22 | 60.00 → 50.00 | 45.17 (12.63..50.00) |
| L3-deep | 24 | 18 | 18 | 60.00 → 50.00 | 44.25 (42.00..50.00) |
| L4-deep | 25 | 17 | 17 | 60.00 → 50.00 | 47.50 (35.24..50.00) |
| cut35-pilot | 22 | 1 | 1 | 29.15 → 35.59 | 35.09 (0.83..50.00) |
| cut35-L1 | 25 | 12 | 12 | 43.03 → 50.00 | 44.99 (26.48..50.00) |
| cut35-L2 | 25 | 15 | 15 | 60.00 → 50.00 | 41.74 (12.63..45.01) |
| cut35-L3 | 24 | 10 | 10 | 43.02 → 48.25 | 41.86 (26.77..50.00) |
| cut35-L4 | 25 | 11 | 11 | 43.02 → 45.00 | 44.25 (33.29..50.00) |
| grain 20kb deep (pooled 37/60) | 60 | 37 | 37 | 60.00 → 50.00 at every locality | 31.93–46.51 per locality |
| grain 50kb deep (pooled 14/23) | 23 | 14 | 14 | 60.00 → 50.00 at every locality | 28.61–48.49 per locality |
| grain 20kb 35x (pooled 22/60) | 60 | 22 | 22 | — | 34.62–46.50 per locality |
| grain 50kb 35x (pooled 9/23) | 23 | 9 | 9 | — | 31.06–44.34 per locality |

(The full 30-row table with every per-family column is the generated artifact
`/home/erikg/hprcv2/qv-official-restatement/aggregate-official-tables.txt`.)

**THE YEAST FLEET RESTATED** (17 components, 1,132 truth-pair-expressible loci;
`realign-sequence-qv-allwave-tables.txt`): rank-1 **880/880 all
official-PERFECT at QV 50**; the 252 non-rank-1 misses keep their ranks and
classes and restate **ALL 252 HIGHER** (delta median **+13.82**, max +46.95;
official QV median 31.94, p10 18.33, min 10.49); **4/252 restate all the way
to official-perfect QV 50** — identical sequence, wrong class call (the
per-hap-best rule prices the truth sequence as perfectly present even when the
called CLASS differs); 115/252 of the misses take a NON-INJECTIVE per-hap-best
selection; 0/504 truth-hap sides carry a partial aligned span at yeast.

## 4. The near-twin QV34–43 band, recomputed — stated honestly

The campaign's named near-twin class (the committed biWFA QVs in [34,44) at
truth ranks 2–4 — the class the human ordinary-genome misses concentrate in,
"indistinguishable material at read level") restated under the official
convention across all 30 families: **51 windows; 47 restate HIGHER, 4 equal,
0 lower; the official band is QV 35.24–50.00 with median 45.00 — and 9 of
the 51 restate all the way to the QV-50 cap.** The biWFA [34,43] digits were
the INJECTIVE-ASSIGNMENT End2End price; the official convention prices the
same windows materially higher — i.e. **the near-twin misses are MORE
sequence-indistinguishable than the committed digits suggested**. The
mechanisms, named per window (the aggregate's table (2) carries all 51 rows):

* **The aligned-span effect (the dominant one)**: allwave measures the
  per-base error over the ALIGNED query span; the biWFA End2End computation
  forced the full span into the denominator. The near-twins differ from their
  best called homolog by a handful of indels flanked by unalignable
  overhangs; allwave trims those overhangs out of the rate, biWFA charged
  them. Per-window receipts: e.g. cut35-pilot locus 5, biWFA 34.57 →
  official 40.78 (both haps' aligned spans 10031/10035 of 10039 End2End
  columns); pilot-deep locus 6, 43.00 → 45.00.
* **The non-injective per-hap-best rule**: both truth haps may take the SAME
  called homolog. Where the near-twin called homolog is the best sequence
  match for BOTH truth haps, the official QV prices both at that alignment
  while the committed yardstick forced a no-replacement assignment across the
  crossed pair. Per-window receipts: L1-deep locus 2, biWFA 43.03 → official
  50.00 (both truth haps take called0 at identity 1.0 — sequence-perfect,
  class call still wrong); likewise L2-deep locus 0, L4-deep loci 13/14/24,
  cut35-L1 loci 2/9, cut35-L3 locus 10, cut35-pilot locus 1.
* **The indel-scoring model**: allwave's dual-piece gap-affine 0,5,8,2,24,1
  places long gaps differently than biWFA's single affine 0,4,6,2 — the
  remaining band rows where neither of the above applies (e.g.
  cut35-L4 locus 2, 40.00 → 43.50).

## 5. The verdict checks — does any campaign verdict change, or only the digits?

* **The yeast residue classes: UNCHANGED IN STRUCTURE.** Rank-1 880/1,132,
  every rank-1 official-perfect; the 252 misses keep every rank and class.
  Only the digits move — all of them UP. The honest nuance: the official
  convention prices the yeast misses at QV 10.5–50 (median 31.9), materially
  higher than the biWFA digits, and 4 windows are sequence-perfect under
  per-hap-best while still calling the wrong class.
* **The human headline numbers: UNCHANGED.** 77/99 generic-genome rank-1
  (all official-perfect), 4/22 segdup hotbed, 49/121 corrected-35x, every
  grain curve point — the rank columns are the gate-proven committed ones;
  the perfect calls restate QV 60 → QV 50 (identity 1.0), the same zero-edit
  property.
* **The 35x datasheet's QUAL band shares: UNCHANGED.** The product QUAL
  columns live in `calls.jsonl`, byte-identical at every family — the
  convention is assessment-side only and does not touch the product's own
  confidence reporting.
* **The grain verdict: STANDS.** The sweep's curve and every rank count are
  unchanged by construction; the QV columns restate (rank-1 QV 60 → 50, the
  misses' digits per the table above). The owner's grain hypothesis remains
  NOT SUPPORTED at either depth; the recommended grain remains the committed
  10kb.

## 6. Receipts and the honest execution notes

* **Receipts of record**: `/home/erikg/hprcv2/qv-official-restatement/` —
  `run-official-cli.sh` (the family runner), `run-all-families.sh` +
  `master-progress.txt` + `master.done` (the 9 non-grain families),
  `run-grain-families.sh` + `master2-progress.txt` + `master2.done` (the 20
  grain families after the fix), per-family `.log/.err/.rss/.wall` +
  `.biwfa-compare.txt` + `.done`, `aggregate-official.py` +
  `aggregate-official-tables.txt` (the generated artifact of record); the
  official out-dirs beside their committed counterparts
  (`cli-product-official-qv`, `cut3-*-cli-official-qv`,
  `grain-sweep/*-w*/cli-deep-official-qv`, `cli-cut35-official-qv` — compact
  receipts kept, bulk sidecars deleted per the disk-floor discipline);
  `/home/erikg/yeast/genome-balanced-diploid-validation-20260930/realign-sequence-qv-allwave-<C>.jsonl`
  (17 components) + `realign-sequence-qv-allwave-tables.txt`.
* **Two runner defects caught and fixed honestly**: (1) the first master's
  grain-family "exit=0" progress entries were FALSE — the progress echo's
  `$?` was clobbered by the `$(date)` command substitution — and the actual
  runs had all failed fast (`axis row count does not match the census window
  span`): the runner never overrode `--partition-graphs` for grain families
  (it fed the committed 10kb graphs with a grain census; the committed grain
  manifests prove the correct input is the grain's own
  `<L>-w<W>/locality-graphs`). Fixed; all 20 grain families re-ran clean.
  (2) The runner's comparability gate demanded `official_allwave` on EVERY
  record and so exit-1'd the five families carrying an inexpressible locus;
  the receipts were always good — the corrected gate (the aggregate's
  section (0)) exempts exactly the committed inexpressible records, and all
  30 families PASS.

Standing rules held: truth assessment-side only (the `--truth-qv` test mode);
no thresholds, no tuning constants (the convention is likegt's verbatim form,
stated); scoreboard machinery unmodified; the 64GiB guard armed in every run
(peaks < 6 GiB — the CLI restatement runs); the disk-floor discipline (bulk
sidecars deleted, 165G free held); NO PR PUSH (the owner pushes).
