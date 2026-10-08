# The diplotype-calling CLI product (the owner's product swap, 2026-10-07)

THE PRODUCT IS THE CLI. The realignment objective is no longer an internal
instrument: the inference CLI's `call-diplotypes` command IS the product —
per-locus diplotype calls with quality, end-to-end over the same proven
machinery (the marginal realignment likelihood, the territory-image extent
normalization, the alignment-induced locality domains, the DP/QUAL bounds of
`src/genome_inference/realign_score.rs` — the scoring layer of record, moved
from `examples/partition_realign_score.rs` into the library so the product
CLI and the instrument example share ONE implementation; the example is now a
thin wrapper preserving the committed surface byte-for-byte).

## The CLI surface

```
impg genome-infer call-diplotypes \
  --panel <syng prefix> \
  --routes <route graph dir> \
  --partition-graphs <locality map dir (window<N>.gfa.map.json)> \
  --component <axis lane, e.g. S288C#0#chrI> \
  --census <multi-matching census receipt (truth-free)> \
  --reads <sample FASTQ> \
  --derive-cache <record-derivation cache> \
  --loci <window ids, comma-separated, distinct and sorted> \
  [--selection locality|territory|default|env] \
  [--rss-budget-gib 64.0] \
  [--truth-qv <haplotype panel path> <haplotype panel path>]   # THE TEST MODE
  [--qv-biwfa <the pinned biWFA helper binary>]                 # test mode only
  --out-dir <NEW output directory>
```

* `--selection` defaults to `locality` in the product CLI — the
  alignment-induced locality domain (rule 9, implying the territory-image
  rule 8): the current best measured selection (880/1,132 = 77.7% truth
  rank-1 on the validation sample, per the committed fleet receipts). The
  instrument example keeps `env` as its default (the committed environment
  switches decide), so every committed runner reproduces byte-for-byte.
* The census input (`--census`) is the sample's own read records and
  placements against the panel — TRUTH-FREE (produced by the committed
  multicensus stage; it contains no truth material). The derive cache is a
  pure function of the FASTQ and the panel. No assessment-side shortcut
  enters the product path: no `--reference-pair`, no truth haplotypes, no
  scoreboard inputs.
* `--out-dir` must not exist; it receives `calls.jsonl` (THE PRODUCT
  OUTPUT), the instrument receipts beside it
  (`instrument-receipt.jsonl` + its sidecars — the checker/translation-gate
  artifacts, stream-written under the same 64 GiB RSS guard and the
  ordered-locus emission discipline), and `manifest.json` (the CLI's
  status envelope with the artifacts' fingerprints).

## THE DEFAULT PRODUCT OUTPUT — `calls.jsonl`

One JSON record per locus, truth-free (no truth key appears anywhere in the
record; verified at the translation gate):

```json
{
  "locus": 2, "partition": 2, "component": "S288C#0#chrMT",
  "model": "anchor-realign-v2-marginal-frame-territory-locality",
  "diplotype": {
    "fold_indices": [0, 102],
    "paths": [
      {"fold_index": 0,
       "members": [{"path_name": "S288C#0#chrMT", "start": 0, "end": 8578}, ...],
       "length": 8578, "territory_image": [0, 8578]},
      ...]
  },
  "class_count": 1081,            // the locus's candidate class universe
  "called_class_count": 1,        // exactly-tied winners included
  "called_classes": [[0, 102]],   // the class identity: unordered fold pairs
  "unit_count": 2543, "record_count": 436,
  "best_log_likelihood": -31884.7,
  "elsewhere_log_prob": -17.87,   // E: the candidate-independent elsewhere branch
  "alternative_log_gap": 212.33,  // winner LL minus the best non-excluded rival
  "qual": 922.14,                 // THE MODEL-INTERNAL QUAL (below)
  "qual_unbounded": true,
  "qual_spectrum_shape": "knee", "qual_knee_distance": 150.0,
  "qual_cluster_size": 3, "qual_cluster_k": 1,
  "qual_alternative_similarity": 6.1e-93,
  "qual_alternative_fold_pairs": [[0, 1]],
  "nearest_rival_distance": 75.0
}
```

* **The two paths**: `diplotype.paths` — the winner class's two folds, each
  spelled as its member rows (panel path name + interval) with its length
  and its territory image (the locus's own comparison coordinate). A
  homozygous call is `i == j` (the class pairs are unordered `i <= j`).
* **The class identity**: `called_classes` — the fold-index pairs of the
  called set (the winner plus any exactly-tied class; every tied class is
  co-winner by the equal-LL rule).
* **E/gap fields**: `elsewhere_log_prob` (the derived generated-elsewhere
  branch E) and `alternative_log_gap` (the winner-vs-best-rival separation,
  the QUAL's own input).

## THE MODEL-INTERNAL QUAL — the exact emitted form (derived)

The owner-ruled cluster semantics, verbatim from the committed machinery
(`spectrum_knee` / `cluster_form_qual` / `qual_p` / `qual_from_p`, the
mirrored pure functions of the instrument of record):

1. **k classes**: the QUAL spectrum is every non-winner class's
   *material distance* to the winner (the differing observed mass over the
   fold node/edge multisets) with its *relative likelihood*
   `e^{LL_c - LL_win}`; the knee cut (`spectrum_knee`) splits the winner's
   material neighborhood (distance `<= cut`: EXCLUDED — a near-identical
   spelling variant of the same call, not a rival hypothesis); `k` is the
   connected-component count of the winner plus every exactly-tied called
   class under the cut (union-find over the tie bands).
2. **The alternative log-gap**: `Δ = LL_win - LL_best-outside` over the
   NON-excluded classes (the nearest genuinely rival hypothesis).
3. **The two-hypothesis flat-prior Phred form**: the winner's normalized
   score `s_win = 1`; the posterior
   `p = s_win / (k·s_win + e^{-Δ}) = 1/(k + e^{-Δ})`;
   `QUAL = -10·log10(1 - p)`.
4. **No meaningless saturation**: where the f64 posterior is representable
   (`p < 1`), the emitted QUAL is the committed cluster value verbatim
   (field-identical to the receipts). Where the committed form's posterior
   rounded to exactly `1.0` (the receipts' `qual: null` +
   `qual_unbounded: true` — an f64 representation limit, NOT a model
   statement), the product evaluates the SAME form in log space:
   `QUAL = 10/ln(10) · [ln(k + e^{-Δ}) − ln(k − 1 + e^{-Δ})]`
   (at `k = 1`: `QUAL = 4.343·Δ + 10·log10(1 + e^{-Δ}) ≈ 4.343·Δ` Phred —
   growing linearly with the realignment separation, never saturating).
   `k > 1` saturates only at the honest cap `10·log10(k/(k−1))` of
   exactly-tied rivals — a real tie, honestly bounded.

## THE TEST MODE — `--truth-qv` (the owner's ruling, locked in)

The called-vs-truth sequence QV is a **VERY SPECIAL TESTING THING**. It
lives ONLY behind `--truth-qv <haplotype> <haplotype>`, which REQUIRES the
truth haplotypes as input (panel path names, e.g.
`S288C#0#chrMT SK1#0#chrMT`):

* The first haplotype must be the component lane — the window's own axis
  row is the truth-first fold (the committed convention, made explicit);
  the second is the input haplotype's fold under the committed path-name
  rule (single-window partitions; exactly one fold).
* The material per fold is the fold's own spelled sequence — the panel
  fetch of the member row `[start, end)` window (the committed
  convention-(a) material; under the locality domain the map rows carry
  the fold intervals, witnessed by the committed constituent gates).
* The alignment is the committed one, through the SAME pinned biWFA
  helper (`--qv-biwfa`, default
  `genome/instrumented/qv-biwfa/target/release/qv-biwfa`): gap-affine
  End2End, match 0 / mismatch 4 / gap-open 6 / gap-extend 2, Medium
  memory, no heuristic.
* The yardstick is the committed assignment rule: both orders of matching
  the two called homologs to the two truth homologs scored by the
  assignment-wide per-base error (summed edits / summed columns); lower
  rate wins, ties to the lower edit count, then the direct order;
  likegt's mean-identity rule computed and any disagreement reported.
* `error = (mismatches + gap columns) / total columns`,
  `QV = -10·log10(error)` with the QV-60 cap at `error <= 1e-6` (a
  perfect, zero-edit call). `perfect = (total edits == 0)`.
* The output is a SEPARATE artifact — `calls.jsonl.truth-qv.jsonl`, one
  record per locus (inexpressible truth pairs reported with
  `truth_pair_expressible: false` and a null QV) — never mixed into
  `calls.jsonl`, never produced without the flag. The default product
  output contains NO truth-referenced field of any kind.

## The translation gate (chrMT first, then the fleet)

The CLI product path and the committed instrument share one
implementation; the gate is measured, not assumed. **chrMT: CLOSED**
(all numbers from the 2026-10-07 session, the regenerated chain anchored
on the surviving receipts of record — the locality receipts and the
census scratch were cleaned from the validation dir, so the chrMT/chrI
censuses were regenerated through the committed multicensus runner and
the chrMT locality receipt regenerated through the committed runner,
then verified against the surviving per-locus QV receipts):

1. the moved library reproduces the committed chain — the regenerated
   chrMT locality receipt matches the surviving per-locus receipts of
   record (`realign-sequence-qv-locality-chrMT.jsonl` called folds +
   truth ranks, 9/9 expressible loci);
2. the CLI's instrument receipt is field-identical to the committed
   runner's receipt on the same binary (0 diffs outside the
   walls/rss timing fields; the one earlier 1-ULP divergence was
   a stale-binary artifact, resolved at the same build);
3. the product `calls.jsonl` is field-identical to the receipt per
   locus — winner folds, class identities, best LL, E, alternative
   log-gap all exact; the QUAL field-identical wherever the committed
   form is representable (1 locus) and the derived log-space extension
   exactly at the 13 f64-saturated loci (the one named, owner-mandated
   divergence: 922.14 at locus 0's Δ=212.33 … 19589.56 at locus 9's
   Δ=4510.66, the linear 4.343·Δ form);
4. the product record mirrors the receipt's serialization convention
   (the territory/locality receipts pass through serde_json's
   from_slice re-emission, whose float parse is up to 1 ULP from the
   correctly-rounded value at some decimals — the product record takes
   the identical round trip under the same conventions, keeping its
   f64 fields bit-identical to the receipt's; under the default
   convention neither record is reparsed);
5. the truth-QV test mode reproduces the committed QV receipts
   field-identically on all 9 expressible chrMT loci (assignment,
   likegt agreement, the four pair-count matrices, chosen pairs,
   error/identity, QV 60.0, perfect);
6. the CLI is deterministic across its own runs (byte-identical
   `calls.jsonl`; the receipt differs only in the walls/rss timing
   fields).

**chrI: CLOSED** — the regenerated chrI receipt anchored 18/18 called
folds + truth ranks vs the committed chrI QV receipts; the CLI
receipt parity TRUE (mod walls/rss); the product calls field-identical
per locus (folds, best LL, E, alternative gap, QUAL); the test mode
field-identical to the committed chrI QV receipts on all 18 expressible
loci — including the four named non-rank-1 residuals (L2/L13/L16/L18)
with their exact error/QV numbers.

**The instrument's own checker, reused**: the committed checker
(`check-realign-scoring.py --territory --locality`) runs over the
regenerated chain the CLI is proven field-identical to — chrMT: **ALL
PHASES PASS, 1,054,875 checks / 0 failures** (the exact committed
gate-2 record; phase-8N aggregate: truth rank-1 9 -> 9 HOLD, CONVERT
[], REGRESS []); chrI: **ALL PHASES PASS, 4,228,180 / 0** (the exact
committed gate-2 record; phase-8N: rank-1 14 -> 14, CONVERT [],
REGRESS [], residual moved [13, 16]). The checker's own cleaned
prerequisites were regenerated through the committed runners first
(the territory 4-wide BEFORE receipts, the locality 4-wide
timered-base receipts, the anchor-projection receipts).

**The full 17-component sample**: driven end-to-end through the CLI by
`genome/instrumented/run-cli-product-fleet.sh` (marker-idempotent per
component: the census regeneration through the committed runner, the
CLI product run + test mode, the per-component gate vs the surviving
committed QV receipts on every field), with the sample aggregate at
`genome/instrumented/cli-product-aggregate.py` (the before/after
receipt vs the committed assessment numbers of record: 880/1,132 =
77.7%, QV median 60.00, p10 14.42). **chrII: CLOSED end-to-end through
the CLI** — 81 loci / 76 expressible / 68 rank-1 (exactly the
committed aggregate table's row), the gate FIELD-IDENTICAL vs the
committed QV receipts on all 76 expressible loci. **chrIII: CLOSED
end-to-end through the CLI** — 38 loci / 28 expressible / 18 rank-1
(exactly the committed aggregate table's row), the gate FIELD-IDENTICAL
vs the committed QV receipts on all 28 expressible loci.
**chrIV: CLOSED end-to-end through the CLI** — 167 loci / 143
expressible / 115 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
143 expressible loci.
**chrV: CLOSED end-to-end through the CLI** — 61 loci / 56
expressible / 39 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
56 expressible loci.
**chrVI: CLOSED end-to-end through the CLI** — 29 loci / 25
expressible / 16 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
25 expressible loci.
**chrVII: CLOSED end-to-end through the CLI** — 117 loci / 100
expressible / 75 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
100 expressible loci.
**chrVIII: CLOSED end-to-end through the CLI** — 60 loci / 50
expressible / 40 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
50 expressible loci.
**chrIX: CLOSED end-to-end through the CLI** — 45 loci / 40
expressible / 29 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
40 expressible loci.
**chrXI: CLOSED end-to-end through the CLI** — 66 loci / 66
expressible / 47 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
66 expressible loci.
**chrX: CLOSED end-to-end through the CLI** — 80 loci / 71
expressible / 54 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
71 expressible loci.
**chrXIII: CLOSED end-to-end through the CLI** — 101 loci / 86
expressible / 68 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
86 expressible loci.
**chrXIV: CLOSED end-to-end through the CLI** — 82 loci / 72
expressible / 60 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
72 expressible loci.
**chrXVI: CLOSED end-to-end through the CLI** — 107 loci / 92
expressible / 71 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
92 expressible loci.
**chrXV: CLOSED end-to-end through the CLI** — 119 loci / 102
expressible / 82 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
102 expressible loci.
**chrXII: CLOSED end-to-end through the CLI** — 119 loci / 98
expressible / 75 rank-1 (exactly the committed aggregate table's
row), the gate FIELD-IDENTICAL vs the committed QV receipts on all
98 expressible loci.

**THE FLEET IS 17/17**: every component closed end-to-end through
the CLI, every per-component row exactly the committed assessment
table's row, every gate FIELD-IDENTICAL vs the committed QV
receipts on every expressible locus — the translation held with
zero divergence. The sample before/after receipt lands at
`cli-product-gate-AGGREGATE.txt`.

The long pole is the per-component
census regeneration (the validation dir was cleaned of the census
scratch; ~6 h of walls across the 15 fleet components); the
record-derivation cache is whole-sample (one cache serves every
component).

## THE SAMPLE BEFORE/AFTER RECEIPT — the fleet 17/17 aggregate

`cli-product-gate-AGGREGATE.txt` (written by
`genome/instrumented/cli-product-aggregate.py`, receipt-side): the
full diploid validation sample measured END-TO-END THROUGH THE CLI
vs the committed assessment numbers of record.

| component | loci | expressible | truth-rank-1 | perfect |
|---|---|---|---|---|
| chrMT | 14 | 9 | 9 | 9 |
| chrI | 21 | 18 | 14 | 14 |
| chrII | 81 | 76 | 68 | 68 |
| chrIII | 38 | 28 | 18 | 18 |
| chrIV | 167 | 143 | 115 | 115 |
| chrV | 61 | 56 | 39 | 39 |
| chrVI | 29 | 25 | 16 | 16 |
| chrVII | 117 | 100 | 75 | 75 |
| chrVIII | 60 | 50 | 40 | 40 |
| chrIX | 45 | 40 | 29 | 29 |
| chrX | 80 | 71 | 54 | 54 |
| chrXI | 66 | 66 | 47 | 47 |
| chrXII | 119 | 98 | 75 | 75 |
| chrXIII | 101 | 86 | 68 | 68 |
| chrXIV | 82 | 72 | 60 | 60 |
| chrXV | 119 | 102 | 82 | 82 |
| chrXVI | 107 | 92 | 71 | 71 |
| **TOTAL** | **1,307** | **1,132** | **880** | **880** |

**THE AGGREGATE: loci 1,307, expressible 1,132, truth rank-1 880,
perfect 880 — 880/1,132 = 77.7%.** The QV distribution over all 1,132
expressible loci: median 60.00, p10 14.42 (nearest-rank), QV>=40 at
78.18%.

**THE TRANSLATION VERDICT: 77.7% SURVIVES THE INSTRUMENT-TO-CLI
TRANSLATION EXACTLY.** The committed assessment numbers of record
(the locality fleet) were 880/1,132 expressible = 77.7%, QV perfect
880, median 60.00, p10 14.42, QV>=40 78.18% — the CLI end-to-end
sample reproduces every one of them, and every per-component row is
the committed table's row exactly. THE
DIVERGENCES FOUND AND DIAGNOSED, all summary-script-side, zero
per-locus divergence anywhere: (1) the aggregate script's first QV
line collected the QV only within the rank-1 branch, printing the
rank-1-only distribution (p10 60.00, >=40 100%) instead of the full
expressible cohort — fixed to collect every expressible locus's QV
(the per-locus values were always field-identical, proven by every
per-component gate); (2) the p10 index convention (int-index 112 vs
the committed tables' nearest-rank 113) produced 14.41 vs the
committed 14.42 — aligned to the committed nearest-rank convention.
No instrument-side divergence exists: the CLI's per-locus outputs
are field-identical to the committed receipts at every expressible
locus of all 17 components (the per-component gates), and the two
gate components' regenerated chains carry the full committed
checker verification (chrMT 1,054,875/0, chrI 4,228,180/0).

## THE TEST-MODE DEMONSTRATION — the truth-QV path, spot-checked

The `--truth-qv` test mode (the ONLY assessment surface, requiring
the truth haplotypes as explicit input, never in the default
output) is demonstrated at chrMT (the dedicated
`cli-product-chrMT-truthqv` run of record): `--truth-qv
S288C#0#chrMT SK1#0#chrMT` — the first haplotype is the component
lane (the axis-row truth-first convention), the material is each
fold's own panel-window sequence, the assignment yardstick the
owner's injective no-replacement rule with likegt's mean-identity
rule computed beside it, biWFA gap-affine End2End at likegt's
penalties. The run's `calls.jsonl.truth-qv.jsonl` reproduces the
committed QV receipts FIELD-IDENTICALLY on all 9 expressible loci:
9 rank-1, 9 perfect, QV 60.0 at each — and the same proof stands
fleet-wide, every component's gate comparing its
`calls.jsonl.truth-qv.jsonl` against the committed QV receipts
field-by-field (truth_rank, rank1, assignment, likegt_assignment,
assignment_disagreement, pairs, chosen_pairs, total_edits,
total_columns, mismatches, gap_columns, error, identity, qv,
perfect, called_folds/truth_folds strain sets) — FIELD-IDENTICAL at
all 17 components, zero diffs. The default `calls.jsonl` carries
NO truth keys anywhere.

## Production constraints

Bounded memory (the 64 GiB RSS guard, the lane cache, the degenerate-locus
wall fix — all inherited from the committed machinery), streaming ordered
emission per locus, no thresholds, no tuning constants, the scoreboard
machinery unmodified, no PR push (the owner pushes). The census stage is
currently an input produced by the committed multicensus runner
(`run-cosine-multicensus.sh`, truth-free); folding the census production
into the CLI itself is named follow-up work, as is the `--loci all`
enumeration (the axis-interval input) and per-component derive-cache
builds for the fleet run.
