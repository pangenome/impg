# Diagnose C4 Query End Artifacts

Task: `diagnose-c4-query-2`
Date: 2026-06-09

## Inputs

Artifacts inspected:

- Baseline one-chunk run: `/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z`
- Fuller comparison run: `/home/erikg/impg/data/run_fuller_c4_20260609T160000Z`

The target BED is identical in both trees:

```text
GRCh38#0#chr6	31982056	32035418	C4_GRCh38_53kb
```

## Finding

The coordinate-level end artifacts are already present in the query-selected
sequence set before AGC extraction, local seed induction, localized polishing,
sorting, rendering, or `gfalook` visualization. The boundary for the exact
left-end dropout and the high-confidence 3-prime tail paths is syng query
interval collection.

The reason expanded localized polish looks worse is separate: moderate and
strong runs apply many exact-path-preserving SPOA replacement chunks. Those
chunks reduce some path steps but add segments, singleton bp, duplicate bp, and
white-space. This is introduced at localized polishing, not at sorting,
rendering, or `gfalook`.

No code was changed in this pass. A small safe fix is not obvious without
adding anchor-provenance tests first: a C4-length outlier filter would be
dataset-specific and could discard real copy-number or structural sequence. The
patch plan below gives the exact code path to fix.

## Command Evidence

Selected sequence spans were extracted from `combined.fa` with:

```sh
awk '/^>/{if (name!="") print name"\t"len; name=substr($0,2); len=0; next} {len+=length($0)} END{if (name!="") print name"\t"len}' .../combined.fa
```

Run-to-run interval-coordinate comparison by stable path key showed one
moderate-only coordinate change:

```text
HG01960#2#JBHIHN010000011.1	diag=32045898-32158486 len=112588	mod=32078636-32158486 len=79850	delta_len=-32738
```

The same comparison against strong printed no differences, so strong matches
the one-chunk diagnostic interval set exactly.

Graph-build summaries show the difference exists before graph construction:

```text
diagnose: sequences=466 total_sequence_bp=53165920 raw_paf_records=541854 filtered_paf_records=431600 paf_length_mismatches=0 seqwish_segments=3850
moderate: sequences=466 total_sequence_bp=53133182 raw_paf_records=541784 filtered_paf_records=431535 paf_length_mismatches=0 seqwish_segments=3852
strong:   sequences=466 total_sequence_bp=53165920 raw_paf_records=541854 filtered_paf_records=431600 paf_length_mismatches=0 seqwish_segments=3850
```

Local seed logs agree:

```text
diagnose: [syng local seed] collected 466 local sequence(s), 53165920 bp
moderate: [syng local seed] collected 466 local sequence(s), 53133182 bp
strong:   [syng local seed] collected 466 local sequence(s), 53165920 bp
```

All three local seeds report `path_validation=pass`; the failures here are not
exact path-spelling corruption inside graph construction.

## 5-prime Left Dropout

The only run-differential 5-prime dropout I found is `HG01960#2#JBHIHN010000011.1`.
Moderate clips the selected interval start 32,738 bp to the right while leaving
the end unchanged. The extracted sequence is therefore 32,738 bp shorter before
AGC extraction or seqwish runs.

| Run | Selected path interval | Length | Delta from one-chunk | Sorted GFA first-last step count | Boundary |
| --- | --- | ---: | ---: | --- | --- |
| one-chunk baseline | `HG01960#2#JBHIHN010000011.1:32045898-32158486` | 112,588 | 0 | `7+` to `3814+`, 3,775 steps | syng query interval collection |
| moderate | `HG01960#2#JBHIHN010000011.1:32078636-32158486` | 79,850 | -32,738 | `518+` to `3963+`, 3,040 steps | syng query interval collection |
| strong | `HG01960#2#JBHIHN010000011.1:32045898-32158486` | 112,588 | 0 | `7+` to `4121+`, 3,725 steps | syng query interval collection |

The moderate first node `518+` is a singleton start in the sorted graph, while
the one-chunk/strong interval starts on the low-depth but shared left-start
class `7+`. That is the visual left dropout. AGC is not the introducing
boundary: it extracts the interval requested by syng, and `path_spellings.tsv`
contains the clipped moderate path name and clipped length.

There is also a stable short-interval class already present in the one-chunk
baseline. These are not run-differential regressions, but if a full C4-copy
span is expected they are the coordinate-level candidates for persistent
low-depth end dropouts:

| Path interval | Length |
| --- | ---: |
| `HG00099#2#CM087361.1:31928455-32008252` | 79,797 |
| `HG00146#1#CM090015.1:31922607-32002404` | 79,797 |
| `HG00272#2#CM094216.1:31931690-32011487` | 79,797 |
| `HG00280#2#CM087160.1:31931045-32010842` | 79,797 |
| `HG00738#1#CM086647.1:31951624-32031421` | 79,797 |
| `HG02258#2#CM085837.1:31934546-32014343` | 79,797 |
| `HG06807#2#CM101434.1:32040804-32120601` | 79,797 |
| `NA19468#2#JBHDXY010000045.1:31915811-31995608` | 79,797 |
| `HG03017#2#CM088928.1:31927766-32007564` | 79,798 |
| `NA18879#1#CM099716.1:31958945-32038743` | 79,798 |
| `HG01952#2#CM087999.1:31995790-32075638` | 79,848 |
| `HG02027#2#JBHDTY010000008.1:5188422-5268271` | 79,849 |
| `HG02132#2#CM086908.1:31997546-32077395` | 79,849 |
| `NA18565#2#CM094342.1:31929266-32009115` | 79,849 |
| `NA18945#2#CM101630.1:32018016-32097865` | 79,849 |
| `NA18960#2#CM101559.1:31867170-31947019` | 79,849 |

Because these short names and lengths are already in the query-selected FASTA,
their boundary is also syng query interval collection. I did not find evidence
that AGC extraction, local seed induction, sorting, rendering, or `gfalook`
shortens them.

## 3-prime Extra Tails

The high-confidence 3-prime/long-tail candidates are the longest selected
intervals. They are already present in `combined.fa` and `path_spellings.tsv`
before localized polishing. The `path_spellings.tsv` entries for representative
paths are identical across one-chunk and strong, and moderate differs only for
the HG01960 left clip above.

The table uses 118,900 bp as the common full-copy reference class only to size
the outlier; this is not a graph-quality gate.

| Path interval | Selected length | Excess over 118,900 | One-chunk sorted first-last step count | Boundary |
| --- | ---: | ---: | --- | --- |
| `NA18952#1#JBHIJS010000065.1:31855636-32040012` | 184,376 | 65,476 | `29+` to `3814+`, 6,589 steps | syng query interval collection |
| `NA18948#2#CM101588.1:31979894-32164269` | 184,375 | 65,475 | `28+` to `3815+`, 6,588 steps | syng query interval collection |
| `HG02392#1#CM099882.1:31884401-32056195` | 171,794 | 52,894 | `2+` to `3831+`, 6,175 steps | syng query interval collection |
| `HG00232#2#JBHIKN010000040.1:31904039-32055683` | 151,644 | 32,744 | `29+` to `3819+`, 5,293 steps | syng query interval collection |
| `HG02083#1#JBHDTV010000048.1:31903697-32055337` | 151,640 | 32,740 | `28+` to `3838+`, 5,281 steps | syng query interval collection |
| `HG01261#2#JBHDVH010000044.1:31939220-32090858` | 151,638 | 32,738 | `29+` to `3814+`, 5,278 steps | syng query interval collection |
| `HG01496#2#CM086388.1:31868198-32019836` | 151,638 | 32,738 | `28+` to `3810+`, 5,278 steps | syng query interval collection |
| `HG01975#1#JBJUVS010000062.1:32044917-32196555` | 151,638 | 32,738 | `28+` to `3814+`, 5,278 steps | syng query interval collection |
| `HG01981#2#JBHDRW010000010.1:31934383-32086021` | 151,638 | 32,738 | `28+` to `3811+`, 5,278 steps | syng query interval collection |
| `HG02004#2#JBHDRU010000054.1:31886084-32037722` | 151,638 | 32,738 | `28+` to `3818+`, 5,278 steps | syng query interval collection |
| `HG01081#1#CM088567.1:31939645-32091282` | 151,637 | 32,737 | `28+` to `3814+`, 5,282 steps | syng query interval collection |
| `HG01978#2#CM089273.1:31881693-32033284` | 151,591 | 32,691 | `28+` to `3814+`, 5,260 steps | syng query interval collection |

Representative `path_spellings.tsv` evidence:

```text
diagnose: HG02392#1#CM099882.1:31884401-32056195	171794
diagnose: NA18948#2#CM101588.1:31979894-32164269	184375
diagnose: NA18952#1#JBHIJS010000065.1:31855636-32040012	184376
moderate: HG02392#1#CM099882.1:31884401-32056195	171794
moderate: NA18948#2#CM101588.1:31979894-32164269	184375
moderate: NA18952#1#JBHIJS010000065.1:31855636-32040012	184376
strong:   HG02392#1#CM099882.1:31884401-32056195	171794
strong:   NA18948#2#CM101588.1:31979894-32164269	184375
strong:   NA18952#1#JBHIJS010000065.1:31855636-32040012	184376
```

This supports the user's suspicion: the long tails are selected as query
intervals by syng. AGC extraction and local seed induction preserve those
intervals exactly; they do not introduce the extra sequence.

## Boundary Matrix

| Artifact | Introduced at | Evidence against later boundaries |
| --- | --- | --- |
| Moderate `HG01960#2` 5-prime/left dropout | syng query interval collection | The clipped coordinate appears in moderate `combined.fa` and `path_spellings.tsv`; AGC and local seed report zero length mismatches and `path_validation=pass`. |
| Stable short-span end-dropout candidates | syng query interval collection | The short coordinate spans are path names in the selected FASTA before graph construction. |
| High-confidence 3-prime/long-tail candidates | syng query interval collection | The long coordinate spans are path names and path spellings in the seed graph before polishing; baseline, moderate, and strong carry the same long intervals. |
| Moderate/strong visual scars, extra segments, extra singleton bp, white-space growth | localized polishing | `localized_polish_summary.tsv` shows exact-preserving SPOA replacements increasing graph size before `gfasort`, render, or `gfalook`. |
| Sorting/rendering/`gfalook` | not introducing exact artifacts | Sorted reports and GFA `P` lines reflect existing path intervals; these tools do not fetch sequence or alter path coordinates. |

## Baseline vs Expanded Polish

The one-chunk run is better because it applies only the small early replacement
`Poa_crush01_candidate_span24bp_med24bp_cov457of466_sites457_steps91-94`.
Moderate and strong expand into repeated left-frontier chunks and late/right
chunks, while the exact-path hard gate stays green.

Localized summary:

| Profile | Applied chunks | Selected chunk bp | Segments before-after | Path steps before-after | White-space bp before-after | Singleton bp before-after |
| --- | ---: | ---: | --- | --- | --- | --- |
| one-chunk | 1 | 2,024 | 3,850-3,851 | 1,785,229-1,784,781 | 3,197,050,792-3,197,892,063 | 2,531-2,531 |
| moderate iter 0 | 11 | 74,979 | 3,852-3,899 | 1,783,604-1,754,112 | 3,192,653,848-3,217,259,508 | 2,722-2,993 |
| moderate iter 1 | 9 | 79,824 | 3,899-4,000 | 1,754,112-1,750,053 | 3,217,259,508-3,323,525,830 | 2,993-3,055 |
| strong iter 0 | 40 | 307,853 | 3,850-4,111 | 1,785,229-1,753,671 | 3,197,050,792-3,342,391,759 | 2,531-4,804 |
| strong iter 1 | 35 | 314,075 | 4,111-4,158 | 1,753,671-1,753,671 | 3,342,391,759-3,350,998,371 | 4,804-4,863 |

Final sorted graph-report diagnostics:

| Profile | Segments | Links | Segment bp | Singleton bp | Segment white-space bp | Sparse coverage bp | Duplicate sequence fraction | Path white-space p99-max |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| one-chunk | 3,851 | 5,220 | 97,960 | 2,531 | 5,656,162 | 11,203 | 0.662685 | 210-26,933 |
| moderate | 4,000 | 5,391 | 103,754 | 3,055 | 8,333,591 | 16,882 | 0.677000 | 262-29,096 |
| strong | 4,158 | 5,560 | 110,122 | 4,863 | 11,269,943 | 23,251 | 0.698172 | 262-23,107 |

Specific chunk families that increase scars:

- Moderate repeats left-frontier chunks after the one-chunk site:
  `span298 steps45-47`, `span130 steps67-69`, `span265 steps77-94`,
  then larger early chunks `span1920 steps188-206`, `span401 steps495-518`,
  and `span2443 steps536-586`.
- Moderate also applies late/right chunks such as `span4 steps3400-3403`,
  `span66 steps3520-3524`, `span18 steps3585-3602`, and
  `span30 steps3753-3775`.
- Strong applies the same families plus a much larger late/right set:
  `span323 steps3243-3245`, `span729 steps3375-3391`,
  `span51 steps3447-3449`, `span95 steps3471-3473`,
  `span505 steps3515-3528`, `span234 steps3557-3563`,
  `span225 steps3632-3648`, `span288 steps3670-3672`, and
  `span2004 steps3707-3748`, then repeats many of those windows in iteration 1.

The replacements are not exact path corruption. Logs show they are accepted
with `hard_gate=exact_path_corruption` and `quality_gate_used=false`. The
metrics above are diagnostic only: they explain the visual scars but should not
be promoted to hard graph-quality gates.

## Patch Plan

1. Preserve syng anchor provenance through the query-output path. Today
   `syng_intervals_to_adjusted` in `src/main.rs` drops anchor data from
   `SyngIndex::query_region_ext`, then `merge_query_adjusted_intervals` merges
   purely by sequence coordinate. Add a debug TSV for GFA/FASTA local seed:
   path key, pre-merge interval, post-merge interval, strand, anchor count,
   query anchor min/max, target anchor min/max, and emitted fetch interval.
2. Add unit tests in `src/syng.rs` around `query_region_ext` and
   `merge_intervals_with_anchors` for a synthetic C4-like case where two
   homologous blocks on one target path are separated by about 32.7 kb or
   65.5 kb. Expected behavior: the emitted interval must not absorb an
   unanchored terminal tail, and a clipped interval like the moderate
   `HG01960#2` start shift must be rejected or split unless its anchor coverage
   supports the clipped boundary.
3. Change syng interval emission so merged intervals are trimmed or split by
   anchor-supported target span plus requested padding/query extension. Do not
   use C4 length cutoffs. Use anchor coverage and chain continuity, otherwise
   real long alleles and copy-number sequence could be lost.
4. Route GFA/FASTA local-seed interval collection through that anchored merge
   path. Keep exact path spelling validation as the hard failure. Use graph
   metrics only to diagnose replacement behavior, not to reject exact paths.
5. Until the syng interval fix lands, keep the C4 release/default diagnostic
   profile at `max-total-chunks=1` for localized polish. Moderate/strong are
   useful stress tests, but they demonstrably add exact-preserving topology and
   should not be treated as better graph output.

## Validation Performed

- Re-read required task criteria and used both required artifact trees.
- Confirmed exact path corruption is not occurring inside graph build or
  localized polish: local seed and polish reports have `path_validation=pass`,
  and graph-build summaries have `paf_length_mismatches=0`.
- No code was changed, so targeted tests and `cargo build`/`cargo test` were
  not run. The exact next code patch is the anchored syng interval emission
  plan above.
