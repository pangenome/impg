# Sequence-k syng filter evaluation

Task: `evaluate-sequence-k`

Lineage:

- Upstream implementation context: `fix-slow-syng`, commit `9b48a46`, which made full syng start-count tests tractable.
- Primitive change: explicit `gfa:syng` raw-materialization filter knob. No wrapper filter is added; all behavior is requested in the engine string and logged by `syng2gfa`.
- Implementation lineage: existing `min-run` weak-shared-syncmer model was generalized with longer exact sequence-context support. Unsupported shared syncmer topology is split into private per-occurrence segments, while existing frequency-filtered syncmers are still removed and bridged.

## Implementation

New syng mask parameter:

```text
gfa:syng:mask,sequence-k=191
```

Aliases accepted by the parser:

```text
sequence-k
seq-k
sequence-span
sequence-context
min-sequence-span
min-sequence-context
min-sequence-bp
mask-sequence-k
mask-seq-k
```

Default is disabled: `min_sequence_span_bp = 0`.

The sequence-context filter runs before bluntification/crush, inside syng raw GFA materialization. For each selected local path walk it records syncmer node id and source bp position. A shared syncmer node remains shared if either condition holds:

- it is in a repeated consecutive-syncmer run satisfying the existing `min-run` model; or
- it is in a repeated exact consecutive-syncmer window whose source-position span is at least `sequence-k` bp.

Weak shared nodes are split into private per-occurrence syncmer segments. This was a correctness pivot during the task: the first drop-node version changed path spellings on C4. Splitting reduced that damage on C4 for `127` and `191`, but it is still not robust after bluntification in all cases; see validation below.

## C4 Commands

Output root:

```text
/home/erikg/impg/data/sequence_k_syng_filter_split_20260529T205140Z
```

Primary region and inputs:

```bash
SYNG=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng
AGC=/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc
REGION='GRCh38#0#chr6:31891045-32123783'
```

Raw sweep engines:

```bash
gfa:syng:mask,min-run=3:nosort
gfa:syng:mask,sequence-k=127:nosort
gfa:syng:mask,sequence-k=191:nosort
gfa:syng:mask,sequence-k=311:nosort
gfa:syng:mask,sequence-k=191,min-run=3:nosort
```

Per raw row:

```bash
/usr/bin/time -v target/release/impg query -t 32 \
  -a "$SYNG" -r "$REGION" --sequence-files "$AGC" -d 50k \
  -o "$ENGINE" -O "$dir/run.nosort" -v 1
gfasort -i "$dir/run.nosort.gfa" -o "$dir/run.Ygs.gfa" -p Ygs -t 32
target/release/impg graph-report -g "$dir/run.Ygs.gfa" \
  -o "$dir/run.Ygs.report.tsv" --format tsv
python3 /tmp/find_stringy_bubbles.py "$dir/run.Ygs.gfa" > "$dir/stringy.txt"
```

Full C4 crush row used the highest path-stable threshold from the raw sweep:

```bash
ENGINE='gfa:syng:mask,sequence-k=191:crush,method=auto,aligner=fastga,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort'
/usr/bin/time -v target/release/impg query -t 32 \
  -a "$SYNG" -r "$REGION" --sequence-files "$AGC" -d 50k \
  -o "$ENGINE" -O "$out/c4_seq191_crush/run.nosort" -v 1
```

## C4 Raw Sweep

`seq weak` is the number of shared local syncmer nodes that did not meet the configured support rule. `seq split` excludes nodes already removed by the frequency filter. `seq clones` is the number of private per-occurrence syncmer segments created.

| Run | Wall | Seq weak | Seq split | Seq clones | Run-supported | Seq-supported | Segments | Segment bp | bp-weighted cov | Singleton bp | WS bridges | Stringy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `min-run=3` | 4:58.80 | 2 | 2 | 4 | 14911 | 0 | 18050 | 389369 | 261.020998 | 70920 | 22080 | 506 |
| `sequence-k=127` | 5:10.11 | 6 | 6 | 23 | 0 | 14907 | 18065 | 389829 | 260.712992 | 71488 | 23170 | 506 |
| `sequence-k=191` | 4:39.59 | 10 | 10 | 31 | 0 | 14903 | 18069 | 389939 | 260.639446 | 71708 | 22447 | 507 |
| `sequence-k=311` | 4:40.67 | 12 | 12 | 35 | 0 | 14901 | 18071 | 389988 | 260.606708 | 71802 | 23810 | 507 |
| `sequence-k=191,min-run=3` | 4:40.82 | 2 | 2 | 4 | 14911 | 14903 | 18050 | 389369 | 261.020998 | 70920 | 22171 | 506 |

Raw interpretation:

- The standalone sequence-k filter has very little effect on C4. It weakens 6, 10, or 12 shared syncmer nodes out of 14,913 shared local nodes.
- Splitting those nodes slightly increases segment count and singleton bp. It does not reduce C4's stringy heuristic; `191` and `311` are one stringy candidate worse than `min-run=3`.
- Combining `sequence-k=191` with `min-run=3` collapses back to the `min-run=3` behavior on the main quality metrics, because the existing consecutive-run model rescues the sequence-k weak nodes that matter here.

## Path Validation

Validator:

```bash
cargo run --release --example compare_gfa_paths -- expected.gfa observed.gfa
```

C4 results:

| Comparison | Expected paths | Observed paths | Missing | Extra | Spelling mismatches |
| --- | ---: | ---: | ---: | ---: | ---: |
| `min-run=3` raw vs `sequence-k=127` raw | 465 | 465 | 0 | 0 | 0 |
| `min-run=3` raw vs `sequence-k=191` raw | 465 | 465 | 0 | 0 | 0 |
| `min-run=3` raw vs `sequence-k=311` raw | 465 | 465 | 0 | 0 | 1 |
| `min-run=3` raw vs `sequence-k=191,min-run=3` raw | 465 | 465 | 0 | 0 | 0 |
| `sequence-k=191` raw vs `sequence-k=191 + crush` | 465 | 465 | 0 | 0 | 0 |

The `sequence-k=311` mismatch is not a dropout: all paths are present. It is a post-blunt spelling change caused by splitting a syncmer into a private copy before bluntification. In the mismatching path, bluntification cut the shared original to a 30 bp segment but left the private copy as 34 bp, shifting the path by 4 bp. That makes `311` unsafe under the current pre-blunt split implementation.

I also added `examples/validate_gfa_path_sources.rs` while investigating source-interval validation. It currently reports that regional syng GFA path names span larger source-coordinate hulls than the path-spelled graph sequence on these HPRC extractions, so the project-standard GFA-to-GFA spelling validator above is the reliable gate for this evaluation. A follow-up should decide whether syng regional path names need tighter interval metadata for direct AGC source validation.

## C4 Crush Comparison

Controls:

- Current best syng/crush: `/home/erikg/impg/data/c4_diagnose_residual_two_20260529T141339Z/outward_highbp_r3/run.Ygs.gfa`
- PGGB control: `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa`

The current-best wall time is the recorded crush-on-existing-graph step, not a fresh syng query.

| Graph | Wall | Segments | Links | Segment bp | bp-weighted cov | Singleton bp | WS bridges | Stringy |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Current best syng/crush | 2:32.56 | 16477 | 23288 | 384319 | 264.450808 | 30925 | 29037 | 100 |
| PGGB control | 13:37.90 | 13288 | 16240 | 234524 | 454.919215 | 2890 | 22851 | 12 |
| `sequence-k=191` raw | 4:39.59 | 18069 | 20955 | 389939 | 260.639446 | 71708 | 22447 | 507 |
| `sequence-k=191 + crush` | 34:40.03 | 19716 | 23417 | 544441 | 186.674929 | 88416 | 29647 | 490 |

Crush log summary for `sequence-k=191 + crush`:

```text
crush: 12 resolved, 0 bailed, 12 candidates seen across 3 rounds
Syng query complete ... in 31m14.027s
Maximum resident set size: 102308324 KB
```

The `sequence-k=191 + crush` row does not fix the residual C4 two-cluster bubbles. It resolves the same number of local crush candidates as the prior `sequence-k=311` trial shape, but the final graph is worse than the current-best syng/crush on segment count, segment bp, bp-weighted coverage, singleton bp, and stringy count. It is much farther from PGGB than the current best.

## Simpler Sanity Check

Region:

```text
GRCh38#0#chr1:1000000-1010000
```

This region was intended as a small sanity check, but local syng extraction expanded to 12,835 paths and a 60 MB GFA, so I used unsorted GFA outputs and skipped `gfasort`.

| Run | Wall | S | L | P | Shared | Weak | Split | Seq-supported |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| `min-run=3` | 5:14.18 | 79823 | 85177 | 12835 | 54377 | 27 | 25 | 0 |
| `sequence-k=191` | 5:14.78 | 92825 | 98161 | 12835 | 54377 | 557 | 553 | 53820 |

Path comparison:

```text
expected_paths        12835
observed_paths        12835
missing_paths         0
extra_paths           0
spelling_mismatches   1148
```

This is the strongest negative result in the evaluation. The option can preserve C4 at `191`, but it is not generally path-stable under pre-blunt splitting.

## PNG

Rendered and uploaded:

```text
https://hypervolu.me/~erik/impg/c4-sequence-k-191-crush.png
```

Upload confirmation:

```text
-rw-r--r-- 1 erik erik 2.2M May 29 22:17 www/impg/c4-sequence-k-191-crush.png
```

## Recommendation

Do not enable sequence-k filtering by default.

Recommended default remains:

```text
sequence-k=0
```

Keep the existing min-run behavior as the safer explicit model for now. If this line of work continues, the next implementation should make sequence-context filtering path-stable before bluntification. Two plausible directions:

1. Split weak topology after a path-preserving blunting step, then validate and document that the filter is post-blunt rather than raw-layer.
2. Teach raw-layer splitting to preserve the exact bluntified path spellings by carrying the shared node's blunt cut coordinates into private copies.

Until one of those exists, `sequence-k` should remain an experimental explicit option only. The current implementation is useful for evaluating the hypothesis, but the chr1 sanity failure and the C4 `311` failure rule it out as a production default.
