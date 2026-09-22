# Validate Flank-Aware Crush On C4

Date: 2026-06-06

Task: `validate-flank-aware-crush-c4`

Implementation under evaluation: `b904c109a46038b63bfe5385b0e9d48a58532529`
(`b904c10`, `impg 0.4.1`). The dependency note named implementation commit
`aa0fcb5`; this WG worktree had the later implementation head `b904c10`, and
the installed binary was rebuilt from this checkout before running the
experiment.

## Result

Flank-aware crush is behaviorally valid, but it is not a clear C4 quality fix in
this bounded run.

The flank-aware arm (`--replacement-flank-bp 64`) preserved all C4 paths exactly
and accepted 23 more replacements than the no-flank arm, with fewer final
segments and path steps. That is a real local compaction improvement. However,
the C4 residual bubble/span and white-space proxy metrics did not materially
improve: the largest residual POVU reference span remained 1312 steps, the
path-white-space maximum remained 25254 bp, and the bridge count was effectively
unchanged. The uploaded Ygs/gfalook renders should therefore be treated as a
diagnostic comparison, not as evidence that flanks fixed the observed C4
underalignment/splitting.

Sanity behavior did pass: the synthetic locus exercised repeated context and
short one-sided path-end flanks, the targeted Rust sanity suite exercised
reverse-orientation trimming/lacing and hard path-corruption rejection, and no
metric or quality-score guard filtered replacements.

## Artifacts

Main run directory:

`/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z`

Important local files:

| Artifact | Path |
|---|---|
| command log | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/commands.sh` |
| run summary | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/run_summary.json` |
| metrics table | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/reports/metrics.tsv` |
| candidate accounting | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/reports/candidate_accounting.tsv` |
| exact path preservation | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/validation/path_preservation.tsv` |
| sanity test table | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/reports/sanity_tests.tsv` |
| artifact table | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/artifacts.tsv` |
| PNG URL table | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/renders/uploaded_urls.tsv` |

Rendered PNGs uploaded to Hypervolume:

| Graph | URL |
|---|---|
| C4 flank disabled | http://hypervolu.me/~erik/impg/flank-aware-crush-c4-c4.flank0-flank_aware_crush_c4_20260606T112300Z.png |
| C4 flank64 | http://hypervolu.me/~erik/impg/flank-aware-crush-c4-c4.flank64-flank_aware_crush_c4_20260606T112300Z.png |
| current best C4 default-selfloop | http://hypervolu.me/~erik/impg/c4-spoa2k-default-selfloop-verify-20260606T065905Z.Ygs.mean-depth.png |
| current best 1:1 no-scaffold diagnostic | http://hypervolu.me/~erik/impg/c4-k311-1to1-noscaffold-spoa2k-default-selfloop-20260606T070148Z.Ygs.mean-depth.png |
| PGGB control | https://hypervolu.me/~erik/impg/c4-pggb-control.png |

Both new PNGs were produced by normal Ygs/gfalook path view commands and are PNG
images around 3200 x 1580 pixels.

## Inputs And Controls

| Role | Path |
|---|---|
| C4 seed input | `/home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa` |
| current best C4 graph | `/home/erikg/impg/data/c4_spoa2k_default_selfloop_verify_20260606T065905Z/sorted/c4.k311.seed.spoa2k.default-selfloop.Ygs.gfa` |
| current best 1:1 no-scaffold diagnostic | `/home/erikg/impg/data/c4_k311_1to1_noscaffold_shorttmp_20260606T070148Z/sorted/c4.k311.1to1.noscaffold.spoa2k.default-selfloop.Ygs.gfa` |
| PGGB control | `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.Ygs.gfa` |
| synthetic sanity GFA | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/sanity/flank_required_repeated_path_end.gfa` |

No expected control was missing. The current best default-selfloop C4 metrics
match the `c4.flank0` arm because that graph is the current non-flank default
under the same resolver/settings.

## Commands

All exact commands are in `commands.sh`. The key before/after C4 commands were:

```bash
/usr/bin/time -v /home/erikg/.cargo/bin/impg crush -g /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa -o /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank0.gfa --method poa --max-iterations 5 --max-traversal-len 2000 --max-median-traversal-len 2000 --max-total-sequence 1m --max-traversals 10k --replacement-flank-bp 0 -t 32 -v 1

/usr/bin/time -v /home/erikg/.cargo/bin/impg crush -g /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/graphs/c4.k311.seed.gfa -o /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank64.gfa --method poa --max-iterations 5 --max-traversal-len 2000 --max-median-traversal-len 2000 --max-total-sequence 1m --max-traversals 10k --replacement-flank-bp 64 -t 32 -v 1
```

The resolver was `poa` for all sanity and C4 crush runs. Unrelated parameters
were held constant between the C4 arms.

The render commands were:

```bash
/usr/bin/time -v /home/erikg/.cargo/bin/gfalook -i /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/sorted/c4.flank0.Ygs.gfa -o /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/renders/c4.flank0.Ygs.gfalook-m.png -m -x 3200 -y 1800 -a 3 -t 32 -v 1

/usr/bin/time -v /home/erikg/.cargo/bin/gfalook -i /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/sorted/c4.flank64.Ygs.gfa -o /home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/renders/c4.flank64.Ygs.gfalook-m.png -m -x 3200 -y 1800 -a 3 -t 32 -v 1
```

The `flank64` stderr log uses compacted flank trim diagnostics. The first
attempt produced an 18 GB unbounded raw stderr file because every trim
diagnostic included complete C4 path lists. The final clean run preserved the
requested/actual flank and trim fields while replacing the full candidate/path
payloads with bounded placeholders; the clean run directory is 297 MB.

## Synthetic Sanity

The external sanity GFA contains two simple indel contexts with repeated
sequence (`AAC ... TTA` appears around both candidate sites). One candidate is
at a path start and one reaches a path end, so the `flank4` run has one-sided
shortened flanks. The run verifies exact path preservation and candidate
accounting for flank-disabled and flank-enabled crush.

The targeted Rust sanity suite also ran:

```bash
/usr/bin/time -v zsh -lc 'source ./env.sh && cargo test --lib flank -- --nocapture'
```

Result: `8 passed; 0 failed; 393 filtered out`. The covered tests include
repeated occurrence context, short path-boundary flanks, reverse occurrence
trim/lacing, global resolver flank trimming, hard trim-boundary rejection, and
diagnostic formatting. Relevant tests are in `src/resolution.rs`:

| Behavior | Test |
|---|---|
| repeated same-occurrence context | `flank_aware_context_is_taken_from_the_same_repeated_occurrence` |
| requested vs actual short path flanks | `flank_aware_context_records_requested_and_short_path_boundary_flanks` |
| reverse orientation | `flank_aware_reverse_occurrence_trims_and_laces_in_original_orientation` |
| global direct resolver trim | `flank_aware_global_direct_resolvers_trim_simple_indel_context` |
| hard boundary/path-corruption rejection | `flank_aware_trim_rejects_inconsistent_boundaries_as_path_corruption_risk` |
| diagnostics | `flank_aware_diagnostics_expose_span_flanks_trim_and_orientation` |

Sanity exact path preservation:

| Run | Expected paths | Observed paths | Path-name mismatches | Sequence mismatches | Full preservation |
|---|---:|---:|---:|---:|---|
| `sanity.flank0` | 2 | 2 | 0 | 0 | true |
| `sanity.flank4` | 2 | 2 | 0 | 0 | true |

Sanity candidate accounting:

| Run | Requested flank | Accepted | Rejected | Flank observations |
|---|---:|---:|---:|---|
| `sanity.flank0` | 0 | 2 | 0 | no flank occurrences logged |
| `sanity.flank4` | 4 | 2 | 0 | actual left/right path flanks were `0:6;4:6`, demonstrating path-end shortening |

The sanity result does not depend on hidden filtering: both candidates were
accepted, hard rejection counts were zero, and `metric_based_filter_used=false`.

## Exact Path Preservation

| Run | Expected paths | Observed paths | Path-name mismatches | Sequence-spelling mismatches | Full preservation |
|---|---:|---:|---:|---:|---|
| `sanity.flank0` | 2 | 2 | 0 | 0 | true |
| `sanity.flank4` | 2 | 2 | 0 | 0 | true |
| `c4.flank0` | 465 | 465 | 0 | 0 | true |
| `c4.flank64` | 465 | 465 | 0 | 0 | true |

Full path-detail TSVs are under:

`/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/validation/*.path-detail.tsv`

## Candidate Accounting

Candidate rejection was allowed only for hard path-corruption checks. None
occurred in these runs. Budget skips are explicit length-budget exclusions, not
metric or quality-score guards.

| Run | Requested flank | Round candidate sum | Attempted/seen | Selected | Accepted | Hard rejected | Budget ineligible | Rounds |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| `sanity.flank0` | 0 | 2 | 2 | 2 | 2 | 0 | 0 | 1 |
| `sanity.flank4` | 4 | 2 | 2 | 2 | 2 | 0 | 0 | 1 |
| `c4.flank0` | 0 | 3025 | 1517 | 1517 | 1517 | 0 | 10 | 4 |
| `c4.flank64` | 64 | 3238 | 1540 | 1540 | 1540 | 0 | 10 | 5 |

For `c4.flank64`, actual logged flank distributions were:

| Field | Distribution |
|---|---|
| orientation | `forward:689341` |
| left actual path bp | `64:689341` |
| right actual path bp | `64:689341` |
| left actual canonical bp | `64:689341` |
| right actual canonical bp | `64:689341` |
| left trim bp | `64:689341` |
| right trim bp | `64:689341` |

For the synthetic `flank4` run, actual left and right path flanks included both
full requested flanks and zero-length path-end flanks:

`0:6;4:6`

Trim-boundary ambiguity/path-corruption counts were zero for every measured
arm. The C4 logs also explicitly report `threshold=disabled (diagnostic only)`
and `no replacement was rejected by ratio` for every replacement-compression
diagnostic round. Candidate-accounting TSV fields record
`quality_guard_rejections=0` and `metric_based_filter_used=false`.

## C4 Metrics

| Graph | Segments | Links | Paths | Path steps | Segment bp | Node cov mean | BP-wtd cov | Cov p10/med/p90 | WS p99 | WS max | WS bridges | POVU sites | Max residual span |
|---|---:|---:|---:|---:|---:|---:|---:|---|---:|---:|---:|---:|---:|
| `c4.flank0` | 8350 | 11289 | 465 | 3121673 | 254167 | 373.853054 | 419.761314 | 1/461/884 | 139 | 25254 | 6740 | 2829 | 1312 |
| `c4.flank64` | 8163 | 11308 | 465 | 2955270 | 256121 | 362.032341 | 416.558869 | 1/458/901 | 140 | 25254 | 6739 | 2904 | 1312 |
| current best default-selfloop | 8350 | 11289 | 465 | 3121673 | 254167 | 373.853054 | 419.761314 | 1/461/884 | 141 | 25254 | 6740 | 2829 | 1312 |
| current best 1:1 no-scaffold | 8301 | 11147 | 465 | 3112253 | 257626 | 374.925069 | 414.125414 | 1/462/884 | 159 | 25250 | 7198 | 2865 | 1304 |
| PGGB control | 13288 | 16240 | 465 | 5538879 | 234524 | 416.833158 | 454.919215 | 1/465/929 | 14 | 219917 | 11876 | 0 | 1710 |

Large residual bubble/reference spans in the new C4 arms are unchanged at the
top: `1312,1298,1298,1297,1297,1297,1297,1297,1296,166` for both `flank0`
and `flank64`.

Interpretation:

- `flank64` did discover and accept more replacement opportunities: 1540 vs
  1517.
- `flank64` reduced final path steps by 166403 (about 5.3%) and segments by
  187 relative to `flank0`.
- The main residual C4 splitting indicators did not improve: residual span max
  stayed 1312, path-white-space max stayed 25254, and WS bridge count changed
  only 6740 to 6739.
- The PGGB control remains structurally different: many more path steps and
  segments, lower path-white-space p99, but much larger path-white-space max.

## Self-Loop And Same-Step Repeats

Graph-report after default normalization:

| Graph | Direct self-loop edges | Adjacent same-step path steps | Self-loop repeat runs |
|---|---:|---:|---:|
| `c4.flank0` | 0 | 0 | 0 |
| `c4.flank64` | 0 | 0 | 0 |
| current best default-selfloop | 0 | 0 | 0 |
| current best 1:1 no-scaffold | 0 | 0 | 0 |

Crush log normalization summaries:

| Graph | Direct self-loop edges before -> after | Adjacent same-step repeats before -> after | Collapsed runs | Created segments |
|---|---:|---:|---:|---:|
| `c4.flank0` | 57 -> 0 | 252715 -> 0 | 28031 | 190 |
| `c4.flank64` | 59 -> 0 | 254090 -> 0 | 28941 | 192 |

## Runtime And Resources

| Command | Wall time | Max RSS |
|---|---:|---:|
| `cargo test --lib flank -- --nocapture` | 0:00.46 | 88064 KB |
| `c4.flank0 crush` | 0:35.35 | 2266792 KB |
| `c4.flank64 crush` | 3:06.53 | 2504012 KB |

The `flank64` wall time includes stderr compaction overhead for bounded
diagnostics. The underlying first attempt without inline compaction finished
faster but produced an unacceptable 18 GB stderr file; the final reported run
is the bounded WG-tracked run.

## Local Graph Artifacts

| Graph | GFA | zst |
|---|---|---|
| `sanity.flank0` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/sanity.flank0.gfa` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/sanity.flank0.gfa.zst` |
| `sanity.flank4` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/sanity.flank4.gfa` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/sanity.flank4.gfa.zst` |
| `c4.flank0` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank0.gfa` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank0.gfa.zst` |
| `c4.flank64` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank64.gfa` | `/home/erikg/impg/data/flank_aware_crush_c4_20260606T112300Z/graphs/c4.flank64.gfa.zst` |

## Conclusion

Flank-aware extraction/trimming works and is not relying on unrelated filtering:
it preserves paths exactly, reports requested/actual flanks, exercises path-end
shortening and reverse-orientation logic in sanity coverage, and rejects no C4
candidates through hidden quality guards.

For C4 graph quality, the answer is no for the observed local
underalignment/splitting fix. `flank64` gives a modest compaction win and more
accepted candidates, but the large residual spans and white-space proxies that
track the visible splitting are essentially unchanged. The next useful
measurement should inspect the specific residual 1312/1298-step C4 bubbles in
the `flank64` graph and determine whether the remaining split is caused by
resolver behavior, candidate discovery boundaries, or a separate normalization
issue.
