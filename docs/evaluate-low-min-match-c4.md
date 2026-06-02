# C4 Low Min-Match Top-Flubble Evaluation

Task: `evaluate-low-min`

Date: 2026-05-28 UTC

Output root:

```text
/home/erikg/impg/data/c4_low_min_match_20260528T163541Z
```

Uploaded PNG:

```text
https://hypervolu.me/~erik/impg/c4-low-min-default-off-poasta-smooth.png
```

Confirmed upload:

```text
ssh erik@hypervolu.me 'ls -lh www/impg/c4-low-min-default-off-poasta-smooth.png'
-rw-r--r-- 1 erik erik 2.0M May 28 17:29 www/impg/c4-low-min-default-off-poasta-smooth.png
```

## Summary

The `min_match_len=299` observed in the earlier top-flubble SweepGA path was not inherited from seqwish itself. It came from impg's wrapper-side adaptive rescue: the replacement path started from `seqwish-k=311`, then lowered the floor to the shortest local traversal length. In the reproduced fixture the shortest traversal was 299 bp, but every indexable exact CIGAR run was shorter than 299, so seqwish received non-empty PAF but indexed no shared runs.

That was an accidental coupling between local replacement induction and the whole-graph `seqwish-k` default. This change removes the hidden coupling by default. Local top-flubble/SweepGA replacement induction now defaults to `min_match_len=1`, which effectively disables exact-run filtering. The old adaptive behavior remains available only when explicitly requested with `min-match-length=adaptive`. Fixed explicit values are also supported, for example `min-match-length=31` or `min-match-length=63`.

On real C4, all threshold variants completed, preserved 465/465 paths, and had 0 path spelling mismatches against the no-crush reference. Lowering the threshold did not cause additional zero-PAF top-level candidates to align; the same 25/484 top-level replacements were applied and 459 zero-PAF candidates were skipped in every threshold run. However, disabling the floor increased seqwish shared-node evidence inside the 25 PAF-bearing blocks.

The second-stage POASTA smoothing path works mechanically after top-level SweepGA: it ran 4 rounds, resolved 1781 smaller bubbles, had 0 bails, preserved all 465 paths, and reduced segment bp from 797,812 to 570,842. It is not yet a clear quality win: the resolver score, whitespace totals, and trivial-stringy heuristic worsened after POASTA smoothing.

## Code Locations

The relevant local wrapper behavior is in these files:

| Area | Location | Notes |
| --- | --- | --- |
| User-facing policy | `src/main.rs:2353-2368`, `src/main.rs:2713-2723`, `src/main.rs:5119-5126` | Parses `--min-match-length` / `min-match-length=` as `off`, `adaptive`, `0`, `1`, or a fixed bp value. |
| Legacy base value | `src/resolution.rs:356` | `DEFAULT_REPLACEMENT_SEQWISH_MIN_MATCH_LEN = 311`; now only used by explicit adaptive policy unless callers set a fixed policy. |
| Replacement policy state | `src/resolution.rs:105-118`, `src/resolution.rs:208-212`, `src/resolution.rs:407-408` | Adds `ReplacementMinMatchLenPolicy`; default is `Fixed(1)` so local exact-run filtering is off. |
| Wrapper to graph config | `src/resolution.rs:4218-4263` | `seqwish_replacement_config()` maps default/off to `min_match_len=1`, `adaptive=false`; maps `adaptive` to `min_match_len=311`, `adaptive=true`; maps fixed values directly. |
| Adaptive clamp | `src/syng_graph.rs:981-1019` | Computes shortest traversal and PAF exact-run stats, but only lowers `min_match_len` when `adaptive_min_match_len=true`. |
| Debug evidence | `src/syng_graph.rs:1063-1099` | Writes per-block `summary.tsv`, now including `adaptive_min_match_len` and `effective_min_match_len`. |
| Final call into seqwish | `src/commands/graph.rs:227-234` | Calls `seqwish::alignments::unpack_paf_alignments(..., config.min_match_len, ...)`. |
| Dependency exact-run filter | `~/.cargo/git/checkouts/seqwish-be7416204a0fe48d/15aee55/src/alignments.rs:141` | Seqwish indexes CIGAR `M` / `=` runs only when `len >= min_match_len`. The wrapper chooses this value. |
| PAF prefilter | `src/commands/graph.rs:675-724` | `filter_generated_paf()` applies SweepGA post-alignment filtering unless `no_filter` is set; replacement `min_map_length=0` follows the effective min-match length. |

The fixture that exposes the 299-vs-311 problem is:

```text
tests/test_data/crush/top_flubble_seqwish_minrun.paf
tests/test_data/crush/top_flubble_seqwish_minrun.fa
```

The CIGARs in that fixture include exact runs like `47=`, `20=`, and `51=`, while the affected path had length 299. An adaptive shortest-traversal clamp from 311 to 299 still blocked those real exact runs. With the new default/off policy, the same block uses `effective_min_match_len=1` and induces shared graph nodes.

## CLI Behavior

For `impg crush`, the new user-facing knob is:

```text
--min-match-length off|adaptive|0|1|<bp>
```

Engine strings can use the same policy:

```text
gfa:syng:mask,min-run=3:crush,method=top-flubble-sweepga,seqwish-k=311,min-match-length=63:nosort
```

The old `--seqwish-k` / `seqwish-k=` option is retained as the adaptive base. It no longer opts the local replacement path into exact-run filtering by itself.

## C4 Commands

The reference no-crush run wrote:

```text
ref/run.nosort.gfa
ref metrics: paths=465 segments=18048 segment_bp=389354
```

The top-level threshold runs used this base engine:

```text
gfa:syng:mask,min-run=3:crush,method=top-flubble-sweepga,aligner=wfmash,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,polish-rounds=0:nosort
```

Threshold variants:

```text
default_off: no min-match-length parameter; default policy is off/effective 1
fixed31:     min-match-length=31
fixed63:     min-match-length=63
adaptive:    min-match-length=adaptive
```

The smoothing run used the `default_off` top-level result:

```text
target/release/impg crush \
  --gfa default_off/run.nosort.gfa \
  --output default_off_poasta_smooth/run.nosort.gfa \
  --method poasta \
  --max-iterations 5 \
  --max-traversal-len 10000 \
  --max-median-traversal-len 1000 \
  --max-total-sequence 1000000 \
  --max-traversals 10000 \
  --poa-scoring 1,4,6,2,26,1 \
  --threads 64 \
  --verbose 1
```

## Threshold Results

All top-level threshold runs:

- completed with exit status 0;
- preserved 465/465 paths;
- had 0 path spelling mismatches against `ref/run.nosort.gfa`;
- had 25 per-block debug summaries;
- had 25/25 nonzero-PAF graph builds;
- applied 25/484 top-level replacements and skipped 459 zero-PAF candidates.

| Run | Effective min | Applied/skipped | Nonzero PAF blocks | Seqwish shared segment sum | Seqwish S sum | Seqwish bp sum | Final S | Final bp | Trivial-stringy | Wall | Post score | ws-total |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| default/off | 1 | 25/459 | 25 | 2635 | 4556 | 433406 | 21540 | 797812 | 491 | 9:55.45 | 198731988 | 20688581799 |
| fixed31 | 31 | 25/459 | 25 | 1550 | 3388 | 440618 | 20372 | 805024 | 478 | 9:45.24 | 205309917 | 19740686157 |
| fixed63 | 63 | 25/459 | 25 | 1420 | 3248 | 444011 | 20232 | 808417 | 476 | 9:34.57 | 207417581 | 19729823343 |
| adaptive | 70-265 | 25/459 | 25 | 1304 | 3134 | 458974 | 20118 | 823380 | 468 | 9:34.55 | 207863512 | 19993037041 |

Interpretation:

- The low/off threshold does not make the 459 zero-PAF candidates align. Those are skipped before seqwish can help.
- Among the 25 blocks that do have PAF, default/off creates the most shared-node evidence (`replacement_shared_segments` sum 2635 vs 1304 adaptive).
- Adaptive remains the most compact by segment count and trivial-stringy count, but it does so through an explicit opt-in exact-run floor now.
- Fixed 63 is a reasonable explicit local value if compactness is weighted more than shared-node evidence, but it is still slightly worse than adaptive on segment count and trivial-stringy.
- Default/off is the correct semantic default because it avoids a hidden exact-run floor and decouples local replacement induction from whole-graph `seqwish-k=311`.

## Smoothing Result

The second-stage POASTA smoothing pass on `default_off` completed:

```text
crush: 1781 resolved, 0 bailed, 1781 candidates seen across 4 rounds
Elapsed wall time: 1:20.58
paths=465 segments=21346 segment_bp=570842
spelling mismatches vs no-crush reference: 0
```

Round summary:

| Round | Selected | Resolved | Segment change | Segment-bp change | Score change |
| ---: | ---: | ---: | --- | --- | --- |
| 1 | 741 | 741 | 21540 -> 21116 | 797812 -> 568634 | 198731988 -> 218854547 |
| 2 | 950 | 950 | 21116 -> 21251 | 568634 -> 569569 | 218854547 -> 223534121 |
| 3 | 81 | 81 | 21251 -> 21334 | 569569 -> 570748 | 223534121 -> 225574699 |
| 4 | 9 | 9 | 21334 -> 21346 | 570748 -> 570842 | 225574699 -> 227333914 |
| 5 | 0 | 0 | no eligible candidates | no eligible candidates | no eligible candidates |

The smoothing path therefore works and converges, but it is mixed:

- Good: segment bp drops by 226,970 bp relative to the default/off top-level run.
- Good: segment count drops slightly, 21,540 -> 21,346.
- Bad: the quality score worsens, 198,731,988 -> 227,333,914.
- Bad: whitespace total worsens, 20,688,581,799 -> 27,762,234,445.
- Bad: trivial-stringy candidates increase, 491 -> 550.

Conclusion: iterative POASTA smoothing is runnable and path preserving after top-level SweepGA, but this parameterization does not yet improve local under-alignment quality.

## Validation

Commands and outcomes:

```text
cargo test --test test_crush_integration c4_fragment_seqwish_regressions_induce_shared_graphs -- --nocapture
PASS

cargo test --test test_crush_integration c4_top_flubble_seqwish_indexes_observed_exact_runs -- --nocapture
PASS

cargo test --test test_crush_integration c4_fragment_lacing_uses_pairwise_induced_graphs -- --nocapture
PASS

cargo test --test test_syng_startcount -- --test-threads=1
PASS

cargo test --all -- --test-threads=1
PASS

cargo install --path .
PASS
```

The full-suite run initially exposed an order-sensitive native crash in `tests/test_syng_startcount.rs`: `test_identical_sequences_get_distinct_start_counts` removed its temp directory while the loaded `SyngIndex` values were still in scope. Under the default Rust test capture mode, the following FFI round-trip test could crash or trip `checkSymRuns`. The test now explicitly drops `loaded` and `index` before removing the directory; the standalone startcount binary and full `cargo test --all -- --test-threads=1` pass afterward.

## Conclusion

Answers to the task questions:

1. The effective 299 threshold was an impg wrapper artifact, not seqwish. `seqwish-k=311` was being lowered by local adaptive wrapper code to the shortest local traversal length, 299 bp in the reproduced block. Seqwish then applied its normal exact-run filter to CIGAR runs, and every relevant exact run was shorter than 299.
2. Removing the hidden threshold does not make additional zero-PAF C4 top-level candidates align. It does make every PAF-bearing C4 top-level block produce shared graph nodes and increases the shared-node evidence within those blocks.
3. The metric effect is mixed. Default/off gives the strongest shared-node evidence and best recorded resolver score for the top-level pass, but worse segment count and trivial-stringy counts than adaptive. Explicit `31` and `63` sit between off and adaptive. All runs preserve paths exactly.
4. The iterative smaller-bubble POASTA smoothing path runs after top-level SweepGA and converges, but this configuration does not yet improve local under-alignment quality. It preserves spelling and reduces segment bp, but worsens whitespace/score and the trivial-stringy heuristic.

The operational fix is still clear: local replacement exact-run filtering must default off. Any nontrivial local min-match floor must be an explicit opt-in policy, not an implicit consequence of `seqwish-k=311`.
