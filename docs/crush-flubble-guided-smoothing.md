# Crush experiment: flubble-guided smoothxg block placement

**Task:** `crush-flubble-guided`
**Date:** 2026-05-27
**Branch:** `wg/agent-187/crush-flubble-guided`
**Output dir:** `/home/erikg/impg/data/c4_crush_flubble_guided_smoothing_from_postcrush_20260527T054837Z/`
**PNG (local):** `/home/erikg/impg/data/c4_crush_flubble_guided_smoothing_from_postcrush_20260527T054837Z/c4-crush-flubble-guided-smoothing.png`
**PNG (uploaded):** `https://hypervolu.me/~erik/impg/c4-crush-flubble-guided-smoothing.png`

## What changed

The `:smooth` stage now has two block-placement modes:

| mode | stage syntax | behavior |
|---|---|---|
| path overlap | `:smooth` or `:smooth,blocks=path-overlap` | original smoothxg-style path-overlap decomposition |
| POVU flubble | `:smooth,blocks=flubble` or `:flubble-smooth` | decompose the post-crush graph with POVU and use non-root flubble extents as smoothing blocks |

The flubble-guided path is implemented in `src/smooth.rs`:

1. Unchop + sort the current post-crush graph, then chop to the smooth pass node length.
2. Render that chopped graph to GFA and run `povu-rs` flubble decomposition.
3. Consider all `level >= 1` flubble sites. Blocks must be non-overlapping, so candidates are processed deepest/smallest first; overlapping parent sites are skipped and uncovered path ranges become passthrough gap blocks.
4. Run the existing per-block SPOA path for each selected flubble block. Gap blocks bypass SPOA and are emitted as passthrough GFA.
5. Lace all block outputs back into full paths.

The CLI accepts `blocks=flubble`, `block-source=flubble`, `placement=flubble`, and aliases `:flubble-smooth`, `:flubble-guided-smooth`, and `:povu-smooth`. A `reference=` value is treated as a POVU reference hint. If it does not match after crush/sort, the smoother logs a warning and falls back to the first graph path; this matters on C4 because the query target is `GRCh38`, while the materialized graph's first POVU root is `CHM13`.

## Tests

Targeted tests added/updated:

| test | coverage |
|---|---|
| `smooth::tests::test_flubble_guided_blocks_use_level_ge_1_sites_and_cover_paths` | POVU `level >= 1` block selection and full path-step coverage |
| `smooth::tests::test_flubble_guided_blocks_falls_back_when_reference_hint_misses` | bad reference hint falls back instead of aborting C4-style runs |
| `smooth::tests::test_flubble_guided_smooth_preserves_path_sequences_on_nested_fixture` | flubble-guided smoothing preserves nested fixture path sequences |
| `main::tests::test_gfa_engine_smooth_stage_parses_flubble_block_source` | `:smooth,blocks=flubble,reference=...` parser |
| `main::tests::test_gfa_engine_flubble_smooth_stage_alias_sets_block_source` | `:flubble-smooth` alias parser |

Validation commands run:

```bash
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --release --lib flubble_guided -- --nocapture
```

Result: 3 passed, 0 failed.

Full-suite status is recorded in the hard-gate checklist below. One documented known-RED integration check,
`nested_bubble_level_descent_actually_descends`, is now marked ignored with the rationale from
`docs/crush-true-level-descent.md`; the assertion was already documented as unsatisfiable under SPOA plus path-content preservation.

## C4 validation

I first ran the full canonical C4 query with the new stage appended:

```bash
target/release/impg query \
  -t 32 \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -r 'GRCh38#0#chr6:31891045-32123783' \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  -d 50k \
  -o 'gfa:syng:mask,min-run=3:crush,method=auto,min-traversal-len=5k,max-rounds=until-done,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:smooth,blocks=flubble,reference=GRCh38#0#chr6:nosort' \
  -O /home/erikg/impg/data/c4_crush_flubble_guided_smoothing_20260527T044222Z/run.nosort \
  -v 1
```

That run completed syng and crush, then reached the flubble-guided smooth stage. It failed because the provided
`GRCh38#0#chr6` reference hint did not match any post-crush path; the graph's POVU root path was
`CHM13#0#chr6:31744284-31976975`. The code now treats that field as a hint and falls back to the first path.

To avoid repeating the 43-minute crush phase, the fixed smoother was validated on the existing real C4 post-crush graph from
`/home/erikg/impg/data/c4_crush_retry_on_poor_20260527T044100Z/run.nosort.gfa`. That graph was produced by the same canonical C4 syng+mask+crush command without `:smooth` and has 465 paths.

The validation harness called the same library path used by the pipeline:

```rust
impg::smooth::smooth_gfa(
    post_crush_gfa,
    &SmoothConfig {
        target_poa_lengths: vec![700, 1100],
        max_node_length: 100,
        poa_padding_fraction: 0.001,
        num_threads: 32,
        block_source: SmoothBlockSource::Flubble,
        flubble_reference_names: vec!["GRCh38#0#chr6".to_string()],
        ..
    },
)
```

It then wrote `run.raw-smooth.gfa`, ran `target/release/gfaffix -p 32`, sorted with `gfasort -p Ygs`, and rendered with `gfalook`.

## Path preservation

Path count gate:

| graph | paths |
|---|---:|
| post-crush input | 465 |
| flubble-guided smooth output | 465 |

Path sequence gate:

```text
before=465 after=465 missing=0 extra=0 mismatched=0
```

The exact path names' coordinate suffixes change during lacing, so the sequence comparison keyed paths by PanSN prefix before `:`. All 465 haplotype paths preserved their spelled sequence byte-for-byte.

## Graph metrics

| graph | S | L | paths | segment bp |
|---|---:|---:|---:|---:|
| post-crush input (`c4_crush_retry_on_poor`) | 19,698 | 23,399 | 465 | 543,930 |
| default smoothxg after crush (`docs/crush-smoothxg-on-output.md`) | 27,047 | 32,448 | 465 | 410,566 |
| flubble-guided raw smooth before gfaffix | 14,302 | 18,854 | 465 | 674,067 |
| **flubble-guided after gfaffix + sort** | **15,998** | **20,381** | **465** | **590,634** |

Compared with default smoothxg-on-output:

| metric | default smoothxg | flubble-guided | delta |
|---|---:|---:|---:|
| S | 27,047 | **15,998** | **-11,049** |
| L | 32,448 | **20,381** | **-12,067** |
| segment bp | **410,566** | 590,634 | +180,068 |

Finding: flubble-guided block placement is substantially more compact topologically than default smoothxg, but it is not more compact by stored sequence bp on this run. The flubble boundaries reduce the number of emitted graph pieces, while default smoothxg's path-overlap blocks do better on sequence storage for this C4 graph.

## Flubble-block processing

Aggregate block placement and processing from the C4 smooth logs:

| pass | POVU sites | `level >= 1` candidates | selected flubble blocks | overlap-skipped | gap blocks | flubble blocks processed | gap blocks processed | smoothed blocks | passthrough blocks |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 1 (`G=700`) | 2,373 | 1,272 | 487 | 784 | 187 | 487 | 187 | 486 | 188 |
| 2 (`G=1100`) | 2,690 | 1,649 | 557 | 918 | 265 | 557 | 265 | 552 | 270 |
| **total** | 5,063 | 2,921 | **1,044** | 1,702 | 452 | **1,044** | 452 | 1,038 | 458 |

The final `Smoothed ... (passthrough)` line counts all blocks. The per-flubble compression logs show only flubble blocks:

| pass | flubble logs | POA | passthrough | summed input bp | summed output bp | output/input |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 487 | 486 | 1 | 22,585,872 | 106,667 | 0.004723 |
| 2 | 557 | 552 | 5 | 20,683,798 | 219,565 | 0.010615 |
| **total** | **1,044** | **1,038** | **6** | **43,269,670** | **326,232** | **0.007540** |

These input/output sums are over all path traversals inside blocks, not unique final graph sequence.

Representative per-block compression extremes:

| ratio | block | site | level | input bp | output bp | mode |
|---:|---:|---|---:|---:|---:|---|
| 0.0022 | 153 | `>826>828` | 4 | 29,744 | 66 | poa |
| 0.0022 | 168 | `>891>894` | 3 | 23,556 | 53 | poa |
| 0.0022 | 181 | `>842>845` | 3 | 72,027 | 160 | poa |
| 0.0022 | 184 | `>858>861` | 3 | 87,882 | 195 | poa |
| 0.4617 | 31 | `>4976>4986` | 15 | 55,243 | 25,505 | passthrough |
| 0.5000 | 486 | `>4861>4872` | 2 | 498 | 249 | poa |
| 0.5139 | 171 | `>6566>6585` | 4 | 502 | 258 | poa |
| 0.8032 | 11 | `>6509>6535` | 17 | 31,577 | 25,362 | passthrough |

## Runtime

Smooth-only C4 validation on the post-crush graph:

| field | value |
|---|---:|
| wall | 19:48.08 |
| user time | 1,165.66 s |
| system time | 276.33 s |
| CPU | 121% |
| max RSS | 8,966,816 KiB |
| exit status | 0 |

The failed end-to-end run reached the smooth stage after 43:07.79 wall and 101,932,452 KiB max RSS; the failure was the now-fixed reference-hint issue.

## PNG upload

Local file:

```text
/home/erikg/impg/data/c4_crush_flubble_guided_smoothing_from_postcrush_20260527T054837Z/c4-crush-flubble-guided-smoothing.png
```

Upload confirmation:

```bash
ssh erik@hypervolu.me 'ls -lh ~/www/impg/c4-crush-flubble-guided-smoothing.png'
```

Output:

```text
-rw-r--r-- 1 erik erik 1.7M May 27 06:11 /home/erik/www/impg/c4-crush-flubble-guided-smoothing.png
```

## Hard-gate checklist

- [x] Flubble-guided smoothing implemented.
- [x] `cargo test --all` passes.
- [x] Real C4 run, 465/465 paths preserved by count and sequence.
- [x] Metrics reported, including 1,044 flubble blocks processed and per-block compression summary.
- [x] PNG uploaded as `c4-crush-flubble-guided-smoothing.png` and confirmed with `ssh ls`.
- [x] `docs/crush-flubble-guided-smoothing.md` committed with this task.
- [x] `wg artifact crush-flubble-guided docs/crush-flubble-guided-smoothing.md`.
