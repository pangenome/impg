# Validate Current C4 Local Graph Quality

Task: `validate-current-c4`

Date: 2026-06-08 UTC

Branch/worktree head: `e3c2494` (`feat: resolve-merge-for (agent-629)`)

## Decision

The current C4 result is **better than prior SYNG-derived C4 outputs on the
specific repeat-artifact failure mode, but still not solved and not yet
clearly PGGB-comparable**.

The important improvement is that the current branch-local reproduction has
zero direct self-loop edges and zero adjacent same-node path steps. Prior C4
k311 reports still had 57 to 60 direct self-loop edges and about 252k to 264k
adjacent same-node path steps after comparable k311 seed or POA cleanup runs.

The remaining blocker is graph quality relative to the PGGB/SmoothXG-style
control. The reproduced SYNG k311 + POA 1kb output still has 2.59x more
singleton bp than the PGGB control and 10x worse path white-space p99, even
though it has fewer path steps, fewer white-space bridges, and a better
white-space max. This is a real improvement over prior SYNG, not convergence.

## Recipe Selection

The resolved local-compression docs identify
`chunk_window_sweepga_seqwish` as the best current fast-profile signal, but
they also say that row is a resolver-distinct synthetic fixture generator and
is not C4-ready. I did not restart the autopoietic loop or run a broad search.

For actual full-C4 validation I used the latest validated real-C4 recipe from
the resolved work lineage:

- Reproduce the full-C4 k311 SweepGA/seqwish seed construction from the
  existing C4 all-vs-all PAF evidence.
- Apply the conservative `--method poa` 1kb cleanup recommended by the k311
  POA threshold sweep.
- Compare against the existing full-C4 PGGB/SmoothXG-style control from
  `/home/erikg/impg/data/c4_pggb_control_20260526T025439Z/`.

No standalone full-C4 SmoothXG-only control was available in the resolved docs.
The PGGB control is the best available PGGB/SmoothXG-style control; its query
log shows the internal PGGB pipeline with smoothing and `gfaffix`.

## Run Root

All generated graph, render, report, and log artifacts were written outside
git:

```text
/home/erikg/impg/data/validate_current_c4_20260608T125053Z
/home/erikg/impg/data/validate_current_c4_latest -> validate_current_c4_20260608T125053Z
```

Only this report and the small metrics TSV/JSON are committed.

## Commands

The branch-local binaries were built first:

```bash
source ./env.sh
cargo build --release --bin impg --example compare_gfa_paths
cargo build --release --bin gfaffix
```

The first seed reproduction attempt reached graph emission but failed at the
final `gfaffix` step because `target/release/gfaffix` was not yet built next
to `target/release/impg`:

```text
Exit status: 1
Wall: 3:14.10
Max RSS: 8,744,752 KB
Blocker: gfaffix binary not found next to impg executable
```

After building `gfaffix`, the same seed construction succeeded:

```bash
target/release/impg graph \
  --sequence-files /home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/c4_whole_region.fa \
  --paf-file /home/erikg/impg/data/c4_low_seqwish_k_20260605T140000Z/debug/k79/raw.paf \
  --gfa-engine seqwish \
  --fastga \
  --num-mappings 1:many \
  --scaffold-filter 1:many \
  --scaffold-jump 0 \
  --min-match-len 311 \
  --temp-dir /home/erikg/impg/data/validate_current_c4_latest/tmp/c4.k311.seed \
  --debug-dir /home/erikg/impg/data/validate_current_c4_latest/debug/c4.k311.seed \
  -g /home/erikg/impg/data/validate_current_c4_latest/graphs/c4.k311.seed.repro.gfa \
  -t 32 -v 1
```

POA cleanup:

```bash
target/release/impg crush \
  -g /home/erikg/impg/data/validate_current_c4_latest/graphs/c4.k311.seed.repro.gfa \
  -o /home/erikg/impg/data/validate_current_c4_latest/graphs/c4.k311.poa1kb.repro.gfa \
  --method poa \
  --max-iterations 5 \
  --max-traversal-len 1k \
  --max-median-traversal-len 1k \
  --max-total-sequence 1m \
  --max-traversals 10k \
  -t 32 -v 1
```

Then I ran `compare_gfa_paths`, `gfasort -p Ygs`, `gfalook -m`, and
`target/release/impg graph-report --format tsv --povu` on the reproduced SYNG
outputs and the existing PGGB control.

## Metrics

Full fields are in:

- `docs/evaluations/validate-current-c4-metrics.tsv`
- `docs/evaluations/validate-current-c4-metrics.json`

Key metrics:

| graph | segments | path steps | segment bp | singleton bp | replay ratio | depth median/p95/max | white-space p99/max | white-space bridges | self-loop edges | adjacent same-node steps | runtime | RSS KB |
| --- | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| SYNG k311 seed repro | 7,346 | 2,865,437 | 253,140 | 5,925 | 421.464304 | 461 / 931 / 16,043 | 199 / 36,423 | 5,958 | 0 | 0 | 3:18.19 | 8,716,128 |
| SYNG k311 + POA 1kb repro | 8,142 | 3,065,010 | 258,628 | 7,492 | 412.520972 | 461 / 931 / 16,043 | 160 / 25,248 | 4,900 | 0 | 0 | 0:33.07 | 2,200,168 |
| PGGB/SmoothXG-style control | 13,288 | 5,538,879 | 234,524 | 2,890 | 454.919215 | 465 / 931 / 7,079 | 16 / 73,181 | 22,851 | 1 | 51 | 13:37.90 | 64,330,340 |

Coverage and long-link proxies:

| graph | node coverage mean | bp-weighted coverage | link jump p99/max | path jump p99/max |
| --- | ---: | ---: | ---: | ---: |
| SYNG k311 seed repro | 390.067656 | 421.464304 | 11 / 2,012 | 5 / 2,012 |
| SYNG k311 + POA 1kb repro | 376.444363 | 412.520972 | 8 / 1,130 | 6 / 1,130 |
| PGGB/SmoothXG-style control | 416.833158 | 454.919215 | 5 / 3,488 | 5 / 3,488 |

Interpretation:

- The reproduced SYNG k311 + POA 1kb graph has much better local repeat
  artifact diagnostics than prior SYNG C4 outputs: 0 direct self-loop edges
  and 0 adjacent same-node path steps.
- The POA cleanup improves white-space max and bridge count relative to the
  reproduced seed, but increases segments, path steps, segment bp, singleton
  bp, and duplicate sequence fraction.
- Relative to PGGB, the SYNG + POA graph is still not a clean win. It has fewer
  segments and path steps and far fewer white-space bridges, but worse singleton
  bp, lower replay ratio, worse p99 white-space, and a much higher path-depth
  max.
- The PGGB control has its own rough edges: higher path steps, more white-space
  bridges, one direct self-loop edge, and 51 adjacent same-node path steps.
  Still, its singleton bp and p99 white-space remain substantially better.

## Path Checks

SYNG POA exact path preservation passed against the reproduced SYNG seed:

```text
expected_paths      465
observed_paths      465
missing_paths       0
extra_paths         0
spelling_mismatches 0
```

The PGGB control could not be checked by direct full-name equality against the
SYNG seed because all 465 path names differ between the two constructions:

```text
expected_paths      465
observed_paths      465
missing_paths       465
extra_paths         465
spelling_mismatches 0
```

This is a name-stability incompatibility between the two C4 control artifacts,
not a shared-name spelling mismatch.

## Uploaded Artifacts

All URLs returned HTTP 200 on 2026-06-08.

| artifact | URL |
| --- | --- |
| SYNG k311 seed PNG | http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.seed.repro.Ygs.png |
| SYNG k311 + POA 1kb PNG | http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.poa1kb.repro.Ygs.png |
| PGGB control PNG | http://hypervolu.me/~erik/impg/c4.validate-current-c4.pggb-control.Ygs.png |
| SYNG k311 seed GFA zst | http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.seed.repro.gfa.zst |
| SYNG k311 + POA 1kb GFA zst | http://hypervolu.me/~erik/impg/c4.validate-current-c4.k311.poa1kb.repro.gfa.zst |

## Validation Notes

No source code was changed. This is validation-only work, so I did not run the
full test suite. The relevant validation was the branch-local release build,
the real C4 reproduction, exact path check for the reproduced SYNG cleanup,
fresh graph reports, render/upload checks, and repository artifact hygiene.

The committed files are small documentation and metrics files only. Large
GFAs, sorted GFAs, PNGs, debug PAFs, and command logs remain outside git under
`/home/erikg/impg/data/validate_current_c4_latest`.
