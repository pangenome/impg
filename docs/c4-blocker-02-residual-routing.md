# C4 Blocker 02: Residual Routing

Task: `c4-blocker-02-route`

## Summary

The blocker-01 largest-mode fix correctly selected the broad full-C4 residual
window:

```text
>272218968..>272219341 root_span=32776bp traversals=143 max=26408 median=38 total=348244 unique_steps=1111
```

Before this change, iterative multi-level routing still treated that candidate
as a small bubble because the auto tier used median traversal length. It routed
to direct SPOA and stalled in the small-POA builder.

`src/resolution.rs` now keeps median-based direct auto routing intact for
ordinary candidates, but upgrades iterative multi-level windows that would have
gone to `Poa` when they also have broad residual scale:

- long max traversal at or above the scalable tier length floor;
- broad bp scale by root span or total traversal bp;
- residual path-step mass by unique step count relative to traversal count.

The default scalable direct route is `Poasta`. If the POASTA auto tier is
disabled, the same overflow path falls through to `Sweepga`. This is a dispatch
rule only. It does not use graph-quality metrics as a safety gate, and candidate
application remains guarded by exact path-sequence preservation.

## Regression Test

Added:

```text
resolution::tests::iterative_multi_level_routes_broad_residual_away_from_small_poa
```

The test uses the representative C4 shape (`median=38`, `max=26408`,
`total=348244`, `root_span=32776`, `unique_steps=1111`) and confirms two
things:

- raw median-only auto routing still demonstrates the old failure by returning
  `Poa`;
- iterative multi-level window routing upgrades the same broad residual to
  `Poasta`.

The test failed before implementation:

```text
assertion `left != right` failed: broad residual windows must route to a scalable builder
left: Poa
right: Poa
```

After implementation it passes.

## Representative C4 Slice

Output directory:

```text
/home/erikg/impg/data/c4_blocker_02_residual_routing_20260603T032503Z
```

Command:

```bash
target/debug/impg crush \
  --gfa /home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfa \
  --output /home/erikg/impg/data/c4_blocker_02_residual_routing_20260603T032503Z/full_c4_largest1_poasta.gfa \
  --method iterative-multi-level \
  --window-mode largest \
  --window-target-bp 30k \
  --max-window-sites 8 \
  --candidate-limit 1 \
  --max-iterations 1 \
  --polish-method poasta \
  --polish-rounds 0 \
  --threads 32
```

Result:

```text
crush: 1 resolved, 0 bailed, 1 candidates seen across 1 rounds
Elapsed wall time: 1:14.62
Max RSS: 1602104 kB
Exit status: 0
```

Routing evidence from `crush.stderr.log`:

```text
replacement-routing=multi-site-or-broad-residual-scalable/single-site-auto, application_gate=exact-path-preservation
generated candidate detail: #1 source=top-level sites=1 level=0 >272218968..>272219341 root_steps=4827..6071 root_span=32776bp traversals=143 max=26408 median=38 total=348244 unique_steps=1111 step_savings=13960
building candidate 1/1 source=top-level sites=1 with Poasta; traversals=143, max-len=26408, median-len=38, total-len=348244, root-span=32776
built candidate detail source=top-level sites=1 method=Poasta ... replacement_segments=153, replacement_bp=26458
```

Local candidate compression from the applied-candidate log:

| metric | input | output |
| --- | ---: | ---: |
| candidate segments | 1111 | 153 |
| candidate segment bp | 29920 | 26458 |
| candidate singleton bp | 2396 | 33 |
| candidate mean bp-weighted coverage | 11.6392 | 13.1621 |

Exact full-graph path-spelling validation:

```text
expected_paths        465
observed_paths        465
missing_paths         0
extra_paths           0
spelling_mismatches   0
```

Path validation artifact:

```text
/home/erikg/impg/data/c4_blocker_02_residual_routing_20260603T032503Z/full_c4_largest1_poasta.path-validation.tsv
```

## Graph Metrics

Graph-report artifacts:

```text
/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/before/c4_full_integrated.graph-report.tsv
/home/erikg/impg/data/c4_blocker_02_residual_routing_20260603T032503Z/full_c4_largest1_poasta.graph-report.tsv
```

| graph | segments | total segment bp | singleton bp | bp-weighted coverage | path white-space p99 | path white-space max | PNG |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| integrated C4 input | 19026 | 462940 | 78435 | 230.460695 | 160890 | 459430 | `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfalook-m.png` |
| largest residual routed to POASTA | 19162 | 488855 | 78142 | 218.243598 | 172012 | 485345 | `/home/erikg/impg/data/c4_blocker_02_residual_routing_20260603T032503Z/full_c4_largest1_poasta.gfalook-m.png` |

The representative slice is meant to validate routing and completion of the
previously stalled build, not to act as a graph-quality acceptance gate. The
global quality deltas remain diagnostic only; exact path preservation is the
hard acceptance invariant.

## Validation

Commands run:

```bash
source ./env.sh && cargo test --lib resolution::tests::iterative_multi_level_routes_broad_residual_away_from_small_poa -- --nocapture
source ./env.sh && cargo test --lib resolution:: -- --nocapture
source ./env.sh && cargo build --bin impg
source ./env.sh && cargo build --bin gfaffix
source ./env.sh && cargo test --lib
source ./env.sh && cargo test --bin impg
target/debug/impg graph-report -g full_c4_largest1_poasta.gfa -o full_c4_largest1_poasta.graph-report.tsv --format tsv
cargo run --quiet --example compare_gfa_paths -- c4_full_integrated.gfa full_c4_largest1_poasta.gfa
gfalook -i full_c4_largest1_poasta.gfa -o full_c4_largest1_poasta.gfalook-m.png -m -x 2200 -y 1200
```

Results:

```text
targeted routing test: passed
cargo test --lib resolution::: 92 passed
cargo build --bin impg: passed
cargo build --bin gfaffix: passed
cargo test --lib: 366 passed
cargo test --bin impg: 81 passed
representative C4 slice: completed, 1 resolved, path spellings preserved
```
