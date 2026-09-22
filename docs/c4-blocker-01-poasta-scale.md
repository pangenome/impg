# C4 Blocker 01: Poasta Residual Scale

Task: `c4-blocker-01-poasta-scale-impl`

## Summary

The C4 integrated run was not failing path correctness. It was selecting the
wrong residual scale under `--window-mode largest --candidate-limit 1`: after
the local name/wfmash fixes, later rounds repeatedly selected a small child
residual whose replacement was path-correct but near identity. The selected
child moved little or no singleton sequence while broader enclosing residuals
were still visible in the same POVU discovery pass.

`src/resolution.rs` now ranks largest-mode iterative multi-level candidates by
a residual-scale score before the old balanced-frontier proxy:

```text
root_span_bp * max(unique_steps, 1)
```

This only changes candidate priority and built-candidate tie-breaking. It does
not add a graph-quality safety gate, and exact path preservation remains the
hard acceptance invariant.

## Evidence

The old integrated C4 log selected a small child in round 8:

```text
>272218192..>272218467 root_span=6733bp traversals=132 max=6734 median=6734 total=825193 unique_steps=111
```

The same discovery pass showed a broader residual that the previous ranking did
not admit under the one-candidate cap:

```text
>272218968..>272219341 root_span=32776bp traversals=143 max=26408 median=38 total=348244 unique_steps=1111
```

The patched one-round representative full-C4 selection run started from the old
integrated graph and selected the broader residual as candidate `#1`:

```text
generated candidate detail: #1 source=top-level sites=1 level=0 >272218968..>272219341 root_steps=4827..6071 root_span=32776bp traversals=143 max=26408 median=38 total=348244 unique_steps=1111 step_savings=13960
```

Artifact:

```text
/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/after/crush.stderr.log
```

That representative run was stopped during the broad POA build after discovery
and selection had completed. It was not used for after-compression metrics.

## Metrics

Baseline and PGGB control were regenerated with:

```bash
target/debug/impg graph-report \
  -g /home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfa \
  -o /home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/before/c4_full_integrated.graph-report.tsv \
  --format tsv

target/debug/impg graph-report \
  -g /home/erikg/impg/data/c4_pggb_control_20260526T025439Z/pggb.nosort.gfa \
  -o /home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/before/pggb_control.graph-report.tsv \
  --format tsv
```

| graph | segments | total segment bp | singleton bp | bp-weighted coverage | path white-space p99 | path white-space max | PNG |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| integrated C4 before fix | 19,026 | 462,940 | 78,435 | 230.460695 | 160,890 | 459,430 | `/home/erikg/impg/data/integrate_c4_local_20260602T2020Z/full_c4/c4_full_integrated.gfalook-m.png` |
| PGGB control | 13,288 | 234,524 | 2,890 | 454.919215 | 14 | 219,917 | n/a |
| saved residual input | 446 | 15,819 | 1,202 | 62.443454 | 3,678 | 9,135 | `/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/before/residual_input.gfalook-m.png` |
| saved residual Poasta output | 102 | 6,763 | 13 | 146.058406 | 6,640 | 6,640 | `/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/after/residual_poasta_output.gfalook-m.png` |

The saved residual slice confirms that selecting the residual as one unit gives
a path-preserving compression win:

```text
segments:       446 -> 102
segment bp:     15,819 -> 6,763
singleton bp:   1,202 -> 13
bp-weighted coverage: 62.443454 -> 146.058406
```

Path-spelling validation for the saved residual:

```text
expected_paths        333
observed_paths        333
missing_paths         0
extra_paths           0
spelling_mismatches   0
```

Artifact:

```text
/home/erikg/impg/data/c4_blocker_01_poasta_scale_20260603T0300Z/after/residual_path_validation.tsv
```

## Validation

Commands run:

```bash
source ./env.sh && cargo test --lib residual_scale -- --nocapture
source ./env.sh && cargo test --lib resolution:: -- --nocapture
source ./env.sh && cargo build
source ./env.sh && cargo test --bin impg -- --nocapture
source ./env.sh && cargo test --lib

source ./env.sh && cargo run --example compare_gfa_paths -- \
  /home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/site_272218192_272218467/input.gfa \
  /home/erikg/impg/data/settle_local_replacement_20260602T190213Z/poasta_residual/site_272218192_272218467/poasta.output.gfa
```

Results:

```text
cargo test --lib residual_scale: 2 passed
cargo test --lib resolution::: 91 passed
cargo build: passed
cargo test --bin impg: 81 passed
cargo test --lib: 365 passed
compare_gfa_paths: 333/333 paths, 0 spelling mismatches
```

Local build notes:

- `vendor/syng` and `vendor/gfaffix` submodules were initialized in the WG
  worktree for validation.
- `wfmash-rs` attempted to rebuild vendored wfmash; the local build used the
  existing `/home/erikg/.cargo/bin/wfmash` binary in Cargo's `wfmash-rs`
  `OUT_DIR`, matching the workaround already documented for prior C4 tasks.

## Next Blocker

The scale selection bug is fixed, but the first broader full-C4 residual chosen
by the patched ranking routes to `Poa` because routing still keys primarily off
median traversal length. That broad residual has `median=38` but
`max=26408` and `total=348244`, and the representative full-C4 run spent
several minutes in the POA build before it was stopped.

The next blocker is therefore residual build routing, not graph-quality
admission: broad high-scale residual windows with long outlier traversals should
route to a scalable builder such as Poasta or SweepGA instead of the small-POA
path, while preserving exact path-sequence validation and semantic
FASTA/PAF/GFA names.
