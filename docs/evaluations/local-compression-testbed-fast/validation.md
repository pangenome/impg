# Local Compression Testbed Fast Validation

Date: 2026-06-07

This artifact records the validation used for the PGGB and SmoothXG control
runner update.

## Commands

- `python3 -m py_compile scripts/local_compression_testbed.py`: passed.
- `python3 scripts/local_compression_testbed.py write-fixtures --root tests/test_data/local_compression`: passed.
- `python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast`: passed after installing the current worktree binary; wrote 143 rows.
- `source ./env.sh && export CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" && cargo test --test test_local_compression_testbed`: passed, 2 tests.
- `source ./env.sh && export CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" && cargo build`: passed.
- `source ./env.sh && export CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" && cargo test`: passed.
- `source ./env.sh && export CMAKE_C_COMPILER="$CC" CMAKE_CXX_COMPILER="$CXX" && cargo install --path .`: passed and replaced the local `impg` and `gfaffix` binaries from this worktree.

The first unconfigured targeted cargo test attempt failed before running tests
because the default build environment could not find the vendored htslib faidx
header. After sourcing `env.sh`, exporting the CMake compiler paths, and
initializing `vendor/syng` and `vendor/gfaffix`, the targeted test, build, full
test suite, install, and profile rerun all passed.

## Scoreboard Summary

- Rows: 143.
- Command status counts: `pass=80`, `path_corrupt=30`, `skipped=33`.
- Topology status counts: `pass=70`, `fail=10`, `not_run=63`.
- Control methods represented: `pggb_control`, `smoothxg_control`, and `pggb_plus_smoothxg_control`.
- Included CI fixtures executed 30 real control commands: 10 rows for each control method.
- Local-tier fixtures skipped 9 control rows by profile policy: 3 rows for each control method.

## Control Interpretation

All included PGGB/SmoothXG controls executed bounded local commands and wrote
command logs, stdout/stderr logs, output GFAs, normalized GFAs, metrics, and
scoreboard rows. The control rows are visible as `path_corrupt`, not hidden or
filtered: exact path preservation failed because the local PGGB path rewrote
PanSN interval-name ends from the expected fixture coordinates. Since exact
path corruption is the only hard rejection, topology checks were left
`not_run` for those rows and the failures remain present in the scoreboard.

Standalone `pggb` and `smoothxg` binaries were not found in the local PATH.
The bounded profile therefore uses the repository-supported local
`impg graph --gfa-engine pggb` path. `smoothxg_control` uses the same local
engine with an explicit empty PAF because the repository does not expose a
standalone smooth-only CLI path.
