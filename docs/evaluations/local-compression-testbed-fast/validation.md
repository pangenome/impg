# Local Compression Testbed Fast Validation

Date: 2026-06-08

## Commands

- `python3 -m py_compile scripts/local_compression_testbed.py`: passed.
- `cargo test --test test_local_compression_testbed local_compression_path_replay_compression_ratio -- --nocapture`: passed, 1 selected test.
- `cargo test --test test_local_compression_testbed local_compression_chunk_window_sweepga_seqwish_nested_top_level_wrong -- --nocapture`: passed, 1 selected test.
- `python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast`: passed, wrote 156 rows.
- TSV assertion script from `docs/evaluations/local-compression-autopoietic/iter-2-synthesis.md`: passed; verified the new `path_replay_compression_ratio` column and the `chunk_window_sweepga_seqwish` row on `nested_top_level_wrong`.
- `cargo build`: passed with existing warnings.
- `cargo test`: passed, including `tests/test_local_compression_testbed.rs` with 6 tests.
- `cargo install --path .`: passed; replaced `/home/erikg/.cargo/bin/impg` and `/home/erikg/.cargo/bin/gfaffix` with this worktree build.

The cargo commands required the project-local native validation environment used
by prior C4 tasks:

```bash
export HTS=/home/erikg/wfmash/build/vendored_htslib
export JEM=/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0
export C_INCLUDE_PATH=$HTS/include:${C_INCLUDE_PATH:-}
export CPLUS_INCLUDE_PATH=$HTS/include:${CPLUS_INCLUDE_PATH:-}
export CPATH=$HTS/include:${CPATH:-}
export LIBRARY_PATH=$HTS/lib:$JEM/lib:${LIBRARY_PATH:-}
export LD_LIBRARY_PATH=$HTS/lib:$JEM/lib:${LD_LIBRARY_PATH:-}
export PKG_CONFIG_PATH=$HTS/lib/pkgconfig:${PKG_CONFIG_PATH:-}
export CFLAGS="-I$HTS/include ${CFLAGS:-}"
export CXXFLAGS="-I$HTS/include ${CXXFLAGS:-}"
export LDFLAGS="-L$HTS/lib -L$JEM/lib ${LDFLAGS:-}"
export CMAKE_PREFIX_PATH=$HTS:${CMAKE_PREFIX_PATH:-}
export CMAKE_INCLUDE_PATH=$HTS/include:${CMAKE_INCLUDE_PATH:-}
export CMAKE_LIBRARY_PATH=$HTS/lib:$JEM/lib:${CMAKE_LIBRARY_PATH:-}
export CARGO_BUILD_JOBS=8
export CMAKE_BUILD_PARALLEL_LEVEL=8
```

`git submodule update --init --recursive vendor/syng vendor/gfaffix` was also
required before cargo validation because this WG worktree initially had
uninitialized submodules.

## Fast Scoreboard Checks

- Rows: 156.
- Methods represented: 12, including `chunk_window_sweepga_seqwish` and all three internal controls.
- Command statuses: 132 `pass`, 24 `skipped`.
- Topology statuses: 82 `pass`, 50 `fail`, 24 `not_run`.
- Exact path preservation: 132 `pass`, 24 `not_run`.
- Hard path corruption rows: 0.
- Internal controls are present and graph-producing for all 11 CI fixtures; the two remaining local-tier fixtures are skipped by fast-profile policy only.
- `path_replay_compression_ratio` is present for graph-producing rows and remains diagnostic only.

## Iteration-Specific Check

`nested_top_level_wrong` remains a fast-profile CI fixture. Its exact-path rows
remain visible:

- `chunk_window_smooth_or_crush`: `pass/pass`, `candidate_count=2`, `bubble_count=2`, `flubble_count=2`, `path_replay_compression_ratio=3.459459`.
- `chunk_window_sweepga_seqwish`: `pass/pass`, `candidate_count=2`, `bubble_count=2`, `flubble_count=2`, `path_replay_compression_ratio=3.459459`.
- Existing one-parent compact placeholders: `pass/fail`, `candidate_count=1`, `bubble_count=1`, `flubble_count=1`, `path_replay_compression_ratio=2.064516`.
- Internal-control rows: `pass/fail`, `candidate_count=1`, `bubble_count=0`, `flubble_count=0`, `path_replay_compression_ratio=1.0`.

No hidden graph-quality acceptance gate was added. The runner still writes rows
for topology failures, control failures, skips, and exact-path failures; exact
path corruption remains the only hard-fail status.
