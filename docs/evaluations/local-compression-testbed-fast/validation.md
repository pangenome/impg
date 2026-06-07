# Local Compression Testbed Fast Validation

Date: 2026-06-07

This artifact records validation for the internal impg control runner update.

## Commands

- `git submodule update --init vendor/syng vendor/gfaffix`: passed.
- `source ./env.sh && cargo build`: passed with existing warnings.
- `python3 -m py_compile scripts/local_compression_testbed.py`: passed.
- `source ./env.sh && cargo test -p impg commands::graph::tests::restores_direct_fasta_coordinate_like_path_names_after_pggb_smoothing`: passed.
- `source ./env.sh && cargo test --test test_graph_poa graph_pggb_cli_preserves_direct_fasta_coordinate_headers_after_smoothing`: passed.
- `source ./env.sh && cargo test --test test_local_compression_testbed`: passed, 3 tests.
- `python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast`: passed, wrote 143 rows.
- `source ./env.sh && cargo test`: passed.
- `source ./env.sh && cargo install --path .`: passed and replaced the global `impg` and `gfaffix` executables with this worktree build.

## Scoreboard Summary

- Rows: 143.
- Control rows: 39 total.
- Included control graph rows: 30 command passes.
- Included control graph rows: 30 exact path-preservation passes.
- Fast-profile local-tier control skips: 9 `profile_excludes_local_fixture` rows.
- Control topology statuses: 30 `fail` and 9 `not_run`; these remain visible topology assertions in the scoreboard and are not path-preservation failures.

## Control Interpretation

- Control rows are generated through `impg graph --gfa-engine pggb`, not standalone `pggb` or `smoothxg`.
- `pggb_control` maps to the internal pggb engine with one smoothing target.
- `pggb_plus_smoothxg_control` maps to the internal pggb engine with two smoothing targets.
- `smoothxg_control` maps to the internal pggb engine's smoothxg-style smoothing stage over an explicit empty PAF, because this CLI does not expose a separate smooth-only graph command.
- The previous PanSN/interval `:1-N` to `:1-(N+1)` path-name rewrite was localized to the internal pggb smoothing/lacing path for direct FASTA headers. The graph builder now restores exact direct FASTA path names after smoothing when the output path spelling matches the original input record, so included control outputs round-trip the input FASTA path names and spellings exactly.
