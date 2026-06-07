# Local Compression Testbed Fast Validation

Task: `follow-up-produce`
Date: 2026-06-07

## Native Build Environment

The Rust build in this worktree needs the local htslib and jemalloc paths used
by earlier C4 validation tasks:

```bash
export C_INCLUDE_PATH=/home/erikg/htslib-local/include
export CPLUS_INCLUDE_PATH=/home/erikg/htslib-local/include
export CPATH=/home/erikg/htslib-local/include
export LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/usr/lib/x86_64-linux-gnu
export LD_LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/home/erikg/.cargo/lib:${LD_LIBRARY_PATH:-}
export PKG_CONFIG_PATH=/home/erikg/htslib-local/lib/pkgconfig:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib/pkgconfig
export CFLAGS=-I/home/erikg/htslib-local/include
export CXXFLAGS=-I/home/erikg/htslib-local/include
export LDFLAGS='-L/home/erikg/htslib-local/lib -Wl,-rpath,/home/erikg/htslib-local/lib -L/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib -Wl,-rpath,/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -lcurl -L/home/erikg/.cargo/lib -Wl,-rpath,/home/erikg/.cargo/lib'
export CMAKE_PREFIX_PATH=/home/erikg/htslib-local
export CMAKE_INCLUDE_PATH=/home/erikg/htslib-local/include
export CMAKE_LIBRARY_PATH=/home/erikg/htslib-local/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/usr/lib/x86_64-linux-gnu
```

## Commands Run

```bash
git submodule update --init --recursive vendor/syng vendor/gfaffix
python3 scripts/local_compression_testbed.py write-fixtures --root tests/test_data/local_compression
python3 scripts/local_compression_testbed.py validate-fixtures --manifest tests/test_data/local_compression/manifest.json
python3 scripts/local_compression_testbed.py run --profile fast --manifest tests/test_data/local_compression/manifest.json --out-dir docs/evaluations/local-compression-testbed-fast
cargo build
cargo test --test test_local_compression_testbed
cargo test
cargo install --path .
```

## Results

- Fixture validation passed: 13 fixtures validated from `tests/test_data/local_compression/manifest.json`.
- Fast profile passed: wrote 143 rows to `docs/evaluations/local-compression-testbed-fast/scoreboard.tsv` and `docs/evaluations/local-compression-testbed-fast/scoreboard.json`.
- Scoreboard coverage: 13 fixtures, 11 methods, 80 produced graph rows, 63 skipped rows.
- Optional control skip reasons are explicit for every optional-control row: `profile_excludes_optional` for CI fixtures and `profile_excludes_local_fixture` for local-only fixtures in the fast profile.
- `cargo build` passed with existing warnings.
- `cargo test --test test_local_compression_testbed` passed: 2 tests passed.
- `cargo test` passed: full unit, integration, and doc-test suite passed; known C4 tests remained ignored.
- `cargo install --path .` passed and replaced `/home/erikg/.cargo/bin/impg` and `/home/erikg/.cargo/bin/gfaffix` from this worktree.
