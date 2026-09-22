# Remove crush replacement quality guards

Task: `remove-crush-replacement`

## Semantics

Crush replacement is now unconditional after graph-validity checks. Once a
bubble, chain, or POVU block is selected, the configured replacement builder
runs and the resulting graph is spliced into the working graph when it passes
the correctness invariants:

- the replacement GFA parses;
- replacement path names/order are preserved for the selected traversals;
- each replacement traversal spells the same sequence as the input region;
- the rewritten whole graph preserves the input path spellings.

Compression ratio, segment count, segment-bp, trivial-stringy count, and visual
quality score are diagnostics only. They are logged so runs are auditable, but
they do not reject an otherwise valid replacement and do not keep the previous
graph.

The legacy compression-retry CLI keys are still accepted for old command
strings. `retry-min-compression-ratio` now only annotates diagnostic logs, and
the deprecated retry counters remain zero.

## Code changes

- Removed the compression-ratio retry path that rebuilt a candidate with an
  alternate aligner and chose between replacements by output segment-bp.
- Changed `DEFAULT_RETRY_MIN_COMPRESSION_RATIO` to `0.0`; nonzero configured
  values are diagnostic only.
- Kept replacement compression logging for accepted plans without using those
  metrics as an acceptance gate.
- Changed chain-POVU smoothxg -> POASTA behavior so a valid smoothed
  replacement is used unconditionally. Direct POASTA is only used as a
  validity path when a traversal is empty or the smoothed path fails to build,
  parse, or preserve traversal spellings.
- Replaced the chain-POVU metric-decision counters with validity-path counters:
  `smooth_used`, `direct_on_empty`, and `direct_on_smooth_failure`.

## Rust validation

The checkout needs the native dependency paths used by the existing C4 work:

```bash
CARGO_BUILD_JOBS=8 \
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo test --all
```

Result: pass.

Additional checks run before the full suite:

```bash
rustfmt --edition 2021 --check src/resolution.rs
cargo test --lib compression_ratio_threshold_is_diagnostic_only -- --nocapture
cargo test --lib chain_povu_resolves_nested_parent_as_one_block -- --nocapture
```

The installed binary was refreshed from this worktree:

```bash
CARGO_BUILD_JOBS=8 \
CFLAGS='-I/home/erikg/wfmash/build/vendored_htslib/include -I/home/erikg/micromamba/include' \
LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -L/home/erikg/micromamba/lib' \
CMAKE_PREFIX_PATH='/home/erikg/wfmash/build/vendored_htslib;/home/erikg/micromamba' \
cargo install --path .
```

Result: pass; `impg` and `gfaffix` were installed from this worktree.

## Real C4 validation

Output directory:

```text
/home/erikg/impg/data/c4_remove_crush_replacement_chain_povu_20260527T141634Z
```

Command:

```bash
out=/home/erikg/impg/data/c4_remove_crush_replacement_chain_povu_20260527T141634Z
LD_LIBRARY_PATH='/home/erikg/wfmash/build/vendored_htslib/lib:/home/erikg/micromamba/lib:'"${LD_LIBRARY_PATH:-}" \
/usr/bin/time -v -o "$out/time.txt" \
  impg query \
    -t 32 \
    -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
    -r 'GRCh38#0#chr6:31891045-32123783' \
    --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
    -d 50k \
    -o 'gfa:syng:mask,min-run=3:crush,method=chain-povu,max-rounds=3,aligner=fastga,min-traversal-len=5k,seqwish-k=311,max-pair-alignments=0,max-paf-bytes=0,no-filter=true,polish-rounds=until-done,polish-method=poasta,polish-max-traversal-len=10k,polish-max-median-traversal-len=1k:nosort' \
    -O "$out/run.nosort" \
    -v 1 \
    > "$out/run.stdout" \
    2> "$out/run.stderr"
```

Result: exit status 0. `/usr/bin/time` wall clock was `7:00.97`; the query log
reported `Syng query complete ... in 3m34.861s` after index loading.

## C4 metrics

Metrics from `run.nosort.gfa`:

| Metric | Value |
|---|---:|
| paths | 465 |
| duplicate path-name groups | 0 |
| segments | 193,318 |
| links | 320,654 |
| segment_bp | 70,358,675 |
| trivial_stringy | 99 |
| block count | 1,571 |
| round count | 3 |
| bailed replacements | 0 |
| candidates seen | 1,571 |
| wall time | 7:00.97 |
| max RSS | 55,705,316 KiB |

The final graph preserved 465/465 paths.

Round summaries from stderr:

```text
round 1: accepted 714/714 block replacement(s); smooth_used=714, direct_on_empty=0, direct_on_smooth_failure=0
round 2: accepted 493/493 block replacement(s); smooth_used=493, direct_on_empty=0, direct_on_smooth_failure=0
round 3: accepted 364/364 block replacement(s); smooth_used=364, direct_on_empty=0, direct_on_smooth_failure=0
crush: 1571 resolved, 0 bailed, 1571 candidates seen across 3 rounds
```

The quality metrics worsened substantially while replacements were still kept,
which is the intended experiment semantics for this task:

```text
round 1 quality: segment-bp 389354 -> 46981500
round 2 quality: segment-bp 46981500 -> 61488140
round 3 quality: segment-bp 61488140 -> 70358675
```

A grep over `run.stderr` found no legacy quality-guard keep/reject markers:
no direct-vs-smooth metric selector, no quality rejection, no keep-old result,
and no compression-retry attempt. Result: no matches.

The trivial-stringy heuristic was run with `/tmp/find_stringy_bubbles.py`:

```text
# trivial-stringy bubble candidates: 99
# parsed 193318 segments, 465 paths
# anchors (visited by >= 372/465 paths): 2016
```

## PNG

The final graph was rendered directly from `run.nosort.gfa`:

```bash
gfalook -i "$out/run.nosort.gfa" \
  -o "$out/c4-remove-crush-replacement-chain-povu.png" \
  -x 2200 -y 1200 -m -O -H
```

Uploaded artifact:

```text
www/impg/c4-remove-crush-replacement-chain-povu.png
```

Remote confirmation:

```text
-rw-r--r-- 1 erik erik 63K May 27 14:37 www/impg/c4-remove-crush-replacement-chain-povu.png
```
