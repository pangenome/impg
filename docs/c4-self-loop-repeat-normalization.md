# C4 Self-Loop Repeat-Unit Normalization

Task: `c4-self-loop-repeat-normalization`

Date: 2026-06-05

## Summary

The C4 direct self-loop edges are already present in the SweepGA seed graph and are not introduced by the 1 kb POA polish. POA slightly reduces them: 152 direct self-loop edges in the seed graph to 144 in the polished graph. Adjacent same-node path repeats likewise drop from 617,670 to 604,675.

These loops are overwhelmingly short repeat-unit representations. In the polished graph, 139/144 self-loop nodes are 1 bp, 2 are 2 bp, 2 are 3 bp, and 1 is 111 bp. All observed adjacent same-node repeat steps are same-signed repeats, so they are path-local run spellings such as `A+,A+,A+` rather than strand-flip artifacts.

I added:

- `impg graph-report` direct self-loop diagnostics: length buckets, top self-loop nodes, node sequence previews, path visits, distinct paths, adjacent same-node/same-step repeats, repeat-run counts, and max run length.
- `impg normalize-self-loops`: an explicit path-preserving normalizer that collapses maximal same-signed self-loop runs into longer run nodes, relaces zero-overlap links from the transformed paths, and rejects output unless `resolution::path_sequences` is exactly unchanged.
- Regression tests for 1 bp and 2 bp self-loop repeat nodes, reverse W-line walks, max unit length scoping, numeric segment ID preservation for C4/gfasort compatibility, and CLI parsing.

Recommendation: keep this as an explicit graph-normalization/reporting tool for now, not as an automatic crush/POA/gfaffix default. It is exact-path preserving and removes a real artifact, but it changes the representation from shared repeat-unit loops to run-length nodes and can worsen order-sensitive graph-report layout metrics unless followed by layout/sort handling. It is suitable as an optional final normalization when downstream consumers or visual reports need loop-free local repeat runs.

## Inputs

- Seed graph:
  `/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa`
- POA-polished graph:
  `/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa`

## Outputs

Output directory:

`/home/erikg/impg/data/c4_self_loop_normalization_20260605T130000Z/`

Key files:

- `seed.graph-report.tsv`
- `seed.graph-report.md`
- `poa1kb.graph-report.tsv`
- `poa1kb.graph-report.md`
- `poa1kb.Ygs.graph-report.tsv`
- `full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.gfa`
- `normalized.graph-report.tsv`
- `normalized.Ygs.graph-report.tsv`
- `normalize-self-loops.stats.json`
- `compare_gfa_paths.seed-vs-normalized.stdout.log`
- `full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.gfa`
- `full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png`

Uploaded PNG:

`https://hypervolu.me/~erik/impg/c4-self-loop-normalized-poa1kb.Ygs.mean-depth.png`

## Commands

Build/test commands used the repo's documented local htslib/jemalloc environment:

```bash
git submodule update --init --recursive vendor/syng vendor/gfaffix

env \
  C_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
  CPLUS_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
  CPATH=/home/erikg/wfmash/build/vendored_htslib/include \
  LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/usr/lib/x86_64-linux-gnu \
  LD_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/home/erikg/.cargo/lib:${LD_LIBRARY_PATH:-} \
  PKG_CONFIG_PATH=/home/erikg/wfmash/build/vendored_htslib/lib/pkgconfig \
  CFLAGS=-I/home/erikg/wfmash/build/vendored_htslib/include \
  CXXFLAGS=-I/home/erikg/wfmash/build/vendored_htslib/include \
  LDFLAGS='-L/home/erikg/wfmash/build/vendored_htslib/lib -Wl,-rpath,/home/erikg/wfmash/build/vendored_htslib/lib -L/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib -Wl,-rpath,/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib -L/usr/lib/x86_64-linux-gnu -Wl,-rpath,/usr/lib/x86_64-linux-gnu -lcurl -L/home/erikg/.cargo/lib -Wl,-rpath,/home/erikg/.cargo/lib' \
  CMAKE_PREFIX_PATH=/home/erikg/wfmash/build/vendored_htslib \
  CMAKE_INCLUDE_PATH=/home/erikg/wfmash/build/vendored_htslib/include \
  CMAKE_LIBRARY_PATH=/home/erikg/wfmash/build/vendored_htslib/lib:/gnu_old/store/wg716bmhd47nxnspc8m0lnmhc1an12n5-jemalloc-5.3.0/lib:/usr/lib/x86_64-linux-gnu \
  cargo test gfa_self_loops --lib

env ... cargo test graph_report --lib
env ... cargo test test_normalize_self_loops_cli_parse --bin impg
env ... cargo build --bin impg
```

Diagnostics and normalization:

```bash
out=/home/erikg/impg/data/c4_self_loop_normalization_20260605T130000Z
seed=/home/erikg/impg/data/c4_whole_region_sweepga_seed_20260603T092921Z/graphs/one_many_minmatch1_scaffold0.initial.gfa
poa=/home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.gfa
norm=$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.gfa

target/debug/impg graph-report -g "$seed" -o "$out/seed.graph-report.tsv" --format tsv --top 20
target/debug/impg graph-report -g "$poa" -o "$out/poa1kb.graph-report.tsv" --format tsv --top 20
target/debug/impg graph-report -g "$seed" -o "$out/seed.graph-report.md" --format markdown --top 20
target/debug/impg graph-report -g "$poa" -o "$out/poa1kb.graph-report.md" --format markdown --top 20

target/debug/impg normalize-self-loops \
  -g "$poa" \
  -o "$norm" \
  --report "$out/normalize-self-loops.stats.json"

target/debug/impg graph-report -g "$norm" -o "$out/normalized.graph-report.tsv" --format tsv --top 20
target/debug/impg graph-report \
  -g /home/erikg/impg/data/c4_full_poa1kb_20260605T104518Z/full_c4.one_many_minmatch1_scaffold0.poa1kb.Ygs.gfa \
  -o "$out/poa1kb.Ygs.graph-report.tsv" \
  --format tsv --top 20

/home/erikg/impg/target/release/examples/compare_gfa_paths "$seed" "$norm" \
  > "$out/compare_gfa_paths.seed-vs-normalized.stdout.log" \
  2> "$out/compare_gfa_paths.seed-vs-normalized.stderr.log"

/home/erikg/.cargo/bin/gfasort \
  -i "$norm" \
  -o "$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.gfa" \
  -p Ygs -t 32

/home/erikg/.cargo/bin/gfalook \
  -i "$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.gfa" \
  -o "$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png" \
  -x 2400 -y 1200 -m -t 32 -v 1

target/debug/impg graph-report \
  -g "$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.gfa" \
  -o "$out/normalized.Ygs.graph-report.tsv" \
  --format tsv --top 20

scp "$out/full_c4.one_many_minmatch1_scaffold0.poa1kb.selfloop-normalized.Ygs.mean-depth.png" \
  erik@hypervolu.me:www/impg/c4-self-loop-normalized-poa1kb.Ygs.mean-depth.png
```

## Results

| Graph | Segments | Links | Paths | Path steps | Segment bp | Direct self-loop edges | Adjacent same-node repeat steps | Repeat runs | Max run | Length buckets |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| SweepGA seed | 8,186 | 11,303 | 465 | 4,632,725 | 232,093 | 152 | 617,670 | 203,597 | 29 | 1 bp: 147; 2 bp: 2; 3 bp: 2; 111 bp: 1 |
| POA 1 kb polished | 9,079 | 12,339 | 465 | 4,843,152 | 240,589 | 144 | 604,675 | 198,015 | 29 | 1 bp: 139; 2 bp: 2; 3 bp: 2; 111 bp: 1 |
| POA 1 kb polished, `gfasort -p Ygs` | 9,079 | 12,339 | 465 | 4,843,152 | 240,589 | 144 | 604,675 | 198,015 | 29 | 1 bp: 139; 2 bp: 2; 3 bp: 2; 111 bp: 1 |
| Self-loop normalized POA | 9,416 | 12,994 | 465 | 4,238,477 | 243,782 | 0 | 0 | 0 | 0 | none |
| Self-loop normalized POA, `gfasort -p Ygs` | 9,416 | 12,994 | 465 | 4,238,477 | 243,782 | 0 | 0 | 0 | 0 | none |

Normalization stats from `normalize-self-loops.stats.json`:

```json
{
  "input_direct_self_loop_edges": 144,
  "output_direct_self_loop_edges": 0,
  "input_adjacent_same_node_path_steps": 604675,
  "output_adjacent_same_node_path_steps": 0,
  "input_adjacent_same_step_path_steps": 604675,
  "output_adjacent_same_step_path_steps": 0,
  "normalized_nodes": 144,
  "collapsed_runs": 198015,
  "created_segments": 337,
  "added_links": 799,
  "removed_self_loop_links": 144,
  "paths_changed": 465,
  "path_steps_before": 4843152,
  "path_steps_after": 4238477
}
```

Required path validation:

```text
expected_paths	465
observed_paths	465
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

## Top Repeat Nodes

The top self-loop nodes are the same dominant short-unit repeats before and after POA. Examples from the polished graph:

| Node | Unit | Visits | Paths | Adjacent same-step repeats | Runs | Max run |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| `1180` | `G` | 161,984 | 465 | 84,442 | 48,924 | 6 |
| `1171` | `C` | 156,947 | 465 | 59,902 | 34,083 | 7 |
| `5655` | `A` | 42,503 | 465 | 25,806 | 13,442 | 5 |
| `8006` | `T` | 21,040 | 459 | 19,207 | 1,833 | 29 |
| `4157` | `A` | 27,707 | 465 | 14,506 | 6,697 | 5 |

The seed graph has the same top nodes and very similar counts. This argues against POA being the source. POA/seqwish replacement did not amplify these loops; it reduced the direct self-loop count by 8 and the adjacent repeat-step count by 12,995.

The persisted C4 POA artifact directory does not include separate unchop/gfaffix intermediate GFAs, so those stages could not be isolated directly from saved files. The saved sorted POA graph does show that `gfasort -p Ygs` is not an amplifier: it has exactly the same 144 direct self-loop edges and 604,675 adjacent same-step repeats as the unsorted POA graph.

## Interpretation

The artifact is a local repeat-run representation: a short segment stores one repeat unit and a direct `L` self-loop plus repeated path steps encode per-path run length. It is path-correct but graph-topology unfriendly, because high-depth single-base nodes acquire huge visitation counts and direct self-loop edges.

The implemented normalizer converts each maximal same-signed run of a targeted self-loop node into a longer run node. For example:

```text
S	a	A
L	a	+	a	+	0M
P	p	a+,a+,a+	*
```

becomes a run-node spelling:

```text
S	a	A
S	a_selfloop_run3	AAA
P	p	a_selfloop_run3+	*
```

On numeric C4 graphs, new run nodes receive numeric IDs above the existing max so `gfasort`/`gfalook` continue to parse the graph.

The transform relaces all adjacent path steps as `0M` links and then hard-gates on exact path spellings. No graph-quality metric is used to accept or reject output. In the C4 polished graph, the exact seed-vs-normalized comparison reports 0 spelling mismatches.

## Placement Recommendation

Do not put this in `gfaffix`. `gfaffix`/sorting are not the source of the loops, and the existing loops are valid graph spellings.

Do not enable it unconditionally in POA replacement or crush finalization yet. It removes the self-loop artifact, but it changes the representation from shared repeat-unit nodes to run-length nodes and may alter order-sensitive layout/report metrics. For example, the unsorted normalized graph has 0 self-loop repeats but still needs sorting for a sensible visual.

Best current placement:

1. Keep the `graph-report` warning and detailed self-loop diagnostics always available.
2. Keep `impg normalize-self-loops` as an explicit, exact-path-preserving post-process.
3. Consider an opt-in crush final normalization stage only if downstream consumers require loop-free local runs. That stage should keep the same hard gate: exact path preservation only.
