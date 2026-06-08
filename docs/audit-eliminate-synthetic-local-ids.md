# Audit: synthetic local IDs in query/crush sequence naming

Date: 2026-06-02
Task: `audit-and-eliminate`

## Required visible-name contract

Visible FASTA headers, PAF query/target names, GFA path names, replacement traversal names, local induction sequence names, and debug artifacts must remain in the same semantic namespace emitted by `impg query -o fasta`, `impg query -o gfa`, and related GFA-style query outputs.

For local haplotype subsets, the only allowed name change is to adjust the range to the selected interval. Source/PanSN/path identity must otherwise be preserved. Sidecar maps may carry extra provenance, but they are not a replacement for correct visible names and there is no hidden alias namespace for downstream tools to recover from.

The local replacement/crush naming convention after this audit is:

- Preserve the first FASTA token of the source path name.
- If the source name already has a trailing `:start-end` interval, compute the selected interval in that coordinate system and emit `source:(start+local_start)-(start+local_end)`.
- If the source name has no trailing interval, emit `source:local_start-local_end`.
- If the same source interval is emitted more than once, append `|duplicate-source-interval-copyN`.
- Do not emit opaque visible names such as `seq_0`, `path_17`, `candidate_3`, `traversal_2`, `local`, `replacement`, or `__impg_bubble_path...` in query/crush/local-induction FASTA, PAF, GFA, or debug artifact paths.

## Findings

| Location | Classification | Result |
| --- | --- | --- |
| `src/resolution.rs:643` `PathRange` | Visible naming input was incomplete. Replacement candidates retained source path text but not selected source coordinates. | Fixed by adding `source_start` and `source_end`, populated from path-position tables for root intervals, discovered candidates, and chained candidate groups. |
| `src/resolution.rs:6859` `candidate_named_sequences` and `src/resolution.rs:6885` `candidate_sequence_name` | Visible synthetic IDs were written into replacement FASTA headers, propagated into PAF names, GFA path names, and `combined.fa` debug artifacts. Old names included `__impg_bubble_path{path_idx}_{ordinal}` and could fall back to path-index-only aliases. | Fixed. Headers now preserve source/PanSN/path identity and adjust only the selected interval range. Missing source names are treated as errors instead of silently inventing aliases. |
| `src/resolution.rs:6931` `candidate_named_sequences_longest_first` | Visible names were sorted and forwarded to induction/lacing tools. | Fixed through the shared semantic-name builder; PAF/GFA consumers receive the same headers as FASTA. |
| `src/resolution.rs:7175` `build_chain_povu_smooth_poasta_replacement` and `src/resolution.rs:8223` `polish_replacement_gfa` | Full-range path names were stripped before smoothing/polishing, collapsing names such as `ref:0-5` to `ref`. This broke the stricter contract because ranges must remain visible. | Fixed by removing the range-stripping step. Polishing now preserves semantic path names. |
| `src/syng_graph.rs:955` `build_gfa_from_paf_and_sequences` and `src/syng_graph.rs:992` debug `combined.fa` | Pass-through visible artifact writer. It did not invent names, but it exposed upstream synthetic replacement names in debug output. | Fixed indirectly by correcting upstream replacement names and active local-induction fixtures. C4 debug `combined.fa` now shows semantic names. |
| `tests/test_data/crush/c4_fragments/*.fa`, `*.paf`, `*.gfa` and `tests/test_data/crush/top_flubble_seqwish_minrun.*` | Active local-induction fixtures passed synthetic names to the seqwish/BIWFA/GFA test path. | Fixed by rewriting fixture-visible names to `C4FIXTURE#0#<fixture>:<start>-<end>`. These are semantic fixture-source names because archived original sample coordinates are not available for these reduced diagnostic fragments. |
| `src/lib.rs:735`, `src/lib.rs:1092`, `src/lib.rs:1110` local syng GFA construction | Visible but already semantically compatible. | No code change. These paths use `SequenceMetadata::path_name()` / selected `seq:start-end` names and therefore stay in the query-output namespace. |
| `src/commands/syng2gfa.rs:2915`, `src/commands/syng2gfa.rs:3130`, `src/commands/syng2gfa.rs:3291` | Visible but already semantically compatible. Split GFA paths use semantic source path names and deterministic `|partN` suffixes for split segments. | No code change. |
| `src/commands/graph.rs:1171` and `src/commands/align.rs` temporary preparation paths | Internal-only temporary file paths; sequence names are parsed/preserved from input first FASTA tokens. | No code change. |
| `src/graph_report.rs` `_povu_flubble_path_{rank}_lv...` report paths | Visible graph-report overlay/debug naming, outside the query/crush/local-induction replacement path covered by this task. | Not changed in this audit. It should be handled separately if graph-report overlays are later brought under the same query-output-name contract. |
| `src/syng.rs` test helpers using `seq_{}` / `region_seq_{}` | Internal test input names only. | No code change. |
| `candidate_index`, `traversal_rank`, and similar fields in genotype/reporting code | Internal/report indices, not FASTA/PAF/GFA sequence or path names. | No code change. |
| Historical diagnostic docs mentioning `__impg_bubble_path...` | Stale documentation of old behavior, not executable output. | Superseded by this report. Production and active fixture outputs were changed. |

## Tests added or updated

- `src/resolution.rs:10344` `duplicate_source_intervals_get_meaningful_unique_suffixes`
- `src/resolution.rs:10394` `candidate_interval_name_matches_query_fasta_contract`
- `src/resolution.rs:10503` `replacement_paf_sequence_names_match_fasta_headers`
- `src/resolution.rs:10545` `local_replacement_gfa_path_names_preserve_semantic_source_names`
- Updated existing replacement header/path-name tests so they reject synthetic local IDs and preserve first-token FASTA semantics.

## Validation

The following commands passed with the htslib/jemalloc validation environment used in this worktree:

- `cargo check --lib`
- `cargo test --lib candidate_ -- --nocapture`
- `cargo test --lib duplicate_source_intervals_get_meaningful_unique_suffixes -- --nocapture`
- `cargo test --lib replacement_paf_sequence_names_match_fasta_headers -- --nocapture`
- `cargo test --lib local_replacement_gfa_path_names_preserve_semantic_source_names -- --nocapture`
- `cargo test --lib chain_povu_resolves_nested_parent_as_one_block -- --nocapture`
- `cargo test --test test_crush_integration c4_fragment_seqwish_regressions_induce_shared_graphs -- --nocapture`
- `cargo test --test test_crush_integration c4_fragment_lacing_uses_pairwise_induced_graphs -- --nocapture`
- `cargo test --bin impg`
- `cargo build`
- `cargo test`
- `cargo install --path .`

C4 diagnostic/local-induction debug artifacts were regenerated under:

`/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/`

Representative semantic `combined.fa` headers:

- `/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/graph_build_0000/combined.fa`: `>C4FIXTURE#0#easy_shared_flank:14-240`
- `/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/graph_build_0001/combined.fa`: `>C4FIXTURE#0#bounded_multi_bubble:2-200`
- `/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/graph_build_0002/combined.fa`: `>C4FIXTURE#0#top_flubble_seqwish_minrun:0-399`
- `/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/graph_build_0003/combined.fa`: `>C4FIXTURE#0#short_floor:0-250`
- `/home/erikg/impg/data/audit-and-eliminate-c4-debug/crush/graph_build_0004/combined.fa`: `>C4FIXTURE#0#duplicated_repeat:3-225`

Synthetic-name scan:

`rg -n "__impg|seq_|path_|candidate_" /home/erikg/impg/data/audit-and-eliminate-c4-debug`

No matches were found.
