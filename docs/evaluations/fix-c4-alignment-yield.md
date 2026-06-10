# Fix C4 Alignment Yield in Localized Polish/Crush

Task: `fix-c4-near`
Date: 2026-06-10

## Diagnosis

The rejected `fix-c4-localized` run did not fix alignment yield. Its expanded
C4 localized-polish run selected 32 dirty chunks and applied 0 replacements.
The exact collapse point in that run was before alignment: every selected
localized chunk was skipped by `localized_candidate_from_dirty_chunk` because
the chunk interval lacked both shared flank anchors on its path. No candidate
sequences, pair generation, PAF records, seqwish induction, or replacement
lacing happened for those chunks.

There was a second code-level hazard in the replacement wrapper: local
replacement defaults still silently requested strict `1:1` SweepGA filtering
through `ResolutionConfig` and the `impg crush` CLI defaults. That is unsafe
for C4-style many-haplotype blocks unless the user explicitly opts into it.
This patch changes default local replacement mapping/scaffold filters to
`many:many` and keeps `replacement_scaffold_mass=0`; explicit CLI or engine
configuration still overrides those settings.

The final C4 validation below shows the useful replacement path does not
collapse at sequence extraction, PanSN/name compatibility, pair generation,
SweepGA/FastGA, wrapper filtering, seqwish induction, or lacing. The hard
safety invariant remains exact path spelling, not graph-quality movement.

## Implementation

- `ResolutionConfig::default()` now uses `replacement_num_mappings=many:many`
  and `replacement_scaffold_filter=many:many`.
- The standalone `impg crush` CLI defaults now match `many:many`.
- `LocalizedPolishConfig::default()` explicitly keeps localized replacements at
  `many:many` and `replacement_scaffold_mass=0`.
- Localized polish selection now prefers chunks whose dirty core interval is
  strictly inside the path when `path_length_bp` is available. Path-edge chunks
  are deferred until after anchorable interior chunks, so the budget is not
  spent entirely on regions that cannot have both flank anchors.
- Replacement evidence now records raw/filtered PAF counts, unique pair counts,
  median/max alignment length, effective min-match/min-map/filter settings,
  seqwish input PAF path, and seqwish graph shape.
- Localized debug output writes
  `localized_polish_alignment_yield.tsv` beside
  `localized_polish_summary.tsv`.
- Direct `impg crush` accepted-replacement logs now include the same
  `alignment_evidence=...` block, not only accepted/failed progress.

Synthetic regressions added:

- `localized_polish_prefers_anchorable_interior_chunks`
- `localized_polish_defaults_do_not_hide_one_to_one_pairwise_filter`
- `replacement_pair_filters_default_to_many_many`
- `test_localized_polish_explicit_replacement_filter_is_honored`

## C4 Validation Run

Generated artifacts are outside git under:

```text
/home/erikg/impg/data/fix_c4_near_20260610T105108Z_crush_manymany_large_final
```

Input graph:

```text
/home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z/graphs_tuned/C4_GRCh38_53kb.gfa
```

Command:

```bash
IMPG_CRUSH_DEBUG_DIR="$out/debug/crush" \
/usr/bin/time -v -o "$out/logs/crush.time.txt" \
  /usr/bin/timeout 600s target/release/impg crush \
    -g /home/erikg/impg/data/diagnose_exact_c4_20260609T091222Z/graphs_tuned/C4_GRCh38_53kb.gfa \
    -o "$out/C4_GRCh38_53kb.crush-sweepga-manymany.large.gfa" \
    --method sweepga \
    --max-iterations 1 \
    --max-span 7k \
    --min-traversal-len 5k \
    --max-traversal-len 7k \
    --max-median-traversal-len 7k \
    --max-total-sequence 2m \
    --max-traversals 10k \
    --polish-rounds 0 \
    -t 16 -v 1
```

Runtime summary:

| metric | value |
|---|---:|
| exit status | 0 |
| wall time | 0:07.35 |
| max RSS | 663,564 KB |
| selected C4 candidates | 2 |
| accepted replacements | 2 |
| bailed candidates | 0 |
| input graph shape | 3,851 segments / 5,220 links / 466 paths |
| output graph shape | 3,990 segments / 5,455 links / 466 paths |

Exact path validation:

```text
expected_paths       466
observed_paths       466
missing_paths        0
extra_paths          0
spelling_mismatches  0
```

## Alignment-Yield Table

`aligned_pair_fraction` is computed as unique unordered non-self aligned
sequence pairs divided by `n * (n - 1) / 2`. PAF records can include multiple
records per unordered pair, so this fraction is more conservative than
record-count ratios. The direct runtime log also records ordered raw/filtered
unique-pair counts in the `alignment_evidence=...` field.

| block / chunk id | step range | path count | candidate sequence count | candidate bp | expected unordered pairs | raw PAF records | filtered PAF records | aligned-pair fraction | median / max aln len | settings | seqwish input PAF | seqwish shape | replacement result |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|---|---|---|
| `Sweepga_crush01_candidate_span6460bp_med6461bp_cov47of466_sites47_steps786-968` | `786..968` | 47 | 47 | 239,974 | 1,081 | 1,469 | 1,422 | 0.657724 | 6,461 / 6,462 | `min_match=1`, `min_map=1`, `num_mappings=many:many`, `scaffold_filter=many:many`, `scaffold_mass=0`, `no_filter=false` | `/home/erikg/impg/data/fix_c4_near_20260610T105108Z_crush_manymany_large_final/debug/crush/graph_build_0000/filtered.paf` | `S=62`, `bp=6573`, `L=81`, `P=47` | applied; exact path spelling passed |
| `Sweepga_crush02_candidate_span23bp_med23bp_cov49of466_sites49_steps2757-2760` | `2757..2760` | 49 | 49 | 80,228 | 1,176 | 9 | 6 | 0.002551 | 26,390 / 26,390 | `min_match=1`, `min_map=1`, `num_mappings=many:many`, `scaffold_filter=many:many`, `scaffold_mass=0`, `no_filter=false` | `/home/erikg/impg/data/fix_c4_near_20260610T105108Z_crush_manymany_large_final/debug/crush/graph_build_0001/filtered.paf` | `S=47`, `bp=27448`, `L=0`, `P=49` | applied; exact path spelling passed |

Key runtime evidence lines from `logs/crush.stderr.log`:

```text
crush round 1: 2 resolved ...; total resolved=2
alignment_evidence=raw_paf_records=1469,filtered_paf_records=1422,...,num_mappings=many:many,scaffold_filter=many:many,scaffold_mass=0,...,seqwish_segments=62,seqwish_links=81,seqwish_paths=47
alignment_evidence=raw_paf_records=9,filtered_paf_records=6,...,num_mappings=many:many,scaffold_filter=many:many,scaffold_mass=0,...,seqwish_segments=47,seqwish_links=0,seqwish_paths=49
```

## Stage Assessment

| stage | C4 result |
|---|---|
| Sequence extraction | passed; 47 and 49 candidate sequences with 239,974 bp and 80,228 bp |
| Naming / length compatibility | passed; `paf_name_mismatches=0`, `paf_length_mismatches=0` in both debug summaries |
| Pair generation / FastGA invocation | passed for selected useful block; 1,469 raw PAF records for the 47-path block |
| Wrapper filtering | passed; 1,469 raw -> 1,422 filtered records with `many:many` filters; no hidden `1:1` default |
| Min exact run / min-match | not a collapse point; effective min-match/min-map are 1 and exact-run scan is skipped |
| Scaffold filtering | not a collapse point; `scaffold_filter=many:many`, `scaffold_mass=0` |
| Seqwish induction | passed; seqwish returned replacement graphs with paths matching candidate count |
| Replacement lacing | passed; 2 replacements applied, `compare_gfa_paths` reports 0 spelling mismatches |
| Graph quality metrics | diagnostic only; replacement was not gated by quality score |

## Artifact Policy

No generated C4 GFA, PNG, PAF, debug directory, or log artifact is committed.
The C4 PNG was not rendered or uploaded because the task invariant is alignment
yield plus exact path preservation; the final graph-quality metrics are
diagnostic and intentionally not used as an acceptance gate.

The committed artifact is this report only.
