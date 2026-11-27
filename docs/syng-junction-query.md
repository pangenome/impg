# `impg junction-query`: junction words versus assembled haplotype paths

`junction-query` is a read-only, batched **panel-context** query, not a variant caller. It distinguishes (a) an exact word somewhere in the selected source sequences, (b) a contiguous word containing both extended arms and the intervening sequence on one named path, and (c) whether the flanks of that occurrence identify unique positions within the reporting and target assemblies. It does not estimate a population frequency or a sample genotype.

```sh
impg junction-query --index /path/to/panel --assemblies /path/to/source.agc \
    --models models.jsonl --out private-results --examples 3 -t 8
```

Supply an immutable `.1gbwt`/`.1khash` prefix with its `.names`, `.meta`, `.spos` and `.pstep`. PanSN `sample#haplotype#contig` names are required. Select complete named assembly sets with `--sets 'REF#0,donor#1'` (default: all in the index). A source archive can include additional sets. Every selected group must include all its named paths, and the reporting reference must be among them. This command does not append, repair or re-stamp an index. Reading old metadata here does **not** enable `genome-infer` on a legacy panel.

## JSONL input

```
{"id":"example1","end1":["REF#0#chr1",5000,"L"],"end2":["REF#0#chr2",6500,"R"],"insert":"AC","coupled":[[5000,6500]]}
```

Cuts are 0-based interbase. `L` retains sequence below the cut; `R` retains sequence at/above the cut. Orient both arms in the declared end1 → end2 molecule direction; for end1 `R` or end2 `L`, the stored plus-reference arm is reverse-complemented. `insert` is in molecule order; `""` means no insertion; `null` is unknown (flanks assessable, junction words not assessable). Supply a unique opaque ID. `coupled` is an optional list of complete equivalent cut pairs including the representative. For a known empty insertion, the source bases transferred between the two arms must agree exactly; unrelated shifts are rejected. Count one model regardless of the number of placements. Do not infer an inserted sequence from its length.

## What is measured

The native command reads source sequence through the AGC library, checks names and lengths against SYNG, and compares signed `walk_path_range` step spellings to source bases at the reporting arms, a midpoint window per path, and displayed target witnesses. A missing legacy start, seedless interval or ambiguous base does not become a negative seed vote. Hard length/step inconsistencies fail the run. Sampled graph checks are not a claim that every graph step was checked.

An overlapping multi-pattern search visits every selected source base on both strands (one forward scan with both word orientations), one contig at a time. It includes both individual 40-base flanks, 22+insert+22 and 33+insert+33 words, the joined 40+insert+40 word, and analogous 96-base individual and joined contexts. Reverse-complement-palindromic words count once per start. It never concatenates contigs, crosses ambiguous bases, or requires a SYNG seed in a short query. A model at a contig edge can have some assessable words and others unassessable. This is exact, not an alignment tolerant of mismatches in a divergent haplotype.

For `verified_same_path_join`, at least one **entire** joined 40- or 96-base-arm word occurs on one target path in order, with both flanks at that width individually unique in both the reporting set and that target set. Multiplicity is per assembly, not across people. `ambiguous_locus` means a joined 40-base-arm word exists but neither width gives a uniquely anchored witness. A short exact word alone does not meet these conditions. `no_verified_join` means no joined 40-base-arm word on these assembled paths; not absence of all homologous structures. Unknown insert/insufficient flanks can produce `unassessable`.

## Output and recovery

- `summary.jsonl`: one row per model in input order; `literal` contains per-word exact status, reason if not assessed, total occurrences, path and set counts, and nonzero counts per set. `placement`, `unique_anchor_sets` and `ambiguous_only_sets` describe the longer-arm test. Display caps never affect these counts.
- `hits.jsonl`: optional-size sample of same-path placements. Target coordinates are **0-based half-open**. `arm1` and `arm2` retain model order even on reverse-strand hits. For a palindrome the representative orientation is shown and labelled `+/-`. Coupled target cuts outside the verified extent are not projected. `arm_graph_evidence` distinguishes checked stored steps from shared signed steps at the same arm-relative offset; zero or null is not a failed literal match.
- `config.json` binds the model file, AGC, six sidecars, program binary (SHA-256), selected sets, seed policy and output options. `units/` contains independently sealed per-set counts; `complete.json` binds them and both result files. Trust only a complete, hash-verified run. Resume by repeating the identical invocation; never combine different configurations. Use a new output directory if any bound input changes. The output root must be 0700. No query sequences are emitted, but coordinates/counts can still be private.
- Scan memory scales with the AGC decompressor slots, active contigs and query automaton, not the product of nodes and haplotypes. `-t` bounds Rayon workers and `IMPG_AGC_DECOMPRESSOR_SLOTS` can bound the number of simultaneous archive readers further.

`python3 tests/junction_query_cli.py /path/to/impg /path/to/agc /new/fixture-dir` exercises a synthetic five-set archive and compares every word/strand/set count against a naive overlapping search, including a held-out positive and corrupted checkpoint rejection.
