# Syng integration — next steps

**PR:** pangenome/impg#162 (feature/syng-integration)
**Status:** Core syng query pipeline is stable (commit `d48afc8`). Three
investigations below are independent of that pipeline's correctness and
chase three distinct questions.

This document is the roadmap for the remaining work. Each section is
self-contained — you can run them in parallel. Recommended order at the
bottom.

---

## 0. Where we are

**Done (committed on feature/syng-integration):**

- `query_region_with_anchors` in `src/syng.rs` — both GBWT orientations
  queried, anchors tagged by `q_orient XOR t_orient`.
- Absolute-coordinate `walk_path` via `GbwtPathStart::first_syncmer_pos`
  (fresh index rebuild required for correctness — `.syng.names` is 7-col).
- `distance_merge_anchored` — bedtools-style `-d`, plus a co-linearity
  signature check (`target_pos − query_pos` for `+`, `target_pos +
  query_pos` for `-`) that prevents merging paralogs with different
  deltas into one super-interval.
- `dedupe_strand_overlaps` — per-path, overlapping `+`/`-` intervals
  resolved by majority anchor count; non-overlapping both kept.
- `refine_boundaries` — BiWFA EndsFree edge realignment anchored at the
  innermost shared syncmer; linear projection fallback.
- CLI: syng-index `-a <prefix>` bed/fasta/gbwt/gfa all route through the
  pipeline. `-d` honoured. `--syng-raw` bypass for debug.

**Observed on yeast235 `S288C#0#chrIV:408000-410000 -d 10000`:**

| Output format | Wall time | Paths | Graph length | Correct? |
| ------------- | --------- | ----- | ------------ | -------- |
| `-o bed`      | 5s    | — | — | Widths reflect real indels, strands clean |
| `-o fasta`    | 5s    | 269 seqs, 460kb | — | Max interval 2146bp (no ragged ones) |
| `-o gfa --gfa-engine seqwish` | 12s | **269** | **23kb** | One path per input, ~20× compression |
| `-o gfa --gfa-engine pggb`    | 2m25s | **14750** | **3.5Mb** | **Pathological inflation** |

The syng query pipeline produces the right input to the GFA engines.
The `pggb` engine's block-POA step is the source of the inflation.

---

## 1. Reference-alignment validation (ground truth)

### Question
Does the syng query produce interval boundaries and strand
annotations that agree with a traditional PAF-based query for
closely-related genomes?

### Why it matters
We've verified the pipeline is *internally consistent* (co-linearity
signatures, strand dedupe, edge realignment). We have not yet compared
it to an independently-derived ground truth. Until we do, any
disagreement between syng-graph builds and PAF-graph builds could be
blamed on either side.

### Inputs
- `yeast235.agc` — the AGC archive.
- `yeast235` — rebuilt 7-column syng index at `~/scrapy/yeast235*`.

### Steps

1. **Extract reference sequence.** Pull S288C#0 only out of the AGC
   (all chromosomes or just chrIV). If the full genome is too much,
   use `AAGCSTART:END` slicing on chrIV to cover `407000-411000` with
   margin. Save as `S288C_chrIV.fa`.

2. **Build sparse all-vs-reference PAF.** Use sweepga with the `-p`
   or equivalent flag to align every haplotype in yeast235 against
   the reference. Example (adjust flags):

   ```bash
   sweepga yeast235.agc \
       --reference S288C_chrIV.fa \
       --threads 8 \
       -o yeast235_vs_S288C.paf
   ```

   Target output: one PAF file with ~235 haplotypes × ~16 chromosomes
   = a few thousand alignments, filtered to the main collinear blocks.

3. **Build an impg index from the PAF.**

   ```bash
   impg index -a yeast235_vs_S288C.paf -s S288C_chrIV.fa
   ```

4. **Query the same region on the PAF-based index.**

   ```bash
   impg query -a yeast235_vs_S288C.paf \
       --sequence-files yeast235.agc \
       -r 'S288C#0#chrIV:408000-410000' \
       --transitive -d 10000 \
       -o bed > paf_query.bed
   ```

5. **Run the syng query on the same region.**

   ```bash
   impg query -a yeast235 \
       --sequence-files yeast235.agc \
       -r 'S288C#0#chrIV:408000-410000' \
       -d 10000 -o bed > syng_query.bed
   ```

6. **Diff the two outputs.** A simple script: for each `(sample, hap,
   contig)` key, report:
   - Present in both? (set diff on keys)
   - Start delta (syng_start - paf_start): distribution; median, max.
   - End delta: same.
   - Strand agreement: rows where one reports `+` and the other `-`.
   - Width ratio: syng_width / paf_width.

   Put this in a small shell/awk script at
   `tests/validation/compare_syng_vs_paf.sh` so it's re-runnable.

### Acceptance criteria

- **Coverage:** every (sample, hap, contig) in the PAF output should
  appear in the syng output (or the omission should be diagnosable —
  e.g. the region on that haplotype lacked any shared syncmer).
- **Boundaries:** start/end deltas should be ≤20bp for ≥90% of
  haplotypes. Outliers should cluster in specific biological contexts
  (e.g. subtelomeric, transposon-rich) rather than being uniformly
  distributed.
- **Strand:** 100% agreement on strand annotation (or documented
  exceptions for known RC regions).

### Deliverables

- `notes/SYNG_VS_PAF_VALIDATION.md` — writeup of the comparison script,
  the generated PAF index, and observed agreement statistics.
- If agreement is poor, the investigation feeds back into the syng
  pipeline (likely the edge realignment or the co-linearity tolerance).

---

## 2. `smoothxg` inflation diagnosis

### Question
Why does impg's pggb engine turn a clean 23kb seqwish graph (269
paths, one per haplotype) into a 3.5Mb smoothed graph (14750 paths)?
Is it a parameter issue, a lacing bug, or a POA misbehaviour?

### Why it matters
The pggb pipeline is supposed to produce a polished version of the
seqwish output. ~150× length expansion and ~55× path fragmentation
imply something fundamental is wrong. Either our embedded smoothxg
is misconfigured or it has a bug the reference implementation doesn't.

### Inputs
- The clean seqwish GFA: `/tmp/yeast_sw.gfa` (23kb, 269 paths). Re-create
  with:

  ```bash
  impg query -a yeast235 --sequence-files yeast235.agc \
      -r 'S288C#0#chrIV:408000-410000' \
      -o gfa --gfa-engine seqwish -O /tmp/yeast_sw -d 10000
  ```

### Steps

1. **Capture smoothxg's intermediate outputs.** Add debug-dir support
   to `src/smooth.rs::smooth_gfa_pass` (or enable existing logging) to
   dump:
   - Post-chop GFA (should have 59437 nodes at `max_node_length=100`).
   - Block decomposition summary (block start/end positions, member paths).
   - Each per-block POA sub-GFA.
   - Pre-lace state (all block GFAs concatenated).
   - Post-lace GFA.

2. **Compare block-by-block.** Each block's POA should produce a mini-
   graph whose total sequence length is comparable to the average
   sequence length within the block. If any single block produces
   disproportionate length (e.g. 100kb output for a block covering
   2kb of input), that's the bug.

3. **Compare against reference smoothxg.** Run the *external* `smoothxg`
   binary on the same `/tmp/yeast_sw.gfa` with matching parameters
   (`-w 700`, `-l 1100`, etc.) and compare output.
   - Same output shape → our embedded smoothxg is correct; the issue is
     parameters / input assumptions.
   - Different output shape → impg has a bug in its smoothxg port.

4. **Inspect path names.** The current pggb output has paths like
   `ANL#0#chrIV:938209-938512` where 938209 is nowhere in the input.
   These names are being GENERATED by smoothxg — find where in the
   code and verify the logic. Likely candidate: the block's path
   range projection is mis-computing absolute coordinates.

5. **Parameter sweep.** If external smoothxg behaves the same, try:
   - `target_poa_length` 1500, 3000 (default 700)
   - `max_node_length` 500 (default 100)
   - `poa_padding_fraction` 0.0 (disable padding)

   Observe output shape per parameter setting. Likely the defaults
   were tuned for larger inputs and produce pathological small-block
   behaviour on 2kb queries.

### Acceptance criteria

- Root cause identified: parameter, lacing bug, or POA bug.
- If parameter: either update defaults or document the right settings
  for small-region queries.
- If code bug: fix and add a regression test.
- Post-fix: graph length on yeast235 `408000-410000` query should be
  ≤ 3× input length (i.e. ≤ ~1.5Mb for 460kb input; realistically
  ~50–100kb for the clean-merged case).

### Deliverables

- Either a one-line parameter fix with a comment explaining why, or
  a targeted code patch to `src/smooth.rs`.
- `notes/SMOOTHXG_INFLATION.md` — writeup of what went wrong and how
  it was fixed.

---

## 3. New partitioned-POA engine (`GfaEngine::PartPoa`)

### Question
Can we produce a clean, fast pangenome graph for closely-related
genomes using `align → partition → per-partition POA → lace`,
bypassing seqwish and smoothxg entirely?

### Why it matters
- For closely-related genomes, the shared-syncmer structure already
  gives us natural partition boundaries (runs of colinear anchors).
- POA on each partition is a single alignment pass per block; simpler
  than seqwish's graph induction + smoothxg's block-POA.
- Gives us an independent engine to compare against seqwish and (once
  fixed) pggb.

### Inputs
- Syng query output: `Vec<HomologousInterval>`. For each interval,
  sequence bytes from the AGC / FASTA.

### Design sketch

```
for each syng query (per target range):
    intervals = query_transitive(...)
    sequences = fetch_sequences(intervals, sequence_index)

    # 1. Align all-vs-all across fetched sequences.
    paf = sweepga_align_in_memory(sequences)

    # 2. Partition by alignment connectivity.
    # Each partition is a connected component of the alignment
    # graph — sequences that share at least one alignment >= threshold.
    partitions = partition_by_alignment(paf, sequences, min_overlap_bp=500)

    # 3. Per-partition POA using SPOA.
    per_partition_gfas = []
    for partition in partitions:
        subset_sequences = [sequences[i] for i in partition.members]
        graph = spoa::build_graph(subset_sequences)
        gfa = graph.generate_gfa(names_for_partition)
        per_partition_gfas.append(gfa)

    # 4. Lace partition GFAs into one.
    # Each input sequence's path threads through multiple partitions;
    # lacing stitches them in correct order.
    final_gfa = lace_subgraphs(per_partition_gfas, sequence_order)

    # 5. Single gfaffix normalization at the end.
    output_gfa = gfaffix(final_gfa)
```

### Steps

1. **Study existing building blocks.**
   - `src/commands/partition.rs` — how partition groups are defined
     for the current partitioned pipeline.
   - `GfaEngine::Poa` (single-shot SPOA) — how SPOA is invoked in-process.
   - `commands::lace::lace_subgraphs` — how per-partition GFAs are
     stitched.
   - `partitioned_gfa_pipeline` (src/lib.rs:430) — the existing
     partitioned-pggb pipeline is structurally very close; we can
     likely fork it and replace the per-partition engine.

2. **Decide partitioning criterion for closely-related inputs.** Two
   options:
   - (a) **Sliding-window** along the reference — partition size `-p`
     in bp. Each window's sequences come from fetching that region
     from every haplotype. Simplest, used by real pggb.
   - (b) **Alignment-connectivity** — connected components of the
     pairwise alignment graph. Handles structural variants better but
     more complex.

   Start with (a) — it matches what `partitioned_gfa_pipeline` already
   does and keeps the scope bounded.

3. **Add `GfaEngine::PartPoa`.**
   - Enum variant in `src/lib.rs`.
   - Parse `partpoa` / `partpoa:N` from the CLI engine spec.
   - Implement `generate_gfa_partpoa_from_intervals` in `src/graph.rs`:
     per-partition SPOA, then lace.
   - Wire it into `dispatch_gfa_engine_with_config` in `src/lib.rs`.

4. **Default `-o gfa --gfa-engine partpoa` for syng queries?** Decide
   once we have output comparison data from steps 1 and 2 (if pggb
   turns out to be fixable, defaulting to pggb still makes sense).

5. **Integration test.** A small fixture: 3 genomes sharing a 3kb
   backbone with small indels, query 2kb region, assert:
   - Graph length ≈ 2kb + indel content (≤ 3kb).
   - Exactly 3 paths (one per input).
   - Each path is continuous (no fragmentation).

### Acceptance criteria

- On yeast235 `408000-410000` query:
  - Wall time ≤ 15s (competitive with `seqwish` engine).
  - 269 paths (one per input), no fragmentation.
  - Graph length on the order of 25–100kb (≤ 5× the 23kb seqwish
    baseline — we'd expect some inflation because POA adds alternate
    alleles explicitly as graph nodes where seqwish collapses them).
- Integration test on synthetic 3-genome fixture passes.
- Output is consumable by odgi / vg without errors.

### Deliverables

- `src/graph.rs`: `generate_gfa_partpoa_from_intervals`.
- `src/lib.rs`: `GfaEngine::PartPoa` wiring.
- `src/main.rs`: CLI parsing for the new engine name.
- Integration test in `tests/test_syng_integration.rs`.
- `notes/SYNG_TRANSITIVE_DESIGN.md` — update the design doc to include
  the new engine option.

---

## Recommended order

1. **Start with #1** (reference-alignment validation). This is the
   ground truth we need before judging whether #2 or #3 fixes the
   right problem. If the syng query disagrees with PAF-based query
   in unexpected ways, #2 and #3 won't help — the syng pipeline has
   more work. If agreement is tight, we can trust the syng output and
   focus on the downstream graph engine.

2. **Then #2 and #3 in parallel.** They're independent tools for the
   same problem (better graph output). #2 fixes the existing pggb
   path; #3 adds a new path. Compare all three (seqwish, fixed pggb,
   new partpoa) on the same query once both are done.

3. **Final comparison writeup.** `notes/GFA_ENGINE_COMPARISON.md` —
   a table of wall time, output quality, and correctness for each
   engine on the yeast235 query. Use this to decide the default engine.

---

## Open assumptions

- We assume the external sweepga CLI is available for #1. If not, an
  in-process library path exists via the existing `sweepga_align`
  function in `src/commands/align.rs`.
- We assume real smoothxg binary is available for #2 comparison.
  If not, we can bisect parameters on the embedded version directly.
- We assume SPOA (via `spoa-rs`) can handle 270 × 2kb sequences per
  partition without memory blowup. Previous benchmarks suggest yes.
