# Read-level syng genotyping: COSIGT/LikeGT baseline and HPRCv2 pilot

## Scope and provenance

This is a source audit, a completed dictionary-mapping smoke test, and a proposed
benchmark. It does not implement or validate a new genotype scorer.

**Design update:** use the compact sample-index specification in
[`genotyping.typ`](genotyping.typ) for new work.
The read-by-candidate table below records an earlier proposal, superseded by
BWT-counted contextual features. The source audit, completed mapping pilot,
and truth/hold-out controls still apply; sample MEM-BWT construction and
matching-coverage queries replace the proposed per-read evidence store.

- impg: `6dd2e73`, installed version 0.5.0.
- Local LikeGT: `/home/erikg/likegt`, clean `main` at `58b91ad`.
- COSIGT source inspected at `821f2a8`, including
  `cmd/cosigt/main.go` and `cosigt_smk/workflow/rules/gafpack.smk`.
- [COSIGT workflow documentation](https://davidebolo1993.github.io/cosigtdoc/workflow/workflow.html)
- [COSIGT paper methods](https://doi.org/10.1186/s13059-026-04242-4)
- [COSIGT repository](https://github.com/davidebolo1993/cosigt)
- [LikeGT repository](https://github.com/pangenome/likegt)

## Keep COSIGT as the baseline

COSIGT compares sample node coverage with sums of candidate path copy-number
vectors, using masked/weighted cosine similarity. Repeated path visits carry
relative copy-number information. It intentionally retains useful multimapping:
its current workflow invokes `gafpack --len-scale --weight-queries`. The paper
specifies aligned bases divided by node length, with read contributions divided
across alternative alignments. Coverage is an effective structural genotype
summary; read-level models should be tested against it, not assumed superior.

LikeGT implements the same cosine objective, but its active data pipelines are
not interchangeable with one another or automatically identical to COSIGT:

- `src/commands/geno.rs:565`: sample coverage invokes `gafpack --len-scale`,
  without `--weight-queries`. Default minimap2 realignment uses `--secondary=no`
  (`:477` onward). Freeze these choices for comparisons.
- `src/commands/geno.rs:597-755`: the active external-read scorer aligns columns
  by node header, applies optional masks/weights, and enumerates all unordered
  genotype combinations with repetition. `--top` limits output, not search.
- `src/commands/hold2out.rs:1653`: reported `correct` is `cosine > 0.95`, with
  rank fixed to 1. This is not truth-genotype concordance. The hold-out path
  also uses older positional matrix parsing and retains the original graph.
- A geometric QV derived from `1-cosine` is not a calibrated genotype posterior.
- Pin odgi/gafpack and verify unequal-length, repeated-node and multimapping
  fixtures before declaring vector equivalence. The Rust scorer preserves the
  units emitted by those external tools; a CLI flag alone is not a complete
  accounting specification.

The existing local `target/release/likegt` binary is older than the source CLI:
its help does not list `geno`. Most external pipeline tools were not on the
current PATH. No LikeGT end-to-end comparison has been run in this session.

## First new model: read-by-haplotype evidence

Edge support remains useful, but need not be a prerequisite to the first
read-level window scorer. Keep the COSIGT/node-pack route as a baseline and
extract a sparse evidence table:

```
fragment/read ID, window ID, candidate source path + interval,
orientation, matched anchors, union of matched syncmer bases,
longest offset-consistent chain, chain breaks/gaps,
alternative occurrence count, unanchored bases, score
```

Existing starting points:

- `src/syng.rs:4022`, `SyngIndex::gbwt_mems_for_walk`, matches signed syncmer
  steps and inter-anchor bp offsets, returning GBWT MEM intervals.
- `src/commands/infer.rs:959`, `add_mem_hits`, compares MEM subruns to candidate
  occurrences; `build_read_walk_evidence` currently aggregates their support
  after local cosine candidate/result retention.

These are exact matches in the **anchor-and-spacing representation**, not a
base-level matching-statistics array. Matching spacing across a gap does not
prove the intervening DNA matches. Define this distinction in exported data.

For each candidate haplotype, measure compatibility within its window and retain
alternative repeated occurrences. A longest match against *any* panel path is
not enough for genotype inference. Do not sum overlapping MEM/anchor spans as
independent matched bases, or reward a candidate merely for having more possible
placements. Include shorter compatible subruns rather than allowing one strong
panel/self match to hide all alternatives.

A subsequent fragment-mixture scorer can use

```
log P(reads | genotype, window)
  = sum_fragments log(sum_haplotype_copies pi_h * P(fragment | h, window)
                      + background/error component)
```

The weights/normalizers must follow the read-generation model, including valid
fragment start positions, length and copy number. Equal 1/2 weights are only a
simplifying assumption, not universally correct for unequal-length alleles.
Matched-anchor statistics initially give scores, not calibrated probabilities;
calibration/error modeling and an explicit unknown/background option are needed.
Retain a depth/count component when evaluating absolute CN. Test fragment pairs
jointly and avoid double-counting the same reads as independent node, edge and
walk likelihoods or across overlapping windows.

Keep qualities or a locator to original FASTQ available for local realignment.
With this index's k=63, reads shorter than 63 bp cannot provide exact anchors;
error/divergence sensitivity must be measured, especially for ancient DNA.

Use homologous source-coordinate intervals from impg for windows. Numeric
coordinates on different haplotypes are not interchangeable. Global-to-local
GFA translation is a separate requirement for graph-node evidence/GAF export;
source-haplotype inference can start without materializing a local graph.

## Available HPRCv2 inputs

Prefix: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng`

- `.1khash`: about 5.1 GiB; 272,213,874-node universe reported by mapper.
- `.1gbwt`: about 8.9 GiB; `.pstep` and `.spos` are also present.
- Metadata: k=63, s=8, seed=7.
- Names table: 38,790 paths, 234 sample names and 466 sample/haplotype pairs,
  including reference entries. Total source sequence length ~1.403 Tbp.
- Sequence archive: `/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc`.
- Small real-read pilot: `hg002_pilot_1k_R{1,2}.fq`, 1,000 pairs of 148bp reads.
- Prior locus smoke outputs exist for C4/RCCX and amylase. Treat these as
  extraction/debug artifacts, not curated genotype truth.

The mapping pilot used the dictionary only. Presence of all full-index sidecars
is not a new validation of full GBWT coordinate queries or AGC source spelling.

## Completed pilot

Outputs: `target/experiments/hprcv2-read-pilot.LqrESX/` (ignored build artifacts).
The original HPRCv2 files were not modified. Mates were combined into a copied
FASTQ with `/1` and `/2` appended to their record IDs to avoid name collisions;
current `map` still treats them as separate reads, not joint fragments.

```
impg map -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -q reads.fq -o proj -O sample.proj --min-anchors 1 -t 8
```

Results (`map.log`, `summary.json`, `sample.proj/`):

| Measurement | Result |
|---|---:|
| Input reads / pairs | 2,000 / 1,000 |
| Reads with >=1 matched anchor | 1,867 (93.35%) |
| Reads with >=2 / >=3 anchors | 1,777 / 1,309 |
| Matched anchors / supported node IDs | 6,347 / 6,288 |
| Median / mean anchors per retained read | 3 / 3.40 |
| Median union of matched syncmer bases | 108 bp |
| Pairs with both / one / neither mate retained | 877 / 113 / 10 |
| Wall time, including loading and both projection passes | 24.27 s |
| Peak RSS | 12.33 GiB |
| Pack pass / GAF pass reported time, excluding loading | 0.583 / 0.037 s |

Pack size was about 31 KiB and compressed GAF about 51 KiB. This is dictionary
membership retention, not graph-placement accuracy, uniqueness, or genotype
accuracy. HG002 is represented in the panel. The tiny first-read subset is not
a representative WGS benchmark; do not extrapolate throughput from it.

## Proposed controlled benchmark

1. Extract and source-spelling-check an ordinary 10-100 kb locus, C4/RCCX, and
   an amylase/CNV locus. Existing coordinate starting points are in
   `scripts/hprcv2-syng-smoke.py`; verify assembly/region identity rather than
   treating approximate locus labels as truth.
2. Freeze candidate paths, homologous interval boundaries, graph versions and
   truth haplotypes. Simulate fixed-seed paired reads from known haplotype pairs
   and retain the FASTQs for every method. Start with 150bp at 1x, 5x and 30x,
   then add longer/error-prone reads, repeat-copy/order decoys and recombinants.
3. Separate evaluation dimensions:
   - same GFA/matrices: LikeGT/COSIGT cosine versus impg cosine, checking exact
     input units/columns first;
   - same syng index/candidates/reads: node support versus optional edge features
     versus per-read haplotype compatibility;
   - end-to-end GFA/COSIGT versus syng: report representation and read-recruitment
     differences, not just scorer differences.
4. Distinguish hold-0, candidate-only hold-out, and strict panel/index/graph
   hold-out. Excluding truth only from the genotype roster still leaves its
   sequence in the dictionary/graph. For strict held-out tests rebuild small
   training-only local indexes/graphs and keep independent truth sequences.
   Do not rewrite the global HPRCv2 index.
5. Report actual unordered genotype/equivalence-class accuracy for hold-0;
   independent sequence distance and attainable panel approximation for hold-out;
   CN error, within-window phase/order accuracy, no-call rate, candidate recall,
   score margins and runtime/RSS. Use switch error for multiwindow experiments.
   Do not use a high cosine threshold as an accuracy label.

First deliverable should be the auditable read-by-candidate table plus a
fixed-truth window benchmark, not a claim of whole-genome inference. Reuse the
existing beam only after validating local read evidence, and retain the
chromosome/coordinate fixes identified in [evidence-notes.md](evidence-notes.md) as
prerequisites for genome-scale sequence output.
