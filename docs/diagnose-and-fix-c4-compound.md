# C4 compound induction diagnosis

Task: `diagnose-and-fix`

Date: 2026-06-02

## Exact f728 region artifacts

The exact first C4 region selected by `implement-true-compound` / `f728b40`
was preserved from:

`/home/erikg/impg/data/implement_true_compound_ranked_20260602T115416Z/debug/graph_build_0000/combined.fa`

Before-fix diagnostic artifacts:

`/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_before_f728/debug/`

After-fix diagnostic artifacts:

`/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_after_fix/debug/`

Before/after summary:

`/home/erikg/impg/data/diagnose_c4_compound_20260602T_exact_metrics_summary.tsv`

Each debug directory contains `combined.fa`, `path_spellings.tsv`, `raw.paf`,
`filtered.paf`, `final.paf`, `raw.metrics.tsv`, `raw.coverage.tsv`,
`filtered.metrics.tsv`, `filtered.coverage.tsv`, `final.metrics.tsv`,
`final.coverage.tsv`, `seqwish.gfa`, and `gfa.metrics.tsv`.

## Stage diagnosis

Filtering was not the failing stage for the f728 run. The saved f728
`raw.paf`, `filtered.paf`, and `final.paf` are equivalent for this no-filter
run and all have the same gap: 9 of 45 path spellings have zero alignment
coverage.

The failure was raw wfmash parameterization for this local replacement
context. The short 3.8 kb allele class was homologous input sequence, but
raw wfmash did not map it under the whole-genome-style identity default. In
the exact diagnostic, keeping the local no-filter many-to-many mapping
multiplicity and using a local replacement identity default of 70% recovered
full coverage.

## Exact f728 metrics

| run | stage | records | pairs | covered seqs | zero seqs | mean coverage |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| before | raw | 689 | 499 | 36/45 | 9 | 0.800000 |
| before | filtered | 689 | 499 | 36/45 | 9 | 0.800000 |
| before | final | 689 | 499 | 36/45 | 9 | 0.800000 |
| after | raw | 880 | 690 | 45/45 | 0 | 1.000000 |
| after | filtered | 880 | 690 | 45/45 | 0 | 1.000000 |
| after | final | 880 | 690 | 45/45 | 0 | 1.000000 |

| run | GFA segments | GFA bp | singleton bp | shared bp |
| --- | ---: | ---: | ---: | ---: |
| before | 144 | 64206 | 34049 | 30157 |
| after | 123 | 30186 | 7 | 30179 |

The final PAF was good after the change and seqwish induction produced a
condensed graph, so seqwish parameterization was not the root failure.

## Code change

Local wfmash replacement induction now:

- writes local FASTA itself but still creates the aligner through SweepGA's
  upstream `create_aligner_adaptive`;
- derives local raw `-n` from SweepGA's upstream `parse_filter_mode`;
- treats no-filter and many-to-many local replacement induction as unbounded
  over the local sequence count;
- uses explicit local replacement identity default `70` when the user did not
  set `--sweepga-map-pct-identity`;
- preserves user identity overrides;
- writes raw, filtered, final PAF coverage metrics and path spelling TSVs in
  replacement debug directories.

No graph-quality rejection gate was added. Path sequence preservation remains
the hard replacement gate.

## Full C4 rerun

Full C4 after-fix artifacts:

`/home/erikg/impg/data/diagnose_c4_compound_20260602T_full_after_fix/`

The full run completed 8/8 rounds and wrote:

`/home/erikg/impg/data/diagnose_c4_compound_20260602T_full_after_fix/c4_full_after_fix.gfa`

Current HEAD's largest-mode first selected region is `source=top-level sites=1`
rather than the older f728 `source=compound-run sites=106` selector. That
selected round-1 SweepGA/seqwish replacement no longer expands before lacing:

- raw/final PAF: 10092 records, 6351 undirected pairs, 117/117 covered;
- seqwish debug GFA: 570 segments, 32994 bp, 93 singleton bp;
- local objective log: input 44848 bp and 4972 singleton bp, output 33278 bp
  and 376 singleton bp.

Remaining C4 quality movement after the first fixed replacement is from later
POASTA rounds and global objective diagnostics, not from pre-lacing expansion
of the selected SweepGA/seqwish replacement graph.
