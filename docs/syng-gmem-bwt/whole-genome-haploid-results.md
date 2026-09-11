# Whole-genome haploid count/thread results

Two complete in-panel haploid controls now exercise the implemented quantitative
pipeline, not only the presence-compatibility bootstrap. These are **source-bundle
and source-interval recovery results**, not calibrated allele accuracy or final
catalog acceptance. Diploid and whole-genome recombinant validation remain open.

## Frozen conditions

- Code: `9e1cc305a0d346b7c06c8f81002e677abdaa49b9`, branch
  `work/genome-mem-bwt-pipeline`; independent review found no issues.
- Parent validation: 433 library + 85 CLI unit + 41 integration tests passed.
- Normal installation: `cargo install --offline --locked --path . --force --bin impg -j4`.
  Installed executable SHA256:
  `d1ef32c9335cc60a482095c110e0e4c23b36e3d4b525a1746adf0f3c054bb303`.
- Same ACGT-only sample64 yeast panel, all 19,421 completed ownership groups and
  383,849 source occurrences. No candidate/seed/frequency/partition-policy change.
- Error-free 150bp reads, independent uniform legal starts and random orientation,
  depth10, seed1729. All17 source paths per assembly, including mitochondria;
  ambiguity-spanning reads retained. This is an artificial one-copy-per-path model,
  not a real-stock ploidy or mitochondrial-copy claim.
- Both controls use the **same declared S288C source-coordinate axis**:
  12,242,942bp, 1,307 intervals, 17 components. No truth-based candidate selection.
- Same model parameters: haploid depth10, background0.1, maximum mean deviance10,
  switch penalty10. No tuning against either truth result.
- Truth is opened after calls and threads are frozen. The reusable sample/catalog,
  axis, truth and installed binary hashes were unchanged over each genotype run.

## Results and denominators

| Measurement | S288C | SK1, fixed S288C axis |
|---|---:|---:|
| Complete sample source bp | 12,242,942 | 12,147,923 |
| Simulated reads | 816,202 | 809,870 |
| Truth-containing ownership groups | 1,218 | 1,323 |
| True source bundle among best candidates, sample-source bp | 11,908,550 (97.27%) | 11,921,183 (98.13%) |
| Unique best bundle truth-compatible, sample-source bp | 1,291,585 | 4,215,469 |
| Thread-resolved reference-axis bp | 10,740,422 (87.73%) | 10,491,799 (85.70%) |
| Resolved reference-axis intervals | 1,085 | 1,055 |
| Resolved intervals with a matching truth group | 1,085 | 1,047 |
| Exact truth-source-compatible resolved intervals | 1,083 | 1,041 |
| Exact source agreement on resolved, truth-assessable reference bp | 10,720,352 / 10,740,422 (99.81%) | 10,355,009 / 10,411,698 (99.46%) |
| Emitted source-path blocks | 50 | 66 |

The source-bundle metric permits tied best candidates. It is not a unique-genotype
accuracy statistic. Reference-axis coverage and sample-source coverage use different
coordinate systems for SK1 and must not be mixed. Eight resolved SK1-control axis
intervals have no corresponding truth group; they are **not counted as correct**
in the assessed agreement denominator. A source-identity disagreement can still
spell similar sequence; sequence-level allele equivalence has not been evaluated.

S288C disagrees at two resolved terminal intervals:
- `S288C#0#chrIV:0-10011`: predicts `DG1768#0#chrIV:6944-16902`.
- `S288C#0#chrXI:0-10059`: predicts `W303#0#chrXI:0-9995`.

SK1 has six assessed resolved disagreements: one interval at reference chrXIV
570043–580059, and five spanning reference chrXVI10044–56717. The dominant resolved
donor is SK1 (10,355,009 reference bp), not the coordinate reference S288C. This
control therefore checks more than successful recovery of the reference itself.

Twenty-two groups recur on the axis, comprising111 intervals /553,177 reference bp.
Their partition genotypes remain, but all appearances are explicitly unthreaded
so evidence is not duplicated. Other unresolved intervals include tied paths,
no positive/local evidence and unknown orientation. Unscaffolded calls are retained.

SK1 truth identifies four assessable **source-path changes** between adjacent
single-occurrence groups on the reference axis; three are recovered. These are
not a synthetic recombination benchmark or exact biological breakpoint claims.
There are nine reported breaks overall, with only three assessable source-path
changes under that restricted truth definition. S288C has no known/reported switches.

## Runtime and storage

Both controls were bounded to four CPUs (252–255), Rayon4, nice10, with no wall-clock
cutoff. Reported stage times include the CLI process; scoring/threading sub-times
include saving their outputs.

| Stage | S288C | SK1 |
|---|---:|---:|
| Sample-index wall time | 66.0s | 983.4s |
| Sample-index user+system CPU time | 65.85s | 92.43s |
| Sample-index size | 39,895,420 bytes | 39,649,676 bytes |
| Genotype/thread/evaluation total | 282.1s | 245.5s |
| Count scoring + save | 22.91s | 20.46s |
| Threading + save | 3.06s | 2.69s |
| Genotype-stage peak RSS | 39,766,016 KiB | 39,763,968 KiB |
| Calls JSON | 212,891,915 bytes | 212,839,260 bytes |

SK1 sample construction received only9% CPU over its measured wall interval. The
large wall-time difference is **not evidence of a15-fold algorithmic slowdown**;
the waiting cause is unestablished. Do not mask reads/repeats or change sampling
based on that elapsed time. Both count indexes peaked near1.28GiB process RSS.

The shared catalog remains41,129,985,798bytes of JSON. Construction took530.3s and
81,445,340KiB peak RSS. This is a scalability limitation, not the compact sample
index. Default quantitative calls omit reproducible per-feature maps; no second
expanded feature catalog is generated. Catalog input/checksum handling dominates
the genotype command's time and memory beyond its26s/23s scoring+threading work.

The original presence-only truth stage was interrupted after a demonstrated I/O
bug: `serde_json::from_reader(File)` issued33.8billion approximately one-byte reads.
Commit `40a7e1d` adds buffering. The failed attempt and completed inputs remain;
no completion/biological result is claimed for that interrupted stage.

## Evidence

S288C root:
`~/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/`
- `inputs/manifest.json`, sample/catalog manifests and resource logs.
- `quantitative-thread-20260911T125231Z/{manifest.json,summary.json,source-discrepancies.json,resources.log,results/}`.

SK1 root:
`~/yeast/genome-mem-bwt-sk1-control-20260911T130448Z/`
- `inputs/manifest.json`, `sample/manifest.json`, `sample.resources.log`.
- `quantitative/{manifest.json,summary.json,source-discrepancies.json,resources.log,results/}`.

Parent scripts and captured reviews/tests:
`target/experiments/genome-mem-bwt-pipeline/` in the original impg worktree.
The original generators/runners were retained when generalized versions were added.
Source and binary hashes, commands, inputs and output checksums are recorded.

## Remaining acceptance gates

The executable now performs count-based partition inference and source-path
threading over a whole genome. It still requires whole-genome recombinant and
diploid controls, independent/rebuilt hold-outs, source-sequence concordance and
callability assessment, calibrated count exposure/background, and joint handling
of shared/nonlocal/boundary evidence and repeated-axis copy states. Current local
factors exclude1,066,487 of11,896,777 catalog features explicitly; catalog acceptance
remains false. No real SK1/Y12 read interpretation is justified by these controls.
