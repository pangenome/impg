# C4 Current-Control Edge Artifact Diagnosis

Run root:
`/home/erikg/impg/data/c4_control_provenance_20260610T154843Z`

Target graph:
`graphs/current_53kb_pggb/C4_GRCh38_53kb.gfa`

Sorted/rendered graph inspected:
`sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa`

## Summary

The left-edge artifact is a source extraction artifact, not graph induction or
rendering. The ten visually cut-off bottom paths all start at final sorted node
`3191`, whose graph layout offset is 33,211 bp. Their final graph spellings
match the source GFA exactly, so the missing 5-prime/left sequence was already
absent from the extracted source spellings before seqwish/smoothing/gfaffix.

The bottom sorting is not caused by short sequence length alone. There are 55
extracted paths shorter than 100 kb, including several shorter than the ten
bottom paths, and most of those start near the left side of the sorted graph.
The bottom group is distinct because its first node is late in the graph, not
only because several members are short.

The right-side condensation localizes to repeated/shared graph intervals near
nodes `5915-6300`, `6312-6374`, and `6685-6820`, with a larger repeated
re-entry at `5639 -> 2966`. This is real graph structure induced from shared
source sequence/PAF alignments and then preserved by the path-preserving graph
pipeline. It is not a gfalook-only rendering artifact, and it is not explained
by the ten truncated-left query intervals.

Recommended next fix:

- For the left cut-off paths, fix query interval collection/extraction. Add a
  pre-build check that flags extracted intervals whose first source base maps to
  a late shared block or whose left flank/anchor coverage is absent.
- For the right-side condensation, do not change interval collection first. If
  the desired output is visually less condensed, change occurrence-aware graph
  construction/smoothing or visualization/unrolling. The evidence points to
  PAF-supported shared/repeated source sequence, not accidental query
  overextension or a sort/render-only problem.

## Commands And Evidence

The original current-control command extracted one interval set from the syng
query and wrote the pggb-style graph:

```bash
/usr/bin/timeout 2400s target/release/impg query \
  -a /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.syng \
  -b /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/c4_53kb.bed \
  -d 100k \
  --sequence-files /home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc \
  --debug-dir /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/debug/current_53kb_pggb_query \
  -o gfa:pggb \
  -O /home/erikg/impg/data/c4_control_provenance_20260610T154843Z/graphs/current_53kb_pggb \
  -t 32 -v 1
```

The query BED row was:

```text
GRCh38#0#chr6	31982056	32035418	C4_GRCh38_53kb
```

Path preservation was checked against the source GFA made from the extracted
FASTA:

```text
expected_paths	466
observed_paths	466
missing_paths	0
extra_paths	0
spelling_mismatches	0
```

This exact spelling check is the key separator: if a final path lacks the left
shared prefix, the source spelling from `combined.fa`/`path_spellings.tsv`
lacked it too.

The pggb-style query log also shows the graph-processing sequence:

```text
Collected 466 merged intervals
Saved combined FASTA to .../debug/current_53kb_pggb_query/combined.fa
Saved raw PAF to .../debug/current_53kb_pggb_query/raw.paf
Saved filtered PAF to .../debug/current_53kb_pggb_query/filtered.paf
[smooth] Pass 1/2 ...
[smooth] Pass 2/2 ...
[gfaffix] found and collapsed 19 shared prefixes ...
```

The graph report supports the presence of repeat/tangle structure in the final
graph:

```text
segments=7252
paths=466
path_steps=3575249
total_segment_bp=89353
path_white_space_bp_max=33593
path_white_space_bridges_ge_threshold=3102
duplicate_sequence_frac=0.820877
local_repeat_context_nodes=355
local_repeat_context_occurrences=559
povu_sites=1618
povu_leaf_sites=1603
```

The one-off parsers used for this report read the sorted GFA `S` line order as
the graph layout coordinate system, then inspected each `P` line's first/last
node and large coordinate jumps. No generated GFA, PAF, PNG, or log artifacts
were created.

Diagnostic command shapes:

```bash
RUN=/home/erikg/impg/data/c4_control_provenance_20260610T154843Z
GFA=$RUN/sorted/current_53kb_pggb.C4_GRCh38_53kb.Ygs.gfa
PAF=$RUN/debug/current_53kb_pggb_query/final.paf

# First/last node and graph-coordinate starts for all paths.
python - "$GFA" <<'PY'
from pathlib import Path
import sys
gfa = Path(sys.argv[1])
seg_len, seg_start, paths = {}, {}, {}
offset = 0
for line in gfa.open():
    if line.startswith("S\t"):
        _, sid, seq, *_ = line.rstrip("\n").split("\t")
        seg_start[sid] = offset
        seg_len[sid] = len(seq)
        offset += len(seq)
    elif line.startswith("P\t"):
        parts = line.rstrip("\n").split("\t")
        paths[parts[1]] = [step[:-1] for step in parts[2].split(",")]
for rank, (name, nodes) in enumerate(paths.items(), 1):
    first, last = nodes[0], nodes[-1]
    if seg_start[first] >= 10000:
        length = sum(seg_len[n] for n in nodes)
        print(rank, name, length, first, seg_start[first], last,
              seg_start[last] + seg_len[last], len(nodes), sep="\t")
PY

# Large path jumps in sorted graph coordinates used the same parser; for each
# adjacent path step a,b, aggregate transitions where:
#   gap = seg_start[b] - (seg_start[a] + seg_len[a])
#   abs(gap) >= 1000
python - <<'PY'
# Summarized here to avoid a long duplicated script in the report.
PY

# PAF support scanned final.paf for qname or tname matching:
#   GRCh38#0#chr6:31949370-32068272
# and overlapping pathbp 105500-109000 or 111500-113200.
python - <<'PY'
# Summarized here to keep the report concise; the counts are reported below.
PY
```

## Left-Edge Cut-Off Paths

These are the only paths whose first node starts at graph bp 33,211 or later.
They are exactly the bottom ten sorted paths, sorted ranks 457-466.

| Sorted rank | Path / extracted source interval | Source length | First node / path bp | Last node / path bp | Source spelling truncated before induction? |
|---:|---|---:|---|---|---|
| 457 | `HG01952#2#CM087999.1:31995790-32075638` | 79,848 | `3191` @ graph 33,211 / path 0 | `7225` @ graph end 86,015 / path end 79,848 | Yes; extracted spelling starts at the late shared block. |
| 458 | `HG01960#2#JBHIHN010000011.1:32078636-32158486` | 79,850 | `3191` @ graph 33,211 / path 0 | `7215` @ graph end 86,005 / path end 79,850 | Yes; extracted spelling starts at the late shared block. |
| 459 | `HG02027#2#JBHDTY010000008.1:5188422-5268271` | 79,849 | `3191` @ graph 33,211 / path 0 | `7225` @ graph end 86,015 / path end 79,849 | Yes; extracted spelling starts at the late shared block. |
| 460 | `HG02132#2#CM086908.1:31997546-32077395` | 79,849 | `3191` @ graph 33,211 / path 0 | `7228` @ graph end 86,019 / path end 79,849 | Yes; extracted spelling starts at the late shared block. |
| 461 | `HG02178#2#CM089938.1:31943310-32055920` | 112,610 | `3191` @ graph 33,211 / path 0 | `7186` @ graph end 85,936 / path end 112,610 | Yes for left-end coverage; total length is not globally shortest, but the left shared block is absent. |
| 462 | `HG03710#2#CM086882.1:31974825-32061063` | 86,238 | `3191` @ graph 33,211 / path 0 | `7191` @ graph end 85,945 / path end 86,238 | Yes; extracted spelling starts at the late shared block. |
| 463 | `HG03874#1#CM089432.1:31915349-32027959` | 112,610 | `3191` @ graph 33,211 / path 0 | `7182` @ graph end 85,932 / path end 112,610 | Yes for left-end coverage; total length is not globally shortest, but the left shared block is absent. |
| 464 | `NA18565#2#CM094342.1:31929266-32009115` | 79,849 | `3191` @ graph 33,211 / path 0 | `7225` @ graph end 86,015 / path end 79,849 | Yes; extracted spelling starts at the late shared block. |
| 465 | `NA18945#2#CM101630.1:32018016-32097865` | 79,849 | `3191` @ graph 33,211 / path 0 | `7225` @ graph end 86,015 / path end 79,849 | Yes; extracted spelling starts at the late shared block. |
| 466 | `NA18960#2#CM101559.1:31867170-31947019` | 79,849 | `3191` @ graph 33,211 / path 0 | `7228` @ graph end 86,019 / path end 79,849 | Yes; extracted spelling starts at the late shared block. |

The extracted length distribution from
`debug/current_53kb_pggb_query/path_spellings.tsv` was:

```text
paths=466
min_length=79797
max_length=184376
mean_length=114019.7
paths_lt_100kb=55
```

The shortest paths are not all in the bottom block. For example:

| Sorted rank | Path | Length | First node @ graph bp |
|---:|---|---:|---|
| 263 | `HG00099#2#CM087361.1:31928455-32008252` | 79,797 | `65` @ 64 |
| 265 | `HG00146#1#CM090015.1:31922607-32002404` | 79,797 | `65` @ 64 |
| 270 | `HG00272#2#CM094216.1:31931690-32011487` | 79,797 | `65` @ 64 |
| 271 | `HG00280#2#CM087160.1:31931045-32010842` | 79,797 | `65` @ 64 |
| 292 | `HG00738#1#CM086647.1:31951624-32031421` | 79,797 | `65` @ 64 |
| 411 | `HG02258#2#CM085837.1:31934546-32014343` | 79,797 | `66` @ 65 |
| 424 | `HG06807#2#CM101434.1:32040804-32120601` | 79,797 | `66` @ 65 |
| 439 | `NA19468#2#JBHDXY010000045.1:31915811-31995608` | 79,797 | `68` @ 67 |
| 190 | `NA18879#1#CM099716.1:31958945-32038743` | 79,798 | `60` @ 59 |
| 338 | `HG03017#2#CM088928.1:31927766-32007564` | 79,798 | `65` @ 64 |

So the visual bottom grouping means "late first node" rather than "short path"
by itself.

## Right-Side Condensation

The strongest final-graph coordinate jumps are path transitions, not renderer
paint artifacts. The table below gives the transition, the graph-coordinate
gap, how many paths contain it, and the equivalent position on the GRCh38
source path `GRCh38#0#chr6:31949370-32068272`.

| Transition | Graph end/start | Gap | Paths affected | GRCh38 path bp | GRCh38 source coord |
|---|---:|---:|---:|---:|---:|
| `5639 -> 2966` | 65,133 -> 31,529 | -33,604 | 419 paths, 467 occurrences | 63,899 | 32,013,269 |
| `6300 -> 5915` | 75,507 -> 71,814 | -3,693 | 465 paths | 106,516 | 32,055,886 |
| `6254 -> 6301` | 74,213 -> 75,507 | +1,294 | 466 paths | 108,435 | 32,057,805 |
| `6817 -> 6312` | 79,567 -> 75,793 | -3,774 | 456 paths | 112,057 | 32,061,427 |
| `6374 -> 6685` | 76,039 -> 78,046 | +2,007 | 466 paths | 112,282 | 32,061,652 |
| `6688 -> 6820` | 78,362 -> 79,569 | +1,207 | 463 paths | 112,597 | 32,061,967 |

For the representative GRCh38 path, the right-side repeated regions are:

```text
nodes 5915-6254: pathbp 103335-105237 and 106516-108435
nodes 6312-6374: pathbp 108718-108943 and 112057-112282
nodes 6685-6817: pathbp 110930-112057 and 112282-112597
nodes 6820-7224: pathbp 112597-118902
```

The final graph spelling for GRCh38 exactly reconstructs the extracted
`combined.fa` sequence:

```text
Combined FASTA GRCh38 length 118902
graph-spelling matches source True
```

The PAF also supports these right-side windows as shared sequence. Scanning
`debug/current_53kb_pggb_query/final.paf` for GRCh38 windows gave:

```text
GRCh38 pathbp 105500-109000: q_records=465, t_records=465
GRCh38 pathbp 111500-113200: q_records=465, t_records=465
```

Representative records include near-full-length or right-tail alignments:

```text
GRCh38:0-118882 -> HG00350#2#JBIREM010000045.1:0-118899
  nmatch/aln=118778/118905 dv=0.0008

GRCh38:0-118886 -> HG02559#2#CM088248.1:0-118902
  nmatch/aln=118877/118902 dv=0.0000

HG03927#1#CM086636.1:111051-118955 -> GRCh38:110930-118825
  nmatch/aln=7885/7904 dv=0.0024
```

This evidence points to graph induction over real shared/repeated input
sequence:

- Not query/extraction: the right-side windows are present in normal-length
  source paths, including GRCh38 and CHM13, and are covered by many PAF
  records.
- Not sorting/rendering only: the jumps are explicit consecutive node
  transitions in the final GFA path lines.
- Not primarily gfaffix: the localized intervals span thousands of graph bp and
  hundreds of paths, while gfaffix reported only 19 shared-prefix collapses.
- Best classification: source/PAF-supported repeat condensation introduced by
  graph induction and retained through smoothing/gfaffix, then made visually
  prominent by sorted graph layout and rendering.

## Validation Checklist

- Suspect truncated paths are named explicitly in the left-edge table.
- Source extraction is distinguished from induction/rendering by
  `compare_gfa_paths` exact spelling validation and `path_spellings.tsv`
  source lengths.
- The right-side condensation is localized by node transition, graph bp, path
  count, GRCh38 path bp, and GRCh38 source coordinate.
- The recommendation separates the fixes: query interval collection for the
  left cut-off group; occurrence-aware graph construction/smoothing or
  visualization for the right-side condensed structure.
