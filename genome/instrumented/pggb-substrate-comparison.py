#!/usr/bin/env python3
"""THE PGGB SUBSTRATE vs THE EXPORT SUBSTRATE, per locus (stage 3, the
comparison measurement): where do the two graph families differ
structurally at the chrMT/chrI loci — the repeat/subtelomeric
localities the export path was weakest at (the foreign/repeat holders,
the seam/repeat-interior residual class)?

Per locus of the component's committed normalized-rule receipt
(realign-exhaustive-territory-<C>):
  * the truth folds and the called winner folds, from the receipt's
    own fold identities (truth_folds / best_fold_indices);
  * the EXPORT sharing between the truth rows and the winner rows:
    the syncmer nodes (signed ids, gap splices excluded) shared by
    their P lines in the raw syng range GFA of the locus's partition
    — exact-substring sharing at syncmer granularity;
  * the PGGB sharing between the same rows: the segments shared by
    their P lines in the pggb build of the SAME partition, bp-weighted
    by the segments' own lengths — alignment-induced sharing at the
    seqwish/smoothed granularity (exact-match runs through the
    wfmash/FastGA alignments, with divergent bases absorbed into
    variant bubbles);
  * the bubble density of each graph at the partition;
  * the co-occurrence verdict: whether the truth and winner rows sit
    in the SAME partition (the seam class's geometry).

Receipt-side only; no thresholds; the sharing fractions are stated
definitions (intersection over union over the rows' own P-line node
sets, bp-weighted for the pggb side where segments carry length).
Usage: pggb-substrate-comparison.py <component>
"""
import json
import sys

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
EXPORT = f"{D}/partition-graphs"
PGGB = f"{D}/pggb-partition-graphs-chrMT-chrI"
NUM_SYNCMER_NODES = 10024605
C = sys.argv[1] if len(sys.argv) > 1 else "chrI"
RECEIPT = f"{D}/realign-exhaustive-territory-{C}.jsonl"


class Gfa:
    def __init__(self, path):
        self.path = path
        self.segs = {}
        self.paths = {}
        with open(path) as f:
            for line in f:
                fields = line.rstrip("\n").split("\t")
                if fields[0] == "S":
                    self.segs[fields[1]] = len(fields[2])
                elif fields[0] == "P":
                    self.paths[fields[1]] = [
                        s.strip() for s in fields[2].split(",") if s.strip()
                    ]

    def node_set(self, row, syncmer_only=False):
        out = set()
        for step in self.paths[row]:
            if step[-1] in "+-":
                ident = step[:-1]
            else:
                ident = step
            if syncmer_only and int(ident) > NUM_SYNCMER_NODES:
                continue
            out.add(ident)
        return out


def main():
    print(f"# the pggb substrate vs the export substrate, component {C}")
    print(
        "# locus partition truth_rank | truth rows winner rows | co-occur | "
        "export shared/union syncmer nodes | pggb shared/union bp | "
        "export bubbles | pggb bubbles | verdict"
    )
    for line in open(RECEIPT):
        d = json.loads(line)
        locus = d["locus"]
        partition = d["partition"]
        fold_ids = d["fold_identities"]
        truth_idx = d["truth_folds"]
        winner_idx = d["best_fold_indices"]
        if truth_idx is None:
            print(
                f"locus {locus}\tpartition {partition}\ttruth INEXPRESSIBLE\t"
                f"winner_rows {len(winner_idx) if winner_idx else 0}"
            )
            continue
        truth_rows = set()
        for i in truth_idx:
            for m in fold_ids[i]["members"]:
                truth_rows.add((m["path_name"], m["start"], m["end"]))
        winner_rows = set()
        for i in winner_idx:
            for m in fold_ids[i]["members"]:
                winner_rows.add((m["path_name"], m["start"], m["end"]))
        export = Gfa(f"{EXPORT}/partition{partition}.gfa")
        pggb = Gfa(f"{PGGB}/partition{partition}/graph.gfa")

        def rows_to_names(rows, gfa_paths):
            names = []
            for (p, s, e) in rows:
                name = f"{p}:{s}-{e}"
                if name in gfa_paths:
                    names.append(name)
            return names

        t_export = rows_to_names(truth_rows, export.paths)
        w_export = rows_to_names(winner_rows, export.paths)
        t_pggb = rows_to_names(truth_rows, pggb.paths)
        w_pggb = rows_to_names(winner_rows, pggb.paths)

        # export sharing: syncmer nodes only, over rows present in the
        # export partition graph
        def export_nodes(names):
            nodes = set()
            for n in names:
                nodes |= export.node_set(n, syncmer_only=True)
            return nodes

        # pggb sharing: bp-weighted segment sets
        def pggb_bp(names):
            bp = set()
            for n in names:
                for ident in pggb.node_set(n):
                    bp.add(ident)
            return bp

        t_in = len(t_export) > 0 and len(w_export) > 0
        if t_in:
            tn, wn = export_nodes(t_export), export_nodes(w_export)
            export_shared, export_union = len(tn & wn), len(tn | wn)
            ts, ws = pggb_bp(t_pggb), pggb_bp(w_pggb)
            pggb_shared_bp = sum(pggb.segs[i] for i in ts & ws)
            pggb_union_bp = sum(pggb.segs[i] for i in ts | ws)
            # bubble counts (both graphs, the partition's own structure)
            def bubbles(gfa):
                outd = {}
                ind = {}
                with open(gfa.path) as f:
                    for line in f:
                        fields = line.rstrip("\n").split("\t")
                        if fields[0] == "L":
                            outd.setdefault(fields[0 + 1], set()).add(fields[2])
                            ind.setdefault(fields[3], set()).add(fields[1])
                return sum(1 for s in outd.values() if len(s) > 1) + sum(
                    1 for s in ind.values() if len(s) > 1
                )

            export_bubbles = bubbles(export)
            pggb_bubbles = bubbles(pggb)
            verdict = ""
            ef = export_shared / export_union if export_union else 0.0
            pf = pggb_shared_bp / pggb_union_bp if pggb_union_bp else 0.0
            if ef < 0.01 and pf >= 0.5:
                verdict = "alignment-merged (export near-disjoint, pggb majority-shared)"
            elif pf > ef + 0.1:
                verdict = "pggb shares more"
            elif abs(pf - ef) <= 0.1:
                verdict = "same-sharing"
            else:
                verdict = "export shares more"
            print(
                f"locus {locus}\tpartition {partition}\ttruth_rank {d['truth_rank']}\t"
                f"truth_rows {len(truth_rows)}\twinner_rows {len(winner_rows)}\t"
                f"co-occur {int(len(t_export) > 0 and len(w_export) > 0)}\t"
                f"export {export_shared}/{export_union} = {ef:.4f}\t"
                f"pggb {pggb_shared_bp}/{pggb_union_bp} bp = {pf:.4f}\t"
                f"bubbles e{export_bubbles} p{pggb_bubbles}\t{verdict}"
            )
        else:
            print(
                f"locus {locus}\tpartition {partition}\ttruth_rank {d['truth_rank']}\t"
                f"truth_rows {len(truth_rows)}\twinner_rows {len(winner_rows)}\t"
                f"co-occur 0\trows not joint in the partition graph (the seam geometry)"
            )


main()
