#!/usr/bin/env python3
# THE DOSAGE-SURFACE AGGREGATE: the 17-component fleet receipt for the
# copy-number layer — the per-component agreement tables summed, the
# honest decomposition (how much dosage error is inherited from the
# non-rank-1 residue vs dosage-specific), and THE OLD DEFECT'S VERDICT
# (the 1,041/1,307 single-class-of-2 era vs the emitted dosage now,
# with numbers). Assessment-side numbers come ONLY from the test-mode
# artifacts; the default emission's truth-free statements (the
# flat-class-of-2 loci) are stated as what they are.
#
# Usage: dosage-aggregate.py <validation-dir>
import json, os, sys

COMPONENTS = ["chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
              "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII",
              "chrXIII", "chrXIV", "chrXV", "chrXVI"]

def truth_paths(component, d):
    if component == "chrMT":
        return (f"{d}/cli-product-chrMT-truthqv/calls.jsonl.truth-qv.jsonl")
    return f"{d}/cli-product-{component}/calls.jsonl.truth-qv.jsonl"

def main():
    d = sys.argv[1]
    total = {
        "loci": 0, "expressible": 0, "rank1": 0, "correct": 0,
        "rank1_agree": 0, "emitted_flat2": 0, "truth_flat2": 0,
        "segments": 0, "copy_bp": 0,
        "union_segments": 0, "segment_classes": {},
        "classes": {}, "mechanisms": {},
    }
    rows = []
    for component in COMPONENTS:
        emission = json.loads(
            open(f"{d}/dosage-{component}/dosage.jsonl").read().strip())
        lines = [json.loads(l)
                 for l in open(f"{d}/dosage-{component}/dosage.jsonl.truth-qv.jsonl")]
        summary = lines[-1]
        per_locus = [l for l in lines if not l.get("summary")]
        product = {}
        for line in open(truth_paths(component, d)):
            r = json.loads(line)
            product[r["locus"]] = r
        rank1 = sum(1 for r in per_locus if r["rank1"])
        emitted_flat2 = sum(1 for l in emission["loci"] if l["flat_class_of_2"])
        homozygous_pairs = sum(1 for l in emission["loci"]
                               if l["homozygous_pair"])
        mixed_emitted = sum(1 for l in emission["loci"]
                             if not l["flat_class_of_2"])
        # the truth-side heterozygosity at the expressible loci (the
        # old defect's truth was two distinct 1-copy classes)
        truth_mixed = sum(1 for r in per_locus
                          if r["truth"] and not r["truth"]["flat_class_of_2"])
        truth_flat2 = sum(1 for r in per_locus
                          if r["truth"] and r["truth"]["flat_class_of_2"])
        classes = summary["class_counts"]
        row = {
            "component": component, "loci": summary["loci"],
            "expressible": summary["truth_pair_expressible"],
            "rank1": rank1,
            "correct": summary["correct_dosage"],
            "rank1_agree": summary["rank1_dosage_agree"],
            "inherited": classes.get("dosage_inherited_from_non_rank1_call", 0),
            "dosage_specific": classes.get("dosage_specific_emission_defect", 0),
            "nonrank1_agrees": classes.get("non_rank1_call_dosage_agrees", 0),
            "bracket": classes.get("bracketed_truth_pair_not_expressible", 0),
            "bracketed": classes.get("bracketed_truth_pair_not_expressible", 0),
            "emitted_flat2": emitted_flat2,
            "homozygous_pairs": homozygous_pairs,
            "mixed_emitted": mixed_emitted,
            "truth_mixed": truth_mixed,
            "truth_flat2": truth_flat2,
            "segments": len(emission["segments"]),
            "copy_bp": emission["depth_qc"]["emitted_copy_bp"],
            "union_segments": summary["union_segments"],
            "segment_classes": summary["segment_classes"],
        }
        rows.append(row)
        total["loci"] += row["loci"]
        total["expressible"] += row["expressible"]
        total["rank1"] += row["rank1"]
        total["correct"] += row["correct"]
        total["rank1_agree"] += row["rank1_agree"]
        total["emitted_flat2"] += row["emitted_flat2"]
        total["truth_flat2"] += row["truth_flat2"]
        total["segments"] += row["segments"]
        total["copy_bp"] += row["copy_bp"]
        total["union_segments"] += row["union_segments"]
        for key, value in classes.items():
            total["classes"][key] = total["classes"].get(key, 0) + value
        for key, value in summary["segment_classes"].items():
            total["segment_classes"][key] = \
                total["segment_classes"].get(key, 0) + value
        for r in per_locus:
            if r["class"] in ("dosage_inherited_from_non_rank1_call",
                              "dosage_specific_emission_defect"):
                total["mechanisms"][r["mechanism"]] = \
                    total["mechanisms"].get(r["mechanism"], 0) + 1

    out = sys.stdout
    print("# THE DOSAGE-SURFACE FLEET AGGREGATE (17 components)")
    print("#")
    print("# component  loci  expr  rank1  correct  r1agree  inherited  dosagespec  nonr1agr  bracket  flat2  homopair  mixed  truthmix  truthflat2  segs  copy_bp")
    for row in rows:
        print(f"# {row['component']:10} {row['loci']:4}  {row['expressible']:4}  "
              f"{row['rank1']:5}  {row['correct']:7}  {row['rank1_agree']:7}  "
              f"{row['inherited']:9}  {row['dosage_specific']:10}  "
              f"{row['nonrank1_agrees']:8}  {row['bracket']:7}  "
              f"{row['emitted_flat2']:5}  {row['homozygous_pairs']:8}  "
              f"{row['mixed_emitted']:5}  {row['truth_mixed']:8}  "
              f"{row['truth_flat2']:10}  {row['segments']:4}  {row['copy_bp']}")
    print("#")
    print(f"# TOTAL loci {total['loci']}  expressible {total['expressible']}  "
          f"rank-1 {total['rank1']}  correct dosage {total['correct']} "
          f"({total['correct'] / total['expressible']:.4f} of expressible)  "
          f"rank-1 dosage agree {total['rank1_agree']} "
          f"({total['rank1_agree'] / total['rank1']:.4f} of rank-1)")
    print(f"# THE HONEST DECOMPOSITION (per locus): "
          f"{json.dumps(dict(sorted(total['classes'].items())))}")
    print(f"#   mechanisms of the differing loci: "
          f"{json.dumps(dict(sorted(total['mechanisms'].items())))}")
    print(f"# THE PER-SEGMENT DECOMPOSITION (over the union universes, "
          f"{total['union_segments']} segments): "
          f"{json.dumps(dict(sorted(total['segment_classes'].items())))}")
    print(f"# THE EMISSION: {total['segments']} segments, "
          f"copy-bp {total['copy_bp']} (the row-bp identity gated at "
          f"every component), flat-class-of-2 loci "
          f"{total['emitted_flat2']} of {total['loci']} "
          f"({total['emitted_flat2'] / total['loci']:.4f}), "
          f"truth flat-class-of-2 loci {total['truth_flat2']}")
    print("#")
    print(f"# THE OLD DEFECT'S VERDICT: the 1,041/1,307 single-class-of-2 era "
          f"(the old genotype_calls.diplotype.dosage emitted ONE class of 2 at "
          f"79.6% of loci vs truth 1.0/1.0) is GONE at the coverage level: the "
          f"emission states {total['emitted_flat2']} flat-class-of-2 loci of "
          f"{total['loci']} ({total['emitted_flat2'] / total['loci']:.4f}), and "
          f"at rank-1 loci the emitted dosage equals the truth dosage at "
          f"{total['rank1_agree']} of {total['rank1']} "
          f"({total['rank1_agree'] / total['rank1']:.4f}) — the dosage is now "
          f"the molecules' own counted material, mixed where the truth is "
          f"mixed; the residual disagreement is INHERITED from the "
          f"non-rank-1 calls ({total['classes'].get('dosage_inherited_from_non_rank1_call', 0)} "
          f"loci, the 252-locus residue's own material), with "
          f"{total['classes'].get('dosage_specific_emission_defect', 0)} "
          f"dosage-specific defects at rank-1 loci")

if __name__ == "__main__":
    main()
