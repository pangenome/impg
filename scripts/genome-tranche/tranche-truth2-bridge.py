#!/usr/bin/env python3
"""THE GENOME TRANCHE - the truth2 bridge (the axis-bridge machinery's
HG002#2 leg, the committed recipe VERBATIM with the maternal target in
place of the paternal one): the owner's truth pair at a territory is
(the HG002#1 axis row - prep/<T>/axis-path.txt, ALREADY the committed
bridge output for hap1 - and the HG002#2 row of the refined bed). At
1,270 of 5,844 accepted territories the refined bed carries MORE THAN
ONE HG002#2 row (the cross-chromosomal segdup-hit class the split's
own axis-row rule receipts), so the first-row convention of the
five-locality A/B would name a WRONG truth at exactly the segdup
territories the tranche most needs - the fix is the bridge itself:
each territory's window (the CHM13#0 bridge row of locus.bed) aligned
asm5 to HG002#2#chr1-22 (extracted from the panel AGC, the same
recipe as hg002-axis-target.fa), and the truth2 row = the refined
bed's same-contig HG002#2 row with MAX OVERLAP vs the bridge interval
(the split's own axis-row rule, applied to the maternal haplotype).
Fail-closed classes, named: no-bridge (no primary same-chrom row),
ambiguous-bridge (>1 row at max overlap), no-truth2-row (the bridge
interval overlaps no refined-bed HG002#2 row).
usage: tranche-truth2-bridge.py <WT> <RV> [nproc]
"""
import glob
import os
import subprocess
import sys

WT = sys.argv[1]
RV = sys.argv[2]
NPROC = int(sys.argv[3] or 16)
AGC = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc"
AGC_BIN = "/home/erikg/agc/bin/agc"
MAPFILE = "/home/erikg/hprcv2/generic-sampling/agc-path-sample-map.txt"
BRIDGE = f"{RV}/rv-axis-bridge.py"


def load_map():
    m, cur = {}, None
    for line in open(MAPFILE):
        if line.startswith(" "):
            m.setdefault(line.strip().split()[0], cur)
        else:
            cur = line.strip().split()[0]
    return m


def extract_target():
    """HG002#2#chr1-22 from the AGC (the hg002-axis-target.fa recipe)."""
    out = f"{WT}/bridge/hg002-2-target.fa"
    if os.path.exists(out) and os.path.getsize(out) > 2.5e9:
        return out
    os.makedirs(f"{WT}/bridge", exist_ok=True)
    sample_of = load_map()
    specs = []
    for c in list(range(1, 23)):
        path = f"HG002#2#chr{c}"
        specs.append(f"{path}@{sample_of[path]}")
    with open(f"{WT}/bridge/target.specs", "w") as f:
        f.write("\n".join(specs) + "\n")
    # the AGC writes the path@sample header; rewrite to the plain path
    subprocess.run(
        f"xargs -a {WT}/bridge/target.specs {AGC_BIN} getctg {AGC} "
        f"> {WT}/bridge/target.raw.fa",
        shell=True, check=True)
    with open(out, "w") as dst:
        for line in open(f"{WT}/bridge/target.raw.fa"):
            if line.startswith(">"):
                line = ">" + line[1:].split("@")[0] + "\n"
            dst.write(line)
    os.remove(f"{WT}/bridge/target.raw.fa")
    return out


def build_windows(tranche):
    """windows.fa = each territory's CHM13#0 bridge row (locus.bed)."""
    out = f"{WT}/bridge/windows.fa"
    os.makedirs(f"{WT}/bridge", exist_ok=True)
    sample_of = load_map()
    specs, names = [], []
    with open(f"{WT}/bridge/window.specs", "w") as f:
        for row in tranche:
            T = row["tid"]
            a = open(f"{WT}/prep/{T}/locus.bed").read().split()
            path, s, e = a[0], int(a[1]), int(a[2])
            f.write(f"{path}@{sample_of[path]}:{s+1}-{e}\n")
            specs.append((path, s, e))
            names.append(T)
    subprocess.run(
        f"xargs -a {WT}/bridge/window.specs {AGC_BIN} getctg {AGC} "
        f"> {WT}/bridge/windows.raw.fa",
        shell=True, check=True)
    recs, buf = [], []
    for line in open(f"{WT}/bridge/windows.raw.fa"):
        if line.startswith(">"):
            if buf:
                recs.append("".join(buf))
            buf = []
        else:
            buf.append(line.strip())
    if buf:
        recs.append("".join(buf))
    assert len(recs) == len(names), (len(recs), len(names))
    with open(out, "w") as f:
        for T, (path, s, e), seq in zip(names, specs, recs):
            assert len(seq) == e - s, (T, len(seq), e - s)
            f.write(f">{T} {path}:{s}-{e}\n")
            for i in range(0, len(seq), 80):
                f.write(seq[i:i + 80] + "\n")
    os.remove(f"{WT}/bridge/windows.raw.fa")
    return out


def main():
    tranche = []
    with open(f"{WT}/tranche.tsv") as f:
        hdr = f.readline().rstrip("\n").split("\t")
        for line in f:
            if line.startswith("#"):
                continue
            tranche.append(dict(zip(hdr, line.rstrip("\n").split("\t"))))
    target = extract_target()
    windows = build_windows(tranche)
    paf = f"{WT}/bridge/truth2-bridge.paf"
    if not (os.path.exists(paf) and os.path.getsize(paf) > 0):
        subprocess.run(["python3", BRIDGE, windows, target, paf, str(NPROC)],
                       check=True)
    # the primary rows by window
    rows = {}
    for line in open(paf):
        a = line.rstrip("\n").split("\t")
        rows.setdefault(a[0], []).append(
            (a[4], a[5], int(a[7]), int(a[8])))
    out = open(f"{WT}/truth2.tsv", "w")
    out.write("tid\ttruth1\ttruth2\tclass\n")
    classes = {}
    for row in tranche:
        T = row["tid"]
        axis = row["axis_path"]
        truth_rows = []
        for line in open(f"{WT}/prep/{T}/locus.refined.bed"):
            a = line.rstrip("\n").split("\t")
            if len(a) >= 3 and a[0].startswith("HG002#2#"):
                truth_rows.append((a[0], int(a[1]), int(a[2])))
        prim = [r for r in rows.get(T, []) if r[0] == "+" or True]
        if not prim:
            cls = "no-bridge"
            pick = ""
        else:
            # max-overlap primary same-contig row (the split's own rule)
            prim.sort(key=lambda r: -(r[3] - r[2]))
            strand, contig, tst, ten = prim[0]
            overlapping = [(p, s, e) for (p, s, e) in truth_rows
                           if p.split(":")[0] == contig
                           and min(e, ten) - max(s, tst) > 0]
            if not overlapping:
                cls, pick = "no-truth2-row", ""
            else:
                overlaps = {r: min(r[2], ten) - max(r[1], tst)
                            for r in overlapping}
                best = max(overlaps.values())
                ties = [r for r, o in overlaps.items() if o == best]
                if len(ties) > 1:
                    cls, pick = "ambiguous-truth2", ""
                else:
                    p, s, e = ties[0]
                    cls, pick = "ok", f"{p}:{s}-{e}"
        classes[cls] = classes.get(cls, 0) + 1
        out.write(f"{T}\t{axis}\t{pick}\t{cls}\n")
    out.close()
    print("truth2 bridge classes:", classes)


if __name__ == "__main__":
    main()
