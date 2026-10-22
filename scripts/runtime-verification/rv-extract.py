#!/usr/bin/env python3
"""THE WHOLE-GENOME RUNTIME VERIFICATION - the territory locus.fa extraction
(the committed extract-territory-locus-fa.py translated to the runtime run's
prep/<TID> layout, VERBATIM conventions): one auditable AGC getctg per member
via ONE xargs invocation per territory; the members = the refined bed's rows
plus the CHM13#0 bridge query row.
usage: rv-extract.py <RV> <TID>
"""
import subprocess, sys, os

RV = sys.argv[1]
TID = sys.argv[2]
AGC = "/home/erikg/hprcv2/HPRC_r2_assemblies_0.6.1.agc"
AGC_BIN = "/home/erikg/agc/bin/agc"
MAPFILE = "/home/erikg/hprcv2/generic-sampling/agc-path-sample-map.txt"

def load_map():
    m, cur = {}, None
    for line in open(MAPFILE):
        if line.startswith(" "):
            m.setdefault(line.strip().split()[0], cur)
        else:
            cur = line.strip().split()[0]
    return m

def main():
    Ldir = f"{RV}/prep/{TID}"
    sample_of = load_map()
    members = []
    with open(f"{Ldir}/locus.refined.bed") as f:
        for line in f:
            a = line.rstrip("\n").split("\t")
            if len(a) >= 4:
                members.append((a[0], int(a[1]), int(a[2]), a[3]))
    with open(f"{Ldir}/locus.bed") as f:
        a = f.read().split()
    members.append((a[0], int(a[1]), int(a[2]), "+"))
    specs, wanted = [], []
    for path, s, e, strand in members:
        sample = sample_of[path]
        specs.append(f"{path}@{sample}:{s+1}-{e}")     # agc: 1-based inclusive
        wanted.append((path, s, e, strand))           # bed: 0-based half-open
    with open(f"{Ldir}/extract.specs.txt", "w") as f:
        f.write("\n".join(specs) + "\n")
    subprocess.run(f"xargs -a {Ldir}/extract.specs.txt {AGC_BIN} getctg {AGC} > {Ldir}/extract.raw.fa",
                   shell=True, check=True, timeout=1800)
    recs, buf = [], []
    for line in open(f"{Ldir}/extract.raw.fa"):
        if line.startswith(">"):
            if buf:
                recs.append("".join(buf))
            buf = []
        elif line.strip():
            buf.append(line.strip())
    if buf:
        recs.append("".join(buf))
    if len(recs) != len(wanted):
        sys.exit(f"EXTRACTION COUNT MISMATCH: {len(recs)} records vs {len(wanted)} members")
    total_bp = 0
    with open(f"{Ldir}/locus.fa", "w") as fa, open(f"{Ldir}/locus.members.bed", "w") as mb:
        for (path, s, e, strand), seq in zip(wanted, recs):
            assert len(seq) == e - s, f"{path}:{s}-{e}: len {len(seq)} != bed {e - s}"
            total_bp += len(seq)
            fa.write(f">{path}:{s}-{e}\n")
            for i in range(0, len(seq), 80):
                fa.write(seq[i:i+80] + "\n")
            mb.write(f"{path}\t{s}\t{e}\t{strand}\n")
    os.remove(f"{Ldir}/extract.raw.fa")
    print(f"{TID}: {len(wanted)} members, {total_bp} bp")

if __name__ == "__main__":
    main()
