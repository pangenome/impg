#!/usr/bin/env python3
"""THE AXIS BRIDGE (the committed bridge recipe, second leg): every grid
window (query, named by TID) aligned to the sample's OWN declared-axis
assembly HG002#1#chr1-22 (the target fasta extracted from the panel AGC),
preset asm5, the SAME acceptance rule as the grid's CHM13 bridge - the axis
interval of a window is its UNIQUE primary same-chromosome row (no-bridge /
ambiguous-bridge the honest fail-closed classes). The committed binary
(minimap2 2.31) is no longer on the box; mappy 2.31 is the same minimap2
core at the same version - the deviation (the python binding, not the CLI
binary) is stated in the datasheet. Output: PAF-shaped rows
(qname qlen qst qen strand tname tlen tst ten nmatch alen mapq tp:A:P|S).
usage: rv-axis-bridge.py <windows.fa> <axis-target.fa> <out.paf> [nproc]
"""
import sys, os
import mappy

windows, target, out, nproc = sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4] or 24)

def rows_for(recs):
    a = ALIGNER
    out_rows = []
    for name, seq in recs:
        for h in a.map(seq, cs=False):
            if not h.m_primary:
                continue
            out_rows.append(
                f"{name}\t{len(seq)}\t{h.q_st}\t{h.q_en}\t{'+' if h.strand == 1 else '-'}\t"
                f"{h.ctg}\t{h.ctg_len}\t{h.r_st}\t{h.r_en}\t{h.mlen}\t{h.blen}\t{h.mapq}\ttp:A:P\n")
    return out_rows

def main():
    global ALIGNER
    ALIGNER = mappy.Aligner(target, preset="asm5")
    if not ALIGNER:
        sys.exit("failed to load the axis target index")
    # window records in grid order; fork-pool shares the index (COW)
    recs = []
    name, buf = None, []
    for line in open(windows):
        if line.startswith(">"):
            if name is not None:
                recs.append((name, "".join(buf)))
            name, buf = line[1:].strip().split()[0], []
        else:
            buf.append(line.strip())
    if name is not None:
        recs.append((name, "".join(buf)))
    print(f"{len(recs)} windows; aligning with {nproc} workers", flush=True)
    import multiprocessing as mp
    with open(out, "w") as f:
        with mp.Pool(nproc) as pool:
            for rows in pool.imap(rows_for, [recs[i::nproc] for i in range(nproc)]):
                for r in rows:
                    f.write(r)
    print("done", flush=True)

if __name__ == "__main__":
    ALIGNER = None
    main()
