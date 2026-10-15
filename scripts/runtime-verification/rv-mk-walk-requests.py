#!/usr/bin/env python3
# THE C-RIDER (measurement only): the walk-request generator. Reads the
# committed partition BEDs (every row of every partition) + the locality
# syng's .names (name -> path_id) and emits the `syng_homology_dump
# walks` request table: path_id, start, end, tag (= the full row name
# PATH:SRCSTART-SRCEND:LOCSTART-LOCEND, the same names the export GFAs'
# P lines and the locality maps use).
import sys
beds, names_path, out_path = sys.argv[1:4]
name2id = {}
with open(names_path) as f:
    for line in f:
        a = line.rstrip("\n").split("\t")
        if len(a) >= 2:
            name2id[a[1]] = int(a[0])
n = 0
missing = 0
with open(out_path, "w") as out:
    import os
    for fn in sorted(os.listdir(beds)):
        if not fn.endswith(".bed"):
            continue
        with open(os.path.join(beds, fn)) as f:
            for line in f:
                if not line.strip():
                    continue
                a = line.rstrip("\n").split("\t")
                row_name, s, e = a[0], int(a[1]), int(a[2])
                pid = name2id.get(row_name)
                if pid is None:
                    missing += 1
                    continue
                out.write(f"{pid}\t{s}\t{e}\t{row_name}:{s}-{e}\n")
                n += 1
print(f"requests: {n} rows, missing names: {missing}")
