#!/usr/bin/env python3
"""THE ORPHAN/CORPSE SWEEP (one-shot): earlier prep waves were killed without
group-kill - their impg children finished or died with NO wrapper alive, so
beds.done markers were never written. For every TID with a partition-w10k-d1k
dir but no beds.done:
  - the partition COMPLETED (the impg-written terminal line 'Partitioned into
    N regions: ... written / ... sequence' in beds.err, >=1 .bed in results):
    write beds.done honestly (the state genuinely exists; the marker-gated
    re-prep then continues at the export step) + beds.wall ESTIMATED from the
    filesystem timestamps (noted - the wrapper's clock is gone). No re-run.
  - else (killed mid-run): rm -rf the partition dir (stale partial beds would
    poison the re-prep's region counts) - the marker-gated re-prep redoes it.
Live (still-running) orphan partitions are skipped here - the watcher handles
them on completion.
usage: rv-orphan-sweep.py <RV>
"""
import os, re, sys, time

RV = sys.argv[1]
PREP = f"{RV}/prep"

# live partition/syng orphans (from ps at sweep time)
live = set()
out = os.popen(
    "ps -eo args | grep 'impg partition' | grep -v grep | "
    "sed -n 's#.*runtime-verification/prep/\\(T[0-9]*\\)/.*#\\1#p'").read()
live |= {l.strip() for l in out.splitlines() if l.strip()}
out = os.popen(
    "ps -eo args | grep 'impg syng' | grep -v grep | "
    "sed -n 's#.*runtime-verification/prep/\\(T[0-9]*\\)/.*#\\1#p'").read()
live |= {l.strip() for l in out.splitlines() if l.strip()}

TERM = re.compile(r"Partitioned into \d+ regions:.*written /")
marked, cleaned, kept_live = [], [], []
for T in sorted(os.listdir(PREP)):
    d = f"{PREP}/{T}"
    pdir = f"{d}/partition-w10k-d1k"
    if not (os.path.isdir(pdir) and not os.path.exists(f"{d}/beds.done")):
        continue
    err = f"{d}/beds.err"
    if T in live:
        kept_live.append(T)
        continue
    term = False
    if os.path.exists(err):
        with open(err, errors="replace") as f:
            for line in f:
                if TERM.search(line):
                    term = True
    if term and os.path.isdir(f"{pdir}/results") \
            and any(x.endswith(".bed") for x in os.listdir(f"{pdir}/results")):
        open(f"{d}/beds.done", "w").close()
        os.utime(f"{d}/beds.done", None)
        # the wall estimate: beds.err mtime (last log line ~ completion)
        # minus the results dir's birth (first partition file ~ start)
        end = os.path.getmtime(err)
        start = os.path.getmtime(f"{pdir}/starting-sequences.txt")
        if start > end:
            start = end - 1
        with open(f"{d}/beds.wall", "w") as f:
            f.write(str(max(1, int(end - start))))
        with open(f"{d}/beds.wall-estimated", "w") as f:
            f.write("orphan-completed; wall estimated from filesystem timestamps\n")
        marked.append(T)
    else:
        os.system(f"rm -rf {pdir}")
        cleaned.append(T)
print(f"sweep: marked-complete {len(marked)} cleaned-partial {len(cleaned)} "
      f"kept-live {len(kept_live)}")
if kept_live:
    print("live:", " ".join(kept_live))
