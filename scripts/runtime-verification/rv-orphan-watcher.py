#!/usr/bin/env python3
"""THE ORPHAN WATCHER: the last 8 orphaned prep-step processes (their bash
wrappers are dead - earlier waves were killed without group-kill) are doing
REAL work (the segdup-tail partitions: 3,000-5,000 regions per territory).
On each one's exit:
  - a partition whose beds.err carries the impg-written terminal line
    ('Partitioned into N regions: ... written / ... sequence'): write
    beds.done + beds.wall measured from /proc's own start time - the
    compute is kept, the marker-gated re-prep continues at the export step.
  - a partition that died WITHOUT the terminal line: rm -rf its partial
    partition dir (stale beds would poison the re-prep's region counts).
  - the syng orphan (its success is unverifiable - the wrapper that would
    capture the exit status is dead): rm -f its partial syng outputs; the
    marker-gated re-prep redoes the step cleanly.
Detached; exits when all tracked pids are gone; touches orphan-watcher.done.
usage: rv-orphan-watcher.py <RV> <pid:kind:TID>...
"""
import os, re, sys, time

RV = sys.argv[1]
CLK = os.sysconf("SC_CLK_TCK")
btime = 0
for line in open("/proc/stat"):
    if line.startswith("btime "):
        btime = int(line.split()[1])
TERM = re.compile(rb"Partitioned into \d+ regions:.*written /")

def start_epoch(pid):
    with open(f"/proc/{pid}/stat") as f:
        parts = f.read().rsplit(") ", 1)[1].split()
    return btime + int(parts[19]) / CLK

def terminal_in(err):
    try:
        with open(err, "rb") as f:
            f.seek(-8192, 2)
            return bool(TERM.search(f.read()))
    except OSError:
        try:
            with open(err, "rb") as f:
                return bool(TERM.search(f.read()))
        except OSError:
            return False

watch = []
for spec in sys.argv[2:]:
    pid, kind, tid = spec.split(":")
    watch.append({"pid": int(pid), "kind": kind, "tid": tid,
                  "start": start_epoch(pid)})
print(f"watching {len(watch)} orphan prep steps: "
      + " ".join(f"{w['tid']}({w['kind']})" for w in watch), flush=True)

while watch:
    for w in list(watch):
        if os.path.exists(f"/proc/{w['pid']}"):
            continue
        d = f"{RV}/prep/{w['tid']}"
        if w["kind"] == "partition":
            if terminal_in(f"{d}/beds.err"):
                open(f"{d}/beds.done", "w").close()
                wall = max(1, int(time.time() - w["start"]))
                with open(f"{d}/beds.wall", "w") as f:
                    f.write(str(wall))
                with open(f"{d}/beds.wall-estimated", "w") as f:
                    f.write(f"orphan-completed; wall from /proc start of pid {w['pid']}\n")
                print(f"{w['tid']}: partition COMPLETE (wall {wall}s) - marked", flush=True)
            else:
                os.system(f"rm -rf {d}/partition-w10k-d1k")
                print(f"{w['tid']}: partition died without terminal line - "
                      "partial cleaned for re-prep", flush=True)
        else:  # syng: exit status unverifiable
            os.system(f"rm -f {d}/locus-pos64.syng {d}/locus-pos64.syng.*")
            print(f"{w['tid']}: syng orphan exited (status lost) - partial "
                  "outputs cleaned for a clean re-run", flush=True)
        watch.remove(w)
    if watch:
        time.sleep(30)
open(f"{RV}/orphan-watcher.done", "w").close()
print("orphan watcher complete", flush=True)
