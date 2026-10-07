#!/usr/bin/env python3
"""THE CLI-PRODUCT SAMPLE AGGREGATE (the owner's product swap, stage 4's
before/after receipt): the full diploid validation sample's numbers as
measured END-TO-END THROUGH THE CLI (the 17 components' calls.jsonl
truth-QV test-mode artifacts - the assessment side lives only there),
compared against the committed assessment numbers of record
(genome-truth-rank-tables-locality.txt: rank-1 880/1,132 expressible =
77.7%, the per-component table).

Receipt-side only; no thresholds; no instrument inputs.
Usage: cli-product-aggregate.py    (writes cli-product-gate-AGGREGATE.txt
at the validation dir)
"""
import json
import os

D = "/home/erikg/yeast/genome-balanced-diploid-validation-20260930"
COMPONENTS = [
    "chrMT", "chrI", "chrII", "chrIII", "chrIV", "chrV", "chrVI",
    "chrVII", "chrVIII", "chrIX", "chrX", "chrXI", "chrXII", "chrXIII",
    "chrXIV", "chrXV", "chrXVI",
]
# chrMT/chrI landed in the stage-1/2 session dirs; the fleet components
# at cli-product-<C>.
DIRS = {
    "chrMT": f"{D}/cli-product-chrMT-truthqv",
    "chrI": f"{D}/cli-product-chrI",
}

out_lines = []
def say(text=""):
    print(text, flush=True)
    out_lines.append(text)

rows = []
total_loci = total_expressible = total_rank1 = total_perfect = 0
qv_errors = []
for c in COMPONENTS:
    base = DIRS.get(c, f"{D}/cli-product-{c}")
    path = f"{base}/calls.jsonl.truth-qv.jsonl"
    if not os.path.exists(path):
        say(f"{c}: NO CLI PRODUCT RUN, pending")
        rows.append((c, None, None, None, None))
        continue
    n = e = r1 = perfect = 0
    for line in open(path):
        t = json.loads(line)
        n += 1
        if t.get("truth_pair_expressible"):
            e += 1
            if t.get("rank1"):
                r1 += 1
                if t.get("perfect"):
                    perfect += 1
                if t.get("qv") is not None:
                    qv_errors.append((t.get("error", 0.0), t.get("qv")))
    rows.append((c, n, e, r1, perfect))
    total_loci += n; total_expressible += e
    total_rank1 += r1; total_perfect += perfect

say("== THE CLI-PRODUCT SAMPLE AGGREGATE (17 components, end-to-end through the CLI)")
say("component  loci  expressible  truth-rank1  perfect")
for (c, n, e, r1, p) in rows:
    say(f"{c:10s} {str(n):>5s}  {str(e):>11s}  {str(r1):>11s}  {str(p):>7s}"
        if n is not None else f"{c:10s} {'PENDING':>5s}")
say()
say(f"THE AGGREGATE: loci {total_loci}, expressible {total_expressible}, "
    f"truth rank-1 {total_rank1}, perfect {total_perfect}")
if total_expressible:
    say(f"  rank-1 rate: {total_rank1}/{total_expressible} = "
        f"{100.0 * total_rank1 / total_expressible:.1f}%")
if qv_errors:
    errors = sorted(x[0] for x in qv_errors)
    qvs = sorted(x[1] for x in qv_errors)
    median = qvs[len(qvs) // 2] if len(qvs) % 2 == 1 else (
        qvs[len(qvs) // 2 - 1] + qvs[len(qvs) // 2]) / 2.0
    p10 = qvs[max(0, int(0.10 * len(qvs)) - 1)]
    say(f"  QV: perfect {total_perfect} (= the rank-1 count by the committed "
        f"measurement), median {median:.2f}, p10 {p10:.2f}, "
        f"QV>=40 at {100.0 * sum(1 for q in qvs if q >= 40) / len(qvs):.2f}%")
say()
say("== THE COMMITTED ASSESSMENT NUMBERS OF RECORD (the locality fleet: "
    "rank-1 880/1,132 expressible = 77.7%, QV perfect 880, median 60.00, "
    "p10 14.42, QV>=40 78.18%)")
say("the before/after receipt: the CLI end-to-end numbers vs the committed "
    "table above; any divergence named per component by the per-component "
    "gate files (cli-product-gate-<C>.txt).")

with open(f"{D}/cli-product-gate-AGGREGATE.txt", "w") as f:
    f.write("\n".join(out_lines) + "\n")
