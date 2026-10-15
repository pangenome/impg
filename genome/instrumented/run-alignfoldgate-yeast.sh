#!/usr/bin/env bash
# STAGE 1 GATE 4: THE YEAST IDENTITY GATE. (a) env-unset reruns of the
# committed chrMT/chrI exhaustive locality receipts through the rebuilt
# binary - byte-identical (all four sidecars + the main receipt); (b) the
# relaxed (env-set) reruns for the spot-check (the folds that should
# stay exact stay exact; no genuinely-distinguishable classes merged).
set -uo pipefail
cd /home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-balanced-diploid-validation-20260930

loci_of () {
  python3 - "$1" << 'EOF'
import json, sys
c = sys.argv[1]
axis = json.load(open(
    "/home/erikg/yeast/genome-mem-bwt-bootstrap-20260911T055203Z/reference-axis-v1.json"
))["intervals"]
n = sum(1 for iv in axis if iv["component"] == f"S288C#0#{c}")
print(",".join(str(i) for i in range(n)))
EOF
}

for C in chrMT chrI; do
  L=$(loci_of "$C")
  # (a) the env-unset identity run
  if [ ! -f "$D/run-alignfoldgate-unset-$C.done" ]; then
    bash genome/instrumented/run-realign-locality.sh 4 "$C" "$L" \
      "realign-locality-alignfoldgate-unset-$C" "alignfoldgate-unset-$C" \
      > "$D/alignfoldgate-unset-$C.log" 2>&1
    echo "$?" > "$D/run-alignfoldgate-unset-$C.exit"
    touch "$D/run-alignfoldgate-unset-$C.done"
  fi
  # (b) the relaxed rerun
  if [ ! -f "$D/run-alignfoldgate-relaxed-$C.done" ]; then
    IMPG_REALIGN_ALIGNMENT_CONNECTED_FOLDS=1 \
    bash genome/instrumented/run-realign-locality.sh 4 "$C" "$L" \
      "realign-locality-alignfoldgate-relaxed-$C" "alignfoldgate-relaxed-$C" \
      > "$D/alignfoldgate-relaxed-$C.log" 2>&1
    echo "$?" > "$D/run-alignfoldgate-relaxed-$C.exit"
    touch "$D/run-alignfoldgate-relaxed-$C.done"
  fi
done

# the byte-identity diffs (env-unset vs the committed receipts of
# record): all four SIDECARS byte-identical; the main receipt
# field-identical on every semantic field with the walls/rss_kb
# exemption (the committed identity-gate convention of slices D/E - a
# rerun's timing fields are its own).
: > "$D/alignfoldgate-identity.txt"
python3 - "$D" << 'EOF' >> "$D/alignfoldgate-identity.txt"
import filecmp, json, sys
D = sys.argv[1]
for C in ("chrMT", "chrI"):
    base_c = f"{D}/realign-locality-exhaustive-{C}"
    base_a = f"{D}/realign-locality-alignfoldgate-unset-{C}"
    for suf in (".exactness.jsonl", ".skeleton.jsonl",
                ".jsonl.records.jsonl", ".jsonl.ingredients.jsonl"):
        same = filecmp.cmp(f"{base_c}{suf}", f"{base_a}{suf}", shallow=False)
        print(f"{'IDENTICAL' if same else 'DIFFERS'} {C}{suf}")
    def load(p):
        out = {}
        for line in open(p):
            d = json.loads(line)
            out[d["locus"]] = d
        return out
    b, a = load(f"{base_c}.jsonl"), load(f"{base_a}.jsonl")
    bad = [
        (l, k)
        for l in b
        for k in b[l]
        if k not in ("walls", "rss_kb") and b[l].get(k) != a[l].get(k)
    ]
    print(
        f"{'IDENTICAL' if not bad else 'DIFFERS'} {C}.jsonl "
        f"(semantic fields, walls/rss_kb exempt)"
        + (f" VIOLATIONS {bad[:3]}" if bad else "")
    )
EOF
echo "=== alignfold yeast gate chain complete $(date -u +%FT%TZ)"
