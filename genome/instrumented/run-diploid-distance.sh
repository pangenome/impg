#!/usr/bin/env bash
# Re-score the saved diploid fleet serially; never starts a product run.
set -euo pipefail
ROOT=/home/erikg/impg-genome-inference
D=/home/erikg/yeast/genome-diploid-validation-20260929
mapfile -t components < <(python3 -B - "$D/private-truth/truth-diploid.json" <<'PY'
import json,sys
for route in json.load(open(sys.argv[1]))['slot1']:
    print(route['component'])
PY
)
for component in "${components[@]}"; do
    tag=fleet
    if [[ "$component" == chrI || "$component" == chrMT ]]; then tag=smoke2; fi
    prefix="$D/distance-$tag-$component"
    rm -f "$prefix.done" "$prefix.exit"
    start=$(date +%s)
    if python3 -B "$ROOT/genome/instrumented/score-diploid-distance.py" "$component" "$tag" \
        > "$prefix.json.tmp" 2> "$prefix.err"; then
        mv "$prefix.json.tmp" "$prefix.json"
        echo 0 > "$prefix.exit"
        echo $(( $(date +%s) - start )) > "$prefix.wall"
        touch "$prefix.done"
        echo "$component: score exit 0 / $(cat "$prefix.wall")s"
    else
        status=$?
        echo "$status" > "$prefix.exit"
        echo $(( $(date +%s) - start )) > "$prefix.wall"
        rm -f "$prefix.json.tmp"
        touch "$prefix.done"
        echo "$component: score exit $status; stopped" >&2
        exit "$status"
    fi
done
python3 -B - "$D" "${components[@]}" <<'PY'
import json,sys
from pathlib import Path
D=Path(sys.argv[1]); components=sys.argv[2:]
rows=[]
for component in components:
    tag='smoke2' if component in ('chrI','chrMT') else 'fleet'
    x=json.loads((D/f'distance-{tag}-{component}.json').read_text())
    rows.append({k:x[k] for k in ('component','distance','error_bp','truth_bp','correct_truth_bp',
                                  'covered_truth_bp','uncovered_truth_bp','inserted_query_bp',
                                  'unaligned_query_bp','duplicated_truth_bp',
                                  'candidate_alignment_blocks','selected_alignment_blocks',
                                  'selected_ploidy','truth_self_error_bp','truth_swapped_error_bp')})
total_error=sum(x['error_bp'] for x in rows)
total_truth=sum(x['truth_bp'] for x in rows)
summary={'distance':total_error/total_truth,'error_bp':total_error,'truth_bp':total_truth,
         'components':rows}
(D/'distance-genome.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Genome aligned diploid distance:', summary['distance'],
      'errors:',total_error,'truth_bp:',total_truth)
PY
