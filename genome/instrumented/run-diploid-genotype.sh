#!/usr/bin/env bash
# Assessment-only: score the existing saved per-locus genotype pairs serially.
set -euo pipefail
ROOT=/home/erikg/impg-genome-inference
D="${IMPG_DIPLOID_VALIDATION_DIR:-/home/erikg/yeast/genome-diploid-validation-20260929}"
mapfile -t components < <(python3 -B - "$D/private-truth/truth-diploid.json" <<'PY'
import json,sys
for route in json.load(open(sys.argv[1]))['slot1']:
    print(route['component'])
PY
)
for component in "${components[@]}"; do
    tag=fleet
    if [[ "$component" == chrI || "$component" == chrMT ]]; then tag=smoke2; fi
    prefix="$D/genotype-distance-$tag-$component"
    rm -f "$prefix.done" "$prefix.exit"
    start=$(date +%s)
    if python3 -B "$ROOT/genome/instrumented/score-diploid-genotype.py" "$component" "$tag" \
        > "$prefix.json.tmp" 2> "$prefix.err"; then
        mv "$prefix.json.tmp" "$prefix.json"
        echo 0 > "$prefix.exit"
        echo $(( $(date +%s) - start )) > "$prefix.wall"
        touch "$prefix.done"
        echo "$component: attested genotype score exit 0 / $(cat "$prefix.wall")s"
    else
        status=$?
        echo "$status" > "$prefix.exit"
        echo $(( $(date +%s) - start )) > "$prefix.wall"
        rm -f "$prefix.json.tmp"
        touch "$prefix.done"
        echo "$component: attested genotype score exit $status; stopped" >&2
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
    x=json.loads((D/f'genotype-distance-{tag}-{component}.json').read_text())
    rows.append({k:x[k] for k in ('component','error_bp','truth_bp','distance','assessable_loci',
                                  'bracketed_loci','bracketed_S288C_axis_bp',
                                  'assignment_transitions','informative_adjacent_boundaries',
                                  'assignment_ties','assignment_track','truth_dosage',
                                  'truth_pair_in_domain_loci',
                                  'exact_truth_pair_when_expressible_loci',
                                  'balanced_called_dosage_class_loci',
                                  'selected_single_class_dosage','selected_ploidy')})
error=sum(row['error_bp'] for row in rows)
truth=sum(row['truth_bp'] for row in rows)
summary={'yardstick':'injective_per_locus_genotype_on_attested_SK1_homologs',
         'distance':error/truth, 'error_bp':error,'truth_bp':truth,
         'assessable_loci':sum(row['assessable_loci'] for row in rows),
         'bracketed_loci':sum(row['bracketed_loci'] for row in rows),
         'assignment_transitions':sum(row['assignment_transitions'] for row in rows),
         'informative_adjacent_boundaries':sum(row['informative_adjacent_boundaries'] for row in rows),
         'assignment_ties':sum(row['assignment_ties'] for row in rows),
         'truth_pair_in_domain_loci':sum(row['truth_pair_in_domain_loci'] for row in rows),
         'exact_truth_pair_when_expressible_loci':sum(row['exact_truth_pair_when_expressible_loci'] for row in rows),
         'balanced_called_dosage_class_loci':sum(row['balanced_called_dosage_class_loci'] for row in rows),
         'bracket_policy':'no score or zero-length SK1 allele at untyped loci',
         'components':rows}
(D/'genotype-distance-genome.json').write_text(json.dumps(summary,indent=2)+'\n')
print('Genome assessable per-locus genotype distance:',summary['distance'],
      'errors:',error,'attested_truth_bp:',truth,
      'assessable_loci:',summary['assessable_loci'],'bracketed_loci:',summary['bracketed_loci'])
PY
