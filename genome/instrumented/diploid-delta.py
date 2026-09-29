#!/usr/bin/env python3
"""Assessment-only haploid-vs-diploid component table (different denominator labels).

Only completed diploid records are shown; never substitute one-slot accuracy
for missing 2x2 accuracy or a bracket for a decision error.
"""
import json
from pathlib import Path

base = Path('/home/erikg/impg/target/experiments/genome-mem-bwt-pipeline/genome-genotyping-v1-item1')
dip = Path('/home/erikg/yeast/genome-diploid-validation-20260929')
truth = json.loads((dip / 'private-truth/truth-diploid.json').read_text())
for row in truth['slot1']:
    comp = row['component']
    hp = base / f'score-rerun-S288C#0#{comp}.json'
    dp = dip / f'score-{"smoke2" if comp in ("chrMT", "chrI") else "fleet"}-{comp}.json'
    if not dp.exists():
        continue
    h = json.loads(hp.read_text())['rows']['selected']
    d = json.loads(dp.read_text())
    if d['exit'] != '0':
        continue
    print(f"{comp:7} H_acc_H={h['agg_haploid_acc']:.4f}/{h['windows']} "
          f"H_switch={h['switches']}/{h['boundaries']} "
          f"D_ploidy={d['selected_ploidy']} "
          f"D_exact={d['exact_homolog_windows']}/{d['windows']} "
          f"bracket={d['bracketed_homolog_windows']} "
          f"D_truth_domain={d['truth_pair_in_domain_windows']}/{d['windows']} "
          f"(slots {d['truth_in_domain_by_slot']}) "
          f"D_pair_vs_selected_m1_nats={d['truth_pair_vs_selected_m1_nats']:.2f} "
          f"D_2x2={d['pair_accuracy_attested_only']} ({d['two_slot_windows']} pairs) "
          f"D_candidate={d['candidate_accuracy_attested_only']} "
          f"({d['candidate_only_windows']} single-slot) "
          f"D_switch={d['switches_per_haplotype']}/{d['switch_boundaries']} "
          f"D_dose={d['selected_single_class_dosage']} / true {d['truth_dosage']} "
          f"walls_H={(base / f'run-rerun-S288C#0#{comp}.wall').read_text().strip()}s,D={d['wall_seconds']}s")
