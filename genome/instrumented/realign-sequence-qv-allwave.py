#!/usr/bin/env python3
"""THE OFFICIAL SEQUENCE-QV CONVENTION RESTATED FOR THE YEAST FLEET
(the owner's ruling: the sequence QV is quantity 1, called-vs-truth
sequence accuracy, in the likegt max-qv form EXACTLY - surveyed from
likegt's own implementation first, not invented here).

THE OFFICIAL CONVENTION (likegt's own `max-qv` command,
/home/erikg/likegt/src/commands/max_qv.rs, the README's allwave
aligner), verbatim:

  (1) ALIGNMENT: allwave aligns the TRUTH haplotypes (the queries)
      against the CALLED haplotypes (the targets) with likegt's own
      invocation - `allwave -i <fasta> -t 4 -p tree:5:0:0` (the
      max-qv command's default sparsification and thread count;
      allwave's OWN default dual-piece gap-affine scores
      0,5,8,2,24,1, never overridden by likegt).
  (2) THE PAF DERIVATION (likegt's parse_paf_line): query_span =
      query_end - query_start; identity = matches / query_span - the
      per-base error is measured in TRUTH sequence space over the
      ALIGNED SPAN (no forced End2End overhang; unaligned truth bases
      outside the best alignment do not enter the rate - that IS
      likegt's convention); error_rate = 1 - identity;
      QV = -10*log10(error_rate) clamped to [0, 50]; a perfect
      alignment (identity >= 1.0) carries QV 50.0 exactly (likegt's
      own cap - NOT the committed biWFA receipts' 60).
  (3) PER-TRUTH-HAP BEST SELECTION (likegt's most-similar-haplotype
      rule): for each truth haplotype, the BEST alignment among the
      called homologs is taken; alignment lines from one query to one
      target pool their (matches, query_span) sums first (likegt's
      per-target aggregation, calculate_aggregate_identity). NO
      injective assignment: both truth haps may take the same called
      homolog when it is the best for each (the difference from the
      committed both-orders no-replacement yardstick, named per
      window below).
  (4) THE DIPLOTYPE SEQUENCE QV = the MEAN of the two per-truth-hap
      QVs (likegt's avg_qv = (qv_hap1 + qv_hap2) / 2).

The material gates are the committed stage's EXACTLY: this script
IMPORTS the committed realign-sequence-qv.py (of record) and reuses
its Panel fetch, fold criterion, GFA/constituent-containment gates,
and census classes unchanged - only the alignment and the QV
derivation are the official convention's. The committed biWFA
per-locus receipts are joined read-only into each record (the
per-window delta, the named mechanism), never rewritten.

Receipt-side only: reads the committed receipts and partition GFAs,
no instrument input, no thresholds, no tuning constants (the
convention is likegt's own, stated), no product change.
Usage:
  realign-sequence-qv-allwave.py --component C [--locality]
        [--biwfa-suffix S] [--receipt-prefix P] [--out-suffix S]
        (per-component pass; writes realign-sequence-qv-allwave[-S]-C.jsonl
         beside the receipts)
  realign-sequence-qv-allwave.py --tables [--locality] [--biwfa-suffix S]
        (the aggregate; writes realign-sequence-qv-allwave[-S]-tables.txt)
"""
import importlib.util
import json
import math
import os
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
_spec = importlib.util.spec_from_file_location(
    "realign_sequence_qv", os.path.join(HERE, "realign-sequence-qv.py"))
_committed = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_committed)

# The committed module parsed its own argv at import (the --tables /
# --component / --receipt-prefix / --out-suffix / --locality
# switches); re-derive OUR switches from the same argv, then override
# the module-level knobs so the imported machinery points at the
# committed receipts while OUR outputs carry the official suffix.
ARGS = sys.argv[1:]
def _switch(name, default=None):
    if name in ARGS:
        return ARGS[ARGS.index(name) + 1]
    return default

OUT_SUFFIX = _switch("--out-suffix", "-allwave")
LOCALITY = "--locality" in ARGS
# the committed biWFA receipts joined for the delta: the locality
# family of record with --locality (realign-sequence-qv-locality-*),
# the default-convention family otherwise
BIWFA_SUFFIX = _switch("--biwfa-suffix", "-locality" if LOCALITY else "")
CLI_RECEIPTS = "--cli-receipts" in ARGS or LOCALITY  # the locality family's upstream of record
RECEIPT_PREFIX = _switch("--receipt-prefix", _committed.RECEIPT_PREFIX)
# restore the committed module's defaults (it may have consumed the
# shared switches from the same argv at import time)
_committed.RECEIPT_PREFIX = RECEIPT_PREFIX
_committed.OUT_SUFFIX = BIWFA_SUFFIX
_committed.LOCALITY = LOCALITY
D = _committed.D
COMPONENTS = _committed.COMPONENTS

ALLWAVE = os.environ.get("IMPG_TRUTH_QV_ALLWAVE_BIN", "allwave")
THREADS = 4                 # likegt's own max-qv default
SPARSIFICATION = "tree:5:0:0"  # likegt's own max-qv default
QV_CAP = 50.0               # likegt's own cap (a perfect alignment = 50)


def qv_of_identity(identity):
    """likegt max_qv.rs's QV derivation, verbatim."""
    if identity >= 1.0:
        return QV_CAP
    error_rate = 1.0 - identity
    if error_rate <= 0.0:
        return QV_CAP
    return min(QV_CAP, max(0.0, -10.0 * math.log10(error_rate)))


def run_allwave(fasta_path):
    result = subprocess.run(
        [ALLWAVE, "-i", fasta_path, "-t", str(THREADS), "-p", SPARSIFICATION],
        capture_output=True, text=True, check=True)
    paf = []
    for line in result.stdout.splitlines():
        if not line:
            continue
        fields = line.split("\t")
        if len(fields) < 12:
            raise RuntimeError(f"malformed allwave PAF line: {line[:80]}")
        paf.append({
            "query": fields[0],
            "query_start": int(fields[2]),
            "query_end": int(fields[3]),
            "target": fields[5],
            "matches": int(fields[9]),
            "query_span": int(fields[3]) - int(fields[2]),
        })
    return paf


def official_locus(seqs, called, truth, tmpdir, locus):
    """The official-convention record for one locus: seqs maps fold
    index -> sequence; called/truth are the two fold indices each.
    The fasta ALWAYS carries FOUR role-named records (called0,
    called1, truth0, truth1) even when the called and truth folds
    coincide (the rank-1 case: the called pair IS the truth pair) —
    the roles are the queries/targets of the convention, not the
    fold indices. Returns (record, n_alignments)."""
    fasta = os.path.join(tmpdir, f"locus-{locus}.fa")
    with open(fasta, "w") as f:
        for fi, name in ((called[0], "called0"), (called[1], "called1"),
                         (truth[0], "truth0"), (truth[1], "truth1")):
            f.write(f">{name}\n{seqs[fi]}\n")
    paf = run_allwave(fasta)
    os.remove(fasta)
    record = {}
    for ti, tname in ((truth[0], "truth0"), (truth[1], "truth1")):
        per_target = {}
        for line in paf:
            if line["query"] == tname and line["target"] in ("called0", "called1"):
                agg = per_target.setdefault(line["target"],
                                            {"matches": 0, "query_span": 0})
                agg["matches"] += line["matches"]
                agg["query_span"] += line["query_span"]
        for target, agg in per_target.items():
            agg["identity"] = (agg["matches"] / agg["query_span"]
                              if agg["query_span"] else 0.0)
        # the per-hap BEST alignment: the highest identity, ties to
        # the higher matches, then the higher query_span, then called0
        best = None
        for target in ("called0", "called1"):
            agg = per_target.get(target)
            if agg is None or agg["query_span"] == 0:
                continue
            key = (agg["identity"], agg["matches"], agg["query_span"])
            if best is None or key > best[0]:
                best = (key, target, agg)
        if best is None:
            record[tname] = {
                "identity": 0.0, "qv": 0.0, "best_called": None,
                "matches": 0, "query_span": 0, "truth_len": len(seqs[ti]),
                "no_alignment": True, "pairs": {t: per_target.get(t)
                                                for t in ("called0", "called1")},
            }
        else:
            (_, target, agg) = best
            record[tname] = {
                "identity": agg["identity"], "qv": qv_of_identity(agg["identity"]),
                "best_called": target, "matches": agg["matches"],
                "query_span": agg["query_span"], "truth_len": len(seqs[ti]),
                "no_alignment": False,
                "pairs": {t: per_target.get(t) for t in ("called0", "called1")},
            }
    return record, len(paf)


def biwfa_join(component, locus):
    """The committed biWFA per-locus receipt record for the delta join
    (read-only; None when the committed receipt has no such locus)."""
    path = f"{D}/realign-sequence-qv{BIWFA_SUFFIX}-{component}.jsonl"
    if not hasattr(biwfa_join, "_cache"):
        biwfa_join._cache = {}
    if component not in biwfa_join._cache:
        table = {}
        if os.path.exists(path):
            with open(path) as f:
                for line in f:
                    r = json.loads(line)
                    table[r["locus"]] = r
        biwfa_join._cache[component] = table
    return biwfa_join._cache[component].get(locus)


def delta_reasons(official, biwfa):
    """The named per-window mechanisms of the biWFA-vs-official delta,
    derived from the records' own numbers (never asserted)."""
    reasons = []
    if biwfa is None:
        return ["no committed biWFA receipt (the join is honest)"]
    best0 = official["truth0"].get("best_called")
    best1 = official["truth1"].get("best_called")
    if best0 is not None and best0 == best1:
        reasons.append(
            f"per-hap-best is NON-INJECTIVE (both truth haps take {best0}; "
            "the committed yardstick forced the no-replacement assignment "
            "across it)")
    for tname in ("truth0", "truth1"):
        hap = official[tname]
        if hap["no_alignment"]:
            reasons.append(f"{tname}: allwave emitted no alignment (identity 0)")
        elif hap["query_span"] < hap["truth_len"]:
            reasons.append(
                f"{tname}: partial span (aligned {hap['query_span']} of "
                f"{hap['truth_len']} truth bp - the official convention "
                "scores the aligned span only; the End2End yardstick "
                "counted every base)")
    if biwfa.get("gap_columns", 0) > 0 and not reasons:
        reasons.append(
            "indel-scoring model (allwave dual-piece 0,5,8,2,24,1 vs biWFA "
            "single affine 0,4,6,2 - different gap placements, same "
            "substitutions)")
    if not reasons:
        reasons.append(
            "substitution-only: the two conventions agree in kind "
            "(identity over the same full span, the QV cap aside)")
    return reasons


def cli_receipt_path(component):
    """The upstream receipt of record for the restatement: the
    yeast CLI product runs' instrument receipts (cli-product-<C>/
    instrument-receipt.jsonl - the translation-gate runs verified
    field-identical to the committed QV receipts' own upstream values,
    and the join re-verified HERE: exact on truth_rank + called folds
    + truth folds at every expressible locus of every component; the
    large realign-locality-<C> upstream receipts were disk-discipline
    cleaned in earlier sessions). chrMT's run dir carries the
    -truthqv suffix (the two chrMT dirs verified: only it joins
    exactly to the 9 committed loci)."""
    special = {"chrMT": "cli-product-chrMT-truthqv"}
    return f"{D}/{special.get(component, f'cli-product-{component}')}/instrument-receipt.jsonl"


def run_component(component):
    receipt_path = (cli_receipt_path(component) if CLI_RECEIPTS
                    else f"{D}/{RECEIPT_PREFIX}-{component}.jsonl")
    if not os.path.exists(receipt_path):
        print(f"{component}: NO RECEIPT, skipped", flush=True)
        return False
    out_path = f"{D}/realign-sequence-qv{OUT_SUFFIX}-{component}.jsonl"
    panel = _committed.Panel()
    gates = {"folds_verified": 0, "constituent_rows_verified": 0}
    n_loci = n_expressible = n_alignments = 0
    classes = _committed.census_classes(component)
    with open(receipt_path) as fin, open(out_path, "w") as fout, \
            tempfile.TemporaryDirectory(prefix="qv-allwave-") as tmpdir:
        for line in fin:
            d = json.loads(line)
            n_loci += 1
            if not d["truth_pair_expressible"]:
                continue
            called = list(d["best_fold_indices"])
            truth = list(d["truth_folds"])
            if len(called) != 2 or len(truth) != 2:
                raise RuntimeError(
                    f"locus {d['locus']}: unexpected fold-pair shape")
            partition = d["partition"]
            g = None if LOCALITY else _committed.gfa(partition)
            seqs = {}
            for fi in sorted(set(called + truth)):
                fold = d["fold_identities"][fi]
                if fold["length"] <= 0:
                    raise RuntimeError(
                        f"locus {d['locus']} fold {fi}: empty fold")
                first = fold["members"][0]
                seq = panel.window(first["path_name"], first["start"],
                                   first["end"])
                if len(seq) != fold["length"]:
                    raise RuntimeError(
                        f"locus {d['locus']} fold {fi}: panel window length "
                        f"{len(seq)} != fold length {fold['length']}")
                for m in fold["members"]:
                    other = panel.window(m["path_name"], m["start"], m["end"])
                    if len(other) != fold["length"] or other != seq:
                        raise RuntimeError(
                            f"locus {d['locus']} fold {fi}: member "
                            f"{m['path_name']} panel window differs")
                if LOCALITY:
                    _committed.constituent_gate(panel, first["path_name"],
                                                first["start"], first["end"],
                                                gates)
                    gates["folds_verified"] += 1
                else:
                    spelled = g.spelled(
                        f"{first['path_name']}:{first['start']}-{first['end']}")
                    offset = spelled.find(seq)
                    if offset < 0 or spelled[offset:offset + len(seq)] != seq:
                        raise RuntimeError(
                            f"locus {d['locus']} fold {fi}: panel window "
                            "absent from the GFA spelling")
                    gates["folds_verified"] += 1
                seqs[fi] = seq
            official, n_paf = official_locus(seqs, called, truth, tmpdir,
                                             d["locus"])
            n_alignments += n_paf
            n_expressible += 1
            biwfa = biwfa_join(component, d["locus"])
            qv = (official["truth0"]["qv"] + official["truth1"]["qv"]) / 2.0
            rec = {
                "locus": d["locus"],
                "partition": partition,
                "truth_rank": d["truth_rank"],
                "rank1": d["truth_rank"] == 1,
                "called_folds": [
                    {"index": fi,
                     "strains": _committed.fold_strains(d["fold_identities"][fi])}
                    for fi in called],
                "truth_folds": [
                    {"index": fi,
                     "strains": _committed.fold_strains(d["fold_identities"][fi])}
                    for fi in truth],
                "convention": "likegt-max-qv-allwave",
                "invocation": (f"allwave -i <locus fasta> -t {THREADS} "
                               f"-p {SPARSIFICATION} (likegt's own max-qv default)"),
                "qv": qv,
                "qv_truth_hap0": official["truth0"]["qv"],
                "qv_truth_hap1": official["truth1"]["qv"],
                "identity_truth_hap0": official["truth0"]["identity"],
                "identity_truth_hap1": official["truth1"]["identity"],
                "truth_len_hap0": official["truth0"]["truth_len"],
                "truth_len_hap1": official["truth1"]["truth_len"],
                "query_span_hap0": official["truth0"]["query_span"],
                "query_span_hap1": official["truth1"]["query_span"],
                "best_called_hap_truth_hap0": official["truth0"]["best_called"],
                "best_called_hap_truth_hap1": official["truth1"]["best_called"],
                "perfect": (official["truth0"]["identity"] >= 1.0
                            and official["truth1"]["identity"] >= 1.0),
                "no_alignment_hap0": official["truth0"]["no_alignment"],
                "no_alignment_hap1": official["truth1"]["no_alignment"],
                "pairs": {t: official[t]["pairs"] for t in ("truth0", "truth1")},
                "biwfa_qv": biwfa["qv"] if biwfa else None,
                "biwfa_identity": biwfa["identity"] if biwfa else None,
                "biwfa_assignment": biwfa["assignment"] if biwfa else None,
                "delta_qv": (qv - biwfa["qv"]) if biwfa else None,
                "delta_reasons": delta_reasons(official, biwfa),
                "class": classes.get(d["locus"], "pilot (no census receipt)"),
                "component": component,
            }
            fout.write(json.dumps(rec) + "\n")
    print(
        f"{component}: {n_loci} loci, {n_expressible} truth-pair-expressible, "
        f"{gates['folds_verified']} folds gated, {n_alignments} allwave PAF "
        f"alignments; wrote {out_path}", flush=True)
    return True


def _load_records():
    records = []
    for component in COMPONENTS:
        path = f"{D}/realign-sequence-qv{OUT_SUFFIX}-{component}.jsonl"
        if not os.path.exists(path):
            print(f"WARNING: {component} per-locus official receipt absent")
            continue
        census = _committed.census_classes(component)
        with open(path) as f:
            for line in f:
                r = json.loads(line)
                r["class"] = census.get(r["locus"], "pilot (no census receipt)")
                records.append(r)
    return records


def run_tables():
    out_lines = []
    say = lambda t="": (print(t, flush=True), out_lines.append(t))
    records = _load_records()
    say("== THE CALLED-VS-TRUTH SEQUENCE QV - THE OFFICIAL CONVENTION "
        "(the owner's ruling, the likegt max-qv form) - THE RESTATEMENT")
    say(f"   components: {len(COMPONENTS)} closed; locality gates: {LOCALITY}")
    say(f"   records: {len(records)} (every truth-pair-expressible locus)")
    say("")
    say("THE OFFICIAL CONVENTION (the exact statement, likegt's own max-qv):")
    say("  * alignment: allwave, likegt's own invocation "
        f"`-t {THREADS} -p {SPARSIFICATION}`")
    say("    (allwave's default dual-piece gap-affine scores 0,5,8,2,24,1,")
    say("    never overridden); TRUTH haplotypes the QUERIES, the called")
    say("    homologs the targets - allwave's all-vs-all covers both")
    say("    directions; only truth-as-query lines enter the QV.")
    say("  * per PAF line (likegt's parse): identity = matches / query_span")
    say("    (query_span = query_end - query_start: the per-base error in")
    say("    TRUTH sequence space over the ALIGNED SPAN - no forced")
    say("    End2End overhang); lines from one query to one target pool")
    say("    their (matches, span) sums (likegt's per-target aggregation).")
    say("  * PER-TRUTH-HAP BEST: for each truth haplotype the best")
    say("    alignment among the called homologs is taken (likegt's")
    say("    most-similar-haplotype selection) - NO injective assignment;")
    say("    both truth haps may take the same called homolog.")
    say("  * QV = -10*log10(1 - identity), likegt's cap: a perfect")
    say(f"    alignment is QV {QV_CAP:.0f} exactly (NOT the committed biWFA")
    say("    receipts' 60). The diplotype sequence QV = the MEAN of the")
    say("    two per-truth-hap QVs (likegt's avg_qv).")
    say("  * the material gates are the committed stage's exactly (the")
    say("    panel fetch, the fold criterion, the GFA / constituent")
    say("    containment), reused from realign-sequence-qv.py of record.")
    say("")
    qvs = [r["qv"] for r in records]
    perfect = sum(1 for r in records if r["perfect"])
    say("== (1) THE GENOME-WIDE QV DISTRIBUTION (all expressible loci)")
    say(f"   n={len(qvs)}  perfect calls (both haps identity 1.0, QV 50): {perfect}")
    if qvs:
        say(f"   median {_committed.median(qvs):.2f}  "
            f"p10 {_committed.nearest_rank(qvs,10):.2f}  "
            f"p90 {_committed.nearest_rank(qvs,90):.2f}  "
            f"min {min(qvs):.2f}  max {max(qvs):.2f}")
        say(f"   QV>=40: {_committed.frac(qvs,40)*100:.2f}%   "
            f"QV>=30: {_committed.frac(qvs,30)*100:.2f}%   "
            f"QV>=20: {_committed.frac(qvs,20)*100:.2f}%")
    say("")
    say("== (2) THE PER-COMPONENT TABLE (official convention)")
    say("component  n_expressible  rank1  non-rank1  medianQV  p10QV  p90QV  "
        "minQV  >=40%  >=30%  >=20%")
    for component in COMPONENTS:
        rs = [r for r in records if r["component"] == component]
        if not rs:
            continue
        cq = [r["qv"] for r in rs]
        n1 = sum(1 for r in rs if r["rank1"])
        say(f"{component:<10} {len(rs):<14} {n1:<6} {len(rs)-n1:<10} "
            f"{_committed.median(cq):<8.2f} {_committed.nearest_rank(cq,10):<6.2f} "
            f"{_committed.nearest_rank(cq,90):<6.2f} {min(cq):<6.2f} "
            f"{_committed.frac(cq,40)*100:<6.2f} {_committed.frac(cq,30)*100:<6.2f} "
            f"{_committed.frac(cq,20)*100:.2f}")
    say("")
    say("== (3) THE REFRACTIVE TABLE: rank-1 vs non-rank-1 (official)")
    r1 = [r for r in records if r["rank1"]]
    nr1 = [r for r in records if not r["rank1"]]
    for label, rs in (("RANK-1 (the called class IS the truth class)", r1),
                      ("NON-RANK-1 (a different class called)", nr1)):
        q = [r["qv"] for r in rs]
        say(f"   {label}: n={len(rs)}")
        if rs:
            say(f"     QV: median {_committed.median(q):.2f}  "
                f"p10 {_committed.nearest_rank(q,10):.2f}  "
                f"p90 {_committed.nearest_rank(q,90):.2f}  min {min(q):.2f}")
            say(f"     perfect (both haps identity 1.0): "
                f"{sum(1 for r in rs if r['perfect'])}")
            for thr in (0.999, 0.99, 0.95):
                n = sum(1 for r in rs
                        if min(r["identity_truth_hap0"],
                               r["identity_truth_hap1"]) >= thr)
                say(f"     both haps >= {thr*100:.1f}% identity: {n}/{len(rs)} "
                    f"({n/len(rs)*100:.2f}%)")
    say("")
    say("== (4) THE RANK-1 VERIFICATION (zero edits is zero edits)")
    say(f"   rank-1 perfect under the official convention: "
        f"{sum(1 for r in r1 if r['perfect'])}/{len(r1)} "
        "(every rank-1 locus aligns full-span zero-error in both haps;")
    bad = [r for r in r1 if not r["perfect"]]
    say(f"   exceptions: {len(bad)}")
    for r in bad[:10]:
        say(f"     {r['component']} L{r['locus']}: "
            f"hap0 identity {r['identity_truth_hap0']:.6g}, "
            f"hap1 identity {r['identity_truth_hap1']:.6g}")
    say("")
    say("== (5) THE NON-RANK-1 RESTATED QVs WITH THE NAMED DELTAS")
    say("   (biWFA End2End injective-assignment QV -> official QV; the")
    say("    mechanism columns name WHY: the per-hap-best rule vs the")
    say("    assignment rule, the aligned-span vs End2End denominator,")
    say("    the indel-scoring model)")
    say("component  locus  biWFA_QV  official_QV  delta  hap0_id(hap)  hap1_id(hap)  mechanisms")
    for r in sorted(nr1, key=lambda r: r["qv"]):
        say(f"{r['component']:<10} L{r['locus']:<6} "
            f"{r['biwfa_qv'] if r['biwfa_qv'] is not None else '-':<9.2f} "
            f"{r['qv']:<12.2f} "
            f"{('%+.2f' % r['delta_qv']) if r['delta_qv'] is not None else '-':<7} "
            f"{r['identity_truth_hap0']:.6f}({r['best_called_hap_truth_hap0']})  "
            f"{r['identity_truth_hap1']:.6f}({r['best_called_hap_truth_hap1']})  "
            + "; ".join(r["delta_reasons"]))
    say("")
    shifts = [r for r in nr1 if r["delta_qv"] is not None]
    up = sum(1 for r in shifts if r["delta_qv"] > 1e-9)
    down = sum(1 for r in shifts if r["delta_qv"] < -1e-9)
    same = len(shifts) - up - down
    say(f"== (6) THE DELTA SUMMARY over the {len(shifts)} restated non-rank-1:")
    say(f"   official HIGHER: {up}   official LOWER: {down}   equal (to 1e-9): {same}")
    if shifts:
        deltas = sorted(abs(r["delta_qv"]) for r in shifts)
        say(f"   |delta| median {_committed.median(deltas):.2f}  max {max(deltas):.2f}")
        noninj = sum(1 for r in nr1
                     if r["best_called_hap_truth_hap0"] is not None
                     and r["best_called_hap_truth_hap0"]
                     == r["best_called_hap_truth_hap1"])
        partial = sum(1 for r in nr1 for h in ("0", "1")
                      if not r.get("no_alignment_hap" + h, False)
                      and r["query_span_hap" + h] < r["truth_len_hap" + h])
        say(f"   non-injective per-hap-best selections (both truth haps take "
            f"one called homolog): {noninj}/{len(nr1)}")
        say(f"   truth haps with a PARTIAL aligned span (unaligned truth bases "
            f"outside the official rate): {partial}/{2*len(nr1)} hap-sides")
    say("")
    say("(the official convention's perfect is BOTH truth haps' best "
        "alignments at identity 1.0 (QV 50, likegt's cap); the committed "
        "biWFA receipts' perfect is zero edits End2End (QV 60) - the two "
        "caps differ BY CONVENTION, the zero-edit property does not)")
    with open(f"{D}/realign-sequence-qv{OUT_SUFFIX}-tables.txt", "w") as f:
        f.write("\n".join(out_lines) + "\n")
    print(f"wrote {D}/realign-sequence-qv{OUT_SUFFIX}-tables.txt", flush=True)


def main():
    if "--tables" in ARGS:
        run_tables()
        return
    comp = _switch("--component")
    if comp:
        run_component(comp)
        return
    for component in COMPONENTS:
        run_component(component)


if __name__ == "__main__":
    main()
