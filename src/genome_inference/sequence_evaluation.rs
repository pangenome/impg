//! Standalone base-level PAF replay. Whole FASTAs remain the denominators.
use super::{atomic_write, read_json, reconstruction::fingerprint, write_json, FORMAT_VERSION};
use crate::sample_mem_bwt::invalid;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::BTreeMap,
    fs::{self, File},
    io::{self, BufRead, BufReader},
    path::Path,
    process::{Command, Stdio},
};
pub const MODEL: &str = "greedy-one-to-one-sequence-evaluation-v1";
pub type Fasta = BTreeMap<String, Vec<u8>>;

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b.to_ascii_uppercase() {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'R' => b'Y',
            b'Y' => b'R',
            b'S' => b'S',
            b'W' => b'W',
            b'K' => b'M',
            b'M' => b'K',
            b'B' => b'V',
            b'V' => b'B',
            b'D' => b'H',
            b'H' => b'D',
            _ => b'N',
        })
        .collect()
}
/// Empty FASTA is a valid no-call result; empty records, duplicate names and
/// non-IUPAC DNA are errors. No header normalization or truth trimming.
pub fn read_fasta(path: &Path) -> io::Result<Fasta> {
    let file = File::open(path)?;
    // niffler's magic detection requires five bytes; empty no-calls and a
    // one-base plain FASTA are both valid inputs shorter than that.
    let reader: Box<dyn io::Read> = if file.metadata()?.len() < 5 {
        Box::new(file)
    } else {
        niffler::get_reader(Box::new(file))
            .map_err(io::Error::other)?
            .0
    };
    let mut result = Fasta::new();
    let mut name = None;
    let mut sequence = Vec::new();
    fn store(result: &mut Fasta, name: Option<String>, seq: &mut Vec<u8>) -> io::Result<()> {
        if let Some(name) = name {
            if seq.is_empty() || result.insert(name, std::mem::take(seq)).is_some() {
                return Err(invalid("empty FASTA record or duplicate name"));
            }
        } else if !seq.is_empty() {
            return Err(invalid("sequence before FASTA header"));
        }
        Ok(())
    }
    for line in BufReader::new(reader).lines() {
        let line = line?;
        if let Some(header) = line.strip_prefix('>') {
            store(&mut result, name.take(), &mut sequence)?;
            name = Some(
                header
                    .split_whitespace()
                    .next()
                    .filter(|s| !s.is_empty())
                    .ok_or_else(|| invalid("empty FASTA name"))?
                    .into(),
            );
        } else {
            for b in line.bytes() {
                let b = b.to_ascii_uppercase();
                if !b"ACGTRYSWKMBDHVN".contains(&b) {
                    return Err(invalid("invalid FASTA DNA"));
                }
                sequence.push(b);
            }
        }
    }
    store(&mut result, name, &mut sequence)?;
    Ok(result)
}
#[derive(Clone, Debug, Default, Serialize, Deserialize)]
pub struct Counts {
    pub matches: u64,
    pub mismatches: u64,
    pub insertion_bp: u64,
    pub deletion_bp: u64,
    pub unassessed_columns: u64,
}
impl Counts {
    fn assessed(&self) -> u64 {
        self.matches + self.mismatches + self.insertion_bp + self.deletion_bp
    }
    fn add(&mut self, c: &Self) {
        self.matches += c.matches;
        self.mismatches += c.mismatches;
        self.insertion_bp += c.insertion_bp;
        self.deletion_bp += c.deletion_bp;
        self.unassessed_columns += c.unassessed_columns;
    }
}
#[derive(Clone, Debug, Serialize)]
struct Alignment {
    line: usize,
    query: String,
    query_start: u64,
    query_end: u64,
    target: String,
    target_start: u64,
    target_end: u64,
    strand: String,
    mapq: u64,
    cigar: String,
    counts: Counts,
    reported_matches: u64,
    reported_columns: u64,
    soft_clipped_bp: u64,
    hard_clipped_bp: u64,
}
fn number(s: &str) -> io::Result<u64> {
    s.parse().map_err(|_| invalid("invalid PAF/CIGAR integer"))
}
fn cigar(s: &str) -> io::Result<Vec<(u64, u8)>> {
    let mut result = Vec::new();
    let mut start = 0;
    for (i, b) in s.bytes().enumerate() {
        if !b.is_ascii_digit() {
            let n = number(&s[start..i])?;
            if n == 0 || !b"M=XIDSH".contains(&b) {
                return Err(invalid(
                    "unsupported/zero-length CIGAR operation (only M=XID and terminal SH)",
                ));
            }
            result.push((n, b));
            start = i + 1;
        }
    }
    if start != s.len() || result.is_empty() {
        return Err(invalid("truncated/empty CIGAR"));
    }
    Ok(result)
}
fn acgt(b: u8) -> bool {
    b"ACGT".contains(&b)
}
fn parse(line: &str, id: usize, query: &Fasta, truth: &Fasta) -> io::Result<Alignment> {
    let f: Vec<_> = line.split('\t').collect();
    if f.len() < 12 {
        return Err(invalid("PAF needs at least 12 tab-separated fields"));
    }
    let q = query
        .get(f[0])
        .ok_or_else(|| invalid("unknown PAF query"))?;
    let t = truth
        .get(f[5])
        .ok_or_else(|| invalid("unknown PAF target"))?;
    let ql = number(f[1])?;
    let qs = number(f[2])?;
    let qe = number(f[3])?;
    let tl = number(f[6])?;
    let ts = number(f[7])?;
    let te = number(f[8])?;
    let reported_matches = number(f[9])?;
    let reported_columns = number(f[10])?;
    let mapq = number(f[11])?;
    if ql != q.len() as u64
        || tl != t.len() as u64
        || qs >= qe
        || qe > ql
        || ts >= te
        || te > tl
        || !matches!(f[4], "+" | "-")
        || mapq > 255
        || reported_matches > reported_columns
    {
        return Err(invalid("PAF length, bounds, strand or count mismatch"));
    }
    let tags: Vec<_> = f[12..]
        .iter()
        .filter_map(|s| s.strip_prefix("cg:Z:"))
        .collect();
    if tags.len() != 1 {
        return Err(invalid("exactly one base-level cg:Z CIGAR required"));
    }
    let ops = cigar(tags[0])?;
    let first = ops
        .iter()
        .position(|(_, b)| !b"SH".contains(b))
        .ok_or_else(|| invalid("clipping-only CIGAR"))?;
    let last = ops.iter().rposition(|(_, b)| !b"SH".contains(b)).unwrap();
    let sum = |slice: &[(u64, u8)]| -> io::Result<u64> {
        slice.iter().try_fold(0u64, |a, (n, _)| {
            a.checked_add(*n).ok_or_else(|| invalid("CIGAR overflow"))
        })
    };
    let lead = sum(&ops[..first])?;
    let tail = sum(&ops[last + 1..])?;
    if (lead > 0 || tail > 0)
        && (lead != if f[4] == "+" { qs } else { ql - qe }
            || tail != if f[4] == "+" { ql - qe } else { qs })
    {
        return Err(invalid(
            "terminal clipping disagrees with full query coordinates",
        ));
    }
    let oriented = if f[4] == "-" {
        reverse_complement(&q[qs as usize..qe as usize])
    } else {
        q[qs as usize..qe as usize].to_vec()
    };
    let target = &t[ts as usize..te as usize];
    let mut qi = 0usize;
    let mut ti = 0usize;
    let mut counts = Counts::default();
    for &(n, op) in &ops[first..=last] {
        if b"SH".contains(&op) {
            return Err(invalid("internal clipping is invalid"));
        }
        let n = usize::try_from(n).map_err(|_| invalid("oversized CIGAR"))?;
        let consumes_q = op != b'D';
        let consumes_t = op != b'I';
        if (consumes_q && n > oriented.len() - qi) || (consumes_t && n > target.len() - ti) {
            return Err(invalid("CIGAR overconsumes interval"));
        }
        for _ in 0..n {
            match op {
                b'I' => {
                    if acgt(oriented[qi]) {
                        counts.insertion_bp += 1;
                    } else {
                        counts.unassessed_columns += 1;
                    }
                    qi += 1;
                }
                b'D' => {
                    if acgt(target[ti]) {
                        counts.deletion_bp += 1;
                    } else {
                        counts.unassessed_columns += 1;
                    }
                    ti += 1;
                }
                _ => {
                    let a = oriented[qi];
                    let b = target[ti];
                    if acgt(a) && acgt(b) {
                        if (op == b'=' && a != b) || (op == b'X' && a == b) {
                            return Err(invalid("CIGAR =/X contradicts actual ACGT bases"));
                        }
                        if a == b {
                            counts.matches += 1;
                        } else {
                            counts.mismatches += 1;
                        }
                    } else {
                        counts.unassessed_columns += 1;
                    }
                    qi += 1;
                    ti += 1;
                }
            }
        }
    }
    if qi != oriented.len()
        || ti != target.len()
        || counts.assessed() + counts.unassessed_columns != reported_columns
    {
        return Err(invalid("CIGAR consumption/alignment-column count mismatch"));
    }
    Ok(Alignment {
        line: id,
        query: f[0].into(),
        query_start: qs,
        query_end: qe,
        target: f[5].into(),
        target_start: ts,
        target_end: te,
        strand: f[4].into(),
        mapq,
        cigar: tags[0].into(),
        counts,
        reported_matches,
        reported_columns,
        soft_clipped_bp: ops.iter().filter(|(_, b)| *b == b'S').map(|(n, _)| n).sum(),
        hard_clipped_bp: ops.iter().filter(|(_, b)| *b == b'H').map(|(n, _)| n).sum(),
    })
}
#[derive(Default)]
struct Coverage(BTreeMap<String, Vec<(u64, u64)>>);
impl Coverage {
    fn overlaps(&self, n: &str, s: u64, e: u64) -> bool {
        self.0
            .get(n)
            .is_some_and(|v| v.iter().any(|&(a, b)| s < b && e > a))
    }
    fn add(&mut self, n: &str, s: u64, e: u64) {
        self.0.entry(n.into()).or_default().push((s, e));
    }
    fn lengths(&self) -> BTreeMap<String, u64> {
        self.0
            .iter()
            .map(|(n, v)| {
                let mut v = v.clone();
                v.sort_unstable();
                let mut end = 0;
                let mut bp = 0;
                for (s, e) in v {
                    bp += e.saturating_sub(s.max(end));
                    end = end.max(e);
                }
                (n.clone(), bp)
            })
            .collect()
    }
    fn bp(&self) -> u64 {
        self.lengths().values().sum()
    }
}
pub fn evaluate(query: &Fasta, truth: &Fasta, paf: impl BufRead) -> io::Result<Value> {
    let mut alignments = Vec::new();
    let mut raw_q = Coverage::default();
    let mut raw_t = Coverage::default();
    for (i, line) in paf.lines().enumerate() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let a = parse(&line, i + 1, query, truth)
            .map_err(|e| invalid(format!("PAF line {}: {e}", i + 1)))?;
        raw_q.add(&a.query, a.query_start, a.query_end);
        raw_t.add(&a.target, a.target_start, a.target_end);
        alignments.push(a);
    }
    // Length before identity: no best-local-identity cherry picking. A conflicting
    // record is rejected in full (no invented truncated CIGAR or partial credit).
    alignments.sort_by(|a, b| {
        let span = |a: &Alignment| (a.query_end - a.query_start).min(a.target_end - a.target_start);
        span(b)
            .cmp(&span(a))
            .then(b.mapq.cmp(&a.mapq))
            .then(a.line.cmp(&b.line))
    });
    let mut selected_q = Coverage::default();
    let mut selected_t = Coverage::default();
    let mut counts = Counts::default();
    let mut selected = Vec::new();
    let mut rejected = Vec::new();
    for a in &alignments {
        let q = selected_q.overlaps(&a.query, a.query_start, a.query_end);
        let t = selected_t.overlaps(&a.target, a.target_start, a.target_end);
        if q || t {
            rejected.push(json!({"alignment":a,"query_overlap":q,"truth_overlap":t}));
        } else {
            selected_q.add(&a.query, a.query_start, a.query_end);
            selected_t.add(&a.target, a.target_start, a.target_end);
            counts.add(&a.counts);
            selected.push(a);
        }
    }
    let mut by_query: BTreeMap<&str, Vec<&Alignment>> = BTreeMap::new();
    for a in &selected {
        by_query.entry(&a.query).or_default().push(a);
    }
    let mut structural = Vec::new();
    for (name, mut rows) in by_query {
        rows.sort_by_key(|a| a.query_start);
        if rows.len() > 1 {
            let mut changes = Vec::new();
            for p in rows.windows(2) {
                let reason = if p[0].target != p[1].target {
                    "target-path-change"
                } else if p[0].strand != p[1].strand {
                    "orientation-change"
                } else if (p[0].strand == "+" && p[1].target_start < p[0].target_end)
                    || (p[0].strand == "-" && p[1].target_end > p[0].target_start)
                {
                    "noncollinear-order"
                } else {
                    "collinear-split"
                };
                changes.push(json!({"left_line":p[0].line,"right_line":p[1].line,"reason":reason,"unaligned_query_gap_bp":p[1].query_start-p[0].query_end}));
            }
            structural
                .push(json!({"query":name,"selected_alignments":rows.len(),"changes":changes}));
        }
    }
    let total = |f: &Fasta| f.values().map(|s| s.len() as u64).sum::<u64>();
    let qb = total(query);
    let tb = total(truth);
    let assessed = counts.assessed();
    let errors = counts.mismatches + counts.insertion_bp + counts.deletion_bp;
    let fraction = |n: u64, d: u64| {
        if d == 0 {
            None
        } else {
            Some(n as f64 / d as f64)
        }
    };
    let per_sequence = |f: &Fasta, raw: &Coverage, selected: &Coverage| {
        let r = raw.lengths();
        let s = selected.lengths();
        f.iter().map(|(n,b)|json!({"name":n,"total_bp":b.len(),"unknown_bp":b.iter().filter(|&&b|!acgt(b)).count(),
            "raw_aligned_union_bp":r.get(n).copied().unwrap_or(0),"selected_aligned_union_bp":s.get(n).copied().unwrap_or(0),"unaligned_bp":b.len() as u64-s.get(n).copied().unwrap_or(0)})).collect::<Vec<_>>()
    };
    Ok(json!({"version":FORMAT_VERSION,"model":MODEL,
        "scope":"Alignment-derived QV, NOT calibrated Merqury or assembly-wide QV. Whole FASTAs are denominators; zero errors on aligned columns is not perfect-genome recovery. Missing mappings may reflect aligner thresholds/repeats/short blocks, not automatically inference error.",
        "selection_policy":"Deterministic greedy whole-record one-to-one: descending min(query span, truth span), descending MAPQ, ascending PAF line. Reject entire record on ANY query or truth coordinate overlap. No identity-based sorting. Conservative, not globally optimal; alternative evidence retained.",
        "query_total_bp":qb,"truth_total_bp":tb,"query_records":query.len(),"truth_records":truth.len(),
        "query_unknown_bp":query.values().map(|s|s.iter().filter(|&&b|!acgt(b)).count() as u64).sum::<u64>(),
        "truth_unknown_bp":truth.values().map(|s|s.iter().filter(|&&b|!acgt(b)).count() as u64).sum::<u64>(),
        "selected_query_unknown_bp":selected.iter().map(|a|query[&a.query][a.query_start as usize..a.query_end as usize].iter().filter(|&&b|!acgt(b)).count() as u64).sum::<u64>(),
        "selected_truth_unknown_bp":selected.iter().map(|a|truth[&a.target][a.target_start as usize..a.target_end as usize].iter().filter(|&&b|!acgt(b)).count() as u64).sum::<u64>(),
        "raw_query_aligned_union_bp":raw_q.bp(),"raw_truth_aligned_union_bp":raw_t.bp(),
        "selected_query_aligned_union_bp":selected_q.bp(),"selected_truth_aligned_union_bp":selected_t.bp(),
        "unaligned_query_bp":qb-selected_q.bp(),"unaligned_truth_bp":tb-selected_t.bp(),
        "query_coverage":fraction(selected_q.bp(),qb),"truth_coverage":fraction(selected_t.bp(),tb),
        "counts":counts,"assessed_columns":assessed,"error_columns":errors,"identity":fraction(counts.matches,assessed),
        "alignment_qv":if assessed>0&&errors>0 {Some(-10.0*(errors as f64/assessed as f64).log10())} else {None},
        "alignment_qv_status":if assessed==0 {"unavailable-no-assessed-columns"} else if errors==0 {"zero-observed-errors"} else {"finite"},
        "raw_alignment_count":alignments.len(),"selected_alignment_count":selected.len(),"rejected_overlap_alignment_count":rejected.len(),
        "raw_query_span_bp_sum":alignments.iter().map(|a|a.query_end-a.query_start).sum::<u64>(),
        "raw_truth_span_bp_sum":alignments.iter().map(|a|a.target_end-a.target_start).sum::<u64>(),
        "selected_alignments":selected,"rejected_alternatives":rejected,"structural_diagnostics":structural,
        "query_sequences":per_sequence(query,&raw_q,&selected_q),"truth_sequences":per_sequence(truth,&raw_t,&selected_t)}))
}

pub fn run(
    query: &Path,
    truth: &Path,
    paf: &Path,
    out: &Path,
    aligner_metadata: Option<&Path>,
) -> io::Result<Value> {
    let mut evaluation = evaluate(
        &read_fasta(query)?,
        &read_fasta(truth)?,
        BufReader::new(File::open(paf)?),
    )?;
    evaluation["inputs"] =
        json!({"query":fingerprint(query)?,"truth":fingerprint(truth)?,"paf":fingerprint(paf)?});
    evaluation["aligner"] = if let Some(path) = aligner_metadata {
        read_json(path)?
    } else {
        json!({"identity":"not supplied; external PAF","parameters":"unknown; no aligner recall claim"})
    };
    // Keep raw evidence even when caller supplied a PAF outside this output directory.
    if paf != out.join("alignments.paf") {
        fs::copy(paf, out.join("alignments.paf"))?;
    }
    write_json(&out.join("evaluation.json"), &evaluation)?;
    Ok(
        json!({"query_coverage":evaluation["query_coverage"],"truth_coverage":evaluation["truth_coverage"],"alignment_qv":evaluation["alignment_qv"],"alignment_qv_status":evaluation["alignment_qv_status"]}),
    )
}

/// Fixed, explicit evaluation settings; no frequency masking, approximate mode,
/// truth-based tuning or internal one-to-one filtering. Raw alternatives survive.
pub fn align(
    query: &Path,
    truth: &Path,
    tool: &Path,
    out: &Path,
    threads: usize,
) -> io::Result<Value> {
    if !(1..=4).contains(&threads) {
        return Err(invalid("alignment threads must be 1..4"));
    }
    let tool = fs::canonicalize(tool)?;
    let version = Command::new(&tool).arg("--version").output()?;
    if !version.status.success() {
        return Err(invalid("aligner version command failed"));
    }
    let identity = format!(
        "{}{}",
        String::from_utf8_lossy(&version.stdout),
        String::from_utf8_lossy(&version.stderr)
    );
    // Pin the inspected option contract, not an arbitrary wfmash-compatible executable.
    if !identity.contains("b55cf75") {
        return Err(invalid("sequence alignment driver requires inspected wfmash b55cf75; other tools can supply external base-level PAF"));
    }
    let args = [
        "-t",
        &threads.to_string(),
        "-H",
        "0",
        "-f",
        "-n",
        "10",
        "-S",
        "10",
        "-s",
        "1000",
        "-l",
        "1000",
        "-p",
        "90",
        "-k",
        "19",
    ];
    let paf = out.join("alignments.paf");
    let log = out.join("aligner.stderr.log");
    // Native tools may create .fai indexes or temporary files. Isolate all such
    // writes in the reserved output, never beside frozen reconstruction/truth.
    let q = read_fasta(query)?;
    let t = read_fasta(truth)?;
    let materialize = |name: &str, records: &Fasta| -> io::Result<std::path::PathBuf> {
        let path = out.join(name);
        let mut bytes = Vec::new();
        for (n, s) in records {
            bytes.extend(format!(">{n}\n").bytes());
            for line in s.chunks(80) {
                bytes.extend(line);
                bytes.push(b'\n');
            }
        }
        atomic_write(&path, &bytes)?;
        fs::canonicalize(path)
    };
    let alignment_query = materialize("alignment-query.fa", &q)?;
    let alignment_truth = materialize("alignment-truth.fa", &t)?;
    let skipped = q.is_empty() || t.is_empty();
    let metadata = json!({"tool":fingerprint(&tool)?,"version":identity,"args":args,"target":alignment_truth,"query":alignment_query,
        "original_query":fingerprint(query)?,"original_truth":fingerprint(truth)?,"working_query":fingerprint(&alignment_query)?,"working_truth":fingerprint(&alignment_truth)?,
        "working_copy_policy":"Uncompressed uppercase IUPAC FASTA, unchanged names/lengths/bases and complete records; no interval trimming or coordinate-suffix normalization.",
        "working_directory":fs::canonicalize(out)?,"skipped_empty_input":skipped,
        "base_level":true,"frequency_masking":"H0; no top-frequency exclusion","mapping_filtering":"-f disables mapping filtering; n10/S10 explicit",
        "limitations":"k19, p90, segment1000/block1000; queries shorter than segment use S10. Very short/repetitive/low-identity blocks can remain unmapped. No recall guarantee. Defaults not listed are pinned by executable fingerprint/version."});
    write_json(&out.join("aligner.json"), &metadata)?;
    // wfmash does not accept an empty query. A no-call is still a valid evaluation.
    if skipped {
        atomic_write(&paf, b"")?;
        atomic_write(&log, b"Empty query or truth: aligner skipped.\n")?;
    } else {
        let status = Command::new(&tool)
            .arg(&alignment_truth)
            .arg(&alignment_query)
            .args(args)
            .current_dir(out)
            .stdout(Stdio::from(File::create(&paf)?))
            .stderr(Stdio::from(File::create(&log)?))
            .status()?;
        if !status.success() {
            return Err(invalid(format!(
                "wfmash failed: {status}; raw PAF/stderr retained"
            )));
        }
    }
    run(query, truth, &paf, out, Some(&out.join("aligner.json")))
}

#[cfg(test)]
mod tests {
    use super::*;
    fn fasta(name: &str, seq: &str) -> Fasta {
        BTreeMap::from([(name.into(), seq.as_bytes().to_vec())])
    }
    fn paf(
        q: &str,
        qs: u64,
        qe: u64,
        strand: &str,
        t: &str,
        ts: u64,
        te: u64,
        cigar: &str,
        cols: u64,
    ) -> String {
        format!(
            "q\t{}\t{qs}\t{qe}\t{strand}\tt\t{}\t{ts}\t{te}\t0\t{cols}\t60\tcg:Z:{cigar}\n",
            q.len(),
            t.len()
        )
    }
    #[test]
    fn exact_mismatch_indels_denominator_qv_and_missing_truth() {
        let q = "ACTGGT";
        let t = "ACCGCTAAAA";
        // AC match, T/C mismatch, G insertion, G match, C deletion, T match.
        let p = paf(q, 0, 6, "+", t, 0, 6, "2=1X1I1M1D1=", 7);
        let r = evaluate(&fasta("q", q), &fasta("t", t), p.as_bytes()).unwrap();
        assert_eq!(
            r["counts"],
            json!({"matches":4,"mismatches":1,"insertion_bp":1,"deletion_bp":1,"unassessed_columns":0})
        );
        assert_eq!(r["assessed_columns"], 7);
        assert_eq!(r["truth_coverage"], 0.6);
        assert_eq!(r["unaligned_truth_bp"], 4);
        assert!(
            (r["alignment_qv"].as_f64().unwrap() + 10.0 * (3.0f64 / 7.0).log10()).abs() < 1e-12
        );
        assert_eq!(r["identity"], json!(4.0 / 7.0));
    }
    #[test]
    fn reverse_clipping_unknowns_and_zero_errors_are_qualified() {
        let q = "TTACGNAA";
        let t = "NCGT";
        let p = paf(q, 2, 6, "-", t, 0, 4, "2S4M2S", 4);
        let r = evaluate(&fasta("q", q), &fasta("t", t), p.as_bytes()).unwrap();
        assert_eq!(r["counts"]["matches"], 3);
        assert_eq!(r["counts"]["unassessed_columns"], 1);
        assert_eq!(r["unaligned_query_bp"], 4);
        assert_eq!(r["alignment_qv_status"], "zero-observed-errors");
        assert!(r["alignment_qv"].is_null());
        assert_eq!(r["selected_alignments"][0]["soft_clipped_bp"], 4);
        let p = paf("NN", 0, 2, "+", "NN", 0, 2, "2M", 2);
        assert_eq!(
            evaluate(&fasta("q", "NN"), &fasta("t", "NN"), p.as_bytes()).unwrap()
                ["alignment_qv_status"],
            "unavailable-no-assessed-columns"
        );
        let p = paf("AN", 0, 2, "+", "NA", 0, 2, "1I1M1D", 3);
        let r = evaluate(&fasta("q", "AN"), &fasta("t", "NA"), p.as_bytes()).unwrap();
        assert_eq!(
            r["counts"],
            json!({"matches":0,"mismatches":0,"insertion_bp":1,"deletion_bp":1,"unassessed_columns":1})
        );
        assert_eq!(reverse_complement(b"ACGTRYSWKMBDHVN"), b"NBDHVKMWSRYACGT");
    }
    #[test]
    fn duplicate_reconstruction_and_alternative_mappings_never_double_credit() {
        let q = BTreeMap::from([
            ("q".into(), b"ACGT".to_vec()),
            ("copy".into(), b"ACGT".to_vec()),
        ]);
        let p = paf("ACGT", 0, 4, "+", "ACGT", 0, 4, "4M", 4);
        let p = format!("{p}{}{}", p.replace("q\t", "copy\t"), p);
        let r = evaluate(&q, &fasta("t", "ACGT"), p.as_bytes()).unwrap();
        assert_eq!(r["counts"]["matches"], 4);
        assert_eq!(r["selected_alignment_count"], 1);
        assert_eq!(r["rejected_overlap_alignment_count"], 2);
        assert_eq!(r["query_coverage"], 0.5);
        assert_eq!(r["raw_query_aligned_union_bp"], 8);
        assert_eq!(r["raw_truth_aligned_union_bp"], 4);
        assert_eq!(r["unaligned_query_bp"], 4);
    }
    #[test]
    fn split_noncollinear_orientation_and_chromosome_diagnostics() {
        let q = fasta("q", "AAAACCCCGGGGTTTT");
        let t = BTreeMap::from([
            ("t".into(), b"CCCCAAAACCCC".to_vec()),
            ("chr2".into(), b"TTTT".to_vec()),
        ]);
        let p = format!(
            "{}{}{}{}",
            paf("AAAACCCCGGGGTTTT", 0, 4, "+", "CCCCAAAACCCC", 4, 8, "4M", 4),
            paf("AAAACCCCGGGGTTTT", 4, 8, "+", "CCCCAAAACCCC", 0, 4, "4M", 4),
            paf(
                "AAAACCCCGGGGTTTT",
                8,
                12,
                "-",
                "CCCCAAAACCCC",
                8,
                12,
                "4M",
                4
            ),
            paf("AAAACCCCGGGGTTTT", 12, 16, "+", "TTTT", 0, 4, "4M", 4)
                .replace("\tt\t", "\tchr2\t")
        );
        let r = evaluate(&q, &t, p.as_bytes()).unwrap();
        let c = &r["structural_diagnostics"][0]["changes"];
        assert_eq!(c[0]["reason"], "noncollinear-order");
        assert_eq!(c[1]["reason"], "orientation-change");
        assert_eq!(c[2]["reason"], "target-path-change");
        assert_eq!(r["counts"]["matches"], 16);
    }
    #[test]
    fn concatenated_xz_fails_closed_instead_of_shortening_truth() {
        // Two valid independent XZ streams, generated with Python lzma.compress.
        // niffler's single-stream XZ decoder cannot consume their concatenation.
        let bytes = include_bytes!("../../tests/fixtures/genome-sequence-two-streams.fa.xz");
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("truth.fa.xz");
        fs::write(&path, &bytes[..76]).unwrap();
        assert_eq!(read_fasta(&path).unwrap(), fasta("truth1", "AAAAAAAAAA"));
        fs::write(&path, &bytes[76..]).unwrap();
        assert_eq!(read_fasta(&path).unwrap(), fasta("truth2", "CCCCCCCCCC"));
        fs::write(&path, bytes).unwrap();
        assert!(read_fasta(&path).is_err());
        // Both records belong to the whole-truth denominator, never just stream 1.
        fs::write(&path, b">truth1\nAAAAAAAAAA\n>truth2\nCCCCCCCCCC\n").unwrap();
        let truth = read_fasta(&path).unwrap();
        assert_eq!(truth.len(), 2);
        assert_eq!(truth.values().map(Vec::len).sum::<usize>(), 20);
    }

    #[test]
    fn no_resolved_output_zero_truth_coverage_and_corrupt_cigar_rejected() {
        let r = evaluate(&Fasta::new(), &fasta("t", "AAAA"), b"".as_slice()).unwrap();
        assert_eq!(r["truth_coverage"], 0.0);
        assert_eq!(r["unaligned_truth_bp"], 4);
        assert_eq!(r["assessed_columns"], 0);
        assert!(r["identity"].is_null());
        for cg in [
            "",
            "4",
            "0M",
            "3M",
            "5M",
            "4X",
            "1M1S2M",
            "4N",
            "18446744073709551616M",
            "1S4M",
        ] {
            let p = paf("AAAA", 0, 4, "+", "AAAA", 0, 4, cg, 4);
            assert!(
                evaluate(&fasta("q", "AAAA"), &fasta("t", "AAAA"), p.as_bytes()).is_err(),
                "{cg}"
            );
        }
        for p in [
            paf("AAAA", 0, 4, "?", "AAAA", 0, 4, "4M", 4),
            paf("AAAA", 0, 4, "+", "AAAA", 0, 4, "4M", 3),
            paf("AAAA", 0, 4, "+", "AAAA", 0, 4, "4M", 4).replace("cg:Z:", "zz:Z:"),
            paf("AAAA", 0, 4, "+", "AAAA", 0, 4, "4M", 4).replace("q\t4", "q\t5"),
            paf("AAAA", 0, 4, "+", "AAAA", 0, 4, "4M", 4).replace("q\t", "foreign\t"),
        ] {
            assert!(evaluate(&fasta("q", "AAAA"), &fasta("t", "AAAA"), p.as_bytes()).is_err());
        }
        let dir = tempfile::tempdir().unwrap();
        let file = dir.path().join("bad.fa");
        fs::write(&file, b"").unwrap();
        assert!(read_fasta(&file).unwrap().is_empty());
        fs::write(&file, b">q\nA").unwrap();
        assert_eq!(read_fasta(&file).unwrap(), fasta("q", "A"));
        for s in [">q\nA\n>q\nA\n", ">q\n", ">q\nZ\n", "A\n", ">\nA\n"] {
            fs::write(&file, s).unwrap();
            assert!(read_fasta(&file).is_err());
        }
    }
}
