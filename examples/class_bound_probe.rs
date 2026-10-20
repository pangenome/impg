//! The admissible-bound prototype probe (receipt-side, NO behavior
//! change anywhere in the product): the frontier-1 branch-and-bound
//! design's pruning-rate measurement at the committed worst-locus
//! shapes, computed OFFLINE from the tranche's committed call
//! receipts — no panel, no cut, no call re-run (docs/syng-gmem-bwt/
//! global-panel-call-mode.md SS10.6.14's autopsy arithmetic is the
//! substrate; the tranche's heavy-tail members are the shapes).
//!
//! THE BOUND BEING PRICED (the design, derived from the receipts —
//! every quantity a committed one, no tuned constants):
//!
//!   THE E-DECOMPOSITION. Every per-(fold, unit) log-likelihood is
//!   either exactly the elsewhere constant E (the no-pin branch — the
//!   committed row-compression exactness argument: an entry equals the
//!   no-pin constant iff the unit is unpinned on the fold) or E plus a
//!   positive excess (a pinned branch; logsumexp2(local, E) > E
//!   strictly). The class pair sum therefore decomposes exactly (real
//!   arithmetic):
//!
//!     S(f,g) = sum_u c_u * mix(ll_f[u], ll_g[u])
//!            = E*C + sum_u c_u * mix0(alpha_f[u], alpha_g[u])
//!
//!   where mix is the committed 1/2-1/2 homolog mixture, mix0(a,b) =
//!   ln(0.5e^a + 0.5e^b), and C = sum_u c_u.
//!
//!   THE PER-FOLD EXCESS-MASS BOUND (the E-per-component additivity
//!   the campaign receipted, SS10.6.14's log_gap_decomposition): for
//!   alpha, beta >= 0, mix0(alpha,beta) <= alpha + beta (proof: WLOG
//!   alpha >= beta, then mix0 = alpha + ln(0.5 + 0.5e^{beta - alpha})
//!   <= alpha <= alpha + beta). Therefore
//!
//!     S(f,g) <= E*C + A_f + A_g,   A_f = sum_u c_u * alpha_f[u]
//!
//!   — ONE SCALAR per fold, computable in one walk over the fold's
//!   PINNED entries (the sparse row holds exactly the alpha > 0
//!   entries): the whole bound index costs one class-pair-walk
//!   equivalent per locus, folds/2 fewer than the exhaustive scan.
//!
//!   THE DERIVED SKIP THRESHOLDS (each from an emitted consumer's own
//!   arithmetic, never tuned):
//!     - the called set / best:   skip iff E*C + A_f + A_g < best
//!       (strict — the pair cannot win or tie the emitted argmax);
//!     - the QUAL posterior terms: exp(ll - best) underflows to
//!       exactly +0.0 iff ll - best <= the f64 exp underflow point
//!       (derived by evaluation below: the largest f64 t with
//!       exp(t) == 0.0); a skipped pair's every QUAL contribution is
//!       then BIT-IDENTICAL to the exact run's (+0.0 additions are
//!       float identity on the non-negative running sums);
//!     - the truth-rank counts (the --truth-qv assessment artifact
//!       emits 1 + #{classes strictly above the truth pair}): skip iff
//!       E*C + A_f + A_g <= truth_ll (the pair can be neither above
//!       nor tied);
//!     - the best-outside fallback (the emitted best divergent rival
//!       and its tied set): the standard branch-and-bound argmax over
//!       the skipped non-excluded pairs, walked in bound-descending
//!       order until the bound falls below the running max — its cost
//!       is bounded by the kept-side counts this probe receipts.
//!
//! THE HONEST ADMISSIBILITY QUESTION THIS PROBE PRICES WITH NUMBERS:
//! the committed gate set (the receipts' own bar) includes the
//! INSTRUMENT RECEIPT's exhaustive `class_log_likelihoods` vector and
//! the ingredients' full `ll_matrix` — a pair-excluding bound changes
//! those bytes BY CONSTRUCTION. This probe measures at the worst
//! shapes: (a) the ANSWER-window pruning rate (the bound's exclusion
//! count at the derived thresholds), (b) the exact-value-preserving
//! memoization rate (bit-equal rows share one sum — the only lever
//! that CAN pass the committed bar), and (c) the bound's own cost —
//! plus the winner pair's bit-identity against the committed receipt
//! at every probed locus (the parse's own end-to-end verification).
//!
//! Usage: class_bound_probe <committed-call-dir> [--loci 0,10] [--out <jsonl>]
//!   <committed-call-dir> = a call output dir carrying
//!   instrument-receipt.jsonl, instrument-receipt.jsonl.ingredients.jsonl,
//!   and calls.jsonl.truth-qv.jsonl (the tranche's callC/<TID> shape).
//!   The monster loci's ingredients lines reach ~12GB; everything is
//!   streamed in chunks, nothing held whole.

use std::borrow::Cow;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use std::path::PathBuf;

/// The committed 1/2-1/2 homolog mixture (realign_score.rs's own
/// closed form, reimplemented here for the offline verification).
fn mix_logsumexp(a: f64, b: f64) -> f64 {
    if a == f64::NEG_INFINITY && b == f64::NEG_INFINITY {
        return f64::NEG_INFINITY;
    }
    let m = a.max(b);
    m + (0.5 * (a - m).exp() + 0.5 * (b - m).exp()).ln()
}

/// The ulp distance between two f64s on the monotone bit axis (0 for
/// equal bits) — realign_score.rs's own receipted delta unit.
fn ulp_distance(a: f64, b: f64) -> u64 {
    let key = |x: f64| {
        let bits = x.to_bits();
        if bits >> 63 == 0 {
            bits | 0x8000_0000_0000_0000
        } else {
            !bits
        }
    };
    key(a).abs_diff(key(b))
}

/// The f64 exp underflow point, DERIVED BY EVALUATION (never a tuned
/// constant): the largest f64 t with exp(t) == +0.0, walked down by
/// ulps from -1075*ln(2) — the real-arithmetic point below which the
/// rounded exponential is the exact +0.0. Every x <= t underflows, so
/// a pair skipped under this threshold contributes exactly +0.0 where
/// the exact run contributes exactly +0.0.
fn exp_underflow_point() -> f64 {
    let mut lo = -744.0f64;
    let mut hi = -1075.0f64 * 2.0f64.ln();
    assert!(lo.exp() > 0.0);
    assert_eq!(hi.exp(), 0.0);
    while hi.to_bits() - lo.to_bits() > 1 {
        // the midpoint on the bit axis WITHOUT the u64 sum overflow
        // (both anchors are negative doubles: their bit patterns sum
        // past u64::MAX).
        let mid = f64::from_bits(lo.to_bits() + (hi.to_bits() - lo.to_bits()) / 2);
        if mid == lo || mid == hi {
            break;
        }
        if mid.exp() == 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    hi
}

fn parse_args() -> (PathBuf, Vec<usize>, Option<PathBuf>) {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!(
            "usage: class_bound_probe <committed-call-dir> [--loci 0,10] [--out <jsonl>]"
        );
        std::process::exit(2);
    }
    let mut dir = None;
    let mut loci = Vec::new();
    let mut out = None;
    let mut i = 1;
    while i < args.len() {
        match args[i].as_str() {
            "--loci" => {
                i += 1;
                for part in args[i].split(',') {
                    loci.push(part.parse().expect("bad locus id"));
                }
            }
            "--out" => {
                i += 1;
                out = Some(PathBuf::from(&args[i]));
            }
            _ => dir = Some(PathBuf::from(&args[i])),
        }
        i += 1;
    }
    (dir.expect("missing call dir"), loci, out)
}

/// One ingredients line, streamed in chunks (the monster loci's lines
/// reach ~12GB; nothing is held whole).
struct LineChunks {
    reader: BufReader<File>,
    buf: Vec<u8>,
    filled: usize,
    pos: usize,
    line_done: bool,
}

impl LineChunks {
    fn open(path: &std::path::Path) -> std::io::Result<LineChunks> {
        Ok(LineChunks {
            reader: BufReader::with_capacity(8 << 20, File::open(path)?),
            buf: vec![0; 16 << 20],
            filled: 0,
            pos: 0,
            line_done: false,
        })
    }

    /// Skip `count` whole lines; false at EOF.
    fn skip_lines(&mut self, count: usize) -> std::io::Result<bool> {
        let mut remaining = count;
        loop {
            if remaining == 0 {
                return Ok(true);
            }
            if self.pos >= self.filled {
                self.filled = self.reader.read(&mut self.buf)?;
                self.pos = 0;
                if self.filled == 0 {
                    return Ok(false);
                }
            }
            let hay = &self.buf[self.pos..self.filled];
            match hay.iter().position(|&b| b == b'\n') {
                Some(idx) => {
                    self.pos += idx + 1;
                    remaining -= 1;
                    if remaining == 0 {
                        return Ok(true);
                    }
                }
                None => self.pos = self.filled,
            }
        }
    }

    /// The next chunk of the current line (empty slice = line ended).
    /// The line is DONE only when a '\n' has actually been seen — a
    /// chunk that ends because the BUFFER filled is mid-line.
    fn next_chunk(&mut self) -> std::io::Result<&[u8]> {
        loop {
            if self.pos < self.filled {
                let start = self.pos;
                match self.buf[self.pos..self.filled]
                    .iter()
                    .position(|&b| b == b'\n')
                {
                    Some(idx) => {
                        self.pos += idx + 1;
                        self.line_done = true;
                        return Ok(&self.buf[start..start + idx]);
                    }
                    None => {
                        self.pos = self.filled;
                        return Ok(&self.buf[start..self.filled]);
                    }
                }
            }
            if self.line_done {
                return Ok(&[]);
            }
            self.filled = self.reader.read(&mut self.buf)?;
            self.pos = 0;
            if self.filled == 0 {
                self.line_done = true;
                return Ok(&[]);
            }
        }
    }
}

/// A marker-aware scanner over ONE ingredients line: `find_marker`
/// leaves the scanner's remaining bytes EXACTLY after the marker
/// (the overshoot after the marker is buffered and served first by
/// `next_bytes`).
struct LineScanner {
    chunks: LineChunks,
    pending: Vec<u8>,
}

impl LineScanner {
    fn open(path: &std::path::Path) -> std::io::Result<LineScanner> {
        Ok(LineScanner {
            chunks: LineChunks::open(path)?,
            pending: Vec::new(),
        })
    }

    fn skip_lines(&mut self, count: usize) -> std::io::Result<bool> {
        self.chunks.skip_lines(count)
    }

    /// The bytes of the current line up to and including the first
    /// `stop` byte; the REMAINDER of the buffer stays in the scanner
    /// (nothing is discarded mid-chunk). Empty result = line ended.
    fn token_until(&mut self, stop: u8) -> std::io::Result<Vec<u8>> {
        let mut out = Vec::new();
        loop {
            if self.pending.is_empty() {
                let chunk = self.chunks.next_chunk()?;
                if chunk.is_empty() {
                    return Ok(out);
                }
                self.pending = chunk.to_vec();
            }
            let mut idx = 0usize;
            while idx < self.pending.len() {
                let b = self.pending[idx];
                idx += 1;
                out.push(b);
                if b == stop {
                    self.pending.drain(..idx);
                    return Ok(out);
                }
            }
            self.pending.clear();
        }
    }

    /// The scanner's remaining bytes, pending-first (empty = the line
    /// ended).
    fn next_bytes(&mut self) -> std::io::Result<Cow<'_, [u8]>> {
        if !self.pending.is_empty() {
            return Ok(Cow::Owned(std::mem::take(&mut self.pending)));
        }
        Ok(Cow::Borrowed(self.chunks.next_chunk()?))
    }

    /// Find the marker in the line's remaining bytes; on success the
    /// scanner's remaining bytes start exactly after the marker. The
    /// scan window is ALWAYS pending ++ the next chunk (the pending
    /// holds either a post-marker leftover or a straddle carry — both
    /// are data to scan, and a new chunk is pulled every iteration so
    /// the window always advances).
    fn find_marker(&mut self, marker: &[u8]) -> std::io::Result<bool> {
        loop {
            let chunk = self.chunks.next_chunk()?;
            if chunk.is_empty() {
                self.pending.clear();
                return Ok(false);
            }
            let mut hay: Vec<u8> = Vec::with_capacity(self.pending.len() + chunk.len());
            hay.extend_from_slice(&self.pending);
            hay.extend_from_slice(chunk);
            let mut from = 0usize;
            let mut found = None;
            while let Some(rel) = hay[from..].iter().position(|&b| b == marker[0]) {
                let at = from + rel;
                if hay[at..].starts_with(marker) {
                    found = Some(at);
                    break;
                }
                from = at + 1;
            }
            if let Some(at) = found {
                self.pending = hay[at + marker.len()..].to_vec();
                return Ok(true);
            }
            // Keep a tail that could straddle the next chunk.
            let keep = marker.len().saturating_sub(1).min(hay.len());
            self.pending = hay[hay.len() - keep..].to_vec();
        }
    }
}

/// The per-row streaming accumulator: the excess mass A_f, the pinned
/// count, the 128-bit token hash (serde_json float tokens are
/// shortest-roundtrip, so token equality IS bit equality), and — for
/// the winner pair's two verification rows — the full dense values.
struct RowAgg {
    unit: usize,
    alpha_sum: f64,
    pinned: u64,
    hash_a: u64,
    hash_b: u64,
    elsewhere: f64,
    elsewhere_token: String,
    values: Option<Vec<f64>>,
}

impl RowAgg {
    fn new(elsewhere: f64, keep_values: bool) -> RowAgg {
        RowAgg {
            unit: 0,
            alpha_sum: 0.0,
            pinned: 0,
            hash_a: 0xcbf2_9ce4_8422_2325,
            hash_b: 0x9e37_79b9_7f4a_7c15,
            elsewhere,
            elsewhere_token: serde_json::to_string(&elsewhere).expect("E token"),
            values: keep_values.then(Vec::new),
        }
    }

    fn hash_byte(hash: u64, byte: u8) -> u64 {
        (hash ^ byte as u64).wrapping_mul(0x100_0000_01b3)
    }

    /// One ll token (the row hash folds a separator byte first, so two
    /// rows with the same value sequence hash equal whichever chunk
    /// boundaries split them).
    fn push_token(&mut self, token: &[u8], counts: &[f64]) {
        self.hash_a = Self::hash_byte(self.hash_a, b',');
        self.hash_b = Self::hash_byte(self.hash_b, b',');
        for &b in token {
            self.hash_a = Self::hash_byte(self.hash_a, b);
            self.hash_b = Self::hash_byte(self.hash_b, b);
        }
        if token == self.elsewhere_token.as_bytes() {
            if let Some(values) = self.values.as_mut() {
                values.push(self.elsewhere);
            }
        } else {
            self.pinned += 1;
            let value: f64 = std::str::from_utf8(token)
                .expect("utf8 ll token")
                .parse()
                .expect("ll float token");
            self.alpha_sum += counts[self.unit] * (value - self.elsewhere);
            if let Some(values) = self.values.as_mut() {
                values.push(value);
            }
        }
        self.unit += 1;
    }
}

/// Count unordered pairs f <= g with excess[f] + excess[g] >= need
/// (inclusive) or > need (strict), by sort + partition point.
fn count_pairs(excess: &[f64], need: f64, strict: bool) -> usize {
    let mut sorted: Vec<f64> = excess.to_vec();
    sorted.sort_by(|a, b| a.partial_cmp(b).expect("finite excess"));
    let mut ordered = 0usize;
    for &a in excess {
        let want = need - a;
        let j = if strict {
            sorted.partition_point(|&x| x <= want)
        } else {
            sorted.partition_point(|&x| x < want)
        };
        ordered += sorted.len() - j;
    }
    // The diagonal (f == f) is counted once in `ordered`; unordered
    // f <= g counts it once too, so: unordered = (ordered + diag) / 2.
    let diag = if strict {
        excess.iter().filter(|&&a| 2.0 * a > need).count()
    } else {
        excess.iter().filter(|&&a| 2.0 * a >= need).count()
    };
    (ordered + diag) / 2
}

fn main() -> std::io::Result<()> {
    let (dir, loci, out_path) = parse_args();
    let exp_underflow = exp_underflow_point();
    eprintln!(
        "[probe] derived exp underflow point: t = {exp_underflow:.17} \
         (exp(t) == 0.0 exactly)"
    );

    // The instrument receipt: the exhaustive class vector + E per locus.
    let receipt_path = dir.join("instrument-receipt.jsonl");
    let mut receipt_by_locus: Vec<serde_json::Value> = Vec::new();
    for line in BufReader::new(File::open(&receipt_path)?).lines() {
        receipt_by_locus.push(serde_json::from_str(&line?)?);
    }
    // The truth-qv assessment artifact: the TEST truth pair + rank.
    let truth_path = dir.join("calls.jsonl.truth-qv.jsonl");
    let mut truth_by_locus: Vec<(Option<(usize, usize)>, Option<u64>)> = Vec::new();
    for line in BufReader::new(File::open(&truth_path)?).lines() {
        let value: serde_json::Value = serde_json::from_str(&line?)?;
        let pair = value["truth_folds"].as_array().and_then(|folds| {
            let a = folds.first()?["index"].as_u64()? as usize;
            let b = folds.get(1)?["index"].as_u64()? as usize;
            Some((a.min(b), a.max(b)))
        });
        truth_by_locus.push((pair, value["truth_rank"].as_u64()));
    }

    let requested: Vec<usize> = if loci.is_empty() {
        receipt_by_locus
            .iter()
            .map(|v| v["locus"].as_u64().expect("locus id") as usize)
            .collect()
    } else {
        loci
    };

    let mut out_file = match &out_path {
        Some(path) => Some(BufWriter::new(File::create(path)?)),
        None => None,
    };
    let stdout = std::io::stdout();
    let mut stdout = stdout.lock();

    let ingredients_path = dir.join("instrument-receipt.jsonl.ingredients.jsonl");
    let member = dir
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();

    let mut weighted: Vec<(u64, f64)> = Vec::new();
    for &locus in &requested {
        let receipt = &receipt_by_locus[locus];
        let elsewhere = receipt["scoring"]["elsewhere_log_prob"]
            .as_f64()
            .expect("elsewhere");
        let folds = receipt["folds"].as_u64().expect("folds") as usize;
        let units = receipt["unit_count"].as_u64().expect("units") as usize;
        let pairs = folds * (folds + 1) / 2;
        let mut class_lls: Vec<f64> = Vec::with_capacity(pairs);
        let mut class_nulls = 0usize;
        for value in receipt["class_log_likelihoods"].as_array().expect("lls") {
            match value.as_f64() {
                Some(ll) => class_lls.push(ll),
                None => {
                    class_nulls += 1;
                    class_lls.push(f64::NEG_INFINITY);
                }
            }
        }
        let finite: Vec<f64> = class_lls
            .iter()
            .copied()
            .filter(|ll| ll.is_finite())
            .collect();
        let best = finite.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        let winner = receipt["best_fold_indices"].as_array().expect("winner");
        let winner_first = winner[0].as_u64().expect("wf") as usize;
        let winner_second = winner[1].as_u64().expect("ws") as usize;

        // ---------------- PASS A: unit counts + the locus verification.
        let mut scanner = LineScanner::open(&ingredients_path)?;
        if !scanner.skip_lines(locus)? {
            return Err(std::io::Error::other("ingredients ended before locus"));
        }
        if !scanner.find_marker(b"\"ll_matrix\":")? {
            return Err(std::io::Error::other("ll_matrix marker missing"));
        }
        if !scanner.find_marker(b",\"locus\":")? {
            return Err(std::io::Error::other("locus marker missing"));
        }
        // The locus value: digits until the next ','. Nothing after the
        // comma is discarded (token_until keeps the remainder).
        let locus_token = scanner.token_until(b',')?;
        let locus_digits: Vec<u8> = locus_token
            .into_iter()
            .filter(|b| b.is_ascii_digit())
            .collect();
        let locus_seen: usize = std::str::from_utf8(&locus_digits)
            .expect("digits")
            .parse()
            .unwrap_or(usize::MAX);
        if locus_seen != locus {
            return Err(std::io::Error::other(format!(
                "ingredients line {locus} carries locus {locus_seen} — \
                 the line-index assumption broke; refusing to guess"
            )));
        }
        if !scanner.find_marker(b"\"units\":")? {
            return Err(std::io::Error::other("units marker missing"));
        }
        // The unit counts: a rolling match on "count": then the
        // integer up to ',' or '}' — chunk-boundary safe.
        const NEEDLE: &[u8] = b"\"count\":";
        let mut state = 0usize;
        let mut value: f64 = 0.0f64;
        let mut saw_digit = false;
        let mut counts: Vec<f64> = Vec::with_capacity(units);
        'units: loop {
            let bytes = scanner.next_bytes()?;
            if bytes.is_empty() {
                break;
            }
            for &b in bytes.iter() {
                if state < NEEDLE.len() {
                    if b == NEEDLE[state] {
                        state += 1;
                    } else {
                        state = usize::from(b == NEEDLE[0]);
                    }
                    continue;
                }
                if b.is_ascii_digit() {
                    value = value * 10.0 + f64::from(b - b'0');
                    saw_digit = true;
                } else if b == b',' || b == b'}' {
                    if saw_digit {
                        counts.push(value);
                        if counts.len() == units {
                            break 'units;
                        }
                    }
                    value = 0.0;
                    saw_digit = false;
                    state = 0;
                }
            }
        }
        if counts.len() != units {
            return Err(std::io::Error::other(format!(
                "unit counts: expected {units}, saw {}",
                counts.len()
            )));
        }
        let c_total: f64 = counts.iter().sum();

        // ---------------- PASS B: the ll rows (aggregates; the winner
        // pair's two rows kept in full for the bit-verification).
        let mut scanner = LineScanner::open(&ingredients_path)?;
        if !scanner.skip_lines(locus)? {
            return Err(std::io::Error::other("second pass skip failed"));
        }
        if !scanner.find_marker(b"\"ll_matrix\":")? {
            return Err(std::io::Error::other("ll_matrix marker missing (pass B)"));
        }
        let keep = |row: usize| row == winner_first || row == winner_second;
        let mut excess: Vec<f64> = Vec::with_capacity(folds);
        let mut pinned: Vec<u64> = Vec::with_capacity(folds);
        let mut hash_a: Vec<u64> = Vec::with_capacity(folds);
        let mut hash_b: Vec<u64> = Vec::with_capacity(folds);
        let mut kept_rows: Vec<Vec<f64>> = Vec::with_capacity(2);
        let mut agg = RowAgg::new(elsewhere, keep(0));
        let mut row = 0usize;
        let mut token: Vec<u8> = Vec::with_capacity(24);
        let mut matrix_ended = false;
        'matrix: loop {
            let bytes = scanner.next_bytes()?;
            if bytes.is_empty() {
                break;
            }
            for &b in bytes.iter() {
                match b {
                    b',' | b']' => {
                        if !token.is_empty() {
                            agg.push_token(&token, &counts);
                            token.clear();
                        }
                        if b == b']' && agg.unit == units {
                            excess.push(agg.alpha_sum);
                            pinned.push(agg.pinned);
                            hash_a.push(agg.hash_a);
                            hash_b.push(agg.hash_b);
                            if let Some(values) = agg.values.take() {
                                kept_rows.push(values);
                            }
                            row += 1;
                            if row == folds {
                                matrix_ended = true;
                                break 'matrix;
                            }
                            agg = RowAgg::new(elsewhere, keep(row));
                        }
                    }
                    b'[' | b' ' => {}
                    _ => token.push(b),
                }
            }
        }
        if !matrix_ended || row != folds {
            return Err(std::io::Error::other(format!(
                "ll_matrix parse incomplete: rows {row} of {folds}"
            )));
        }

        // The winner-pair verification: the sum rebuilt with the
        // committed accumulation from the streamed rows, vs the
        // receipt's own class text (a hom pair uses the one kept row
        // twice — exactly the committed accumulation's own input).
        // THE NAMED RECEIPT-SERIALIZATION CLASS, stated for the
        // interpretation: this repo's serde_json carries NO
        // float_roundtrip feature, so ~31.5% of arbitrary f64s emit
        // JSON text that parses back ONE ULP OFF (unit-receipted);
        // the committed ll_matrix text therefore does not always
        // carry the scan's exact in-memory bits, and an offline
        // re-accumulation can differ from the emitted class sum by
        // 1-2 ulps at pairs whose accumulation crosses affected
        // values — the committed BINARY reproduces its own receipt
        // bit-identically (the rerun receipt), so the class is the
        // emission's, never the scan's or the parse's. The
        // verification therefore reports the ULP DISTANCE and its
        // class: 0 = bit-identical; <= 2 = the receipt-serialization
        // class (the parse is validated, the answer-window analysis
        // is unaffected — 2 ulps against 745-nat thresholds); larger
        // = a genuine parse failure (none observed).
        let expected_kept = 2 - usize::from(winner_first == winner_second);
        let verify_ulp: Option<u64> = (kept_rows.len() == expected_kept).then(|| {
            let (a_values, b_values) = if winner_first == winner_second {
                (&kept_rows[0], &kept_rows[0])
            } else {
                (&kept_rows[0], &kept_rows[1])
            };
            let mut total = 0.0f64;
            for (u, &count) in counts.iter().enumerate() {
                total += count * mix_logsumexp(a_values[u], b_values[u]);
            }
            let flat = winner_second * (winner_second + 1) / 2 + winner_first;
            ulp_distance(total, class_lls[flat])
        });
        let verify = match verify_ulp {
            Some(0) => "bit-identical",
            Some(1..=2) => "receipt-serialization-ulps",
            Some(ulps) => {
                return Err(std::io::Error::other(format!(
                    "winner-pair recomputation off by {ulps} ulps — a genuine                      parse failure, refusing to emit an unvalidated receipt"
                )));
            }
            None => {
                return Err(std::io::Error::other(
                    "the winner pair's rows were not both kept — internal error",
                ));
            }
        };
        let verify_ulp = verify_ulp.unwrap_or(u64::MAX);

        // ---------------- The counterfactual arithmetic.
        let floor = elsewhere * c_total;
        let winner_excess = best - floor;
        let within_underflow = finite
            .iter()
            .filter(|&&ll| ll > best + exp_underflow)
            .count();
        let second_gap = {
            let mut second = f64::NEG_INFINITY;
            for &ll in &finite {
                if ll < best && ll > second {
                    second = ll;
                }
            }
            best - second
        };
        // Production answer window: keep iff E*C + A_f + A_g >= best + t.
        let kept_production = count_pairs(&excess, best + exp_underflow - floor, false);
        // Assessment window: keep iff the bound is strictly above the
        // truth TEST pair's own class sum (the rank count needs every
        // pair that could sit above it).
        let truth = truth_by_locus.get(locus).and_then(|(p, _)| *p);
        let truth_rank_test = truth_by_locus.get(locus).and_then(|(_, r)| *r);
        let truth_ll_test = truth.map(|(a, b)| class_lls[b * (b + 1) / 2 + a]);
        let kept_assessment = truth_ll_test
            .filter(|ll| ll.is_finite())
            .map(|ll| count_pairs(&excess, ll - floor, true));
        // The exact-preserving memoization rate: bit-equal rows share
        // one sum per unordered (row-class, row-class) pair.
        let mut classes: Vec<(u64, u64)> = hash_a
            .iter()
            .zip(hash_b.iter())
            .map(|(&a, &b)| (a, b))
            .collect();
        classes.sort();
        classes.dedup();
        let distinct_rows = classes.len();
        let distinct_pairs = distinct_rows * (distinct_rows + 1) / 2;

        let pinned_entries: u64 = pinned.iter().sum();
        let fold_unit_pairs = folds as u64 * units as u64;
        let scan_unit_pairs = pairs as u64 * units as u64;
        let pruned = pairs - kept_production.min(pairs);
        let rate = pruned as f64 / pairs as f64;
        let mut excess_desc = excess.clone();
        excess_desc.sort_by(|a, b| b.partial_cmp(a).expect("finite excess"));
        let record = serde_json::json!({
            "member": member,
            "locus": locus,
            "folds": folds,
            "units": units,
            "pairs": pairs,
            "class_lls_nulls": class_nulls,
            "elsewhere": elsewhere,
            "c_total": c_total,
            "winner_pair": [winner_first, winner_second],
            "landscape": {
                "best": best,
                "second_gap": second_gap,
                "within_exp_underflow": within_underflow,
                "winner_excess_over_E_floor": winner_excess,
            },
            "pinned": {
                "fold_unit_pairs": fold_unit_pairs,
                "pinned_entries": pinned_entries,
                "pinned_share": pinned_entries as f64 / fold_unit_pairs as f64,
                "per_row_max": pinned.iter().copied().max().unwrap_or(0),
            },
            "excess_mass": {
                "a_max": excess.iter().copied().fold(0.0f64, f64::max),
                "a_median": excess_desc[excess_desc.len() / 2],
                "a_p10": excess_desc[excess_desc.len() / 10],
                "a_p90": excess_desc[excess_desc.len() - 1 - excess_desc.len() / 10],
            },
            "counterfactual_bound": {
                "kept_production": kept_production,
                "pruned_production": pruned,
                "pruning_rate_production": rate,
                "kept_assessment": kept_assessment,
                "truth_ll_test": truth_ll_test,
                "truth_rank_test": truth_rank_test,
                "bound_walk_pinned_entries": pinned_entries,
                "scan_unit_pairs": scan_unit_pairs,
                "bound_walk_over_scan":
                    pinned_entries as f64 / scan_unit_pairs as f64,
            },
            "memoization_exact_preserving": {
                "distinct_rows": distinct_rows,
                "distinct_pairs": distinct_pairs,
                "pruning_rate": 1.0 - distinct_pairs as f64 / pairs as f64,
            },
            "verification": {
                "winner_pair_recomputation_ulp_distance": verify_ulp,
                "winner_pair_recomputation_class": verify,
            },
        });
        let line = serde_json::to_string(&record)?;
        writeln!(stdout, "{line}")?;
        stdout.flush()?;
        if let Some(writer) = out_file.as_mut() {
            writeln!(writer, "{line}")?;
        }
        eprintln!(
            "[probe] {member} locus {locus}: pairs {pairs}, kept {kept_production} \
             (rate {rate:.4}), memo {:.4}, pinned {:.4}, verify {verify} ({verify_ulp} ulp)",
            1.0 - distinct_pairs as f64 / pairs as f64,
            pinned_entries as f64 / fold_unit_pairs as f64,
        );
        weighted.push((scan_unit_pairs, rate));
    }
    if let Some(writer) = out_file.as_mut() {
        writer.flush()?;
    }
    let total_weight: u64 = weighted.iter().map(|(w, _)| *w).sum();
    let weighted_rate: f64 = weighted
        .iter()
        .map(|(w, r)| *w as f64 * r)
        .sum::<f64>()
        / total_weight.max(1) as f64;
    let summary = serde_json::json!({
        "member": member,
        "loci_measured": weighted.len(),
        "scan_unit_pairs_total": total_weight,
        "wall_weighted_pruning_rate_production": weighted_rate,
        "exp_underflow_point": exp_underflow,
    });
    let line = serde_json::to_string(&summary)?;
    writeln!(stdout, "{line}")?;
    if let Some(mut writer) = out_file.take() {
        writeln!(writer, "{line}")?;
        writer.flush()?;
    }
    Ok(())
}
