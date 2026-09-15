//! Exact sparse GLOBAL deltas, not donor-local likelihood or uncharged pricing.
use super::*;
pub type Counts = Vec<BTreeMap<[u64; 3], u64>>;
#[derive(Debug)]
pub struct Baseline {
    pub assignment: routes::Assignment,
    pub score: f64,
    pub counts: Counts,
    pub epoch: u8,
}
pub fn baseline(e: routes::Evaluation, epoch: u8, lengths: usize) -> Baseline {
    let mut counts = vec![BTreeMap::new(); lengths];
    for f in e.factors {
        for (li, n) in f.counts_by_length.into_iter().enumerate() {
            if n > 0 {
                counts[li].insert(f.tokens, n);
            }
        }
    }
    Baseline {
        assignment: e.assignment,
        score: e.relative_objective,
        counts,
        epoch,
    }
}
#[derive(Serialize)]
pub struct Change {
    pub tokens: [u64; 3],
    pub observed: u64,
    pub old_counts: Vec<u64>,
    pub new_counts: Vec<u64>,
    pub old_signal: f64,
    pub new_signal: f64,
    pub exact_delta: f64,
    pub linear_delta: f64,
}
#[derive(Serialize)]
pub struct Delta {
    pub objective: f64,
    pub exact_delta: f64,
    pub linear_delta: f64,
    pub changes: Vec<Change>,
}
pub fn signal(q: &[u64], histogram: &[u64], denominator: u64) -> io::Result<f64> {
    let sum = q.iter().zip(histogram).try_fold(0u64, |n, (&q, &h)| {
        n.checked_add(
            q.checked_mul(h)
                .ok_or_else(|| invalid("count product overflow"))?,
        )
        .ok_or_else(|| invalid("count sum overflow"))
    })?;
    ensure(
        sum <= 1 << 53 && denominator <= 1 << 53,
        "integer exposure conversion precision limit",
    )?;
    Ok(10.0 * sum as f64 / denominator as f64)
}
pub fn delta(
    e: &mut routes::Evaluator<'_>,
    base: &Baseline,
    candidate: &routes::Assignment,
    observed: &BTreeMap<[u64; 3], u64>,
    histogram: &[u64],
    denominator: u64,
    meter: &mut Meter,
) -> io::Result<Delta> {
    let mut changes: BTreeMap<[u64; 3], Vec<u64>> = BTreeMap::new();
    for (old, new) in base.assignment.routes.iter().zip(&candidate.routes) {
        if old == new {
            continue;
        }
        for (route, subtract) in [(old, true), (new, false)] {
            meter.profile(route.length()?, histogram.len(), "explicit_changed_route")?;
            meter.explicit_profile_attempts += 1;
            let hits = e.cache_hits;
            let attempted = e.route_counts(route);
            meter.route_cache_hits += e.cache_hits - hits;
            let (counts, corrected) = attempted.map_err(public_feature_limit)?;
            meter.completed_route_profiles_observed += 1;
            meter.corrected_starts += corrected;
            for (li, row) in counts.iter().enumerate() {
                for (&tokens, &n) in row {
                    meter.work(1)?;
                    ensure(
                        changes.contains_key(&tokens) || changes.len() < meter.limits.max_features,
                        "budget: feature_terms",
                    )?;
                    let entry = changes.entry(tokens).or_insert_with(|| {
                        base.counts
                            .iter()
                            .map(|q| q.get(&tokens).copied().unwrap_or(0))
                            .collect()
                    });
                    entry[li] = if subtract {
                        entry[li].checked_sub(n)
                    } else {
                        entry[li].checked_add(n)
                    }
                    .ok_or_else(|| invalid("checked physical delta count overflow/underflow"))?;
                }
            }
        }
    }
    ensure(
        changes.len() <= meter.limits.max_features,
        "budget: feature_terms",
    )?;
    // Admit the whole known touched pass before charging an exact-score attempt.
    // Successful total work is unchanged; an interrupted pass is not an objective.
    meter.work(changes.len() as u64)?;
    meter.admit_score()?;
    let mut rows = Vec::with_capacity(changes.len());
    let (mut exact, mut linear) = (0.0, 0.0);
    for (tokens, new_counts) in changes {
        let old_counts: Vec<_> = base
            .counts
            .iter()
            .map(|q| q.get(&tokens).copied().unwrap_or(0))
            .collect();
        let old_signal = signal(&old_counts, histogram, denominator)?;
        let new_signal = signal(&new_counts, histogram, denominator)?;
        let observed = observed.get(&tokens).copied().unwrap_or(0);
        let d = new_signal - old_signal;
        let exact_delta = d - observed as f64 * (d / (0.1 + old_signal)).ln_1p();
        let linear_delta = (1.0 - observed as f64 / (0.1 + old_signal)) * d;
        ensure(
            exact_delta.is_finite() && linear_delta.is_finite(),
            "nonfinite sparse full-change loss",
        )?;
        exact += exact_delta;
        linear += linear_delta;
        rows.push(Change {
            tokens,
            observed,
            old_counts,
            new_counts,
            old_signal,
            new_signal,
            exact_delta,
            linear_delta,
        });
    }
    meter.exact_scores += 1;
    meter.exact_delta_scores += 1;
    Ok(Delta {
        objective: base.score + exact,
        exact_delta: exact,
        linear_delta: linear,
        changes: rows,
    })
}
pub fn parity(
    base: &Baseline,
    delta: &Delta,
    actual: &routes::Evaluation,
    observed: &BTreeMap<[u64; 3], u64>,
    meter: &mut Meter,
) -> io::Result<()> {
    close(delta.objective, actual.relative_objective)?;
    let changes: BTreeMap<_, _> = delta.changes.iter().map(|c| (c.tokens, c)).collect();
    let mut actual_keys = BTreeSet::new();
    for f in &actual.factors {
        meter.work(1)?;
        actual_keys.insert(f.tokens);
        let expected: Vec<_> = changes
            .get(&f.tokens)
            .map(|c| c.new_counts.clone())
            .unwrap_or_else(|| {
                base.counts
                    .iter()
                    .map(|q| q.get(&f.tokens).copied().unwrap_or(0))
                    .collect()
            });
        ensure(
            expected == f.counts_by_length
                && observed.get(&f.tokens).copied().unwrap_or(0) == f.observed,
            "public complete integer count parity failed",
        )?;
    }
    for c in &delta.changes {
        ensure(
            actual_keys.contains(&c.tokens)
                || (c.observed == 0 && c.new_counts.iter().all(|&n| n == 0)),
            "missing public realized/positive factor",
        )?;
    }
    Ok(())
}
