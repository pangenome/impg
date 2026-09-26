//! Fixed-parameter finite-hull mathematics, not a genome representation.
use serde::Serialize;
use std::io;
pub const BETA: f64 = 0.1;
pub const DEPTH: f64 = 10.0;
pub const TIE: f64 = 1e-9;
pub const CHECK: f64 = 1e-7;
pub const GAP_STOP: f64 = 1e-9;
pub const MAX_PRODUCTS: usize = 4096;
pub const MAX_COLUMNS: usize = 1024;
pub const MAX_FEATURES: usize = 4096;
pub const MAX_CELLS: usize = 4_194_304;
pub const MAX_WORK: usize = 600_000_000;
pub const MAX_BYTES: usize = 128 * 1024 * 1024;
pub const ITERATIONS: usize = 128;
pub const BISECTIONS: usize = 64;
pub fn require(ok: bool, message: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(io::Error::other(message))
    }
}
pub fn close(a: f64, b: f64) -> io::Result<()> {
    require(
        a.is_finite() && b.is_finite() && (a - b).abs() <= CHECK * (1.0 + a.abs().max(b.abs())),
        "numerical identity/bound check failed",
    )
}
pub fn loss(c: &[f64], s: &[f64]) -> f64 {
    c.iter()
        .zip(s)
        .map(|(&c, &s)| s - c * (s / BETA).ln_1p())
        .sum()
}
pub fn gradient(c: &[f64], s: &[f64]) -> Vec<f64> {
    c.iter()
        .zip(s)
        .map(|(&c, &s)| 1.0 - c / (BETA + s))
        .collect()
}
pub fn delta(c: &[f64], s: &[f64], t: &[f64]) -> f64 {
    c.iter()
        .zip(s)
        .zip(t)
        .map(|((&c, &s), &t)| (t - s) - c * ((t - s) / (BETA + s)).ln_1p())
        .sum()
}
pub fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
#[derive(Debug, Serialize)]
pub struct Pricing {
    pub chosen: usize,
    pub prices: Vec<f64>,
    pub raw_gap: Option<f64>,
    pub numerical_tangent_lower_bound: Option<f64>,
    pub full_explicit_universe: bool,
}
/// Incomplete candidate pricing supplies an upper bound on the pricing minimum,
/// so even its apparent zero gap must not be exposed as a full-domain bound.
pub fn price(c: &[f64], s: &[f64], columns: &[Vec<f64>], complete: bool) -> Pricing {
    let g = gradient(c, s);
    let prices: Vec<_> = columns.iter().map(|v| dot(&g, v)).collect();
    let chosen = (0..prices.len())
        .min_by(|&a, &b| prices[a].total_cmp(&prices[b]).then(a.cmp(&b)))
        .unwrap();
    let gap = dot(&g, s) - prices[chosen];
    Pricing {
        chosen,
        prices,
        raw_gap: complete.then_some(gap),
        numerical_tangent_lower_bound: complete.then_some(loss(c, s) - gap),
        full_explicit_universe: complete,
    }
}
#[derive(Debug, Serialize)]
pub struct Trace {
    pub iteration: usize,
    pub signal: Vec<f64>,
    pub residual: Vec<f64>,
    pub convex_weights_not_genomes: Vec<f64>,
    pub objective: f64,
    pub pricing: Pricing,
    pub lower_bound_minus_discrete_minimum: Option<f64>,
    pub step_to_next: f64,
}
#[derive(Debug, Serialize)]
pub struct Relaxation {
    pub status: String,
    pub numerical_scope: String,
    pub rigorous_real_certificate: bool,
    pub discrete_minimum: f64,
    pub correlated_optimum_column_indices: Vec<usize>,
    pub scalar_work: usize,
    pub logical_numeric_bytes: usize,
    pub trace: Vec<Trace>,
}
pub fn optimize(
    c: &[f64],
    columns: &[Vec<f64>],
    public_scores: &[f64],
    complete: bool,
) -> io::Result<Relaxation> {
    require(
        !c.is_empty()
            && c.len() <= MAX_FEATURES
            && !columns.is_empty()
            && columns.len() <= MAX_COLUMNS,
        "feature/column cap or empty universe",
    )?;
    let cells = c.len() * columns.len();
    // Includes integer profiles, signals, and all retained trace numeric vectors.
    let bytes = 24 * cells + 8 * ITERATIONS * (4 * c.len() + 2 * columns.len());
    require(
        cells <= MAX_CELLS && bytes <= MAX_BYTES,
        "numeric cell/memory cap; no full-universe result",
    )?;
    require(
        public_scores.len() == columns.len() && c.iter().all(|v| v.is_finite() && *v >= 0.0),
        "invalid observations/scores",
    )?;
    for (v, &score) in columns.iter().zip(public_scores) {
        require(
            v.len() == c.len() && v.iter().all(|s| s.is_finite() && *s >= 0.0),
            "invalid signal column",
        )?;
        close(loss(c, v), score)?;
    }
    let best = public_scores
        .iter()
        .copied()
        .min_by(f64::total_cmp)
        .unwrap();
    let support = if complete {
        public_scores
            .iter()
            .enumerate()
            .filter_map(|(i, &v)| ((v - best).abs() <= TIE).then_some(i))
            .collect()
    } else {
        vec![]
    };
    let mut s = columns[0].clone();
    let mut weights = vec![0.0; columns.len()];
    weights[0] = 1.0;
    let mut trace = Vec::new();
    let mut work = 0;
    let mut status = "iteration_limit";
    for iteration in 0..ITERATIONS {
        let charge = cells + (2 * BISECTIONS + 12) * c.len();
        require(
            work + charge <= MAX_WORK,
            "scalar work cap; no full-universe result",
        )?;
        work += charge;
        let objective = loss(c, &s);
        let residual = gradient(c, &s);
        let pricing = price(c, &s, columns, complete);
        if let Some(lb) = pricing.numerical_tangent_lower_bound {
            require(
                lb <= best + CHECK * (1.0 + best.abs()),
                "numerical tangent bound exceeds exhausted discrete minimum",
            )?;
        }
        let stop = pricing
            .raw_gap
            .is_some_and(|gap| gap <= GAP_STOP * (1.0 + objective.abs()));
        let target = &columns[pricing.chosen];
        let direction: Vec<_> = target.iter().zip(&s).map(|(t, s)| t - s).collect();
        let slope = |a: f64| {
            c.iter()
                .zip(&s)
                .zip(&direction)
                .map(|((&c, &s), &d)| (1.0 - c / (BETA + s + a * d)) * d)
                .sum::<f64>()
        };
        let step = if stop || iteration + 1 == ITERATIONS || slope(0.0) >= 0.0 {
            0.0
        } else if slope(1.0) <= 0.0 {
            1.0
        } else {
            let (mut lo, mut hi) = (0.0, 1.0);
            for _ in 0..BISECTIONS {
                let mid = (lo + hi) / 2.0;
                if slope(mid) < 0.0 {
                    lo = mid;
                } else {
                    hi = mid;
                }
            }
            (lo + hi) / 2.0
        };
        let discrepancy = pricing.numerical_tangent_lower_bound.map(|lb| lb - best);
        trace.push(Trace {
            iteration,
            signal: s.clone(),
            residual,
            convex_weights_not_genomes: weights.clone(),
            objective,
            pricing,
            lower_bound_minus_discrete_minimum: discrepancy,
            step_to_next: step,
        });
        if stop {
            status = "numerical_gap_stop";
            break;
        }
        if iteration + 1 == ITERATIONS {
            break;
        }
        let next: Vec<_> = s
            .iter()
            .zip(target)
            .map(|(s, t)| (1.0 - step) * s + step * t)
            .collect();
        require(
            loss(c, &next) <= objective + CHECK * (1.0 + objective.abs()),
            "line search increased loss",
        )?;
        for w in &mut weights {
            *w *= 1.0 - step;
        }
        weights[trace.last().unwrap().pricing.chosen] += step;
        s = next;
    }
    Ok(Relaxation {
        status: status.into(),
        numerical_scope: if complete {
            "f64 only; exhausted supplied finite layout, not genome-wide"
        } else {
            "incomplete pricing; no full-universe bound or support"
        }
        .into(),
        rigorous_real_certificate: false,
        discrete_minimum: best,
        correlated_optimum_column_indices: support,
        scalar_work: work,
        logical_numeric_bytes: bytes,
        trace,
    })
}
