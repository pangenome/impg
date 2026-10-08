//! Complete pair loss; public haploid losses are charged but NEVER added.
use super::*;
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct Pair {
    pub copies: [routes::Assignment; 2],
}
impl Pair {
    /// A comparison key only: never reorder an immutable baseline's slots.
    pub fn key(&self) -> Self {
        let mut key = self.clone();
        key.copies.sort();
        key
    }
}
#[derive(Serialize)]
pub struct Factor {
    pub tokens: [u64; 3],
    pub counts: u64,
    pub observed: u64,
    pub signal: f64,
    pub loss: f64,
}
#[derive(Clone, Serialize)]
pub struct Scored {
    pub pair: Pair,
    pub objective: f64,
}
pub fn add(a: u64, b: u64) -> io::Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| invalid("paired integer count overflow"))
}
pub fn loss(q: u64, observed: u64, histogram: u64) -> io::Result<(f64, f64)> {
    ensure(
        observed <= 1 << 53 && histogram > 0,
        "paired count precision/histogram",
    )?;
    let denominator = histogram
        .checked_mul(150)
        .ok_or_else(|| invalid("histogram overflow"))?;
    let signal = haploid::score::signal(&[q], &[histogram], denominator)?;
    let loss = signal - observed as f64 * (signal / 0.1).ln_1p();
    ensure(loss.is_finite(), "nonfinite paired loss")?;
    Ok((signal, loss))
}
pub struct Scorer<'a, 'b> {
    pub e: &'a mut routes::Evaluator<'b>,
    pub meter: &'a mut haploid::Meter,
    pub ledger: &'a mut File,
    pub pair_attempts: &'a mut u64,
    pub pair_completed: &'a mut u64,
    pub pair_cap: u64,
    pub histogram: u64,
}
impl Scorer<'_, '_> {
    pub fn evaluate(&mut self, pair: &Pair, native: bool) -> io::Result<Scored> {
        for copy in &pair.copies {
            ensure(copy.routes.len() <= 4, "budget: paired_molecules")?;
            for route in &copy.routes {
                ensure(route.length()? <= 4096, "budget: paired_molecule_bp")?;
            }
        }
        ensure(*self.pair_attempts < self.pair_cap, "budget: paired_scores")?;
        *self.pair_attempts += 1;
        let ordinal = *self.pair_attempts;
        line(
            self.ledger,
            &json!({"kind":"pair_admission","ordinal":ordinal,"pair":pair,"native":native}),
        )?;
        let mut counts = BTreeMap::new();
        let mut normalized = pair.clone();
        for (slot, copy) in pair.copies.iter().enumerate() {
            self.meter.score_room(1)?;
            self.meter
                .validation(copy, self.e.graph.k, self.e.graph.port_count)?;
            let length = copy
                .routes
                .iter()
                .try_fold(0u64, |s, r| add(s, r.length()?))?;
            self.meter.profile(length, 1, "paired_public_haploid")?;
            self.meter.admit_score()?;
            self.meter.public_score_attempts += 1;
            if native {
                self.meter.native_score_attempts += 1;
            }
            self.meter.public_evaluation_route_reservations += copy.routes.len() as u64;
            let hits = self.e.cache_hits;
            let result = self.e.evaluate(copy);
            self.meter.route_cache_hits += self.e.cache_hits - hits;
            let ev = result.map_err(haploid::public_feature_limit)?;
            self.meter.exact_scores += 1;
            self.meter.public_scores += 1;
            if native {
                self.meter.native_scores += 1;
            }
            self.meter.completed_route_profiles_observed += copy.routes.len() as u64;
            self.meter.corrected_starts += ev.corrected_starts;
            line(
                self.ledger,
                &json!({"kind":"public_haploid_return","pair_ordinal":ordinal,"copy_slot":slot,"public":ev,"is_paired_score":false}),
            )?;
            normalized.copies[slot] = ev.assignment;
            for f in ev.factors {
                self.meter.work(1)?;
                ensure(f.counts_by_length.len() == 1, "paired L150 factor shape")?;
                ensure(
                    counts.contains_key(&f.tokens) || counts.len() < self.meter.limits.max_features,
                    "budget: paired_feature_union",
                )?;
                let entry = counts.entry(f.tokens).or_insert((0, f.observed));
                ensure(entry.1 == f.observed, "observations changed across copies")?;
                entry.0 = add(entry.0, f.counts_by_length[0])?;
            }
        }
        // Admit the complete known pass; failed public/union work is not an objective.
        self.meter.work(counts.len() as u64)?;
        let mut factors = Vec::with_capacity(counts.len());
        let mut objective = 0.0;
        for (tokens, (q, c)) in counts {
            let (signal, loss) = loss(q, c, self.histogram)?;
            objective += loss;
            factors.push(Factor {
                tokens,
                counts: q,
                observed: c,
                signal,
                loss,
            });
        }
        ensure(objective.is_finite(), "nonfinite complete paired objective")?;
        *self.pair_completed += 1;
        let scored = Scored {
            pair: normalized,
            objective,
        };
        line(
            self.ledger,
            &json!({"kind":"paired_score","ordinal":ordinal,"native":native,"scored":scored,"factors":factors,"background_once":0.1,"per_copy_depth":10,"universe":routes::UNIVERSE,"symbolic_zero_terms":"unmaterialized zero/zero registry/rule terms contribute zero"}),
        )?;
        Ok(scored)
    }
}
