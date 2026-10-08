//! Complete-assignment objective and deliberately bounded exhaustive search.
//! Components are diagnostic, not a claim of scalable factorized optimization.
use super::*;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Assignment {
    pub version: u32,
    pub model: String,
    pub compiled_checksum: String,
    /// Ordered by layout.slots; indexes refer to the frozen compiled alternatives.
    pub choices: Vec<usize>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Term {
    pub slot: usize,
    pub alternative: usize,
    pub exposure: f64,
}
#[derive(Clone, Debug, Serialize)]
pub struct Factor {
    pub tokens: [u64; 3],
    pub original_ids: Vec<usize>,
    pub observed: u64,
    pub terms: Vec<Term>,
    pub scope: Vec<usize>,
    pub constant: bool,
    pub background_only: bool,
}
pub struct Problem<'a> {
    pub compiled: &'a Compiled,
    pub factors: Vec<Factor>,
    pub components: Vec<Vec<usize>>,
    pub depth: f64,
    pub background: f64,
    rates: Vec<Vec<Vec<(usize, f64)>>>,
}
impl<'a> Problem<'a> {
    pub fn new(
        compiled: &'a Compiled,
        sample: &sample::SampleIndex,
        depth: f64,
        background: f64,
    ) -> io::Result<Self> {
        compiled.validate()?;
        require(
            sample.version == super::super::FORMAT_VERSION
                && sample.panel == compiled.layout.panel
                && sample.count_policy == COUNT_POLICY,
            "joint sample identity/count policy mismatch",
        )?;
        require(
            depth.is_finite() && depth > 0.0 && background.is_finite() && background > 0.0,
            "depth and background must be finite positive",
        )?;
        require(
            !sample.stats.read_lengths.is_empty()
                && sample
                    .stats
                    .read_lengths
                    .iter()
                    .all(|(&l, &n)| n > 0 && compiled.read_lengths.contains(&(l as u64))),
            "sample contains empty/uncompiled read histogram",
        )?;
        let denominator: f64 = sample
            .stats
            .read_lengths
            .iter()
            .map(|(&l, &n)| l as f64 * n as f64)
            .sum();
        require(
            denominator.is_finite() && denominator > 0.0,
            "invalid read histogram denominator",
        )?;
        // Enumeration must succeed in full before a problem can be scored. Sample
        // only features remain explicit residual factors with mu=background.
        let observed = sample.counts.observed_pairs()?;
        let mut rows: BTreeMap<_, _> = compiled
            .definitions
            .iter()
            .map(|d| {
                (
                    d.tokens,
                    Factor {
                        tokens: d.tokens,
                        original_ids: d.original_ids.clone(),
                        observed: 0,
                        terms: Vec::new(),
                        scope: Vec::new(),
                        constant: false,
                        background_only: false,
                    },
                )
            })
            .collect();
        for (&tokens, &count) in &observed {
            let f = rows.entry(tokens).or_insert(Factor {
                tokens,
                original_ids: Vec::new(),
                observed: 0,
                terms: Vec::new(),
                scope: Vec::new(),
                constant: false,
                background_only: false,
            });
            f.observed = count;
        }
        let mut factors: Vec<_> = rows.into_values().collect();
        let lookup: BTreeMap<_, _> = factors
            .iter()
            .enumerate()
            .map(|(i, f)| (f.tokens, i))
            .collect();
        let mut rates = Vec::new();
        for (s, profiles) in compiled.profiles.iter().enumerate() {
            let mut alternatives = Vec::new();
            for (a, p) in profiles.iter().enumerate() {
                let mut terms = Vec::new();
                for c in &p.contributions {
                    let numerator: f64 = compiled
                        .read_lengths
                        .iter()
                        .zip(&c.totals)
                        .map(|(&l, q)| {
                            sample
                                .stats
                                .read_lengths
                                .get(&(l as usize))
                                .copied()
                                .unwrap_or(0) as f64
                                * (q[0] as f64 / 2.0 + q[1] as f64 / 2.0)
                        })
                        .sum();
                    let exposure = numerator / denominator;
                    require(
                        exposure.is_finite() && (depth * exposure).is_finite(),
                        "nonfinite physical exposure",
                    )?;
                    if exposure > 0.0 {
                        let f = lookup[&c.tokens];
                        factors[f].terms.push(Term {
                            slot: s,
                            alternative: a,
                            exposure,
                        });
                        terms.push((f, exposure));
                    }
                }
                alternatives.push(terms);
            }
            rates.push(alternatives);
        }
        for f in &mut factors {
            f.scope = f
                .terms
                .iter()
                .map(|t| t.slot)
                .collect::<BTreeSet<_>>()
                .into_iter()
                .collect();
            f.background_only = f.terms.is_empty();
            f.constant = f.scope.iter().all(|&s| {
                let values: Vec<_> = (0..compiled.layout.slots[s].alternatives.len())
                    .map(|a| {
                        f.terms
                            .iter()
                            .find(|t| t.slot == s && t.alternative == a)
                            .map_or(0.0, |t| t.exposure)
                    })
                    .collect();
                values.iter().all(|&v| v == values[0])
            });
        }
        // Union shared factors AND physical-resource conflicts. No independent
        // solve uses components yet; this keeps the future factorization seam honest.
        let n = compiled.layout.slots.len();
        let mut labels: Vec<_> = (0..n).collect();
        let mut connect = |scope: &[usize]| {
            if let Some(&first) = scope.first() {
                let label = labels[first];
                for &s in scope.iter().skip(1) {
                    let old = labels[s];
                    for l in &mut labels {
                        if *l == old {
                            *l = label;
                        }
                    }
                }
            }
        };
        for f in &factors {
            connect(&f.scope);
        }
        let mut resource_slots: BTreeMap<String, BTreeSet<usize>> = BTreeMap::new();
        for (s, slot) in compiled.layout.slots.iter().enumerate() {
            for a in 0..slot.alternatives.len() {
                for resource in compiled.layout.resources(s, a) {
                    resource_slots.entry(resource).or_default().insert(s);
                }
            }
        }
        for scope in resource_slots.values() {
            connect(&scope.iter().copied().collect::<Vec<_>>());
        }
        let mut components: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
        for (s, label) in labels.into_iter().enumerate() {
            components.entry(label).or_default().push(s);
        }
        Ok(Self {
            compiled,
            factors,
            components: components.into_values().collect(),
            depth,
            background,
            rates,
        })
    }
    fn score(&self, choices: &[usize]) -> io::Result<(f64, Vec<f64>, Vec<f64>)> {
        require(
            self.compiled.layout.feasible(choices)?,
            "infeasible simultaneous physical-instance reuse",
        )?;
        let mut mus = vec![self.background; self.factors.len()];
        for (s, &a) in choices.iter().enumerate() {
            for &(f, e) in &self.rates[s][a] {
                mus[f] += self.depth * e;
            }
        }
        let losses: Vec<_> = self
            .factors
            .iter()
            .zip(&mus)
            .map(|(f, &mu)| mu - f.observed as f64 * mu.ln())
            .collect();
        let objective: f64 = losses.iter().sum();
        require(
            objective.is_finite() && mus.iter().all(|m| m.is_finite()),
            "nonfinite joint objective",
        )?;
        Ok((objective, mus, losses))
    }
}
#[derive(Clone, Debug, Serialize)]
pub struct Evaluation {
    pub version: u32,
    pub model: String,
    pub choices: Vec<usize>,
    pub objective: f64,
    pub rates: Vec<f64>,
    pub losses: Vec<f64>,
    pub complete_assignment: bool,
    pub context_complete: bool,
    pub candidate_universe: String,
    pub sequence_emission_authorized: bool,
}
pub fn evaluate(problem: &Problem<'_>, assignment: &Assignment) -> io::Result<Evaluation> {
    require(
        assignment.version == VERSION
            && assignment.model == MODEL
            && assignment.compiled_checksum == problem.compiled.digest()?,
        "incompatible joint assignment/provenance",
    )?;
    let (objective, rates, losses) = problem.score(&assignment.choices)?;
    Ok(Evaluation {
        version: VERSION,
        model: MODEL.into(),
        choices: assignment.choices.clone(),
        objective,
        rates,
        losses,
        complete_assignment: true,
        context_complete: true,
        candidate_universe: "finite-explicit-layout-only".into(),
        sequence_emission_authorized: false,
    })
}
#[derive(Clone, Debug, Serialize)]
pub struct OptimalAssignment {
    pub choices: Vec<usize>,
    pub objective: f64,
}
#[derive(Debug, Serialize)]
pub struct SearchResult {
    pub version: u32,
    pub model: String,
    pub method: String,
    pub status: String,
    pub assignments_examined: u64,
    pub feasible_assignments: u64,
    pub max_assignments: u64,
    pub max_optima: usize,
    pub tie_epsilon: f64,
    pub incumbent: Option<OptimalAssignment>,
    pub lower_bound: Option<f64>,
    pub search_exhausted: bool,
    pub global_objective_certified: bool,
    pub correlated_optima_complete: bool,
    pub correlated_optima: Vec<OptimalAssignment>,
    pub components: Vec<Vec<usize>>,
    pub context_complete: bool,
    pub candidate_universe_complete: bool,
    pub sequence_emission_authorized: bool,
}
/// Exhaustive lexicographic assignments, including infeasible assignments in the
/// budget. No unary pruning or floating branch bounds. Certification is within
/// this explicit universe and the stated f64/tie-epsilon objective, not a posterior.
pub fn solve(
    problem: &Problem<'_>,
    max_assignments: u64,
    max_optima: usize,
    epsilon: f64,
) -> io::Result<SearchResult> {
    require(
        max_assignments > 0 && max_optima > 0 && epsilon.is_finite() && epsilon >= 0.0,
        "invalid exact-search budget/tolerance",
    )?;
    let mut choices = vec![0; problem.compiled.layout.slots.len()];
    let mut examined = 0;
    let mut feasible = 0;
    let mut exhausted = false;
    let mut incumbent: Option<OptimalAssignment> = None;
    let mut optima: Vec<OptimalAssignment> = Vec::new();
    let mut discarded = false;
    while examined < max_assignments {
        examined += 1;
        if problem.compiled.layout.feasible(&choices)? {
            feasible += 1;
            let (objective, _, _) = problem.score(&choices)?;
            if incumbent.as_ref().is_none_or(|b| objective < b.objective) {
                incumbent = Some(OptimalAssignment {
                    choices: choices.clone(),
                    objective,
                });
                optima.retain(|o| o.objective - objective <= epsilon);
            }
            if objective - incumbent.as_ref().unwrap().objective <= epsilon {
                if optima.len() < max_optima {
                    optima.push(OptimalAssignment {
                        choices: choices.clone(),
                        objective,
                    });
                } else {
                    discarded = true;
                }
            }
        }
        let mut s = choices.len();
        loop {
            if s == 0 {
                exhausted = true;
                break;
            }
            s -= 1;
            choices[s] += 1;
            if choices[s] < problem.compiled.layout.slots[s].alternatives.len() {
                break;
            }
            choices[s] = 0;
        }
        if exhausted {
            break;
        }
    }
    // If storage was ever exhausted, conservatively withhold complete support,
    // even if a later lower incumbent might have made discarded ties irrelevant.
    let complete = exhausted && !discarded;
    let status = if !exhausted {
        "budget-exhausted"
    } else if incumbent.is_none() {
        "infeasible"
    } else if !complete {
        "optimum-certified-alternatives-incomplete"
    } else {
        "optimum-certified"
    };
    Ok(SearchResult {
        version: VERSION,
        model: MODEL.into(),
        method: "bounded-exhaustive-f64-v1".into(),
        status: status.into(),
        assignments_examined: examined,
        feasible_assignments: feasible,
        max_assignments,
        max_optima,
        tie_epsilon: epsilon,
        lower_bound: if exhausted {
            incumbent.as_ref().map(|i| i.objective)
        } else {
            None
        },
        incumbent,
        search_exhausted: exhausted,
        global_objective_certified: exhausted,
        correlated_optima_complete: complete,
        correlated_optima: optima,
        components: problem.components.clone(),
        context_complete: true,
        candidate_universe_complete: false,
        sequence_emission_authorized: false,
    })
}
