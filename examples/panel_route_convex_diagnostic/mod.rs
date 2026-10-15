//! Standalone bounded finite diagnostic. All profiles originate in public replay.
pub mod math;
use impg::{
    genome_inference::{self as genome, joint, panel_routes as routes, sample},
    syng::SyngIndex,
};
use math::*;
use serde_json::{json, Value};
use std::{collections::BTreeMap, fs, io, path::Path, time::Instant};
pub fn write(path: &Path, value: &impl serde::Serialize) -> io::Result<()> {
    genome::write_json(path, value)
}
fn read(path: &Path) -> io::Result<Value> {
    require(
        fs::metadata(path)?.len() <= 64 * 1024 * 1024,
        "input JSON cap",
    )?;
    genome::read_json(path)
}
fn integer(q: u64) -> io::Result<f64> {
    require(q <= 1u64 << 53, "integer not exactly representable in f64")?;
    Ok(q as f64)
}
/// Independent all-pairs span check, deliberately not the backend implementation.
pub fn canonical_feasible(segments: impl IntoIterator<Item = routes::Segment>) -> bool {
    let spans: Vec<_> = segments.into_iter().collect();
    spans.iter().enumerate().all(|(i, a)| {
        spans[i + 1..]
            .iter()
            .all(|b| a.source != b.source || a.end <= b.start || b.end <= a.start)
    })
}
fn assignment_routes(compiled: &joint::Compiled, choices: &[usize]) -> Vec<routes::Route> {
    let instances: BTreeMap<_, _> = compiled
        .layout
        .instances
        .iter()
        .map(|i| (&i.id, i.source))
        .collect();
    choices
        .iter()
        .enumerate()
        .map(|(s, &a)| {
            routes::Route {
                segments: compiled.layout.slots[s].alternatives[a]
                    .pieces
                    .iter()
                    .map(|p| routes::Segment {
                        source: instances[&p.instance],
                        start: p.start,
                        end: p.end,
                        reverse: p.strand == "-",
                    })
                    .collect(),
            }
            .normalize()
        })
        .collect()
}
fn panel_assignment(
    evaluator: &routes::Evaluator<'_>,
    physical: &[routes::Route],
) -> io::Result<routes::Assignment> {
    let family = evaluator
        .graph
        .families
        .iter()
        .position(|f| {
            f.paths.len() == physical.len()
                && f.paths.iter().zip(physical).all(|(&target, r)| {
                    let first = &r.segments[0];
                    let last = r.segments.last().unwrap();
                    first.source == target
                        && !first.reverse
                        && first.start == 0
                        && last.source == target
                        && !last.reverse
                        && last.end == evaluator.graph.lanes[target].length
                })
        })
        .ok_or_else(|| io::Error::other("no matching native endpoint family"))?;
    Ok(routes::Assignment {
        version: 1,
        model: routes::MODEL.into(),
        graph_checksum: evaluator.graph.digest()?,
        family,
        routes: physical.to_vec(),
    })
}
pub fn analyze(
    compiled: &joint::Compiled,
    sample: &sample::SampleIndex,
    mut route_evaluator: Option<&mut routes::Evaluator<'_>>,
    out: &Path,
) -> io::Result<Value> {
    let start = Instant::now();
    compiled.validate()?;
    let product = compiled.layout.slots.iter().try_fold(1usize, |n, s| {
        n.checked_mul(s.alternatives.len())
            .ok_or_else(|| io::Error::other("product overflow"))
    })?;
    require(
        product <= MAX_PRODUCTS,
        "product cap; enumeration incomplete",
    )?;
    let problem = joint::Problem::new(compiled, sample, DEPTH, BETA)?;
    require(
        problem.factors.len() <= MAX_FEATURES,
        "feature cap; closure incomplete",
    )?;
    let tokens: Vec<_> = problem.factors.iter().map(|f| f.tokens).collect();
    let lookup: BTreeMap<_, _> = tokens.iter().enumerate().map(|(i, &t)| (t, i)).collect();
    let c: Vec<_> = problem
        .factors
        .iter()
        .map(|f| integer(f.observed))
        .collect::<io::Result<_>>()?;
    let histogram: Vec<_> = compiled
        .read_lengths
        .iter()
        .map(|l| {
            sample
                .stats
                .read_lengths
                .get(&(*l as usize))
                .copied()
                .unwrap_or(0)
        })
        .collect();
    let denominator =
        compiled
            .read_lengths
            .iter()
            .zip(&histogram)
            .try_fold(0u64, |s, (&l, &h)| {
                s.checked_add(
                    l.checked_mul(h)
                        .ok_or_else(|| io::Error::other("histogram product overflow"))?,
                )
                .ok_or_else(|| io::Error::other("histogram sum overflow"))
            })?;
    let denominator = integer(denominator)?;
    let background: f64 = c.iter().map(|c| BETA - c * BETA.ln()).sum();
    let mut columns = Vec::new();
    let mut scores = Vec::new();
    let mut rows = Vec::new();
    let mut feasibility = Vec::new();
    let mut count_groups: BTreeMap<Vec<Vec<u64>>, Vec<usize>> = BTreeMap::new();
    let mut generated_zero = false;
    let mut unsupported_positive = false;
    for ordinal in 0..product {
        let mut rest = ordinal;
        let mut choices = vec![0; compiled.layout.slots.len()];
        for s in (0..choices.len()).rev() {
            choices[s] = rest % compiled.layout.slots[s].alternatives.len();
            rest /= compiled.layout.slots[s].alternatives.len();
        }
        let physical = assignment_routes(compiled, &choices);
        let canonical = canonical_feasible(physical.iter().flat_map(|r| r.segments.clone()));
        let feasible = compiled.layout.feasible(&choices)?;
        require(
            canonical == feasible,
            "explicit instance and canonical-span feasibility differ; stop",
        )?;
        let mut route_assignment = None;
        let mut public_validation_error = None;
        if let Some(evaluator) = route_evaluator.as_deref_mut() {
            let a = panel_assignment(evaluator, &physical)?;
            let validation = evaluator.validate_assignment(&a);
            require(
                validation.is_ok() == feasible,
                "public assignment feasibility disagrees with independent enumeration",
            )?;
            public_validation_error = validation.err().map(|e| e.to_string());
            route_assignment = Some(a);
        }
        feasibility.push(json!({"ordinal":ordinal,"choices":choices,"layout_feasible":feasible,"canonical_feasible":canonical,"public_validation_error":public_validation_error}));
        if !feasible {
            continue;
        }
        require(
            columns.len() < MAX_COLUMNS && (columns.len() + 1) * tokens.len() <= MAX_CELLS,
            "column/cell cap; enumeration incomplete",
        )?;
        let retained_numeric_bytes =
            8 * (columns.len() + 1) * tokens.len() * (2 * compiled.read_lengths.len() + 3)
                + 8 * ITERATIONS * (4 * tokens.len() + 2 * (columns.len() + 1));
        require(
            retained_numeric_bytes <= MAX_BYTES,
            "retained numeric state cap; enumeration incomplete",
        )?;
        let mut q = vec![vec![0u64; compiled.read_lengths.len()]; tokens.len()];
        for (s, &a) in choices.iter().enumerate() {
            for contribution in &compiled.profiles[s][a].contributions {
                let row = &mut q[lookup[&contribution.tokens]];
                for (target, totals) in row.iter_mut().zip(&contribution.totals) {
                    require(totals[0] == totals[1], "raw orientation count mismatch")?;
                    *target = target
                        .checked_add(totals[0])
                        .ok_or_else(|| io::Error::other("physical count sum overflow"))?;
                }
            }
        }
        let signal: Vec<f64> = q
            .iter()
            .map(|q| {
                let numerator = q.iter().zip(&histogram).try_fold(0u64, |s, (&q, &h)| {
                    s.checked_add(
                        q.checked_mul(h)
                            .ok_or_else(|| io::Error::other("weighted count product overflow"))?,
                    )
                    .ok_or_else(|| io::Error::other("weighted count sum overflow"))
                })?;
                Ok(DEPTH * integer(numerator)? / denominator)
            })
            .collect::<io::Result<_>>()?;
        let full = joint::evaluate(
            &problem,
            &joint::Assignment {
                version: 1,
                model: joint::MODEL.into(),
                compiled_checksum: compiled.digest()?,
                choices: choices.clone(),
            },
        )?;
        let independent = loss(&c, &signal);
        let joint_relative = full.objective - background;
        close(independent, joint_relative)?;
        for (&s, &mu) in signal.iter().zip(&full.rates) {
            close(BETA + s, mu)?;
        }
        let mut public_score = joint_relative;
        if let Some(evaluator) = route_evaluator.as_deref_mut() {
            let actual = evaluator.evaluate(route_assignment.as_ref().unwrap())?;
            for f in &actual.factors {
                let i = *lookup.get(&f.tokens).ok_or_else(|| {
                    io::Error::other("realized feature absent from complete union")
                })?;
                require(
                    q[i] == f.counts_by_length && integer(f.observed)? == c[i],
                    "public exact counts disagree",
                )?;
                close(signal[i], f.signal)?;
            }
            let actual_tokens: BTreeMap<_, _> =
                actual.factors.iter().map(|f| (f.tokens, f)).collect();
            for (i, t) in tokens.iter().enumerate() {
                require(
                    actual_tokens.contains_key(t) || (c[i] == 0.0 && signal[i] == 0.0),
                    "public factor missing positive or realized signal",
                )?;
            }
            public_score = actual.relative_objective;
            close(independent, public_score)?;
        }
        generated_zero |= c.iter().zip(&signal).any(|(&c, &s)| c == 0.0 && s > 0.0);
        unsupported_positive |= c.iter().zip(&signal).any(|(&c, &s)| c > 0.0 && s == 0.0);
        let index = columns.len();
        count_groups.entry(q.clone()).or_default().push(index);
        rows.push(json!({"column":index,"ordinal":ordinal,"choices":choices,"physical_routes":physical,"panel_assignment":route_assignment,"counts_by_feature_length":q,"signal":signal,"public_relative_objective":public_score,"joint_absolute_objective":full.objective,"absolute_background_shift":background,"independent_relative_objective":independent,"independent_minus_public":independent-public_score,"independent_minus_joint_relative":independent-joint_relative}));
        columns.push(signal);
        scores.push(public_score);
    }
    require(!columns.is_empty(), "exhausted infeasible universe")?;
    let mut identities = Vec::new();
    for (i, t) in columns.iter().enumerate() {
        let s = &columns[0];
        let exact = delta(&c, s, t);
        let difference = scores[i] - scores[0];
        close(exact, difference)?;
        let middle: Vec<_> = s.iter().zip(t).map(|(s, t)| (s + t) / 2.0).collect();
        let direction: Vec<_> = s.iter().zip(t).map(|(s, t)| t - s).collect();
        let plus: Vec<_> = middle
            .iter()
            .zip(&direction)
            .map(|(s, d)| s + 1e-6 * d)
            .collect();
        let minus: Vec<_> = middle
            .iter()
            .zip(&direction)
            .map(|(s, d)| s - 1e-6 * d)
            .collect();
        let finite_difference = (loss(&c, &plus) - loss(&c, &minus)) / 2e-6;
        let derivative = dot(&gradient(&c, &middle), &direction);
        close(derivative, finite_difference)?;
        let tangent = dot(&gradient(&c, s), &direction);
        require(
            tangent <= exact + CHECK * (1.0 + exact.abs()),
            "finite-change convex tangent violated",
        )?;
        identities.push(json!({"from":0,"to":i,"exact_delta":exact,"public_score_difference":difference,"delta_discrepancy":exact-difference,"initial_directional_derivative":tangent,"midpoint_derivative":derivative,"central_difference":finite_difference,"derivative_discrepancy":derivative-finite_difference}));
    }
    let relaxation = optimize(&c, &columns, &scores, true)?;
    let equivalent: Vec<_> = count_groups.values().filter(|g| g.len() > 1).collect();
    write(&out.join("compiled-public.json"), compiled)?;
    write(
        &out.join("features.json"),
        &json!({"tokens":tokens,"observed":c,"read_lengths":compiled.read_lengths,"histogram":histogram,"denominator":denominator,"background_shift":background,"joint_factor_ledger":problem.factors}),
    )?;
    write(&out.join("hypotheses.json"), &rows)?;
    write(&out.join("feasibility.json"), &feasibility)?;
    write(&out.join("identities.json"), &identities)?;
    write(&out.join("relaxation.json"), &relaxation)?;
    let summary = json!({"status":"complete_explicit_finite_universe","products_examined":product,"feasible_columns":columns.len(),"features":tokens.len(),"depth":DEPTH,"background":BETA,"read_lengths":compiled.read_lengths,"histogram":histogram,"public_panel_evaluations":if route_evaluator.is_some(){columns.len()}else{0},"public_joint_evaluations":columns.len(),"generated_zero":generated_zero,"unsupported_positive":unsupported_positive,"exact_count_equivalent_column_groups":equivalent,"equivalence_scope":"integer count vectors only; inspect normalized physical_routes for distinct assignments; score equality does not establish count or DNA equality","correlated_optimum_columns":relaxation.correlated_optimum_column_indices,"best_integral_public_objective":relaxation.discrete_minimum,"elapsed_seconds":start.elapsed().as_secs_f64(),"sequence_emission_authorized":false,"rigorous_real_certificate":false,"genome_wide_bound":false,"gate_b":"NO-GO pending independent review and parent decision"});
    write(&out.join("summary.json"), &summary)?;
    Ok(summary)
}
pub fn run_fixture(root: &Path, out: &Path) -> io::Result<Value> {
    fs::create_dir(out)?;
    write(
        &out.join("status.json"),
        &json!({"status":"running","full_universe_bound":null}),
    )?;
    let result = (|| {
        let prefix = if root.join("panel.syng.meta").exists() {
            root.join("panel.syng")
        } else {
            root.join("p.syng")
        };
        let frozen = read(&root.join("finite/joint-profiles.json"))?;
        let historical: joint::Compiled = serde_json::from_value(frozen["payload"].clone())?;
        require(
            frozen["version"] == 1
                && frozen["model"] == joint::MODEL
                && frozen["checksum"].as_str() == Some(historical.digest()?.as_str()),
            "historical envelope checksum/model mismatch",
        )?;
        let sizes: Vec<_> = historical
            .layout
            .slots
            .iter()
            .map(|s| s.alternatives.len())
            .collect();
        let source_lengths: Vec<_> = historical.layout.sources.iter().map(|s| s.length).collect();
        require(
            (sizes == [513]
                && source_lengths == [400, 400, 120]
                && historical.read_lengths == [150, 500])
                || (sizes == [64, 64]
                    && source_lengths == [220, 220]
                    && historical.read_lengths == [220]),
            "standalone input is not a declared finite513/coupled64x64 synthetic fixture",
        )?;
        let identity = genome::PanelIdentity::read(prefix.to_str().unwrap())?;
        let panel = SyngIndex::load(prefix.to_str().unwrap(), Default::default())?;
        let graph = routes::Graph::load(&root.join("routes"), &identity)?;
        // Do not spoof historical compiler identity: recompile via current public
        // operator, and compare complete integer profiles to the frozen oracle.
        let layout = historical.layout.clone();
        let lengths = historical.read_lengths.clone();
        let catalog = root.join("catalog/catalog.json");
        let compiled = joint::compile(
            &panel,
            identity.clone(),
            layout,
            &graph.source_paths,
            lengths,
            Some(&catalog),
        )?;
        let previous = historical.profiles;
        require(
            previous.len() == compiled.profiles.len(),
            "frozen slot count mismatch",
        )?;
        for (a, b) in previous.iter().zip(&compiled.profiles) {
            require(a.len() == b.len(), "frozen alternative count mismatch")?;
            for (a, b) in a.iter().zip(b) {
                require(
                    serde_json::to_value(&a.contributions)?
                        == serde_json::to_value(&b.contributions)?
                        && a.admitted_starts == b.admitted_starts,
                    "frozen public integer profile mismatch",
                )?;
            }
        }
        write(
            &out.join("frozen-profile-comparison.json"),
            &json!({"all_profiles_exactly_equal":true,"historical_compiler":frozen["payload"]["compiler_identity"],"current_compiler":compiled.compiler_identity,"frozen_envelope_checksum":frozen["checksum"],"input":root}),
        )?;
        let sample = sample::SampleIndex::load(&root.join("sample.membwt"), &identity)?;
        let mut evaluator = routes::Evaluator::new(
            &root.join("routes"),
            &graph,
            &panel,
            &sample,
            DEPTH,
            BETA,
            MAX_FEATURES,
            100_000,
        )?;
        analyze(&compiled, &sample, Some(&mut evaluator), out)
    })();
    write(
        &out.join("status.json"),
        &match &result {
            Ok(_) => {
                json!({"status":"succeeded","numerical_scope":"exhausted explicit finite layout only"})
            }
            Err(e) => {
                json!({"status":"failed_or_incomplete","error":e.to_string(),"full_universe_bound":null,"correlated_support_complete":false})
            }
        },
    )?;
    result
}
