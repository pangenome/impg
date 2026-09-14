//! New finite copy model. The production assignment/capacity policy is unchanged.
use impg::{
    genome_inference::{self as genome, panel_routes as routes, sample},
    syng::SyngIndex,
};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::{self, File},
    io::{self, BufWriter, Write},
    path::Path,
    time::Instant,
};
pub const MAX_FEATURES: usize = 100_000;
pub const MAX_CELLS: usize = 2_000_000;
pub const MAX_BYTES: usize = 128 * 1024 * 1024;
pub const TIE: f64 = 1e-9;
pub type Counts = BTreeMap<[u64; 3], u64>;
pub fn require(ok: bool, text: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(io::Error::other(text))
    }
}
pub fn add(a: u64, b: u64) -> io::Result<u64> {
    a.checked_add(b)
        .ok_or_else(|| io::Error::other("integer count overflow"))
}
fn exact(x: u64) -> io::Result<f64> {
    require(x <= 1 << 53, "integer exceeds exact f64 range")?;
    Ok(x as f64)
}
pub fn close(a: f64, b: f64) -> io::Result<()> {
    require(
        a.is_finite() && b.is_finite() && (a - b).abs() <= 1e-8 * (1.0 + a.abs()),
        "independent formula mismatch",
    )
}
/// Exactly two whole-hypothesis resources; there is no caller-supplied copy ID.
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct FixedDiploid {
    pub copies: [routes::Assignment; 2],
}
impl FixedDiploid {
    pub fn validate(&self, e: &mut routes::Evaluator<'_>) -> io::Result<Self> {
        let mut copies = [
            e.validate_assignment(&self.copies[0])?,
            e.validate_assignment(&self.copies[1])?,
        ];
        copies.sort(); // Exchange labels, never remove the second physical copy.
        Ok(Self { copies })
    }
}
/// Conservative logical charge, not allocator/RSS or public-operator scratch.
/// 256 bytes/cell covers sparse tree entries, dense ledgers and serialization
/// buffers; a fixed 1MiB covers bounded assignment/route/row metadata.
pub fn state_charge(cells: usize) -> io::Result<usize> {
    require(cells <= MAX_CELLS, "retained count-cell cap")?;
    let bytes = cells
        .checked_mul(256)
        .and_then(|x| x.checked_add(1024 * 1024))
        .ok_or_else(|| io::Error::other("state charge overflow"))?;
    require(bytes <= MAX_BYTES, "retained state cap")?;
    Ok(bytes)
}
pub fn joint_loss(q: &[u64], c: &[u64], histogram: u64) -> io::Result<f64> {
    require(
        q.len() == c.len() && histogram > 0,
        "invalid factors/histogram",
    )?;
    let denominator = exact(
        histogram
            .checked_mul(150)
            .ok_or_else(|| io::Error::other("histogram overflow"))?,
    )?;
    let mut relative = 0.0;
    let mut independent = 0.0;
    for (&q, &c) in q.iter().zip(c) {
        let numerator = exact(
            q.checked_mul(histogram)
                .ok_or_else(|| io::Error::other("weighted count overflow"))?,
        )?;
        let s = 10.0 * numerator / denominator;
        let c = exact(c)?;
        relative += s - c * (s / 0.1).ln_1p();
        independent += (0.1 + s) - c * (0.1 + s).ln() - (0.1 - c * 0.1f64.ln());
    }
    close(relative, independent)?;
    Ok(relative)
}
#[derive(Serialize, Deserialize)]
pub struct PairRow {
    pub pair: [usize; 2],
    pub counts_by_length_150: Vec<u64>,
    pub relative_objective: f64,
    pub realized_zero_features: usize,
    pub unsupported_positive_features: usize,
    pub shared_signal_features: usize,
}
/// Exhaustive read replay uses the public sample BWT, not route event/cache code.
/// Each window is an anonymous singleton record; full reverse replay checks RC orbits.
pub fn replay(
    e: &mut routes::Evaluator<'_>,
    a: &routes::Assignment,
    out: &Path,
    ordinal: usize,
    replay_reads: &mut usize,
) -> io::Result<Counts> {
    let mut expected = Counts::new();
    for r in &a.routes {
        let (q, _) = e.route_counts(r)?;
        for (&t, &n) in &q[0] {
            let old = expected.get(&t).copied().unwrap_or(0);
            expected.insert(t, add(old, n)?);
        }
    }
    for reverse in [false, true] {
        let path = out.join(format!("replay-{ordinal}-{reverse}.fa"));
        let mut f = BufWriter::new(File::create(&path)?);
        let before = *replay_reads;
        for r in &a.routes {
            let length = r.length()?;
            let n = length.saturating_sub(149) as usize;
            require(*replay_reads + n <= 100_000, "all-window replay read cap")?;
            *replay_reads += n;
            let mut dna = e.sources.spell(r, 0, length)?;
            if reverse {
                dna = impg::graph::reverse_complement(&dna);
            }
            for window in dna.windows(150) {
                writeln!(
                    f,
                    ">anonymous\n{}",
                    std::str::from_utf8(window).map_err(io::Error::other)?
                )?;
            }
        }
        f.flush()?;
        if *replay_reads == before {
            // Empty FASTA is rejected by the public reader. No admitted physical
            // windows means an empty count vector, not an invented sentinel read.
            require(expected.is_empty(), "zero-start route has nonzero profile")?;
            continue;
        }
        let s = sample::build(e.panel, e.graph.panel.clone(), &[path])?;
        require(
            s.counts.observed_pairs()? == expected,
            "all-window public BWT / complete-route integer mismatch",
        )?;
    }
    Ok(expected)
}
pub fn analyze(
    e: &mut routes::Evaluator<'_>,
    sample: &sample::SampleIndex,
    supplied: Vec<routes::Assignment>,
    out: &Path,
) -> io::Result<Value> {
    let start = Instant::now();
    require(
        e.graph.read_lengths == [150] && e.depth == 10.0 && e.background == 0.1,
        "frozen diagnostic requires L150/per-copy depth10/background0.1",
    )?;
    require(
        !supplied.is_empty() && supplied.len() <= 8,
        "haploid domain cap/empty domain",
    )?;
    let mut domain = BTreeSet::new();
    for a in supplied {
        require(
            a.routes.len() <= 4 && a.routes.iter().map(|r| r.segments.len()).sum::<usize>() <= 32,
            "molecule/segment cap",
        )?;
        for r in &a.routes {
            require(r.length()? <= 4096, "molecule length cap")?;
        }
        domain.insert(e.validate_assignment(&a)?);
    }
    let domain: Vec<_> = domain.into_iter().collect();
    let pairs = domain.len() * (domain.len() + 1) / 2;
    require(pairs <= 36, "paired domain cap")?;
    genome::write_json(&out.join("normalized-domain.json"), &domain)?;
    let observed = sample.counts.observed_pairs()?;
    require(observed.len() <= MAX_FEATURES, "observed feature cap")?;
    let mut keys: BTreeSet<_> = observed.keys().copied().collect();
    let mut profiles = Vec::new();
    let mut replay_reads = 0;
    let mut cells = observed.len();
    let mut peak_charge = state_charge(cells)?;
    let mut haploid_losses = Vec::new();
    for (i, a) in domain.iter().enumerate() {
        // Admit a worst-case next public profile before allocating it. Operator
        // evaluation and replay internals have their own max_terms/scratch scope.
        peak_charge = peak_charge.max(state_charge(cells + MAX_FEATURES)?);
        let ev = e.evaluate(a)?;
        let expected = replay(e, a, out, i, &mut replay_reads)?;
        let mut profile = Counts::new();
        for f in &ev.factors {
            require(
                f.observed == observed.get(&f.tokens).copied().unwrap_or(0),
                "observed count altered",
            )?;
            require(
                f.counts_by_length.len() == 1
                    && f.counts_by_length[0] == expected.get(&f.tokens).copied().unwrap_or(0),
                "public profile/replay mismatch",
            )?;
            if f.counts_by_length[0] > 0 {
                profile.insert(f.tokens, f.counts_by_length[0]);
            }
        }
        require(
            profile == expected,
            "complete public profile missing replay term",
        )?;
        keys.extend(profile.keys().copied());
        require(keys.len() <= MAX_FEATURES, "finite feature union cap")?;
        cells += profile.len();
        peak_charge = peak_charge.max(state_charge(cells)?);
        profiles.push(profile);
        haploid_losses.push(ev.relative_objective);
    }
    let tokens: Vec<_> = keys.into_iter().collect();
    let c: Vec<_> = tokens
        .iter()
        .map(|t| observed.get(t).copied().unwrap_or(0))
        .collect();
    // Includes all retained pair rows, working vectors, token/observed arrays and
    // an independent comparison profile, before their allocation.
    let charged_cells = cells + (pairs + 8) * tokens.len();
    peak_charge = peak_charge.max(state_charge(charged_cells)?);
    genome::write_json(
        &out.join("features.json"),
        &json!({"tokens":tokens,"observed_once":c,"universe":routes::UNIVERSE,"registry_definitions":e.graph.registry_count,"symbolic_zero_terms":"unmaterialized registry/rule zero-count zero-signal contributes zero"}),
    )?;
    let h = *sample
        .stats
        .read_lengths
        .get(&150)
        .ok_or_else(|| io::Error::other("missing histogram"))?;
    let mut rows = Vec::new();
    let cache_before = e.cache_hits;
    let mut aa_double_checks = 0;
    for i in 0..domain.len() {
        for j in i..domain.len() {
            let pair = FixedDiploid {
                copies: [domain[i].clone(), domain[j].clone()],
            }
            .validate(e)?;
            let mut independent = Counts::new();
            for copy in &pair.copies {
                for route in &copy.routes {
                    let (q, _) = e.route_counts(route)?;
                    for (&t, &v) in &q[0] {
                        let old = independent.get(&t).copied().unwrap_or(0);
                        independent.insert(t, add(old, v)?);
                    }
                }
            }
            let mut q = Vec::with_capacity(tokens.len());
            let mut shared = 0;
            for t in &tokens {
                let a = profiles[i].get(t).copied().unwrap_or(0);
                let b = profiles[j].get(t).copied().unwrap_or(0);
                let sum = add(a, b)?;
                require(
                    sum == independent.get(t).copied().unwrap_or(0),
                    "copy profile sum/cache multiplicity mismatch",
                )?;
                if a > 0 && b > 0 {
                    shared += 1;
                }
                if i == j {
                    require(
                        sum == a
                            .checked_mul(2)
                            .ok_or_else(|| io::Error::other("dosage overflow"))?,
                        "A/A dosage mismatch",
                    )?;
                }
                q.push(sum);
            }
            if i == j {
                aa_double_checks += 1;
            }
            let loss = joint_loss(&q, &c, h)?;
            rows.push(PairRow {
                pair: [i, j],
                relative_objective: loss,
                realized_zero_features: q
                    .iter()
                    .zip(&c)
                    .filter(|(q, c)| **q > 0 && **c == 0)
                    .count(),
                unsupported_positive_features: q
                    .iter()
                    .zip(&c)
                    .filter(|(q, c)| **q == 0 && **c > 0)
                    .count(),
                shared_signal_features: shared,
                counts_by_length_150: q,
            });
        }
    }
    let minimum = rows
        .iter()
        .map(|r| r.relative_objective)
        .fold(f64::INFINITY, f64::min);
    let support: Vec<_> = rows
        .iter()
        .filter(|r| (r.relative_objective - minimum).abs() <= TIE)
        .map(|r| r.pair)
        .collect();
    let mut equivalent = BTreeMap::<&Vec<u64>, Vec<[usize; 2]>>::new();
    for row in &rows {
        equivalent
            .entry(&row.counts_by_length_150)
            .or_default()
            .push(row.pair);
    }
    let groups: Vec<_> = equivalent.values().filter(|g| g.len() > 1).collect();
    let mut old = domain[0].clone();
    old.routes.extend(domain[0].routes.clone());
    let old_error = e.validate_assignment(&old).err().map(|x| x.to_string());
    require(
        old_error.is_some(),
        "old concatenated assignment unexpectedly valid",
    )?;
    // The selection artifact is finalized before any fixture truth assessment.
    genome::write_json(&out.join("pairs.json"), &rows)?;
    genome::write_json(
        &out.join("selection.json"),
        &json!({"minimum":minimum,"correlated_pairs":support,"absolute_tie_epsilon":TIE,"finite_domain_complete":true,"support_complete_within_supplied_domain":true,"automatic_generation":false}),
    )?;
    let summary = json!({"model":"fixed-two-whole-haploid-resources-finite-v1","haploid_assignments":domain.len(),"paired_hypotheses":pairs,"minimum":minimum,"correlated_optima":support,"exact_count_equivalent_pair_groups":groups,"aa_exact_double_checks":aa_double_checks,"pair_route_cache_hits":e.cache_hits-cache_before,"old_concatenation_actual_rejection":old_error,"haploid_losses_not_added":haploid_losses,"observed_features":observed.len(),"realized_union_features":tokens.len(),"replayed_reads":replay_reads,"retained_count_cells_charged":charged_cells,"peak_logical_state_charge_bytes":peak_charge,"state_scope":"conservative diagnostic numerical/metadata charge; public evaluator/replay scratch and measured RSS separate","nominal_per_copy_depth":10,"nominal_haploid_reference_depth":20,"background_once":0.1,"read_length":150,"actual_reads":sample.stats.reads,"actual_bases":sample.stats.bases,"elapsed_seconds":start.elapsed().as_secs_f64(),"representation_operator_checks_passed":true,"general_diploid_support":false,"sequence_emission_authorized":false});
    genome::write_json(&out.join("summary.json"), &summary)?;
    Ok(summary)
}
pub fn run(
    panel: &Path,
    graph_root: &Path,
    sample_path: &Path,
    domain_path: &Path,
    out: &Path,
) -> io::Result<Value> {
    fs::create_dir(out)?;
    let result = (|| {
        require(
            fs::metadata(domain_path)?.len() <= 1024 * 1024,
            "domain JSON input cap",
        )?;
        let identity = genome::PanelIdentity::read(
            panel
                .to_str()
                .ok_or_else(|| io::Error::other("non UTF8 panel"))?,
        )?;
        let graph = routes::Graph::load(graph_root, &identity)?;
        let panel = SyngIndex::load(panel.to_str().unwrap(), Default::default())?;
        let sample = sample::SampleIndex::load(sample_path, &identity)?;
        let domain = genome::read_json(domain_path)?;
        let mut e = routes::Evaluator::new(
            graph_root,
            &graph,
            &panel,
            &sample,
            10.0,
            0.1,
            MAX_FEATURES,
            MAX_FEATURES,
        )?;
        analyze(&mut e, &sample, domain, out)
    })();
    if let Err(error) = &result {
        genome::write_json(
            &out.join("failure.json"),
            &json!({"error":error.to_string(),"support_complete":false}),
        )?;
    }
    result
}
