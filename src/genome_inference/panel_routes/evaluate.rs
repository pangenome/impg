use super::*;
use profiles::Counts;
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct Segment {
    pub source: usize,
    pub start: u64,
    pub end: u64,
    pub reverse: bool,
}
impl Segment {
    pub fn entry(&self) -> u64 {
        if self.reverse {
            self.end
        } else {
            self.start
        }
    }
    pub fn exit(&self) -> u64 {
        if self.reverse {
            self.start
        } else {
            self.end
        }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct Route {
    pub segments: Vec<Segment>,
}
impl Route {
    pub fn length(&self) -> io::Result<u64> {
        self.segments.iter().try_fold(0u64, |n, p| {
            n.checked_add(
                p.end
                    .checked_sub(p.start)
                    .ok_or_else(|| invalid("reversed physical interval"))?,
            )
            .ok_or_else(|| invalid("route length overflow"))
        })
    }
    pub fn normalize(&self) -> Self {
        let mut out: Vec<Segment> = Vec::new();
        for p in &self.segments {
            if let Some(last) = out.last_mut() {
                if last.source == p.source && last.reverse == p.reverse && last.exit() == p.entry()
                {
                    last.start = last.start.min(p.start);
                    last.end = last.end.max(p.end);
                    continue;
                }
            }
            out.push(p.clone());
        }
        Self { segments: out }
    }
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq, PartialOrd, Ord)]
#[serde(deny_unknown_fields)]
pub struct Assignment {
    pub version: u32,
    pub model: String,
    pub graph_checksum: String,
    pub family: usize,
    pub routes: Vec<Route>,
}
#[derive(Clone, Debug, Serialize)]
pub struct Factor {
    pub tokens: [u64; 3],
    pub counts_by_length: Vec<u64>,
    pub observed: u64,
    pub signal: f64,
    pub relative_loss: f64,
    pub positive_residual: bool,
}
#[derive(Debug, Serialize)]
pub struct Evaluation {
    pub version: u32,
    pub model: String,
    pub assignment: Assignment,
    pub objective_kind: String,
    pub relative_objective: f64,
    pub universe: String,
    pub factors: Vec<Factor>,
    pub sample_positive_features: usize,
    pub registry_definitions: u64,
    pub symbolic_zero_terms: String,
    pub admitted_starts: Vec<u64>,
    pub corrected_starts: u64,
    pub native_range_starts: u64,
    pub route_cache_hits: u64,
    pub native_evaluation_complete: bool,
    pub generation_rule_complete: bool,
    pub topology_inventory_complete_for_panel: bool,
    pub biological_topology_complete: bool,
    pub sequence_emission_authorized: bool,
}
pub struct Evaluator<'a> {
    pub graph: &'a Graph,
    pub panel: &'a SyngIndex,
    pub ports: graph::Ports,
    pub sources: graph::Sources,
    store: profiles::Store,
    observed: Counts,
    histogram: Vec<u64>,
    denominator: f64,
    pub depth: f64,
    pub background: f64,
    pub max_terms: usize,
    cache: BTreeMap<Route, Vec<Counts>>,
    cache_terms: usize,
    max_cache_terms: usize,
    pub cache_hits: u64,
    pub(super) graph_checksum: String,
}
/// Canonical source-forward span capacity, shared by both orientations and all
/// storage/ownership descriptions. Disjoint uses are admitted even across slots.
pub fn resources_feasible(segments: impl IntoIterator<Item = Segment>) -> bool {
    let mut spans: BTreeMap<usize, Vec<(u64, u64)>> = BTreeMap::new();
    for p in segments {
        spans.entry(p.source).or_default().push((p.start, p.end));
    }
    spans.values_mut().all(|s| {
        s.sort_unstable();
        s.windows(2).all(|w| w[0].1 <= w[1].0)
    })
}
impl<'a> Evaluator<'a> {
    pub fn new(
        root: &Path,
        graph: &'a Graph,
        panel: &'a SyngIndex,
        sample: &super::super::sample::SampleIndex,
        depth: f64,
        background: f64,
        max_terms: usize,
        max_cache_terms: usize,
    ) -> io::Result<Self> {
        graph.validate()?;
        require(
            sample.version == 1
                && sample.panel == graph.panel
                && sample.count_policy == COUNT_POLICY,
            "sample/route operator or dictionary mismatch",
        )?;
        require(
            depth.is_finite()
                && depth > 0.0
                && background.is_finite()
                && background > 0.0
                && max_terms > 0,
            "invalid route scoring parameters",
        )?;
        require(
            !sample.stats.read_lengths.is_empty()
                && sample
                    .stats
                    .read_lengths
                    .iter()
                    .all(|(&l, &n)| n > 0 && graph.read_lengths.contains(&(l as u64))),
            "unsupported/incomplete route read histogram",
        )?;
        require(
            panel.syncmer_length_bp() as u64 == graph.k
                && panel.name_map.path_to_name.len() == graph.lanes.len(),
            "native panel shape mismatch",
        )?;
        for lane in &graph.lanes {
            require(
                panel.name_map.path_to_name[lane.id] == lane.name
                    && panel.name_map.path_to_length[lane.id] == lane.length,
                "native panel lane mismatch",
            )?;
        }
        for (p, f) in graph.source_paths.iter().zip(&graph.source_files) {
            require(
                *f == super::super::reconstruction::fingerprint(Path::new(p))?,
                "route source provenance mismatch",
            )?;
        }
        // Verify auxiliary lookup files before any native reader can recreate or
        // consume them. Sequence bytes alone do not bind FAI/GZI address translation.
        graph::verify_source_access(&graph.source_paths, &graph.source_access)?;
        let sources = graph::Sources::open(
            &graph.source_paths,
            graph
                .lanes
                .iter()
                .map(|l| (l.name.clone(), l.length))
                .collect(),
        )?;
        let histogram: Vec<_> = graph
            .read_lengths
            .iter()
            .map(|&l| {
                sample
                    .stats
                    .read_lengths
                    .get(&(l as usize))
                    .copied()
                    .unwrap_or(0)
            })
            .collect();
        let denominator: f64 = graph
            .read_lengths
            .iter()
            .zip(&histogram)
            .map(|(&l, &n)| l as f64 * n as f64)
            .sum();
        require(
            denominator.is_finite() && denominator > 0.0,
            "invalid histogram normalization",
        )?;
        let observed = sample.counts.observed_pairs()?;
        require(
            observed.len() <= max_terms,
            "sample support resource cap exhausted; no partial scoring",
        )?;
        Ok(Self {
            graph,
            panel,
            ports: graph::Ports::open(root, graph)?,
            sources,
            store: profiles::Store::new(root),
            observed,
            histogram,
            denominator,
            depth,
            background,
            max_terms,
            cache: BTreeMap::new(),
            cache_terms: 0,
            max_cache_terms,
            cache_hits: 0,
            graph_checksum: graph.digest()?,
        })
    }
    pub fn validate_assignment(&mut self, a: &Assignment) -> io::Result<Assignment> {
        require(
            a.version == VERSION
                && a.model == MODEL
                && a.graph_checksum == self.graph_checksum
                && a.family < self.graph.families.len(),
            "route assignment identity mismatch",
        )?;
        let family = &self.graph.families[a.family];
        require(
            a.routes.len() == family.paths.len(),
            "incomplete topology assignment",
        )?;
        let mut normalized = a.clone();
        for (route, &target) in normalized.routes.iter_mut().zip(&family.paths) {
            require(!route.segments.is_empty(), "empty route")?;
            for p in &route.segments {
                require(
                    p.source < self.graph.lanes.len()
                        && p.start < p.end
                        && p.end <= self.graph.lanes[p.source].length,
                    "invalid/zero-length physical traversal",
                )?;
            }
            *route = route.normalize();
            route.length()?;
            let first = &route.segments[0];
            let last = route.segments.last().unwrap();
            require(
                first.source == target
                    && !first.reverse
                    && first.start == 0
                    && last.source == target
                    && !last.reverse
                    && last.end == self.graph.lanes[target].length,
                "route violates native endpoint pairing",
            )?;
            for pair in route.segments.windows(2) {
                let left = self
                    .ports
                    .at_cut(self.graph, pair[0].source, pair[0].exit(), pair[0].reverse)?
                    .ok_or_else(|| invalid("undeclared outgoing switch"))?;
                let right = self
                    .ports
                    .at_cut(self.graph, pair[1].source, pair[1].entry(), pair[1].reverse)?
                    .ok_or_else(|| invalid("undeclared incoming switch"))?;
                require(
                    left.word == right.word,
                    "switch words/oriented cuts disagree",
                )?;
                for p in [left, right] {
                    let mut dna =
                        self.sources
                            .fetch(p.source, p.anchor, p.anchor + self.graph.k)?;
                    if p.reverse {
                        dna = crate::graph::reverse_complement(&dna);
                    }
                    require(dna == p.word, "switch DNA provenance mismatch")?;
                }
            }
        }
        require(
            resources_feasible(normalized.routes.iter().flat_map(|r| r.segments.clone())),
            "conflicting canonical source-span reuse across genome",
        )?;
        Ok(normalized)
    }
    pub fn route_counts(&mut self, route: &Route) -> io::Result<(Vec<Counts>, u64)> {
        if let Some(c) = self.cache.get(route) {
            self.cache_hits += 1;
            return Ok((c.clone(), 0));
        }
        require(
            !route.segments.is_empty()
                && route.segments.iter().all(|p| {
                    p.source < self.graph.lanes.len()
                        && p.start < p.end
                        && p.end <= self.graph.lanes[p.source].length
                }),
            "invalid route profile traversal",
        )?;
        require(
            resources_feasible(route.segments.clone()),
            "overlapping route profile resources",
        )?;
        let molecule = route.length()?;
        let mut result = Vec::new();
        let mut corrected = 0;
        for (li, &length) in self.graph.read_lengths.iter().enumerate() {
            let mut counts = Counts::new();
            let mut offset = 0;
            // Baseline owns a start at EVERY traversed base where a native read
            // exists, including starts later removed by the declared molecule end.
            for p in &route.segments {
                let (lo, hi) = native_domain(
                    p,
                    0,
                    p.end - p.start,
                    length,
                    self.graph.lanes[p.source].length,
                );
                self.store.integrate(
                    self.graph,
                    p.source,
                    li,
                    lo,
                    hi,
                    &mut counts,
                    false,
                    self.max_terms,
                )?;
            }
            let admitted = if length <= molecule {
                molecule - length + 1
            } else {
                0
            };
            let mut windows = Vec::new();
            for p in route.segments.iter().take(route.segments.len() - 1) {
                offset += p.end - p.start;
                windows.push((offset.saturating_sub(length - 1), offset));
            }
            windows.push((admitted, molecule));
            let windows = union(windows);
            for &(a, b) in &windows {
                let mut offset = 0;
                for p in &route.segments {
                    let end = offset + p.end - p.start;
                    let start = a.max(offset);
                    let stop = b.min(end);
                    if start < stop {
                        let (lo, hi) = native_domain(
                            p,
                            start - offset,
                            stop - offset,
                            length,
                            self.graph.lanes[p.source].length,
                        );
                        self.store.integrate(
                            self.graph,
                            p.source,
                            li,
                            lo,
                            hi,
                            &mut counts,
                            true,
                            self.max_terms,
                        )?;
                    }
                    offset = end;
                }
                let end = b.min(admitted);
                if a < end {
                    add(&mut corrected, end - a)?;
                    profiles::events(
                        self.panel,
                        molecule,
                        length,
                        a,
                        end,
                        self.graph.core_bp,
                        |lo, hi| self.sources.spell(route, lo, hi),
                        |e| {
                            for (t, q) in e.counts {
                                profiles::accumulate(
                                    &mut counts,
                                    t,
                                    q.checked_mul(e.end - e.start)
                                        .ok_or_else(|| invalid("corrected count overflow"))?,
                                    false,
                                    self.max_terms,
                                )?;
                            }
                            Ok(())
                        },
                    )?;
                }
            }
            result.push(counts);
        }
        // Charge route bindings as well: zero-anchor/zero-feature routes must not
        // create an unbounded cache of keys with zero charged feature terms.
        let terms: usize = result
            .iter()
            .map(|c| c.len())
            .sum::<usize>()
            .saturating_add(route.segments.len());
        if terms <= self.max_cache_terms {
            if self.cache_terms.saturating_add(terms) > self.max_cache_terms {
                self.cache.clear();
                self.cache_terms = 0;
            }
            self.cache.insert(route.clone(), result.clone());
            self.cache_terms += terms;
        }
        Ok((result, corrected))
    }
    pub fn evaluate(&mut self, assignment: &Assignment) -> io::Result<Evaluation> {
        let assignment = self.validate_assignment(assignment)?;
        let mut counts = vec![Counts::new(); self.graph.read_lengths.len()];
        let mut admitted_starts = vec![0; counts.len()];
        let mut corrected_starts = 0;
        let before = self.store.queried_starts;
        let cache_before = self.cache_hits;
        for route in &assignment.routes {
            let (q, corrected) = self.route_counts(route)?;
            add(&mut corrected_starts, corrected)?;
            let m = route.length()?;
            for (li, c) in q.into_iter().enumerate() {
                let l = self.graph.read_lengths[li];
                if l <= m {
                    add(&mut admitted_starts[li], m - l + 1)?;
                }
                for (t, n) in c {
                    profiles::accumulate(&mut counts[li], t, n, false, self.max_terms)?;
                }
            }
        }
        let mut keys: BTreeSet<_> = self.observed.keys().copied().collect();
        for c in &counts {
            keys.extend(c.keys().copied());
            require(
                keys.len() <= self.max_terms,
                "realized feature union exceeds resource cap; no partial objective",
            )?;
        }
        let mut factors = Vec::new();
        let mut relative_objective = 0.0;
        for tokens in keys {
            let q: Vec<_> = counts
                .iter()
                .map(|c| c.get(&tokens).copied().unwrap_or(0))
                .collect();
            let observed = self.observed.get(&tokens).copied().unwrap_or(0);
            let exposure: f64 = q
                .iter()
                .zip(&self.histogram)
                .map(|(&q, &n)| q as f64 * n as f64)
                .sum::<f64>()
                / self.denominator;
            let signal = self.depth * exposure;
            let loss = if signal == 0.0 {
                0.0
            } else {
                signal - observed as f64 * (signal / self.background).ln_1p()
            };
            require(
                signal.is_finite() && loss.is_finite(),
                "nonfinite relative objective",
            )?;
            relative_objective += loss;
            factors.push(Factor {
                tokens,
                counts_by_length: q,
                observed,
                signal,
                relative_loss: loss,
                positive_residual: observed > 0 && signal == 0.0,
            });
        }
        require(
            relative_objective.is_finite(),
            "nonfinite relative objective sum",
        )?;
        Ok(Evaluation{version:VERSION,model:MODEL.into(),assignment,objective_kind:"fixed-universe-background-relative-NLL-only".into(),relative_objective,universe:UNIVERSE.into(),factors,sample_positive_features:self.observed.len(),registry_definitions:self.graph.registry_count,symbolic_zero_terms:"all unmaterialized registry/route-rule features have zero observation and zero selected signal; exact relative contribution zero".into(),admitted_starts,corrected_starts,native_range_starts:self.store.queried_starts-before,route_cache_hits:self.cache_hits-cache_before,native_evaluation_complete:true,generation_rule_complete:true,topology_inventory_complete_for_panel:true,biological_topology_complete:false,sequence_emission_authorized:false})
    }
    /// Analytic conservative version of the independent nonnegative-signal
    /// relaxation: ln(r)<=r-1 gives -(c-beta)^2/beta. Outward operations also
    /// cover integer->f64 conversion. None denotes a valid -infinity overflow.
    pub fn lower_bound(&self) -> Option<f64> {
        let mut bound = 0.0f64;
        for &c in self.observed.values() {
            let upper = (c as f64).next_up();
            if upper > self.background {
                let delta = (upper - self.background).next_up();
                let term = ((delta * delta).next_up() / self.background).next_up();
                bound = (bound - term).next_down();
                if !bound.is_finite() {
                    return None;
                }
            }
        }
        Some(bound)
    }
}
/// Native forward start coordinates for a traversal-local half-open start range.
/// Reverse read starting at local t corresponds to x=end-L-t, not end-1-t.
pub(super) fn native_domain(
    p: &Segment,
    lo: u64,
    hi: u64,
    length: u64,
    source_length: u64,
) -> (u64, u64) {
    if length > source_length {
        return (0, 0);
    }
    let cap = i128::from(source_length - length + 1);
    let (a, b) = if p.reverse {
        (
            i128::from(p.end) - i128::from(length) - i128::from(hi) + 1,
            i128::from(p.end) - i128::from(length) - i128::from(lo) + 1,
        )
    } else {
        (i128::from(p.start + lo), i128::from(p.start + hi))
    };
    let a = a.clamp(0, cap) as u64;
    let b = b.clamp(0, cap) as u64;
    (a, b)
}
fn union(mut windows: Vec<(u64, u64)>) -> Vec<(u64, u64)> {
    windows.sort_unstable();
    let mut out: Vec<(u64, u64)> = Vec::new();
    for (a, b) in windows {
        if a >= b {
            continue;
        }
        if let Some(last) = out.last_mut() {
            if a <= last.1 {
                last.1 = last.1.max(b);
                continue;
            }
        }
        out.push((a, b));
    }
    out
}
