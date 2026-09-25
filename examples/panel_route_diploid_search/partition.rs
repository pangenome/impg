//! Experimental catalog-partition diplotypes over the complete MEM-subwalk space.
//!
//! A partition is a catalog group of homologous source intervals. Physical
//! occurrences remain distinct even when their sequences or profiles are equal.
//!
//! The production route registry currently projects evidence to adjacent
//! `[node,gap,node]` keys. That projection—not L150 or `WeightedBwt`—loses longer
//! linkage. This isolated scorer exercises the retained arbitrary odd-length
//! query space without changing the production scorer.
use super::{ensure, invalid};
use impg::{
    genome_inference::{catalog, mem_records, panel_routes as routes},
    graph::reverse_complement,
    sample_mem_bwt::{canonical, WeightedBwt},
    syng::SyngIndex,
};
use serde::Serialize;
use std::{collections::BTreeMap, io};

pub type FeatureKey = Vec<u64>;
pub type Profile = BTreeMap<FeatureKey, u64>;

#[derive(Clone, Debug, Serialize, PartialEq, Eq, PartialOrd, Ord)]
pub struct SourceRange {
    pub partition: usize,
    pub occurrence: usize,
    pub source: usize,
    pub start: u64,
    pub end: u64,
    pub reverse: bool,
}

#[derive(Clone, Debug)]
pub struct Allele {
    pub range: SourceRange,
    pub sequence: Vec<u8>,
    pub interior: Profile,
}

#[derive(Clone, Debug)]
pub struct Seam {
    pub left: SourceRange,
    pub right: SourceRange,
    pub profile: Profile,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub enum FeatureClass {
    Owned(usize),
    Boundary(usize),
    Deferred(String),
}

#[derive(Clone, Debug, Serialize)]
pub struct UniverseEntry {
    pub key: FeatureKey,
    pub observed: u64,
    pub class: FeatureClass,
}

#[derive(Clone, Debug)]
pub struct Universe {
    pub entries: Vec<UniverseEntry>,
}

#[derive(Clone, Debug, Serialize, PartialEq, Eq)]
pub struct FeatureSpectrumBin {
    pub node_count: usize,
    pub token_length: usize,
    pub features: u64,
    pub observed_multiplicity: u64,
    pub predicted_multiplicity: u64,
    pub local_features: u64,
    pub shared_features: u64,
}

/// Preserve every catalog occurrence. An explicit strand yields one oriented
/// allele; an unoriented BED3 occurrence yields both orientations rather than
/// silently inventing one.
pub fn catalog_ranges(catalog: &catalog::Catalog) -> io::Result<Vec<Vec<SourceRange>>> {
    let mut out = Vec::with_capacity(catalog.groups.len());
    for (partition, group) in catalog.groups.iter().enumerate() {
        let mut ranges = Vec::new();
        for &occurrence_id in &group.occurrences {
            let occurrence = catalog
                .occurrences
                .get(occurrence_id)
                .ok_or_else(|| invalid("catalog group occurrence is missing"))?;
            ensure(
                occurrence.id == occurrence_id
                    && occurrence.group == partition
                    && occurrence.source < catalog.sources.len()
                    && occurrence.interval.path == catalog.sources[occurrence.source].path
                    && occurrence.interval.start < occurrence.interval.end
                    && occurrence.interval.end <= catalog.sources[occurrence.source].length,
                "invalid catalog source interval binding",
            )?;
            let orientations: &[bool] = match occurrence.interval.strand.as_deref() {
                Some("+") => &[false],
                Some("-") => &[true],
                None => &[false, true],
                Some(_) => return Err(invalid("invalid catalog occurrence strand")),
            };
            for &reverse in orientations {
                ranges.push(SourceRange {
                    partition,
                    occurrence: occurrence_id,
                    source: occurrence.source,
                    start: occurrence.interval.start,
                    end: occurrence.interval.end,
                    reverse,
                });
            }
        }
        ensure(!ranges.is_empty(), "empty catalog partition")?;
        out.push(ranges);
    }
    Ok(out)
}

fn feature_cap(max_features: usize) -> io::Result<()> {
    ensure(
        max_features > 0 && max_features <= 50_000,
        "invalid maximal MEM subwalk feature cap",
    )
}

fn add(profile: &mut Profile, key: FeatureKey, n: u64, max_features: usize) -> io::Result<()> {
    feature_cap(max_features)?;
    if !profile.contains_key(&key) {
        ensure(
            profile.len() < max_features,
            "budget: maximal_mem_subwalk_features",
        )?;
    }
    let value = profile.entry(key).or_default();
    *value = value
        .checked_add(n)
        .ok_or_else(|| invalid("maximal MEM subwalk count overflow"))?;
    Ok(())
}

/// Add every node-to-node contiguous subwalk of every retained maximal record.
/// Each position contributes, so repeats and maximal-record multiplicities are
/// retained. Keys are RC-canonical and include single-node and full-record keys.
pub fn add_record_subwalks(
    profile: &mut Profile,
    record: &[u64],
    max_features: usize,
) -> io::Result<()> {
    add_record_subwalks_weighted(profile, record, 1, max_features)
}

/// Weighted form used by event-compressed replay. `multiplicity` is the exact
/// number of consecutive read starts represented by one invariant MEM query.
pub fn add_record_subwalks_weighted(
    profile: &mut Profile,
    record: &[u64],
    multiplicity: u64,
    max_features: usize,
) -> io::Result<()> {
    ensure(
        !record.is_empty() && record.len() % 2 == 1 && multiplicity > 0,
        "invalid weighted maximal MEM token record",
    )?;
    for start in (0..record.len()).step_by(2) {
        for end in (start..record.len()).step_by(2) {
            add(
                profile,
                canonical(&record[start..=end]),
                multiplicity,
                max_features,
            )?;
        }
    }
    Ok(())
}

pub fn profile_read(panel: &SyngIndex, read: &[u8], max_features: usize) -> io::Result<Profile> {
    feature_cap(max_features)?;
    let mut profile = Profile::new();
    for record in mem_records::canonical_mem_records(panel, read)? {
        add_record_subwalks(&mut profile, &record, max_features)?;
    }
    Ok(profile)
}

fn merge(target: &mut Profile, source: &Profile, max_features: usize) -> io::Result<()> {
    for (key, &count) in source {
        add(target, key.clone(), count, max_features)?;
    }
    Ok(())
}

/// Profile complete read starts wholly contained by one exact allele subrange.
pub fn profile_interior(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    max_features: usize,
) -> io::Result<Profile> {
    ensure(read_length > 0, "zero partition read length")?;
    let mut profile = Profile::new();
    if sequence.len() < read_length {
        return Ok(profile);
    }
    for start in 0..=sequence.len() - read_length {
        merge(
            &mut profile,
            &profile_read(panel, &sequence[start..start + read_length], max_features)?,
            max_features,
        )?;
    }
    Ok(profile)
}

/// Profile only complete windows that cross a compatible physical join. The
/// partition cut is never passed to the matcher as a molecule terminal.
pub fn profile_seam(
    panel: &SyngIndex,
    left: &[u8],
    right: &[u8],
    read_length: usize,
    max_features: usize,
) -> io::Result<Profile> {
    ensure(read_length > 0, "zero partition read length")?;
    if read_length == 1 || left.is_empty() || right.is_empty() {
        return Ok(Profile::new());
    }
    let left_take = left.len().min(read_length - 1);
    let right_take = right.len().min(read_length - 1);
    let mut context = left[left.len() - left_take..].to_vec();
    let boundary = context.len();
    context.extend_from_slice(&right[..right_take]);
    let mut profile = Profile::new();
    if context.len() >= read_length {
        for start in 0..=context.len() - read_length {
            if start < boundary && start + read_length > boundary {
                merge(
                    &mut profile,
                    &profile_read(panel, &context[start..start + read_length], max_features)?,
                    max_features,
                )?;
            }
        }
    }
    Ok(profile)
}

pub fn allele(
    range: SourceRange,
    sequence: Vec<u8>,
    panel: &SyngIndex,
    read_length: usize,
    max_features: usize,
) -> io::Result<Allele> {
    ensure(
        sequence.len() as u64 == range.end - range.start && sequence.len() >= read_length,
        "partition allele is not a complete read-spanning source interval",
    )?;
    let sequence = if range.reverse {
        reverse_complement(&sequence)
    } else {
        sequence
    };
    let interior = profile_interior(panel, &sequence, read_length, max_features)?;
    Ok(Allele {
        range,
        sequence,
        interior,
    })
}

pub fn catalog_alleles(
    catalog: &catalog::Catalog,
    panel: &SyngIndex,
    read_length: usize,
    max_features: usize,
    mut fetch: impl FnMut(usize, u64, u64) -> io::Result<Vec<u8>>,
) -> io::Result<Vec<Vec<Allele>>> {
    catalog_ranges(catalog)?
        .into_iter()
        .map(|partition| {
            partition
                .into_iter()
                .map(|range| {
                    let sequence = fetch(range.source, range.start, range.end)?;
                    allele(range, sequence, panel, read_length, max_features)
                })
                .collect()
        })
        .collect()
}

pub fn seams_from_public_links(
    partitions: &[Vec<Allele>],
    links: &[catalog::Link],
    panel: &SyngIndex,
    read_length: usize,
    max_features: usize,
) -> io::Result<Vec<Vec<Seam>>> {
    ensure(!partitions.is_empty(), "empty linked partition domain")?;
    let mut out = vec![Vec::new(); partitions.len() - 1];
    for link in links.iter().filter(|link| link.oriented_continuation) {
        let left_candidates = partitions
            .iter()
            .flatten()
            .filter(|allele| allele.range.occurrence == link.from)
            .collect::<Vec<_>>();
        let right_candidates = partitions
            .iter()
            .flatten()
            .filter(|allele| allele.range.occurrence == link.to)
            .collect::<Vec<_>>();
        for left in &left_candidates {
            for right in &right_candidates {
                if left.range.reverse != right.range.reverse {
                    continue;
                }
                ensure(
                    left.range.partition + 1 == right.range.partition,
                    "public link skips partition adjacency",
                )?;
                out[left.range.partition].push(seam(
                    left,
                    right,
                    panel,
                    read_length,
                    max_features,
                )?);
            }
        }
    }
    for boundary in &out {
        ensure(
            !boundary.is_empty(),
            "partition boundary has no public legal seam",
        )?;
    }
    Ok(out)
}

pub fn seam(
    left: &Allele,
    right: &Allele,
    panel: &SyngIndex,
    read_length: usize,
    max_features: usize,
) -> io::Result<Seam> {
    ensure(
        left.range.partition + 1 == right.range.partition,
        "nonadjacent partition seam",
    )?;
    Ok(Seam {
        left: left.range.clone(),
        right: right.range.clone(),
        profile: profile_seam(
            panel,
            &left.sequence,
            &right.sequence,
            read_length,
            max_features,
        )?,
    })
}

/// Freeze the disjoint scoring ledger from public allele interiors and every
/// public legal seam. Sample reads are never enumerated here; observations are
/// point queries against the aggregate `WeightedBwt` only.
pub fn universe(
    alleles: &[Vec<Allele>],
    seams: &[Vec<Seam>],
    sample: &WeightedBwt,
    max_features: usize,
) -> io::Result<Universe> {
    feature_cap(max_features)?;
    ensure(
        !alleles.is_empty() && seams.len() + 1 == alleles.len(),
        "invalid feature classification domain",
    )?;
    let mut incidence: BTreeMap<FeatureKey, (Vec<usize>, Vec<usize>)> = BTreeMap::new();
    for (partition, candidates) in alleles.iter().enumerate() {
        for allele in candidates {
            for key in allele.interior.keys() {
                let partitions = &mut incidence.entry(key.clone()).or_default().0;
                if !partitions.contains(&partition) {
                    partitions.push(partition);
                }
            }
        }
    }
    for (boundary, candidates) in seams.iter().enumerate() {
        for seam in candidates {
            for key in seam.profile.keys() {
                let boundaries = &mut incidence.entry(key.clone()).or_default().1;
                if !boundaries.contains(&boundary) {
                    boundaries.push(boundary);
                }
            }
        }
    }
    ensure(
        incidence.len() <= max_features,
        "budget: fixed_feature_universe",
    )?;
    let entries = incidence
        .into_iter()
        .map(|(key, (partitions, boundaries))| {
            let class = match (partitions.as_slice(), boundaries.as_slice()) {
                ([partition], []) => FeatureClass::Owned(*partition),
                ([], [boundary]) => FeatureClass::Boundary(*boundary),
                (p, b) => FeatureClass::Deferred(if !p.is_empty() && !b.is_empty() {
                    "interior_and_boundary_contributions".into()
                } else if p.len() > 1 {
                    "multiple_partition_interiors".into()
                } else if b.len() > 1 {
                    "multiple_adjacencies".into()
                } else {
                    "nonlocal_candidate_feature".into()
                }),
            };
            Ok(UniverseEntry {
                observed: sample.count(&key)?,
                key,
                class,
            })
        })
        .collect::<io::Result<Vec<_>>>()?;
    let universe = Universe { entries };
    let classified = universe.entries.len();
    let unique = universe
        .entries
        .iter()
        .map(|entry| &entry.key)
        .collect::<std::collections::BTreeSet<_>>()
        .len();
    ensure(classified == unique, "non-disjoint feature ledger")?;
    Ok(universe)
}

/// Summarize the fixed universe without filtering any feature length. A feature
/// is local only when its candidate incidence is confined to exactly one supplied
/// partition; seam-only, multi-partition, and sample-only keys are reported as
/// shared/nonlocal. Multiplicities are descriptive evidence, not definitive calls.
pub fn feature_spectrum(
    universe: &Universe,
    predicted: &Profile,
) -> io::Result<Vec<FeatureSpectrumBin>> {
    let mut bins: BTreeMap<(usize, usize), FeatureSpectrumBin> = BTreeMap::new();
    for entry in &universe.entries {
        ensure(
            !entry.key.is_empty() && entry.key.len() % 2 == 1,
            "invalid feature in spectrum audit",
        )?;
        let node_count = (entry.key.len() + 1) / 2;
        let local = matches!(entry.class, FeatureClass::Owned(_));
        let bin = bins
            .entry((node_count, entry.key.len()))
            .or_insert(FeatureSpectrumBin {
                node_count,
                token_length: entry.key.len(),
                features: 0,
                observed_multiplicity: 0,
                predicted_multiplicity: 0,
                local_features: 0,
                shared_features: 0,
            });
        bin.features = bin
            .features
            .checked_add(1)
            .ok_or_else(|| invalid("feature spectrum overflow"))?;
        bin.observed_multiplicity = bin
            .observed_multiplicity
            .checked_add(entry.observed)
            .ok_or_else(|| invalid("feature spectrum overflow"))?;
        bin.predicted_multiplicity = bin
            .predicted_multiplicity
            .checked_add(predicted.get(&entry.key).copied().unwrap_or(0))
            .ok_or_else(|| invalid("feature spectrum overflow"))?;
        if local {
            bin.local_features = bin
                .local_features
                .checked_add(1)
                .ok_or_else(|| invalid("feature spectrum overflow"))?;
        } else {
            bin.shared_features = bin
                .shared_features
                .checked_add(1)
                .ok_or_else(|| invalid("feature spectrum overflow"))?;
        }
    }
    Ok(bins.into_values().collect())
}

pub fn sum_profiles<'a>(
    profiles: impl IntoIterator<Item = &'a Profile>,
    max_features: usize,
) -> io::Result<Profile> {
    let mut total = Profile::new();
    for profile in profiles {
        merge(&mut total, profile, max_features)?;
    }
    Ok(total)
}

#[derive(Clone, Debug, Serialize)]
pub struct ScoreModel {
    pub read_length: u64,
    pub histogram: u64,
    pub denominator: f64,
    pub depth: f64,
    pub background: f64,
}

impl ScoreModel {
    pub fn validate(&self) -> io::Result<()> {
        ensure(
            self.read_length == 150
                && self.histogram > 0
                && self.denominator == self.read_length as f64 * self.histogram as f64
                && self.depth.is_finite()
                && self.depth > 0.0
                && self.background.is_finite()
                && self.background > 0.0,
            "invalid partition score normalization",
        )
    }
    pub fn loss(&self, q: u64, observed: u64) -> io::Result<f64> {
        self.validate()?;
        ensure(observed <= 1 << 53, "partition observation precision")?;
        // Mirror the exact scorer's precision bound on the q*histogram product
        // before any floating-point exposure conversion.
        let product = q
            .checked_mul(self.histogram)
            .ok_or_else(|| invalid("partition exposure product overflow"))?;
        ensure(
            product <= 1 << 53,
            "partition exposure conversion precision limit",
        )?;
        let signal = product as f64 * self.depth / self.denominator;
        let loss = signal - observed as f64 * (signal / self.background).ln_1p();
        ensure(loss.is_finite(), "nonfinite partition factor loss")?;
        Ok(loss)
    }
}

fn score_where(
    a: &Profile,
    b: &Profile,
    universe: &Universe,
    model: &ScoreModel,
    mut include: impl FnMut(&FeatureClass) -> bool,
) -> io::Result<f64> {
    let mut objective = 0.0;
    for entry in universe
        .entries
        .iter()
        .filter(|entry| include(&entry.class))
    {
        let q = a
            .get(&entry.key)
            .copied()
            .unwrap_or(0)
            .checked_add(b.get(&entry.key).copied().unwrap_or(0))
            .ok_or_else(|| invalid("diploid MEM subwalk count overflow"))?;
        objective += model.loss(q, entry.observed)?;
    }
    ensure(objective.is_finite(), "nonfinite partition diplotype loss")?;
    Ok(objective)
}

/// Authoritative complete score includes owned, boundary, and deferred keys.
pub fn score_pair(
    a: &Profile,
    b: &Profile,
    universe: &Universe,
    model: &ScoreModel,
) -> io::Result<f64> {
    score_where(a, b, universe, model, |_| true)
}

#[derive(Clone, Debug, Serialize)]
pub struct Diplotype {
    pub alleles: [usize; 2],
    pub loss: f64,
}

pub fn local_diplotypes(
    partition: usize,
    alleles: &[Allele],
    universe: &Universe,
    model: &ScoreModel,
) -> io::Result<Vec<Diplotype>> {
    let mut out = Vec::new();
    for a in 0..alleles.len() {
        for b in a..alleles.len() {
            out.push(Diplotype {
                alleles: [a, b],
                loss: score_where(
                    &alleles[a].interior,
                    &alleles[b].interior,
                    universe,
                    model,
                    |class| matches!(class, FeatureClass::Owned(owner) if *owner == partition),
                )?,
            });
        }
    }
    out.sort_by(|a, b| {
        a.loss
            .total_cmp(&b.loss)
            .then_with(|| a.alleles.cmp(&b.alleles))
    });
    Ok(out)
}

/// Return both homolog matchings when physically compatible. Equal profiles are
/// not collapsed: source coordinates remain part of every matching.
pub fn boundary_matchings(
    left: &Diplotype,
    right: &Diplotype,
    left_alleles: &[Allele],
    right_alleles: &[Allele],
    seams: &[Seam],
    max_features: usize,
) -> io::Result<Vec<([[usize; 2]; 2], Profile)>> {
    ensure(
        left.alleles.iter().all(|&i| i < left_alleles.len())
            && right.alleles.iter().all(|&i| i < right_alleles.len()),
        "diplotype allele index out of range",
    )?;
    let mut out = Vec::new();
    for matching in [
        [
            [left.alleles[0], right.alleles[0]],
            [left.alleles[1], right.alleles[1]],
        ],
        [
            [left.alleles[0], right.alleles[1]],
            [left.alleles[1], right.alleles[0]],
        ],
    ] {
        let mut selected = Vec::new();
        for [left_index, right_index] in matching {
            let left_range = &left_alleles[left_index].range;
            let right_range = &right_alleles[right_index].range;
            let Some(profile) = seams
                .iter()
                .find(|s| &s.left == left_range && &s.right == right_range)
                .map(|s| &s.profile)
            else {
                selected.clear();
                break;
            };
            selected.push(profile);
        }
        if selected.len() == 2 {
            out.push((matching, sum_profiles(selected, max_features)?));
        }
    }
    Ok(out)
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct Accounting {
    pub work: u64,
    pub profile_work: u64,
    pub paired_scores: u64,
    pub haploid_profiles: u64,
    pub state_bytes: u64,
    pub retained_states: u64,
    pub matchings: u64,
    pub reconstructions: u64,
}

impl Accounting {
    fn charge(value: &mut u64, amount: u64, cap: u64, reason: &str) -> io::Result<()> {
        let next = value.checked_add(amount).ok_or_else(|| invalid(reason))?;
        ensure(next <= cap, reason)?;
        *value = next;
        Ok(())
    }
    fn work(&mut self, n: u64) -> io::Result<()> {
        self.work_bounded(n, 4_000_000)
    }
    fn work_bounded(&mut self, n: u64, max_work: u64) -> io::Result<()> {
        ensure(max_work <= 4_000_000, "invalid partition work cap")?;
        Self::charge(&mut self.work, n, max_work, "budget: partition_work")
    }
    fn profile(&mut self, n: u64) -> io::Result<()> {
        Self::charge(
            &mut self.profile_work,
            n,
            8_000_000,
            "budget: partition_profile_work",
        )
    }
    fn state(&mut self, n: u64) -> io::Result<()> {
        Self::charge(
            &mut self.state_bytes,
            n,
            134_217_728,
            "budget: partition_state",
        )
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct Backpointer {
    pub predecessor: Option<usize>,
    pub diplotype: [usize; 2],
    pub matching: Option<u8>,
}

#[derive(Clone, Debug)]
pub struct ChainState {
    pub choices: [Vec<usize>; 2],
    pub spans: [Vec<SourceRange>; 2],
    pub proposal_loss: f64,
    pub backpointers: Vec<Backpointer>,
}

#[derive(Clone, Debug)]
pub struct ChainSearch {
    pub complete: bool,
    pub stop_reason: Option<String>,
    pub states: Vec<ChainState>,
    pub accounting: Accounting,
}

fn logical_state_bytes(partitions: usize) -> u64 {
    // Covers Vec headers plus allocated choices, two SourceRange ledgers,
    // backpointers, and scalar fields with margin; intentionally conservative.
    1024 + 256 * partitions as u64
}

fn class_terms(universe: &Universe, class: &FeatureClass) -> u64 {
    universe
        .entries
        .iter()
        .filter(|entry| &entry.class == class)
        .count() as u64
}

/// Partition-major local diplotype chaining. It never enumerates complete
/// haplotypes independently: each retained state is created by one of the two
/// homolog matchings and carries independent source-span ledgers and backpointers.
pub fn chain_diplotypes(
    partitions: &[Vec<Allele>],
    seams: &[Vec<Seam>],
    universe: &Universe,
    model: &ScoreModel,
) -> io::Result<ChainSearch> {
    chain_diplotypes_bounded(partitions, seams, universe, model, 50_000, 4_000_000)
}

pub fn chain_diplotypes_bounded(
    partitions: &[Vec<Allele>],
    seams: &[Vec<Seam>],
    universe: &Universe,
    model: &ScoreModel,
    max_features: usize,
    max_work: u64,
) -> io::Result<ChainSearch> {
    feature_cap(max_features)?;
    ensure(
        max_work > 0 && max_work <= 4_000_000,
        "invalid partition work cap",
    )?;
    ensure(
        !partitions.is_empty() && partitions.len() <= 64 && seams.len() + 1 == partitions.len(),
        "invalid partition chain",
    )?;
    let mut accounting = Accounting::default();
    for alleles in partitions {
        for allele in alleles {
            if let Err(error) =
                accounting.profile((allele.sequence.len() - model.read_length as usize + 1) as u64)
            {
                return Ok(ChainSearch {
                    complete: false,
                    stop_reason: Some(error.to_string()),
                    states: Vec::new(),
                    accounting,
                });
            }
        }
    }
    for boundary in seams {
        if let Err(error) = accounting.profile(boundary.len() as u64 * (model.read_length - 1)) {
            return Ok(ChainSearch {
                complete: false,
                stop_reason: Some(error.to_string()),
                states: Vec::new(),
                accounting,
            });
        }
    }
    let mut locals = Vec::new();
    for (partition, alleles) in partitions.iter().enumerate() {
        let pairs = (alleles.len() * (alleles.len() + 1) / 2) as u64;
        let Some(local_work) =
            pairs.checked_mul(class_terms(universe, &FeatureClass::Owned(partition)))
        else {
            return Ok(ChainSearch {
                complete: false,
                stop_reason: Some("partition work overflow".into()),
                states: Vec::new(),
                accounting,
            });
        };
        if let Err(error) = accounting.work_bounded(local_work, max_work) {
            return Ok(ChainSearch {
                complete: false,
                stop_reason: Some(error.to_string()),
                states: Vec::new(),
                accounting,
            });
        }
        match local_diplotypes(partition, alleles, universe, model) {
            Ok(diplotypes) => locals.push(diplotypes),
            Err(error) => {
                return Ok(ChainSearch {
                    complete: false,
                    stop_reason: Some(error.to_string()),
                    states: Vec::new(),
                    accounting,
                });
            }
        }
    }
    let mut states = Vec::new();
    for dip in &locals[0] {
        let spans = [
            vec![partitions[0][dip.alleles[0]].range.clone()],
            vec![partitions[0][dip.alleles[1]].range.clone()],
        ];
        if let Err(error) = accounting.state(logical_state_bytes(1)) {
            return Ok(ChainSearch {
                complete: false,
                stop_reason: Some(error.to_string()),
                states,
                accounting,
            });
        }
        accounting.retained_states += 1;
        states.push(ChainState {
            choices: [vec![dip.alleles[0]], vec![dip.alleles[1]]],
            spans,
            proposal_loss: dip.loss,
            backpointers: vec![Backpointer {
                predecessor: None,
                diplotype: dip.alleles,
                matching: None,
            }],
        });
    }
    for partition in 1..partitions.len() {
        let previous = states;
        let mut next = Vec::new();
        for (predecessor, state) in previous.iter().enumerate() {
            for dip in &locals[partition] {
                let matching_count = if dip.alleles[0] == dip.alleles[1] {
                    1
                } else {
                    2
                };
                for matching_id in 0..matching_count {
                    let next_indices = if matching_id == 0 {
                        dip.alleles
                    } else {
                        [dip.alleles[1], dip.alleles[0]]
                    };
                    accounting.matchings += 1;
                    if let Err(error) = accounting.work_bounded(1, max_work) {
                        return Ok(ChainSearch {
                            complete: false,
                            stop_reason: Some(error.to_string()),
                            states: previous.clone(),
                            accounting,
                        });
                    }
                    let mut seam_profiles = Vec::new();
                    let mut legal = true;
                    for copy in 0..2 {
                        let left =
                            &partitions[partition - 1][*state.choices[copy].last().unwrap()].range;
                        let right = &partitions[partition][next_indices[copy]].range;
                        let Some(found) = seams[partition - 1]
                            .iter()
                            .find(|seam| &seam.left == left && &seam.right == right)
                        else {
                            legal = false;
                            break;
                        };
                        let mut spans = state.spans[copy].clone();
                        spans.push(right.clone());
                        if !source_spans_feasible(&spans) {
                            legal = false;
                            break;
                        }
                        seam_profiles.push(&found.profile);
                    }
                    if !legal {
                        continue;
                    }
                    let boundary_profile = match sum_profiles(seam_profiles, max_features) {
                        Ok(profile) => profile,
                        Err(error) => {
                            return Ok(ChainSearch {
                                complete: false,
                                stop_reason: Some(error.to_string()),
                                states: previous.clone(),
                                accounting,
                            });
                        }
                    };
                    if let Err(error) = accounting.work_bounded(
                        class_terms(universe, &FeatureClass::Boundary(partition - 1)),
                        max_work,
                    ) {
                        return Ok(ChainSearch {
                            complete: false,
                            stop_reason: Some(error.to_string()),
                            states: previous.clone(),
                            accounting,
                        });
                    }
                    let transition_loss = match score_where(
                        &boundary_profile,
                        &Profile::new(),
                        universe,
                        model,
                        |class| matches!(class, FeatureClass::Boundary(owner) if *owner == partition - 1),
                    ) {
                        Ok(loss) => loss,
                        Err(error) => {
                            return Ok(ChainSearch {
                                complete: false,
                                stop_reason: Some(error.to_string()),
                                states: previous.clone(),
                                accounting,
                            });
                        }
                    };
                    if let Err(error) = accounting.state(logical_state_bytes(partition + 1)) {
                        return Ok(ChainSearch {
                            complete: false,
                            stop_reason: Some(error.to_string()),
                            states: previous.clone(),
                            accounting,
                        });
                    }
                    let mut child = state.clone();
                    for copy in 0..2 {
                        child.choices[copy].push(next_indices[copy]);
                        child.spans[copy]
                            .push(partitions[partition][next_indices[copy]].range.clone());
                    }
                    child.proposal_loss += dip.loss + transition_loss;
                    child.backpointers.push(Backpointer {
                        predecessor: Some(predecessor),
                        diplotype: dip.alleles,
                        matching: Some(matching_id as u8),
                    });
                    accounting.retained_states += 1;
                    next.push(child);
                }
            }
        }
        states = next;
        if states.is_empty() {
            return Ok(ChainSearch {
                complete: true,
                stop_reason: None,
                states,
                accounting,
            });
        }
    }
    Ok(ChainSearch {
        complete: true,
        stop_reason: None,
        states,
        accounting,
    })
}

#[derive(Clone, Debug, Serialize)]
pub struct ChainedCandidate {
    pub state: usize,
    pub objective: f64,
}

pub fn rescore_chained(
    search: &mut ChainSearch,
    panel: &SyngIndex,
    partitions: &[Vec<Allele>],
    universe: &Universe,
    model: &ScoreModel,
    fixed_suffix: &[u8],
    max_features: usize,
) -> io::Result<Vec<ChainedCandidate>> {
    let mut out = Vec::new();
    'states: for (state_index, state) in search.states.iter().enumerate() {
        if search.accounting.paired_scores >= 2048 || search.accounting.haploid_profiles + 2 > 4096
        {
            search.complete = false;
            search.stop_reason = Some("budget: complete_partition_rescore".into());
            break;
        }
        search.accounting.paired_scores += 1;
        search.accounting.haploid_profiles += 2;
        search.accounting.reconstructions += 2;
        if let Err(error) = search.accounting.work(universe.entries.len() as u64) {
            search.complete = false;
            search.stop_reason = Some(error.to_string());
            break;
        }
        let mut profiles = Vec::with_capacity(2);
        for copy in 0..2 {
            let choices = &state.choices[copy];
            if choices.len() != partitions.len() {
                // A cap-stopped chain retains partial-depth states; rescoring
                // them stays an incomplete-search result, never a hard error.
                search.complete = false;
                search.stop_reason = Some("incomplete chained reconstruction".into());
                break 'states;
            }
            let sequence_len = choices.iter().enumerate().try_fold(
                fixed_suffix.len(),
                |length, (partition, &choice)| {
                    ensure(
                        choice < partitions[partition].len(),
                        "invalid chained allele",
                    )?;
                    length
                        .checked_add(partitions[partition][choice].sequence.len())
                        .ok_or_else(|| invalid("chained sequence length overflow"))
                },
            )?;
            ensure(
                sequence_len >= model.read_length as usize,
                "short chained sequence",
            )?;
            if let Err(error) = search
                .accounting
                .profile((sequence_len - model.read_length as usize + 1) as u64)
            {
                search.complete = false;
                search.stop_reason = Some(error.to_string());
                break 'states;
            }
            let mut sequence = Vec::with_capacity(sequence_len);
            for (partition, &choice) in choices.iter().enumerate() {
                sequence.extend_from_slice(&partitions[partition][choice].sequence);
            }
            sequence.extend_from_slice(fixed_suffix);
            if sequence.len() != sequence_len {
                return Err(invalid("chained sequence allocation mismatch"));
            }
            profiles.push(profile_complete_sequence(
                panel,
                &sequence,
                model.read_length as usize,
                max_features,
            )?);
        }
        out.push(ChainedCandidate {
            state: state_index,
            objective: score_pair(&profiles[0], &profiles[1], universe, model)?,
        });
    }
    out.sort_by(|a, b| {
        a.objective
            .total_cmp(&b.objective)
            .then_with(|| a.state.cmp(&b.state))
    });
    Ok(out)
}

pub fn reconstruct_sequence(chain: &[&Allele]) -> Vec<u8> {
    chain
        .iter()
        .flat_map(|a| a.sequence.iter().copied())
        .collect()
}

pub fn reconstruct_route(chain: &[&Allele]) -> routes::Route {
    routes::Route {
        segments: chain
            .iter()
            .map(|a| routes::Segment {
                source: a.range.source,
                start: a.range.start,
                end: a.range.end,
                reverse: a.range.reverse,
            })
            .collect(),
    }
}

/// Replace disjoint intervals in one native route while retaining fixed terminal
/// and other-molecule routes. Public validation remains authoritative.
pub fn reconstruct_assignment(
    template: &routes::Assignment,
    route_slot: usize,
    donor_ranges: &[SourceRange],
    target_source: usize,
    target_length: u64,
) -> io::Result<routes::Assignment> {
    ensure(
        route_slot < template.routes.len(),
        "invalid assignment route slot",
    )?;
    let mut donors = donor_ranges.to_vec();
    donors.sort_by_key(|range| (range.start, range.end));
    ensure(
        source_spans_feasible(&donors),
        "overlapping assignment donor ranges",
    )?;
    let mut segments = Vec::new();
    let mut cursor = 0;
    for donor in donors {
        ensure(
            !donor.reverse && cursor <= donor.start && donor.end <= target_length,
            "invalid assignment donor replacement",
        )?;
        if cursor < donor.start {
            segments.push(routes::Segment {
                source: target_source,
                start: cursor,
                end: donor.start,
                reverse: false,
            });
        }
        segments.push(routes::Segment {
            source: donor.source,
            start: donor.start,
            end: donor.end,
            reverse: false,
        });
        cursor = donor.end;
    }
    if cursor < target_length {
        segments.push(routes::Segment {
            source: target_source,
            start: cursor,
            end: target_length,
            reverse: false,
        });
    }
    ensure(segments.len() <= 64, "budget: partition_segments")?;
    let mut assignment = template.clone();
    assignment.routes[route_slot] = routes::Route { segments };
    Ok(assignment)
}

pub fn profile_complete_sequence(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    max_features: usize,
) -> io::Result<Profile> {
    profile_interior(panel, sequence, read_length, max_features)
}

#[derive(Clone, Debug)]
pub struct Haplotype {
    pub choices: Vec<usize>,
    pub ranges: Vec<SourceRange>,
    pub sequence: Vec<u8>,
    pub route: routes::Route,
    pub profile: Profile,
}

fn source_spans_feasible(ranges: &[SourceRange]) -> bool {
    for (i, a) in ranges.iter().enumerate() {
        if ranges[i + 1..]
            .iter()
            .any(|b| a.source == b.source && a.start < b.end && b.start < a.end)
        {
            return false;
        }
    }
    true
}

/// Exhaustively enumerate the declared public finite chain domain. Compatibility
/// comes only from the supplied seam table, never sample or truth support.
pub fn enumerate_haplotypes_assessment_oracle(
    panel: &SyngIndex,
    partitions: &[Vec<Allele>],
    seams: &[Vec<Seam>],
    read_length: usize,
    max_features: usize,
    max_haplotypes: usize,
) -> io::Result<Vec<Haplotype>> {
    ensure(
        !partitions.is_empty()
            && partitions.len() <= 64
            && seams.len() + 1 == partitions.len()
            && max_haplotypes > 0
            && max_haplotypes <= 4096,
        "invalid partition chain domain",
    )?;
    let mut choices = vec![Vec::<usize>::new()];
    for (partition, alleles) in partitions.iter().enumerate() {
        ensure(!alleles.is_empty(), "empty partition allele domain")?;
        let mut next = Vec::new();
        for prefix in choices {
            for index in 0..alleles.len() {
                let compatible = partition == 0
                    || seams[partition - 1].iter().any(|seam| {
                        seam.left == partitions[partition - 1][*prefix.last().unwrap()].range
                            && seam.right == alleles[index].range
                    });
                if compatible {
                    ensure(next.len() < max_haplotypes, "budget: partition_haplotypes")?;
                    let mut child = prefix.clone();
                    child.push(index);
                    next.push(child);
                }
            }
        }
        choices = next;
    }
    let mut out = Vec::with_capacity(choices.len());
    for choices in choices {
        let chain = choices
            .iter()
            .enumerate()
            .map(|(partition, &allele)| &partitions[partition][allele])
            .collect::<Vec<_>>();
        let ranges = chain.iter().map(|a| a.range.clone()).collect::<Vec<_>>();
        if !source_spans_feasible(&ranges) {
            continue;
        }
        let sequence = reconstruct_sequence(&chain);
        let route = reconstruct_route(&chain);
        let profile = profile_complete_sequence(panel, &sequence, read_length, max_features)?;
        out.push(Haplotype {
            choices,
            ranges,
            sequence,
            route,
            profile,
        });
    }
    ensure(!out.is_empty(), "no complete physical partition haplotype")?;
    Ok(out)
}

#[derive(Clone, Debug, Serialize)]
pub struct RankedPair {
    pub haplotypes: [usize; 2],
    pub objective: f64,
}

pub fn rank_pairs_assessment_oracle(
    haplotypes: &[Haplotype],
    universe: &Universe,
    model: &ScoreModel,
    max_pairs: usize,
) -> io::Result<Vec<RankedPair>> {
    let required = haplotypes
        .len()
        .checked_mul(haplotypes.len() + 1)
        .and_then(|n| n.checked_div(2))
        .ok_or_else(|| invalid("partition pair domain overflow"))?;
    ensure(
        max_pairs > 0 && max_pairs <= 2048 && required <= max_pairs,
        "budget: partition_pair_scores",
    )?;
    let mut out = Vec::with_capacity(required);
    for a in 0..haplotypes.len() {
        for b in a..haplotypes.len() {
            out.push(RankedPair {
                haplotypes: [a, b],
                objective: score_pair(
                    &haplotypes[a].profile,
                    &haplotypes[b].profile,
                    universe,
                    model,
                )?,
            });
        }
    }
    out.sort_by(|a, b| {
        a.objective
            .total_cmp(&b.objective)
            .then_with(|| a.haplotypes.cmp(&b.haplotypes))
    });
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use impg::{
        sample_mem_bwt::{encode_walk, WeightedBwt},
        syng::SyncmerParams,
    };

    fn dna(n: usize, mut state: u64) -> Vec<u8> {
        (0..n)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state & 3) as usize]
            })
            .collect()
    }

    fn range(
        partition: usize,
        occurrence: usize,
        source: usize,
        start: usize,
        end: usize,
    ) -> SourceRange {
        SourceRange {
            partition,
            occurrence,
            source,
            start: start as u64,
            end: end as u64,
            reverse: false,
        }
    }

    fn sample_bwt(panel: &SyngIndex, reads: &[Vec<u8>]) -> WeightedBwt {
        let mut records = BTreeMap::new();
        for read in reads {
            for record in mem_records::canonical_mem_records(panel, read).unwrap() {
                *records.entry(record).or_insert(0u64) += 1;
            }
        }
        WeightedBwt::build(&records).unwrap()
    }

    #[test]
    fn catalog_groups_preserve_occurrence_coordinates_multiplicity_and_orientation() {
        let catalog: catalog::Catalog = serde_json::from_value(serde_json::json!({
            "version": 1,
            "panel": {"checksum_algorithm":"test","sidecars":[]},
            "catalog_accepted": false,
            "validation": "test",
            "sources": [
                {"id":0,"path":"A#0#chr","length":500},
                {"id":1,"path":"B#0#chr","length":500}
            ],
            "groups": [{"id":"p0","scaffold":null,"occurrences":[0,1]}],
            "occurrences": [
                {"id":0,"group":0,"source":0,"interval":{"path":"A#0#chr","start":10,"end":210,"strand":"-"},"fully_contained_anchors":0,"anchor_status":"test","feature_multiplicities":{}},
                {"id":1,"group":0,"source":1,"interval":{"path":"B#0#chr","start":20,"end":220,"strand":null},"fully_contained_anchors":0,"anchor_status":"test","feature_multiplicities":{}}
            ],
            "features": [],
            "links": []
        }))
        .unwrap();
        catalog.validate().unwrap();
        let ranges = catalog_ranges(&catalog).unwrap();
        assert_eq!(ranges.len(), 1);
        assert_eq!(ranges[0].len(), 3);
        assert_eq!(
            (
                ranges[0][0].occurrence,
                ranges[0][0].start,
                ranges[0][0].end,
                ranges[0][0].reverse
            ),
            (0, 10, 210, true)
        );
        assert_eq!(ranges[0][1].occurrence, 1);
        assert!(!ranges[0][1].reverse);
        assert_eq!(ranges[0][2].occurrence, 1);
        assert!(ranges[0][2].reverse);

        let sequences = vec![dna(500, 91), dna(500, 92)];
        let panel = SyngIndex::build(
            SyncmerParams::default(),
            vec![
                ("A#0#chr".into(), sequences[0].clone()),
                ("B#0#chr".into(), sequences[1].clone()),
            ]
            .into_iter(),
        );
        let typed = catalog_alleles(&catalog, &panel, 150, 50_000, |source, start, end| {
            Ok(sequences[source][start as usize..end as usize].to_vec())
        })
        .unwrap();
        assert_eq!(typed[0].len(), 3);
        assert_eq!(
            typed[0][0].sequence,
            reverse_complement(&sequences[0][10..210])
        );
        assert_eq!(typed[0][1].sequence, sequences[1][20..220]);
        assert_eq!(
            typed[0][2].sequence,
            reverse_complement(&sequences[1][20..220])
        );

        let mut forward_left = typed[0][1].clone();
        forward_left.range.partition = 0;
        let mut reverse_left = typed[0][2].clone();
        reverse_left.range.partition = 0;
        let mut forward_right = typed[0][1].clone();
        forward_right.range.partition = 1;
        forward_right.range.occurrence = 2;
        let mut reverse_right = typed[0][2].clone();
        reverse_right.range.partition = 1;
        reverse_right.range.occurrence = 2;
        let linked = seams_from_public_links(
            &[
                vec![forward_left, reverse_left],
                vec![forward_right, reverse_right],
            ],
            &[catalog::Link {
                from: 1,
                to: 2,
                source_gap_bp: 0,
                relation: "both-orientations".into(),
                oriented_continuation: true,
            }],
            &panel,
            150,
            50_000,
        )
        .unwrap();
        assert_eq!(linked[0].len(), 2);
        assert!(linked[0]
            .iter()
            .any(|seam| !seam.left.reverse && !seam.right.reverse));
        assert!(linked[0]
            .iter()
            .any(|seam| seam.left.reverse && seam.right.reverse));
    }

    #[test]
    fn subwalk_accumulation_equals_weighted_bwt_orbit_counts() {
        let records = [
            (
                canonical(&encode_walk(&[(1, 0), (-1, 17), (1, 34)]).unwrap()),
                3,
            ),
            (
                canonical(&encode_walk(&[(2, 0), (2, 11), (2, 22), (2, 33)]).unwrap()),
                2,
            ),
        ];
        let stored = records.iter().cloned().collect::<BTreeMap<_, _>>();
        let bwt = WeightedBwt::build(&stored).unwrap();
        let mut direct = Profile::new();
        for (record, multiplicity) in records {
            for _ in 0..multiplicity {
                add_record_subwalks(&mut direct, &record, 1_000).unwrap();
            }
        }
        let singleton_counts = direct
            .iter()
            .filter(|(key, _)| key.len() == 1)
            .collect::<Vec<_>>();
        assert!(!singleton_counts.is_empty());
        for (key, &count) in singleton_counts {
            assert_eq!(bwt.count(key).unwrap(), count, "singleton {key:?}");
        }
        assert!(direct.keys().any(|key| key.len() > 3));
        for (key, count) in direct {
            assert_eq!(bwt.count(&key).unwrap(), count, "key {key:?}");
        }
        let long = canonical(&encode_walk(&[(3, 0), (4, 9)]).unwrap());
        assert_eq!(
            add_record_subwalks(&mut Profile::new(), &long, 1)
                .unwrap_err()
                .to_string(),
            "budget: maximal_mem_subwalk_features"
        );
    }

    #[test]
    fn five_token_linkage_survives_weighted_bwt_but_not_pair_projection() {
        let encoded = |nodes: [i32; 3]| {
            canonical(&encode_walk(&[(nodes[0], 0), (nodes[1], 10), (nodes[2], 20)]).unwrap())
        };
        // Both hypotheses have exactly the same singleton and adjacent-pair
        // marginals around the shared middle node. Only complete 5-token records
        // preserve which outer nodes co-occurred.
        let cis_records = [encoded([1, 2, 3]), encoded([4, 2, 5])];
        let trans_records = [encoded([1, 2, 5]), encoded([4, 2, 3])];
        let sample = WeightedBwt::build(
            &cis_records
                .iter()
                .cloned()
                .map(|record| (record, 1))
                .collect(),
        )
        .unwrap();
        let profile = |record: &FeatureKey| {
            let mut profile = Profile::new();
            add_record_subwalks(&mut profile, record, 100).unwrap();
            profile
        };
        let cis = cis_records.iter().map(profile).collect::<Vec<_>>();
        let trans = trans_records.iter().map(profile).collect::<Vec<_>>();
        let candidate_alleles = cis
            .iter()
            .chain(&trans)
            .enumerate()
            .map(|(index, profile)| Allele {
                range: range(0, index, index, 0, 150),
                sequence: vec![b'A'; 150],
                interior: profile.clone(),
            })
            .collect::<Vec<_>>();
        let universe = universe(&[candidate_alleles], &[], &sample, 100).unwrap();
        let cis_sum = sum_profiles(cis.iter(), 100).unwrap();
        let trans_sum = sum_profiles(trans.iter(), 100).unwrap();
        assert!(universe
            .entries
            .iter()
            .filter(|entry| entry.key.len() <= 3)
            .all(|entry| cis_sum.get(&entry.key) == trans_sum.get(&entry.key)));
        assert!(universe.entries.iter().any(|entry| {
            entry.key.len() == 5
                && cis_sum.get(&entry.key).copied().unwrap_or(0)
                    != trans_sum.get(&entry.key).copied().unwrap_or(0)
        }));
        let model = ScoreModel {
            read_length: 150,
            histogram: 1,
            denominator: 150.0,
            depth: 10.0,
            background: 0.1,
        };
        assert!(
            score_pair(&cis[0], &cis[1], &universe, &model).unwrap()
                < score_pair(&trans[0], &trans[1], &universe, &model).unwrap()
        );
    }

    #[test]
    fn four_catalog_partitions_retain_long_mem_boundary_phase() {
        const PART: usize = 180;
        const READ: usize = 150;
        const MAX: usize = 50_000;
        let a = [dna(PART, 11), dna(PART, 12), dna(PART, 13), dna(PART, 14)];
        let mut b = a.clone();
        for (partition, allele) in b.iter_mut().enumerate() {
            let (lo, hi) = match partition {
                1 => (PART - 42, PART - 18),
                2 => (18, 42),
                _ => (72, 96),
            };
            for base in &mut allele[lo..hi] {
                *base = match *base {
                    b'A' => b'C',
                    b'C' => b'G',
                    b'G' => b'T',
                    _ => b'A',
                };
            }
        }
        let hap = |mask: usize| {
            (0..4)
                .flat_map(|partition| {
                    if mask & (1 << partition) == 0 {
                        a[partition].clone()
                    } else {
                        b[partition].clone()
                    }
                })
                .collect::<Vec<_>>()
        };
        // All public local combinations exist, but the two truth mosaics are
        // absent as native whole paths.
        let source_sequences = (0usize..16)
            .filter(|mask| ![10, 12].contains(mask))
            .map(|mask| (format!("H{mask}#0#chr"), hap(mask)))
            .collect::<Vec<_>>();
        let b_source = source_sequences
            .iter()
            .position(|(name, _)| name == "H15#0#chr")
            .unwrap();
        let panel = SyngIndex::build(
            SyncmerParams::default(),
            source_sequences.clone().into_iter(),
        );
        let mut partitions = Vec::new();
        for partition in 0..4 {
            let start = partition * PART;
            let end = start + PART;
            partitions.push(vec![
                allele(
                    range(partition, 2 * partition, 0, start, end),
                    a[partition].clone(),
                    &panel,
                    READ,
                    MAX,
                )
                .unwrap(),
                allele(
                    range(partition, 2 * partition + 1, b_source, start, end),
                    b[partition].clone(),
                    &panel,
                    READ,
                    MAX,
                )
                .unwrap(),
            ]);
        }
        let mut public_links = Vec::new();
        for boundary in 0..3 {
            for left in &partitions[boundary] {
                for right in &partitions[boundary + 1] {
                    public_links.push(catalog::Link {
                        from: left.range.occurrence,
                        to: right.range.occurrence,
                        source_gap_bp: 0,
                        relation: "public-port-compatible".into(),
                        oriented_continuation: true,
                    });
                }
            }
        }
        let boundaries =
            seams_from_public_links(&partitions, &public_links, &panel, READ, MAX).unwrap();
        assert!(boundaries.iter().all(|boundary| boundary.len() == 4));

        let truth_sequences = [hap(12), hap(10)];
        assert!(truth_sequences
            .iter()
            .all(|truth| source_sequences.iter().all(|(_, native)| native != truth)));
        let reads = truth_sequences
            .iter()
            .flat_map(|sequence| {
                (0..=sequence.len() - READ).map(move |start| sequence[start..start + READ].to_vec())
            })
            .collect::<Vec<_>>();
        let sample = sample_bwt(&panel, &reads);
        let feature_ledger = universe(&partitions, &boundaries, &sample, MAX).unwrap();
        assert!(feature_ledger
            .entries
            .iter()
            .any(|entry| entry.key.len() > 3));
        assert_eq!(
            feature_ledger.entries.len(),
            feature_ledger
                .entries
                .iter()
                .map(|entry| &entry.key)
                .collect::<std::collections::BTreeSet<_>>()
                .len()
        );
        assert!(feature_ledger
            .entries
            .iter()
            .any(|entry| matches!(entry.class, FeatureClass::Owned(_))));
        assert!(feature_ledger
            .entries
            .iter()
            .any(|entry| matches!(entry.class, FeatureClass::Boundary(_))));
        assert!(feature_ledger.entries.iter().any(|entry| matches!(
            &entry.class,
            FeatureClass::Deferred(reason) if reason == "interior_and_boundary_contributions"
        )));

        let aa_chain = [
            &partitions[0][0],
            &partitions[1][0],
            &partitions[2][1],
            &partitions[3][1],
        ];
        let bb_chain = [
            &partitions[0][0],
            &partitions[1][1],
            &partitions[2][0],
            &partitions[3][1],
        ];
        let ab_chain = [
            &partitions[0][0],
            &partitions[1][0],
            &partitions[2][0],
            &partitions[3][1],
        ];
        let ba_chain = [
            &partitions[0][0],
            &partitions[1][1],
            &partitions[2][1],
            &partitions[3][1],
        ];
        let profile_chain = |chain: &[&Allele]| {
            profile_complete_sequence(&panel, &reconstruct_sequence(chain), READ, MAX).unwrap()
        };
        let ledger_chain = |chain: &[&Allele]| {
            let mut profiles = chain
                .iter()
                .map(|allele| &allele.interior)
                .collect::<Vec<_>>();
            for (boundary, adjacent) in chain.windows(2).enumerate() {
                profiles.push(
                    &boundaries[boundary]
                        .iter()
                        .find(|seam| {
                            seam.left == adjacent[0].range && seam.right == adjacent[1].range
                        })
                        .unwrap()
                        .profile,
                );
            }
            sum_profiles(profiles, MAX).unwrap()
        };
        let aa = profile_chain(&aa_chain);
        let bb = profile_chain(&bb_chain);
        let ab = profile_chain(&ab_chain);
        let ba = profile_chain(&ba_chain);
        assert_eq!(ledger_chain(&aa_chain), aa);
        assert_eq!(ledger_chain(&bb_chain), bb);
        assert_eq!(ledger_chain(&ab_chain), ab);
        assert_eq!(ledger_chain(&ba_chain), ba);
        let predicted = sum_profiles([&aa, &bb], MAX).unwrap();
        let spectrum = feature_spectrum(&feature_ledger, &predicted).unwrap();
        assert!(spectrum.iter().any(|bin| {
            bin.node_count == 1
                && bin.token_length == 1
                && bin.features > 0
                && bin.observed_multiplicity > 0
                && bin.predicted_multiplicity > 0
                && bin.local_features > 0
        }));
        assert!(spectrum
            .iter()
            .any(|bin| bin.node_count == 2 && bin.token_length == 3 && bin.features > 0));
        assert!(spectrum.iter().any(|bin| {
            bin.node_count >= 3
                && bin.token_length >= 5
                && bin.features > 0
                && bin.shared_features > 0
        }));
        let histogram = reads.len() as u64;
        let model = ScoreModel {
            read_length: READ as u64,
            histogram,
            denominator: (READ as u64 * histogram) as f64,
            depth: 150.0,
            background: 0.1,
        };
        let capped =
            chain_diplotypes_bounded(&partitions, &boundaries, &feature_ledger, &model, MAX, 1)
                .unwrap();
        assert!(!capped.complete);
        assert_eq!(
            capped.stop_reason.as_deref(),
            Some("budget: partition_work")
        );
        let local_work = partitions
            .iter()
            .enumerate()
            .map(|(partition, alleles)| {
                (alleles.len() * (alleles.len() + 1) / 2) as u64
                    * class_terms(&feature_ledger, &FeatureClass::Owned(partition))
            })
            .sum::<u64>();
        let mut mid_chain = chain_diplotypes_bounded(
            &partitions,
            &boundaries,
            &feature_ledger,
            &model,
            MAX,
            local_work + 1,
        )
        .unwrap();
        assert!(!mid_chain.complete);
        assert_eq!(
            mid_chain.stop_reason.as_deref(),
            Some("budget: partition_work")
        );
        assert!(!mid_chain.states.is_empty());
        // Rescoring a cap-stopped partial-depth chain reports an incomplete
        // search instead of a hard error.
        let incomplete_rescore = rescore_chained(
            &mut mid_chain,
            &panel,
            &partitions,
            &feature_ledger,
            &model,
            &[],
            MAX,
        )
        .unwrap();
        assert!(!mid_chain.complete);
        assert_eq!(
            mid_chain.stop_reason.as_deref(),
            Some("incomplete chained reconstruction")
        );
        assert!(incomplete_rescore.is_empty());
        let cis = score_pair(&aa, &bb, &feature_ledger, &model).unwrap();
        let trans = score_pair(&ab, &ba, &feature_ledger, &model).unwrap();
        let independent = feature_ledger
            .entries
            .iter()
            .map(|entry| {
                let q = aa.get(&entry.key).copied().unwrap_or(0)
                    + bb.get(&entry.key).copied().unwrap_or(0);
                let signal = q as f64 * model.depth / model.read_length as f64;
                signal - entry.observed as f64 * (signal / model.background).ln_1p()
            })
            .sum::<f64>();
        assert!((cis - independent).abs() < 1e-9);
        let separately_summed_haploid =
            score_where(&aa, &Profile::new(), &feature_ledger, &model, |_| true).unwrap()
                + score_where(&bb, &Profile::new(), &feature_ledger, &model, |_| true).unwrap();
        assert!((cis - separately_summed_haploid).abs() > 1e-6);
        assert!(
            cis < trans,
            "linked truth {cis} must beat alternate phase {trans}"
        );
        assert!(feature_ledger.entries.iter().any(|entry| {
            entry.key.len() > 3
                && aa.get(&entry.key).copied().unwrap_or(0)
                    + bb.get(&entry.key).copied().unwrap_or(0)
                    != ab.get(&entry.key).copied().unwrap_or(0)
                        + ba.get(&entry.key).copied().unwrap_or(0)
        }));
        assert!(feature_ledger
            .entries
            .iter()
            .any(|entry| entry.observed == 0));
        assert_eq!(
            local_diplotypes(1, &partitions[1], &feature_ledger, &model)
                .unwrap()
                .len(),
            3
        );

        let left = Diplotype {
            alleles: [0, 1],
            loss: 0.0,
        };
        let right = Diplotype {
            alleles: [0, 1],
            loss: 0.0,
        };
        let matchings = boundary_matchings(
            &left,
            &right,
            &partitions[1],
            &partitions[2],
            &boundaries[1],
            MAX,
        )
        .unwrap();
        assert_eq!(matchings.len(), 2);
        assert_ne!(matchings[0].0, matchings[1].0);
        assert_ne!(matchings[0].1, matchings[1].1);
        assert_eq!(reconstruct_route(&aa_chain).segments.len(), 4);
        assert_eq!(reconstruct_sequence(&aa_chain), truth_sequences[0]);

        let chained = chain_diplotypes(&partitions, &boundaries, &feature_ledger, &model).unwrap();
        assert!(chained.complete);
        assert!(!chained.states.is_empty());
        assert_eq!(chained.accounting.retained_states, 255);
        assert_eq!(chained.accounting.matchings, 252);
        let mut chained_pairs = chained
            .states
            .iter()
            .map(|state| {
                let profiles = state.choices.each_ref().map(|choices| {
                    let chain = choices
                        .iter()
                        .enumerate()
                        .map(|(partition, &allele)| &partitions[partition][allele])
                        .collect::<Vec<_>>();
                    profile_chain(&chain)
                });
                (
                    score_pair(&profiles[0], &profiles[1], &feature_ledger, &model).unwrap(),
                    state,
                )
            })
            .collect::<Vec<_>>();
        chained_pairs.sort_by(|a, b| a.0.total_cmp(&b.0));
        assert_eq!(
            chained_pairs[0].1.choices,
            [vec![0, 0, 1, 1], vec![0, 1, 0, 1]]
        );
        let dosage = (0..4)
            .map(|partition| {
                chained_pairs[0].1.choices[0][partition] + chained_pairs[0].1.choices[1][partition]
            })
            .collect::<Vec<_>>();
        assert_eq!(dosage, vec![0, 1, 1, 2]);
        assert!(chained.accounting.work > 0);
        assert!(chained.accounting.profile_work > 0);
        assert!(chained.accounting.state_bytes <= 134_217_728);
        assert!(chained.stop_reason.is_none());
        eprintln!(
            "linked_partition_accounting={}",
            serde_json::to_string(&chained.accounting).unwrap()
        );
        assert!((chained_pairs[0].0 - cis).abs() < 1e-9);
        assert!(
            (chained_pairs[0].1.proposal_loss - chained_pairs[0].0).abs() > 1e-6,
            "deferred terms must be absent locally and charged once globally"
        );

        // Remove only central spanning evidence. Both phase classes remain
        // physically distinct, have identical retained count vectors, and tie.
        let mut unlinked_boundaries = boundaries.clone();
        for seam in &mut unlinked_boundaries[1] {
            seam.profile.clear();
        }
        let unlinked_universe = universe(&partitions, &unlinked_boundaries, &sample, MAX).unwrap();
        let unlinked = chain_diplotypes(
            &partitions,
            &unlinked_boundaries,
            &unlinked_universe,
            &model,
        )
        .unwrap();
        let phase = |state: &ChainState| {
            let mut copies = state.choices.clone();
            copies.sort();
            copies
        };
        let truth_phase = [vec![0, 0, 1, 1], vec![0, 1, 0, 1]];
        let alternate_phase = [vec![0, 0, 0, 1], vec![0, 1, 1, 1]];
        let truth_state = unlinked
            .states
            .iter()
            .find(|state| phase(state) == truth_phase)
            .unwrap();
        let alternate_state = unlinked
            .states
            .iter()
            .find(|state| phase(state) == alternate_phase)
            .unwrap();
        assert!((truth_state.proposal_loss - alternate_state.proposal_loss).abs() < 1e-9);
        let unlinked_profile = |state: &ChainState| {
            let copy = |slot: usize| {
                let chain = state.choices[slot]
                    .iter()
                    .enumerate()
                    .map(|(partition, &index)| &partitions[partition][index])
                    .collect::<Vec<_>>();
                let mut pieces = chain
                    .iter()
                    .map(|allele| &allele.interior)
                    .collect::<Vec<_>>();
                for (boundary, adjacent) in chain.windows(2).enumerate() {
                    pieces.push(
                        &unlinked_boundaries[boundary]
                            .iter()
                            .find(|seam| {
                                seam.left == adjacent[0].range && seam.right == adjacent[1].range
                            })
                            .unwrap()
                            .profile,
                    );
                }
                sum_profiles(pieces, MAX).unwrap()
            };
            let first = copy(0);
            let second = copy(1);
            sum_profiles([&first, &second], MAX).unwrap()
        };
        assert_eq!(
            unlinked_profile(truth_state),
            unlinked_profile(alternate_state)
        );

        // Public-link removal makes that continuation unavailable rather than
        // being repaired by a Cartesian fallback.
        let illegal_links = public_links
            .iter()
            .filter(|link| {
                !(link.from == partitions[1][0].range.occurrence
                    && link.to == partitions[2][1].range.occurrence)
            })
            .map(|link| catalog::Link {
                from: link.from,
                to: link.to,
                source_gap_bp: link.source_gap_bp,
                relation: link.relation.clone(),
                oriented_continuation: link.oriented_continuation,
            })
            .collect::<Vec<_>>();
        let legal_only =
            seams_from_public_links(&partitions, &illegal_links, &panel, READ, MAX).unwrap();
        let legal_search = chain_diplotypes(
            &partitions,
            &legal_only,
            &universe(&partitions, &legal_only, &sample, MAX).unwrap(),
            &model,
        )
        .unwrap();
        assert!(legal_search.states.iter().all(|state| {
            (0..2).all(|copy| !(state.choices[copy][1] == 0 && state.choices[copy][2] == 1))
        }));

        let haplotypes =
            enumerate_haplotypes_assessment_oracle(&panel, &partitions, &boundaries, READ, MAX, 32)
                .unwrap();
        assert_eq!(haplotypes.len(), 16);
        assert!(haplotypes.iter().all(|haplotype| {
            haplotype.ranges.len() == 4
                && haplotype.route.segments.len() == 4
                && haplotype.sequence.len() == 4 * PART
        }));
        let ranked =
            rank_pairs_assessment_oracle(&haplotypes, &feature_ledger, &model, 256).unwrap();
        let best_choices = ranked[0]
            .haplotypes
            .map(|index| haplotypes[index].choices.clone());
        assert_eq!(best_choices, [vec![0, 0, 1, 1], vec![0, 1, 0, 1]]);
        assert!((ranked[0].objective - cis).abs() < 1e-9);
        assert!(rank_pairs_assessment_oracle(&haplotypes, &feature_ledger, &model, 135).is_err());
    }
}
