//! Shared vocabulary for backend-neutral genotyping and inference.
//!
//! The current implementation is syng-backed, but the genotyping pipeline is
//! intended to support multiple graph backends, evidence projections, and
//! scoring methods. These enums keep user-visible metadata stable while giving
//! future backends a common language.

use rustc_hash::FxHashMap;
use std::io;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FeatureSpace {
    SyngSyncmerNode,
    GfaSegment,
    VariationGraphNode,
    LocalHaplotypeSegment,
    HaplotypeSequence,
    MemHit,
}

impl FeatureSpace {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::SyngSyncmerNode => "syng-syncmer-node",
            Self::GfaSegment => "gfa-segment",
            Self::VariationGraphNode => "variation-graph-node",
            Self::LocalHaplotypeSegment => "local-haplotype-segment",
            Self::HaplotypeSequence => "haplotype-sequence",
            Self::MemHit => "mem-hit",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EvidenceBackend {
    Pack,
    Projection,
    SyncmerWalk,
    VariationGraphAlignment,
    HaplotypeAlignment,
    MemProjection,
}

impl EvidenceBackend {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pack => "pack",
            Self::Projection => "projection",
            Self::SyncmerWalk => "syncmer-walk",
            Self::VariationGraphAlignment => "variation-graph-alignment",
            Self::HaplotypeAlignment => "haplotype-alignment",
            Self::MemProjection => "mem-projection",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ScoringMethod {
    Cos,
    CountLikelihood,
    AlignmentLikelihood,
    LearnedLocalEmission,
    GbwtImputation,
}

impl ScoringMethod {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Cos => "cos",
            Self::CountLikelihood => "count-likelihood",
            Self::AlignmentLikelihood => "alignment-likelihood",
            Self::LearnedLocalEmission => "learned-local-emission",
            Self::GbwtImputation => "gbwt-imputation",
        }
    }

    pub fn metric_str(self) -> &'static str {
        match self {
            Self::Cos => "cosine",
            Self::CountLikelihood => "count-likelihood",
            Self::AlignmentLikelihood => "alignment-likelihood",
            Self::LearnedLocalEmission => "learned-local-emission",
            Self::GbwtImputation => "gbwt-imputation",
        }
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct CombinationScore {
    pub combination: Vec<usize>,
    pub similarity: f64,
    pub qv: f64,
    pub dot: f64,
    pub sample_norm: f64,
    pub genotype_norm: f64,
}

pub fn feature_universe<'a, I>(candidate_features: I) -> Vec<u32>
where
    I: IntoIterator<Item = &'a [(u32, u64)]>,
{
    let mut seen: FxHashMap<u32, ()> = FxHashMap::default();
    for features in candidate_features {
        for &(feature_id, _) in features {
            seen.insert(feature_id, ());
        }
    }
    let mut features: Vec<u32> = seen.into_keys().collect();
    features.sort_unstable();
    features
}

pub fn sample_norm_sq_for_features(sample_counts: &FxHashMap<u32, u64>, features: &[u32]) -> f64 {
    features
        .iter()
        .map(|feature_id| sample_counts.get(feature_id).copied().unwrap_or(0) as f64)
        .map(|v| v * v)
        .sum()
}

pub fn cosine_for_feature_counts(
    candidate_features: &[(u32, u64)],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
) -> f64 {
    if sample_norm_sq == 0.0 {
        return 0.0;
    }
    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for &(feature_id, count) in candidate_features {
        let g = count as f64;
        genotype_norm_sq += g * g;
        dot += g * sample_counts.get(&feature_id).copied().unwrap_or(0) as f64;
    }
    if genotype_norm_sq == 0.0 {
        0.0
    } else {
        dot / (sample_norm_sq.sqrt() * genotype_norm_sq.sqrt())
    }
}

pub fn feature_universe_f64<'a, I>(candidate_features: I) -> Vec<u32>
where
    I: IntoIterator<Item = &'a [(u32, f64)]>,
{
    let mut seen: FxHashMap<u32, ()> = FxHashMap::default();
    for features in candidate_features {
        for &(feature_id, _) in features {
            seen.insert(feature_id, ());
        }
    }
    let mut features: Vec<u32> = seen.into_keys().collect();
    features.sort_unstable();
    features
}

pub fn sample_norm_sq_for_weight_features(
    sample_weights: &FxHashMap<u32, f64>,
    features: &[u32],
) -> f64 {
    features
        .iter()
        .map(|feature_id| sample_weights.get(feature_id).copied().unwrap_or(0.0))
        .map(|v| v * v)
        .sum()
}

pub fn cosine_for_feature_weights(
    candidate_features: &[(u32, f64)],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
) -> f64 {
    if sample_norm_sq == 0.0 {
        return 0.0;
    }
    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for &(feature_id, weight) in candidate_features {
        genotype_norm_sq += weight * weight;
        dot += weight * sample_weights.get(&feature_id).copied().unwrap_or(0.0);
    }
    if genotype_norm_sq == 0.0 {
        0.0
    } else {
        dot / (sample_norm_sq.sqrt() * genotype_norm_sq.sqrt())
    }
}

pub fn score_cosine_combination(
    combination: &[usize],
    candidate_features: &[&[(u32, u64)]],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
) -> CombinationScore {
    let mut genotype_counts: FxHashMap<u32, u64> = FxHashMap::default();
    for &idx in combination {
        for &(feature_id, count) in candidate_features[idx] {
            *genotype_counts.entry(feature_id).or_insert(0) += count;
        }
    }

    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for (feature_id, count) in genotype_counts {
        let g = count as f64;
        genotype_norm_sq += g * g;
        dot += g * sample_counts.get(&feature_id).copied().unwrap_or(0) as f64;
    }

    let sample_norm = sample_norm_sq.sqrt();
    let genotype_norm = genotype_norm_sq.sqrt();
    let similarity = if sample_norm == 0.0 || genotype_norm == 0.0 {
        0.0
    } else {
        dot / (sample_norm * genotype_norm)
    };
    let qv = if similarity >= 1.0 {
        999.0
    } else if similarity <= 0.0 {
        0.0
    } else {
        -10.0 * (1.0 - similarity).log10()
    };

    CombinationScore {
        combination: combination.to_vec(),
        similarity,
        qv,
        dot,
        sample_norm,
        genotype_norm,
    }
}

pub fn score_cosine_combination_f64(
    combination: &[usize],
    candidate_features: &[&[(u32, f64)]],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
) -> CombinationScore {
    let mut genotype_weights: FxHashMap<u32, f64> = FxHashMap::default();
    for &idx in combination {
        for &(feature_id, weight) in candidate_features[idx] {
            *genotype_weights.entry(feature_id).or_insert(0.0) += weight;
        }
    }

    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for (feature_id, weight) in genotype_weights {
        genotype_norm_sq += weight * weight;
        dot += weight * sample_weights.get(&feature_id).copied().unwrap_or(0.0);
    }

    let sample_norm = sample_norm_sq.sqrt();
    let genotype_norm = genotype_norm_sq.sqrt();
    let similarity = if sample_norm == 0.0 || genotype_norm == 0.0 {
        0.0
    } else {
        dot / (sample_norm * genotype_norm)
    };
    let qv = if similarity >= 1.0 {
        999.0
    } else if similarity <= 0.0 {
        0.0
    } else {
        -10.0 * (1.0 - similarity).log10()
    };

    CombinationScore {
        combination: combination.to_vec(),
        similarity,
        qv,
        dot,
        sample_norm,
        genotype_norm,
    }
}

struct CosineCombinationSearch<'a> {
    candidate_features: &'a [&'a [(u32, u64)]],
    sample_counts: &'a FxHashMap<u32, u64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
    visited: u64,
    results: Vec<CombinationScore>,
}

impl CosineCombinationSearch<'_> {
    fn run(mut self) -> io::Result<Vec<CombinationScore>> {
        let mut current = Vec::with_capacity(self.ploidy);
        self.visit(0, &mut current)?;
        self.results.sort_by(|a, b| {
            b.similarity
                .total_cmp(&a.similarity)
                .then_with(|| b.dot.total_cmp(&a.dot))
                .then_with(|| a.combination.cmp(&b.combination))
        });
        Ok(self.results)
    }

    fn visit(&mut self, start: usize, current: &mut Vec<usize>) -> io::Result<()> {
        if current.len() == self.ploidy {
            self.visited += 1;
            if self.visited > self.max_combinations {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "genotype combination search exceeded --max-combinations ({})",
                        self.max_combinations
                    ),
                ));
            }
            self.results.push(score_cosine_combination(
                current,
                self.candidate_features,
                self.sample_counts,
                self.sample_norm_sq,
            ));
            return Ok(());
        }

        for idx in start..self.candidate_features.len() {
            current.push(idx);
            self.visit(idx, current)?;
            current.pop();
        }
        Ok(())
    }
}

pub fn run_cosine_combination_search(
    candidate_features: &[&[(u32, u64)]],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
) -> io::Result<Vec<CombinationScore>> {
    CosineCombinationSearch {
        candidate_features,
        sample_counts,
        sample_norm_sq,
        ploidy,
        max_combinations,
        visited: 0,
        results: Vec::new(),
    }
    .run()
}

struct CosineCombinationSearchF64<'a> {
    candidate_features: &'a [&'a [(u32, f64)]],
    sample_weights: &'a FxHashMap<u32, f64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
    visited: u64,
    results: Vec<CombinationScore>,
}

impl CosineCombinationSearchF64<'_> {
    fn run(mut self) -> io::Result<Vec<CombinationScore>> {
        let mut current = Vec::with_capacity(self.ploidy);
        self.visit(0, &mut current)?;
        self.results.sort_by(|a, b| {
            b.similarity
                .total_cmp(&a.similarity)
                .then_with(|| b.dot.total_cmp(&a.dot))
                .then_with(|| a.combination.cmp(&b.combination))
        });
        Ok(self.results)
    }

    fn visit(&mut self, start: usize, current: &mut Vec<usize>) -> io::Result<()> {
        if current.len() == self.ploidy {
            self.visited += 1;
            if self.visited > self.max_combinations {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "genotype combination search exceeded --max-combinations ({})",
                        self.max_combinations
                    ),
                ));
            }
            self.results.push(score_cosine_combination_f64(
                current,
                self.candidate_features,
                self.sample_weights,
                self.sample_norm_sq,
            ));
            return Ok(());
        }

        for idx in start..self.candidate_features.len() {
            current.push(idx);
            self.visit(idx, current)?;
            current.pop();
        }
        Ok(())
    }
}

pub fn run_cosine_combination_search_f64(
    candidate_features: &[&[(u32, f64)]],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
) -> io::Result<Vec<CombinationScore>> {
    CosineCombinationSearchF64 {
        candidate_features,
        sample_weights,
        sample_norm_sq,
        ploidy,
        max_combinations,
        visited: 0,
        results: Vec::new(),
    }
    .run()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn stable_user_visible_names() {
        assert_eq!(FeatureSpace::SyngSyncmerNode.as_str(), "syng-syncmer-node");
        assert_eq!(FeatureSpace::GfaSegment.as_str(), "gfa-segment");
        assert_eq!(EvidenceBackend::Pack.as_str(), "pack");
        assert_eq!(ScoringMethod::Cos.as_str(), "cos");
        assert_eq!(ScoringMethod::Cos.metric_str(), "cosine");
    }

    #[test]
    fn cosine_combination_search_ranks_with_replacement_genotypes() {
        let candidate_a = vec![(1, 1)];
        let candidate_b = vec![(2, 1)];
        let candidates = vec![candidate_a.as_slice(), candidate_b.as_slice()];
        let mut sample = FxHashMap::default();
        sample.insert(1, 2);
        let features = feature_universe(candidates.iter().copied());
        assert_eq!(features, vec![1, 2]);
        let sample_norm_sq = sample_norm_sq_for_features(&sample, &features);

        let results =
            run_cosine_combination_search(&candidates, &sample, sample_norm_sq, 2, 10).unwrap();

        assert_eq!(results[0].combination, vec![0, 0]);
        assert!((results[0].similarity - 1.0).abs() < 1e-9);
        assert!(results[1].similarity < results[0].similarity);
    }

    #[test]
    fn cosine_combination_search_enforces_budget() {
        let candidate_a = vec![(1, 1)];
        let candidate_b = vec![(2, 1)];
        let candidates = vec![candidate_a.as_slice(), candidate_b.as_slice()];
        let mut sample = FxHashMap::default();
        sample.insert(1, 1);
        sample.insert(2, 1);
        let features = feature_universe(candidates.iter().copied());
        let sample_norm_sq = sample_norm_sq_for_features(&sample, &features);

        let err =
            run_cosine_combination_search(&candidates, &sample, sample_norm_sq, 2, 2).unwrap_err();

        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(err
            .to_string()
            .contains("genotype combination search exceeded --max-combinations"));
    }

    #[test]
    fn cosine_combination_search_f64_ranks_weighted_genotypes() {
        let candidate_a = vec![(1, 0.5), (2, 1.0)];
        let candidate_b = vec![(1, 1.0)];
        let candidates = vec![candidate_a.as_slice(), candidate_b.as_slice()];
        let mut sample = FxHashMap::default();
        sample.insert(1, 1.0);
        sample.insert(2, 1.0);
        let features = feature_universe_f64(candidates.iter().copied());
        let sample_norm_sq = sample_norm_sq_for_weight_features(&sample, &features);

        let results =
            run_cosine_combination_search_f64(&candidates, &sample, sample_norm_sq, 1, 10).unwrap();

        assert_eq!(results[0].combination, vec![0]);
        assert!(results[0].similarity > results[1].similarity);
    }
}
