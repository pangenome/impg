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

fn feature_universe_with<'a, T: 'a, I>(candidate_features: I) -> Vec<u32>
where
    I: IntoIterator<Item = &'a [(u32, T)]>,
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

fn sample_norm_sq_with<T: Copy>(
    sample: &FxHashMap<u32, T>,
    features: &[u32],
    to_f64: impl Fn(T) -> f64,
) -> f64 {
    features
        .iter()
        .filter_map(|feature_id| sample.get(feature_id).copied())
        .map(|v| {
            let x = to_f64(v);
            x * x
        })
        .sum()
}

fn cosine_for_features_with<T: Copy>(
    candidate_features: &[(u32, T)],
    sample: &FxHashMap<u32, T>,
    sample_norm_sq: f64,
    to_f64: impl Fn(T) -> f64,
) -> f64 {
    if sample_norm_sq <= 0.0 || !sample_norm_sq.is_finite() {
        return 0.0;
    }
    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for &(feature_id, val) in candidate_features {
        let g = to_f64(val);
        genotype_norm_sq += g * g;
        if let Some(&s) = sample.get(&feature_id) {
            dot += g * to_f64(s);
        }
    }
    if genotype_norm_sq <= 0.0 || !genotype_norm_sq.is_finite() || !dot.is_finite() {
        return 0.0;
    }
    let sim = dot / (sample_norm_sq.sqrt() * genotype_norm_sq.sqrt());
    if sim.is_finite() {
        sim
    } else {
        0.0
    }
}

fn finalize_combination_score(
    combination: &[usize],
    dot: f64,
    sample_norm: f64,
    genotype_norm_sq: f64,
) -> CombinationScore {
    let genotype_norm = if genotype_norm_sq > 0.0 && genotype_norm_sq.is_finite() {
        genotype_norm_sq.sqrt()
    } else {
        0.0
    };
    let similarity = if sample_norm <= 0.0
        || !sample_norm.is_finite()
        || genotype_norm <= 0.0
        || !dot.is_finite()
    {
        0.0
    } else {
        let sim = dot / (sample_norm * genotype_norm);
        if sim.is_finite() {
            sim
        } else {
            0.0
        }
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

pub fn feature_universe<'a, I>(candidate_features: I) -> Vec<u32>
where
    I: IntoIterator<Item = &'a [(u32, u64)]>,
{
    feature_universe_with(candidate_features)
}

pub fn sample_norm_sq_for_features(sample_counts: &FxHashMap<u32, u64>, features: &[u32]) -> f64 {
    sample_norm_sq_with(sample_counts, features, |v| v as f64)
}

pub fn cosine_for_feature_counts(
    candidate_features: &[(u32, u64)],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
) -> f64 {
    cosine_for_features_with(candidate_features, sample_counts, sample_norm_sq, |v| {
        v as f64
    })
}

pub fn feature_universe_f64<'a, I>(candidate_features: I) -> Vec<u32>
where
    I: IntoIterator<Item = &'a [(u32, f64)]>,
{
    feature_universe_with(candidate_features)
}

pub fn sample_norm_sq_for_weight_features(
    sample_weights: &FxHashMap<u32, f64>,
    features: &[u32],
) -> f64 {
    sample_norm_sq_with(sample_weights, features, |v| v)
}

pub fn cosine_for_feature_weights(
    candidate_features: &[(u32, f64)],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
) -> f64 {
    cosine_for_features_with(candidate_features, sample_weights, sample_norm_sq, |v| v)
}

trait FeatureWeight: Copy + Default {
    fn add_saturating(self, rhs: Self) -> Self;
    fn into_f64(self) -> f64;
}

impl FeatureWeight for u64 {
    #[inline]
    fn add_saturating(self, rhs: Self) -> Self {
        self.saturating_add(rhs)
    }

    #[inline]
    fn into_f64(self) -> f64 {
        self as f64
    }
}

impl FeatureWeight for f64 {
    #[inline]
    fn add_saturating(self, rhs: Self) -> Self {
        self + rhs
    }

    #[inline]
    fn into_f64(self) -> f64 {
        self
    }
}

fn score_combination_with<W: FeatureWeight>(
    combination: &[usize],
    candidate_features: &[&[(u32, W)]],
    sample: &FxHashMap<u32, W>,
    sample_norm_sq: f64,
) -> CombinationScore {
    let mut genotype: FxHashMap<u32, W> = FxHashMap::default();
    for &idx in combination {
        for &(feature_id, value) in candidate_features[idx] {
            let entry = genotype.entry(feature_id).or_default();
            *entry = entry.add_saturating(value);
        }
    }

    let mut dot = 0.0;
    let mut genotype_norm_sq = 0.0;
    for (feature_id, value) in genotype {
        let g = value.into_f64();
        genotype_norm_sq += g * g;
        if let Some(&s) = sample.get(&feature_id) {
            dot += g * s.into_f64();
        }
    }

    let sample_norm = if sample_norm_sq > 0.0 && sample_norm_sq.is_finite() {
        sample_norm_sq.sqrt()
    } else {
        0.0
    };
    finalize_combination_score(combination, dot, sample_norm, genotype_norm_sq)
}

pub fn score_cosine_combination(
    combination: &[usize],
    candidate_features: &[&[(u32, u64)]],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
) -> CombinationScore {
    score_combination_with(combination, candidate_features, sample_counts, sample_norm_sq)
}

pub fn score_cosine_combination_f64(
    combination: &[usize],
    candidate_features: &[&[(u32, f64)]],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
) -> CombinationScore {
    score_combination_with(combination, candidate_features, sample_weights, sample_norm_sq)
}

/// Dense-indexed accumulator avoiding per-combination hash map allocations during DFS.
fn run_dense_cosine_search<W: FeatureWeight>(
    candidate_features: &[&[(u32, W)]],
    sample: &FxHashMap<u32, W>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
) -> io::Result<Vec<CombinationScore>> {
    let mut feature_to_dense: FxHashMap<u32, usize> = FxHashMap::default();
    let mut dense_candidates: Vec<Vec<(usize, W)>> = Vec::with_capacity(candidate_features.len());
    for &features in candidate_features {
        let mut mapped = Vec::with_capacity(features.len());
        for &(feature_id, value) in features {
            let next_idx = feature_to_dense.len();
            let dense_idx = *feature_to_dense.entry(feature_id).or_insert(next_idx);
            mapped.push((dense_idx, value));
        }
        dense_candidates.push(mapped);
    }

    let mut dense_sample = vec![0.0f64; feature_to_dense.len()];
    for (&feature_id, &dense_idx) in &feature_to_dense {
        if let Some(&val) = sample.get(&feature_id) {
            dense_sample[dense_idx] = val.into_f64();
        }
    }

    let sample_norm = if sample_norm_sq > 0.0 && sample_norm_sq.is_finite() {
        sample_norm_sq.sqrt()
    } else {
        0.0
    };

    struct Search<'a, W: FeatureWeight> {
        dense_candidates: &'a [Vec<(usize, W)>],
        dense_sample: &'a [f64],
        sample_norm: f64,
        ploidy: usize,
        max_combinations: u64,
        visited: u64,
        acc: Vec<W>,
        seen_gen: Vec<u32>,
        touched: Vec<usize>,
        cur_gen: u32,
        results: Vec<CombinationScore>,
    }

    impl<W: FeatureWeight> Search<'_, W> {
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
                self.cur_gen = self.cur_gen.wrapping_add(1);
                if self.cur_gen == 0 {
                    self.seen_gen.fill(0);
                    self.cur_gen = 1;
                }
                self.touched.clear();
                for &idx in current.iter() {
                    for &(feature_idx, value) in &self.dense_candidates[idx] {
                        if self.seen_gen[feature_idx] != self.cur_gen {
                            self.seen_gen[feature_idx] = self.cur_gen;
                            self.acc[feature_idx] = value;
                            self.touched.push(feature_idx);
                        } else {
                            self.acc[feature_idx] = self.acc[feature_idx].add_saturating(value);
                        }
                    }
                }
                let mut dot = 0.0;
                let mut genotype_norm_sq = 0.0;
                for &feature_idx in &self.touched {
                    let g = self.acc[feature_idx].into_f64();
                    genotype_norm_sq += g * g;
                    dot += g * self.dense_sample[feature_idx];
                }
                self.results.push(finalize_combination_score(
                    current,
                    dot,
                    self.sample_norm,
                    genotype_norm_sq,
                ));
                return Ok(());
            }

            for idx in start..self.dense_candidates.len() {
                current.push(idx);
                self.visit(idx, current)?;
                current.pop();
            }
            Ok(())
        }
    }

    let num_features = dense_sample.len();
    let mut search = Search {
        dense_candidates: &dense_candidates,
        dense_sample: &dense_sample,
        sample_norm,
        ploidy,
        max_combinations,
        visited: 0,
        acc: vec![W::default(); num_features],
        seen_gen: vec![0u32; num_features],
        touched: Vec::new(),
        cur_gen: 0,
        results: Vec::new(),
    };
    let mut current = Vec::with_capacity(ploidy);
    search.visit(0, &mut current)?;
    search.results.sort_by(|a, b| {
        b.similarity
            .total_cmp(&a.similarity)
            .then_with(|| b.dot.total_cmp(&a.dot))
            .then_with(|| a.combination.cmp(&b.combination))
    });
    Ok(search.results)
}

pub fn run_cosine_combination_search(
    candidate_features: &[&[(u32, u64)]],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
) -> io::Result<Vec<CombinationScore>> {
    run_dense_cosine_search(
        candidate_features,
        sample_counts,
        sample_norm_sq,
        ploidy,
        max_combinations,
    )
}

pub fn run_cosine_combination_search_f64(
    candidate_features: &[&[(u32, f64)]],
    sample_weights: &FxHashMap<u32, f64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
) -> io::Result<Vec<CombinationScore>> {
    run_dense_cosine_search(
        candidate_features,
        sample_weights,
        sample_norm_sq,
        ploidy,
        max_combinations,
    )
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

    #[test]
    fn cosine_guards_zero_norm_and_non_finite_weights() {
        let candidate_a = vec![(1, 1.0), (2, 2.0)];
        let candidate_b = vec![(1, 0.0)];
        let candidates = vec![candidate_a.as_slice(), candidate_b.as_slice()];
        let mut nan_sample = FxHashMap::default();
        nan_sample.insert(1, f64::NAN);
        assert_eq!(
            cosine_for_feature_weights(&candidate_a, &nan_sample, 1.0),
            0.0
        );
        let zero_score = score_cosine_combination_f64(&[0], &candidates, &nan_sample, f64::NAN);
        assert_eq!(zero_score.similarity, 0.0);
        assert_eq!(zero_score.qv, 0.0);

        let empty_sample = FxHashMap::default();
        let results =
            run_cosine_combination_search_f64(&candidates, &empty_sample, 0.0, 2, 10).unwrap();
        assert_eq!(results.len(), 3);
        assert!(results.iter().all(|r| r.similarity == 0.0 && r.qv == 0.0));
    }
}
