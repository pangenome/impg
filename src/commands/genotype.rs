use clap::ValueEnum;
use log::info;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::io::{self, Write};

use crate::commands::partition;
use crate::pack;
use crate::syng;

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum CandidateMode {
    /// Keep path candidates whose shared anchors span the requested reference range.
    Spanning,
    /// Score every query-gathered candidate interval independently.
    Overlapping,
}

#[derive(Debug, Clone)]
pub struct HaplotypeCandidate {
    pub path_name: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub anchors: usize,
    pub query_span_fraction: f64,
    pub features: Vec<(u32, u64)>,
    pub oriented_walk: Vec<i32>,
    single_similarity: f64,
}

struct CandidateFeatures {
    unoriented: Vec<(u32, u64)>,
    oriented_walk: Vec<i32>,
}

#[derive(Debug)]
struct SyngCandidateGroup {
    start: u64,
    end: u64,
    anchors: Vec<syng::Anchor>,
}

#[derive(Debug)]
pub struct CosigtResult {
    pub combination: Vec<usize>,
    pub similarity: f64,
    pub qv: f64,
    pub dot: f64,
    pub sample_norm: f64,
    pub genotype_norm: f64,
}

pub struct SyngCosigtConfig<'a> {
    pub syng_prefix: &'a str,
    pub pack_path: &'a str,
    pub target_range: &'a str,
    pub candidate_mode: CandidateMode,
    pub ploidy: usize,
    pub top_n: usize,
    pub candidate_top_k: usize,
    pub max_combinations: u64,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub min_anchors: usize,
    pub min_span_fraction: f64,
}

pub struct SyngCosigtQuery<'a> {
    pub target_range: &'a str,
    pub candidate_mode: CandidateMode,
    pub ploidy: usize,
    pub top_n: usize,
    pub candidate_top_k: usize,
    pub max_combinations: u64,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub min_anchors: usize,
    pub min_span_fraction: f64,
}

#[derive(Debug)]
pub struct SyngCosigtOutput {
    pub region_name: String,
    pub candidate_mode: CandidateMode,
    pub ploidy: usize,
    pub candidates: Vec<HaplotypeCandidate>,
    pub selected_features: Vec<u32>,
    pub pack_nonzero_nodes: u64,
    pub pack_universe_nodes: Option<u64>,
    pub pack_retained_records: Option<u64>,
    pub pack_syncmer_anchors: Option<u64>,
    pub results: Vec<CosigtResult>,
}

fn anchor_span_fraction(
    anchors: &[syng::Anchor],
    region_start: u64,
    region_end: u64,
    syncmer_len: u64,
) -> f64 {
    if anchors.is_empty() || region_end <= region_start {
        return 0.0;
    }
    let anchor_start = anchors
        .iter()
        .map(|a| a.query_pos)
        .min()
        .unwrap_or(region_start);
    let anchor_end = anchors
        .iter()
        .map(|a| a.query_pos.saturating_add(syncmer_len))
        .max()
        .unwrap_or(anchor_start);
    let clipped_start = anchor_start.max(region_start);
    let clipped_end = anchor_end.min(region_end);
    if clipped_end <= clipped_start {
        return 0.0;
    }
    (clipped_end - clipped_start) as f64 / (region_end - region_start) as f64
}

fn syng_candidate_features(
    syng_index: &syng::SyngIndex,
    path_name: &str,
    start: u64,
    end: u64,
) -> io::Result<CandidateFeatures> {
    let path_idx = syng_index
        .name_map
        .name_to_path
        .get(path_name)
        .copied()
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "candidate path '{}' is missing from syng name map",
                    path_name
                ),
            )
        })? as usize;
    let mut counts: FxHashMap<u32, u64> = FxHashMap::default();
    let mut oriented_walk = Vec::new();
    for (signed_node, _) in syng_index.walk_path_range(path_idx, start, end)? {
        *counts.entry(signed_node.unsigned_abs()).or_insert(0) += 1;
        oriented_walk.push(signed_node);
    }
    let mut unoriented: Vec<(u32, u64)> = counts.into_iter().collect();
    unoriented.sort_unstable_by_key(|&(feature_id, _)| feature_id);
    Ok(CandidateFeatures {
        unoriented,
        oriented_walk,
    })
}

#[allow(clippy::too_many_arguments)]
fn collect_syng_candidates(
    syng_index: &syng::SyngIndex,
    target_name: &str,
    region_start: u64,
    region_end: u64,
    padding: u64,
    extension: u64,
    mode: CandidateMode,
    min_anchors: usize,
    min_span_fraction: f64,
) -> io::Result<Vec<HaplotypeCandidate>> {
    let syncmer_len = syng_index.syncmer_length_bp() as u64;
    let hits = syng_index.query_region_with_anchors_ext(
        target_name,
        region_start,
        region_end,
        padding,
        extension,
    )?;
    let mut candidates = Vec::new();

    match mode {
        CandidateMode::Overlapping => {
            for hit in hits {
                if hit.anchors.len() < min_anchors || hit.end <= hit.start {
                    continue;
                }
                let query_span_fraction =
                    anchor_span_fraction(&hit.anchors, region_start, region_end, syncmer_len);
                let features =
                    syng_candidate_features(syng_index, &hit.genome, hit.start, hit.end)?;
                if features.unoriented.is_empty() {
                    continue;
                }
                candidates.push(HaplotypeCandidate {
                    path_name: hit.genome,
                    start: hit.start,
                    end: hit.end,
                    strand: hit.strand,
                    anchors: hit.anchors.len(),
                    query_span_fraction,
                    features: features.unoriented,
                    oriented_walk: features.oriented_walk,
                    single_similarity: 0.0,
                });
            }
        }
        CandidateMode::Spanning => {
            let mut groups: BTreeMap<(String, char), SyngCandidateGroup> = BTreeMap::new();
            for hit in hits {
                if hit.end <= hit.start {
                    continue;
                }
                let entry =
                    groups
                        .entry((hit.genome, hit.strand))
                        .or_insert_with(|| SyngCandidateGroup {
                            start: u64::MAX,
                            end: 0,
                            anchors: Vec::new(),
                        });
                entry.start = entry.start.min(hit.start);
                entry.end = entry.end.max(hit.end);
                entry.anchors.extend(hit.anchors);
            }

            for ((path_name, strand), mut group) in groups {
                group.anchors.sort_by(|a, b| {
                    a.query_pos
                        .cmp(&b.query_pos)
                        .then(a.target_pos.cmp(&b.target_pos))
                });
                group.anchors.dedup_by(|a, b| {
                    a.query_pos == b.query_pos
                        && a.target_pos == b.target_pos
                        && a.node_id == b.node_id
                });
                let query_span_fraction =
                    anchor_span_fraction(&group.anchors, region_start, region_end, syncmer_len);
                if group.anchors.len() < min_anchors || query_span_fraction < min_span_fraction {
                    continue;
                }
                let features =
                    syng_candidate_features(syng_index, &path_name, group.start, group.end)?;
                if features.unoriented.is_empty() {
                    continue;
                }
                candidates.push(HaplotypeCandidate {
                    path_name,
                    start: group.start,
                    end: group.end,
                    strand,
                    anchors: group.anchors.len(),
                    query_span_fraction,
                    features: features.unoriented,
                    oriented_walk: features.oriented_walk,
                    single_similarity: 0.0,
                });
            }
        }
    }

    candidates.sort_by(|a, b| {
        a.path_name
            .cmp(&b.path_name)
            .then(a.strand.cmp(&b.strand))
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });
    Ok(candidates)
}

fn locus_features(candidates: &[HaplotypeCandidate]) -> Vec<u32> {
    let mut seen: FxHashMap<u32, ()> = FxHashMap::default();
    for candidate in candidates {
        for &(feature_id, _) in &candidate.features {
            seen.insert(feature_id, ());
        }
    }
    let mut features: Vec<u32> = seen.into_keys().collect();
    features.sort_unstable();
    features
}

fn sample_norm_sq_for_features(sample_counts: &FxHashMap<u32, u64>, features: &[u32]) -> f64 {
    features
        .iter()
        .map(|feature_id| sample_counts.get(feature_id).copied().unwrap_or(0) as f64)
        .map(|v| v * v)
        .sum()
}

fn cosine_for_candidate_features(
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

fn score_combination(
    combination: &[usize],
    candidates: &[HaplotypeCandidate],
    sample_counts: &FxHashMap<u32, u64>,
    sample_norm_sq: f64,
) -> CosigtResult {
    let mut genotype_counts: FxHashMap<u32, u64> = FxHashMap::default();
    for &idx in combination {
        for &(feature_id, count) in &candidates[idx].features {
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

    CosigtResult {
        combination: combination.to_vec(),
        similarity,
        qv,
        dot,
        sample_norm,
        genotype_norm,
    }
}

struct CosigtSearch<'a> {
    candidates: &'a [HaplotypeCandidate],
    sample_counts: &'a FxHashMap<u32, u64>,
    sample_norm_sq: f64,
    ploidy: usize,
    max_combinations: u64,
    visited: u64,
    results: Vec<CosigtResult>,
}

impl CosigtSearch<'_> {
    fn run(mut self) -> io::Result<Vec<CosigtResult>> {
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
            self.results.push(score_combination(
                current,
                self.candidates,
                self.sample_counts,
                self.sample_norm_sq,
            ));
            return Ok(());
        }

        for idx in start..self.candidates.len() {
            current.push(idx);
            self.visit(idx, current)?;
            current.pop();
        }
        Ok(())
    }
}

fn rank_cosigt(
    candidates: &mut Vec<HaplotypeCandidate>,
    sample_counts: &FxHashMap<u32, u64>,
    ploidy: usize,
    top_n: usize,
    candidate_top_k: usize,
    max_combinations: u64,
) -> io::Result<(Vec<CosigtResult>, Vec<u32>)> {
    let all_locus_features = locus_features(candidates);
    let all_sample_norm_sq = sample_norm_sq_for_features(sample_counts, &all_locus_features);
    for candidate in candidates.iter_mut() {
        candidate.single_similarity =
            cosine_for_candidate_features(&candidate.features, sample_counts, all_sample_norm_sq);
    }
    candidates.sort_by(|a, b| {
        b.single_similarity
            .total_cmp(&a.single_similarity)
            .then(b.anchors.cmp(&a.anchors))
            .then(a.path_name.cmp(&b.path_name))
            .then(a.start.cmp(&b.start))
    });
    if candidate_top_k > 0 && candidates.len() > candidate_top_k {
        candidates.truncate(candidate_top_k);
    }

    let selected_features = locus_features(candidates);
    let sample_norm_sq = sample_norm_sq_for_features(sample_counts, &selected_features);
    if sample_norm_sq == 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sample pack has zero coverage over candidate graph features",
        ));
    }

    let search = CosigtSearch {
        candidates,
        sample_counts,
        sample_norm_sq,
        ploidy,
        max_combinations,
        visited: 0,
        results: Vec::new(),
    };
    let mut results = search.run()?;
    results.truncate(top_n);
    Ok((results, selected_features))
}

fn format_candidate_region(candidate: &HaplotypeCandidate) -> String {
    format!(
        "{}:{}-{}({})",
        candidate.path_name, candidate.start, candidate.end, candidate.strand
    )
}

pub fn run_syng_cosigt<W: Write>(out: &mut W, config: &SyngCosigtConfig<'_>) -> io::Result<()> {
    info!("Loading syng index from prefix: {}", config.syng_prefix);
    let syng_index = syng::SyngIndex::load(config.syng_prefix, syng::SyncmerParams::default())?;
    info!("Loading sample pack coverage from {}", config.pack_path);
    let coverage = pack::read(config.pack_path)?;
    let query = SyngCosigtQuery {
        target_range: config.target_range,
        candidate_mode: config.candidate_mode,
        ploidy: config.ploidy,
        top_n: config.top_n,
        candidate_top_k: config.candidate_top_k,
        max_combinations: config.max_combinations,
        syng_padding: config.syng_padding,
        syng_extension: config.syng_extension,
        min_anchors: config.min_anchors,
        min_span_fraction: config.min_span_fraction,
    };
    let result = compute_syng_cosigt(&syng_index, &coverage, &query)?;
    write_syng_cosigt_output(out, &result)
}

pub fn compute_syng_cosigt(
    syng_index: &syng::SyngIndex,
    coverage: &pack::Coverage,
    query: &SyngCosigtQuery<'_>,
) -> io::Result<SyngCosigtOutput> {
    if query.ploidy == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--ploidy must be greater than 0",
        ));
    }
    if query.top_n == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--top-n must be greater than 0",
        ));
    }
    if !(0.0..=1.0).contains(&query.min_span_fraction) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--min-span-fraction must be between 0 and 1",
        ));
    }

    let (target_name, (range_start, range_end), region_name) =
        partition::parse_target_range(query.target_range)?;
    if range_start < 0 || range_end <= range_start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "target range must be non-empty with non-negative coordinates: {}",
                query.target_range
            ),
        ));
    }
    let region_start = range_start as u64;
    let region_end = range_end as u64;

    let mut candidates = collect_syng_candidates(
        syng_index,
        &target_name,
        region_start,
        region_end,
        query.syng_padding,
        query.syng_extension,
        query.candidate_mode,
        query.min_anchors,
        query.min_span_fraction,
    )?;
    if candidates.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no syng genotype candidates were found for the requested range",
        ));
    }

    let (results, selected_features) = rank_cosigt(
        &mut candidates,
        &coverage.counts,
        query.ploidy,
        query.top_n,
        query.candidate_top_k,
        query.max_combinations,
    )?;

    Ok(SyngCosigtOutput {
        region_name,
        candidate_mode: query.candidate_mode,
        ploidy: query.ploidy,
        candidates,
        selected_features,
        pack_nonzero_nodes: coverage.nonzero_nodes,
        pack_universe_nodes: coverage.universe_nodes,
        pack_retained_records: coverage.retained_records,
        pack_syncmer_anchors: coverage.syncmer_anchors,
        results,
    })
}

pub fn write_syng_cosigt_output<W: Write>(
    out: &mut W,
    result: &SyngCosigtOutput,
) -> io::Result<()> {
    writeln!(out, "#impg genotype cos")?;
    writeln!(out, "#region\t{}", result.region_name)?;
    writeln!(out, "#method\tcos")?;
    writeln!(out, "#metric\tcosine")?;
    writeln!(out, "#alias\tcosigt")?;
    writeln!(out, "#feature_space\tsyng-syncmer-node")?;
    writeln!(out, "#candidate_mode\t{:?}", result.candidate_mode)?;
    writeln!(out, "#ploidy\t{}", result.ploidy)?;
    writeln!(out, "#candidates\t{}", result.candidates.len())?;
    writeln!(out, "#locus_features\t{}", result.selected_features.len())?;
    writeln!(out, "#pack_nonzero_nodes\t{}", result.pack_nonzero_nodes)?;
    if let Some(universe_nodes) = result.pack_universe_nodes {
        writeln!(out, "#pack_universe_nodes\t{}", universe_nodes)?;
    }
    if let Some(retained_records) = result.pack_retained_records {
        writeln!(out, "#pack_retained_records\t{}", retained_records)?;
    }
    if let Some(syncmer_anchors) = result.pack_syncmer_anchors {
        writeln!(out, "#pack_syncmer_anchors\t{}", syncmer_anchors)?;
    }
    writeln!(
        out,
        "#rank\tmethod\tploidy\tsimilarity\tqv\tdot\tsample_norm\tgenotype_norm\thaplotypes\tregions\tcandidate_anchors\tcandidate_span_fractions"
    )?;
    for (rank, genotype_result) in result.results.iter().enumerate() {
        let haplotypes: Vec<&str> = genotype_result
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].path_name.as_str())
            .collect();
        let regions: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|&idx| format_candidate_region(&result.candidates[idx]))
            .collect();
        let anchors: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].anchors.to_string())
            .collect();
        let spans: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|&idx| format!("{:.6}", result.candidates[idx].query_span_fraction))
            .collect();
        writeln!(
            out,
            "{}\tcos\t{}\t{:.9}\t{:.3}\t{:.3}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            rank + 1,
            result.ploidy,
            genotype_result.similarity,
            genotype_result.qv,
            genotype_result.dot,
            genotype_result.sample_norm,
            genotype_result.genotype_norm,
            haplotypes.join(","),
            regions.join(","),
            anchors.join(","),
            spans.join(","),
        )?;
    }

    Ok(())
}
