use clap::ValueEnum;
use log::info;
use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, BTreeSet};
use std::fs::File;
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

use crate::commands::partition;
use crate::genotyping::{self, FeatureSpace, ScoringMethod};
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
    pub oriented_walk: Vec<syng::SyngWalkStep>,
    single_similarity: f64,
}

struct CandidateFeatures {
    unoriented: Vec<(u32, u64)>,
    oriented_walk: Vec<syng::SyngWalkStep>,
}

#[derive(Debug)]
struct SyngCandidateGroup {
    start: u64,
    end: u64,
    anchors: Vec<syng::Anchor>,
}

pub type CosigtResult = genotyping::CombinationScore;

pub struct SyngCosigtConfig<'a> {
    pub syng_prefix: &'a str,
    pub pack_path: &'a str,
    pub target_range: &'a str,
    pub emit_report_path: Option<&'a str>,
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

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum GraphFeatureIdMode {
    /// Use numeric GFA segment names when all segment names are positive u32s; otherwise use dense import order.
    Auto,
    /// Assign feature IDs 1..N in S-line import order.
    Dense,
    /// Require every S-line segment name to be a positive u32 and use that as the feature ID.
    SegmentName,
}

impl GraphFeatureIdMode {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Auto => "auto",
            Self::Dense => "dense",
            Self::SegmentName => "segment-name",
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum GraphContributionModel {
    /// Score raw graph-node traversal counts.
    Raw,
    /// Divide sample node counts by segment length and count candidate partial traversals by covered-bp / segment-length.
    LengthNormalized,
}

impl GraphContributionModel {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Raw => "raw",
            Self::LengthNormalized => "length-normalized",
        }
    }
}

pub enum GraphSource<'a> {
    Gfa {
        path: &'a str,
    },
    RenderBundle {
        path: &'a str,
    },
    Query {
        gfa: &'a str,
        description: &'a str,
        build_command: Option<&'a str>,
        default_target_path: Option<&'a str>,
    },
}

pub struct GraphCosigtConfig<'a> {
    pub graph_source: GraphSource<'a>,
    pub pack_path: &'a str,
    pub target_path: Option<&'a str>,
    pub pack_feature_space: Option<&'a str>,
    pub pack_graph_id: Option<&'a str>,
    pub feature_id_mode: GraphFeatureIdMode,
    pub contribution_model: GraphContributionModel,
    pub emit_report_path: Option<&'a str>,
    pub ploidy: usize,
    pub top_n: usize,
    pub candidate_top_k: usize,
    pub max_combinations: u64,
}

#[derive(Debug, Clone)]
pub struct GraphSegment {
    pub name: String,
    pub feature_id: u32,
    pub length: u64,
    sequence: Option<String>,
}

#[derive(Debug, Clone)]
pub struct GraphPathStep {
    pub segment_index: usize,
    pub feature_id: u32,
    pub orientation: char,
    pub start: u64,
    pub end: u64,
}

#[derive(Debug, Clone)]
pub struct GraphPath {
    pub name: String,
    pub source_record: char,
    pub steps: Vec<GraphPathStep>,
    pub length: u64,
}

#[derive(Debug, Clone)]
pub struct NormalizedGraph {
    pub source_label: String,
    pub source_path: Option<PathBuf>,
    pub build_command: Option<String>,
    pub feature_space: String,
    pub graph_id: String,
    pub requested_feature_id_mode: GraphFeatureIdMode,
    pub effective_feature_id_mode: GraphFeatureIdMode,
    pub segments: Vec<GraphSegment>,
    pub paths: Vec<GraphPath>,
    feature_ids: BTreeSet<u32>,
    feature_to_segment: FxHashMap<u32, usize>,
    path_name_to_index: FxHashMap<String, usize>,
}

#[derive(Debug, Clone)]
pub struct GraphHaplotypeCandidate {
    pub path_name: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
    pub anchors: usize,
    pub query_span_fraction: f64,
    pub features: Vec<(u32, u64)>,
    pub scoring_features: Vec<(u32, f64)>,
    single_similarity: f64,
}

#[derive(Debug, Clone, Default)]
pub struct GraphPackMetadata {
    pub feature_space: Option<String>,
    pub graph_id: Option<String>,
    pub feature_id_mode: Option<String>,
    pub contribution_model: Option<String>,
    pub source: Option<String>,
}

#[derive(Debug)]
pub struct GraphCosigtOutput {
    pub region_name: String,
    pub target_path: Option<String>,
    pub ploidy: usize,
    pub contribution_model: GraphContributionModel,
    pub graph: NormalizedGraph,
    pub candidates: Vec<GraphHaplotypeCandidate>,
    pub selected_features: Vec<u32>,
    pub sample_raw_counts: FxHashMap<u32, u64>,
    pub sample_weights: FxHashMap<u32, f64>,
    pub pack_metadata: GraphPackMetadata,
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
    for (signed_node, bp_pos) in syng_index.walk_path_range(path_idx, start, end)? {
        *counts.entry(signed_node.unsigned_abs()).or_insert(0) += 1;
        oriented_walk.push(syng::SyngWalkStep {
            signed_node,
            bp_pos,
        });
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

fn rank_cosigt(
    candidates: &mut Vec<HaplotypeCandidate>,
    sample_counts: &FxHashMap<u32, u64>,
    ploidy: usize,
    top_n: usize,
    candidate_top_k: usize,
    max_combinations: u64,
) -> io::Result<(Vec<CosigtResult>, Vec<u32>)> {
    let all_locus_features = genotyping::feature_universe(
        candidates
            .iter()
            .map(|candidate| candidate.features.as_slice()),
    );
    let all_sample_norm_sq =
        genotyping::sample_norm_sq_for_features(sample_counts, &all_locus_features);
    for candidate in candidates.iter_mut() {
        candidate.single_similarity = genotyping::cosine_for_feature_counts(
            &candidate.features,
            sample_counts,
            all_sample_norm_sq,
        );
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

    let selected_features = genotyping::feature_universe(
        candidates
            .iter()
            .map(|candidate| candidate.features.as_slice()),
    );
    let sample_norm_sq = genotyping::sample_norm_sq_for_features(sample_counts, &selected_features);
    if sample_norm_sq == 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sample pack has zero coverage over candidate graph features",
        ));
    }

    let candidate_features: Vec<&[(u32, u64)]> = candidates
        .iter()
        .map(|candidate| candidate.features.as_slice())
        .collect();
    let mut results = genotyping::run_cosine_combination_search(
        &candidate_features,
        sample_counts,
        sample_norm_sq,
        ploidy,
        max_combinations,
    )?;
    results.truncate(top_n);
    Ok((results, selected_features))
}

fn format_candidate_region(candidate: &HaplotypeCandidate) -> String {
    format!(
        "{}:{}-{}({})",
        candidate.path_name, candidate.start, candidate.end, candidate.strand
    )
}

fn candidate_count_for_feature(candidate: &HaplotypeCandidate, node_id: u32) -> u64 {
    candidate
        .features
        .binary_search_by_key(&node_id, |&(feature_id, _)| feature_id)
        .map(|idx| candidate.features[idx].1)
        .unwrap_or(0)
}

fn candidate_count_mass(candidate: &HaplotypeCandidate) -> u64 {
    candidate.features.iter().map(|&(_, count)| count).sum()
}

fn candidate_repeated_stats(candidate: &HaplotypeCandidate) -> (usize, u64, u64) {
    let repeated_nodes = candidate
        .features
        .iter()
        .filter(|&&(_, count)| count > 1)
        .count();
    let repeated_extra_count = candidate
        .features
        .iter()
        .map(|&(_, count)| count.saturating_sub(1))
        .sum();
    let max_node_count = candidate
        .features
        .iter()
        .map(|&(_, count)| count)
        .max()
        .unwrap_or(0);
    (repeated_nodes, repeated_extra_count, max_node_count)
}

fn candidate_sample_overlap(
    candidate: &HaplotypeCandidate,
    sample_counts: &FxHashMap<u32, u64>,
) -> (usize, u64, f64) {
    let mut overlap_unique_nodes = 0usize;
    let mut overlap_sample_count_mass = 0u64;
    let mut dot = 0.0;
    for &(node_id, candidate_count) in &candidate.features {
        let sample_count = sample_counts.get(&node_id).copied().unwrap_or(0);
        if sample_count > 0 {
            overlap_unique_nodes += 1;
            overlap_sample_count_mass += sample_count;
            dot += sample_count as f64 * candidate_count as f64;
        }
    }
    (overlap_unique_nodes, overlap_sample_count_mass, dot)
}

fn format_optional_u64(value: Option<u64>) -> String {
    value
        .map(|value| value.to_string())
        .unwrap_or_else(|| "NA".to_string())
}

fn write_report_kv<W: Write>(
    out: &mut W,
    key: &str,
    value: impl std::fmt::Display,
) -> io::Result<()> {
    writeln!(out, "{key}\t{value}")
}

pub fn write_syng_cosigt_report<W: Write>(
    out: &mut W,
    config: &SyngCosigtConfig<'_>,
    coverage: &pack::Coverage,
    result: &SyngCosigtOutput,
) -> io::Result<()> {
    writeln!(out, "#impg genotype cos report")?;
    writeln!(out, "#format\tsectioned-tsv-v1")?;

    writeln!(out, "#section\tinput_metadata")?;
    writeln!(out, "key\tvalue")?;
    write_report_kv(out, "syng_prefix", config.syng_prefix)?;
    write_report_kv(out, "pack_path", config.pack_path)?;
    write_report_kv(out, "target_range", config.target_range)?;
    write_report_kv(out, "resolved_region", &result.region_name)?;
    write_report_kv(out, "feature_space", FeatureSpace::SyngSyncmerNode.as_str())?;
    write_report_kv(out, "method", ScoringMethod::Cos.as_str())?;
    write_report_kv(out, "metric", ScoringMethod::Cos.metric_str())?;
    write_report_kv(
        out,
        "candidate_mode",
        format!("{:?}", config.candidate_mode),
    )?;
    write_report_kv(out, "ploidy", config.ploidy)?;
    write_report_kv(out, "top_n", config.top_n)?;
    write_report_kv(out, "candidate_top_k", config.candidate_top_k)?;
    write_report_kv(out, "max_combinations", config.max_combinations)?;
    write_report_kv(out, "syng_padding", config.syng_padding)?;
    write_report_kv(out, "syng_extension", config.syng_extension)?;
    write_report_kv(out, "min_anchors", config.min_anchors)?;
    write_report_kv(
        out,
        "min_span_fraction",
        format!("{:.6}", config.min_span_fraction),
    )?;
    write_report_kv(out, "pack_nonzero_nodes", coverage.nonzero_nodes)?;
    write_report_kv(
        out,
        "pack_universe_nodes",
        format_optional_u64(coverage.universe_nodes),
    )?;
    write_report_kv(
        out,
        "pack_retained_records",
        format_optional_u64(coverage.retained_records),
    )?;
    write_report_kv(
        out,
        "pack_syncmer_anchors",
        format_optional_u64(coverage.syncmer_anchors),
    )?;
    write_report_kv(
        out,
        "sample_pack_counting_semantics",
        "distinct_nodes_per_read",
    )?;
    write_report_kv(
        out,
        "sample_pack_counting_detail",
        "impg map pack counts each distinct syng node at most once per retained read; repeated node occurrences in one read do not increase that node count",
    )?;

    let sample_locus_count_mass: u64 = result
        .selected_features
        .iter()
        .map(|node_id| coverage.counts.get(node_id).copied().unwrap_or(0))
        .sum();
    let sample_locus_overlap_features = result
        .selected_features
        .iter()
        .filter(|node_id| coverage.counts.get(node_id).copied().unwrap_or(0) > 0)
        .count();
    let sample_norm_sq =
        genotyping::sample_norm_sq_for_features(&coverage.counts, &result.selected_features);
    writeln!(out, "#section\tpack_evidence_summary")?;
    writeln!(out, "metric\tvalue")?;
    write_report_kv(out, "pack_nonzero_nodes", coverage.nonzero_nodes)?;
    write_report_kv(
        out,
        "pack_universe_nodes",
        format_optional_u64(coverage.universe_nodes),
    )?;
    write_report_kv(
        out,
        "pack_retained_records",
        format_optional_u64(coverage.retained_records),
    )?;
    write_report_kv(
        out,
        "pack_syncmer_anchors",
        format_optional_u64(coverage.syncmer_anchors),
    )?;
    write_report_kv(
        out,
        "selected_locus_features",
        result.selected_features.len(),
    )?;
    write_report_kv(
        out,
        "locus_feature_overlap_nonzero_nodes",
        sample_locus_overlap_features,
    )?;
    write_report_kv(
        out,
        "locus_feature_overlap_sample_count_mass",
        sample_locus_count_mass,
    )?;
    write_report_kv(
        out,
        "sample_norm_over_selected_locus_features",
        format!("{:.6}", sample_norm_sq.sqrt()),
    )?;

    writeln!(out, "#section\tsample_locus_features")?;
    writeln!(out, "node_id\tsample_count")?;
    for &node_id in &result.selected_features {
        writeln!(
            out,
            "{}\t{}",
            node_id,
            coverage.counts.get(&node_id).copied().unwrap_or(0)
        )?;
    }

    writeln!(out, "#section\tcandidates")?;
    writeln!(
        out,
        "candidate_index\tpath\tinterval\tstart\tend\tstrand\tanchors\tspan_fraction\tfeature_count\ttotal_candidate_node_count_mass\tunique_nodes\trepeated_nodes\trepeated_extra_count\tmax_node_count\tsingle_haplotype_cosine\tsample_overlap_unique_nodes\tsample_overlap_sample_count_mass\tsample_overlap_dot_contribution"
    )?;
    for (idx, candidate) in result.candidates.iter().enumerate() {
        let total_mass = candidate_count_mass(candidate);
        let (repeated_nodes, repeated_extra_count, max_node_count) =
            candidate_repeated_stats(candidate);
        let (overlap_unique_nodes, overlap_sample_count_mass, overlap_dot) =
            candidate_sample_overlap(candidate, &coverage.counts);
        writeln!(
            out,
            "{}\t{}\t{}:{}-{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{}\t{}\t{:.3}",
            idx,
            candidate.path_name,
            candidate.path_name,
            candidate.start,
            candidate.end,
            candidate.start,
            candidate.end,
            candidate.strand,
            candidate.anchors,
            candidate.query_span_fraction,
            candidate.features.len(),
            total_mass,
            candidate.features.len(),
            repeated_nodes,
            repeated_extra_count,
            max_node_count,
            candidate.single_similarity,
            overlap_unique_nodes,
            overlap_sample_count_mass,
            overlap_dot,
        )?;
    }

    writeln!(out, "#section\tcandidate_features")?;
    writeln!(
        out,
        "candidate_index\tpath\tnode_id\tsample_count\tcandidate_count\tdot_contribution"
    )?;
    for (idx, candidate) in result.candidates.iter().enumerate() {
        for &(node_id, candidate_count) in &candidate.features {
            let sample_count = coverage.counts.get(&node_id).copied().unwrap_or(0);
            let dot_contribution = sample_count as f64 * candidate_count as f64;
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{:.3}",
                idx, candidate.path_name, node_id, sample_count, candidate_count, dot_contribution,
            )?;
        }
    }

    writeln!(out, "#section\tresult_scores")?;
    writeln!(
        out,
        "rank\tmethod\tploidy\tsimilarity\tqv\tdot\tsample_norm\tgenotype_norm\tcandidate_indices\thaplotypes\tregions\tcandidate_anchors\tcandidate_span_fractions"
    )?;
    for (rank, genotype_result) in result.results.iter().enumerate() {
        let candidate_indices: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|idx| idx.to_string())
            .collect();
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
            "{}\t{}\t{}\t{:.9}\t{:.3}\t{:.3}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            rank + 1,
            ScoringMethod::Cos.as_str(),
            result.ploidy,
            genotype_result.similarity,
            genotype_result.qv,
            genotype_result.dot,
            genotype_result.sample_norm,
            genotype_result.genotype_norm,
            candidate_indices.join(","),
            haplotypes.join(","),
            regions.join(","),
            anchors.join(","),
            spans.join(","),
        )?;
    }

    writeln!(out, "#section\tresult_features")?;
    writeln!(
        out,
        "rank\tnode_id\tsample_count\tcandidate_counts\tgenotype_count\tdot_contribution"
    )?;
    for (rank, genotype_result) in result.results.iter().enumerate() {
        let mut genotype_counts = BTreeMap::new();
        for &candidate_idx in &genotype_result.combination {
            for &(node_id, count) in &result.candidates[candidate_idx].features {
                *genotype_counts.entry(node_id).or_insert(0u64) += count;
            }
        }
        for (node_id, genotype_count) in genotype_counts {
            let sample_count = coverage.counts.get(&node_id).copied().unwrap_or(0);
            let candidate_counts: Vec<String> = genotype_result
                .combination
                .iter()
                .map(|&candidate_idx| {
                    candidate_count_for_feature(&result.candidates[candidate_idx], node_id)
                        .to_string()
                })
                .collect();
            let dot_contribution = sample_count as f64 * genotype_count as f64;
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{:.3}",
                rank + 1,
                node_id,
                sample_count,
                candidate_counts.join(","),
                genotype_count,
                dot_contribution,
            )?;
        }
    }

    Ok(())
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
    write_syng_cosigt_output(out, &result)?;
    if let Some(path) = config.emit_report_path {
        let mut report = BufWriter::new(File::create(path)?);
        write_syng_cosigt_report(&mut report, config, &coverage, &result)?;
        report.flush()?;
    }
    Ok(())
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
    writeln!(out, "#method\t{}", ScoringMethod::Cos.as_str())?;
    writeln!(out, "#metric\t{}", ScoringMethod::Cos.metric_str())?;
    writeln!(out, "#alias\tcosigt")?;
    writeln!(
        out,
        "#feature_space\t{}",
        FeatureSpace::SyngSyncmerNode.as_str()
    )?;
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
            "{}\t{}\t{}\t{:.9}\t{:.3}\t{:.3}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            rank + 1,
            ScoringMethod::Cos.as_str(),
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

struct LoadedGraphSource {
    source_label: String,
    source_path: Option<PathBuf>,
    build_command: Option<String>,
    default_target_path: Option<String>,
    feature_space: String,
    gfa_text: String,
}

fn load_graph_source(source: &GraphSource<'_>) -> io::Result<LoadedGraphSource> {
    match source {
        GraphSource::Gfa { path } => Ok(LoadedGraphSource {
            source_label: "gfa".to_string(),
            source_path: Some(PathBuf::from(path)),
            build_command: None,
            default_target_path: None,
            feature_space: FeatureSpace::GfaSegment.as_str().to_string(),
            gfa_text: std::fs::read_to_string(path)?,
        }),
        GraphSource::RenderBundle { path } => {
            let bundle = crate::render_bundle::load_bundle(path)?;
            let feature_space = bundle.manifest.feature_space.clone();
            if !matches!(
                feature_space.as_str(),
                "gfa-segment" | "variation-graph-node"
            ) {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "render bundle feature_space '{}' is not a graph-node genotype space",
                        feature_space
                    ),
                ));
            }
            let graph_gfa = bundle.manifest.graph_gfa.as_deref().ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "render bundle does not contain graph_gfa",
                )
            })?;
            let graph_path = crate::render_bundle::resolve_bundle_path(&bundle.root, graph_gfa);
            let default_target_path = bundle.rendered_target_range(None).ok();
            Ok(LoadedGraphSource {
                source_label: format!("render-bundle:{}", bundle.manifest.engine),
                source_path: Some(PathBuf::from(path)),
                build_command: None,
                default_target_path,
                feature_space,
                gfa_text: std::fs::read_to_string(graph_path)?,
            })
        }
        GraphSource::Query {
            gfa,
            description,
            build_command,
            default_target_path,
        } => Ok(LoadedGraphSource {
            source_label: format!("query:{description}"),
            source_path: None,
            build_command: build_command.map(ToString::to_string),
            default_target_path: default_target_path.map(ToString::to_string),
            feature_space: FeatureSpace::GfaSegment.as_str().to_string(),
            gfa_text: (*gfa).to_string(),
        }),
    }
}

fn parse_segment_length(fields: &[&str], line_no: usize) -> io::Result<(u64, Option<String>)> {
    let seq = fields.get(2).copied().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("GFA S-line {} is missing sequence", line_no),
        )
    })?;
    if seq != "*" {
        return Ok((seq.len() as u64, Some(seq.to_string())));
    }
    for field in fields.iter().skip(3) {
        if let Some(value) = field.strip_prefix("LN:i:") {
            let len = value.parse::<u64>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GFA S-line {} has invalid LN tag '{}'", line_no, field),
                )
            })?;
            return Ok((len, None));
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "GFA S-line {} has '*' sequence but no LN:i:length tag",
            line_no
        ),
    ))
}

fn parse_p_step_token(token: &str, line_no: usize) -> io::Result<(String, char)> {
    if let Some(name) = token.strip_suffix('+') {
        if !name.is_empty() {
            return Ok((name.to_string(), '+'));
        }
    }
    if let Some(name) = token.strip_suffix('-') {
        if !name.is_empty() {
            return Ok((name.to_string(), '-'));
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!("GFA P-line {} has invalid step token '{}'", line_no, token),
    ))
}

fn parse_w_walk(walk: &str, line_no: usize) -> io::Result<Vec<(String, char)>> {
    let mut parsed = Vec::new();
    let mut orientation: Option<char> = None;
    let mut name = String::new();
    for ch in walk.chars() {
        if ch == '>' || ch == '<' {
            if let Some(prev) = orientation {
                if name.is_empty() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA W-line {} has an empty walk step", line_no),
                    ));
                }
                parsed.push((std::mem::take(&mut name), prev));
            }
            orientation = Some(if ch == '>' { '+' } else { '-' });
        } else {
            if orientation.is_none() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GFA W-line {} walk must start with '>' or '<'", line_no),
                ));
            }
            name.push(ch);
        }
    }
    if let Some(prev) = orientation {
        if name.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GFA W-line {} has an empty trailing walk step", line_no),
            ));
        }
        parsed.push((name, prev));
    }
    if parsed.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("GFA W-line {} has an empty walk", line_no),
        ));
    }
    Ok(parsed)
}

fn parse_overlap_bases(raw: &str) -> u64 {
    if raw == "*" {
        return 0;
    }
    raw.strip_suffix('M')
        .and_then(|value| value.parse::<u64>().ok())
        .unwrap_or(0)
}

fn build_graph_path(
    name: String,
    source_record: char,
    parsed_steps: Vec<(String, char)>,
    overlaps: &[&str],
    line_no: usize,
    segment_name_to_index: &FxHashMap<String, usize>,
    segments: &[GraphSegment],
) -> io::Result<GraphPath> {
    let mut steps = Vec::with_capacity(parsed_steps.len());
    let mut pos = 0u64;
    for (idx, (segment_name, orientation)) in parsed_steps.into_iter().enumerate() {
        let segment_index = segment_name_to_index
            .get(&segment_name)
            .copied()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GFA {}-line {} references unknown segment '{}'",
                        source_record, line_no, segment_name
                    ),
                )
            })?;
        let segment = &segments[segment_index];
        let start = pos;
        let end = start.checked_add(segment.length).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GFA path '{}' length overflows u64", name),
            )
        })?;
        steps.push(GraphPathStep {
            segment_index,
            feature_id: segment.feature_id,
            orientation,
            start,
            end,
        });
        pos = end;
        if idx < overlaps.len() {
            let overlap = parse_overlap_bases(overlaps[idx]);
            pos = pos.saturating_sub(overlap.min(segment.length));
        }
    }
    Ok(GraphPath {
        name,
        source_record,
        steps,
        length: pos,
    })
}

fn resolve_feature_id_mode(
    segments: &[GraphSegment],
    requested: GraphFeatureIdMode,
) -> io::Result<GraphFeatureIdMode> {
    let segment_name_mode_valid = || {
        let mut seen = BTreeSet::new();
        segments.iter().all(|segment| {
            segment
                .name
                .parse::<u32>()
                .ok()
                .filter(|&id| id > 0)
                .is_some_and(|id| seen.insert(id))
        })
    };
    match requested {
        GraphFeatureIdMode::Auto => {
            if segment_name_mode_valid() {
                Ok(GraphFeatureIdMode::SegmentName)
            } else {
                Ok(GraphFeatureIdMode::Dense)
            }
        }
        GraphFeatureIdMode::Dense => Ok(GraphFeatureIdMode::Dense),
        GraphFeatureIdMode::SegmentName => {
            if segment_name_mode_valid() {
                Ok(GraphFeatureIdMode::SegmentName)
            } else {
                Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    "--graph-feature-id-mode segment-name requires every GFA S-line name to be a unique positive u32",
                ))
            }
        }
    }
}

fn fnv1a_update(hash: &mut u64, bytes: &[u8]) {
    const FNV_PRIME: u64 = 0x0000_0100_0000_01b3;
    for byte in bytes {
        *hash ^= u64::from(*byte);
        *hash = hash.wrapping_mul(FNV_PRIME);
    }
    *hash ^= 0xff;
    *hash = hash.wrapping_mul(FNV_PRIME);
}

fn graph_id_for(graph: &NormalizedGraph) -> String {
    let mut hash = 0xcbf2_9ce4_8422_2325u64;
    fnv1a_update(&mut hash, graph.feature_space.as_bytes());
    fnv1a_update(
        &mut hash,
        graph.effective_feature_id_mode.as_str().as_bytes(),
    );
    for segment in &graph.segments {
        fnv1a_update(&mut hash, b"S");
        fnv1a_update(&mut hash, segment.name.as_bytes());
        fnv1a_update(&mut hash, segment.feature_id.to_string().as_bytes());
        fnv1a_update(&mut hash, segment.length.to_string().as_bytes());
        if let Some(sequence) = &segment.sequence {
            fnv1a_update(&mut hash, sequence.as_bytes());
        }
    }
    for path in &graph.paths {
        fnv1a_update(&mut hash, b"P");
        fnv1a_update(&mut hash, path.name.as_bytes());
        fnv1a_update(&mut hash, &[path.source_record as u8]);
        for step in &path.steps {
            fnv1a_update(
                &mut hash,
                graph.segments[step.segment_index].name.as_bytes(),
            );
            fnv1a_update(&mut hash, &[step.orientation as u8]);
        }
    }
    format!("gfa-fnv1a64:{hash:016x}")
}

pub fn parse_normalized_gfa(
    gfa_text: &str,
    source_label: impl Into<String>,
    source_path: Option<PathBuf>,
    build_command: Option<String>,
    feature_space: impl Into<String>,
    requested_feature_id_mode: GraphFeatureIdMode,
) -> io::Result<NormalizedGraph> {
    let feature_space = feature_space.into();
    let mut segments = Vec::new();
    let mut segment_name_to_index = FxHashMap::default();

    for (line_idx, line) in gfa_text.lines().enumerate() {
        let line_no = line_idx + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.first().copied() != Some("S") {
            continue;
        }
        let name = fields.get(1).copied().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GFA S-line {} is missing segment name", line_no),
            )
        })?;
        if segment_name_to_index.contains_key(name) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GFA repeats segment name '{}'", name),
            ));
        }
        let (length, sequence) = parse_segment_length(&fields, line_no)?;
        if length == 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GFA segment '{}' has zero length", name),
            ));
        }
        segment_name_to_index.insert(name.to_string(), segments.len());
        segments.push(GraphSegment {
            name: name.to_string(),
            feature_id: 0,
            length,
            sequence,
        });
    }
    if segments.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GFA contains no S-lines",
        ));
    }

    let effective_feature_id_mode = resolve_feature_id_mode(&segments, requested_feature_id_mode)?;
    let mut feature_ids = BTreeSet::new();
    let mut feature_to_segment = FxHashMap::default();
    for (idx, segment) in segments.iter_mut().enumerate() {
        segment.feature_id = match effective_feature_id_mode {
            GraphFeatureIdMode::Dense => (idx + 1) as u32,
            GraphFeatureIdMode::SegmentName => segment.name.parse::<u32>().map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("segment '{}' is not a positive u32", segment.name),
                )
            })?,
            GraphFeatureIdMode::Auto => unreachable!("auto was resolved"),
        };
        if !feature_ids.insert(segment.feature_id) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate graph feature ID {}", segment.feature_id),
            ));
        }
        feature_to_segment.insert(segment.feature_id, idx);
    }

    let mut paths = Vec::new();
    let mut path_name_to_index = FxHashMap::default();
    for (line_idx, line) in gfa_text.lines().enumerate() {
        let line_no = line_idx + 1;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("P") => {
                let name = fields.get(1).copied().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA P-line {} is missing path name", line_no),
                    )
                })?;
                let walk = fields.get(2).copied().ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA P-line {} is missing segment walk", line_no),
                    )
                })?;
                let parsed_steps = if walk == "*" || walk.is_empty() {
                    Vec::new()
                } else {
                    walk.split(',')
                        .map(|token| parse_p_step_token(token, line_no))
                        .collect::<io::Result<Vec<_>>>()?
                };
                if parsed_steps.is_empty() {
                    continue;
                }
                let overlaps: Vec<&str> = fields
                    .get(3)
                    .map(|value| value.split(',').collect())
                    .unwrap_or_default();
                let path = build_graph_path(
                    name.to_string(),
                    'P',
                    parsed_steps,
                    &overlaps,
                    line_no,
                    &segment_name_to_index,
                    &segments,
                )?;
                if path_name_to_index
                    .insert(path.name.clone(), paths.len())
                    .is_some()
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA repeats path/walk name '{}'", path.name),
                    ));
                }
                paths.push(path);
            }
            Some("W") => {
                if fields.len() < 7 {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA W-line {} has fewer than 7 fields", line_no),
                    ));
                }
                let sample = fields[1];
                let hap = fields[2];
                let seqid = fields[3];
                let name = if hap == "*" {
                    format!("{sample}#{seqid}")
                } else {
                    format!("{sample}#{hap}#{seqid}")
                };
                let parsed_steps = parse_w_walk(fields[6], line_no)?;
                let path = build_graph_path(
                    name,
                    'W',
                    parsed_steps,
                    &[],
                    line_no,
                    &segment_name_to_index,
                    &segments,
                )?;
                if path_name_to_index
                    .insert(path.name.clone(), paths.len())
                    .is_some()
                {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GFA repeats path/walk name '{}'", path.name),
                    ));
                }
                paths.push(path);
            }
            _ => {}
        }
    }
    if paths.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "GFA contains no P-lines or W-lines",
        ));
    }

    let mut graph = NormalizedGraph {
        source_label: source_label.into(),
        source_path,
        build_command,
        feature_space,
        graph_id: String::new(),
        requested_feature_id_mode,
        effective_feature_id_mode,
        segments,
        paths,
        feature_ids,
        feature_to_segment,
        path_name_to_index,
    };
    graph.graph_id = graph_id_for(&graph);
    Ok(graph)
}

pub fn parse_graph_path_spec(spec: &str) -> io::Result<(String, Option<(u64, u64)>)> {
    if let Some((name, coords)) = spec.rsplit_once(':') {
        if let Some((start, end)) = coords.split_once('-') {
            if let (Ok(start), Ok(end)) = (start.parse::<u64>(), end.parse::<u64>()) {
                if end <= start {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("target path interval '{}' is empty", spec),
                    ));
                }
                return Ok((name.to_string(), Some((start, end))));
            }
        }
    }
    Ok((spec.to_string(), None))
}

fn build_graph_candidate_for_path(
    graph: &NormalizedGraph,
    path: &GraphPath,
    start: u64,
    end: u64,
    target_span: u64,
    contribution_model: GraphContributionModel,
) -> Option<GraphHaplotypeCandidate> {
    if end <= start {
        return None;
    }
    let mut raw_counts: BTreeMap<u32, u64> = BTreeMap::new();
    let mut scoring_counts: BTreeMap<u32, f64> = BTreeMap::new();
    let mut anchors = 0usize;
    let mut covered_bp = 0u64;
    for step in &path.steps {
        let overlap_start = start.max(step.start);
        let overlap_end = end.min(step.end);
        if overlap_end <= overlap_start {
            continue;
        }
        anchors += 1;
        let overlap_bp = overlap_end - overlap_start;
        covered_bp += overlap_bp;
        *raw_counts.entry(step.feature_id).or_insert(0) += 1;
        let segment_len = graph.segments[step.segment_index].length.max(1) as f64;
        let weight = match contribution_model {
            GraphContributionModel::Raw => 1.0,
            GraphContributionModel::LengthNormalized => overlap_bp as f64 / segment_len,
        };
        *scoring_counts.entry(step.feature_id).or_insert(0.0) += weight;
    }
    if raw_counts.is_empty() {
        return None;
    }
    let query_span_fraction = if target_span == 0 {
        1.0
    } else {
        (covered_bp.min(target_span)) as f64 / target_span as f64
    };
    Some(GraphHaplotypeCandidate {
        path_name: path.name.clone(),
        start,
        end,
        strand: '+',
        anchors,
        query_span_fraction,
        features: raw_counts.into_iter().collect(),
        scoring_features: scoring_counts.into_iter().collect(),
        single_similarity: 0.0,
    })
}

pub fn collect_graph_candidates(
    graph: &NormalizedGraph,
    target_path: Option<&str>,
    contribution_model: GraphContributionModel,
) -> io::Result<(Vec<GraphHaplotypeCandidate>, String, Option<String>)> {
    let (target_name, target_interval) = if let Some(target_path) = target_path {
        let (name, interval) = parse_graph_path_spec(target_path)?;
        let target_idx = graph.path_name_to_index.get(&name).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("target path '{}' is not present in graph", name),
            )
        })?;
        let target_len = graph.paths[*target_idx].length;
        if let Some((start, end)) = interval {
            if end > target_len {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!(
                        "target path interval '{}:{}-{}' exceeds path length {}",
                        name, start, end, target_len
                    ),
                ));
            }
            (Some(name), Some((start, end)))
        } else {
            (Some(name), None)
        }
    } else {
        (None, None)
    };

    let target_span = target_interval.map(|(start, end)| end - start).unwrap_or(0);
    let mut candidates = Vec::new();
    for path in &graph.paths {
        let (start, end) = if let Some((start, end)) = target_interval {
            (start.min(path.length), end.min(path.length))
        } else {
            (0, path.length)
        };
        if let Some(candidate) =
            build_graph_candidate_for_path(graph, path, start, end, target_span, contribution_model)
        {
            candidates.push(candidate);
        }
    }
    candidates.sort_by(|a, b| {
        a.path_name
            .cmp(&b.path_name)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
    });

    let region_name = match (&target_name, target_interval) {
        (Some(name), Some((start, end))) => format!("{name}:{start}-{end}"),
        (Some(name), None) => name.clone(),
        (None, None) => "all-graph-paths".to_string(),
        (None, Some(_)) => unreachable!(),
    };
    Ok((candidates, region_name, target_name))
}

fn sample_weights_for_graph(
    graph: &NormalizedGraph,
    coverage: &pack::Coverage,
    contribution_model: GraphContributionModel,
) -> io::Result<FxHashMap<u32, f64>> {
    let mut weights = FxHashMap::default();
    for (&feature_id, &count) in &coverage.counts {
        let segment_idx = graph
            .feature_to_segment
            .get(&feature_id)
            .copied()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "sample graph pack contains feature ID {} that is not present in graph {}",
                        feature_id, graph.graph_id
                    ),
                )
            })?;
        let segment_len = graph.segments[segment_idx].length.max(1) as f64;
        let weight = match contribution_model {
            GraphContributionModel::Raw => count as f64,
            GraphContributionModel::LengthNormalized => count as f64 / segment_len,
        };
        if weight > 0.0 {
            weights.insert(feature_id, weight);
        }
    }
    Ok(weights)
}

fn rank_graph_cosigt(
    candidates: &mut Vec<GraphHaplotypeCandidate>,
    sample_weights: &FxHashMap<u32, f64>,
    ploidy: usize,
    top_n: usize,
    candidate_top_k: usize,
    max_combinations: u64,
) -> io::Result<(Vec<CosigtResult>, Vec<u32>)> {
    let all_locus_features = genotyping::feature_universe_f64(
        candidates
            .iter()
            .map(|candidate| candidate.scoring_features.as_slice()),
    );
    let all_sample_norm_sq =
        genotyping::sample_norm_sq_for_weight_features(sample_weights, &all_locus_features);
    for candidate in candidates.iter_mut() {
        candidate.single_similarity = genotyping::cosine_for_feature_weights(
            &candidate.scoring_features,
            sample_weights,
            all_sample_norm_sq,
        );
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

    let selected_features = genotyping::feature_universe_f64(
        candidates
            .iter()
            .map(|candidate| candidate.scoring_features.as_slice()),
    );
    let sample_norm_sq =
        genotyping::sample_norm_sq_for_weight_features(sample_weights, &selected_features);
    if sample_norm_sq == 0.0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "sample graph pack has zero coverage over candidate graph features",
        ));
    }

    let candidate_features: Vec<&[(u32, f64)]> = candidates
        .iter()
        .map(|candidate| candidate.scoring_features.as_slice())
        .collect();
    let mut results = genotyping::run_cosine_combination_search_f64(
        &candidate_features,
        sample_weights,
        sample_norm_sq,
        ploidy,
        max_combinations,
    )?;
    results.truncate(top_n);
    Ok((results, selected_features))
}

fn merge_pack_metadata(target: &mut GraphPackMetadata, key: &str, value: &str, source: &str) {
    let value = value.trim().to_string();
    match key.trim() {
        "feature_space" | "pack_feature_space" => {
            target.feature_space.get_or_insert(value);
        }
        "graph_id" | "pack_graph_id" => {
            target.graph_id.get_or_insert(value);
        }
        "feature_id_mode" | "graph_feature_id_mode" | "pack_feature_id_mode" => {
            target.feature_id_mode.get_or_insert(value);
        }
        "contribution_model" | "graph_contribution_model" | "pack_contribution_model" => {
            target.contribution_model.get_or_insert(value);
        }
        _ => return,
    }
    target.source.get_or_insert_with(|| source.to_string());
}

fn read_pack_metadata_sidecar(path: &Path, metadata: &mut GraphPackMetadata) -> io::Result<()> {
    let path_str = path.to_string_lossy();
    let candidates = [
        PathBuf::from(format!("{path_str}.meta.tsv")),
        PathBuf::from(format!("{path_str}.metadata.tsv")),
    ];
    for candidate in candidates {
        if !candidate.exists() {
            continue;
        }
        let file = File::open(&candidate)?;
        for line in BufReader::new(file).lines() {
            let line = line?;
            let line = line.trim();
            if line.is_empty() || line.starts_with('#') {
                continue;
            }
            if let Some((key, value)) = line.split_once('\t') {
                merge_pack_metadata(metadata, key, value, &candidate.to_string_lossy());
            }
        }
    }
    Ok(())
}

fn read_pack_metadata_inline(path: &Path, metadata: &mut GraphPackMetadata) -> io::Result<()> {
    let mut file = File::open(path)?;
    let mut buf = [0u8; 8192];
    let n = file.read(&mut buf)?;
    if n >= pack::BINARY_MAGIC.len() && &buf[..pack::BINARY_MAGIC.len()] == pack::BINARY_MAGIC {
        return Ok(());
    }
    let Ok(prefix) = std::str::from_utf8(&buf[..n]) else {
        return Ok(());
    };
    for line in prefix.lines() {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        if !line.starts_with('#') {
            break;
        }
        let line = line.trim_start_matches('#');
        if let Some((key, value)) = line.split_once('\t') {
            merge_pack_metadata(metadata, key, value, &path.to_string_lossy());
        }
    }
    Ok(())
}

fn read_graph_pack_metadata(path: &str) -> io::Result<GraphPackMetadata> {
    let path = Path::new(path);
    let mut metadata = GraphPackMetadata::default();
    read_pack_metadata_inline(path, &mut metadata)?;
    read_pack_metadata_sidecar(path, &mut metadata)?;
    Ok(metadata)
}

fn validate_graph_pack_compatibility(
    graph: &NormalizedGraph,
    coverage: &pack::Coverage,
    pack_path: &str,
    pack_feature_space: Option<&str>,
    pack_graph_id: Option<&str>,
    contribution_model: GraphContributionModel,
) -> io::Result<GraphPackMetadata> {
    let mut metadata = read_graph_pack_metadata(pack_path)?;
    if let Some(feature_space) = pack_feature_space {
        metadata.feature_space = Some(feature_space.to_string());
        metadata.source = Some("cli".to_string());
    }
    if let Some(graph_id) = pack_graph_id {
        metadata.graph_id = Some(graph_id.to_string());
        metadata.source = Some("cli".to_string());
    }

    let feature_space = metadata.feature_space.as_deref().ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            "graph genotype packs must declare feature_space=gfa-segment or variation-graph-node in #feature_space TSV metadata, a .meta.tsv sidecar, or --pack-feature-space",
        )
    })?;
    if !matches!(feature_space, "gfa-segment" | "variation-graph-node") {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "pack feature_space '{}' is incompatible with graph genotype feature_space '{}'",
                feature_space, graph.feature_space
            ),
        ));
    }
    if let Some(graph_id) = metadata.graph_id.as_deref() {
        if graph_id != graph.graph_id {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "pack graph_id '{}' does not match loaded graph_id '{}'",
                    graph_id, graph.graph_id
                ),
            ));
        }
    }
    if let Some(mode) = metadata.feature_id_mode.as_deref() {
        if mode != graph.effective_feature_id_mode.as_str() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "pack feature_id_mode '{}' does not match graph feature_id_mode '{}'",
                    mode,
                    graph.effective_feature_id_mode.as_str()
                ),
            ));
        }
    }
    if let Some(model) = metadata.contribution_model.as_deref() {
        if model != contribution_model.as_str() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "pack graph_contribution_model '{}' does not match requested graph_contribution_model '{}'",
                    model,
                    contribution_model.as_str()
                ),
            ));
        }
    }
    for feature_id in coverage.counts.keys() {
        if !graph.feature_ids.contains(feature_id) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "pack feature ID {} is not present in the loaded graph feature universe",
                    feature_id
                ),
            ));
        }
    }
    Ok(metadata)
}

pub fn run_graph_cosigt<W: Write>(out: &mut W, config: &GraphCosigtConfig<'_>) -> io::Result<()> {
    let output = compute_graph_cosigt(config)?;
    write_graph_cosigt_output(out, &output)?;
    if let Some(path) = config.emit_report_path {
        let mut report = BufWriter::new(File::create(path)?);
        write_graph_cosigt_report(&mut report, config, &output)?;
        report.flush()?;
    }
    Ok(())
}

pub fn compute_graph_cosigt(config: &GraphCosigtConfig<'_>) -> io::Result<GraphCosigtOutput> {
    if config.ploidy == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--ploidy must be greater than 0",
        ));
    }
    if config.top_n == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--top-n must be greater than 0",
        ));
    }
    let loaded = load_graph_source(&config.graph_source)?;
    let target_path = config
        .target_path
        .map(ToString::to_string)
        .or(loaded.default_target_path.clone());
    let graph = parse_normalized_gfa(
        &loaded.gfa_text,
        loaded.source_label,
        loaded.source_path,
        loaded.build_command,
        loaded.feature_space,
        config.feature_id_mode,
    )?;
    info!(
        "Loaded graph {} with {} segment(s), {} path(s), feature_id_mode={}",
        graph.graph_id,
        graph.segments.len(),
        graph.paths.len(),
        graph.effective_feature_id_mode.as_str()
    );
    let coverage = pack::read(config.pack_path)?;
    let pack_metadata = validate_graph_pack_compatibility(
        &graph,
        &coverage,
        config.pack_path,
        config.pack_feature_space,
        config.pack_graph_id,
        config.contribution_model,
    )?;
    let sample_weights = sample_weights_for_graph(&graph, &coverage, config.contribution_model)?;
    let (mut candidates, region_name, resolved_target_path) =
        collect_graph_candidates(&graph, target_path.as_deref(), config.contribution_model)?;
    if candidates.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no graph genotype candidates were found",
        ));
    }
    let (results, selected_features) = rank_graph_cosigt(
        &mut candidates,
        &sample_weights,
        config.ploidy,
        config.top_n,
        config.candidate_top_k,
        config.max_combinations,
    )?;
    Ok(GraphCosigtOutput {
        region_name,
        target_path: resolved_target_path,
        ploidy: config.ploidy,
        contribution_model: config.contribution_model,
        graph,
        candidates,
        selected_features,
        sample_raw_counts: coverage.counts.clone(),
        sample_weights,
        pack_metadata,
        pack_nonzero_nodes: coverage.nonzero_nodes,
        pack_universe_nodes: coverage.universe_nodes,
        pack_retained_records: coverage.retained_records,
        pack_syncmer_anchors: coverage.syncmer_anchors,
        results,
    })
}

fn format_graph_candidate_region(candidate: &GraphHaplotypeCandidate) -> String {
    format!(
        "{}:{}-{}({})",
        candidate.path_name, candidate.start, candidate.end, candidate.strand
    )
}

fn graph_candidate_count_for_feature(candidate: &GraphHaplotypeCandidate, node_id: u32) -> u64 {
    candidate
        .features
        .binary_search_by_key(&node_id, |&(feature_id, _)| feature_id)
        .map(|idx| candidate.features[idx].1)
        .unwrap_or(0)
}

fn graph_candidate_weight_for_feature(candidate: &GraphHaplotypeCandidate, node_id: u32) -> f64 {
    candidate
        .scoring_features
        .binary_search_by_key(&node_id, |&(feature_id, _)| feature_id)
        .map(|idx| candidate.scoring_features[idx].1)
        .unwrap_or(0.0)
}

fn graph_segment_by_feature(graph: &NormalizedGraph, feature_id: u32) -> io::Result<&GraphSegment> {
    graph
        .feature_to_segment
        .get(&feature_id)
        .map(|&idx| &graph.segments[idx])
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("feature ID {} is not present in graph", feature_id),
            )
        })
}

pub fn write_graph_cosigt_output<W: Write>(
    out: &mut W,
    result: &GraphCosigtOutput,
) -> io::Result<()> {
    writeln!(out, "#impg genotype cos")?;
    writeln!(out, "#region\t{}", result.region_name)?;
    writeln!(out, "#method\t{}", ScoringMethod::Cos.as_str())?;
    writeln!(out, "#metric\t{}", ScoringMethod::Cos.metric_str())?;
    writeln!(out, "#alias\tcosigt")?;
    writeln!(out, "#feature_space\t{}", result.graph.feature_space)?;
    writeln!(out, "#graph_source\t{}", result.graph.source_label)?;
    writeln!(out, "#graph_id\t{}", result.graph.graph_id)?;
    writeln!(
        out,
        "#graph_feature_id_mode\t{}",
        result.graph.effective_feature_id_mode.as_str()
    )?;
    writeln!(
        out,
        "#graph_contribution_model\t{}",
        result.contribution_model.as_str()
    )?;
    if let Some(path) = &result.target_path {
        writeln!(out, "#target_path\t{}", path)?;
    }
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
            .map(|&idx| format_graph_candidate_region(&result.candidates[idx]))
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
            "{}\t{}\t{}\t{:.9}\t{:.3}\t{:.3}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}",
            rank + 1,
            ScoringMethod::Cos.as_str(),
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

pub fn write_graph_cosigt_report<W: Write>(
    out: &mut W,
    config: &GraphCosigtConfig<'_>,
    result: &GraphCosigtOutput,
) -> io::Result<()> {
    writeln!(out, "#impg genotype cos report")?;
    writeln!(out, "#format\tsectioned-tsv-v1")?;

    writeln!(out, "#section\tinput_metadata")?;
    writeln!(out, "key\tvalue")?;
    write_report_kv(out, "graph_source", &result.graph.source_label)?;
    write_report_kv(
        out,
        "graph_source_path",
        result
            .graph
            .source_path
            .as_ref()
            .map(|p| p.to_string_lossy().to_string())
            .unwrap_or_else(|| "NA".to_string()),
    )?;
    write_report_kv(
        out,
        "graph_build_command",
        result.graph.build_command.as_deref().unwrap_or("NA"),
    )?;
    write_report_kv(out, "graph_id", &result.graph.graph_id)?;
    write_report_kv(out, "feature_space", &result.graph.feature_space)?;
    write_report_kv(
        out,
        "requested_feature_id_mode",
        result.graph.requested_feature_id_mode.as_str(),
    )?;
    write_report_kv(
        out,
        "effective_feature_id_mode",
        result.graph.effective_feature_id_mode.as_str(),
    )?;
    write_report_kv(out, "pack_path", config.pack_path)?;
    write_report_kv(
        out,
        "pack_declared_feature_space",
        result
            .pack_metadata
            .feature_space
            .as_deref()
            .unwrap_or("NA"),
    )?;
    write_report_kv(
        out,
        "pack_declared_graph_id",
        result.pack_metadata.graph_id.as_deref().unwrap_or("NA"),
    )?;
    write_report_kv(
        out,
        "pack_declared_contribution_model",
        result
            .pack_metadata
            .contribution_model
            .as_deref()
            .unwrap_or("NA"),
    )?;
    write_report_kv(
        out,
        "pack_metadata_source",
        result.pack_metadata.source.as_deref().unwrap_or("NA"),
    )?;
    write_report_kv(out, "resolved_region", &result.region_name)?;
    write_report_kv(
        out,
        "target_path",
        result.target_path.as_deref().unwrap_or("NA"),
    )?;
    write_report_kv(out, "method", ScoringMethod::Cos.as_str())?;
    write_report_kv(out, "metric", ScoringMethod::Cos.metric_str())?;
    write_report_kv(
        out,
        "contribution_model",
        result.contribution_model.as_str(),
    )?;
    write_report_kv(out, "ploidy", result.ploidy)?;
    write_report_kv(out, "top_n", config.top_n)?;
    write_report_kv(out, "candidate_top_k", config.candidate_top_k)?;
    write_report_kv(out, "max_combinations", config.max_combinations)?;
    write_report_kv(out, "graph_segments", result.graph.segments.len())?;
    write_report_kv(out, "graph_paths", result.graph.paths.len())?;
    write_report_kv(out, "candidate_count", result.candidates.len())?;
    write_report_kv(out, "pack_nonzero_nodes", result.pack_nonzero_nodes)?;
    write_report_kv(
        out,
        "pack_universe_nodes",
        format_optional_u64(result.pack_universe_nodes),
    )?;
    write_report_kv(
        out,
        "pack_retained_records",
        format_optional_u64(result.pack_retained_records),
    )?;
    write_report_kv(
        out,
        "pack_syncmer_anchors",
        format_optional_u64(result.pack_syncmer_anchors),
    )?;

    writeln!(out, "#section\tgraph_feature_universe")?;
    writeln!(
        out,
        "feature_id\tsegment_name\tsegment_length\tsample_raw_count\tsample_weight"
    )?;
    for segment in &result.graph.segments {
        let raw = result
            .sample_raw_counts
            .get(&segment.feature_id)
            .copied()
            .unwrap_or(0);
        let weight = result
            .sample_weights
            .get(&segment.feature_id)
            .copied()
            .unwrap_or(0.0);
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{:.9}",
            segment.feature_id, segment.name, segment.length, raw, weight
        )?;
    }

    writeln!(out, "#section\tpack_evidence_summary")?;
    writeln!(out, "metric\tvalue")?;
    write_report_kv(
        out,
        "selected_locus_features",
        result.selected_features.len(),
    )?;
    let sample_locus_weight_mass: f64 = result
        .selected_features
        .iter()
        .map(|node_id| result.sample_weights.get(node_id).copied().unwrap_or(0.0))
        .sum();
    let sample_locus_overlap_features = result
        .selected_features
        .iter()
        .filter(|node_id| result.sample_weights.get(node_id).copied().unwrap_or(0.0) > 0.0)
        .count();
    let sample_norm_sq = genotyping::sample_norm_sq_for_weight_features(
        &result.sample_weights,
        &result.selected_features,
    );
    write_report_kv(
        out,
        "locus_feature_overlap_nonzero_nodes",
        sample_locus_overlap_features,
    )?;
    write_report_kv(
        out,
        "locus_feature_overlap_sample_weight_mass",
        format!("{sample_locus_weight_mass:.9}"),
    )?;
    write_report_kv(
        out,
        "sample_norm_over_selected_locus_features",
        format!("{:.6}", sample_norm_sq.sqrt()),
    )?;

    writeln!(out, "#section\tsample_locus_features")?;
    writeln!(
        out,
        "feature_id\tsegment_name\tsegment_length\tsample_weight"
    )?;
    for &feature_id in &result.selected_features {
        let segment = graph_segment_by_feature(&result.graph, feature_id)?;
        writeln!(
            out,
            "{}\t{}\t{}\t{:.9}",
            feature_id,
            segment.name,
            segment.length,
            result
                .sample_weights
                .get(&feature_id)
                .copied()
                .unwrap_or(0.0)
        )?;
    }

    writeln!(out, "#section\tcandidates")?;
    writeln!(
        out,
        "candidate_index\tpath\tinterval\tstart\tend\tstrand\tpath_steps\tspan_fraction\tfeature_count\ttotal_candidate_node_count_mass\tunique_nodes\trepeated_nodes\trepeated_extra_count\tmax_node_count\tsingle_haplotype_cosine\tsample_overlap_unique_nodes\tsample_overlap_weight_mass\tsample_overlap_dot_contribution"
    )?;
    for (idx, candidate) in result.candidates.iter().enumerate() {
        let total_mass: u64 = candidate.features.iter().map(|&(_, count)| count).sum();
        let repeated_nodes = candidate
            .features
            .iter()
            .filter(|&&(_, count)| count > 1)
            .count();
        let repeated_extra_count: u64 = candidate
            .features
            .iter()
            .map(|&(_, count)| count.saturating_sub(1))
            .sum();
        let max_node_count = candidate
            .features
            .iter()
            .map(|&(_, count)| count)
            .max()
            .unwrap_or(0);
        let mut overlap_unique_nodes = 0usize;
        let mut overlap_weight_mass = 0.0;
        let mut overlap_dot = 0.0;
        for &(feature_id, candidate_weight) in &candidate.scoring_features {
            let sample_weight = result
                .sample_weights
                .get(&feature_id)
                .copied()
                .unwrap_or(0.0);
            if sample_weight > 0.0 {
                overlap_unique_nodes += 1;
                overlap_weight_mass += sample_weight;
                overlap_dot += sample_weight * candidate_weight;
            }
        }
        writeln!(
            out,
            "{}\t{}\t{}:{}-{}\t{}\t{}\t{}\t{}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{}\t{:.9}\t{:.9}",
            idx,
            candidate.path_name,
            candidate.path_name,
            candidate.start,
            candidate.end,
            candidate.start,
            candidate.end,
            candidate.strand,
            candidate.anchors,
            candidate.query_span_fraction,
            candidate.features.len(),
            total_mass,
            candidate.features.len(),
            repeated_nodes,
            repeated_extra_count,
            max_node_count,
            candidate.single_similarity,
            overlap_unique_nodes,
            overlap_weight_mass,
            overlap_dot,
        )?;
    }

    writeln!(out, "#section\tcandidate_features")?;
    writeln!(
        out,
        "candidate_index\tpath\tfeature_id\tsegment_name\tsegment_length\tsample_weight\tcandidate_count\tcandidate_weight\tdot_contribution"
    )?;
    for (idx, candidate) in result.candidates.iter().enumerate() {
        for &(feature_id, candidate_weight) in &candidate.scoring_features {
            let segment = graph_segment_by_feature(&result.graph, feature_id)?;
            let sample_weight = result
                .sample_weights
                .get(&feature_id)
                .copied()
                .unwrap_or(0.0);
            let candidate_count = graph_candidate_count_for_feature(candidate, feature_id);
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{:.9}\t{}\t{:.9}\t{:.9}",
                idx,
                candidate.path_name,
                feature_id,
                segment.name,
                segment.length,
                sample_weight,
                candidate_count,
                candidate_weight,
                sample_weight * candidate_weight,
            )?;
        }
    }

    writeln!(out, "#section\tresult_scores")?;
    writeln!(
        out,
        "rank\tmethod\tploidy\tsimilarity\tqv\tdot\tsample_norm\tgenotype_norm\tcandidate_indices\thaplotypes\tregions\tcandidate_anchors\tcandidate_span_fractions"
    )?;
    for (rank, genotype_result) in result.results.iter().enumerate() {
        let candidate_indices: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|idx| idx.to_string())
            .collect();
        let haplotypes: Vec<&str> = genotype_result
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].path_name.as_str())
            .collect();
        let regions: Vec<String> = genotype_result
            .combination
            .iter()
            .map(|&idx| format_graph_candidate_region(&result.candidates[idx]))
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
            "{}\t{}\t{}\t{:.9}\t{:.3}\t{:.3}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}",
            rank + 1,
            ScoringMethod::Cos.as_str(),
            result.ploidy,
            genotype_result.similarity,
            genotype_result.qv,
            genotype_result.dot,
            genotype_result.sample_norm,
            genotype_result.genotype_norm,
            candidate_indices.join(","),
            haplotypes.join(","),
            regions.join(","),
            anchors.join(","),
            spans.join(","),
        )?;
    }

    writeln!(out, "#section\tresult_features")?;
    writeln!(
        out,
        "rank\tfeature_id\tsegment_name\tsegment_length\tsample_weight\tcandidate_weights\tgenotype_weight\tdot_contribution"
    )?;
    for (rank, genotype_result) in result.results.iter().enumerate() {
        let mut genotype_weights = BTreeMap::new();
        for &candidate_idx in &genotype_result.combination {
            for &(feature_id, weight) in &result.candidates[candidate_idx].scoring_features {
                *genotype_weights.entry(feature_id).or_insert(0.0) += weight;
            }
        }
        for (feature_id, genotype_weight) in genotype_weights {
            let segment = graph_segment_by_feature(&result.graph, feature_id)?;
            let sample_weight = result
                .sample_weights
                .get(&feature_id)
                .copied()
                .unwrap_or(0.0);
            let candidate_weights: Vec<String> = genotype_result
                .combination
                .iter()
                .map(|&candidate_idx| {
                    format!(
                        "{:.9}",
                        graph_candidate_weight_for_feature(
                            &result.candidates[candidate_idx],
                            feature_id
                        )
                    )
                })
                .collect();
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{:.9}\t{}\t{:.9}\t{:.9}",
                rank + 1,
                feature_id,
                segment.name,
                segment.length,
                sample_weight,
                candidate_weights.join(","),
                genotype_weight,
                sample_weight * genotype_weight,
            )?;
        }
    }

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    const TWO_HAP_GFA: &str = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
P\th1\t1+,2+,4+\t*
P\th2\t1+,3+,4+\t*
";

    fn write_pack(dir: &Path, rows: &[(u32, u64)], with_metadata: bool) -> PathBuf {
        let path = dir.join("sample.pack.tsv");
        let mut text = String::new();
        if with_metadata {
            text.push_str("#feature_space\tgfa-segment\n");
            text.push_str("#feature_id_mode\tsegment-name\n");
        }
        text.push_str("#node_id\tcount\n");
        for (node_id, count) in rows {
            text.push_str(&format!("{node_id}\t{count}\n"));
        }
        std::fs::write(&path, text).unwrap();
        path
    }

    fn graph_config<'a>(
        gfa: &'a str,
        pack_path: &'a str,
        contribution_model: GraphContributionModel,
    ) -> GraphCosigtConfig<'a> {
        GraphCosigtConfig {
            graph_source: GraphSource::Query {
                gfa,
                description: "unit-test",
                build_command: Some("impg query -o gfa ..."),
                default_target_path: None,
            },
            pack_path,
            target_path: None,
            pack_feature_space: None,
            pack_graph_id: None,
            feature_id_mode: GraphFeatureIdMode::SegmentName,
            contribution_model,
            emit_report_path: None,
            ploidy: 2,
            top_n: 5,
            candidate_top_k: 0,
            max_combinations: 100,
        }
    }

    #[test]
    fn parses_p_line_paths_and_extracts_intervals() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCC
S\t3\tG
P\tref\t1+,2+,3+\t*
P\talt\t1+,3+\t*
";
        let graph = parse_normalized_gfa(
            gfa,
            "test",
            None,
            None,
            FeatureSpace::GfaSegment.as_str(),
            GraphFeatureIdMode::SegmentName,
        )
        .unwrap();
        assert_eq!(graph.segments.len(), 3);
        assert_eq!(graph.paths[graph.path_name_to_index["ref"]].length, 7);

        let (candidates, region, target) =
            collect_graph_candidates(&graph, Some("ref:2-6"), GraphContributionModel::Raw).unwrap();
        assert_eq!(region, "ref:2-6");
        assert_eq!(target.as_deref(), Some("ref"));
        let ref_candidate = candidates
            .iter()
            .find(|candidate| candidate.path_name == "ref")
            .unwrap();
        assert_eq!(ref_candidate.features, vec![(1, 1), (2, 1)]);
        assert_eq!(ref_candidate.anchors, 2);
        assert!((ref_candidate.query_span_fraction - 1.0).abs() < 1e-9);
    }

    #[test]
    fn star_segments_use_ln_tags_for_candidate_intervals_and_weights() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\t*\tLN:i:4
S\t2\t*\tLN:i:6
S\t3\tA
P\tref\t1+,2+,3+\t*
P\talt\t1+,3+\t*
";
        let graph = parse_normalized_gfa(
            gfa,
            "test",
            None,
            None,
            FeatureSpace::GfaSegment.as_str(),
            GraphFeatureIdMode::SegmentName,
        )
        .unwrap();
        assert_eq!(graph.segments[0].length, 4);
        assert_eq!(graph.segments[1].length, 6);
        assert_eq!(graph.paths[graph.path_name_to_index["ref"]].length, 11);

        let (candidates, region, target) = collect_graph_candidates(
            &graph,
            Some("ref:3-7"),
            GraphContributionModel::LengthNormalized,
        )
        .unwrap();
        assert_eq!(region, "ref:3-7");
        assert_eq!(target.as_deref(), Some("ref"));

        let ref_candidate = candidates
            .iter()
            .find(|candidate| candidate.path_name == "ref")
            .unwrap();
        assert_eq!(ref_candidate.features, vec![(1, 1), (2, 1)]);
        assert_eq!(ref_candidate.anchors, 2);
        assert_eq!(ref_candidate.scoring_features.len(), 2);
        assert!((ref_candidate.scoring_features[0].1 - 0.25).abs() < 1e-12);
        assert!((ref_candidate.scoring_features[1].1 - 0.5).abs() < 1e-12);
    }

    #[test]
    fn repeated_gfa_path_steps_are_counted_in_candidate_vectors() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCC
S\t3\tG
P\tsingle\t1+,2+,3+\t*
P\tduplicated\t1+,2+,2+,3+\t*
";
        let graph = parse_normalized_gfa(
            gfa,
            "test",
            None,
            None,
            FeatureSpace::GfaSegment.as_str(),
            GraphFeatureIdMode::SegmentName,
        )
        .unwrap();

        let (candidates, _, _) =
            collect_graph_candidates(&graph, None, GraphContributionModel::Raw).unwrap();
        let duplicated = candidates
            .iter()
            .find(|candidate| candidate.path_name == "duplicated")
            .unwrap();

        assert_eq!(duplicated.features, vec![(1, 1), (2, 2), (3, 1)]);
        assert_eq!(
            duplicated.scoring_features,
            vec![(1, 1.0), (2, 2.0), (3, 1.0)]
        );
        assert_eq!(duplicated.anchors, 4);
    }

    #[test]
    fn parses_w_line_walks_as_pansn_paths() {
        let gfa = "\
H\tVN:Z:1.1
S\t1\tAC
S\t2\tGGG
S\t3\tTA
W\tsample\t0\tref\t0\t4\t>1>3
W\tsample\t0\tins\t0\t7\t>1>2>3
";
        let graph = parse_normalized_gfa(
            gfa,
            "test",
            None,
            None,
            FeatureSpace::GfaSegment.as_str(),
            GraphFeatureIdMode::SegmentName,
        )
        .unwrap();
        assert!(graph.path_name_to_index.contains_key("sample#0#ref"));
        assert!(graph.path_name_to_index.contains_key("sample#0#ins"));

        let (candidates, region, _) = collect_graph_candidates(
            &graph,
            Some("sample#0#ins:2-5"),
            GraphContributionModel::Raw,
        )
        .unwrap();
        assert_eq!(region, "sample#0#ins:2-5");
        let ins = candidates
            .iter()
            .find(|candidate| candidate.path_name == "sample#0#ins")
            .unwrap();
        assert_eq!(ins.features, vec![(2, 1)]);
        assert_eq!(ins.anchors, 1);
    }

    #[test]
    fn synthetic_graph_calls_expected_heterozygote() {
        let dir = tempfile::tempdir().unwrap();
        let pack = write_pack(dir.path(), &[(1, 2), (2, 1), (3, 1), (4, 2)], true);
        let config = graph_config(
            TWO_HAP_GFA,
            pack.to_str().unwrap(),
            GraphContributionModel::Raw,
        );
        let result = compute_graph_cosigt(&config).unwrap();
        let best = &result.results[0];
        let best_haps: Vec<&str> = best
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].path_name.as_str())
            .collect();
        assert_eq!(best_haps, vec!["h1", "h2"]);
        assert!((best.similarity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn synthetic_graph_calls_expected_homozygote() {
        let dir = tempfile::tempdir().unwrap();
        let pack = write_pack(dir.path(), &[(1, 2), (3, 2), (4, 2)], true);
        let config = graph_config(
            TWO_HAP_GFA,
            pack.to_str().unwrap(),
            GraphContributionModel::Raw,
        );
        let result = compute_graph_cosigt(&config).unwrap();
        let best = &result.results[0];
        let best_haps: Vec<&str> = best
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].path_name.as_str())
            .collect();
        assert_eq!(best_haps, vec!["h2", "h2"]);
        assert!((best.similarity - 1.0).abs() < 1e-9);
    }

    #[test]
    fn prebuilt_gfa_and_query_graph_sources_are_equivalent() {
        let dir = tempfile::tempdir().unwrap();
        let gfa_path = dir.path().join("locus.gfa");
        std::fs::write(&gfa_path, TWO_HAP_GFA).unwrap();
        let pack = write_pack(dir.path(), &[(1, 2), (2, 1), (3, 1), (4, 2)], true);
        let pack_path = pack.to_str().unwrap();
        let prebuilt = GraphCosigtConfig {
            graph_source: GraphSource::Gfa {
                path: gfa_path.to_str().unwrap(),
            },
            pack_path,
            target_path: None,
            pack_feature_space: None,
            pack_graph_id: None,
            feature_id_mode: GraphFeatureIdMode::SegmentName,
            contribution_model: GraphContributionModel::Raw,
            emit_report_path: None,
            ploidy: 2,
            top_n: 5,
            candidate_top_k: 0,
            max_combinations: 100,
        };
        let query = graph_config(TWO_HAP_GFA, pack_path, GraphContributionModel::Raw);

        let prebuilt_result = compute_graph_cosigt(&prebuilt).unwrap();
        let query_result = compute_graph_cosigt(&query).unwrap();

        assert_eq!(
            prebuilt_result.results[0].combination,
            query_result.results[0].combination
        );
        assert!(
            (prebuilt_result.results[0].similarity - query_result.results[0].similarity).abs()
                < 1e-12
        );
    }

    #[test]
    fn render_bundle_graph_source_calls_expected_heterozygote() {
        use crate::render_bundle::{
            write_manifest, write_translation_binary, RenderBundlePaths, RenderManifest,
            RenderTranslationTables, RenderedPathRecord, StepTranslationRecord,
        };
        use crate::sequence_namespace::{SequenceNamespace, SourceInterval};

        let dir = tempfile::tempdir().unwrap();
        let paths = RenderBundlePaths::new(dir.path());
        std::fs::write(&paths.graph_gfa, TWO_HAP_GFA).unwrap();
        let mut namespace = SequenceNamespace::new();
        let source_id = namespace.add_sequence("source", 3);
        let tables = RenderTranslationTables {
            namespace,
            rendered_paths: vec![RenderedPathRecord {
                rendered_path_id: 0,
                rendered_name: "h1".to_string(),
                source_interval: SourceInterval::new(source_id, 0, 3, '+').unwrap(),
                gbwt_path_id: None,
            }],
            step_samples: vec![
                StepTranslationRecord {
                    rendered_path_id: 0,
                    rendered_step: 0,
                    source_bp: 0,
                    feature_id: 1,
                    orientation: '+',
                },
                StepTranslationRecord {
                    rendered_path_id: 0,
                    rendered_step: 1,
                    source_bp: 1,
                    feature_id: 2,
                    orientation: '+',
                },
                StepTranslationRecord {
                    rendered_path_id: 0,
                    rendered_step: 2,
                    source_bp: 2,
                    feature_id: 4,
                    orientation: '+',
                },
            ],
        };
        let manifest = RenderManifest::new_local_graph(
            "unit-test".to_string(),
            "source-index".to_string(),
            "source:0-3".to_string(),
            &paths,
            1,
            1,
            tables.step_samples.len(),
        );
        write_manifest(&paths.manifest, &manifest).unwrap();
        write_translation_binary(&paths.translation, &tables).unwrap();

        let pack = write_pack(dir.path(), &[(1, 2), (2, 1), (3, 1), (4, 2)], true);
        let config = GraphCosigtConfig {
            graph_source: GraphSource::RenderBundle {
                path: dir.path().to_str().unwrap(),
            },
            pack_path: pack.to_str().unwrap(),
            target_path: None,
            pack_feature_space: None,
            pack_graph_id: None,
            feature_id_mode: GraphFeatureIdMode::SegmentName,
            contribution_model: GraphContributionModel::Raw,
            emit_report_path: None,
            ploidy: 2,
            top_n: 5,
            candidate_top_k: 0,
            max_combinations: 100,
        };

        let result = compute_graph_cosigt(&config).unwrap();
        assert_eq!(result.region_name, "h1:0-3");
        let best_haps: Vec<&str> = result.results[0]
            .combination
            .iter()
            .map(|&idx| result.candidates[idx].path_name.as_str())
            .collect();
        assert_eq!(best_haps, vec!["h1", "h2"]);
    }

    #[test]
    fn length_normalized_model_changes_variable_length_node_weights() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tAAAAAAAAAA
S\t3\tC
P\tlong\t2+,3+\t*
P\tshort\t1+,3+\t*
";
        let dir = tempfile::tempdir().unwrap();
        let pack = write_pack(dir.path(), &[(1, 1), (2, 1), (3, 1)], true);
        let mut raw = graph_config(gfa, pack.to_str().unwrap(), GraphContributionModel::Raw);
        raw.ploidy = 1;
        let raw_result = compute_graph_cosigt(&raw).unwrap();
        let mut normalized = graph_config(
            gfa,
            pack.to_str().unwrap(),
            GraphContributionModel::LengthNormalized,
        );
        normalized.ploidy = 1;
        let normalized_result = compute_graph_cosigt(&normalized).unwrap();

        assert_eq!(raw_result.sample_weights[&2], 1.0);
        assert!((normalized_result.sample_weights[&2] - 0.1).abs() < 1e-12);
        let raw_best = raw_result.candidates[raw_result.results[0].combination[0]]
            .path_name
            .as_str();
        let normalized_best = normalized_result.candidates
            [normalized_result.results[0].combination[0]]
            .path_name
            .as_str();
        assert_eq!(raw_best, "long");
        assert_eq!(normalized_best, "short");
    }

    #[test]
    fn graph_pack_requires_feature_space_declaration() {
        let dir = tempfile::tempdir().unwrap();
        let pack = write_pack(dir.path(), &[(1, 2), (2, 1), (3, 1), (4, 2)], false);
        let config = graph_config(
            TWO_HAP_GFA,
            pack.to_str().unwrap(),
            GraphContributionModel::Raw,
        );
        let err = compute_graph_cosigt(&config).unwrap_err();
        assert_eq!(err.kind(), io::ErrorKind::InvalidInput);
        assert!(err
            .to_string()
            .contains("graph genotype packs must declare"));
    }
}
