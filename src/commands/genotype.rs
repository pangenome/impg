use clap::ValueEnum;
use log::info;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::fs::File;
use std::io::{self, BufWriter, Write};

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
