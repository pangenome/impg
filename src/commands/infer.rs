use clap::ValueEnum;
use log::info;
use rustc_hash::FxHashMap;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::commands::{genotype, graph, partition};
use crate::genotyping::{EvidenceBackend, FeatureSpace, ScoringMethod};
use crate::graph::reverse_complement;
use crate::pack;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::syng;
use crate::{EngineOpts, GfaEngine, SyngImpgWrapper};

#[derive(Debug, Clone)]
pub struct InferTarget {
    pub partition: String,
    pub chrom: String,
    pub start: i32,
    pub end: i32,
    pub name: String,
}

#[derive(Debug, Clone)]
pub struct PartitionDiscoveryConfig<'a> {
    pub window_size: usize,
    pub merge_distance: i32,
    pub starting_sequences_file: Option<&'a str>,
    pub selection_mode: &'a str,
    pub min_missing_size: i32,
    pub min_boundary_distance: i32,
    pub syng_padding: u64,
    pub syng_min_chain_anchors: usize,
    pub syng_min_chain_fraction: f64,
    pub rehome_singletons: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, ValueEnum)]
pub enum StitchMode {
    /// Emit independent local calls only.
    None,
    /// Stitch retained local calls into phased mosaics with a beam search.
    Beam,
}

pub struct InferConfig<'a> {
    pub syng_prefix: &'a str,
    pub pack_path: &'a str,
    pub targets: Vec<InferTarget>,
    pub discovery: Option<PartitionDiscoveryConfig<'a>>,
    pub candidate_mode: genotype::CandidateMode,
    pub ploidy: usize,
    pub top_n: usize,
    pub candidate_top_k: usize,
    pub max_combinations: u64,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub min_anchors: usize,
    pub min_span_fraction: f64,
    pub gaf_path: Option<&'a str>,
    pub stitch: StitchMode,
    pub phase_block_size: usize,
    pub stitch_beam: usize,
    pub switch_penalty: f64,
    pub stitch_gap: u64,
    pub read_link_weight: f64,
    pub min_read_link_anchors: usize,
    pub strict_stitch: bool,
    pub emit_mosaic: Option<&'a str>,
    pub emit_fasta: Option<&'a str>,
    pub emit_gfa: Option<&'a str>,
    pub sequence_index: Option<&'a UnifiedSequenceIndex>,
}

#[derive(Debug)]
pub struct LocalCallSet {
    pub target: InferTarget,
    pub output: Option<genotype::SyngCosigtOutput>,
    pub error: Option<String>,
}

#[derive(Debug, Clone)]
struct OrderedState {
    result_idx: usize,
    candidates: Vec<usize>,
    emission: f64,
    read_emissions: Vec<CandidateEmissionScore>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
struct CandidateKey {
    call_idx: usize,
    candidate_idx: usize,
}

#[derive(Debug, Clone, Copy)]
struct WalkOccurrence {
    key: CandidateKey,
    position: u32,
}

#[derive(Debug, Default, Clone)]
struct LinkEvidence {
    read_weight: f64,
    anchor_weight: f64,
}

#[derive(Debug, Default)]
struct ReadWalkEvidence {
    candidate_support: FxHashMap<CandidateKey, LinkEvidence>,
    links: FxHashMap<(CandidateKey, CandidateKey), LinkEvidence>,
    records: u64,
    records_with_candidate_hits: u64,
    candidate_hits: u64,
    linked_read_pairs: u64,
}

#[derive(Debug, Default, Clone)]
struct CandidateEmissionScore {
    read_weight: f64,
    anchor_weight: f64,
    reward: f64,
}

#[derive(Debug, Clone)]
struct PhaseTransitionScore {
    kind: TransitionKind,
    cost: f64,
    read_link_reads: f64,
    read_link_anchors: f64,
    read_link_reward: f64,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum TransitionKind {
    Start,
    Continue,
    Recombine,
    Uncertain,
}

impl TransitionKind {
    fn as_str(self) -> &'static str {
        match self {
            TransitionKind::Start => "start",
            TransitionKind::Continue => "continue",
            TransitionKind::Recombine => "recombine",
            TransitionKind::Uncertain => "uncertain",
        }
    }
}

#[derive(Debug, Clone)]
struct MosaicStep {
    call_idx: usize,
    state: OrderedState,
    transitions: Vec<PhaseTransitionScore>,
}

#[derive(Debug, Clone)]
struct MosaicPath {
    score: f64,
    steps: Vec<MosaicStep>,
}

fn parse_bed_like(path: &str, group_partitions: bool) -> io::Result<Vec<InferTarget>> {
    let file = File::open(path)?;
    let reader = BufReader::new(file);
    let mut targets = Vec::new();
    let mut first_by_partition: BTreeMap<String, InferTarget> = BTreeMap::new();

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("BED line {} has fewer than 3 fields", line_no + 1),
            ));
        }
        let start = fields[1].parse::<i32>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid BED start on line {}: {}", line_no + 1, e),
            )
        })?;
        let end = fields[2].parse::<i32>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("invalid BED end on line {}: {}", line_no + 1, e),
            )
        })?;
        if start < 0 || end <= start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "invalid BED range on line {}: {}-{}",
                    line_no + 1,
                    start,
                    end
                ),
            ));
        }
        let name = fields
            .get(3)
            .map(|s| s.trim())
            .filter(|s| !s.is_empty() && *s != ".")
            .map(str::to_string)
            .unwrap_or_else(|| format!("{}:{}-{}", fields[0], start, end));
        let partition = if group_partitions {
            name.clone()
        } else {
            format!("{}", targets.len())
        };
        let target = InferTarget {
            partition: partition.clone(),
            chrom: fields[0].to_string(),
            start,
            end,
            name,
        };
        if group_partitions {
            first_by_partition.entry(partition).or_insert(target);
        } else {
            targets.push(target);
        }
    }

    if group_partitions {
        targets.extend(first_by_partition.into_values());
    }
    Ok(targets)
}

pub fn parse_target_range(target_range: &str) -> io::Result<InferTarget> {
    let (chrom, (start, end), name) = partition::parse_target_range(target_range)?;
    Ok(InferTarget {
        partition: "0".to_string(),
        chrom,
        start,
        end,
        name,
    })
}

pub fn parse_target_bed(path: &str) -> io::Result<Vec<InferTarget>> {
    parse_bed_like(path, false)
}

pub fn parse_partitions(path: &str) -> io::Result<Vec<InferTarget>> {
    parse_bed_like(path, true)
}

fn temp_partition_dir() -> io::Result<PathBuf> {
    let stamp = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .map_err(io::Error::other)?
        .as_nanos();
    let dir = std::env::temp_dir().join(format!(
        "impg-infer-partitions-{}-{stamp}",
        std::process::id()
    ));
    fs::create_dir_all(&dir)?;
    Ok(dir)
}

fn discover_partitions(
    syng_prefix: &str,
    discovery: &PartitionDiscoveryConfig<'_>,
) -> io::Result<Vec<InferTarget>> {
    let temp_dir = temp_partition_dir()?;
    let temp_dir_str = temp_dir.to_string_lossy().to_string();
    let syng_index = syng::SyngIndex::load(syng_prefix, syng::SyncmerParams::default())?;
    let seq_index = syng_index.build_seq_index();
    let wrapper = {
        let base = SyngImpgWrapper::new(syng_index, seq_index, discovery.syng_padding);
        if discovery.syng_min_chain_anchors > 0 {
            base.with_chain_filter(
                discovery.syng_min_chain_anchors,
                discovery.syng_min_chain_fraction,
            )
        } else {
            base
        }
    };
    let engine_config = EngineOpts {
        engine: GfaEngine::Poa,
        syng_gfa_mode: None,
        syng_params: None,
        pipeline: graph::GraphBuildConfig::default(),
        target_poa_lengths: vec![700, 1100],
        max_node_length: 100,
        poa_padding_fraction: 0.001,
        partition_size: None,
        crush_config: None,
    };

    info!(
        "Discovering inference partitions with window_size={} merge_distance={}",
        discovery.window_size, discovery.merge_distance
    );
    partition::partition_alignments(
        &wrapper,
        discovery.window_size,
        discovery.starting_sequences_file,
        discovery.selection_mode,
        discovery.merge_distance,
        None,
        discovery.min_missing_size,
        discovery.min_boundary_distance,
        false,
        1,
        0,
        0,
        "bed",
        Some(&temp_dir_str),
        None,
        None,
        false,
        false,
        false,
        false,
        &engine_config,
        discovery.rehome_singletons,
    )?;

    let partitions_path = temp_dir.join("partitions.bed");
    let result = parse_partitions(partitions_path.to_str().ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidInput, "temporary path is not UTF-8")
    })?);
    fs::remove_dir_all(&temp_dir).ok();
    result
}

fn target_range_string(target: &InferTarget) -> String {
    format!("{}:{}-{}", target.chrom, target.start, target.end)
}

fn candidate_region(candidate: &genotype::HaplotypeCandidate) -> String {
    format!(
        "{}:{}-{}({})",
        candidate.path_name, candidate.start, candidate.end, candidate.strand
    )
}

fn resolve_targets(config: &InferConfig<'_>) -> io::Result<Vec<InferTarget>> {
    let targets = if config.targets.is_empty() {
        let discovery = config.discovery.as_ref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                "infer requires --target-range, --target-bed, --partitions, or discovery settings",
            )
        })?;
        discover_partitions(config.syng_prefix, discovery)?
    } else {
        config.targets.clone()
    };

    if targets.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "no inference targets were found",
        ));
    }
    Ok(targets)
}

fn split_phase_blocks(targets: Vec<InferTarget>, block_size: usize) -> Vec<InferTarget> {
    if block_size == 0 {
        return targets;
    }

    let mut blocks = Vec::new();
    for target in targets {
        let len = target.end.saturating_sub(target.start) as usize;
        if len <= block_size {
            blocks.push(target);
            continue;
        }

        let mut block_start = target.start;
        let mut block_idx = 0usize;
        while block_start < target.end {
            let block_end = block_start
                .saturating_add(block_size.min(i32::MAX as usize) as i32)
                .min(target.end);
            blocks.push(InferTarget {
                partition: format!("{}#block{}", target.partition, block_idx),
                chrom: target.chrom.clone(),
                start: block_start,
                end: block_end,
                name: format!("{}#block{}", target.name, block_idx),
            });
            block_start = block_end;
            block_idx += 1;
        }
    }
    blocks
}

fn compute_local_call_sets(
    syng_index: &syng::SyngIndex,
    coverage: &pack::Coverage,
    targets: &[InferTarget],
    config: &InferConfig<'_>,
) -> Vec<LocalCallSet> {
    let mut call_sets = Vec::with_capacity(targets.len());
    for target in targets {
        let target_range = target_range_string(target);
        let query = genotype::SyngCosigtQuery {
            target_range: &target_range,
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
        match genotype::compute_syng_cosigt(syng_index, coverage, &query) {
            Ok(output) => call_sets.push(LocalCallSet {
                target: target.clone(),
                output: Some(output),
                error: None,
            }),
            Err(e) => call_sets.push(LocalCallSet {
                target: target.clone(),
                output: None,
                error: Some(e.to_string().replace('\t', " ")),
            }),
        }
    }
    call_sets
}

fn write_local_infer_output<W: Write>(
    out: &mut W,
    call_sets: &[LocalCallSet],
    config: &InferConfig<'_>,
) -> io::Result<()> {
    writeln!(out, "#impg infer")?;
    writeln!(out, "#evidence_backend\t{}", EvidenceBackend::Pack.as_str())?;
    if config.gaf_path.is_some() {
        writeln!(out, "#read_walks\tgaf")?;
        writeln!(
            out,
            "#read_link_model\tadjacent-candidate-walk-links;min_anchors={};weight={:.6}",
            config.min_read_link_anchors, config.read_link_weight
        )?;
    }
    writeln!(out, "#score\t{}", ScoringMethod::Cos.as_str())?;
    writeln!(
        out,
        "#feature_space\t{}",
        FeatureSpace::SyngSyncmerNode.as_str()
    )?;
    writeln!(out, "#targets\t{}", call_sets.len())?;
    if config.phase_block_size > 0 {
        writeln!(out, "#phase_block_size\t{}", config.phase_block_size)?;
    }
    writeln!(out, "#candidate_mode\t{:?}", config.candidate_mode)?;
    writeln!(out, "#ploidy\t{}", config.ploidy)?;
    writeln!(
        out,
        "#rank\tpartition\tchrom\tstart\tend\tmethod\tploidy\tsimilarity\tqv\thaplotypes\tregions\tcandidate_anchors\tcandidate_span_fractions\tstatus"
    )?;

    for call_set in call_sets {
        let target = &call_set.target;
        match &call_set.output {
            Some(result) => {
                for (rank, genotype_result) in result.results.iter().enumerate() {
                    let haplotypes: Vec<&str> = genotype_result
                        .combination
                        .iter()
                        .map(|&idx| result.candidates[idx].path_name.as_str())
                        .collect();
                    let regions: Vec<String> = genotype_result
                        .combination
                        .iter()
                        .map(|&idx| candidate_region(&result.candidates[idx]))
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
                        "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{:.3}\t{}\t{}\t{}\t{}\tPASS",
                        rank + 1,
                        target.partition,
                        target.chrom,
                        target.start,
                        target.end,
                        ScoringMethod::Cos.as_str(),
                        result.ploidy,
                        genotype_result.similarity,
                        genotype_result.qv,
                        haplotypes.join(","),
                        regions.join(","),
                        anchors.join(","),
                        spans.join(","),
                    )?;
                }
            }
            None => {
                writeln!(
                    out,
                    "1\t{}\t{}\t{}\t{}\t{}\t{}\t0.000000000\t0.000\t.\t.\t.\t.\tNO_CALL:{}",
                    target.partition,
                    target.chrom,
                    target.start,
                    target.end,
                    ScoringMethod::Cos.as_str(),
                    config.ploidy,
                    call_set.error.as_deref().unwrap_or("unknown error")
                )?;
            }
        }
    }

    Ok(())
}

fn unique_ordered_permutations(values: &[usize]) -> Vec<Vec<usize>> {
    fn visit(
        values: &[usize],
        used: &mut [bool],
        current: &mut Vec<usize>,
        out: &mut Vec<Vec<usize>>,
    ) {
        if current.len() == values.len() {
            if !out.contains(current) {
                out.push(current.clone());
            }
            return;
        }
        for idx in 0..values.len() {
            if used[idx] {
                continue;
            }
            used[idx] = true;
            current.push(values[idx]);
            visit(values, used, current, out);
            current.pop();
            used[idx] = false;
        }
    }

    let mut out = Vec::new();
    let mut used = vec![false; values.len()];
    let mut current = Vec::with_capacity(values.len());
    visit(values, &mut used, &mut current, &mut out);
    out
}

fn candidate_read_emission(
    read_evidence: Option<&ReadWalkEvidence>,
    key: CandidateKey,
    read_link_weight: f64,
) -> CandidateEmissionScore {
    let support = read_evidence
        .and_then(|evidence| evidence.candidate_support.get(&key))
        .cloned()
        .unwrap_or_default();
    CandidateEmissionScore {
        read_weight: support.read_weight,
        anchor_weight: support.anchor_weight,
        reward: read_link_reward(support.anchor_weight, read_link_weight),
    }
}

fn ordered_states(
    call_idx: usize,
    call_set: &LocalCallSet,
    config: &InferConfig<'_>,
    read_evidence: Option<&ReadWalkEvidence>,
) -> Vec<OrderedState> {
    let Some(output) = &call_set.output else {
        return Vec::new();
    };
    let mut states = Vec::new();
    for (result_idx, result) in output.results.iter().enumerate() {
        for candidates in unique_ordered_permutations(&result.combination) {
            let mut read_emissions = Vec::with_capacity(candidates.len());
            let mut read_emission_reward = 0.0;
            let mut rewarded_candidates = Vec::new();
            for &candidate_idx in &candidates {
                let emission = candidate_read_emission(
                    read_evidence,
                    CandidateKey {
                        call_idx,
                        candidate_idx,
                    },
                    config.read_link_weight,
                );
                if !rewarded_candidates.contains(&candidate_idx) {
                    read_emission_reward += emission.reward;
                    rewarded_candidates.push(candidate_idx);
                }
                read_emissions.push(emission);
            }
            states.push(OrderedState {
                result_idx,
                candidates,
                emission: result.qv + read_emission_reward,
                read_emissions,
            });
        }
    }
    states
}

fn start_transition() -> PhaseTransitionScore {
    PhaseTransitionScore {
        kind: TransitionKind::Start,
        cost: 0.0,
        read_link_reads: 0.0,
        read_link_anchors: 0.0,
        read_link_reward: 0.0,
    }
}

fn parse_gaf_path_nodes(path: &str, nodes: &mut Vec<i32>) -> io::Result<()> {
    nodes.clear();
    let bytes = path.as_bytes();
    let mut i = 0usize;
    while i < bytes.len() {
        let orientation = bytes[i];
        if orientation != b'>' && orientation != b'<' {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GAF path contains non-orientation byte at offset {i}: {path}"),
            ));
        }
        i += 1;
        let start = i;
        while i < bytes.len() && bytes[i].is_ascii_digit() {
            i += 1;
        }
        if i == start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GAF path step is missing numeric syncmer node: {path}"),
            ));
        }
        let node = path[start..i].parse::<i32>().map_err(|e| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("invalid syncmer node in GAF path '{path}': {e}"),
            )
        })?;
        if node <= 0 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GAF path step must be a positive syncmer node: {path}"),
            ));
        }
        nodes.push(if orientation == b'<' { -node } else { node });
    }
    Ok(())
}

fn parse_gaf_query_positions(fields: &[&str], anchors: usize) -> io::Result<Option<Vec<u64>>> {
    for field in fields.iter().skip(12) {
        if let Some(values) = field.strip_prefix("qp:B:I") {
            let mut positions = Vec::with_capacity(anchors);
            let values = values.strip_prefix(',').unwrap_or(values);
            if !values.is_empty() {
                for value in values.split(',') {
                    positions.push(value.parse::<u64>().map_err(|e| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("invalid qp:B:I query position in GAF tag '{field}': {e}"),
                        )
                    })?);
                }
            }
            if positions.len() != anchors {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GAF qp:B:I tag has {} positions but path has {} syncmer nodes",
                        positions.len(),
                        anchors
                    ),
                ));
            }
            return Ok(Some(positions));
        }
    }
    Ok(None)
}

fn build_read_walk_steps(
    nodes: &[i32],
    positions: &[u64],
    out: &mut Vec<syng::SyngWalkStep>,
) {
    out.clear();
    out.extend(
        nodes
            .iter()
            .zip(positions.iter())
            .map(|(&signed_node, &bp_pos)| syng::SyngWalkStep {
                signed_node,
                bp_pos,
            }),
    );
    out.sort_unstable_by(|a, b| a.bp_pos.cmp(&b.bp_pos).then(a.signed_node.cmp(&b.signed_node)));
}

fn reverse_read_walk_steps(
    steps: &[syng::SyngWalkStep],
    query_len: u64,
    syncmer_len: u64,
    out: &mut Vec<syng::SyngWalkStep>,
) {
    out.clear();
    out.extend(steps.iter().rev().map(|step| syng::SyngWalkStep {
        signed_node: -step.signed_node,
        bp_pos: query_len.saturating_sub(step.bp_pos.saturating_add(syncmer_len)),
    }));
    out.sort_unstable_by(|a, b| a.bp_pos.cmp(&b.bp_pos).then(a.signed_node.cmp(&b.signed_node)));
}

fn sorted_candidate_hits(
    counts: FxHashMap<CandidateKey, u32>,
    min_read_link_anchors: usize,
) -> Vec<(usize, Vec<(usize, u32)>)> {
    let mut by_call: BTreeMap<usize, Vec<(usize, u32)>> = BTreeMap::new();
    let min_read_link_anchors = min_read_link_anchors.max(1) as u32;
    for (key, count) in counts {
        if count >= min_read_link_anchors {
            by_call
                .entry(key.call_idx)
                .or_default()
                .push((key.candidate_idx, count));
        }
    }
    let mut hits: Vec<(usize, Vec<(usize, u32)>)> = by_call.into_iter().collect();
    for (_, candidates) in &mut hits {
        candidates.sort_unstable_by(|a, b| b.1.cmp(&a.1).then(a.0.cmp(&b.0)));
    }
    hits
}

fn add_read_links(evidence: &mut ReadWalkEvidence, hits_by_call: &[(usize, Vec<(usize, u32)>)]) {
    if hits_by_call.len() < 2 {
        return;
    }

    for adjacent in hits_by_call.windows(2) {
        let (prev_call_idx, prev_hits) = &adjacent[0];
        let (curr_call_idx, curr_hits) = &adjacent[1];
        if prev_hits.is_empty() || curr_hits.is_empty() {
            continue;
        }
        let denominator = (prev_hits.len() * curr_hits.len()) as f64;
        for &(prev_candidate_idx, prev_count) in prev_hits {
            for &(curr_candidate_idx, curr_count) in curr_hits {
                let prev_key = CandidateKey {
                    call_idx: *prev_call_idx,
                    candidate_idx: prev_candidate_idx,
                };
                let curr_key = CandidateKey {
                    call_idx: *curr_call_idx,
                    candidate_idx: curr_candidate_idx,
                };
                let entry = evidence.links.entry((prev_key, curr_key)).or_default();
                entry.read_weight += 1.0 / denominator;
                entry.anchor_weight += f64::from(prev_count.min(curr_count)) / denominator;
                evidence.linked_read_pairs += 1;
            }
        }
    }
}

fn add_candidate_support(
    evidence: &mut ReadWalkEvidence,
    hits_by_call: &[(usize, Vec<(usize, u32)>)],
) {
    if hits_by_call.is_empty() {
        return;
    }
    evidence.records_with_candidate_hits += 1;
    evidence.candidate_hits += hits_by_call
        .iter()
        .map(|(_, hits)| hits.len() as u64)
        .sum::<u64>();

    for (call_idx, hits) in hits_by_call {
        if hits.is_empty() {
            continue;
        }
        let denominator = hits.len() as f64;
        for &(candidate_idx, count) in hits {
            let key = CandidateKey {
                call_idx: *call_idx,
                candidate_idx,
            };
            let entry = evidence.candidate_support.entry(key).or_default();
            entry.read_weight += 1.0 / denominator;
            entry.anchor_weight += f64::from(count) / denominator;
        }
    }
}

fn lis_len(values: &[u32]) -> u32 {
    let mut tails: Vec<u32> = Vec::new();
    for &value in values {
        let idx = tails.partition_point(|&tail| tail < value);
        if idx == tails.len() {
            tails.push(value);
        } else {
            tails[idx] = value;
        }
    }
    tails.len() as u32
}

fn add_whole_walk_orientation_hits(
    counts: &mut FxHashMap<CandidateKey, u32>,
    walk_index: &FxHashMap<i32, Vec<WalkOccurrence>>,
    nodes: &[i32],
) {
    let mut positions_by_candidate: FxHashMap<CandidateKey, Vec<u32>> = FxHashMap::default();
    for &node in nodes {
        if let Some(occurrences) = walk_index.get(&node) {
            for &occurrence in occurrences {
                positions_by_candidate
                    .entry(occurrence.key)
                    .or_default()
                    .push(occurrence.position);
            }
        }
    }

    for (key, positions) in positions_by_candidate {
        let matched = lis_len(&positions);
        if matched == 0 {
            continue;
        }
        let entry = counts.entry(key).or_insert(0);
        *entry = (*entry).max(matched);
    }
}

fn add_whole_walk_hits(
    counts: &mut FxHashMap<CandidateKey, u32>,
    walk_index: &FxHashMap<i32, Vec<WalkOccurrence>>,
    nodes: &[i32],
    reverse_nodes: &mut Vec<i32>,
) {
    if nodes.is_empty() {
        return;
    }
    add_whole_walk_orientation_hits(counts, walk_index, nodes);
    reverse_nodes.clear();
    reverse_nodes.extend(nodes.iter().rev().map(|node| -*node));
    add_whole_walk_orientation_hits(counts, walk_index, reverse_nodes);
}

fn candidate_mem_overlap_len(
    read_steps: &[syng::SyngWalkStep],
    mem_start: usize,
    mem_end: usize,
    read_start: usize,
    candidate_walk: &[syng::SyngWalkStep],
    candidate_start: usize,
) -> u32 {
    if read_start >= mem_end
        || candidate_start >= candidate_walk.len()
        || read_steps[read_start].signed_node != candidate_walk[candidate_start].signed_node
    {
        return 0;
    }

    let mut read_left = read_start;
    let mut candidate_left = candidate_start;
    while read_left > mem_start && candidate_left > 0 {
        let prev_read = read_left - 1;
        let prev_candidate = candidate_left - 1;
        if read_steps[prev_read].signed_node != candidate_walk[prev_candidate].signed_node {
            break;
        }
        let Some(read_offset) = read_steps[read_left]
            .bp_pos
            .checked_sub(read_steps[prev_read].bp_pos)
        else {
            break;
        };
        let Some(candidate_offset) = candidate_walk[candidate_left]
            .bp_pos
            .checked_sub(candidate_walk[prev_candidate].bp_pos)
        else {
            break;
        };
        if read_offset != candidate_offset {
            break;
        }
        read_left = prev_read;
        candidate_left = prev_candidate;
    }

    let mut read_right = read_start;
    let mut candidate_right = candidate_start;
    while read_right + 1 < mem_end && candidate_right + 1 < candidate_walk.len() {
        let next_read = read_right + 1;
        let next_candidate = candidate_right + 1;
        if read_steps[next_read].signed_node != candidate_walk[next_candidate].signed_node {
            break;
        }
        let Some(read_offset) = read_steps[next_read]
            .bp_pos
            .checked_sub(read_steps[read_right].bp_pos)
        else {
            break;
        };
        let Some(candidate_offset) = candidate_walk[next_candidate]
            .bp_pos
            .checked_sub(candidate_walk[candidate_right].bp_pos)
        else {
            break;
        };
        if read_offset != candidate_offset {
            break;
        }
        read_right = next_read;
        candidate_right = next_candidate;
    }

    (read_right - read_left + 1) as u32
}

fn add_mem_hits(
    counts: &mut FxHashMap<CandidateKey, u32>,
    call_sets: &[LocalCallSet],
    walk_index: &FxHashMap<i32, Vec<WalkOccurrence>>,
    read_steps: &[syng::SyngWalkStep],
    mems: &[syng::SyngGbwtMem],
) {
    for mem in mems {
        let mut best_by_candidate: FxHashMap<CandidateKey, u32> = FxHashMap::default();
        for read_idx in mem.step_start..mem.step_end {
            let Some(occurrences) = walk_index.get(&read_steps[read_idx].signed_node) else {
                continue;
            };
            for &occurrence in occurrences {
                let Some(output) = &call_sets[occurrence.key.call_idx].output else {
                    continue;
                };
                let candidate = &output.candidates[occurrence.key.candidate_idx];
                let overlap = candidate_mem_overlap_len(
                    read_steps,
                    mem.step_start,
                    mem.step_end,
                    read_idx,
                    &candidate.oriented_walk,
                    occurrence.position as usize,
                );
                let entry = best_by_candidate.entry(occurrence.key).or_insert(0);
                *entry = (*entry).max(overlap);
            }
        }
        for (key, anchors) in best_by_candidate {
            if anchors > 0 {
                *counts.entry(key).or_insert(0) += anchors;
            }
        }
    }
}

fn build_read_walk_evidence(
    syng_index: &syng::SyngIndex,
    call_sets: &[LocalCallSet],
    gaf_path: &str,
    min_read_link_anchors: usize,
) -> io::Result<ReadWalkEvidence> {
    let mut walk_index: FxHashMap<i32, Vec<WalkOccurrence>> = FxHashMap::default();
    for (call_idx, call_set) in call_sets.iter().enumerate() {
        let Some(output) = &call_set.output else {
            continue;
        };
        for (candidate_idx, candidate) in output.candidates.iter().enumerate() {
            let key = CandidateKey {
                call_idx,
                candidate_idx,
            };
            for (position, &node) in candidate.oriented_walk.iter().enumerate() {
                walk_index
                    .entry(node.signed_node)
                    .or_default()
                    .push(WalkOccurrence {
                        key,
                        position: position as u32,
                    });
            }
        }
    }
    for occurrences in walk_index.values_mut() {
        occurrences.sort_unstable_by(|a, b| {
            a.key
                .call_idx
                .cmp(&b.key.call_idx)
                .then(a.key.candidate_idx.cmp(&b.key.candidate_idx))
                .then(b.position.cmp(&a.position))
        });
    }

    if walk_index.is_empty() {
        return Ok(ReadWalkEvidence::default());
    }

    let file = File::open(gaf_path)?;
    let (reader, _) = niffler::get_reader(Box::new(file)).map_err(|e| {
        io::Error::other(format!("failed to open read-walk GAF '{}': {e}", gaf_path))
    })?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, reader);
    let mut evidence = ReadWalkEvidence::default();
    let mut nodes = Vec::new();
    let mut reverse_nodes = Vec::new();
    let mut read_steps = Vec::new();
    let mut reverse_steps = Vec::new();
    let syncmer_len = syng_index.syncmer_length_bp() as u64;

    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        evidence.records += 1;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 6 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GAF line {} has fewer than 6 fields", line_no + 1),
            ));
        }
        parse_gaf_path_nodes(fields[5], &mut nodes)?;
        let mut walk_counts: FxHashMap<CandidateKey, u32> = FxHashMap::default();
        if let Some(positions) = parse_gaf_query_positions(&fields, nodes.len())? {
            let query_len = fields[1].parse::<u64>().map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid GAF query length on line {}: {e}", line_no + 1),
                )
            })?;
            build_read_walk_steps(&nodes, &positions, &mut read_steps);
            let mems = syng_index.gbwt_mems_for_walk(&read_steps)?;
            add_mem_hits(&mut walk_counts, call_sets, &walk_index, &read_steps, &mems);
            reverse_read_walk_steps(&read_steps, query_len, syncmer_len, &mut reverse_steps);
            let reverse_mems = syng_index.gbwt_mems_for_walk(&reverse_steps)?;
            add_mem_hits(
                &mut walk_counts,
                call_sets,
                &walk_index,
                &reverse_steps,
                &reverse_mems,
            );
        } else {
            add_whole_walk_hits(&mut walk_counts, &walk_index, &nodes, &mut reverse_nodes);
        }
        let walk_hits_by_call = sorted_candidate_hits(walk_counts, min_read_link_anchors);
        add_candidate_support(&mut evidence, &walk_hits_by_call);
        add_read_links(&mut evidence, &walk_hits_by_call);
    }

    Ok(evidence)
}

fn read_link_reward(anchor_weight: f64, read_link_weight: f64) -> f64 {
    if anchor_weight <= 0.0 || read_link_weight <= 0.0 {
        0.0
    } else {
        read_link_weight * 10.0 * (1.0 + anchor_weight).log10()
    }
}

fn phase_transition(
    prev: &genotype::HaplotypeCandidate,
    curr: &genotype::HaplotypeCandidate,
    switch_penalty: f64,
    stitch_gap: u64,
) -> (TransitionKind, f64) {
    if prev.path_name == curr.path_name && prev.strand == curr.strand {
        let gap = if curr.start >= prev.end {
            curr.start - prev.end
        } else {
            prev.start.saturating_sub(curr.end)
        };
        if gap <= stitch_gap || curr.start <= prev.end {
            return (TransitionKind::Continue, 0.0);
        }
        return (TransitionKind::Uncertain, switch_penalty * 0.5);
    }
    (TransitionKind::Recombine, switch_penalty)
}

fn transition_between(
    prev_call: &LocalCallSet,
    prev_call_idx: usize,
    prev_state: &OrderedState,
    curr_call: &LocalCallSet,
    curr_call_idx: usize,
    curr_state: &OrderedState,
    config: &InferConfig<'_>,
    read_evidence: Option<&ReadWalkEvidence>,
) -> (Vec<PhaseTransitionScore>, f64, f64) {
    let Some(prev_output) = &prev_call.output else {
        return (Vec::new(), 0.0, 0.0);
    };
    let Some(curr_output) = &curr_call.output else {
        return (Vec::new(), 0.0, 0.0);
    };
    let mut transitions = Vec::with_capacity(prev_state.candidates.len());
    let mut total_cost = 0.0;
    let mut total_reward = 0.0;
    for (&prev_idx, &curr_idx) in prev_state
        .candidates
        .iter()
        .zip(curr_state.candidates.iter())
    {
        let (kind, phase_cost) = phase_transition(
            &prev_output.candidates[prev_idx],
            &curr_output.candidates[curr_idx],
            config.switch_penalty,
            config.stitch_gap,
        );
        let link = read_evidence
            .and_then(|evidence| {
                evidence.links.get(&(
                    CandidateKey {
                        call_idx: prev_call_idx,
                        candidate_idx: prev_idx,
                    },
                    CandidateKey {
                        call_idx: curr_call_idx,
                        candidate_idx: curr_idx,
                    },
                ))
            })
            .cloned()
            .unwrap_or_default();
        let reward = read_link_reward(link.anchor_weight, config.read_link_weight);
        total_cost += phase_cost;
        total_reward += reward;
        transitions.push(PhaseTransitionScore {
            kind,
            cost: phase_cost,
            read_link_reads: link.read_weight,
            read_link_anchors: link.anchor_weight,
            read_link_reward: reward,
        });
    }
    (transitions, total_cost, total_reward)
}

fn stitch_mosaic(
    call_sets: &[LocalCallSet],
    config: &InferConfig<'_>,
    read_evidence: Option<&ReadWalkEvidence>,
) -> io::Result<Option<MosaicPath>> {
    let mut indexed_states: Vec<(usize, Vec<OrderedState>)> = Vec::new();
    for (idx, call_set) in call_sets.iter().enumerate() {
        let states = ordered_states(idx, call_set, config, read_evidence);
        if !states.is_empty() {
            indexed_states.push((idx, states));
        }
    }
    if indexed_states.is_empty() {
        return Ok(None);
    }

    let (first_idx, first_states) = &indexed_states[0];
    let mut beam: Vec<MosaicPath> = first_states
        .iter()
        .map(|state| MosaicPath {
            score: state.emission,
            steps: vec![MosaicStep {
                call_idx: *first_idx,
                state: state.clone(),
                transitions: vec![start_transition(); state.candidates.len()],
            }],
        })
        .collect();
    beam.sort_by(|a, b| b.score.total_cmp(&a.score));
    beam.truncate(config.stitch_beam);

    for (call_idx, states) in indexed_states.iter().skip(1) {
        let mut next = Vec::new();
        for path in &beam {
            let prev_step = path.steps.last().expect("beam path has at least one step");
            let prev_call = &call_sets[prev_step.call_idx];
            for state in states {
                let (transitions, transition_cost, read_link_reward) = transition_between(
                    prev_call,
                    prev_step.call_idx,
                    &prev_step.state,
                    &call_sets[*call_idx],
                    *call_idx,
                    state,
                    config,
                    read_evidence,
                );
                let mut candidate = path.clone();
                candidate.score += state.emission - transition_cost + read_link_reward;
                candidate.steps.push(MosaicStep {
                    call_idx: *call_idx,
                    state: state.clone(),
                    transitions,
                });
                next.push(candidate);
            }
        }
        next.sort_by(|a, b| b.score.total_cmp(&a.score));
        next.truncate(config.stitch_beam);
        beam = next;
    }

    Ok(beam.into_iter().max_by(|a, b| a.score.total_cmp(&b.score)))
}

fn transition_for_phase(step: &MosaicStep, phase: usize) -> PhaseTransitionScore {
    step.transitions
        .get(phase)
        .cloned()
        .unwrap_or_else(|| PhaseTransitionScore {
            kind: TransitionKind::Uncertain,
            cost: 0.0,
            read_link_reads: 0.0,
            read_link_anchors: 0.0,
            read_link_reward: 0.0,
        })
}

fn read_emission_for_phase(step: &MosaicStep, phase: usize) -> CandidateEmissionScore {
    step.state
        .read_emissions
        .get(phase)
        .cloned()
        .unwrap_or_default()
}

fn write_mosaic_tsv(path: &str, mosaic: &MosaicPath, call_sets: &[LocalCallSet]) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "#impg infer mosaic")?;
    writeln!(out, "#score\t{:.6}", mosaic.score)?;
    writeln!(
        out,
        "#phase\tpartition\tchrom\tstart\tend\tpath\tpath_start\tpath_end\tstrand\tlocal_rank\tlocal_similarity\tlocal_qv\ttransition_kind\ttransition_cost\tread_link_reads\tread_link_anchors\tread_link_reward\tread_emission_reads\tread_emission_anchors\tread_emission_reward"
    )?;
    for step in &mosaic.steps {
        let call_set = &call_sets[step.call_idx];
        let output = call_set.output.as_ref().expect("mosaic step has output");
        let result = &output.results[step.state.result_idx];
        for (phase, &candidate_idx) in step.state.candidates.iter().enumerate() {
            let candidate = &output.candidates[candidate_idx];
            let transition = transition_for_phase(step, phase);
            let read_emission = read_emission_for_phase(step, phase);
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{:.3}\t{}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}\t{:.3}",
                phase,
                call_set.target.partition,
                call_set.target.chrom,
                call_set.target.start,
                call_set.target.end,
                candidate.path_name,
                candidate.start,
                candidate.end,
                candidate.strand,
                step.state.result_idx + 1,
                result.similarity,
                result.qv,
                transition.kind.as_str(),
                transition.cost,
                transition.read_link_reads,
                transition.read_link_anchors,
                transition.read_link_reward,
                read_emission.read_weight,
                read_emission.anchor_weight,
                read_emission.reward,
            )?;
        }
    }
    out.flush()
}

fn fetch_candidate_sequence(
    sequence_index: &UnifiedSequenceIndex,
    candidate: &genotype::HaplotypeCandidate,
) -> io::Result<Vec<u8>> {
    let mut seq = sequence_index.fetch_sequence(
        &candidate.path_name,
        candidate.start as i32,
        candidate.end as i32,
    )?;
    if candidate.strand == '-' {
        seq = reverse_complement(&seq);
    }
    Ok(seq)
}

fn write_wrapped_fasta<W: Write>(out: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writeln!(out, ">{name}")?;
    for chunk in seq.chunks(80) {
        out.write_all(chunk)?;
        out.write_all(b"\n")?;
    }
    Ok(())
}

fn write_mosaic_fasta(
    path: &str,
    mosaic: &MosaicPath,
    call_sets: &[LocalCallSet],
    sequence_index: &UnifiedSequenceIndex,
    strict: bool,
) -> io::Result<()> {
    let ploidy = mosaic
        .steps
        .first()
        .map(|step| step.state.candidates.len())
        .unwrap_or(0);
    let mut phase_sequences = vec![Vec::<u8>::new(); ploidy];
    for step in &mosaic.steps {
        let call_set = &call_sets[step.call_idx];
        let output = call_set.output.as_ref().expect("mosaic step has output");
        for (phase, &candidate_idx) in step.state.candidates.iter().enumerate() {
            let kind = transition_for_phase(step, phase).kind;
            if matches!(kind, TransitionKind::Recombine | TransitionKind::Uncertain)
                && !phase_sequences[phase].is_empty()
            {
                if strict {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "strict stitch rejected {} join in partition {}",
                            kind.as_str(),
                            call_set.target.partition
                        ),
                    ));
                }
                phase_sequences[phase].extend_from_slice(b"NNNNNNNNNN");
            }
            let seq = fetch_candidate_sequence(sequence_index, &output.candidates[candidate_idx])?;
            phase_sequences[phase].extend_from_slice(&seq);
        }
    }
    let mut out = BufWriter::new(File::create(path)?);
    for (phase, seq) in phase_sequences.iter().enumerate() {
        write_wrapped_fasta(&mut out, &format!("impg_phase_{phase}"), seq)?;
    }
    out.flush()
}

fn write_mosaic_gfa(
    path: &str,
    mosaic: &MosaicPath,
    call_sets: &[LocalCallSet],
    sequence_index: &UnifiedSequenceIndex,
    strict: bool,
) -> io::Result<()> {
    let ploidy = mosaic
        .steps
        .first()
        .map(|step| step.state.candidates.len())
        .unwrap_or(0);
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "H\tVN:Z:1.0")?;
    let mut phase_segments: Vec<Vec<String>> = vec![Vec::new(); ploidy];
    for (step_idx, step) in mosaic.steps.iter().enumerate() {
        let call_set = &call_sets[step.call_idx];
        let output = call_set.output.as_ref().expect("mosaic step has output");
        for (phase, &candidate_idx) in step.state.candidates.iter().enumerate() {
            let kind = transition_for_phase(step, phase).kind;
            if strict
                && matches!(kind, TransitionKind::Recombine | TransitionKind::Uncertain)
                && !phase_segments[phase].is_empty()
            {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "strict stitch rejected {} join in partition {}",
                        kind.as_str(),
                        call_set.target.partition
                    ),
                ));
            }
            let candidate = &output.candidates[candidate_idx];
            let segment = format!("phase{phase}_{step_idx}");
            let seq = fetch_candidate_sequence(sequence_index, candidate)?;
            writeln!(out, "S\t{}\t{}", segment, String::from_utf8_lossy(&seq))?;
            if let Some(prev) = phase_segments[phase].last() {
                writeln!(
                    out,
                    "L\t{}\t+\t{}\t+\t0M\tTK:Z:{}",
                    prev,
                    segment,
                    kind.as_str()
                )?;
            }
            phase_segments[phase].push(segment);
        }
    }
    for (phase, segments) in phase_segments.iter().enumerate() {
        let walk = segments
            .iter()
            .map(|seg| format!("{seg}+"))
            .collect::<Vec<_>>()
            .join(",");
        writeln!(out, "P\timpg_phase_{phase}\t{}\t*", walk)?;
    }
    out.flush()
}

pub fn run_syng_pack_infer<W: Write>(out: &mut W, config: &InferConfig<'_>) -> io::Result<()> {
    let targets = split_phase_blocks(resolve_targets(config)?, config.phase_block_size);

    info!("Loading syng index from prefix: {}", config.syng_prefix);
    let syng_index = syng::SyngIndex::load(config.syng_prefix, syng::SyncmerParams::default())?;
    info!("Loading sample pack coverage from {}", config.pack_path);
    let coverage = pack::read(config.pack_path)?;
    if let Some(gaf_path) = config.gaf_path {
        info!("Using read-walk projection from {}", gaf_path);
    }

    let call_sets = compute_local_call_sets(&syng_index, &coverage, &targets, config);
    write_local_infer_output(out, &call_sets, config)?;

    let wants_stitch = config.stitch == StitchMode::Beam
        || config.emit_mosaic.is_some()
        || config.emit_fasta.is_some()
        || config.emit_gfa.is_some();
    if wants_stitch {
        let read_evidence = if let Some(gaf_path) = config.gaf_path {
            let evidence = build_read_walk_evidence(
                &syng_index,
                &call_sets,
                gaf_path,
                config.min_read_link_anchors,
            )?;
            info!(
                "Loaded read-walk evidence from {}: {} GAF records, {} records with candidate hits, {} candidate hits, {} local candidate scores, {} linked candidate pairs, {} scored transitions",
                gaf_path,
                evidence.records,
                evidence.records_with_candidate_hits,
                evidence.candidate_hits,
                evidence.candidate_support.len(),
                evidence.linked_read_pairs,
                evidence.links.len()
            );
            Some(evidence)
        } else {
            None
        };
        let mosaic =
            stitch_mosaic(&call_sets, config, read_evidence.as_ref())?.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "no PASS local calls available to stitch into a mosaic",
                )
            })?;
        if let Some(path) = config.emit_mosaic {
            write_mosaic_tsv(path, &mosaic, &call_sets)?;
        }
        if let Some(path) = config.emit_fasta {
            let sequence_index = config.sequence_index.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--emit-fasta requires --sequence-files or --sequence-list",
                )
            })?;
            write_mosaic_fasta(
                path,
                &mosaic,
                &call_sets,
                sequence_index,
                config.strict_stitch,
            )?;
        }
        if let Some(path) = config.emit_gfa {
            let sequence_index = config.sequence_index.ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    "--emit-gfa requires --sequence-files or --sequence-list",
                )
            })?;
            write_mosaic_gfa(
                path,
                &mosaic,
                &call_sets,
                sequence_index,
                config.strict_stitch,
            )?;
        }
    }

    Ok(())
}
