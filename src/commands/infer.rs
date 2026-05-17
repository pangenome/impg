use clap::ValueEnum;
use log::info;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::commands::{genotype, graph, partition};
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
    pub stitch_beam: usize,
    pub switch_penalty: f64,
    pub stitch_gap: u64,
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
    transition_kinds: Vec<TransitionKind>,
    transition_cost: f64,
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
        pipeline: graph::GraphBuildConfig::default(),
        target_poa_lengths: vec![700, 1100],
        max_node_length: 100,
        poa_padding_fraction: 0.001,
        partition_size: None,
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
    writeln!(out, "#evidence_backend\tpack")?;
    if config.gaf_path.is_some() {
        writeln!(out, "#read_walks\tgaf")?;
    }
    writeln!(out, "#score\tcos")?;
    writeln!(out, "#feature_space\tsyng-syncmer-node")?;
    writeln!(out, "#targets\t{}", call_sets.len())?;
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
                        "{}\t{}\t{}\t{}\t{}\tcos\t{}\t{:.9}\t{:.3}\t{}\t{}\t{}\t{}\tPASS",
                        rank + 1,
                        target.partition,
                        target.chrom,
                        target.start,
                        target.end,
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
                    "1\t{}\t{}\t{}\t{}\tcos\t{}\t0.000000000\t0.000\t.\t.\t.\t.\tNO_CALL:{}",
                    target.partition,
                    target.chrom,
                    target.start,
                    target.end,
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

fn ordered_states(call_set: &LocalCallSet) -> Vec<OrderedState> {
    let Some(output) = &call_set.output else {
        return Vec::new();
    };
    let mut states = Vec::new();
    for (result_idx, result) in output.results.iter().enumerate() {
        for candidates in unique_ordered_permutations(&result.combination) {
            states.push(OrderedState {
                result_idx,
                candidates,
                emission: result.qv,
            });
        }
    }
    states
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
    prev_state: &OrderedState,
    curr_call: &LocalCallSet,
    curr_state: &OrderedState,
    config: &InferConfig<'_>,
) -> (Vec<TransitionKind>, f64) {
    let Some(prev_output) = &prev_call.output else {
        return (Vec::new(), 0.0);
    };
    let Some(curr_output) = &curr_call.output else {
        return (Vec::new(), 0.0);
    };
    let mut kinds = Vec::with_capacity(prev_state.candidates.len());
    let mut cost = 0.0;
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
        kinds.push(kind);
        cost += phase_cost;
    }
    (kinds, cost)
}

fn stitch_mosaic(
    call_sets: &[LocalCallSet],
    config: &InferConfig<'_>,
) -> io::Result<Option<MosaicPath>> {
    let indexed_states: Vec<(usize, Vec<OrderedState>)> = call_sets
        .iter()
        .enumerate()
        .filter_map(|(idx, call_set)| {
            let states = ordered_states(call_set);
            if states.is_empty() {
                None
            } else {
                Some((idx, states))
            }
        })
        .collect();
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
                transition_kinds: vec![TransitionKind::Start; state.candidates.len()],
                transition_cost: 0.0,
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
                let (kinds, transition_cost) = transition_between(
                    prev_call,
                    &prev_step.state,
                    &call_sets[*call_idx],
                    state,
                    config,
                );
                let mut candidate = path.clone();
                candidate.score += state.emission - transition_cost;
                candidate.steps.push(MosaicStep {
                    call_idx: *call_idx,
                    state: state.clone(),
                    transition_kinds: kinds,
                    transition_cost,
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

fn write_mosaic_tsv(path: &str, mosaic: &MosaicPath, call_sets: &[LocalCallSet]) -> io::Result<()> {
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "#impg infer mosaic")?;
    writeln!(out, "#score\t{:.6}", mosaic.score)?;
    writeln!(
        out,
        "#phase\tpartition\tchrom\tstart\tend\tpath\tpath_start\tpath_end\tstrand\tlocal_rank\tlocal_similarity\tlocal_qv\ttransition_kind\ttransition_cost"
    )?;
    for step in &mosaic.steps {
        let call_set = &call_sets[step.call_idx];
        let output = call_set.output.as_ref().expect("mosaic step has output");
        let result = &output.results[step.state.result_idx];
        for (phase, &candidate_idx) in step.state.candidates.iter().enumerate() {
            let candidate = &output.candidates[candidate_idx];
            writeln!(
                out,
                "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.9}\t{:.3}\t{}\t{:.3}",
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
                step.transition_kinds
                    .get(phase)
                    .copied()
                    .unwrap_or(TransitionKind::Uncertain)
                    .as_str(),
                step.transition_cost,
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
            let kind = step
                .transition_kinds
                .get(phase)
                .copied()
                .unwrap_or(TransitionKind::Uncertain);
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
            let kind = step
                .transition_kinds
                .get(phase)
                .copied()
                .unwrap_or(TransitionKind::Uncertain);
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
    let targets = resolve_targets(config)?;

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
        let mosaic = stitch_mosaic(&call_sets, config)?.ok_or_else(|| {
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
