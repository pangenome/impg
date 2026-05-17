use log::info;
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

use crate::commands::{genotype, graph, partition};
use crate::pack;
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

pub fn run_syng_pack_infer<W: Write>(out: &mut W, config: &InferConfig<'_>) -> io::Result<()> {
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

    info!("Loading syng index from prefix: {}", config.syng_prefix);
    let syng_index = syng::SyngIndex::load(config.syng_prefix, syng::SyncmerParams::default())?;
    info!("Loading sample pack coverage from {}", config.pack_path);
    let coverage = pack::read(config.pack_path)?;

    writeln!(out, "#impg infer")?;
    writeln!(out, "#evidence_backend\tpack")?;
    writeln!(out, "#score\tcos")?;
    writeln!(out, "#feature_space\tsyng-syncmer-node")?;
    writeln!(out, "#targets\t{}", targets.len())?;
    writeln!(out, "#candidate_mode\t{:?}", config.candidate_mode)?;
    writeln!(out, "#ploidy\t{}", config.ploidy)?;
    writeln!(
        out,
        "#rank\tpartition\tchrom\tstart\tend\tmethod\tploidy\tsimilarity\tqv\thaplotypes\tregions\tcandidate_anchors\tcandidate_span_fractions\tstatus"
    )?;

    for target in &targets {
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

        match genotype::compute_syng_cosigt(&syng_index, &coverage, &query) {
            Ok(result) => {
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
            Err(e) => {
                writeln!(
                    out,
                    "1\t{}\t{}\t{}\t{}\tcos\t{}\t0.000000000\t0.000\t.\t.\t.\t.\tNO_CALL:{}",
                    target.partition,
                    target.chrom,
                    target.start,
                    target.end,
                    config.ploidy,
                    e.to_string().replace('\t', " ")
                )?;
            }
        }
    }

    Ok(())
}
