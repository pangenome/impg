//! Explicit local sequence-set to seed-graph induction.
//!
//! This module is the adapter boundary for localized graph construction:
//! callers collect exact local haplotype sequences once, then choose a seed
//! graph route over that same named sequence set.  The implementation keeps
//! user-visible path names as the aligner/seqwish sequence IDs so no synthetic
//! local namespace can leak into output paths.

use crate::commands;
use crate::graph;
use crate::impg_index;
use crate::resolution;
use crate::sequence_index;
use crate::syng;
use crate::syng_graph;
use coitrees::Interval;
use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::io::{self, Read};
use std::path::Path;
use std::time::Instant;

#[derive(Clone, Debug)]
pub struct LocalSequenceRecord {
    pub path_name: String,
    pub source_name: String,
    pub source_start: u64,
    pub source_end: u64,
    pub strand: char,
    pub sequence: Vec<u8>,
}

impl LocalSequenceRecord {
    pub fn len(&self) -> usize {
        self.sequence.len()
    }

    pub fn is_empty(&self) -> bool {
        self.sequence.is_empty()
    }
}

#[derive(Clone, Debug, Default)]
pub struct LocalSequenceSet {
    records: Vec<LocalSequenceRecord>,
}

impl LocalSequenceSet {
    pub fn new(records: Vec<LocalSequenceRecord>) -> Self {
        Self { records }
    }

    pub fn collect_from_intervals(
        impg: &impl impg_index::ImpgIndex,
        query_intervals: &[Interval<u32>],
        sequence_index: &sequence_index::UnifiedSequenceIndex,
    ) -> io::Result<Self> {
        let sequences = graph::prepare_sequences(impg, query_intervals, sequence_index)?;
        let records = sequences
            .into_iter()
            .map(|(seq, meta)| {
                let (source_start, source_end) = if meta.strand == '+' {
                    (meta.start as u64, (meta.start + meta.size) as u64)
                } else {
                    (
                        (meta.total_length as i32 - meta.start - meta.size) as u64,
                        (meta.total_length as i32 - meta.start) as u64,
                    )
                };
                LocalSequenceRecord {
                    path_name: meta.path_name(),
                    source_name: meta.name,
                    source_start,
                    source_end,
                    strand: meta.strand,
                    sequence: seq.into_bytes(),
                }
            })
            .collect();
        Ok(Self { records })
    }

    pub fn records(&self) -> &[LocalSequenceRecord] {
        &self.records
    }

    pub fn is_empty(&self) -> bool {
        self.records.is_empty()
    }

    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn total_bp(&self) -> usize {
        self.records.iter().map(LocalSequenceRecord::len).sum()
    }

    pub fn named_sequences(&self) -> Vec<(String, Vec<u8>)> {
        self.records
            .iter()
            .map(|record| (record.path_name.clone(), record.sequence.clone()))
            .collect()
    }

    pub fn sweepga_named_refs(&self) -> Vec<(String, &[u8])> {
        self.records
            .iter()
            .map(|record| (record.path_name.clone(), record.sequence.as_slice()))
            .collect()
    }

    pub fn syng_anchor_members(&self) -> Vec<(String, Vec<u8>, syng_graph::Member)> {
        self.records
            .iter()
            .map(|record| {
                (
                    record.path_name.clone(),
                    record.sequence.clone(),
                    syng_graph::Member {
                        chrom: record.source_name.clone(),
                        fwd_start: record.source_start,
                        fwd_end: record.source_end,
                        strand: record.strand,
                    },
                )
            })
            .collect()
    }

    fn expected_path_sequences(&self) -> BTreeMap<String, String> {
        self.records
            .iter()
            .map(|record| {
                (
                    record.path_name.clone(),
                    String::from_utf8_lossy(&record.sequence).into_owned(),
                )
            })
            .collect()
    }
}

#[derive(Clone, Debug)]
pub struct SyngAnchorSeedConfig {
    pub syng_padding: u64,
    pub max_depth: u16,
    pub k_near: usize,
    pub k_far: usize,
    pub random_fraction: f64,
}

impl Default for SyngAnchorSeedConfig {
    fn default() -> Self {
        Self {
            syng_padding: 120,
            max_depth: 1,
            k_near: 3,
            k_far: 1,
            random_fraction: 0.01,
        }
    }
}

#[derive(Clone, Debug)]
pub enum LocalSeedRoute {
    WholeRegionSweepgaSeqwish,
    SyngAnchorSeededSeqwish(SyngAnchorSeedConfig),
    SyngDerived {
        mode: commands::syng2gfa::SyngGfaMode,
        params: syng::SyncmerParams,
        frequency_mask: commands::syng2gfa::SyngGfaFrequencyMask,
    },
}

impl LocalSeedRoute {
    pub fn label(&self) -> &'static str {
        match self {
            Self::WholeRegionSweepgaSeqwish => "whole-region-sweepga-seqwish",
            Self::SyngAnchorSeededSeqwish(_) => "syng-anchor-seeded-seqwish",
            Self::SyngDerived { .. } => "syng-derived",
        }
    }
}

#[derive(Clone)]
pub struct LocalSeedInductionConfig {
    pub route: LocalSeedRoute,
    pub graph_config: commands::graph::GraphBuildConfig,
}

impl LocalSeedInductionConfig {
    pub fn whole_region_sweepga_seqwish(graph_config: commands::graph::GraphBuildConfig) -> Self {
        Self {
            route: LocalSeedRoute::WholeRegionSweepgaSeqwish,
            graph_config,
        }
    }
}

#[derive(Clone, Debug)]
pub struct LocalSeedReport {
    pub seed_source: String,
    pub command_config: String,
    pub sequence_count: usize,
    pub total_bp: usize,
    pub raw_paf_records: Option<usize>,
    pub path_count: usize,
    pub path_validation: String,
    pub elapsed_secs: f64,
}

#[derive(Clone, Debug)]
pub struct LocalSeedGraph {
    pub gfa: String,
    pub report: LocalSeedReport,
}

pub fn induce_seed_graph(
    sequence_set: &LocalSequenceSet,
    config: &LocalSeedInductionConfig,
    syng_index: Option<&syng::SyngIndex>,
) -> io::Result<LocalSeedGraph> {
    let total_start = Instant::now();
    let sequence_count = sequence_set.len();
    let total_bp = sequence_set.total_bp();
    log::info!(
        "[local seed] starting route={} sequences={} total_bp={} config={}",
        config.route.label(),
        sequence_count,
        total_bp,
        command_config(&config.route, &config.graph_config)
    );

    let (gfa, raw_paf_records) = match &config.route {
        LocalSeedRoute::WholeRegionSweepgaSeqwish => {
            build_whole_region_sweepga_seqwish_seed(sequence_set, &config.graph_config)?
        }
        LocalSeedRoute::SyngAnchorSeededSeqwish(anchor_config) => {
            build_syng_anchor_seeded_seqwish_seed(
                sequence_set,
                &config.graph_config,
                syng_index.ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        "syng-anchor-seeded seed route requires a syng index",
                    )
                })?,
                anchor_config,
            )?
        }
        LocalSeedRoute::SyngDerived {
            mode,
            params,
            frequency_mask,
        } => {
            let mut out = Vec::new();
            crate::write_syng_region_gfa_from_sequences_with_params(
                &sequence_set.named_sequences(),
                &mut out,
                *params,
                *mode,
                frequency_mask.clone(),
            )?;
            let gfa = String::from_utf8(out).map_err(|err| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("local syng-derived seed GFA is not UTF-8: {err}"),
                )
            })?;
            (gfa, None)
        }
    };

    let path_count = validate_seed_path_sequences(sequence_set, &gfa)?;
    let elapsed_secs = total_start.elapsed().as_secs_f64();
    let report = LocalSeedReport {
        seed_source: config.route.label().to_string(),
        command_config: command_config(&config.route, &config.graph_config),
        sequence_count,
        total_bp,
        raw_paf_records,
        path_count,
        path_validation: "pass".to_string(),
        elapsed_secs,
    };
    log::info!(
        "[local seed] completed route={} sequences={} total_bp={} raw_paf_records={} paths={} path_validation={} elapsed={:.3}s",
        report.seed_source,
        report.sequence_count,
        report.total_bp,
        report
            .raw_paf_records
            .map(|value| value.to_string())
            .unwrap_or_else(|| "n/a".to_string()),
        report.path_count,
        report.path_validation,
        report.elapsed_secs,
    );
    write_debug_report(&config.graph_config, &report, &gfa)?;
    Ok(LocalSeedGraph { gfa, report })
}

fn build_whole_region_sweepga_seqwish_seed(
    sequence_set: &LocalSequenceSet,
    graph_config: &commands::graph::GraphBuildConfig,
) -> io::Result<(String, Option<usize>)> {
    if sequence_set.is_empty() {
        return Ok((String::from("H\tVN:Z:1.0\n"), Some(0)));
    }

    let named = sequence_set.sweepga_named_refs();
    let align_config = sweepga_align_config(sequence_set, graph_config);
    let align_start = Instant::now();
    let paf_file = sweepga::library_api::sweepga_align(&named, &align_config)
        .map_err(|err| io::Error::other(format!("local seed SweepGA alignment failed: {err}")))?;
    let align_elapsed = align_start.elapsed().as_secs_f64();
    log::info!(
        "[local seed] SweepGA/FastGA invocation returned PAF path={} elapsed={:.3}s",
        paf_file.path().display(),
        align_elapsed
    );
    let paf_load_start = Instant::now();
    let mut paf = String::new();
    std::fs::File::open(paf_file.path())?.read_to_string(&mut paf)?;
    let paf_bytes = paf.len();
    let mut paf_lines = paf
        .lines()
        .filter(|line| !line.trim().is_empty())
        .map(str::to_owned)
        .collect::<Vec<_>>();
    paf_lines.sort_unstable();
    paf_lines.dedup();
    let raw_paf_records = paf_lines.len();
    let mut paf = paf_lines.join("\n");
    if !paf.is_empty() {
        paf.push('\n');
    }
    let inter_sequence_paf_records = count_inter_sequence_paf_records(&paf);
    log::info!(
        "[local seed] loaded/sorted/deduped raw PAF bytes={} unique_records={} inter_sequence_records={} elapsed={:.3}s",
        paf_bytes,
        raw_paf_records,
        inter_sequence_paf_records,
        paf_load_start.elapsed().as_secs_f64()
    );
    log::info!(
        "[local seed] SweepGA/FastGA alignment emitted {} unique raw PAF record(s), {} inter-sequence record(s) in {:.3}s",
        raw_paf_records,
        inter_sequence_paf_records,
        align_elapsed + paf_load_start.elapsed().as_secs_f64()
    );
    if inter_sequence_paf_records == 0 {
        log::info!(
            "[local seed] no inter-sequence SweepGA/FastGA alignments; emitting disconnected exact seed graph"
        );
        return Ok((
            build_disconnected_sequence_gfa(sequence_set),
            Some(raw_paf_records),
        ));
    }

    let induce_start = Instant::now();
    let gfa = syng_graph::build_gfa_from_paf_and_sequences(
        &sequence_set.named_sequences(),
        &paf,
        graph_config,
    )?;
    log::info!(
        "[local seed] seqwish induction from local SweepGA/FastGA PAF completed in {:.3}s",
        induce_start.elapsed().as_secs_f64()
    );
    Ok((gfa, Some(raw_paf_records)))
}

fn build_syng_anchor_seeded_seqwish_seed(
    sequence_set: &LocalSequenceSet,
    graph_config: &commands::graph::GraphBuildConfig,
    syng_index: &syng::SyngIndex,
    anchor_config: &SyngAnchorSeedConfig,
) -> io::Result<(String, Option<usize>)> {
    if sequence_set.is_empty() {
        return Ok((String::from("H\tVN:Z:1.0\n"), Some(0)));
    }
    let members = sequence_set.syng_anchor_members();
    let seed = members
        .first()
        .map(|(_, _, member)| member)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "empty local seed members"))?;
    let paf = syng_graph::build_paf_anchor_seeded(
        &members,
        &seed.chrom,
        seed.fwd_start,
        seed.fwd_end,
        syng_index,
        anchor_config.syng_padding,
        anchor_config.max_depth,
        anchor_config.k_near,
        anchor_config.k_far,
        anchor_config.random_fraction,
    );
    let raw_paf_records = count_paf_records(&paf);
    let gfa = syng_graph::build_gfa_from_paf_and_sequences(
        &sequence_set.named_sequences(),
        &paf,
        graph_config,
    )?;
    Ok((gfa, Some(raw_paf_records)))
}

fn sweepga_align_config(
    sequence_set: &LocalSequenceSet,
    graph_config: &commands::graph::GraphBuildConfig,
) -> sweepga::library_api::SweepgaAlignConfig {
    let aligner = graph_config.aligner.clone();
    sweepga::library_api::SweepgaAlignConfig {
        num_threads: graph_config.num_threads,
        kmer_frequency: local_seed_kmer_frequency(sequence_set, graph_config),
        min_aln_length: graph_config.min_aln_length,
        // Keep this first alignment pass raw; the shared seqwish tail below
        // applies the documented GraphBuildConfig filter exactly once.
        no_filter: true,
        num_mappings: "many:many".to_string(),
        scaffold_filter: "many:many".to_string(),
        scaffold_mass: 0,
        sparsify: graph_config.sparsify.clone(),
        mash_params: graph_config.mash_params.clone(),
        aligner: aligner.clone(),
        temp_dir: graph_config.temp_dir.clone(),
        map_pct_identity: if aligner.eq_ignore_ascii_case("wfmash") {
            graph_config.map_pct_identity.clone()
        } else {
            None
        },
        batch_bytes: graph_config.batch_bytes.clone(),
        ..sweepga::library_api::SweepgaAlignConfig::default()
    }
}

fn local_seed_kmer_frequency(
    sequence_set: &LocalSequenceSet,
    graph_config: &commands::graph::GraphBuildConfig,
) -> usize {
    if let Some(value) = graph_config.frequency {
        return value.max(1);
    }
    let haplotypes = sweepga::pansn::count_pansn_keys(
        sequence_set
            .records()
            .iter()
            .map(|record| record.path_name.as_str()),
        sweepga::pansn::PanSnLevel::Haplotype,
    );
    haplotypes
        .saturating_mul(graph_config.frequency_multiplier.max(1))
        .max(1)
}

pub fn validate_seed_path_sequences(
    sequence_set: &LocalSequenceSet,
    gfa: &str,
) -> io::Result<usize> {
    let expected = sequence_set.expected_path_sequences();
    let observed = resolution::path_sequences(gfa)?
        .into_iter()
        .collect::<BTreeMap<_, _>>();

    let missing = expected
        .keys()
        .filter(|name| !observed.contains_key(*name))
        .cloned()
        .collect::<Vec<_>>();
    let extra = observed
        .keys()
        .filter(|name| !expected.contains_key(*name))
        .cloned()
        .collect::<Vec<_>>();
    let mismatched = expected
        .iter()
        .filter_map(|(name, expected_seq)| {
            observed
                .get(name)
                .filter(|observed_seq| *observed_seq != expected_seq)
                .map(|observed_seq| (name.clone(), expected_seq.len(), observed_seq.len()))
        })
        .collect::<Vec<_>>();

    if !missing.is_empty() || !extra.is_empty() || !mismatched.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "local seed path corruption: missing={missing:?} extra={extra:?} mismatched={mismatched:?}"
            ),
        ));
    }

    Ok(observed.len())
}

fn count_paf_records(paf: &str) -> usize {
    paf.lines().filter(|line| !line.trim().is_empty()).count()
}

fn count_inter_sequence_paf_records(paf: &str) -> usize {
    paf.lines()
        .filter(|line| {
            let mut fields = line.split('\t');
            let Some(query_name) = fields.next() else {
                return false;
            };
            let Some(target_name) = fields.nth(4) else {
                return false;
            };
            query_name != target_name
        })
        .count()
}

fn build_disconnected_sequence_gfa(sequence_set: &LocalSequenceSet) -> String {
    let mut out = String::from("H\tVN:Z:1.0\n");
    for (idx, record) in sequence_set.records().iter().enumerate() {
        let node_id = idx + 1;
        let sequence = String::from_utf8_lossy(&record.sequence);
        let _ = writeln!(out, "S\t{node_id}\t{sequence}");
        let _ = writeln!(out, "P\t{}\t{}+\t*", record.path_name, node_id);
    }
    out
}

fn command_config(
    route: &LocalSeedRoute,
    graph_config: &commands::graph::GraphBuildConfig,
) -> String {
    let mut out = String::new();
    let _ = write!(
        out,
        "impg-local-seed route={} aligner={} threads={} min-aln-len={} min-match-len={} no-filter={} num-mappings={} scaffold-jump={} scaffold-mass={} scaffold-filter={} min-map-length={} min-identity={} sparsify={} repeat-max={} min-repeat-dist={} transclose-batch={} disk-backed={}",
        route.label(),
        graph_config.aligner,
        graph_config.num_threads,
        graph_config.min_aln_length,
        graph_config.min_match_len,
        graph_config.no_filter,
        graph_config.num_mappings,
        graph_config.scaffold_jump,
        graph_config.scaffold_mass,
        graph_config.scaffold_filter,
        graph_config.min_map_length,
        graph_config.min_identity,
        graph_config.sparsify.description(),
        graph_config.repeat_max,
        graph_config.min_repeat_dist,
        graph_config.transclose_batch,
        graph_config.disk_backed,
    );
    match route {
        LocalSeedRoute::WholeRegionSweepgaSeqwish => {}
        LocalSeedRoute::SyngAnchorSeededSeqwish(config) => {
            let _ = write!(
                out,
                " syng-padding={} max-depth={} k-near={} k-far={} random-fraction={}",
                config.syng_padding,
                config.max_depth,
                config.k_near,
                config.k_far,
                config.random_fraction
            );
        }
        LocalSeedRoute::SyngDerived {
            mode,
            params,
            frequency_mask,
        } => {
            let _ = write!(
                out,
                " syng-mode={} syncmer-k={} smer-k={} seed={} mask-enabled={} mask-top={} mask-min-shared-run={}",
                mode.label(),
                params.k + params.w,
                params.k,
                params.seed,
                frequency_mask.enabled(),
                frequency_mask.drop_top_fraction,
                frequency_mask.min_shared_run,
            );
        }
    }
    out
}

fn write_debug_report(
    graph_config: &commands::graph::GraphBuildConfig,
    report: &LocalSeedReport,
    gfa: &str,
) -> io::Result<()> {
    let Some(debug_dir) = graph_config.debug_dir.as_deref() else {
        return Ok(());
    };
    std::fs::create_dir_all(debug_dir)?;
    let path = Path::new(debug_dir).join("local_seed_report.tsv");
    let content = format!(
        concat!(
            "seed_source\t{}\n",
            "command_config\t{}\n",
            "sequence_count\t{}\n",
            "total_bp\t{}\n",
            "raw_paf_records\t{}\n",
            "path_count\t{}\n",
            "path_validation\t{}\n",
            "elapsed_secs\t{:.6}\n"
        ),
        report.seed_source,
        report.command_config,
        report.sequence_count,
        report.total_bp,
        report
            .raw_paf_records
            .map(|value| value.to_string())
            .unwrap_or_else(|| "n/a".to_string()),
        report.path_count,
        report.path_validation,
        report.elapsed_secs,
    );
    std::fs::write(path, content)?;
    write_dirty_region_debug_reports(debug_dir, gfa);
    Ok(())
}

fn write_dirty_region_debug_reports(debug_dir: &str, gfa: &str) {
    let options = crate::graph_badness::DirtyRegionOptions::default();
    let report = match crate::graph_badness::analyze_gfa("local_seed", gfa, &options) {
        Ok(report) => report,
        Err(err) => {
            log::warn!("[local seed] dirty-region diagnostics failed: {err}");
            return;
        }
    };
    let debug_dir = Path::new(debug_dir);
    for format in ["json", "tsv"] {
        let text = match crate::graph_badness::format_dirty_report(&report, format) {
            Ok(text) => text,
            Err(err) => {
                log::warn!("[local seed] dirty-region {format} formatting failed: {err}");
                continue;
            }
        };
        let path = debug_dir.join(format!("local_seed_dirty_regions.{format}"));
        if let Err(err) = std::fs::write(&path, text) {
            log::warn!(
                "[local seed] failed to write dirty-region report {}: {err}",
                path.display()
            );
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use sweepga::knn_graph::SparsificationStrategy;

    fn graph_config() -> commands::graph::GraphBuildConfig {
        commands::graph::GraphBuildConfig {
            num_threads: 1,
            aligner: "fastga".to_string(),
            min_match_len: 1,
            no_filter: true,
            min_map_length: 1,
            sparsify: SparsificationStrategy::None,
            show_progress: false,
            ..commands::graph::GraphBuildConfig::default()
        }
    }

    fn make_sequence(len: usize, seed: u32) -> Vec<u8> {
        let mut state = seed;
        (0..len)
            .map(|_| {
                state = state.wrapping_mul(1664525).wrapping_add(1013904223);
                b"ACGT"[((state >> 24) & 3) as usize]
            })
            .collect()
    }

    fn small_sequence_set() -> LocalSequenceSet {
        let seq_a = make_sequence(512, 17);
        let mut seq_b = seq_a.clone();
        seq_b[128] = if seq_b[128] == b'A' { b'C' } else { b'A' };
        seq_b[320] = if seq_b[320] == b'G' { b'T' } else { b'G' };
        LocalSequenceSet::new(vec![
            LocalSequenceRecord {
                path_name: "HG001#1#chr6:10-522".to_string(),
                source_name: "HG001#1#chr6".to_string(),
                source_start: 10,
                source_end: 522,
                strand: '+',
                sequence: seq_a,
            },
            LocalSequenceRecord {
                path_name: "HG002#1#chr6:10-522".to_string(),
                source_name: "HG002#1#chr6".to_string(),
                source_start: 10,
                source_end: 522,
                strand: '+',
                sequence: seq_b,
            },
        ])
    }

    #[test]
    fn whole_region_sweepga_seqwish_seed_preserves_path_names_and_sequences() {
        let sequence_set = small_sequence_set();
        let config = LocalSeedInductionConfig::whole_region_sweepga_seqwish(graph_config());

        let seed = induce_seed_graph(&sequence_set, &config, None)
            .expect("whole-region SweepGA/seqwish local seed should build");
        let paths = resolution::path_sequences(&seed.gfa)
            .expect("seed GFA paths should parse")
            .into_iter()
            .collect::<BTreeMap<_, _>>();

        assert_eq!(
            paths.keys().cloned().collect::<Vec<_>>(),
            vec![
                "HG001#1#chr6:10-522".to_string(),
                "HG002#1#chr6:10-522".to_string()
            ]
        );
        assert_eq!(
            paths["HG001#1#chr6:10-522"],
            String::from_utf8(sequence_set.records()[0].sequence.clone()).unwrap()
        );
        assert_eq!(
            paths["HG002#1#chr6:10-522"],
            String::from_utf8(sequence_set.records()[1].sequence.clone()).unwrap()
        );
        assert!(paths
            .keys()
            .all(|name| !name.starts_with("local_") && !name.starts_with("__impg")));
        assert_eq!(seed.report.seed_source, "whole-region-sweepga-seqwish");
        assert!(
            seed.report.raw_paf_records.unwrap_or(0) > 0,
            "varied local seed fixture should exercise the SweepGA/FastGA PAF route"
        );
        assert!(seed.report.command_config.contains("aligner=fastga"));
        assert_eq!(seed.report.path_validation, "pass");
        assert!(seed.report.elapsed_secs >= 0.0);
    }

    #[test]
    fn exact_path_corruption_is_the_seed_driver_hard_gate() {
        let sequence_set = small_sequence_set();
        let corrupt = concat!(
            "H\tVN:Z:1.0\n",
            "S\t1\tACGT\n",
            "P\tHG001#1#chr6:10-74\t1+\t*\n",
            "P\tlocal_1\t1+\t*\n"
        );

        let err = validate_seed_path_sequences(&sequence_set, corrupt).unwrap_err();
        assert!(err.to_string().contains("local seed path corruption"));
    }
}
