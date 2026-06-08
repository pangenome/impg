//! Diagnostic graph-quality signals and dirty-region extraction.
//!
//! This layer is intentionally diagnostic-only.  It converts local seed GFA
//! signals into path-coordinate candidate chunks for later local polishing, but
//! it does not accept or reject graph replacement decisions.

use crate::graph_report::{describe_gfa, GraphReportOptions};
use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fmt::Write as _;
use std::io;

#[derive(Clone, Debug)]
pub struct DirtyRegionOptions {
    pub graph_report: GraphReportOptions,
    pub top_n: usize,
    pub min_white_space_gap_bp: usize,
    pub min_sparse_run_bp: usize,
    pub max_sparse_coverage_path_fraction: f64,
    pub min_path_jump: usize,
    pub min_link_jump: usize,
    pub low_depth_max_visits: usize,
    pub min_low_depth_run_bp: usize,
    pub small_loop_max_steps: usize,
    pub merge_distance_bp: usize,
    pub flank_bp: usize,
    pub max_chunk_bp: Option<usize>,
}

impl Default for DirtyRegionOptions {
    fn default() -> Self {
        let graph_report = GraphReportOptions {
            top_n: 64,
            min_white_space_gap_bp: 1_000,
            max_sparse_coverage_path_fraction: 0.25,
            ..GraphReportOptions::default()
        };
        Self {
            graph_report,
            top_n: 64,
            min_white_space_gap_bp: 1_000,
            min_sparse_run_bp: 1_000,
            max_sparse_coverage_path_fraction: 0.25,
            min_path_jump: 1_000,
            min_link_jump: 1_000,
            low_depth_max_visits: 1,
            min_low_depth_run_bp: 1_000,
            small_loop_max_steps: 4,
            merge_distance_bp: 500,
            flank_bp: 1_000,
            max_chunk_bp: Some(50_000),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct DirtyRegionConfigSummary {
    pub top_n: usize,
    pub min_white_space_gap_bp: usize,
    pub min_sparse_run_bp: usize,
    pub max_sparse_coverage_path_fraction: f64,
    pub min_path_jump: usize,
    pub min_link_jump: usize,
    pub low_depth_max_visits: usize,
    pub min_low_depth_run_bp: usize,
    pub small_loop_max_steps: usize,
    pub merge_distance_bp: usize,
    pub flank_bp: usize,
    pub max_chunk_bp: Option<usize>,
}

impl From<&DirtyRegionOptions> for DirtyRegionConfigSummary {
    fn from(options: &DirtyRegionOptions) -> Self {
        Self {
            top_n: options.top_n,
            min_white_space_gap_bp: options.min_white_space_gap_bp,
            min_sparse_run_bp: options.min_sparse_run_bp,
            max_sparse_coverage_path_fraction: options.max_sparse_coverage_path_fraction,
            min_path_jump: options.min_path_jump,
            min_link_jump: options.min_link_jump,
            low_depth_max_visits: options.low_depth_max_visits,
            min_low_depth_run_bp: options.min_low_depth_run_bp,
            small_loop_max_steps: options.small_loop_max_steps,
            merge_distance_bp: options.merge_distance_bp,
            flank_bp: options.flank_bp,
            max_chunk_bp: options.max_chunk_bp,
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct DirtyRegionReport {
    pub source: String,
    pub replacement_decision: String,
    pub metrics_are_diagnostic_only: bool,
    pub graph_report_status: String,
    pub graph_report_failures: Vec<String>,
    pub graph_report_warnings: Vec<String>,
    pub options: DirtyRegionConfigSummary,
    pub metrics: GraphBadnessMetrics,
    pub diagnostics: Vec<DirtyDiagnostic>,
    pub dirty_sites: Vec<DirtySite>,
    pub candidate_chunks: Vec<DirtyChunk>,
}

#[derive(Clone, Debug, Serialize)]
pub struct GraphBadnessMetrics {
    pub segment_count: usize,
    pub segment_bp: usize,
    pub link_count: usize,
    pub path_count: usize,
    pub total_path_steps: usize,
    pub spelled_path_bp: u64,
    pub path_replay_compression_ratio: Option<f64>,
    pub node_path_step_visits: usize,
    pub bp_weighted_path_depth_mean: f64,
    pub singleton_nodes: usize,
    pub singleton_bp: usize,
    pub low_depth_nodes: usize,
    pub low_depth_bp: usize,
    pub bp_weighted_path_depth_distribution: Vec<PathDepthBucket>,
    pub link_jump_p99: usize,
    pub link_jump_max: usize,
    pub path_jump_p99: usize,
    pub path_jump_max: usize,
    pub white_space_proxy_bp_total: u64,
    pub white_space_proxy_bp_p99: usize,
    pub white_space_proxy_bp_max: usize,
    pub white_space_bridges_ge_threshold: usize,
    pub sparse_coverage_segment_bp: usize,
    pub direct_self_loop_edges: usize,
    pub direct_self_loop_nodes: usize,
    pub adjacent_same_node_path_steps: usize,
    pub small_loop_sites: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct PathDepthBucket {
    pub label: String,
    pub min: usize,
    pub max: Option<usize>,
    pub segments: usize,
    pub bp: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct DirtyDiagnostic {
    pub id: String,
    pub kind: String,
    pub source_metric: String,
    pub metric_value: f64,
    pub severity: f64,
    pub message: String,
    pub path_name: Option<String>,
    pub path_start_bp: Option<usize>,
    pub path_end_bp: Option<usize>,
    pub graph_start_bp: Option<usize>,
    pub graph_end_bp: Option<usize>,
}

#[derive(Clone, Debug, Serialize)]
pub struct DirtySite {
    pub id: String,
    pub kind: String,
    pub source_metric: String,
    pub reason: String,
    pub metric_value: f64,
    pub severity: f64,
    pub path_name: Option<String>,
    pub path_start_bp: Option<usize>,
    pub path_end_bp: Option<usize>,
    pub path_start_step: Option<usize>,
    pub path_end_step: Option<usize>,
    pub path_length_bp: Option<usize>,
    pub graph_start_bp: Option<usize>,
    pub graph_end_bp: Option<usize>,
    pub start_order: Option<usize>,
    pub end_order: Option<usize>,
    pub nodes: Vec<String>,
}

#[derive(Clone, Debug, Serialize)]
pub struct DirtyChunk {
    pub id: String,
    pub path_name: Option<String>,
    pub path_start_bp: Option<usize>,
    pub path_end_bp: Option<usize>,
    pub path_length_bp: Option<usize>,
    pub graph_start_bp: Option<usize>,
    pub graph_end_bp: Option<usize>,
    pub start_order: Option<usize>,
    pub end_order: Option<usize>,
    pub site_count: usize,
    pub site_ids: Vec<String>,
    pub kinds: Vec<String>,
    pub nodes: Vec<String>,
    pub max_severity: f64,
    pub flank_bp: usize,
    pub capped_by_budget: bool,
}

#[derive(Clone, Debug)]
struct Step {
    node: String,
}

#[derive(Clone, Debug)]
struct Link {
    from: String,
    to: String,
}

#[derive(Clone, Debug)]
struct PathWalk {
    name: String,
    steps: Vec<Step>,
    step_starts: Vec<usize>,
    length_bp: usize,
}

impl PathWalk {
    fn new(name: String, steps: Vec<Step>, segment_lengths: &HashMap<String, usize>) -> Self {
        let mut step_starts = Vec::with_capacity(steps.len());
        let mut offset = 0usize;
        for step in &steps {
            step_starts.push(offset);
            offset = offset.saturating_add(segment_lengths.get(&step.node).copied().unwrap_or(0));
        }
        Self {
            name,
            steps,
            step_starts,
            length_bp: offset,
        }
    }

    fn step_end_bp(&self, step_exclusive: usize) -> usize {
        if step_exclusive >= self.steps.len() {
            self.length_bp
        } else {
            self.step_starts[step_exclusive]
        }
    }
}

#[derive(Clone, Debug, Default)]
struct ParsedGfa {
    order: HashMap<String, usize>,
    ordered_nodes: Vec<String>,
    segment_lengths: HashMap<String, usize>,
    links: Vec<Link>,
    paths: Vec<PathWalk>,
}

#[derive(Clone, Debug)]
struct OrderedNodeRun {
    nodes: Vec<String>,
    start_order: usize,
    end_order: usize,
    bp: usize,
}

#[derive(Default)]
struct ChunkAcc {
    path_name: Option<String>,
    raw_path_start_bp: Option<usize>,
    raw_path_end_bp: Option<usize>,
    path_length_bp: Option<usize>,
    graph_start_bp: Option<usize>,
    graph_end_bp: Option<usize>,
    start_order: Option<usize>,
    end_order: Option<usize>,
    path_start_step: Option<usize>,
    path_end_step: Option<usize>,
    site_ids: Vec<String>,
    kinds: BTreeSet<String>,
    nodes: BTreeSet<String>,
    max_severity: f64,
}

pub fn analyze_gfa(
    source: impl Into<String>,
    gfa_text: &str,
    options: &DirtyRegionOptions,
) -> io::Result<DirtyRegionReport> {
    let source = source.into();
    let parsed = parse_gfa(gfa_text);
    let graph_report = describe_gfa(source.clone(), gfa_text, &options.graph_report)?;
    let visits_by_node = total_visits_by_node(&parsed);
    let distinct_paths_by_node = distinct_paths_by_node(&parsed);
    let bp_starts = segment_bp_starts(&parsed);
    let mut sites = Vec::new();

    detect_white_space_sites(&parsed, &bp_starts, options, &mut sites);
    detect_path_jump_sites(&parsed, &bp_starts, options, &mut sites);
    detect_long_link_sites(&parsed, &bp_starts, options, &mut sites);
    detect_low_depth_sites(
        &parsed,
        &bp_starts,
        &visits_by_node,
        &distinct_paths_by_node,
        options,
        &mut sites,
    );
    detect_loop_sites(&parsed, &bp_starts, options, &mut sites);

    assign_site_ids(&mut sites);
    let candidate_chunks = merge_dirty_sites(&parsed, &sites, options);
    let diagnostics = sites.iter().map(diagnostic_from_site).collect::<Vec<_>>();
    let metrics = compute_badness_metrics(
        &parsed,
        &graph_report.metrics,
        &visits_by_node,
        &sites,
        options,
    );

    Ok(DirtyRegionReport {
        source,
        replacement_decision: "diagnostic-only".to_string(),
        metrics_are_diagnostic_only: true,
        graph_report_status: graph_report.status,
        graph_report_failures: graph_report.failures,
        graph_report_warnings: graph_report.warnings,
        options: DirtyRegionConfigSummary::from(options),
        metrics,
        diagnostics,
        dirty_sites: sites,
        candidate_chunks,
    })
}

pub fn format_dirty_report(report: &DirtyRegionReport, format: &str) -> io::Result<String> {
    match format {
        "json" => serde_json::to_string_pretty(report)
            .map(|mut out| {
                out.push('\n');
                out
            })
            .map_err(|err| io::Error::other(format!("failed to serialize dirty report: {err}"))),
        "tsv" => Ok(report.to_tsv()),
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported dirty report format '{other}'"),
        )),
    }
}

impl DirtyRegionReport {
    pub fn to_tsv(&self) -> String {
        let mut out = String::new();
        out.push_str(
            "record_type\tid\tsource\tkind\tsource_metric\tmetric_value\tseverity\tpath_name\tpath_start_bp\tpath_end_bp\tpath_start_step\tpath_end_step\tpath_length_bp\tgraph_start_bp\tgraph_end_bp\tstart_order\tend_order\tsite_count\tnodes\tevidence\n",
        );
        for (key, value) in self.metric_rows() {
            write_tsv_row(
                &mut out,
                &[
                    "metric",
                    &key,
                    &self.source,
                    "metric",
                    &key,
                    &value,
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                    "",
                ],
            );
        }
        for site in &self.dirty_sites {
            let nodes = site.nodes.join(",");
            write_tsv_row(
                &mut out,
                &[
                    "site",
                    &site.id,
                    &self.source,
                    &site.kind,
                    &site.source_metric,
                    &format_float(site.metric_value),
                    &format_float(site.severity),
                    site.path_name.as_deref().unwrap_or(""),
                    &option_usize(site.path_start_bp),
                    &option_usize(site.path_end_bp),
                    &option_usize(site.path_start_step),
                    &option_usize(site.path_end_step),
                    &option_usize(site.path_length_bp),
                    &option_usize(site.graph_start_bp),
                    &option_usize(site.graph_end_bp),
                    &option_usize(site.start_order),
                    &option_usize(site.end_order),
                    "1",
                    &nodes,
                    &site.reason,
                ],
            );
        }
        for chunk in &self.candidate_chunks {
            let nodes = chunk.nodes.join(",");
            let evidence = chunk.site_ids.join(",");
            write_tsv_row(
                &mut out,
                &[
                    "chunk",
                    &chunk.id,
                    &self.source,
                    &chunk.kinds.join(","),
                    "merged_dirty_sites",
                    "",
                    &format_float(chunk.max_severity),
                    chunk.path_name.as_deref().unwrap_or(""),
                    &option_usize(chunk.path_start_bp),
                    &option_usize(chunk.path_end_bp),
                    "",
                    "",
                    &option_usize(chunk.path_length_bp),
                    &option_usize(chunk.graph_start_bp),
                    &option_usize(chunk.graph_end_bp),
                    &option_usize(chunk.start_order),
                    &option_usize(chunk.end_order),
                    &chunk.site_count.to_string(),
                    &nodes,
                    &evidence,
                ],
            );
        }
        out
    }

    fn metric_rows(&self) -> Vec<(String, String)> {
        let m = &self.metrics;
        let mut rows = vec![
            (
                "replacement_decision".to_string(),
                self.replacement_decision.clone(),
            ),
            (
                "metrics_are_diagnostic_only".to_string(),
                self.metrics_are_diagnostic_only.to_string(),
            ),
            (
                "graph_report_status".to_string(),
                self.graph_report_status.clone(),
            ),
            ("segment_count".to_string(), m.segment_count.to_string()),
            ("segment_bp".to_string(), m.segment_bp.to_string()),
            ("link_count".to_string(), m.link_count.to_string()),
            ("path_count".to_string(), m.path_count.to_string()),
            (
                "total_path_steps".to_string(),
                m.total_path_steps.to_string(),
            ),
            ("spelled_path_bp".to_string(), m.spelled_path_bp.to_string()),
            (
                "path_replay_compression_ratio".to_string(),
                m.path_replay_compression_ratio
                    .map(format_float)
                    .unwrap_or_default(),
            ),
            (
                "bp_weighted_path_depth_mean".to_string(),
                format_float(m.bp_weighted_path_depth_mean),
            ),
            ("singleton_nodes".to_string(), m.singleton_nodes.to_string()),
            ("singleton_bp".to_string(), m.singleton_bp.to_string()),
            ("low_depth_nodes".to_string(), m.low_depth_nodes.to_string()),
            ("low_depth_bp".to_string(), m.low_depth_bp.to_string()),
            ("link_jump_p99".to_string(), m.link_jump_p99.to_string()),
            ("link_jump_max".to_string(), m.link_jump_max.to_string()),
            ("path_jump_p99".to_string(), m.path_jump_p99.to_string()),
            ("path_jump_max".to_string(), m.path_jump_max.to_string()),
            (
                "white_space_proxy_bp_total".to_string(),
                m.white_space_proxy_bp_total.to_string(),
            ),
            (
                "white_space_proxy_bp_p99".to_string(),
                m.white_space_proxy_bp_p99.to_string(),
            ),
            (
                "white_space_proxy_bp_max".to_string(),
                m.white_space_proxy_bp_max.to_string(),
            ),
            (
                "white_space_bridges_ge_threshold".to_string(),
                m.white_space_bridges_ge_threshold.to_string(),
            ),
            (
                "sparse_coverage_segment_bp".to_string(),
                m.sparse_coverage_segment_bp.to_string(),
            ),
            (
                "direct_self_loop_edges".to_string(),
                m.direct_self_loop_edges.to_string(),
            ),
            (
                "direct_self_loop_nodes".to_string(),
                m.direct_self_loop_nodes.to_string(),
            ),
            (
                "adjacent_same_node_path_steps".to_string(),
                m.adjacent_same_node_path_steps.to_string(),
            ),
            (
                "small_loop_sites".to_string(),
                m.small_loop_sites.to_string(),
            ),
        ];
        for bucket in &m.bp_weighted_path_depth_distribution {
            rows.push((
                format!("path_depth_bp_bucket:{}", bucket.label),
                bucket.bp.to_string(),
            ));
        }
        rows
    }
}

fn parse_gfa(gfa_text: &str) -> ParsedGfa {
    let mut order = HashMap::new();
    let mut ordered_nodes = Vec::new();
    let mut segment_lengths = HashMap::new();
    let mut raw_paths = Vec::<(String, Vec<Step>)>::new();
    let mut links = Vec::new();

    for line in gfa_text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        if fields.is_empty() {
            continue;
        }
        match fields[0] {
            "S" if fields.len() >= 3 => {
                let node = fields[1].to_string();
                let next = ordered_nodes.len();
                if order.insert(node.clone(), next).is_none() {
                    ordered_nodes.push(node.clone());
                }
                segment_lengths.insert(node, sequence_len(fields[2]));
            }
            "L" if fields.len() >= 5 => {
                links.push(Link {
                    from: fields[1].to_string(),
                    to: fields[3].to_string(),
                });
            }
            "P" if fields.len() >= 3 => {
                let steps = split_p_steps(fields[2]);
                if !steps.is_empty() {
                    raw_paths.push((fields[1].to_string(), steps));
                }
            }
            "W" if fields.len() >= 7 => {
                let steps = split_w_steps(fields[6]);
                if !steps.is_empty() {
                    raw_paths.push((w_line_name(&fields), steps));
                }
            }
            _ => {}
        }
    }

    let paths = raw_paths
        .into_iter()
        .map(|(name, steps)| PathWalk::new(name, steps, &segment_lengths))
        .collect();

    ParsedGfa {
        order,
        ordered_nodes,
        segment_lengths,
        links,
        paths,
    }
}

fn sequence_len(sequence: &str) -> usize {
    if sequence == "*" {
        0
    } else {
        sequence.len()
    }
}

fn split_p_steps(raw: &str) -> Vec<Step> {
    if raw == "*" {
        return Vec::new();
    }
    raw.split(',')
        .filter_map(|token| {
            let (node, orientation) = token.split_at(token.len().checked_sub(1)?);
            let orientation = orientation.chars().next()?;
            matches!(orientation, '+' | '-').then(|| Step {
                node: node.to_string(),
            })
        })
        .collect()
}

fn split_w_steps(raw: &str) -> Vec<Step> {
    if raw == "*" {
        return Vec::new();
    }
    let mut steps = Vec::new();
    let mut current_orientation: Option<char> = None;
    let mut current_start = 0usize;
    for (idx, ch) in raw.char_indices() {
        if !matches!(ch, '>' | '<') {
            continue;
        }
        if current_orientation.is_some() {
            let node = &raw[current_start..idx];
            if !node.is_empty() {
                steps.push(Step {
                    node: node.to_string(),
                });
            }
        }
        current_orientation = Some(if ch == '>' { '+' } else { '-' });
        current_start = idx + ch.len_utf8();
    }
    if current_orientation.is_some() {
        let node = &raw[current_start..];
        if !node.is_empty() {
            steps.push(Step {
                node: node.to_string(),
            });
        }
    }
    steps
}

fn w_line_name(fields: &[&str]) -> String {
    let sample = fields.get(1).copied().unwrap_or("*");
    let hap = fields.get(2).copied().unwrap_or("*");
    let seq = fields.get(3).copied().unwrap_or("*");
    let start = fields.get(4).copied().unwrap_or("*");
    let end = fields.get(5).copied().unwrap_or("*");
    format!("{sample}#{hap}#{seq}:{start}-{end}")
}

fn segment_bp_starts(parsed: &ParsedGfa) -> HashMap<String, usize> {
    let mut starts = HashMap::with_capacity(parsed.ordered_nodes.len());
    let mut offset = 0usize;
    for node in &parsed.ordered_nodes {
        starts.insert(node.clone(), offset);
        offset = offset.saturating_add(parsed.segment_lengths.get(node).copied().unwrap_or(0));
    }
    starts
}

fn total_visits_by_node(parsed: &ParsedGfa) -> HashMap<String, usize> {
    let mut visits = HashMap::new();
    for path in &parsed.paths {
        for step in &path.steps {
            *visits.entry(step.node.clone()).or_insert(0) += 1;
        }
    }
    visits
}

fn distinct_paths_by_node(parsed: &ParsedGfa) -> HashMap<String, usize> {
    let mut counts = HashMap::new();
    for path in &parsed.paths {
        let mut seen = HashSet::new();
        for step in &path.steps {
            if seen.insert(step.node.as_str()) {
                *counts.entry(step.node.clone()).or_insert(0) += 1;
            }
        }
    }
    counts
}

fn compute_badness_metrics(
    parsed: &ParsedGfa,
    graph_metrics: &crate::graph_report::GraphMetrics,
    visits_by_node: &HashMap<String, usize>,
    sites: &[DirtySite],
    options: &DirtyRegionOptions,
) -> GraphBadnessMetrics {
    let segment_bp = parsed.segment_lengths.values().copied().sum::<usize>();
    let spelled_path_bp = parsed
        .paths
        .iter()
        .map(|path| path.length_bp as u64)
        .sum::<u64>();
    let path_replay_compression_ratio =
        (segment_bp > 0).then_some(spelled_path_bp as f64 / segment_bp as f64);
    let mut low_depth_nodes = 0usize;
    let mut low_depth_bp = 0usize;
    for node in &parsed.ordered_nodes {
        let visits = visits_by_node.get(node).copied().unwrap_or(0);
        if visits <= options.low_depth_max_visits {
            low_depth_nodes += 1;
            low_depth_bp =
                low_depth_bp.saturating_add(parsed.segment_lengths.get(node).copied().unwrap_or(0));
        }
    }
    let small_loop_sites = sites
        .iter()
        .filter(|site| site.kind == "small_loop")
        .count();

    GraphBadnessMetrics {
        segment_count: graph_metrics.segments,
        segment_bp: graph_metrics.total_segment_bp,
        link_count: graph_metrics.links,
        path_count: graph_metrics.paths,
        total_path_steps: graph_metrics.path_steps,
        spelled_path_bp,
        path_replay_compression_ratio,
        node_path_step_visits: graph_metrics.node_path_step_visits,
        bp_weighted_path_depth_mean: graph_metrics.node_coverage_bp_weighted_mean,
        singleton_nodes: graph_metrics.singleton_nodes,
        singleton_bp: graph_metrics.singleton_bp,
        low_depth_nodes,
        low_depth_bp,
        bp_weighted_path_depth_distribution: path_depth_distribution(parsed, visits_by_node),
        link_jump_p99: graph_metrics.link_jump_p99,
        link_jump_max: graph_metrics.link_jump_max,
        path_jump_p99: graph_metrics.path_jump_p99,
        path_jump_max: graph_metrics.path_jump_max,
        white_space_proxy_bp_total: graph_metrics.path_white_space_bp_total,
        white_space_proxy_bp_p99: graph_metrics.path_white_space_bp_p99,
        white_space_proxy_bp_max: graph_metrics.path_white_space_bp_max,
        white_space_bridges_ge_threshold: graph_metrics.path_white_space_bridges_ge_threshold,
        sparse_coverage_segment_bp: graph_metrics.sparse_coverage_segment_bp,
        direct_self_loop_edges: graph_metrics.direct_self_loop_edges,
        direct_self_loop_nodes: graph_metrics.direct_self_loop_nodes,
        adjacent_same_node_path_steps: graph_metrics.adjacent_same_node_path_steps,
        small_loop_sites,
    }
}

fn path_depth_distribution(
    parsed: &ParsedGfa,
    visits_by_node: &HashMap<String, usize>,
) -> Vec<PathDepthBucket> {
    let mut buckets = [
        ("0", 0, Some(0)),
        ("1", 1, Some(1)),
        ("2", 2, Some(2)),
        ("3-4", 3, Some(4)),
        ("5-8", 5, Some(8)),
        ("9-16", 9, Some(16)),
        ("17-32", 17, Some(32)),
        ("33-64", 33, Some(64)),
        ("65-128", 65, Some(128)),
        ("129+", 129, None),
    ]
    .into_iter()
    .map(|(label, min, max)| PathDepthBucket {
        label: label.to_string(),
        min,
        max,
        segments: 0,
        bp: 0,
    })
    .collect::<Vec<_>>();

    for node in &parsed.ordered_nodes {
        let visits = visits_by_node.get(node).copied().unwrap_or(0);
        let bp = parsed.segment_lengths.get(node).copied().unwrap_or(0);
        let idx = match visits {
            0 => 0,
            1 => 1,
            2 => 2,
            3..=4 => 3,
            5..=8 => 4,
            9..=16 => 5,
            17..=32 => 6,
            33..=64 => 7,
            65..=128 => 8,
            _ => 9,
        };
        if let Some(bucket) = buckets.get_mut(idx) {
            bucket.segments += 1;
            bucket.bp = bucket.bp.saturating_add(bp);
        }
    }
    buckets
}

fn detect_white_space_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    options: &DirtyRegionOptions,
    sites: &mut Vec<DirtySite>,
) {
    if options.min_white_space_gap_bp == 0 {
        return;
    }
    for path in &parsed.paths {
        for (idx, pair) in path.steps.windows(2).enumerate() {
            let Some((gap_bp, _, _)) =
                ordered_gap_bp(parsed, bp_starts, &pair[0].node, &pair[1].node)
            else {
                continue;
            };
            if gap_bp < options.min_white_space_gap_bp {
                continue;
            }
            sites.push(path_site(
                parsed,
                bp_starts,
                "underaligned_white_space",
                "path_white_space_bp",
                format!("path step bridges {gap_bp} bp of segment-order white space after sorting"),
                gap_bp as f64,
                gap_bp as f64 / options.min_white_space_gap_bp.max(1) as f64,
                path,
                idx,
                idx + 2,
                vec![pair[0].node.clone(), pair[1].node.clone()],
            ));
        }
    }
}

fn detect_path_jump_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    options: &DirtyRegionOptions,
    sites: &mut Vec<DirtySite>,
) {
    if options.min_path_jump == 0 {
        return;
    }
    for path in &parsed.paths {
        for (idx, pair) in path.steps.windows(2).enumerate() {
            let Some(&from_order) = parsed.order.get(&pair[0].node) else {
                continue;
            };
            let Some(&to_order) = parsed.order.get(&pair[1].node) else {
                continue;
            };
            let jump = from_order.abs_diff(to_order);
            if jump < options.min_path_jump {
                continue;
            }
            sites.push(path_site(
                parsed,
                bp_starts,
                "long_path_jump",
                "path_jump",
                format!("consecutive path steps jump {jump} segment-order positions"),
                jump as f64,
                jump as f64 / options.min_path_jump.max(1) as f64,
                path,
                idx,
                idx + 2,
                vec![pair[0].node.clone(), pair[1].node.clone()],
            ));
        }
    }
}

fn detect_long_link_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    options: &DirtyRegionOptions,
    sites: &mut Vec<DirtySite>,
) {
    if options.min_link_jump == 0 {
        return;
    }
    for link in &parsed.links {
        let Some(&from_order) = parsed.order.get(&link.from) else {
            continue;
        };
        let Some(&to_order) = parsed.order.get(&link.to) else {
            continue;
        };
        let jump = from_order.abs_diff(to_order);
        if jump < options.min_link_jump {
            continue;
        }

        let mut emitted_path_site = false;
        for path in &parsed.paths {
            for (idx, pair) in path.steps.windows(2).enumerate() {
                if (pair[0].node == link.from && pair[1].node == link.to)
                    || (pair[0].node == link.to && pair[1].node == link.from)
                {
                    emitted_path_site = true;
                    sites.push(path_site(
                        parsed,
                        bp_starts,
                        "long_link",
                        "link_jump",
                        format!("GFA link spans {jump} segment-order positions"),
                        jump as f64,
                        jump as f64 / options.min_link_jump.max(1) as f64,
                        path,
                        idx,
                        idx + 2,
                        vec![pair[0].node.clone(), pair[1].node.clone()],
                    ));
                }
            }
        }
        if !emitted_path_site {
            sites.push(graph_site(
                parsed,
                bp_starts,
                "long_link",
                "link_jump",
                format!("unsupported GFA link spans {jump} segment-order positions"),
                jump as f64,
                jump as f64 / options.min_link_jump.max(1) as f64,
                vec![link.from.clone(), link.to.clone()],
            ));
        }
    }
}

fn detect_low_depth_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    visits_by_node: &HashMap<String, usize>,
    distinct_paths_by_node: &HashMap<String, usize>,
    options: &DirtyRegionOptions,
    sites: &mut Vec<DirtySite>,
) {
    let low_depth_runs = ordered_runs(parsed, |node| {
        visits_by_node.get(node).copied().unwrap_or(0) <= options.low_depth_max_visits
    });
    for run in low_depth_runs {
        if run.bp < options.min_low_depth_run_bp {
            continue;
        }
        emit_run_sites(
            parsed,
            bp_starts,
            sites,
            &run,
            "low_depth_singleton",
            "path_depth",
            format!(
                "ordered run has {} bp at <= {} path-step visit(s)",
                run.bp, options.low_depth_max_visits
            ),
            run.bp as f64,
            run.bp as f64 / options.min_low_depth_run_bp.max(1) as f64,
        );
    }

    let sparse_fraction = if options.max_sparse_coverage_path_fraction.is_finite() {
        options.max_sparse_coverage_path_fraction.clamp(0.0, 1.0)
    } else {
        0.25
    };
    let sparse_runs = ordered_runs(parsed, |node| {
        if parsed.paths.is_empty() {
            return false;
        }
        let distinct = distinct_paths_by_node.get(node).copied().unwrap_or(0);
        distinct as f64 / parsed.paths.len() as f64 <= sparse_fraction
    });
    for run in sparse_runs {
        if run.bp < options.min_sparse_run_bp {
            continue;
        }
        emit_run_sites(
            parsed,
            bp_starts,
            sites,
            &run,
            "sparse_path_depth",
            "sparse_coverage_path_fraction",
            format!(
                "ordered run has {} bp covered by <= {:.3} of paths",
                run.bp, sparse_fraction
            ),
            run.bp as f64,
            run.bp as f64 / options.min_sparse_run_bp.max(1) as f64,
        );
    }
}

fn detect_loop_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    options: &DirtyRegionOptions,
    sites: &mut Vec<DirtySite>,
) {
    let direct_loop_nodes = parsed
        .links
        .iter()
        .filter(|link| link.from == link.to)
        .map(|link| link.from.clone())
        .collect::<BTreeSet<_>>();
    for node in &direct_loop_nodes {
        let mut emitted = false;
        for path in &parsed.paths {
            for (idx, step) in path.steps.iter().enumerate() {
                if &step.node == node {
                    emitted = true;
                    sites.push(path_site(
                        parsed,
                        bp_starts,
                        "self_loop",
                        "direct_self_loop_edges",
                        format!("node {node} has a direct self-loop link"),
                        1.0,
                        1.0,
                        path,
                        idx,
                        idx + 1,
                        vec![node.clone()],
                    ));
                }
            }
        }
        if !emitted {
            sites.push(graph_site(
                parsed,
                bp_starts,
                "self_loop",
                "direct_self_loop_edges",
                format!("node {node} has a direct self-loop link but no path visit"),
                1.0,
                1.0,
                vec![node.clone()],
            ));
        }
    }

    for path in &parsed.paths {
        for (idx, pair) in path.steps.windows(2).enumerate() {
            if pair[0].node == pair[1].node {
                sites.push(path_site(
                    parsed,
                    bp_starts,
                    "adjacent_repeat",
                    "adjacent_same_node_path_steps",
                    format!("path repeats node {} in adjacent steps", pair[0].node),
                    1.0,
                    1.0,
                    path,
                    idx,
                    idx + 2,
                    vec![pair[0].node.clone()],
                ));
            }
        }
    }

    if options.small_loop_max_steps < 3 {
        return;
    }
    let mut seen = BTreeSet::new();
    for path in &parsed.paths {
        for start in 0..path.steps.len() {
            let max_end = (start + options.small_loop_max_steps - 1).min(path.steps.len() - 1);
            for end in (start + 2)..=max_end {
                if path.steps[start].node != path.steps[end].node {
                    continue;
                }
                let nodes = path.steps[start..=end]
                    .iter()
                    .map(|step| step.node.clone())
                    .collect::<Vec<_>>();
                let key = format!("{}:{start}:{end}:{}", path.name, nodes.join(","));
                if !seen.insert(key) {
                    continue;
                }
                sites.push(path_site(
                    parsed,
                    bp_starts,
                    "small_loop",
                    "small_loop_path_revisit",
                    format!(
                        "path revisits node {} within {} steps",
                        path.steps[start].node,
                        end - start + 1
                    ),
                    (end - start + 1) as f64,
                    (end - start + 1) as f64 / options.small_loop_max_steps.max(1) as f64,
                    path,
                    start,
                    end + 1,
                    nodes,
                ));
            }
        }
    }
}

fn ordered_runs(
    parsed: &ParsedGfa,
    mut predicate: impl FnMut(&str) -> bool,
) -> Vec<OrderedNodeRun> {
    let mut runs = Vec::new();
    let mut current_nodes = Vec::new();
    let mut current_start = 0usize;
    let mut current_bp = 0usize;

    for (order, node) in parsed.ordered_nodes.iter().enumerate() {
        if predicate(node) {
            if current_nodes.is_empty() {
                current_start = order;
            }
            current_bp =
                current_bp.saturating_add(parsed.segment_lengths.get(node).copied().unwrap_or(0));
            current_nodes.push(node.clone());
        } else if !current_nodes.is_empty() {
            runs.push(OrderedNodeRun {
                end_order: order - 1,
                start_order: current_start,
                nodes: std::mem::take(&mut current_nodes),
                bp: current_bp,
            });
            current_bp = 0;
        }
    }
    if !current_nodes.is_empty() {
        runs.push(OrderedNodeRun {
            start_order: current_start,
            end_order: parsed.ordered_nodes.len().saturating_sub(1),
            nodes: current_nodes,
            bp: current_bp,
        });
    }
    runs
}

fn emit_run_sites(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    sites: &mut Vec<DirtySite>,
    run: &OrderedNodeRun,
    kind: &str,
    source_metric: &str,
    reason: String,
    metric_value: f64,
    severity: f64,
) {
    let run_nodes = run.nodes.iter().map(String::as_str).collect::<HashSet<_>>();
    let mut emitted = false;
    for path in &parsed.paths {
        let mut idx = 0usize;
        while idx < path.steps.len() {
            if !run_nodes.contains(path.steps[idx].node.as_str()) {
                idx += 1;
                continue;
            }
            let start = idx;
            while idx < path.steps.len() && run_nodes.contains(path.steps[idx].node.as_str()) {
                idx += 1;
            }
            emitted = true;
            let nodes = path.steps[start..idx]
                .iter()
                .map(|step| step.node.clone())
                .collect::<Vec<_>>();
            sites.push(path_site(
                parsed,
                bp_starts,
                kind,
                source_metric,
                reason.clone(),
                metric_value,
                severity,
                path,
                start,
                idx,
                nodes,
            ));
        }
    }
    if !emitted {
        sites.push(graph_site(
            parsed,
            bp_starts,
            kind,
            source_metric,
            format!(
                "{}; no path traversal found for ordered run {}-{}",
                reason, run.start_order, run.end_order
            ),
            metric_value,
            severity,
            run.nodes.clone(),
        ));
    }
}

#[allow(clippy::too_many_arguments)]
fn path_site(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    kind: &str,
    source_metric: &str,
    reason: String,
    metric_value: f64,
    severity: f64,
    path: &PathWalk,
    start_step: usize,
    end_step: usize,
    nodes: Vec<String>,
) -> DirtySite {
    let start_step = start_step.min(path.steps.len());
    let end_step = end_step.min(path.steps.len()).max(start_step);
    let path_start_bp = path
        .step_starts
        .get(start_step)
        .copied()
        .unwrap_or(path.length_bp);
    let path_end_bp = path.step_end_bp(end_step);
    let (graph_start_bp, graph_end_bp, start_order, end_order) =
        graph_span_for_nodes(parsed, bp_starts, &nodes);
    DirtySite {
        id: String::new(),
        kind: kind.to_string(),
        source_metric: source_metric.to_string(),
        reason,
        metric_value,
        severity,
        path_name: Some(path.name.clone()),
        path_start_bp: Some(path_start_bp),
        path_end_bp: Some(path_end_bp.max(path_start_bp)),
        path_start_step: Some(start_step),
        path_end_step: Some(end_step),
        path_length_bp: Some(path.length_bp),
        graph_start_bp,
        graph_end_bp,
        start_order,
        end_order,
        nodes: canonical_nodes(nodes),
    }
}

#[allow(clippy::too_many_arguments)]
fn graph_site(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    kind: &str,
    source_metric: &str,
    reason: String,
    metric_value: f64,
    severity: f64,
    nodes: Vec<String>,
) -> DirtySite {
    let (graph_start_bp, graph_end_bp, start_order, end_order) =
        graph_span_for_nodes(parsed, bp_starts, &nodes);
    DirtySite {
        id: String::new(),
        kind: kind.to_string(),
        source_metric: source_metric.to_string(),
        reason,
        metric_value,
        severity,
        path_name: None,
        path_start_bp: None,
        path_end_bp: None,
        path_start_step: None,
        path_end_step: None,
        path_length_bp: None,
        graph_start_bp,
        graph_end_bp,
        start_order,
        end_order,
        nodes: canonical_nodes(nodes),
    }
}

fn canonical_nodes(nodes: Vec<String>) -> Vec<String> {
    nodes
        .into_iter()
        .collect::<BTreeSet<_>>()
        .into_iter()
        .collect()
}

fn graph_span_for_nodes(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    nodes: &[String],
) -> (Option<usize>, Option<usize>, Option<usize>, Option<usize>) {
    let mut graph_start_bp = None::<usize>;
    let mut graph_end_bp = None::<usize>;
    let mut start_order = None::<usize>;
    let mut end_order = None::<usize>;
    for node in nodes {
        if let Some(&order) = parsed.order.get(node) {
            start_order = Some(start_order.map_or(order, |current| current.min(order)));
            end_order = Some(end_order.map_or(order, |current| current.max(order)));
        }
        if let Some(&start_bp) = bp_starts.get(node) {
            let end_bp =
                start_bp.saturating_add(parsed.segment_lengths.get(node).copied().unwrap_or(0));
            graph_start_bp = Some(graph_start_bp.map_or(start_bp, |current| current.min(start_bp)));
            graph_end_bp = Some(graph_end_bp.map_or(end_bp, |current| current.max(end_bp)));
        }
    }
    (graph_start_bp, graph_end_bp, start_order, end_order)
}

fn ordered_gap_bp(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    from: &str,
    to: &str,
) -> Option<(usize, usize, usize)> {
    let from_start = *bp_starts.get(from)?;
    let to_start = *bp_starts.get(to)?;
    let from_end = from_start.saturating_add(parsed.segment_lengths.get(from).copied()?);
    let to_end = to_start.saturating_add(parsed.segment_lengths.get(to).copied()?);
    if from_end <= to_start {
        Some((to_start.saturating_sub(from_end), from_end, to_start))
    } else if to_end <= from_start {
        Some((from_start.saturating_sub(to_end), to_end, from_start))
    } else {
        Some((0, from_start.min(to_start), from_end.max(to_end)))
    }
}

fn assign_site_ids(sites: &mut Vec<DirtySite>) {
    let mut seen = BTreeSet::new();
    sites.retain(|site| {
        let key = format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            site.kind,
            site.path_name.as_deref().unwrap_or(""),
            option_usize(site.path_start_bp),
            option_usize(site.path_end_bp),
            option_usize(site.graph_start_bp),
            option_usize(site.graph_end_bp),
            site.source_metric,
            site.nodes.join(","),
        );
        seen.insert(key)
    });
    sites.sort_by(|a, b| {
        a.path_name
            .cmp(&b.path_name)
            .then(a.path_start_bp.cmp(&b.path_start_bp))
            .then(a.graph_start_bp.cmp(&b.graph_start_bp))
            .then(a.kind.cmp(&b.kind))
    });
    for (idx, site) in sites.iter_mut().enumerate() {
        site.id = format!("site_{idx:04}");
    }
}

fn diagnostic_from_site(site: &DirtySite) -> DirtyDiagnostic {
    DirtyDiagnostic {
        id: site.id.clone(),
        kind: site.kind.clone(),
        source_metric: site.source_metric.clone(),
        metric_value: site.metric_value,
        severity: site.severity,
        message: site.reason.clone(),
        path_name: site.path_name.clone(),
        path_start_bp: site.path_start_bp,
        path_end_bp: site.path_end_bp,
        graph_start_bp: site.graph_start_bp,
        graph_end_bp: site.graph_end_bp,
    }
}

fn merge_dirty_sites(
    parsed: &ParsedGfa,
    sites: &[DirtySite],
    options: &DirtyRegionOptions,
) -> Vec<DirtyChunk> {
    let mut grouped = BTreeMap::<Option<String>, Vec<&DirtySite>>::new();
    for site in sites {
        grouped
            .entry(site.path_name.clone())
            .or_default()
            .push(site);
    }

    let path_lengths = parsed
        .paths
        .iter()
        .map(|path| (path.name.clone(), path.length_bp))
        .collect::<HashMap<_, _>>();
    let graph_total_bp = parsed
        .ordered_nodes
        .iter()
        .map(|node| parsed.segment_lengths.get(node).copied().unwrap_or(0))
        .sum::<usize>();
    let mut chunks = Vec::new();

    for (path_name, mut group) in grouped {
        if group.is_empty() {
            continue;
        }
        group.sort_by_key(|site| {
            (
                site.path_start_bp.or(site.graph_start_bp).unwrap_or(0),
                site.path_end_bp.or(site.graph_end_bp).unwrap_or(0),
                site.id.clone(),
            )
        });

        let mut current = ChunkAcc::from_site(group[0]);
        for site in group.into_iter().skip(1) {
            if can_merge(&current, site, options.merge_distance_bp) {
                current.add_site(site);
            } else {
                chunks.push(finalize_chunk(
                    current,
                    path_name
                        .as_ref()
                        .and_then(|name| path_lengths.get(name).copied()),
                    graph_total_bp,
                    options,
                ));
                current = ChunkAcc::from_site(site);
            }
        }
        chunks.push(finalize_chunk(
            current,
            path_name
                .as_ref()
                .and_then(|name| path_lengths.get(name).copied()),
            graph_total_bp,
            options,
        ));
    }

    chunks.sort_by(|a, b| {
        a.path_name
            .cmp(&b.path_name)
            .then(a.path_start_bp.cmp(&b.path_start_bp))
            .then(a.graph_start_bp.cmp(&b.graph_start_bp))
            .then(Reverse(a.site_count).cmp(&Reverse(b.site_count)))
    });
    for (idx, chunk) in chunks.iter_mut().enumerate() {
        chunk.id = format!("chunk_{idx:04}");
    }
    chunks
}

impl ChunkAcc {
    fn from_site(site: &DirtySite) -> Self {
        let mut acc = Self {
            path_name: site.path_name.clone(),
            raw_path_start_bp: site.path_start_bp.or(site.graph_start_bp),
            raw_path_end_bp: site.path_end_bp.or(site.graph_end_bp),
            path_length_bp: site.path_length_bp,
            graph_start_bp: site.graph_start_bp,
            graph_end_bp: site.graph_end_bp,
            start_order: site.start_order,
            end_order: site.end_order,
            path_start_step: site.path_start_step,
            path_end_step: site.path_end_step,
            site_ids: Vec::new(),
            kinds: BTreeSet::new(),
            nodes: BTreeSet::new(),
            max_severity: 0.0,
        };
        acc.add_site(site);
        acc
    }

    fn add_site(&mut self, site: &DirtySite) {
        self.raw_path_start_bp = min_option(
            self.raw_path_start_bp,
            site.path_start_bp.or(site.graph_start_bp),
        );
        self.raw_path_end_bp =
            max_option(self.raw_path_end_bp, site.path_end_bp.or(site.graph_end_bp));
        self.path_length_bp = self.path_length_bp.or(site.path_length_bp);
        self.graph_start_bp = min_option(self.graph_start_bp, site.graph_start_bp);
        self.graph_end_bp = max_option(self.graph_end_bp, site.graph_end_bp);
        self.start_order = min_option(self.start_order, site.start_order);
        self.end_order = max_option(self.end_order, site.end_order);
        self.path_start_step = min_option(self.path_start_step, site.path_start_step);
        self.path_end_step = max_option(self.path_end_step, site.path_end_step);
        self.site_ids.push(site.id.clone());
        self.kinds.insert(site.kind.clone());
        self.nodes.extend(site.nodes.iter().cloned());
        self.max_severity = self.max_severity.max(site.severity);
    }
}

fn can_merge(current: &ChunkAcc, site: &DirtySite, merge_distance_bp: usize) -> bool {
    if current.path_name != site.path_name {
        return false;
    }
    let Some(current_end) = current.raw_path_end_bp else {
        return false;
    };
    let Some(next_start) = site.path_start_bp.or(site.graph_start_bp) else {
        return false;
    };
    next_start <= current_end.saturating_add(merge_distance_bp)
}

fn finalize_chunk(
    acc: ChunkAcc,
    path_length_from_graph: Option<usize>,
    graph_total_bp: usize,
    options: &DirtyRegionOptions,
) -> DirtyChunk {
    let coordinate_limit = acc
        .path_length_bp
        .or(path_length_from_graph)
        .or_else(|| (graph_total_bp > 0).then_some(graph_total_bp));
    let mut start = acc.raw_path_start_bp.unwrap_or(0);
    let mut end = acc.raw_path_end_bp.unwrap_or(start);
    start = start.saturating_sub(options.flank_bp);
    end = end.saturating_add(options.flank_bp);
    if let Some(limit) = coordinate_limit {
        end = end.min(limit);
        start = start.min(end);
    }
    let mut capped_by_budget = false;
    if let Some(max_chunk_bp) = options.max_chunk_bp.filter(|value| *value > 0) {
        if end.saturating_sub(start) > max_chunk_bp {
            capped_by_budget = true;
            let raw_start = acc.raw_path_start_bp.unwrap_or(start);
            let raw_end = acc.raw_path_end_bp.unwrap_or(end);
            let center = raw_start.saturating_add(raw_end).saturating_div(2);
            start = center.saturating_sub(max_chunk_bp / 2);
            end = start.saturating_add(max_chunk_bp);
            if let Some(limit) = coordinate_limit {
                if end > limit {
                    end = limit;
                    start = end.saturating_sub(max_chunk_bp);
                }
            }
        }
    }

    DirtyChunk {
        id: String::new(),
        path_name: acc.path_name,
        path_start_bp: acc.path_length_bp.or(path_length_from_graph).map(|_| start),
        path_end_bp: acc.path_length_bp.or(path_length_from_graph).map(|_| end),
        path_length_bp: acc.path_length_bp.or(path_length_from_graph),
        graph_start_bp: acc.graph_start_bp,
        graph_end_bp: acc.graph_end_bp,
        start_order: acc.start_order,
        end_order: acc.end_order,
        site_count: acc.site_ids.len(),
        site_ids: acc.site_ids,
        kinds: acc.kinds.into_iter().collect(),
        nodes: acc.nodes.into_iter().collect(),
        max_severity: acc.max_severity,
        flank_bp: options.flank_bp,
        capped_by_budget,
    }
}

fn min_option<T: Ord>(a: Option<T>, b: Option<T>) -> Option<T> {
    match (a, b) {
        (Some(a), Some(b)) => Some(a.min(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn max_option<T: Ord>(a: Option<T>, b: Option<T>) -> Option<T> {
    match (a, b) {
        (Some(a), Some(b)) => Some(a.max(b)),
        (Some(a), None) => Some(a),
        (None, Some(b)) => Some(b),
        (None, None) => None,
    }
}

fn write_tsv_row(out: &mut String, fields: &[&str]) {
    let escaped = fields
        .iter()
        .map(|field| tsv_field(field))
        .collect::<Vec<_>>();
    let _ = writeln!(out, "{}", escaped.join("\t"));
}

fn tsv_field(raw: &str) -> String {
    raw.replace(['\t', '\n', '\r'], " ")
}

fn option_usize(value: Option<usize>) -> String {
    value.map(|value| value.to_string()).unwrap_or_default()
}

fn format_float(value: f64) -> String {
    if value.is_finite() {
        format!("{value:.6}")
    } else {
        String::new()
    }
}
