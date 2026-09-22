use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;

#[derive(Clone, Debug)]
pub struct GraphReportOptions {
    pub top_n: usize,
    pub max_link_jump_frac: f64,
    pub max_link_jump_p99: usize,
    pub max_path_jump_p99: usize,
    pub min_largest_component_frac: f64,
    pub min_common_start_frac: f64,
    pub min_common_end_frac: f64,
    pub max_internal_tips: usize,
    pub warn_duplicate_sequence_frac: f64,
    pub repeat_context_max_minor: usize,
    pub repeat_context_min_dominance: f64,
    pub min_white_space_gap_bp: usize,
    pub min_white_space_region_support: usize,
    pub max_path_white_space_p99: usize,
    pub max_sparse_coverage_path_fraction: f64,
    pub high_coverage_min_path_fraction: f64,
    pub include_povu: bool,
    pub povu_reference_names: Vec<String>,
}

impl Default for GraphReportOptions {
    fn default() -> Self {
        Self {
            top_n: 10,
            max_link_jump_frac: 0.25,
            max_link_jump_p99: 5_000,
            max_path_jump_p99: 5_000,
            min_largest_component_frac: 0.98,
            min_common_start_frac: 0.70,
            min_common_end_frac: 0.70,
            max_internal_tips: 0,
            warn_duplicate_sequence_frac: 0.10,
            repeat_context_max_minor: 2,
            repeat_context_min_dominance: 0.80,
            min_white_space_gap_bp: 1_000,
            min_white_space_region_support: 1,
            max_path_white_space_p99: 5_000,
            max_sparse_coverage_path_fraction: 0.25,
            high_coverage_min_path_fraction: 0.50,
            include_povu: false,
            povu_reference_names: Vec::new(),
        }
    }
}

#[derive(Clone, Debug, Serialize)]
pub struct GraphReport {
    pub source: String,
    pub status: String,
    pub failures: Vec<String>,
    pub warnings: Vec<String>,
    pub white_space_gap_threshold_bp: usize,
    pub sparse_coverage_path_fraction_threshold: f64,
    pub metrics: GraphMetrics,
    pub top_long_links: Vec<LinkJump>,
    pub top_path_jumps: Vec<PathJump>,
    pub top_white_space_jumps: Vec<WhiteSpaceJump>,
    pub top_white_space_regions: Vec<WhiteSpaceRegion>,
    pub sparse_coverage_runs: Vec<SparseCoverageRun>,
    pub top_depth_nodes: Vec<DepthNode>,
    pub depth_runs: Vec<DepthRun>,
    pub local_repeat_contexts: Vec<LocalRepeatContext>,
    pub self_loops: crate::gfa_self_loops::SelfLoopReport,
    pub povu: Option<PovuArchitecture>,
    pub povu_error: Option<String>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct GraphMetrics {
    pub segments: usize,
    pub links: usize,
    pub paths: usize,
    pub path_steps: usize,
    pub total_segment_bp: usize,
    pub node_path_step_visits: usize,
    pub node_coverage_mean: f64,
    pub node_coverage_bp_weighted_mean: f64,
    pub node_coverage_p10: usize,
    pub node_coverage_median: usize,
    pub node_coverage_p90: usize,
    pub singleton_nodes: usize,
    pub singleton_bp: usize,
    pub high_coverage_threshold: usize,
    pub high_coverage_nodes: usize,
    pub high_coverage_bp: usize,
    pub node_coverage_histogram: Vec<NodeCoverageBucket>,
    pub components: usize,
    pub largest_component_nodes: usize,
    pub largest_component_frac: f64,
    pub tips: usize,
    pub internal_tips: usize,
    pub common_start: Option<NodeSupport>,
    pub common_end: Option<NodeSupport>,
    pub link_jump_p95: usize,
    pub link_jump_p99: usize,
    pub link_jump_max: usize,
    pub path_jump_p95: usize,
    pub path_jump_p99: usize,
    pub path_jump_max: usize,
    pub path_white_space_bp_p95: usize,
    pub path_white_space_bp_p99: usize,
    pub path_white_space_bp_max: usize,
    pub path_white_space_bp_total: u64,
    pub path_white_space_bp_mean: f64,
    pub path_white_space_bridges: usize,
    pub path_white_space_bridges_ge_threshold: usize,
    pub segment_occupancy_bp_fraction: f64,
    pub segment_white_space_bp_fraction: f64,
    pub segment_white_space_bp_total: u64,
    pub segment_occupancy_path_fraction_p05: f64,
    pub segment_occupancy_path_fraction_median: f64,
    pub segment_occupancy_path_fraction_p95: f64,
    pub sparse_coverage_segment_bp: usize,
    pub duplicate_sequence_groups: usize,
    pub duplicate_sequence_nodes: usize,
    pub duplicate_sequence_frac: f64,
    pub max_duplicate_count: usize,
    pub path_depth_min: usize,
    pub path_depth_median: usize,
    pub path_depth_p95: usize,
    pub path_depth_max: usize,
    pub reused_nodes: usize,
    pub local_repeat_context_nodes: usize,
    pub local_repeat_context_occurrences: usize,
    pub direct_self_loop_edges: usize,
    pub direct_self_loop_nodes: usize,
    pub adjacent_same_node_path_steps: usize,
    pub adjacent_same_step_path_steps: usize,
    pub self_loop_repeat_runs: usize,
    pub self_loop_max_repeat_run_len: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct NodeCoverageBucket {
    pub label: String,
    pub min: usize,
    pub max: Option<usize>,
    pub nodes: usize,
    pub bp: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct NodeSupport {
    pub node: String,
    pub count: usize,
    pub fraction: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LinkJump {
    pub from: String,
    pub to: String,
    pub jump: usize,
    pub path_support: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct PathJump {
    pub path: String,
    pub step: usize,
    pub from: String,
    pub to: String,
    pub jump: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct WhiteSpaceJump {
    pub path: String,
    pub step: usize,
    pub from: String,
    pub to: String,
    pub gap_bp: usize,
    pub from_order: usize,
    pub to_order: usize,
    pub gap_start_bp: usize,
    pub gap_end_bp: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct WhiteSpaceRegion {
    pub start_bp: usize,
    pub end_bp: usize,
    pub length_bp: usize,
    pub crossing_path_steps: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct SparseCoverageRun {
    pub start_order: usize,
    pub end_order: usize,
    pub nodes: usize,
    pub start_node: String,
    pub end_node: String,
    pub bp: usize,
    pub missing_path_bp: u64,
    pub min_path_fraction: f64,
    pub max_path_fraction: f64,
    pub mean_path_fraction: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct DepthNode {
    pub node: String,
    pub order: usize,
    pub sequence_len: usize,
    pub sequence_preview: String,
    pub total_occurrences: usize,
    pub distinct_paths: usize,
    pub max_occurrences_in_one_path: usize,
    pub max_path: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct DepthRun {
    pub start_order: usize,
    pub end_order: usize,
    pub nodes: usize,
    pub start_node: String,
    pub end_node: String,
    pub min_depth: usize,
    pub max_depth: usize,
    pub mean_depth: f64,
}

#[derive(Clone, Debug, Serialize)]
pub struct LocalRepeatContext {
    pub node: String,
    pub total_occurrences: usize,
    pub dominant_count: usize,
    pub minor_occurrences: usize,
    pub dominant_fraction: f64,
    pub dominant_left: String,
    pub dominant_right: String,
}

#[derive(Clone, Debug, Serialize)]
pub struct PovuArchitecture {
    pub reference_path: String,
    pub reference_path_index: usize,
    pub sites: usize,
    pub leaf_sites: usize,
    pub level_counts: Vec<LevelCount>,
    pub top_sites: Vec<PovuSiteSummary>,
}

#[derive(Clone, Debug, Serialize)]
pub struct LevelCount {
    pub level: usize,
    pub count: usize,
}

#[derive(Clone, Debug, Serialize)]
pub struct PovuSiteSummary {
    pub id: String,
    pub parent_id: Option<String>,
    pub level: usize,
    pub is_leaf: bool,
    pub reference_start_step: usize,
    pub reference_end_step: usize,
    pub reference_span_steps: usize,
    pub start: String,
    pub end: String,
}

#[derive(Clone, Debug, Default, Eq, PartialEq)]
pub struct FlubblePathEmbeddingStats {
    pub sites: usize,
    pub paths_added: usize,
    pub sites_without_traversal: usize,
}

#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum FlubblePathEmbeddingMode {
    Dominant,
    PerPath,
}

#[derive(Clone, Debug)]
struct Step {
    node: String,
    orientation: char,
}

#[derive(Clone, Debug)]
struct PathWalk {
    name: String,
    steps: Vec<Step>,
}

#[derive(Clone, Debug)]
struct ParsedGfa {
    order: HashMap<String, usize>,
    segment_lengths: HashMap<String, usize>,
    segment_sequences: HashMap<String, String>,
    sequence_counts: HashMap<String, usize>,
    links: Vec<(String, String)>,
    paths: Vec<PathWalk>,
}

pub fn describe_gfa(
    source: impl Into<String>,
    gfa_text: &str,
    options: &GraphReportOptions,
) -> io::Result<GraphReport> {
    let parsed = parse_gfa(gfa_text);
    let self_loops = crate::gfa_self_loops::diagnose_gfa_with_top_n(gfa_text, options.top_n)?;
    let metrics = compute_metrics(&parsed, options, &self_loops);
    let bp_starts = segment_bp_starts(&parsed);
    let path_support = path_pair_support(&parsed.paths);
    let top_long_links = top_link_jumps(&parsed, &path_support, options.top_n);
    let top_path_jumps = top_path_jumps(&parsed, options.top_n);
    let top_white_space_jumps = top_white_space_jumps(&parsed, &bp_starts, options.top_n);
    let top_white_space_regions =
        top_white_space_regions(&parsed, &bp_starts, options.top_n, options);
    let sparse_coverage_runs = sparse_coverage_runs(&parsed, options.top_n, options);
    let top_depth_nodes = top_depth_nodes(&parsed, options.top_n);
    let depth_runs = depth_runs(&parsed, options.top_n);
    let mut local_repeat_contexts = local_repeat_contexts(&parsed.paths, options);
    local_repeat_contexts.truncate(options.top_n);
    let (povu, povu_error) = if options.include_povu {
        match povu_architecture(gfa_text, &options.povu_reference_names, options.top_n) {
            Ok(summary) => (Some(summary), None),
            Err(err) => (None, Some(err.to_string())),
        }
    } else {
        (None, None)
    };

    let mut failures = Vec::new();
    let mut warnings = Vec::new();
    if metrics.components > 1 {
        failures.push("components>1".to_string());
    }
    if metrics.largest_component_frac < options.min_largest_component_frac {
        failures.push("largest_component_frac".to_string());
    }
    if metrics.internal_tips > options.max_internal_tips {
        failures.push(format!("internal_tips>{}", options.max_internal_tips));
    }
    if metrics
        .common_start
        .as_ref()
        .is_some_and(|s| s.fraction < options.min_common_start_frac)
    {
        failures.push("common_start_frac".to_string());
    }
    if metrics
        .common_end
        .as_ref()
        .is_some_and(|s| s.fraction < options.min_common_end_frac)
    {
        failures.push("common_end_frac".to_string());
    }
    if metrics.link_jump_p99 > options.max_link_jump_p99 {
        failures.push("link_jump_p99".to_string());
    }
    if metrics.path_jump_p99 > options.max_path_jump_p99 {
        failures.push("path_jump_p99".to_string());
    }
    if metrics.path_white_space_bp_p99 > options.max_path_white_space_p99 {
        failures.push("path_white_space_bp_p99".to_string());
    }
    if metrics.segments > 0
        && metrics.link_jump_max > (metrics.segments as f64 * options.max_link_jump_frac) as usize
    {
        failures.push("link_jump_max_frac".to_string());
    }
    if metrics.duplicate_sequence_frac > options.warn_duplicate_sequence_frac {
        warnings.push("duplicate_sequence_frac".to_string());
    }
    if metrics.local_repeat_context_occurrences > 0 {
        warnings.push("local_repeat_contexts".to_string());
    }
    if self_loops.direct_self_loop_edges > 0 {
        warnings.push("direct_self_loop_edges".to_string());
    }
    if self_loops.adjacent_same_step_path_steps > 0 {
        warnings.push("adjacent_same_step_path_repeats".to_string());
    }
    if povu_error.is_some() {
        warnings.push("povu_decomposition_failed".to_string());
    }

    let status = if failures.is_empty() {
        "PASS"
    } else {
        "REVIEW"
    }
    .to_string();

    Ok(GraphReport {
        source: source.into(),
        status,
        failures,
        warnings,
        white_space_gap_threshold_bp: options.min_white_space_gap_bp,
        sparse_coverage_path_fraction_threshold: options.max_sparse_coverage_path_fraction,
        metrics,
        top_long_links,
        top_path_jumps,
        top_white_space_jumps,
        top_white_space_regions,
        sparse_coverage_runs,
        top_depth_nodes,
        depth_runs,
        local_repeat_contexts,
        self_loops,
        povu,
        povu_error,
    })
}

pub fn format_report(report: &GraphReport, format: &str) -> io::Result<String> {
    match format {
        "markdown" | "md" => Ok(report.to_markdown()),
        "json" => serde_json::to_string_pretty(report)
            .map(|mut s| {
                s.push('\n');
                s
            })
            .map_err(|e| io::Error::other(format!("failed to serialize graph report: {e}"))),
        "tsv" => Ok(report.to_tsv()),
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported graph report format '{other}'"),
        )),
    }
}

impl GraphReport {
    pub fn to_markdown(&self) -> String {
        let m = &self.metrics;
        let mut out = String::new();
        out.push_str("# Graph Characterization\n\n");
        out.push_str(&format!("- Source: `{}`\n", self.source));
        out.push_str(&format!("- Status: `{}`\n", self.status));
        out.push_str(&format!(
            "- Size: {} segments, {} links, {} paths, {} path steps, {} segment bp\n",
            m.segments, m.links, m.paths, m.path_steps, m.total_segment_bp
        ));
        if !self.failures.is_empty() {
            out.push_str(&format!("- Review flags: `{}`\n", self.failures.join(", ")));
        }
        if !self.warnings.is_empty() {
            out.push_str(&format!("- Warnings: `{}`\n", self.warnings.join(", ")));
        }

        out.push_str("\n## Large-Scale Architecture\n\n");
        out.push_str(&format!(
            "- Path-depth profile across segment order: min {}, median {}, p95 {}, max {}; reused nodes above path count: {}\n",
            m.path_depth_min,
            m.path_depth_median,
            m.path_depth_p95,
            m.path_depth_max,
            m.reused_nodes
        ));
        out.push_str(&format!(
            "- Node path coverage: {} total path-step visits over {} nodes; mean {:.2}, bp-weighted mean {:.2}; p10/median/p90 {}/{}/{}\n",
            m.node_path_step_visits,
            m.segments,
            m.node_coverage_mean,
            m.node_coverage_bp_weighted_mean,
            m.node_coverage_p10,
            m.node_coverage_median,
            m.node_coverage_p90
        ));
        out.push_str(&format!(
            "- Singleton nodes (coverage=1): {} nodes, {} bp; high-coverage shared nodes (coverage>={}): {} nodes, {} bp\n",
            m.singleton_nodes,
            m.singleton_bp,
            m.high_coverage_threshold,
            m.high_coverage_nodes,
            m.high_coverage_bp
        ));
        if !m.node_coverage_histogram.is_empty() {
            let histogram = m
                .node_coverage_histogram
                .iter()
                .filter(|bucket| bucket.nodes > 0)
                .map(|bucket| format!("{}:{} nodes/{} bp", bucket.label, bucket.nodes, bucket.bp))
                .collect::<Vec<_>>()
                .join(", ");
            if !histogram.is_empty() {
                out.push_str(&format!("- Node coverage histogram: {histogram}\n"));
            }
        }
        out.push_str(&format!(
            "- Path-by-segment occupancy: {:.3} covered, {:.3} white-space cells ({} missing path-bp); node occupancy fractions p05 {:.3}, median {:.3}, p95 {:.3}\n",
            m.segment_occupancy_bp_fraction,
            m.segment_white_space_bp_fraction,
            m.segment_white_space_bp_total,
            m.segment_occupancy_path_fraction_p05,
            m.segment_occupancy_path_fraction_median,
            m.segment_occupancy_path_fraction_p95
        ));
        if !self.top_depth_nodes.is_empty() {
            out.push_str("- Highest path-depth nodes:\n");
            for node in &self.top_depth_nodes {
                out.push_str(&format!(
                    "  - `{}` order {}: total {}, paths {}, max/path {} in `{}`, seq_len {}, seq `{}`\n",
                    node.node,
                    node.order,
                    node.total_occurrences,
                    node.distinct_paths,
                    node.max_occurrences_in_one_path,
                    node.max_path,
                    node.sequence_len,
                    node.sequence_preview
                ));
            }
        }
        if !self.depth_runs.is_empty() {
            out.push_str("- Highest-depth order runs:\n");
            for run in &self.depth_runs {
                out.push_str(&format!(
                    "  - orders {}-{} ({} nodes, `{}` to `{}`): depth min {}, max {}, mean {:.2}\n",
                    run.start_order,
                    run.end_order,
                    run.nodes,
                    run.start_node,
                    run.end_node,
                    run.min_depth,
                    run.max_depth,
                    run.mean_depth
                ));
            }
        }
        if !self.sparse_coverage_runs.is_empty() {
            out.push_str(&format!(
                "- Largest sparse-coverage runs by missing path-bp (segment path fraction <= {:.3}):\n",
                self.sparse_coverage_path_fraction_threshold
            ));
            for run in &self.sparse_coverage_runs {
                out.push_str(&format!(
                    "  - orders {}-{} ({} nodes, {} bp, `{}` to `{}`): missing path-bp {}, path fraction min/mean/max {:.3}/{:.3}/{:.3}\n",
                    run.start_order,
                    run.end_order,
                    run.nodes,
                    run.bp,
                    run.start_node,
                    run.end_node,
                    run.missing_path_bp,
                    run.min_path_fraction,
                    run.mean_path_fraction,
                    run.max_path_fraction
                ));
            }
        }
        if let Some(povu) = &self.povu {
            out.push_str(&format!(
                "- POVU flubble decomposition on reference `{}`: {} sites, {} leaves\n",
                povu.reference_path, povu.sites, povu.leaf_sites
            ));
            if !povu.level_counts.is_empty() {
                let levels = povu
                    .level_counts
                    .iter()
                    .map(|lc| format!("{}:{}", lc.level, lc.count))
                    .collect::<Vec<_>>()
                    .join(", ");
                out.push_str(&format!("- POVU site levels: {levels}\n"));
            }
            if !povu.top_sites.is_empty() {
                out.push_str("- Largest POVU sites by reference-step span:\n");
                for site in &povu.top_sites {
                    out.push_str(&format!(
                        "  - `{}` level {} leaf {} span {} steps, parent `{}`, boundaries `{}` -> `{}`\n",
                        site.id,
                        site.level,
                        site.is_leaf,
                        site.reference_span_steps,
                        site.parent_id.as_deref().unwrap_or("."),
                        site.start,
                        site.end
                    ));
                }
            }
        } else if let Some(err) = &self.povu_error {
            out.push_str(&format!("- POVU flubble decomposition failed: {err}\n"));
        }

        out.push_str("\n## Topology\n\n");
        out.push_str(&format!(
            "- Components: {} (largest: {} nodes, {:.3})\n",
            m.components, m.largest_component_nodes, m.largest_component_frac
        ));
        out.push_str(&format!(
            "- Tips: {} total, {} internal/non-endpoint\n",
            m.tips, m.internal_tips
        ));
        if let Some(start) = &m.common_start {
            out.push_str(&format!(
                "- Most common path start: `{}` ({}/{}, {:.3})\n",
                start.node, start.count, m.paths, start.fraction
            ));
        }
        if let Some(end) = &m.common_end {
            out.push_str(&format!(
                "- Most common path end: `{}` ({}/{}, {:.3})\n",
                end.node, end.count, m.paths, end.fraction
            ));
        }

        out.push_str("\n## Layout And Adjacency\n\n");
        out.push_str(&format!(
            "- Link jumps in current segment order: p95 {}, p99 {}, max {}\n",
            m.link_jump_p95, m.link_jump_p99, m.link_jump_max
        ));
        out.push_str(&format!(
            "- Consecutive path-step jumps: p95 {}, p99 {}, max {}\n",
            m.path_jump_p95, m.path_jump_p99, m.path_jump_max
        ));
        out.push_str(&format!(
            "- Path white-space bridges in current segment order: {} nonzero, {} >= {} bp; total {} bp, mean {:.1}, p95 {}, p99 {}, max {}\n",
            m.path_white_space_bridges,
            m.path_white_space_bridges_ge_threshold,
            self.white_space_gap_threshold_bp,
            m.path_white_space_bp_total,
            m.path_white_space_bp_mean,
            m.path_white_space_bp_p95,
            m.path_white_space_bp_p99,
            m.path_white_space_bp_max
        ));
        if !self.top_long_links.is_empty() {
            out.push_str("- Top long links:\n");
            for link in &self.top_long_links {
                out.push_str(&format!(
                    "  - `{}` -> `{}`: jump {}, path_support {}\n",
                    link.from, link.to, link.jump, link.path_support
                ));
            }
        }
        if !self.top_white_space_jumps.is_empty() {
            out.push_str("- Top path white-space bridges:\n");
            for jump in &self.top_white_space_jumps {
                out.push_str(&format!(
                    "  - `{}` step {}: `{}` -> `{}` bridges {} bp (orders {} -> {}, gap bp {}-{})\n",
                    jump.path,
                    jump.step,
                    jump.from,
                    jump.to,
                    jump.gap_bp,
                    jump.from_order,
                    jump.to_order,
                    jump.gap_start_bp,
                    jump.gap_end_bp
                ));
            }
        }
        if !self.top_white_space_regions.is_empty() {
            out.push_str("- Most-crossed white-space intervals:\n");
            for region in &self.top_white_space_regions {
                out.push_str(&format!(
                    "  - bp {}-{} ({} bp): crossed by {} consecutive path step(s)\n",
                    region.start_bp, region.end_bp, region.length_bp, region.crossing_path_steps
                ));
            }
        }
        if !self.top_path_jumps.is_empty() {
            out.push_str("- Top path jumps:\n");
            for jump in &self.top_path_jumps {
                out.push_str(&format!(
                    "  - `{}` step {}: `{}` -> `{}` jump {}\n",
                    jump.path, jump.step, jump.from, jump.to, jump.jump
                ));
            }
        }

        out.push_str("\n## Repeat And Glue Signals\n\n");
        out.push_str(&format!(
            "- Duplicate segment sequences: {} groups, {} nodes, max copy count {}, node fraction {:.3}\n",
            m.duplicate_sequence_groups,
            m.duplicate_sequence_nodes,
            m.max_duplicate_count,
            m.duplicate_sequence_frac
        ));
        out.push_str(&format!(
            "- Rare repeated local contexts: {} nodes, {} occurrences\n",
            m.local_repeat_context_nodes, m.local_repeat_context_occurrences
        ));
        if !self.local_repeat_contexts.is_empty() {
            out.push_str("- Top rare repeated-context nodes:\n");
            for ctx in &self.local_repeat_contexts {
                out.push_str(&format!(
                    "  - `{}`: total {}, dominant {} ({:.3}) between `{}` and `{}`, minor {}\n",
                    ctx.node,
                    ctx.total_occurrences,
                    ctx.dominant_count,
                    ctx.dominant_fraction,
                    ctx.dominant_left,
                    ctx.dominant_right,
                    ctx.minor_occurrences
                ));
            }
        }
        out.push_str(&format!(
            "- Direct self-loop L edges: {} edge(s) on {} node(s); adjacent same-node path steps {}, same-signed repeat steps {}, repeat runs {}, max run {}\n",
            m.direct_self_loop_edges,
            m.direct_self_loop_nodes,
            m.adjacent_same_node_path_steps,
            m.adjacent_same_step_path_steps,
            m.self_loop_repeat_runs,
            m.self_loop_max_repeat_run_len
        ));
        if !self.self_loops.length_buckets.is_empty() {
            let buckets = self
                .self_loops
                .length_buckets
                .iter()
                .map(|bucket| {
                    format!(
                        "{}bp:{} nodes/{} edges",
                        bucket.sequence_len, bucket.nodes, bucket.direct_self_loop_edges
                    )
                })
                .collect::<Vec<_>>()
                .join(", ");
            out.push_str(&format!("- Direct self-loop node lengths: {buckets}\n"));
        }
        if !self.self_loops.top_nodes.is_empty() {
            out.push_str("- Top direct self-loop nodes:\n");
            for node in &self.self_loops.top_nodes {
                out.push_str(&format!(
                    "  - `{}` order {}: seq_len {}, seq `{}`, self-loop edges {} ({:?}), visits {}, paths {}, adjacent repeats same-node/same-step {}/{}, runs {}, max run {} in `{}`\n",
                    node.node,
                    node.order,
                    node.sequence_len,
                    node.sequence_preview,
                    node.direct_self_loop_edges,
                    node.link_orientations,
                    node.total_path_visits,
                    node.distinct_paths,
                    node.adjacent_same_node_path_steps,
                    node.adjacent_same_step_path_steps,
                    node.repeated_path_runs,
                    node.max_repeat_run_len,
                    node.max_repeat_run_path
                ));
            }
        }

        out.push_str("\n## Interpretation\n\n");
        for line in self.interpretation_lines() {
            out.push_str("- ");
            out.push_str(&line);
            out.push('\n');
        }
        out
    }

    pub fn to_tsv(&self) -> String {
        let m = &self.metrics;
        let common_start = m
            .common_start
            .as_ref()
            .map(|s| s.node.as_str())
            .unwrap_or("");
        let common_end = m.common_end.as_ref().map(|s| s.node.as_str()).unwrap_or("");
        let start_frac = m.common_start.as_ref().map(|s| s.fraction).unwrap_or(0.0);
        let end_frac = m.common_end.as_ref().map(|s| s.fraction).unwrap_or(0.0);
        let header = [
            "source",
            "status",
            "failures",
            "warnings",
            "segments",
            "links",
            "paths",
            "path_steps",
            "total_segment_bp",
            "node_path_step_visits",
            "node_coverage_mean",
            "node_coverage_bp_weighted_mean",
            "node_coverage_p10",
            "node_coverage_median",
            "node_coverage_p90",
            "singleton_nodes",
            "singleton_bp",
            "high_coverage_threshold",
            "high_coverage_nodes",
            "high_coverage_bp",
            "components",
            "largest_component_frac",
            "internal_tips",
            "common_start",
            "common_start_frac",
            "common_end",
            "common_end_frac",
            "link_jump_p99",
            "link_jump_max",
            "path_jump_p99",
            "path_jump_max",
            "path_white_space_bp_p99",
            "path_white_space_bp_max",
            "path_white_space_bridges_ge_threshold",
            "segment_occupancy_bp_fraction",
            "segment_white_space_bp_fraction",
            "segment_white_space_bp_total",
            "sparse_coverage_segment_bp",
            "path_depth_median",
            "path_depth_p95",
            "path_depth_max",
            "reused_nodes",
            "duplicate_sequence_frac",
            "local_repeat_context_nodes",
            "local_repeat_context_occurrences",
            "direct_self_loop_edges",
            "direct_self_loop_nodes",
            "adjacent_same_node_path_steps",
            "adjacent_same_step_path_steps",
            "self_loop_repeat_runs",
            "self_loop_max_repeat_run_len",
            "self_loop_length_buckets",
            "povu_sites",
            "povu_leaf_sites",
        ];
        let values = vec![
            self.source.clone(),
            self.status.clone(),
            self.failures.join(","),
            self.warnings.join(","),
            m.segments.to_string(),
            m.links.to_string(),
            m.paths.to_string(),
            m.path_steps.to_string(),
            m.total_segment_bp.to_string(),
            m.node_path_step_visits.to_string(),
            format!("{:.6}", m.node_coverage_mean),
            format!("{:.6}", m.node_coverage_bp_weighted_mean),
            m.node_coverage_p10.to_string(),
            m.node_coverage_median.to_string(),
            m.node_coverage_p90.to_string(),
            m.singleton_nodes.to_string(),
            m.singleton_bp.to_string(),
            m.high_coverage_threshold.to_string(),
            m.high_coverage_nodes.to_string(),
            m.high_coverage_bp.to_string(),
            m.components.to_string(),
            format!("{:.6}", m.largest_component_frac),
            m.internal_tips.to_string(),
            common_start.to_string(),
            format!("{:.6}", start_frac),
            common_end.to_string(),
            format!("{:.6}", end_frac),
            m.link_jump_p99.to_string(),
            m.link_jump_max.to_string(),
            m.path_jump_p99.to_string(),
            m.path_jump_max.to_string(),
            m.path_white_space_bp_p99.to_string(),
            m.path_white_space_bp_max.to_string(),
            m.path_white_space_bridges_ge_threshold.to_string(),
            format!("{:.6}", m.segment_occupancy_bp_fraction),
            format!("{:.6}", m.segment_white_space_bp_fraction),
            m.segment_white_space_bp_total.to_string(),
            m.sparse_coverage_segment_bp.to_string(),
            m.path_depth_median.to_string(),
            m.path_depth_p95.to_string(),
            m.path_depth_max.to_string(),
            m.reused_nodes.to_string(),
            format!("{:.6}", m.duplicate_sequence_frac),
            m.local_repeat_context_nodes.to_string(),
            m.local_repeat_context_occurrences.to_string(),
            m.direct_self_loop_edges.to_string(),
            m.direct_self_loop_nodes.to_string(),
            m.adjacent_same_node_path_steps.to_string(),
            m.adjacent_same_step_path_steps.to_string(),
            m.self_loop_repeat_runs.to_string(),
            m.self_loop_max_repeat_run_len.to_string(),
            self.self_loop_length_bucket_tsv(),
            self.povu.as_ref().map(|p| p.sites).unwrap_or(0).to_string(),
            self.povu
                .as_ref()
                .map(|p| p.leaf_sites)
                .unwrap_or(0)
                .to_string(),
        ];
        debug_assert_eq!(header.len(), values.len());
        format!("{}\n{}\n", header.join("\t"), values.join("\t"))
    }

    fn self_loop_length_bucket_tsv(&self) -> String {
        self.self_loops
            .length_buckets
            .iter()
            .map(|bucket| {
                format!(
                    "{}bp:{}nodes:{}edges",
                    bucket.sequence_len, bucket.nodes, bucket.direct_self_loop_edges
                )
            })
            .collect::<Vec<_>>()
            .join(",")
    }

    fn interpretation_lines(&self) -> Vec<String> {
        let m = &self.metrics;
        let mut lines = Vec::new();
        if m.segments == 0 {
            lines.push("The graph has no segment records; graph construction likely produced an empty result.".to_string());
            return lines;
        }
        if m.components == 1 {
            lines.push("The graph is connected in the undirected S/L topology.".to_string());
        } else {
            lines.push(format!(
                "The graph has {} connected components; disconnected debris or over-split intervals may be present.",
                m.components
            ));
        }
        if let (Some(start), Some(end)) = (&m.common_start, &m.common_end) {
            if start.fraction >= 0.70 && end.fraction >= 0.70 {
                lines.push("Most paths share common entry and exit nodes, consistent with a focused local locus.".to_string());
            } else {
                lines.push("Path starts or ends are dispersed; this can indicate contig breaks, incomplete interval merging, or uncollapsed duplicated boundary nodes.".to_string());
            }
        }
        if m.link_jump_max > 0 && m.segments > 0 && m.link_jump_max > m.segments / 4 {
            let support = self
                .top_long_links
                .first()
                .map(|l| l.path_support)
                .unwrap_or(0);
            if support > 0 {
                lines.push("At least one long sorted-order link is also supported by path adjacency, so it is not just a renderer artifact.".to_string());
            } else {
                lines.push("The largest long sorted-order links are not observed as consecutive path steps in the parsed paths; inspect sorting/layout separately from biological adjacency.".to_string());
            }
        } else {
            lines.push(
                "Current segment-order link jumps are bounded relative to graph size.".to_string(),
            );
        }
        if m.segment_white_space_bp_fraction > 0.50 {
            lines.push("More than half of the path-by-segment matrix is white space; in a homologous core locus this usually means underalignment, repeat/paralog splitting, or large rare insertions rather than a compact pangenome representation.".to_string());
        } else if m.segment_white_space_bp_fraction > 0.20 {
            lines.push("The path-by-segment matrix has substantial white space; inspect the sparse-coverage runs to distinguish real insertion alleles from repeat-driven underalignment.".to_string());
        }
        lines.push(format!(
            "Mean node coverage is {:.2} path-step visits per node; bp-weighted mean path-depth is {:.2} haplotype bases per stored segment base. Prefer the bp-weighted value when comparing graph condensation because it discounts many tiny high-depth nodes and emphasizes shared representation across longer sequence.",
            m.node_coverage_mean,
            m.node_coverage_bp_weighted_mean
        ));
        if m.singleton_nodes > 0 {
            lines.push(format!(
                "{} singleton node(s) covering {} bp have coverage=1; persistent singleton bp usually marks private sequence, underalignment, or unshared replacement fragments.",
                m.singleton_nodes, m.singleton_bp
            ));
        }
        if m.path_white_space_bp_p99 > self.white_space_gap_threshold_bp {
            lines.push("Many consecutive path steps bridge long sorted-order gaps; these are the horizontal white-space spans visible in path renders and are direct targets for additional local alignment/crushing.".to_string());
        }
        if m.local_repeat_context_occurrences > 0 {
            lines.push("Rare repeated local contexts remain; these are candidates for single-syncmer repeat-loop glue or local copy-specific boundary artifacts.".to_string());
        } else {
            lines.push("No rare repeated local-context signal was detected under the configured dominance threshold.".to_string());
        }
        if m.direct_self_loop_edges > 0 {
            lines.push(format!(
                "{} direct self-loop edge(s) are present; if they are dominated by short repeat-unit nodes with same-signed adjacent path repeats, explicit run-node normalization can remove the loop without changing path spellings.",
                m.direct_self_loop_edges
            ));
        }
        if m.reused_nodes > 0 {
            lines.push("Some nodes are traversed more times than there are paths, which is a direct signal of collapsed copy-reuse in the local architecture.".to_string());
        }
        if let Some(povu) = &self.povu {
            if povu.sites > 0 {
                lines.push(format!(
                    "POVU sees {} nested flubble/site records; the largest listed sites are the best first targets for deciding whether the graph has one central SV/CNV region or scattered artifacts.",
                    povu.sites
                ));
            } else {
                lines.push(
                    "POVU did not find variable flubble sites on the selected reference path."
                        .to_string(),
                );
            }
        }
        if m.duplicate_sequence_frac > 0.10 {
            lines.push("Many segment sequences are duplicated; this is expected in syng-derived graphs but should be checked if tips or long jumps are also elevated.".to_string());
        }
        lines
    }
}

fn parse_gfa(gfa_text: &str) -> ParsedGfa {
    let mut order = HashMap::new();
    let mut segment_lengths = HashMap::new();
    let mut segment_sequences = HashMap::new();
    let mut sequence_counts = HashMap::new();
    let mut links = Vec::new();
    let mut paths = Vec::new();

    for line in gfa_text.lines() {
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.is_empty() {
            continue;
        }
        match fields[0] {
            "S" if fields.len() >= 3 => {
                let node = fields[1].to_string();
                let next = order.len();
                order.entry(node.clone()).or_insert(next);
                let sequence = fields[2];
                segment_lengths.insert(node, sequence_len(sequence));
                segment_sequences.insert(fields[1].to_string(), sequence.to_string());
                *sequence_counts.entry(sequence.to_string()).or_insert(0) += 1;
            }
            "L" if fields.len() >= 5 => {
                links.push((fields[1].to_string(), fields[3].to_string()));
            }
            "P" if fields.len() >= 3 => {
                let steps = split_p_steps(fields[2]);
                if !steps.is_empty() {
                    paths.push(PathWalk {
                        name: fields[1].to_string(),
                        steps,
                    });
                }
            }
            "W" if fields.len() >= 7 => {
                let steps = split_w_steps(fields[6]);
                if !steps.is_empty() {
                    paths.push(PathWalk {
                        name: w_line_name(&fields),
                        steps,
                    });
                }
            }
            _ => {}
        }
    }

    ParsedGfa {
        order,
        segment_lengths,
        segment_sequences,
        sequence_counts,
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
                orientation,
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
        if let Some(orientation) = current_orientation {
            let node = &raw[current_start..idx];
            if !node.is_empty() {
                steps.push(Step {
                    node: node.to_string(),
                    orientation,
                });
            }
        }
        current_orientation = Some(if ch == '>' { '+' } else { '-' });
        current_start = idx + ch.len_utf8();
    }
    if let Some(orientation) = current_orientation {
        let node = &raw[current_start..];
        if !node.is_empty() {
            steps.push(Step {
                node: node.to_string(),
                orientation,
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

fn compute_metrics(
    parsed: &ParsedGfa,
    options: &GraphReportOptions,
    self_loops: &crate::gfa_self_loops::SelfLoopReport,
) -> GraphMetrics {
    let segments = parsed.order.len();
    let path_steps = parsed.paths.iter().map(|p| p.steps.len()).sum();
    let component_sizes = component_sizes(&parsed.order, &parsed.links);
    let largest_component_nodes = component_sizes.first().copied().unwrap_or(0);
    let largest_component_frac = if segments == 0 {
        0.0
    } else {
        largest_component_nodes as f64 / segments as f64
    };
    let mut degree: HashMap<&str, usize> = HashMap::new();
    for (from, to) in &parsed.links {
        *degree.entry(from.as_str()).or_insert(0) += 1;
        *degree.entry(to.as_str()).or_insert(0) += 1;
    }
    let (starts, ends) = endpoint_counts(&parsed.paths);
    let endpoint_nodes: HashSet<&str> = starts
        .keys()
        .chain(ends.keys())
        .map(|node| node.as_str())
        .collect();
    let mut tips = 0usize;
    let mut internal_tips = 0usize;
    for node in parsed.order.keys() {
        if degree.get(node.as_str()).copied().unwrap_or(0) <= 1 {
            tips += 1;
            if !endpoint_nodes.contains(node.as_str()) {
                internal_tips += 1;
            }
        }
    }
    let link_jumps = sorted_jumps(
        &parsed.order,
        parsed
            .links
            .iter()
            .map(|(from, to)| (from.as_str(), to.as_str())),
    );
    let path_pairs: Vec<(&str, &str)> = parsed
        .paths
        .iter()
        .flat_map(|path| {
            path.steps
                .windows(2)
                .map(|w| (w[0].node.as_str(), w[1].node.as_str()))
        })
        .collect();
    let path_jumps = sorted_jumps(&parsed.order, path_pairs.into_iter());
    let bp_starts = segment_bp_starts(parsed);
    let mut path_white_space_bps = path_white_space_gaps(parsed, &bp_starts);
    path_white_space_bps.sort_unstable();
    let path_white_space_bp_total_u128 = path_white_space_bps
        .iter()
        .map(|&gap| gap as u128)
        .sum::<u128>();
    let path_white_space_bp_total = path_white_space_bp_total_u128.min(u64::MAX as u128) as u64;
    let path_white_space_bridges = path_white_space_bps.iter().filter(|&&gap| gap > 0).count();
    let path_white_space_bridges_ge_threshold = path_white_space_bps
        .iter()
        .filter(|&&gap| gap >= options.min_white_space_gap_bp)
        .count();
    let path_white_space_bp_mean = if path_white_space_bps.is_empty() {
        0.0
    } else {
        path_white_space_bp_total_u128 as f64 / path_white_space_bps.len() as f64
    };
    let distinct_depths = distinct_path_depths(parsed);
    let mut occupancy_path_fracs = Vec::with_capacity(parsed.order.len());
    let mut observed_path_bp = 0u128;
    let mut sparse_coverage_segment_bp = 0usize;
    for node in parsed.order.keys() {
        let len = parsed.segment_lengths.get(node).copied().unwrap_or(0);
        let distinct_paths = distinct_depths.get(node.as_str()).copied().unwrap_or(0);
        observed_path_bp += len as u128 * distinct_paths as u128;
        let frac = if parsed.paths.is_empty() {
            0.0
        } else {
            distinct_paths as f64 / parsed.paths.len() as f64
        };
        occupancy_path_fracs.push(frac);
        if frac <= options.max_sparse_coverage_path_fraction {
            sparse_coverage_segment_bp += len;
        }
    }
    occupancy_path_fracs.sort_unstable_by(f64::total_cmp);
    let possible_path_bp = parsed.paths.len() as u128
        * parsed.segment_lengths.values().copied().sum::<usize>() as u128;
    let segment_white_space_bp = possible_path_bp.saturating_sub(observed_path_bp);
    let segment_occupancy_bp_fraction = if possible_path_bp == 0 {
        0.0
    } else {
        observed_path_bp as f64 / possible_path_bp as f64
    };
    let segment_white_space_bp_fraction = if possible_path_bp == 0 {
        0.0
    } else {
        segment_white_space_bp as f64 / possible_path_bp as f64
    };
    let common_start = top_support(&starts, parsed.paths.len());
    let common_end = top_support(&ends, parsed.paths.len());
    let duplicate_sequence_groups = parsed
        .sequence_counts
        .values()
        .filter(|&&count| count > 1)
        .count();
    let duplicate_sequence_nodes = parsed
        .sequence_counts
        .values()
        .filter(|&&count| count > 1)
        .sum();
    let duplicate_sequence_frac = if segments == 0 {
        0.0
    } else {
        duplicate_sequence_nodes as f64 / segments as f64
    };
    let max_duplicate_count = parsed.sequence_counts.values().copied().max().unwrap_or(0);
    let repeat_contexts = local_repeat_contexts(&parsed.paths, options);
    let local_repeat_context_occurrences = repeat_contexts
        .iter()
        .map(|ctx| ctx.minor_occurrences)
        .sum::<usize>();
    let mut path_depths = path_depths(parsed);
    path_depths.sort_unstable();
    let path_depth_min = path_depths.first().copied().unwrap_or(0);
    let path_depth_median = quantile(&path_depths, 0.50);
    let path_depth_p95 = quantile(&path_depths, 0.95);
    let path_depth_max = path_depths.last().copied().unwrap_or(0);
    let node_coverage = node_coverage_summary(parsed, options.high_coverage_min_path_fraction);
    let reused_nodes = if parsed.paths.is_empty() {
        0
    } else {
        path_depths
            .iter()
            .filter(|&&depth| depth > parsed.paths.len())
            .count()
    };

    GraphMetrics {
        segments,
        links: parsed.links.len(),
        paths: parsed.paths.len(),
        path_steps,
        total_segment_bp: parsed.segment_lengths.values().sum(),
        node_path_step_visits: node_coverage.total_visits,
        node_coverage_mean: node_coverage.mean,
        node_coverage_bp_weighted_mean: node_coverage.bp_weighted_mean,
        node_coverage_p10: node_coverage.p10,
        node_coverage_median: node_coverage.median,
        node_coverage_p90: node_coverage.p90,
        singleton_nodes: node_coverage.singleton_nodes,
        singleton_bp: node_coverage.singleton_bp,
        high_coverage_threshold: node_coverage.high_coverage_threshold,
        high_coverage_nodes: node_coverage.high_coverage_nodes,
        high_coverage_bp: node_coverage.high_coverage_bp,
        node_coverage_histogram: node_coverage.histogram,
        components: component_sizes.len(),
        largest_component_nodes,
        largest_component_frac,
        tips,
        internal_tips,
        common_start,
        common_end,
        link_jump_p95: quantile(&link_jumps, 0.95),
        link_jump_p99: quantile(&link_jumps, 0.99),
        link_jump_max: link_jumps.last().copied().unwrap_or(0),
        path_jump_p95: quantile(&path_jumps, 0.95),
        path_jump_p99: quantile(&path_jumps, 0.99),
        path_jump_max: path_jumps.last().copied().unwrap_or(0),
        path_white_space_bp_p95: quantile(&path_white_space_bps, 0.95),
        path_white_space_bp_p99: quantile(&path_white_space_bps, 0.99),
        path_white_space_bp_max: path_white_space_bps.last().copied().unwrap_or(0),
        path_white_space_bp_total,
        path_white_space_bp_mean,
        path_white_space_bridges,
        path_white_space_bridges_ge_threshold,
        segment_occupancy_bp_fraction,
        segment_white_space_bp_fraction,
        segment_white_space_bp_total: segment_white_space_bp.min(u64::MAX as u128) as u64,
        segment_occupancy_path_fraction_p05: quantile_f64(&occupancy_path_fracs, 0.05),
        segment_occupancy_path_fraction_median: quantile_f64(&occupancy_path_fracs, 0.50),
        segment_occupancy_path_fraction_p95: quantile_f64(&occupancy_path_fracs, 0.95),
        sparse_coverage_segment_bp,
        duplicate_sequence_groups,
        duplicate_sequence_nodes,
        duplicate_sequence_frac,
        max_duplicate_count,
        path_depth_min,
        path_depth_median,
        path_depth_p95,
        path_depth_max,
        reused_nodes,
        local_repeat_context_nodes: repeat_contexts.len(),
        local_repeat_context_occurrences,
        direct_self_loop_edges: self_loops.direct_self_loop_edges,
        direct_self_loop_nodes: self_loops.direct_self_loop_nodes,
        adjacent_same_node_path_steps: self_loops.adjacent_same_node_path_steps,
        adjacent_same_step_path_steps: self_loops.adjacent_same_step_path_steps,
        self_loop_repeat_runs: self_loops.repeated_path_runs,
        self_loop_max_repeat_run_len: self_loops.max_repeat_run_len,
    }
}

#[derive(Clone, Debug, Default)]
struct NodeCoverageSummary {
    total_visits: usize,
    mean: f64,
    bp_weighted_mean: f64,
    p10: usize,
    median: usize,
    p90: usize,
    singleton_nodes: usize,
    singleton_bp: usize,
    high_coverage_threshold: usize,
    high_coverage_nodes: usize,
    high_coverage_bp: usize,
    histogram: Vec<NodeCoverageBucket>,
}

fn node_coverage_summary(parsed: &ParsedGfa, high_min_path_fraction: f64) -> NodeCoverageSummary {
    let mut visits_by_node: HashMap<&str, usize> = HashMap::new();
    for path in &parsed.paths {
        for step in &path.steps {
            *visits_by_node.entry(step.node.as_str()).or_insert(0) += 1;
        }
    }

    let high_threshold = high_coverage_threshold(parsed.paths.len(), high_min_path_fraction);
    let mut coverages = Vec::with_capacity(parsed.order.len());
    let mut total_visits = 0usize;
    let mut total_bp = 0usize;
    let mut weighted_visit_bp = 0u128;
    let mut singleton_nodes = 0usize;
    let mut singleton_bp = 0usize;
    let mut high_coverage_nodes = 0usize;
    let mut high_coverage_bp = 0usize;
    let mut histogram = empty_node_coverage_histogram();

    for node in parsed.order.keys() {
        let coverage = visits_by_node.get(node.as_str()).copied().unwrap_or(0);
        let len = parsed.segment_lengths.get(node).copied().unwrap_or(0);
        coverages.push(coverage);
        total_visits = total_visits.saturating_add(coverage);
        total_bp = total_bp.saturating_add(len);
        weighted_visit_bp =
            weighted_visit_bp.saturating_add((coverage as u128).saturating_mul(len as u128));
        if coverage == 1 {
            singleton_nodes += 1;
            singleton_bp = singleton_bp.saturating_add(len);
        }
        if coverage >= high_threshold {
            high_coverage_nodes += 1;
            high_coverage_bp = high_coverage_bp.saturating_add(len);
        }
        add_to_node_coverage_histogram(&mut histogram, coverage, len);
    }

    coverages.sort_unstable();
    let mean = if coverages.is_empty() {
        0.0
    } else {
        total_visits as f64 / coverages.len() as f64
    };
    let bp_weighted_mean = if total_bp == 0 {
        0.0
    } else {
        weighted_visit_bp as f64 / total_bp as f64
    };

    NodeCoverageSummary {
        total_visits,
        mean,
        bp_weighted_mean,
        p10: quantile(&coverages, 0.10),
        median: quantile(&coverages, 0.50),
        p90: quantile(&coverages, 0.90),
        singleton_nodes,
        singleton_bp,
        high_coverage_threshold: high_threshold,
        high_coverage_nodes,
        high_coverage_bp,
        histogram,
    }
}

fn high_coverage_threshold(paths: usize, min_path_fraction: f64) -> usize {
    let fraction = if min_path_fraction.is_finite() {
        min_path_fraction.clamp(0.0, 1.0)
    } else {
        0.50
    };
    ((paths as f64 * fraction).ceil() as usize).max(2)
}

fn empty_node_coverage_histogram() -> Vec<NodeCoverageBucket> {
    [
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
    .map(|(label, min, max)| NodeCoverageBucket {
        label: label.to_string(),
        min,
        max,
        nodes: 0,
        bp: 0,
    })
    .collect()
}

fn add_to_node_coverage_histogram(
    histogram: &mut [NodeCoverageBucket],
    coverage: usize,
    bp: usize,
) {
    let idx = match coverage {
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
    if let Some(bucket) = histogram.get_mut(idx) {
        bucket.nodes += 1;
        bucket.bp = bucket.bp.saturating_add(bp);
    }
}

fn component_sizes(order: &HashMap<String, usize>, links: &[(String, String)]) -> Vec<usize> {
    let mut parent: Vec<usize> = (0..order.len()).collect();
    let mut size = vec![1usize; order.len()];

    fn find(parent: &mut [usize], mut x: usize) -> usize {
        while parent[x] != x {
            parent[x] = parent[parent[x]];
            x = parent[x];
        }
        x
    }

    for (from, to) in links {
        let Some(&a) = order.get(from) else {
            continue;
        };
        let Some(&b) = order.get(to) else {
            continue;
        };
        let mut ra = find(&mut parent, a);
        let mut rb = find(&mut parent, b);
        if ra == rb {
            continue;
        }
        if size[ra] < size[rb] {
            std::mem::swap(&mut ra, &mut rb);
        }
        parent[rb] = ra;
        size[ra] += size[rb];
    }

    let mut counts: HashMap<usize, usize> = HashMap::new();
    for idx in 0..order.len() {
        let root = find(&mut parent, idx);
        *counts.entry(root).or_insert(0) += 1;
    }
    let mut sizes: Vec<usize> = counts.into_values().collect();
    sizes.sort_unstable_by_key(|&n| Reverse(n));
    sizes
}

fn sorted_jumps<'a>(
    order: &HashMap<String, usize>,
    pairs: impl Iterator<Item = (&'a str, &'a str)>,
) -> Vec<usize> {
    let mut jumps: Vec<usize> = pairs
        .filter_map(|(a, b)| Some(order.get(a)?.abs_diff(*order.get(b)?)))
        .collect();
    jumps.sort_unstable();
    jumps
}

fn quantile(values: &[usize], frac: f64) -> usize {
    if values.is_empty() {
        return 0;
    }
    let idx = ((frac * values.len() as f64).ceil() as usize)
        .saturating_sub(1)
        .min(values.len() - 1);
    values[idx]
}

fn quantile_f64(values: &[f64], frac: f64) -> f64 {
    if values.is_empty() {
        return 0.0;
    }
    let idx = ((frac * values.len() as f64).ceil() as usize)
        .saturating_sub(1)
        .min(values.len() - 1);
    values[idx]
}

fn segment_bp_starts(parsed: &ParsedGfa) -> HashMap<String, usize> {
    let mut nodes_by_order = vec![String::new(); parsed.order.len()];
    for (node, &order) in &parsed.order {
        if order < nodes_by_order.len() {
            nodes_by_order[order] = node.clone();
        }
    }
    let mut starts = HashMap::with_capacity(parsed.order.len());
    let mut offset = 0usize;
    for node in nodes_by_order {
        starts.insert(node.clone(), offset);
        offset = offset.saturating_add(parsed.segment_lengths.get(&node).copied().unwrap_or(0));
    }
    starts
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

fn path_white_space_gaps(parsed: &ParsedGfa, bp_starts: &HashMap<String, usize>) -> Vec<usize> {
    parsed
        .paths
        .iter()
        .flat_map(|path| {
            path.steps.windows(2).filter_map(|pair| {
                ordered_gap_bp(parsed, bp_starts, &pair[0].node, &pair[1].node)
                    .map(|(gap, _, _)| gap)
            })
        })
        .collect()
}

fn endpoint_counts(paths: &[PathWalk]) -> (HashMap<String, usize>, HashMap<String, usize>) {
    let mut starts = HashMap::new();
    let mut ends = HashMap::new();
    for path in paths {
        if let Some(first) = path.steps.first() {
            *starts.entry(step_label(first)).or_insert(0) += 1;
        }
        if let Some(last) = path.steps.last() {
            *ends.entry(step_label(last)).or_insert(0) += 1;
        }
    }
    (starts, ends)
}

fn top_support(counts: &HashMap<String, usize>, total: usize) -> Option<NodeSupport> {
    let (node, count) = counts
        .iter()
        .max_by(|(node_a, count_a), (node_b, count_b)| {
            count_a.cmp(count_b).then(node_b.cmp(node_a))
        })?;
    Some(NodeSupport {
        node: node.clone(),
        count: *count,
        fraction: if total == 0 {
            0.0
        } else {
            *count as f64 / total as f64
        },
    })
}

fn path_pair_support(paths: &[PathWalk]) -> HashMap<(String, String), usize> {
    let mut support = HashMap::new();
    for path in paths {
        for pair in path.steps.windows(2) {
            let a = step_label(&pair[0]);
            let b = step_label(&pair[1]);
            *support.entry((a.clone(), b.clone())).or_insert(0) += 1;
            *support.entry((b, a)).or_insert(0) += 1;
        }
    }
    support
}

fn top_link_jumps(
    parsed: &ParsedGfa,
    path_support: &HashMap<(String, String), usize>,
    top_n: usize,
) -> Vec<LinkJump> {
    let mut jumps: Vec<LinkJump> = parsed
        .links
        .iter()
        .filter_map(|(from, to)| {
            let jump = parsed.order.get(from)?.abs_diff(*parsed.order.get(to)?);
            let support = oriented_link_support(path_support, from, to);
            Some(LinkJump {
                from: from.clone(),
                to: to.clone(),
                jump,
                path_support: support,
            })
        })
        .collect();
    jumps.sort_unstable_by_key(|link| Reverse(link.jump));
    jumps.truncate(top_n);
    jumps
}

fn oriented_link_support(
    path_support: &HashMap<(String, String), usize>,
    from: &str,
    to: &str,
) -> usize {
    ['+', '-']
        .into_iter()
        .flat_map(|from_strand| {
            ['+', '-']
                .into_iter()
                .map(move |to_strand| (format!("{from}{from_strand}"), format!("{to}{to_strand}")))
        })
        .filter_map(|key| path_support.get(&key).copied())
        .sum()
}

fn top_path_jumps(parsed: &ParsedGfa, top_n: usize) -> Vec<PathJump> {
    let mut jumps = Vec::new();
    for path in &parsed.paths {
        for (idx, pair) in path.steps.windows(2).enumerate() {
            let Some(&from_order) = parsed.order.get(&pair[0].node) else {
                continue;
            };
            let Some(&to_order) = parsed.order.get(&pair[1].node) else {
                continue;
            };
            jumps.push(PathJump {
                path: path.name.clone(),
                step: idx,
                from: step_label(&pair[0]),
                to: step_label(&pair[1]),
                jump: from_order.abs_diff(to_order),
            });
        }
    }
    jumps.sort_unstable_by_key(|jump| Reverse(jump.jump));
    jumps.truncate(top_n);
    jumps
}

fn top_white_space_jumps(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    top_n: usize,
) -> Vec<WhiteSpaceJump> {
    let mut jumps = Vec::new();
    for path in &parsed.paths {
        for (idx, pair) in path.steps.windows(2).enumerate() {
            let Some((gap_bp, gap_start_bp, gap_end_bp)) =
                ordered_gap_bp(parsed, bp_starts, &pair[0].node, &pair[1].node)
            else {
                continue;
            };
            if gap_bp == 0 {
                continue;
            }
            let Some(&from_order) = parsed.order.get(&pair[0].node) else {
                continue;
            };
            let Some(&to_order) = parsed.order.get(&pair[1].node) else {
                continue;
            };
            jumps.push(WhiteSpaceJump {
                path: path.name.clone(),
                step: idx,
                from: step_label(&pair[0]),
                to: step_label(&pair[1]),
                gap_bp,
                from_order,
                to_order,
                gap_start_bp,
                gap_end_bp,
            });
        }
    }
    jumps.sort_unstable_by_key(|jump| {
        Reverse((jump.gap_bp, jump.from_order.abs_diff(jump.to_order)))
    });
    jumps.truncate(top_n);
    jumps
}

fn top_white_space_regions(
    parsed: &ParsedGfa,
    bp_starts: &HashMap<String, usize>,
    top_n: usize,
    options: &GraphReportOptions,
) -> Vec<WhiteSpaceRegion> {
    let mut events: BTreeMap<usize, isize> = BTreeMap::new();
    for path in &parsed.paths {
        for pair in path.steps.windows(2) {
            let Some((gap_bp, gap_start_bp, gap_end_bp)) =
                ordered_gap_bp(parsed, bp_starts, &pair[0].node, &pair[1].node)
            else {
                continue;
            };
            if gap_bp < options.min_white_space_gap_bp || gap_start_bp >= gap_end_bp {
                continue;
            }
            *events.entry(gap_start_bp).or_default() += 1;
            *events.entry(gap_end_bp).or_default() -= 1;
        }
    }

    let mut regions = Vec::new();
    let mut current_support = 0isize;
    let mut previous_pos = None;
    for (pos, delta) in events {
        if let Some(start) = previous_pos {
            if pos > start && current_support >= options.min_white_space_region_support as isize {
                regions.push(WhiteSpaceRegion {
                    start_bp: start,
                    end_bp: pos,
                    length_bp: pos.saturating_sub(start),
                    crossing_path_steps: current_support as usize,
                });
            }
        }
        current_support += delta;
        previous_pos = Some(pos);
    }

    regions.sort_unstable_by_key(|region| {
        Reverse((
            region.crossing_path_steps,
            region.length_bp,
            region.start_bp,
        ))
    });
    regions.truncate(top_n);
    regions
}

fn top_depth_nodes(parsed: &ParsedGfa, top_n: usize) -> Vec<DepthNode> {
    #[derive(Default)]
    struct Acc {
        total: usize,
        distinct_paths: usize,
        max_in_path: usize,
        max_path: String,
    }

    let mut totals: HashMap<&str, Acc> = HashMap::new();
    for path in &parsed.paths {
        let mut counts: HashMap<&str, usize> = HashMap::new();
        for step in &path.steps {
            *counts.entry(step.node.as_str()).or_insert(0) += 1;
        }
        for (node, count) in counts {
            let acc = totals.entry(node).or_default();
            acc.total += count;
            acc.distinct_paths += 1;
            if count > acc.max_in_path {
                acc.max_in_path = count;
                acc.max_path = path.name.clone();
            }
        }
    }

    let mut nodes: Vec<DepthNode> = totals
        .into_iter()
        .filter_map(|(node, acc)| {
            let order = *parsed.order.get(node)?;
            let sequence = parsed
                .segment_sequences
                .get(node)
                .map(String::as_str)
                .unwrap_or("");
            Some(DepthNode {
                node: node.to_string(),
                order,
                sequence_len: sequence_len(sequence),
                sequence_preview: sequence_preview(sequence, 48),
                total_occurrences: acc.total,
                distinct_paths: acc.distinct_paths,
                max_occurrences_in_one_path: acc.max_in_path,
                max_path: acc.max_path,
            })
        })
        .collect();
    nodes.sort_unstable_by_key(|node| {
        Reverse((
            node.total_occurrences,
            node.max_occurrences_in_one_path,
            node.distinct_paths,
        ))
    });
    nodes.truncate(top_n);
    nodes
}

fn sequence_preview(sequence: &str, max_len: usize) -> String {
    if sequence.len() <= max_len {
        sequence.to_string()
    } else {
        format!("{}...", &sequence[..max_len])
    }
}

fn path_depths(parsed: &ParsedGfa) -> Vec<usize> {
    let mut depths: HashMap<&str, usize> = HashMap::new();
    for path in &parsed.paths {
        for step in &path.steps {
            *depths.entry(step.node.as_str()).or_insert(0) += 1;
        }
    }
    parsed
        .order
        .keys()
        .map(|node| depths.get(node.as_str()).copied().unwrap_or(0))
        .collect()
}

fn distinct_path_depths(parsed: &ParsedGfa) -> HashMap<&str, usize> {
    let mut depths: HashMap<&str, usize> = HashMap::new();
    for path in &parsed.paths {
        let mut seen_in_path: HashSet<&str> = HashSet::new();
        for step in &path.steps {
            if seen_in_path.insert(step.node.as_str()) {
                *depths.entry(step.node.as_str()).or_insert(0) += 1;
            }
        }
    }
    depths
}

fn depth_runs(parsed: &ParsedGfa, top_n: usize) -> Vec<DepthRun> {
    if parsed.order.is_empty() || parsed.paths.is_empty() {
        return Vec::new();
    }
    let mut nodes_by_order = vec![String::new(); parsed.order.len()];
    for (node, &order) in &parsed.order {
        nodes_by_order[order] = node.clone();
    }
    let mut depths_by_order = vec![0usize; parsed.order.len()];
    for path in &parsed.paths {
        for step in &path.steps {
            if let Some(&order) = parsed.order.get(&step.node) {
                depths_by_order[order] += 1;
            }
        }
    }
    let max_depth = depths_by_order.iter().copied().max().unwrap_or(0);
    if max_depth == 0 {
        return Vec::new();
    }
    let threshold = ((max_depth as f64 * 0.80).ceil() as usize)
        .max(parsed.paths.len().saturating_div(2))
        .max(2);
    let mut runs = Vec::new();
    let mut start = None;
    for (idx, &depth) in depths_by_order.iter().enumerate() {
        if depth >= threshold {
            start.get_or_insert(idx);
        } else if let Some(run_start) = start.take() {
            runs.push(make_depth_run(
                run_start,
                idx.saturating_sub(1),
                &nodes_by_order,
                &depths_by_order,
            ));
        }
    }
    if let Some(run_start) = start {
        runs.push(make_depth_run(
            run_start,
            depths_by_order.len().saturating_sub(1),
            &nodes_by_order,
            &depths_by_order,
        ));
    }
    runs.sort_unstable_by_key(|run| {
        Reverse((run.max_depth, (run.mean_depth * 1000.0) as usize, run.nodes))
    });
    runs.truncate(top_n);
    runs
}

fn make_depth_run(
    start_order: usize,
    end_order: usize,
    nodes_by_order: &[String],
    depths_by_order: &[usize],
) -> DepthRun {
    let slice = &depths_by_order[start_order..=end_order];
    let min_depth = slice.iter().copied().min().unwrap_or(0);
    let max_depth = slice.iter().copied().max().unwrap_or(0);
    let sum = slice.iter().sum::<usize>();
    DepthRun {
        start_order,
        end_order,
        nodes: end_order.saturating_sub(start_order) + 1,
        start_node: nodes_by_order[start_order].clone(),
        end_node: nodes_by_order[end_order].clone(),
        min_depth,
        max_depth,
        mean_depth: if slice.is_empty() {
            0.0
        } else {
            sum as f64 / slice.len() as f64
        },
    }
}

fn sparse_coverage_runs(
    parsed: &ParsedGfa,
    top_n: usize,
    options: &GraphReportOptions,
) -> Vec<SparseCoverageRun> {
    if parsed.order.is_empty() || parsed.paths.is_empty() {
        return Vec::new();
    }
    let mut nodes_by_order = vec![String::new(); parsed.order.len()];
    for (node, &order) in &parsed.order {
        if order < nodes_by_order.len() {
            nodes_by_order[order] = node.clone();
        }
    }
    let distinct_depths = distinct_path_depths(parsed);
    let fractions_by_order = nodes_by_order
        .iter()
        .map(|node| {
            distinct_depths.get(node.as_str()).copied().unwrap_or(0) as f64
                / parsed.paths.len() as f64
        })
        .collect::<Vec<_>>();

    let mut runs = Vec::new();
    let mut start = None;
    for (idx, &fraction) in fractions_by_order.iter().enumerate() {
        if fraction <= options.max_sparse_coverage_path_fraction {
            start.get_or_insert(idx);
        } else if let Some(run_start) = start.take() {
            runs.push(make_sparse_coverage_run(
                run_start,
                idx.saturating_sub(1),
                &nodes_by_order,
                &fractions_by_order,
                parsed,
            ));
        }
    }
    if let Some(run_start) = start {
        runs.push(make_sparse_coverage_run(
            run_start,
            fractions_by_order.len().saturating_sub(1),
            &nodes_by_order,
            &fractions_by_order,
            parsed,
        ));
    }

    runs.sort_unstable_by_key(|run| {
        Reverse((run.missing_path_bp, run.bp as u64, run.nodes as u64))
    });
    runs.truncate(top_n);
    runs
}

fn make_sparse_coverage_run(
    start_order: usize,
    end_order: usize,
    nodes_by_order: &[String],
    fractions_by_order: &[f64],
    parsed: &ParsedGfa,
) -> SparseCoverageRun {
    let mut bp = 0usize;
    let mut missing_path_bp = 0u128;
    let mut observed_path_bp = 0u128;
    let mut min_path_fraction = f64::INFINITY;
    let mut max_path_fraction = 0.0f64;
    for idx in start_order..=end_order {
        let node = &nodes_by_order[idx];
        let len = parsed.segment_lengths.get(node).copied().unwrap_or(0);
        let fraction = fractions_by_order[idx];
        let distinct_paths = (fraction * parsed.paths.len() as f64).round() as usize;
        bp = bp.saturating_add(len);
        observed_path_bp += len as u128 * distinct_paths as u128;
        missing_path_bp += len as u128 * parsed.paths.len().saturating_sub(distinct_paths) as u128;
        min_path_fraction = min_path_fraction.min(fraction);
        max_path_fraction = max_path_fraction.max(fraction);
    }
    SparseCoverageRun {
        start_order,
        end_order,
        nodes: end_order.saturating_sub(start_order) + 1,
        start_node: nodes_by_order[start_order].clone(),
        end_node: nodes_by_order[end_order].clone(),
        bp,
        missing_path_bp: missing_path_bp.min(u64::MAX as u128) as u64,
        min_path_fraction: if min_path_fraction.is_finite() {
            min_path_fraction
        } else {
            0.0
        },
        max_path_fraction,
        mean_path_fraction: if bp == 0 || parsed.paths.is_empty() {
            0.0
        } else {
            observed_path_bp as f64 / (bp as f64 * parsed.paths.len() as f64)
        },
    }
}

fn local_repeat_contexts(
    paths: &[PathWalk],
    options: &GraphReportOptions,
) -> Vec<LocalRepeatContext> {
    if options.repeat_context_max_minor == 0 || paths.is_empty() {
        return Vec::new();
    }

    let mut context_counts: HashMap<String, HashMap<(String, String), usize>> = HashMap::new();
    for path in paths {
        for (idx, step) in path.steps.iter().enumerate() {
            let left = path
                .steps
                .get(idx.wrapping_sub(1))
                .filter(|_| idx > 0)
                .map(step_label)
                .unwrap_or_else(|| "^".to_string());
            let right = path
                .steps
                .get(idx + 1)
                .map(step_label)
                .unwrap_or_else(|| "$".to_string());
            *context_counts
                .entry(step_label(step))
                .or_default()
                .entry((left, right))
                .or_insert(0) += 1;
        }
    }

    let mut contexts = Vec::new();
    for (node, counts) in context_counts {
        if counts.len() <= 1 {
            continue;
        }
        let mut ranked: Vec<((String, String), usize)> = counts.into_iter().collect();
        ranked.sort_unstable_by(|a, b| b.1.cmp(&a.1).then_with(|| a.0.cmp(&b.0)));
        if ranked.len() > 1 && ranked[0].1 == ranked[1].1 {
            continue;
        }
        let total = ranked.iter().map(|(_, count)| *count).sum::<usize>();
        let (dominant_context, dominant_count) = ranked[0].clone();
        let minor_occurrences = total.saturating_sub(dominant_count);
        let dominant_fraction = if total == 0 {
            0.0
        } else {
            dominant_count as f64 / total as f64
        };
        if minor_occurrences <= options.repeat_context_max_minor
            && dominant_fraction >= options.repeat_context_min_dominance
        {
            contexts.push(LocalRepeatContext {
                node,
                total_occurrences: total,
                dominant_count,
                minor_occurrences,
                dominant_fraction,
                dominant_left: dominant_context.0,
                dominant_right: dominant_context.1,
            });
        }
    }
    contexts.sort_unstable_by_key(|ctx| Reverse((ctx.minor_occurrences, ctx.total_occurrences)));
    contexts
}

fn step_label(step: &Step) -> String {
    let mut label = step.node.clone();
    label.push(step.orientation);
    label
}

fn povu_architecture(
    gfa_text: &str,
    reference_names: &[String],
    top_n: usize,
) -> io::Result<PovuArchitecture> {
    let graph = povu::NativeGfa::parse(gfa_text)
        .map_err(|err| io::Error::other(format!("POVU parse failed: {err}")))?;
    let decomposition = graph
        .decompose_flubbles(reference_names)
        .map_err(|err| io::Error::other(format!("POVU decomposition failed: {err}")))?;
    let mut level_counts = BTreeMap::<usize, usize>::new();
    let mut leaf_sites = 0usize;
    for site in &decomposition.sites {
        *level_counts.entry(site.level).or_default() += 1;
        if site.is_leaf {
            leaf_sites += 1;
        }
    }
    let mut top_sites: Vec<PovuSiteSummary> = decomposition
        .sites
        .iter()
        .map(|site| PovuSiteSummary {
            id: site.id.clone(),
            parent_id: site.parent_id.clone(),
            level: site.level,
            is_leaf: site.is_leaf,
            reference_start_step: site.reference_start_step,
            reference_end_step: site.reference_end_step,
            reference_span_steps: site
                .reference_end_step
                .saturating_sub(site.reference_start_step),
            start: site.start.token(),
            end: site.end.token(),
        })
        .collect();
    top_sites.sort_unstable_by_key(|site| {
        Reverse((
            site.reference_span_steps,
            site.reference_start_step,
            site.reference_end_step,
        ))
    });
    top_sites.truncate(top_n);

    Ok(PovuArchitecture {
        reference_path: decomposition.reference_path,
        reference_path_index: decomposition.reference_path_index,
        sites: decomposition.sites.len(),
        leaf_sites,
        level_counts: level_counts
            .into_iter()
            .map(|(level, count)| LevelCount { level, count })
            .collect(),
        top_sites,
    })
}

pub fn gfa_with_dominant_flubble_paths(
    gfa_text: &str,
    reference_names: &[String],
) -> io::Result<(String, FlubblePathEmbeddingStats)> {
    gfa_with_flubble_paths(
        gfa_text,
        reference_names,
        FlubblePathEmbeddingMode::Dominant,
    )
}

pub fn gfa_with_flubble_paths(
    gfa_text: &str,
    reference_names: &[String],
    mode: FlubblePathEmbeddingMode,
) -> io::Result<(String, FlubblePathEmbeddingStats)> {
    let graph = povu::NativeGfa::parse(gfa_text)
        .map_err(|err| io::Error::other(format!("POVU parse failed: {err}")))?;
    let decomposition = graph
        .decompose_flubbles(reference_names)
        .map_err(|err| io::Error::other(format!("POVU decomposition failed: {err}")))?;

    let mut existing_path_names = graph
        .paths
        .iter()
        .map(|path| path.name.clone())
        .collect::<HashSet<_>>();
    let mut embedded = String::with_capacity(gfa_text.len() + decomposition.sites.len() * 128);
    embedded.push_str(gfa_text);
    if !embedded.ends_with('\n') {
        embedded.push('\n');
    }

    let mut stats = FlubblePathEmbeddingStats {
        sites: decomposition.sites.len(),
        ..Default::default()
    };

    for (rank, site) in decomposition.sites.iter().enumerate() {
        let span_steps = site
            .reference_end_step
            .saturating_sub(site.reference_start_step);
        let added = match mode {
            FlubblePathEmbeddingMode::Dominant => append_dominant_flubble_path(
                &mut embedded,
                &mut existing_path_names,
                &graph.paths,
                site,
                rank,
                span_steps,
            ),
            FlubblePathEmbeddingMode::PerPath => append_per_path_flubble_paths(
                &mut embedded,
                &mut existing_path_names,
                &graph.paths,
                site,
                rank,
                span_steps,
            ),
        };
        if added == 0 {
            stats.sites_without_traversal += 1;
        } else {
            stats.paths_added += added;
        }
    }

    Ok((embedded, stats))
}

#[derive(Clone, Debug)]
struct DominantFlubbleTraversal {
    steps: Vec<povu::native_gfa::Step>,
    support: usize,
    total_hits: usize,
    paths: usize,
}

fn dominant_flubble_traversal(
    paths: &[povu::native_gfa::PathRecord],
    site: &povu::native_gfa::FlubbleSite,
) -> Option<DominantFlubbleTraversal> {
    let mut counts = HashMap::<Vec<povu::native_gfa::Step>, usize>::new();
    let mut total_hits = 0usize;
    let mut paths_with_hits = 0usize;

    for path in paths {
        let mut path_hits = 0usize;
        for (start_idx, step) in path.steps.iter().enumerate() {
            if step != &site.start {
                continue;
            }
            let Some(end_offset) = path.steps[start_idx + 1..]
                .iter()
                .position(|candidate| candidate == &site.end)
            else {
                continue;
            };
            let end_idx = start_idx + 1 + end_offset;
            let traversal = path.steps[start_idx..=end_idx].to_vec();
            *counts.entry(traversal).or_default() += 1;
            total_hits += 1;
            path_hits += 1;
        }
        if path_hits > 0 {
            paths_with_hits += 1;
        }
    }

    counts
        .into_iter()
        .max_by(|(left_steps, left_count), (right_steps, right_count)| {
            left_count
                .cmp(right_count)
                .then_with(|| left_steps.len().cmp(&right_steps.len()))
                .then_with(|| left_steps.cmp(right_steps))
        })
        .map(|(steps, support)| DominantFlubbleTraversal {
            steps,
            support,
            total_hits,
            paths: paths_with_hits,
        })
}

fn append_dominant_flubble_path(
    embedded: &mut String,
    existing: &mut HashSet<String>,
    paths: &[povu::native_gfa::PathRecord],
    site: &povu::native_gfa::FlubbleSite,
    rank: usize,
    span_steps: usize,
) -> usize {
    let Some(dominant) = dominant_flubble_traversal(paths, site) else {
        return 0;
    };
    let name = unique_flubble_path_name(
        rank,
        site,
        span_steps,
        dominant.support,
        dominant.total_hits,
        dominant.paths,
        existing,
    );
    append_gfa_p_line(embedded, &name, &dominant.steps);
    1
}

fn append_per_path_flubble_paths(
    embedded: &mut String,
    existing: &mut HashSet<String>,
    paths: &[povu::native_gfa::PathRecord],
    site: &povu::native_gfa::FlubbleSite,
    rank: usize,
    span_steps: usize,
) -> usize {
    let mut added = 0usize;
    for (path_idx, path) in paths.iter().enumerate() {
        let mut path_hits = 0usize;
        for (start_idx, step) in path.steps.iter().enumerate() {
            if step != &site.start {
                continue;
            }
            let Some(end_offset) = path.steps[start_idx + 1..]
                .iter()
                .position(|candidate| candidate == &site.end)
            else {
                continue;
            };
            let end_idx = start_idx + 1 + end_offset;
            let name = unique_per_path_flubble_path_name(
                rank, site, span_steps, path_idx, path_hits, &path.name, existing,
            );
            append_gfa_p_line(embedded, &name, &path.steps[start_idx..=end_idx]);
            added += 1;
            path_hits += 1;
        }
    }
    added
}

fn unique_flubble_path_name(
    rank: usize,
    site: &povu::native_gfa::FlubbleSite,
    span_steps: usize,
    support: usize,
    total_hits: usize,
    paths: usize,
    existing: &mut HashSet<String>,
) -> String {
    let leaf = if site.is_leaf { "leaf" } else { "internal" };
    let id = sanitize_gfa_path_name(&site.id, 80);
    let base = format!(
        "_povu_flubble_{rank}_lv{}_{}_span{}_support{}_hits{}_paths{}_{}",
        site.level, leaf, span_steps, support, total_hits, paths, id
    );
    if existing.insert(base.clone()) {
        return base;
    }
    for suffix in 1usize.. {
        let candidate = format!("{base}_{suffix}");
        if existing.insert(candidate.clone()) {
            return candidate;
        }
    }
    unreachable!("unbounded suffix loop should always return")
}

fn unique_per_path_flubble_path_name(
    rank: usize,
    site: &povu::native_gfa::FlubbleSite,
    span_steps: usize,
    path_idx: usize,
    path_hit: usize,
    path_name: &str,
    existing: &mut HashSet<String>,
) -> String {
    let leaf = if site.is_leaf { "leaf" } else { "internal" };
    let id = sanitize_gfa_path_name(&site.id, 80);
    let source = sanitize_gfa_path_name(path_name, 96);
    let base = format!(
        "_povu_flubble_path_{rank}_lv{}_{}_span{}_path{}_hit{}_{}_{}",
        site.level, leaf, span_steps, path_idx, path_hit, source, id
    );
    if existing.insert(base.clone()) {
        return base;
    }
    for suffix in 1usize.. {
        let candidate = format!("{base}_{suffix}");
        if existing.insert(candidate.clone()) {
            return candidate;
        }
    }
    unreachable!("unbounded suffix loop should always return")
}

fn append_gfa_p_line(embedded: &mut String, name: &str, steps: &[povu::native_gfa::Step]) {
    let walk = steps
        .iter()
        .map(povu_step_to_p_line_token)
        .collect::<Vec<_>>()
        .join(",");
    embedded.push_str("P\t");
    embedded.push_str(name);
    embedded.push('\t');
    embedded.push_str(&walk);
    embedded.push_str("\t*\n");
}

fn sanitize_gfa_path_name(value: &str, max_len: usize) -> String {
    let mut out = String::new();
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '-' | '.' | '#') {
            out.push(ch);
        } else {
            out.push('_');
        }
        if out.len() >= max_len {
            break;
        }
    }
    if out.is_empty() {
        "site".to_string()
    } else {
        out
    }
}

fn povu_step_to_p_line_token(step: &povu::native_gfa::Step) -> String {
    let orientation = match step.strand {
        povu::native_gfa::Strand::Forward => '+',
        povu::native_gfa::Strand::Reverse => '-',
    };
    format!("{}{}", step.segment, orientation)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reports_basic_topology_and_long_links() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tTTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t1\t+\t4\t+\t0M
P\tp1\t1+,2+,3+,4+\t*
P\tp2\t1+,4+\t*
";
        let report = describe_gfa("tiny.gfa", gfa, &GraphReportOptions::default()).unwrap();
        assert_eq!(report.metrics.segments, 4);
        assert_eq!(report.metrics.paths, 2);
        assert_eq!(report.metrics.components, 1);
        assert_eq!(report.metrics.link_jump_max, 3);
        assert_eq!(report.metrics.path_white_space_bp_max, 8);
        assert!(report.metrics.segment_white_space_bp_fraction > 0.0);
        assert_eq!(report.top_depth_nodes[0].total_occurrences, 2);
        assert_eq!(report.top_long_links[0].from, "1");
        assert_eq!(report.top_long_links[0].to, "4");
        assert_eq!(report.top_long_links[0].path_support, 1);
    }

    #[test]
    fn reports_bp_white_space_and_sparse_coverage_runs() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCCCCCC
S\t3\tGGGG
S\t4\tTTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\twith_insertion\t1+,2+,4+\t*
P\twithout_insertion\t1+,3+,4+\t*
";
        let opts = GraphReportOptions {
            min_white_space_gap_bp: 1,
            max_sparse_coverage_path_fraction: 0.50,
            ..Default::default()
        };
        let report = describe_gfa("white.gfa", gfa, &opts).unwrap();
        assert_eq!(report.metrics.path_white_space_bp_max, 8);
        assert_eq!(report.metrics.path_white_space_bridges_ge_threshold, 2);
        assert_eq!(report.top_white_space_jumps[0].gap_bp, 8);
        assert!(!report.sparse_coverage_runs.is_empty());
        assert!(report.to_markdown().contains("Path-by-segment occupancy"));
    }

    #[test]
    fn reports_w_lines_and_local_repeat_contexts() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tTTTT
S\t5\tACAC
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t2\t+\t0M
L\t2\t+\t5\t+\t0M
W\ts\t0\tchr1\t0\t10\t>1>2>3>4>2>5
W\tt\t0\tchr1\t0\t10\t>1>2>3>4>5
";
        let opts = GraphReportOptions {
            repeat_context_max_minor: 1,
            repeat_context_min_dominance: 0.60,
            ..Default::default()
        };
        let report = describe_gfa("w.gfa", gfa, &opts).unwrap();
        assert_eq!(report.metrics.paths, 2);
        assert!(report.metrics.local_repeat_context_nodes >= 1);
        assert!(report
            .to_markdown()
            .contains("Rare repeated local contexts"));
    }

    #[test]
    fn embeds_dominant_flubble_traversal_paths() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAA
S\t2\tC
S\t3\tG
S\t4\tTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt1\t1+,3+,4+\t*
P\talt2\t1+,3+,4+\t*
";
        let (annotated, stats) =
            gfa_with_dominant_flubble_paths(gfa, &["ref".to_string()]).unwrap();
        assert_eq!(stats.sites, 1);
        assert_eq!(stats.paths_added, 1);
        assert_eq!(stats.sites_without_traversal, 0);
        assert!(
            annotated.contains("_povu_flubble_0_lv0_leaf_span2_support2_hits3_paths3"),
            "{annotated}"
        );
        assert!(
            annotated.contains("\t1+,3+,4+\t*"),
            "dominant alt traversal should be embedded:\n{annotated}"
        );
    }

    #[test]
    fn embeds_per_path_flubble_traversal_paths() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAA
S\t2\tC
S\t3\tG
S\t4\tTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tref\t1+,2+,4+\t*
P\talt1\t1+,3+,4+\t*
P\talt2\t1+,3+,4+\t*
";
        let (annotated, stats) =
            gfa_with_flubble_paths(gfa, &["ref".to_string()], FlubblePathEmbeddingMode::PerPath)
                .unwrap();
        assert_eq!(stats.sites, 1);
        assert_eq!(stats.paths_added, 3);
        assert_eq!(stats.sites_without_traversal, 0);
        assert!(
            annotated.contains("_povu_flubble_path_0_lv0_leaf_span2_path0_hit0_ref"),
            "{annotated}"
        );
        assert!(
            annotated.contains("_povu_flubble_path_0_lv0_leaf_span2_path1_hit0_alt1"),
            "{annotated}"
        );
        assert!(
            annotated.contains("_povu_flubble_path_0_lv0_leaf_span2_path2_hit0_alt2"),
            "{annotated}"
        );
    }

    #[test]
    fn reports_node_path_coverage_metrics_on_known_fixture() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCC
S\t3\tGGGGGG
S\t4\tTTTTTTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t1\t+\t4\t+\t0M
L\t4\t+\t3\t+\t0M
P\tp1\t1+,2+,3+\t*
P\tp2\t1+,2+,3+\t*
P\tp3\t1+,4+,3+\t*
W\ts\t0\tchr1\t0\t14\t>1>2>2>3
";
        let report = describe_gfa("coverage.gfa", gfa, &GraphReportOptions::default()).unwrap();
        let m = &report.metrics;
        assert_eq!(m.segments, 4);
        assert_eq!(m.paths, 4);
        assert_eq!(m.path_steps, 13);
        assert_eq!(m.node_path_step_visits, 13);
        assert!((m.node_coverage_mean - 3.25).abs() < 1e-9);
        assert!((m.node_coverage_bp_weighted_mean - 2.80).abs() < 1e-9);
        assert_eq!(m.node_coverage_p10, 1);
        assert_eq!(m.node_coverage_median, 4);
        assert_eq!(m.node_coverage_p90, 4);
        assert_eq!(m.singleton_nodes, 1);
        assert_eq!(m.singleton_bp, 8);
        assert_eq!(m.high_coverage_threshold, 2);
        assert_eq!(m.high_coverage_nodes, 3);
        assert_eq!(m.high_coverage_bp, 12);

        let singleton_bucket = m
            .node_coverage_histogram
            .iter()
            .find(|bucket| bucket.label == "1")
            .unwrap();
        assert_eq!(singleton_bucket.nodes, 1);
        assert_eq!(singleton_bucket.bp, 8);
        let shared_bucket = m
            .node_coverage_histogram
            .iter()
            .find(|bucket| bucket.label == "3-4")
            .unwrap();
        assert_eq!(shared_bucket.nodes, 3);
        assert_eq!(shared_bucket.bp, 12);
        assert!(report
            .to_markdown()
            .contains("Node path coverage: 13 total path-step visits"));
        assert!(report
            .to_markdown()
            .contains("bp-weighted mean path-depth is 2.80"));
    }

    #[test]
    fn formats_json_and_tsv() {
        let gfa = "H\tVN:Z:1.0\nS\t1\tA\nP\tp\t1+\t*\n";
        let report = describe_gfa("one.gfa", gfa, &GraphReportOptions::default()).unwrap();
        assert!(format_report(&report, "json")
            .unwrap()
            .contains("\"segments\": 1"));
        assert!(format_report(&report, "tsv")
            .unwrap()
            .starts_with("source\tstatus"));
    }
}
