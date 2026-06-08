use impg::graph_badness::DirtyRegionOptions;
use impg::graph_report::GraphReportOptions;
use impg::local_seed::{LocalSeedGraph, LocalSeedReport, LocalSequenceRecord, LocalSequenceSet};
use impg::localized_polish::{polish_seed_graph, LocalizedPolishConfig};
use impg::resolution::{
    path_sequences, LocalizedResolutionStatus, LocalizedResolverConfig, ResolutionConfig,
    ResolutionMethod,
};
use std::collections::BTreeMap;

fn path_map(gfa: &str) -> BTreeMap<String, String> {
    path_sequences(gfa)
        .expect("GFA path spellings should parse")
        .into_iter()
        .collect()
}

fn detector_options() -> DirtyRegionOptions {
    let graph_report = GraphReportOptions {
        top_n: 64,
        min_white_space_gap_bp: 1,
        max_sparse_coverage_path_fraction: 0.34,
        ..GraphReportOptions::default()
    };
    DirtyRegionOptions {
        graph_report,
        top_n: 64,
        min_white_space_gap_bp: 1,
        min_sparse_run_bp: 1,
        max_sparse_coverage_path_fraction: 0.34,
        min_path_jump: 3,
        min_link_jump: 3,
        low_depth_max_visits: 1,
        min_low_depth_run_bp: 1,
        small_loop_max_steps: 4,
        merge_distance_bp: 5,
        flank_bp: 1,
        max_chunk_bp: Some(100),
    }
}

fn local_config() -> LocalizedPolishConfig {
    let mut dirty_options = detector_options();
    dirty_options.flank_bp = 4;
    let mut resolver_resolution = ResolutionConfig {
        method: ResolutionMethod::Poa,
        replacement_flank_bp: 4,
        polish_iterations: 0,
        ..ResolutionConfig::default()
    };
    resolver_resolution.auto_spoa_max_traversal_len = 1_000;
    LocalizedPolishConfig {
        max_iterations: 3,
        max_chunks_per_iteration: 3,
        max_total_chunks: 8,
        max_total_chunk_bp: Some(1_000),
        max_runtime_secs: None,
        dirty_options,
        resolver: LocalizedResolverConfig {
            resolution: resolver_resolution,
            max_chunk_bp: Some(100),
            min_shared_paths: 2,
            skip_capped_chunks: true,
        },
        debug_dir: None,
    }
}

fn seed_graph(gfa: &str) -> LocalSeedGraph {
    LocalSeedGraph {
        gfa: gfa.to_string(),
        report: LocalSeedReport {
            seed_source: "synthetic-local-seed".to_string(),
            command_config: "synthetic".to_string(),
            sequence_count: path_map(gfa).len(),
            total_bp: path_map(gfa).values().map(String::len).sum(),
            raw_paf_records: None,
            path_count: path_map(gfa).len(),
            path_validation: "pass".to_string(),
            elapsed_secs: 0.0,
        },
    }
}

fn sequence_set(paths: &[(&str, &str)]) -> LocalSequenceSet {
    LocalSequenceSet::new(
        paths
            .iter()
            .map(|(path_name, sequence)| {
                let source_name = path_name
                    .rsplit_once(':')
                    .map(|(source, _)| source)
                    .unwrap_or(path_name)
                    .to_string();
                LocalSequenceRecord {
                    path_name: (*path_name).to_string(),
                    source_name,
                    source_start: 0,
                    source_end: sequence.len() as u64,
                    strand: '+',
                    sequence: sequence.as_bytes().to_vec(),
                }
            })
            .collect(),
    )
}

fn assert_paths_preserved(before_gfa: &str, after_gfa: &str) {
    let before = path_map(before_gfa);
    let after = path_map(after_gfa);
    assert_eq!(
        before, after,
        "localized polish must preserve every path name and spelling"
    );
    assert!(after
        .keys()
        .all(|name| !name.starts_with("local_") && !name.starts_with("__impg")));
}

#[test]
fn localized_polish_clean_graph_unchanged() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tHG001#1#chr1:0-12\t1+,2+,3+\t*
P\tHG002#1#chr1:0-12\t1+,2+,3+\t*
P\tHG003#1#chr1:0-12\t1+,2+,3+\t*
";
    let sequences = sequence_set(&[
        ("HG001#1#chr1:0-12", "AAAACCCCGGGG"),
        ("HG002#1#chr1:0-12", "AAAACCCCGGGG"),
        ("HG003#1#chr1:0-12", "AAAACCCCGGGG"),
    ]);

    let result =
        polish_seed_graph(&sequences, seed_graph(gfa), &local_config()).expect("polish runs");

    assert_eq!(result.gfa, gfa);
    assert_eq!(result.report.final_status, "converged-clean");
    assert_eq!(result.report.final_path_validation, "pass");
    assert!(result.report.metrics_are_diagnostic_only);
    assert_eq!(result.report.total_chunks_attempted, 0);
    assert_eq!(result.report.iterations[0].candidate_chunks, 0);
}

#[test]
fn localized_polish_underaligned_graph_improves_diagnostics_and_reports_dirty_chunks() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tCA
S\t3\tGG
L\t1\t+\t2\t+\t0M
L\t2\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tHG001#1#chr1:0-8\t1+,2+,2+,3+\t*
P\tHG002#1#chr1:0-10\t1+,2+,2+,2+,3+\t*
";
    let sequences = sequence_set(&[
        ("HG001#1#chr1:0-8", "AACACAGG"),
        ("HG002#1#chr1:0-10", "AACACACAGG"),
    ]);

    let result = polish_seed_graph(&sequences, seed_graph(gfa), &local_config())
        .expect("localized dirty-region polishing should run");

    assert_paths_preserved(gfa, &result.gfa);
    assert_eq!(result.report.hard_rejection_gate, "exact_path_corruption");
    assert!(result.report.metrics_are_diagnostic_only);
    assert!(
        result.report.total_chunks_attempted > 0,
        "{:#?}",
        result.report
    );
    assert!(
        result.report.total_chunks_applied > 0,
        "{:#?}",
        result.report
    );

    let first = &result.report.iterations[0];
    assert!(!first.selected_chunk_ids.is_empty(), "{first:#?}");
    assert!(first.candidate_chunks > 0, "{first:#?}");
    let after = first
        .metrics_after
        .as_ref()
        .expect("applied iteration should record after metrics");
    assert!(
        after.direct_self_loop_edges < first.metrics_before.direct_self_loop_edges
            || after.small_loop_sites < first.metrics_before.small_loop_sites
            || after.white_space_proxy_bp_total < first.metrics_before.white_space_proxy_bp_total
            || after.path_jump_max < first.metrics_before.path_jump_max
            || after.segment_count < first.metrics_before.segment_count,
        "expected at least one diagnostic metric to improve\nbefore={:#?}\nafter={:#?}",
        first.metrics_before,
        after
    );
    assert!(first
        .resolver_reports
        .iter()
        .all(|report| !report.quality_gate_used));
}

#[test]
fn localized_polish_path_corruption_is_hard_failure_before_replacement() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tHG001#1#chr1:0-12\t1+,2+,3+\t*
P\tHG002#1#chr1:0-12\t1+,2+,3+\t*
";
    let sequences = sequence_set(&[
        ("HG001#1#chr1:0-12", "AAAACCCCGGGG"),
        ("HG002#1#chr1:0-12", "AAAAGGGGGGGG"),
    ]);

    let err = polish_seed_graph(&sequences, seed_graph(gfa), &local_config()).unwrap_err();
    let message = err.to_string();
    assert!(
        message.contains("path validation failed") || message.contains("path corruption"),
        "{message}"
    );
}

#[test]
fn localized_polish_budget_skip_keeps_metrics_diagnostic() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tC
S\t3\tG
S\t4\tTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr1:0-5\t1+,2+,4+\t*
P\tHG002#1#chr1:0-5\t1+,3+,4+\t*
";
    let sequences = sequence_set(&[("HG001#1#chr1:0-5", "AACTT"), ("HG002#1#chr1:0-5", "AAGTT")]);
    let mut config = local_config();
    config.resolver.max_chunk_bp = Some(1);

    let result = polish_seed_graph(&sequences, seed_graph(gfa), &config)
        .expect("budget skip should report without failing path validation");

    assert_eq!(result.gfa, gfa);
    assert!(result.report.metrics_are_diagnostic_only);
    assert!(result.report.iterations[0].candidate_chunks > 0);
    assert_eq!(result.report.total_chunks_applied, 0);
    assert!(result.report.iterations[0]
        .resolver_reports
        .iter()
        .all(|report| {
            report.status == LocalizedResolutionStatus::Skipped && !report.quality_gate_used
        }));
}
