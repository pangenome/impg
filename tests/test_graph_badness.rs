use impg::graph_badness::{analyze_gfa, format_dirty_report, DirtyRegionOptions};
use impg::graph_report::GraphReportOptions;

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
        max_chunk_bp: Some(20),
    }
}

fn has_site_kind(report: &impg::graph_badness::DirtyRegionReport, kind: &str) -> bool {
    report.dirty_sites.iter().any(|site| site.kind == kind)
}

fn assert_diagnostic_only(report: &impg::graph_badness::DirtyRegionReport) {
    assert_eq!(report.replacement_decision, "diagnostic-only");
    assert!(report.metrics_are_diagnostic_only);
}

#[test]
fn dirty_region_detector_leaves_clean_graph_without_chunks() {
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

    let report = analyze_gfa("clean", gfa, &detector_options()).unwrap();

    assert_diagnostic_only(&report);
    assert!(report.dirty_sites.is_empty(), "{:#?}", report.dirty_sites);
    assert!(report.candidate_chunks.is_empty());
    assert_eq!(report.metrics.segment_count, 3);
    assert_eq!(report.metrics.singleton_bp, 0);
    assert_eq!(report.metrics.spelled_path_bp, 36);
    assert_eq!(report.metrics.path_replay_compression_ratio, Some(3.0));
}

#[test]
fn dirty_region_detector_finds_underaligned_bubble_whitespace() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAAAA
S\t2\tCCCC
S\t3\tGGGG
S\t4\tTTTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr1:0-12\t1+,2+,4+\t*
P\tHG002#1#chr1:0-12\t1+,3+,4+\t*
P\tHG003#1#chr1:0-12\t1+,2+,4+\t*
";

    let report = analyze_gfa("bubble", gfa, &detector_options()).unwrap();

    assert_diagnostic_only(&report);
    assert!(has_site_kind(&report, "underaligned_white_space"));
    assert!(!report.diagnostics.is_empty());
    assert!(report
        .candidate_chunks
        .iter()
        .any(|chunk| chunk.path_name.as_deref() == Some("HG001#1#chr1:0-12")));
    assert!(report.metrics.white_space_proxy_bp_max >= 4);
}

#[test]
fn dirty_region_detector_finds_self_loop_microtangle() {
    let mut options = detector_options();
    options.min_white_space_gap_bp = 99;
    options.min_path_jump = 99;
    options.min_link_jump = 99;
    options.min_sparse_run_bp = 99;
    options.min_low_depth_run_bp = 99;
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tC
S\t3\tGG
L\t1\t+\t2\t+\t0M
L\t2\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tHG001#1#chr1:0-5\t1+,2+,2+,2+,3+\t*
P\tHG002#1#chr1:0-5\t1+,2+,3+\t*
";

    let report = analyze_gfa("loop", gfa, &options).unwrap();

    assert_diagnostic_only(&report);
    assert!(has_site_kind(&report, "self_loop"));
    assert!(has_site_kind(&report, "adjacent_repeat"));
    assert!(has_site_kind(&report, "small_loop"));
    assert_eq!(report.metrics.direct_self_loop_edges, 1);
    assert!(report.metrics.adjacent_same_node_path_steps >= 2);
    assert!(!report.candidate_chunks.is_empty());
}

#[test]
fn dirty_region_detector_finds_long_link_like_layout_jump() {
    let mut options = detector_options();
    options.min_sparse_run_bp = 99;
    options.min_low_depth_run_bp = 99;
    let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
S\t3\tG
S\t4\tT
S\t5\tA
S\t6\tC
L\t1\t+\t6\t+\t0M
P\tHG001#1#chr1:0-2\t1+,6+\t*
";

    let report = analyze_gfa("jump", gfa, &options).unwrap();

    assert_diagnostic_only(&report);
    assert!(has_site_kind(&report, "long_path_jump"));
    assert!(has_site_kind(&report, "long_link"));
    assert!(report.metrics.path_jump_max >= 5);
    assert!(report.metrics.link_jump_max >= 5);
    assert!(report
        .candidate_chunks
        .iter()
        .any(|chunk| chunk.path_name.as_deref() == Some("HG001#1#chr1:0-2")));
}

#[test]
fn dirty_region_detector_finds_low_depth_singleton_chunk_and_formats_reports() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAAA
S\t2\tCC
S\t3\tGGGGGG
S\t4\tTT
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr1:0-7\t1+,2+,4+\t*
P\tHG002#1#chr1:0-7\t1+,2+,4+\t*
P\tHG003#1#chr1:0-13\t1+,2+,3+,4+\t*
";

    let report = analyze_gfa("singleton", gfa, &detector_options()).unwrap();

    assert_diagnostic_only(&report);
    assert!(has_site_kind(&report, "low_depth_singleton"));
    assert!(has_site_kind(&report, "sparse_path_depth"));
    assert!(report.metrics.singleton_bp >= 6);
    assert!(report.metrics.low_depth_bp >= 6);
    assert!(report
        .candidate_chunks
        .iter()
        .any(|chunk| chunk.path_name.as_deref() == Some("HG003#1#chr1:0-13")));

    let json = format_dirty_report(&report, "json").unwrap();
    assert!(json.contains("\"candidate_chunks\""));
    assert!(json.contains("\"path_replay_compression_ratio\""));
    let tsv = format_dirty_report(&report, "tsv").unwrap();
    assert!(tsv.starts_with("record_type\tid\tsource"));
    assert!(tsv.contains("site\tsite_"));
    assert!(tsv.contains("chunk\tchunk_"));
    assert!(tsv.contains("path_replay_compression_ratio"));
}

#[test]
fn dirty_region_detector_merges_nearby_sites_with_flank_and_budget() {
    let mut options = detector_options();
    options.min_white_space_gap_bp = 99;
    options.min_path_jump = 99;
    options.min_link_jump = 99;
    options.min_sparse_run_bp = 99;
    options.min_low_depth_run_bp = 1;
    options.merge_distance_bp = 5;
    options.flank_bp = 1;
    options.max_chunk_bp = Some(5);
    let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tCC
S\t3\tG
S\t4\tTT
S\t5\tA
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
L\t4\t+\t5\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t5\t+\t0M
P\tHG001#1#chr1:0-3\t1+,3+,5+\t*
P\tHG002#1#chr1:0-3\t1+,3+,5+\t*
P\tHG003#1#chr1:0-7\t1+,2+,3+,4+,5+\t*
";

    let report = analyze_gfa("merge", gfa, &options).unwrap();
    let chunk = report
        .candidate_chunks
        .iter()
        .find(|chunk| chunk.path_name.as_deref() == Some("HG003#1#chr1:0-7"))
        .expect("singleton sites on HG003 should merge into one chunk");

    assert_diagnostic_only(&report);
    assert!(chunk.site_count >= 2, "{chunk:#?}");
    assert_eq!(chunk.flank_bp, 1);
    assert!(chunk.capped_by_budget, "{chunk:#?}");
    assert!(
        chunk.path_end_bp.unwrap() - chunk.path_start_bp.unwrap() <= 5,
        "{chunk:#?}"
    );
}
