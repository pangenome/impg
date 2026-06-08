use impg::graph_badness::DirtyChunk;
use impg::resolution::{
    path_sequences, resolve_gfa_dirty_chunks, LocalizedResolutionStatus, LocalizedResolverConfig,
    ResolutionConfig, ResolutionMethod,
};
use std::collections::BTreeMap;

fn path_map(gfa: &str) -> BTreeMap<String, String> {
    path_sequences(gfa)
        .expect("GFA path spellings should parse")
        .into_iter()
        .collect()
}

fn dirty_chunk(id: &str, path_name: &str, start: usize, end: usize, flank_bp: usize) -> DirtyChunk {
    DirtyChunk {
        id: id.to_string(),
        path_name: Some(path_name.to_string()),
        path_start_bp: Some(start.saturating_sub(flank_bp)),
        path_end_bp: Some(end.saturating_add(flank_bp)),
        core_path_start_bp: Some(start),
        core_path_end_bp: Some(end),
        path_length_bp: None,
        graph_start_bp: None,
        graph_end_bp: None,
        start_order: None,
        end_order: None,
        site_count: 1,
        site_ids: vec![format!("{id}_site")],
        kinds: vec!["test".to_string()],
        nodes: Vec::new(),
        max_severity: 1.0,
        flank_bp,
        capped_by_budget: false,
    }
}

fn local_config(method: ResolutionMethod, flank_bp: usize) -> LocalizedResolverConfig {
    let mut resolution = ResolutionConfig {
        method,
        replacement_flank_bp: flank_bp,
        polish_iterations: 0,
        ..ResolutionConfig::default()
    };
    resolution.auto_spoa_max_traversal_len = 1_000;
    LocalizedResolverConfig {
        resolution,
        max_chunk_bp: None,
        min_shared_paths: 2,
        skip_capped_chunks: true,
    }
}

fn assert_paths_preserved(before_gfa: &str, after_gfa: &str) {
    let before = path_map(before_gfa);
    let after = path_map(after_gfa);
    assert_eq!(
        before.keys().cloned().collect::<Vec<_>>(),
        after.keys().cloned().collect::<Vec<_>>(),
        "path names should be preserved exactly"
    );
    assert_eq!(before, after, "path spellings should be preserved exactly");
    assert!(after
        .keys()
        .all(|name| !name.starts_with("local_") && !name.starts_with("__impg")));
}

#[test]
fn localized_resolver_flanked_indel_trims_back_and_preserves_paths() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tG
S\t3\tGT
S\t4\tCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr1:0-5\t1+,2+,4+\t*
P\tHG002#1#chr1:0-6\t1+,3+,4+\t*
";
    let chunk = dirty_chunk("indel", "HG002#1#chr1:0-6", 2, 4, 2);
    let result = resolve_gfa_dirty_chunks(gfa, &[chunk], &local_config(ResolutionMethod::Poa, 2))
        .expect("localized POA indel replacement should run");

    assert_eq!(result.reports.len(), 1);
    let report = &result.reports[0];
    assert_eq!(report.status, LocalizedResolutionStatus::Applied);
    assert_eq!(report.method, Some(ResolutionMethod::Poa));
    assert_eq!(report.path_validation, "pass");
    assert!(report.trimmed_back, "{report:?}");
    assert_eq!(report.requested_flank_bp, 2);
    assert!(!report.quality_gate_used);
    assert_eq!(report.hard_rejection_gate, "exact_path_corruption");
    assert!(report.after.is_some(), "{report:?}");
    assert_paths_preserved(gfa, &result.gfa);
}

#[test]
fn localized_resolver_rebuilds_repeat_microtangle_without_path_name_drift() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tCA
S\t3\tGG
L\t1\t+\t2\t+\t0M
L\t2\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tHG001#1#chr2:0-8\t1+,2+,2+,3+\t*
P\tHG002#1#chr2:0-10\t1+,2+,2+,2+,3+\t*
";
    let chunk = dirty_chunk("repeat-microtangle", "HG002#1#chr2:0-10", 2, 8, 2);
    let result = resolve_gfa_dirty_chunks(gfa, &[chunk], &local_config(ResolutionMethod::Poa, 2))
        .expect("localized resolver should handle repeat microtangle");

    let report = &result.reports[0];
    assert_eq!(report.status, LocalizedResolutionStatus::Applied);
    assert_eq!(report.input_paths, 2);
    assert!(report.replacement_segments.unwrap_or(0) > 0);
    assert!(!report.quality_gate_used);
    assert_paths_preserved(gfa, &result.gfa);
}

#[test]
fn localized_resolver_handles_underaligned_bubble_as_diagnostic_not_gate() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tACGT
S\t2\tTTTTGGGG
S\t3\tTTTTAGGG
S\t4\tCCAA
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr3:0-16\t1+,2+,4+\t*
P\tHG002#1#chr3:0-16\t1+,3+,4+\t*
";
    let chunk = dirty_chunk("underaligned-bubble", "HG001#1#chr3:0-16", 4, 12, 4);
    let result = resolve_gfa_dirty_chunks(gfa, &[chunk], &local_config(ResolutionMethod::Poa, 4))
        .expect("localized resolver should rebuild underaligned bubble");

    let report = &result.reports[0];
    assert_eq!(report.status, LocalizedResolutionStatus::Applied);
    assert!(report.before.segments >= 4);
    assert!(report.after.is_some(), "{report:?}");
    assert_eq!(report.path_validation, "pass");
    assert!(!report.quality_gate_used);
    assert!(
        !report.reason.contains("quality"),
        "graph-quality metrics must not be a hidden rejection gate: {report:?}"
    );
    assert_paths_preserved(gfa, &result.gfa);
}

#[test]
fn localized_resolver_reports_budget_skip_without_quality_gate() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tG
S\t3\tGT
S\t4\tCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr4:0-5\t1+,2+,4+\t*
P\tHG002#1#chr4:0-6\t1+,3+,4+\t*
";
    let chunk = dirty_chunk("budget", "HG002#1#chr4:0-6", 2, 4, 2);
    let mut config = local_config(ResolutionMethod::Poa, 2);
    config.max_chunk_bp = Some(1);
    let result =
        resolve_gfa_dirty_chunks(gfa, &[chunk], &config).expect("budget skip should report");

    assert_eq!(result.gfa, gfa);
    let report = &result.reports[0];
    assert_eq!(report.status, LocalizedResolutionStatus::Skipped);
    assert!(report.reason.contains("budget"), "{report:?}");
    assert_eq!(report.path_validation, "not-run");
    assert!(!report.quality_gate_used);
    assert_eq!(report.hard_rejection_gate, "exact_path_corruption");
}

#[test]
fn localized_resolver_reports_sweepga_failure_without_spoa_fallback() {
    let gfa = "\
H\tVN:Z:1.0
S\t1\tAA
S\t2\tTTTT
S\t3\tTTTA
S\t4\tCC
L\t1\t+\t2\t+\t0M
L\t2\t+\t4\t+\t0M
L\t1\t+\t3\t+\t0M
L\t3\t+\t4\t+\t0M
P\tHG001#1#chr5:0-8\t1+,2+,4+\t*
P\tHG002#1#chr5:0-8\t1+,3+,4+\t*
";
    let chunk = dirty_chunk("sweepga-failure", "HG001#1#chr5:0-8", 2, 6, 2);
    let mut config = local_config(ResolutionMethod::Auto, 2);
    config.resolution.auto_spoa_max_traversal_len = 1;
    config.resolution.sweepga_aligner = "definitely-not-a-real-aligner".to_string();

    let result =
        resolve_gfa_dirty_chunks(gfa, &[chunk], &config).expect("failure should be reported");
    let report = &result.reports[0];
    assert_eq!(report.method, Some(ResolutionMethod::Sweepga));
    assert_eq!(report.status, LocalizedResolutionStatus::Failed);
    assert!(
        report.reason.contains("SweepGA") || report.reason.contains("align"),
        "{report:?}"
    );
    assert!(!report.quality_gate_used);
    assert_eq!(
        result.gfa, gfa,
        "failed SweepGA tier must not silently fall back to a different replacement"
    );
}
