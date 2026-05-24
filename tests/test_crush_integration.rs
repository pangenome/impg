//! End-to-end integration tests for the crush public API.
//!
//! These tests call ONLY the publicly-exported `resolve_gfa_bubbles` and
//! `path_sequences` functions — no direct calls to any inner resolution,
//! SyngGraph, or bubble functions.  They use real GFA fixture files from
//! `tests/test_data/crush/`.

use impg::resolution::{path_sequences, resolve_gfa_bubbles, ResolutionConfig};
use std::collections::HashMap;

fn seq_map(gfa: &str) -> HashMap<String, String> {
    path_sequences(gfa)
        .expect("path_sequences failed")
        .into_iter()
        .collect()
}

/// GFA 1.0 (P-line) crush: the public API resolves the insertion bubble and
/// preserves every path sequence exactly.
#[test]
fn crush_public_api_preserves_paths_p_line_gfa() {
    let gfa =
        include_str!("test_data/crush/small_insertion.gfa");
    let before = seq_map(gfa);
    assert!(!before.is_empty(), "fixture must have at least one path");

    let config = ResolutionConfig::default();
    let result = resolve_gfa_bubbles(gfa, &config)
        .expect("resolve_gfa_bubbles failed on P-line GFA");

    assert!(
        result.stats.resolved >= 1,
        "expected at least 1 resolved bubble, got {:?}",
        result.stats,
    );

    let after = seq_map(&result.gfa);
    assert_eq!(
        before.len(),
        after.len(),
        "path count must be preserved after crush",
    );
    for (name, before_seq) in &before {
        let after_seq = after
            .get(name)
            .unwrap_or_else(|| panic!("path '{}' missing from crushed output", name));
        assert_eq!(
            before_seq, after_seq,
            "path '{}' sequence changed after crush",
            name,
        );
    }
}

/// GFA 1.1 (W-line / walk) crush: the same public entry point dispatches to
/// the walk parser and resolves the insertion bubble while preserving sequences.
#[test]
fn crush_public_api_preserves_paths_w_line_gfa() {
    let gfa =
        include_str!("test_data/crush/small_insertion_walks.gfa");
    let before = seq_map(gfa);
    assert!(!before.is_empty(), "fixture must have at least one walk");

    let config = ResolutionConfig::default();
    let result = resolve_gfa_bubbles(gfa, &config)
        .expect("resolve_gfa_bubbles failed on W-line GFA");

    assert!(
        result.stats.resolved >= 1,
        "expected at least 1 resolved bubble in W-line GFA, got {:?}",
        result.stats,
    );

    let after = seq_map(&result.gfa);
    assert_eq!(
        before.len(),
        after.len(),
        "walk count must be preserved after crush",
    );
    for (name, before_seq) in &before {
        let after_seq = after
            .get(name)
            .unwrap_or_else(|| panic!("walk '{}' missing from crushed output", name));
        assert_eq!(
            before_seq, after_seq,
            "walk '{}' sequence changed after crush",
            name,
        );
    }
}

/// Both GFA variants produce the same per-path sequences after crush,
/// confirming the unified dispatch path is consistent.
#[test]
fn crush_public_api_p_and_w_line_produce_same_sequences() {
    let gfa_p = include_str!("test_data/crush/small_insertion.gfa");
    let gfa_w = include_str!("test_data/crush/small_insertion_walks.gfa");

    let config = ResolutionConfig::default();
    let result_p = resolve_gfa_bubbles(gfa_p, &config)
        .expect("resolve_gfa_bubbles failed on P-line GFA");
    let result_w = resolve_gfa_bubbles(gfa_w, &config)
        .expect("resolve_gfa_bubbles failed on W-line GFA");

    assert_eq!(
        result_p.stats.resolved, result_w.stats.resolved,
        "both GFA variants should resolve the same number of bubbles",
    );

    let seqs_p = seq_map(&result_p.gfa);
    let seqs_w = seq_map(&result_w.gfa);

    assert_eq!(
        seqs_p.len(),
        seqs_w.len(),
        "path count must match between P-line and W-line resolved output",
    );
    for (name, seq_p) in &seqs_p {
        let seq_w = seqs_w
            .get(name)
            .unwrap_or_else(|| panic!("path '{}' present in P-line but missing in W-line output", name));
        assert_eq!(
            seq_p, seq_w,
            "path '{}' sequence differs between P-line and W-line resolved output",
            name,
        );
    }
}
