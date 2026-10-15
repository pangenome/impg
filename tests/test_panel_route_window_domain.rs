//! The window-domain extension (owner-ruled Policy A + minimal universe
//! extension): a component window's candidate rows are the anchor group's
//! full BED plus every component-family row from any other group whose
//! interval numerically overlaps the window's coordinate range, each row
//! keeping its OWNING group (dual-role groups keep their existing axis
//! slot; pure-new groups get appended component-universe slots). The
//! extension's structural census (zero duplicate row identities, zero
//! overlapping rows per family sequence) is the no-duplication red line's
//! build-time gate.
#[path = "../examples/panel_route_diploid_search/mod.rs"]
mod search;

use search::genome_wide::{build_window_domain_extension, AxisFile, AxisInterval};
use std::{
    collections::BTreeSet,
    fs,
    path::PathBuf,
};

fn axis_interval(component: &str, start: u64, end: u64, group: &str) -> AxisInterval {
    AxisInterval {
        component: component.to_string(),
        start,
        end,
        group: group.to_string(),
        reference_occurrence: 0,
        reference_strand: "+".to_string(),
        orientations: Default::default(),
    }
}

fn write_bed(dir: &PathBuf, group: &str, rows: &[(&str, u64, u64)]) {
    let text = rows
        .iter()
        .map(|(source, start, end)| format!("{source}\t{start}\t{end}"))
        .collect::<Vec<_>>()
        .join("\n");
    fs::write(dir.join(format!("{group}.bed")), text).unwrap();
}

#[test]
fn window_domain_extension_policy_a_rows_and_ownership() {
    let temp = tempfile::tempdir().unwrap();
    let bed_directory = temp.path().to_path_buf();
    // Family: the two chrIII-suffix sources. chrV rows are skipped (foreign
    // family). Lengths only bound-validate the rows.
    let lanes = vec![
        ("S288C#0#chrIII".to_string(), 2_000u64),
        ("SK1#0#chrIII".to_string(), 2_000),
        ("S288C#0#chrV".to_string(), 1_000),
    ];
    let component_sources: BTreeSet<usize> = [0usize, 1].into_iter().collect();
    let axis = AxisFile {
        version: 1,
        coordinate_system: "declared-S288C#0-assembly-source-axis-v1".to_string(),
        intervals: vec![
            axis_interval("S288C#0#chrIII", 0, 500, "gA"),
            axis_interval("S288C#0#chrIII", 500, 1000, "gB"),
            // A foreign-component window: its group's rows must never enter
            // the chrIII extension.
            axis_interval("S288C#0#chrV", 0, 400, "gC"),
        ],
    };
    // gA: W0's anchor group. Its full BED also holds an SK1 row straddling
    // into W1's range (a dual-role group's row added at the OTHER window).
    write_bed(
        &bed_directory,
        "gA",
        &[
            ("S288C#0#chrIII", 0, 300),
            ("S288C#0#chrIII", 300, 500),
            ("SK1#0#chrIII", 400, 700),
        ],
    );
    // gB: W1's anchor.
    write_bed(&bed_directory, "gB", &[("S288C#0#chrIII", 500, 1000)]);
    // gD: a pure-new group whose SK1 row overlaps W0 only.
    write_bed(&bed_directory, "gD", &[("SK1#0#chrIII", 100, 400)]);
    // gC: the chrV axis group — a chrV row (foreign family, skipped).
    write_bed(&bed_directory, "gC", &[("S288C#0#chrV", 0, 100)]);
    // gZ: a group whose only chrIII row lies outside every window.
    write_bed(&bed_directory, "gZ", &[("SK1#0#chrIII", 1500, 1600)]);

    let extension = build_window_domain_extension(
        &axis,
        &bed_directory,
        &lanes,
        "S288C#0#chrIII",
        &component_sources,
    )
    .unwrap();

    // Ownership: gA anchors W0 (slot 0), gB anchors W1 (slot 1); gD is the
    // single pure-new group (appended slot component_loci + 0 = 2).
    assert_eq!(extension.component_loci, 2);
    let dual_role: BTreeSet<_> = extension.dual_role.iter().map(|(k, v)| (k.clone(), *v)).collect();
    let expected: BTreeSet<_> = [("gA".to_string(), 0usize), ("gB".to_string(), 1)].into_iter().collect();
    assert_eq!(dual_role, expected);
    assert_eq!(extension.pure_new.len(), 1);
    assert_eq!(extension.pure_new[0].0, "gD");
    assert_eq!(extension.pure_new[0].1, vec![(1usize, 100u64, 400u64)]);

    // W0's added rows: gD's SK1 row (owner 2), forward+reverse pairs.
    let w0 = &extension.added[0];
    assert_eq!(w0.len(), 2);
    assert!(w0.iter().all(|range| {
        range.source == 1 && range.start == 100 && range.end == 400 && range.partition == 2
    }));
    assert!(w0.iter().any(|range| !range.reverse));
    assert!(w0.iter().any(|range| range.reverse));
    // W1's added rows: gA's SK1 row (owner 0 — the dual-role group's
    // EXISTING axis slot), forward+reverse.
    let w1 = &extension.added[1];
    assert_eq!(w1.len(), 2);
    assert!(w1.iter().all(|range| {
        range.source == 1 && range.start == 400 && range.end == 700 && range.partition == 0
    }));
    // The anchor rows themselves are never re-added at their own window.
    assert!(w0.iter().all(|range| range.start != 0));
    assert!(w1.iter().all(|range| range.source != 0));

    // The no-duplication red line's structural census: clean tiling.
    assert_eq!(extension.duplicate_row_identities, 0);
    assert_eq!(extension.overlapping_row_pairs, 0);
}

#[test]
fn window_domain_extension_counts_tiling_violations() {
    let temp = tempfile::tempdir().unwrap();
    let bed_directory = temp.path().to_path_buf();
    let lanes = vec![
        ("S288C#0#chrIII".to_string(), 2_000u64),
        ("SK1#0#chrIII".to_string(), 2_000),
    ];
    let component_sources: BTreeSet<usize> = [0usize, 1].into_iter().collect();
    let axis = AxisFile {
        version: 1,
        coordinate_system: "declared-S288C#0-assembly-source-axis-v1".to_string(),
        intervals: vec![axis_interval("S288C#0#chrIII", 0, 1000, "gA")],
    };
    write_bed(&bed_directory, "gA", &[("S288C#0#chrIII", 0, 1000)]);
    // Two groups claiming the SAME row identity (a tiling violation the
    // census must count; the first sorted group keeps ownership).
    write_bed(&bed_directory, "gX", &[("SK1#0#chrIII", 100, 400)]);
    write_bed(&bed_directory, "gY", &[("SK1#0#chrIII", 100, 400)]);
    let extension = build_window_domain_extension(
        &axis,
        &bed_directory,
        &lanes,
        "S288C#0#chrIII",
        &component_sources,
    )
    .unwrap();
    assert_eq!(extension.duplicate_row_identities, 1);
    // The duplicate row is DROPPED from the second group (first sorted group
    // keeps ownership), so the kept rows are again disjoint: the adjacent
    // overlap census stays 0 while the identity census carries the
    // violation.
    assert_eq!(extension.overlapping_row_pairs, 0);
    // Only the first sorted group (gX) owns the row.
    assert_eq!(extension.pure_new.len(), 1);
    assert_eq!(extension.pure_new[0].0, "gX");
}
