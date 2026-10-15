//! The partition-embedded graph export (the owner's ruling, 2026-11-06):
//! per-partition embedded pangenome graphs for the partitions covering
//! chrMT and chrI, exported as GFA with syncmer node ids CONTINUOUS with
//! the global syng interning (frame-independent identity) so the proven
//! likelihood instrument and truth index work unchanged.
//!
//! Construction path (stated exactly — the panel's own alignments only,
//! no new truth information, no caching of derived state):
//!   * the partition's MATERIAL is the panel's own alignment-induced
//!     membership, the completed `impg partition` BED
//!     (`-w 10000 -d 1000 --selection-mode haplotype
//!     --syng-min-chain-anchors 5 --syng-min-chain-fraction 0.5`,
//!     source commit 295bca9, 19,421 partitions, 100% of panel bp) —
//!     material is in the partition because it GENUINELY ALIGNS into the
//!     partition locality (the syng's own chain criteria). No
//!     content-homing, no window territories enter this stage.
//!   * the GRAPH is the impg graph machinery's syng-native range GFA:
//!     every member row becomes one `SyngGfaPathRange` and
//!     `write_range_gfa_with_mode` (the same writer `impg syng2gfa` and
//!     `impg query -o gfa --gfa-engine syng` use) materializes the raw
//!     syng overlap graph over those ranges. The frequency mask is
//!     DISABLED: no private splits, no clones, no frequency removal —
//!     every syncmer node keeps its global interning abs id, which is
//!     the ruling's key property (any masking re-interns nodes and
//!     breaks instrument/truth-index continuity). Inter-syncmer gap
//!     segments (offset > syncmer length) are spliced from real panel
//!     DNA via the AGC and interned by the writer AFTER
//!     `num_syncmer_nodes()`, so S names never alias the syncmer ids.
//!
//! Output per partition:
//!   `<out_dir>/partition<N>.gfa`       — H, S (syncmer + gap), L, P
//!                                         lines (GFA 1.0; P lines spell
//!                                         each member range with signed
//!                                         global node ids)
//!   `<out_dir>/partition<N>.gfa.map.json` — the node-id mapping back to
//!                                         the global syng: the
//!                                         n_syncmer_nodes boundary,
//!                                         the member list, and the
//!                                         segment/link/path counts.
//!
//! Usage:
//!   partition_graph_export <syng_prefix> <sequences.agc> <partitions_root> <out_dir> <partition_id>...
//!
//! Assessment-side only: reads the panel index and the partition BEDs,
//! writes GFA receipts; no threshold enters anything, no product file is
//! touched, no inference path is modified.

use impg::commands::syng2gfa::{
    write_range_gfa_with_mode, GfaVersion, SyngGfaFrequencyMask, SyngGfaMode, SyngGfaPathRange,
};
use impg::sequence_index::UnifiedSequenceIndex;
use impg::syng::{SyncmerParams, SyngIndex};
use std::env;
use std::fs;
use std::io::{BufWriter, Write};
use std::path::Path;

/// One member row of a partition BED: the alignment-induced membership.
/// The partition run's BEDs are source-forward (`impg partition` emits no
/// strand column), so the range strand is always `+`.
#[derive(Clone, Debug, PartialEq, Eq)]
struct MemberRow {
    path_name: String,
    start: u64,
    end: u64,
}

/// Parse a partition membership BED into member rows. Rows are
/// `name<TAB>start<TAB>end` (extra columns are tolerated and ignored);
/// blank lines are skipped; coordinates must be non-degenerate.
fn parse_member_bed(text: &str) -> Result<Vec<MemberRow>, String> {
    let mut rows = Vec::new();
    for (number, line) in text.lines().enumerate() {
        let line = line.trim_end_matches('\r');
        if line.trim().is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 3 {
            return Err(format!(
                "BED line {}: need >= 3 tab-separated fields, got {}",
                number + 1,
                fields.len()
            ));
        }
        let start: u64 = fields[1]
            .parse()
            .map_err(|_| format!("BED line {}: start is not an integer", number + 1))?;
        let end: u64 = fields[2]
            .parse()
            .map_err(|_| format!("BED line {}: end is not an integer", number + 1))?;
        if end <= start {
            return Err(format!(
                "BED line {}: degenerate interval [{start},{end})",
                number + 1
            ));
        }
        rows.push(MemberRow {
            path_name: fields[0].to_string(),
            start,
            end,
        });
    }
    if rows.is_empty() {
        return Err("partition BED has no member rows".to_string());
    }
    Ok(rows)
}

fn usage(program: &str) -> ! {
    eprintln!(
        "usage: {program} export <syng_prefix> <sequences.agc> <partitions_root> <out_dir> <partition_id>..."
    );
    std::process::exit(1);
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 6 || args[1] != "export" {
        usage(&args[0]);
    }
    let prefix = &args[2];
    let sequences = &args[3];
    let partitions_root = &args[4];
    let out_dir = &args[5];
    let partition_ids: Vec<u32> = args[6..]
        .iter()
        .map(|value| {
            value
                .parse()
                .unwrap_or_else(|_| usage(&args[0]))
        })
        .collect();
    if partition_ids.is_empty() {
        usage(&args[0]);
    }

    eprintln!("[partition-graph-export] loading syng index from {prefix}");
    let index = SyngIndex::load(prefix, SyncmerParams::default())
        .expect("failed to load syng index");
    let n_syncmer_nodes = index.num_syncmer_nodes() as u64;
    eprintln!(
        "[partition-graph-export] index loaded: {} syncmer nodes, {} paths",
        n_syncmer_nodes,
        index.name_map.name_to_path.len()
    );
    eprintln!("[partition-graph-export] loading sequences from {sequences}");
    let sequence_files = vec![sequences.clone()];
    let sequence_index =
        UnifiedSequenceIndex::from_files(&sequence_files).expect("failed to load sequences");

    fs::create_dir_all(out_dir).expect("create output directory");

    for &partition in &partition_ids {
        let bed_path = Path::new(partitions_root).join(format!("partition{partition}.bed"));
        let bed_text = fs::read_to_string(&bed_path)
            .unwrap_or_else(|_| panic!("read {}", bed_path.display()));
        let members = parse_member_bed(&bed_text)
            .unwrap_or_else(|error| panic!("parse {}: {error}", bed_path.display()));

        // One range per member row; the P line name mirrors the query
        // machinery's `name:start-end` convention so every path in the
        // GFA is traceable to its BED row.
        let mut ranges = Vec::with_capacity(members.len());
        for member in &members {
            let path_idx = index
                .name_map
                .name_to_path
                .get(&member.path_name)
                .unwrap_or_else(|| {
                    panic!(
                        "partition {partition}: member path '{}' not in the panel name map",
                        member.path_name
                    )
                });
            ranges.push(SyngGfaPathRange {
                path_idx: *path_idx as usize,
                name: format!("{}:{}-{}", member.path_name, member.start, member.end),
                start: member.start,
                end: member.end,
                strand: '+',
            });
        }

        let gfa_path = Path::new(out_dir).join(format!("partition{partition}.gfa"));
        let start = std::time::Instant::now();
        let file = fs::File::create(&gfa_path)
            .unwrap_or_else(|_| panic!("create {}", gfa_path.display()));
        let mut writer = BufWriter::with_capacity(4 * 1024 * 1024, file);
        let (segments, links, paths, skipped, gaps, gap_bp) =
            write_range_gfa_with_mode(
                &index,
                &ranges,
                &mut writer,
                GfaVersion::V1_0,
                Some(&sequence_index),
                SyngGfaMode::Raw,
                SyngGfaFrequencyMask::disabled(),
            )
            .unwrap_or_else(|error| panic!("render partition {partition}: {error}"));
        writer.flush().expect("flush gfa");
        drop(writer);

        // The node-id mapping sidecar: the essential statement is the
        // boundary between global syncmer ids and render-local gap
        // segments, plus the member list and the counts (the checker
        // re-derives everything else from the GFA itself).
        let map = serde_json::json!({
            "partition": partition,
            "construction": {
                "material": "impg partition membership BED (alignment-induced; syng chain criteria, commit 295bca9)",
                "graph": "raw syng range GFA via impg::commands::syng2gfa::write_range_gfa_with_mode",
                "mode": "raw",
                "frequency_mask": "disabled (no private splits, no clones, no removal)",
                "gap_fill": sequences,
            },
            "n_syncmer_nodes": n_syncmer_nodes,
            "node_id_mapping": {
                "syncmer_nodes": "S names in [1, n_syncmer_nodes] are the GLOBAL syng syncmer interning abs ids (identity, not a re-interning)",
                "gap_nodes": "S names > n_syncmer_nodes are render-local inter-syncmer gap segments (interned by the writer after n_syncmer_nodes; never alias syncmer ids)",
                "path_steps": "P lines spell signed global syncmer ids (sign = storage strand)",
            },
            "counts": {
                "members": members.len(),
                "segments": segments,
                "gap_segments": gaps,
                "gap_bp": gap_bp,
                "links": links,
                "paths": paths,
                "skipped_paths": skipped,
            },
            "members": members
                .iter()
                .map(|m| serde_json::json!({
                    "path_name": m.path_name,
                    "start": m.start,
                    "end": m.end,
                }))
                .collect::<Vec<_>>(),
        });
        let map_path = Path::new(out_dir).join(format!("partition{partition}.gfa.map.json"));
        fs::write(&map_path, serde_json::to_string_pretty(&map).expect("serialize map"))
            .unwrap_or_else(|_| panic!("write {}", map_path.display()));

        eprintln!(
            "[partition-graph-export] partition {partition}: {} members -> {} S ({} gap, {} gap bp), {} L, {} P ({} skipped) in {:.1}s -> {}",
            members.len(),
            segments,
            gaps,
            gap_bp,
            links,
            paths,
            skipped,
            start.elapsed().as_secs_f64(),
            gfa_path.display()
        );
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// The BED parser must accept the partition run's 3-column rows, skip
/// blank lines, and reject degenerate/malformed rows with a named error.
#[test]
fn member_bed_rows_parse_and_validate() {
    let text = "AAA#0#chrI\t0\t10006\n\nAAC#0#chrI\t9042\t15392\r\nBTE#3#block28_contig1\t40000\t50000\textra\tcolumns\n";
    let rows = parse_member_bed(text).expect("parse");
    assert_eq!(rows.len(), 3);
    assert_eq!(
        rows[0],
        MemberRow {
            path_name: "AAA#0#chrI".to_string(),
            start: 0,
            end: 10006
        }
    );
    assert_eq!(rows[2].path_name, "BTE#3#block28_contig1");
    assert_eq!((rows[2].start, rows[2].end), (40000, 50000));

    assert!(parse_member_bed("name\t100\t50\n")
        .unwrap_err()
        .contains("degenerate"));
    assert!(parse_member_bed("name\t0\n")
        .unwrap_err()
        .contains("need >= 3"));
    assert!(parse_member_bed("name\tx\t50\n")
        .unwrap_err()
        .contains("start is not an integer"));
    assert!(parse_member_bed("\n\n").unwrap_err().contains("no member rows"));
}
