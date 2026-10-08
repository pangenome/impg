//! Standalone router/spelled phantom-span diagnosis (one-off, the
//! placement-gap ruling's deliverable 1): reproduces oracle_segment_spans
//! on the truth route's S288C chrIII piece, isolates the tagged record that
//! emits the flood-feature spans at path-9564 coordinates, and compares the
//! actual ACGT at those coordinates against the flood k-mer's verified BTE
//! home. Read-only; no model code is touched.
//!
//! Usage: router_spell_diag <panel.syng> <routes-dir>

use impg::genome_inference::{mem_records, panel_routes as routes, read_json};
use impg::sample_mem_bwt::{canonical, encode_walk};
use impg::syng::{SyncmerParams, SyngIndex};
use std::io;

const READ_LENGTH: usize = 150;

/// Copy of `examples/panel_route_mem_routed.rs::oracle_segment_spans`'s
/// event-run machinery (the spelled side's span builder) — the diagnostic
/// reproduces it exactly; the copy is pinned to this diagnosis.
fn event_boundaries(
    panel: &SyngIndex,
    sequence: &[u8],
    read_length: usize,
    start_lo: usize,
    start_hi: usize,
) -> io::Result<Vec<u64>> {
    let views = mem_records::raw_mem_anchor_positions(panel, sequence)?;
    let length = read_length as u64;
    let k = panel.syncmer_length_bp() as u64;
    let lo = start_lo as u64;
    let hi = start_hi as u64;
    let mut events = vec![lo, hi];
    for positions in &views {
        for &position in positions {
            for event in [
                (position + k).saturating_sub(length),
                position.saturating_add(1),
            ] {
                if event > lo && event < hi {
                    events.push(event);
                }
            }
        }
    }
    events.sort_unstable();
    events.dedup();
    Ok(events)
}

fn oracle_segment_spans(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<std::collections::BTreeMap<Vec<u64>, Vec<(u64, u64)>>> {
    let read_length = READ_LENGTH;
    let mut spans: std::collections::BTreeMap<Vec<u64>, Vec<(u64, u64)>> = Default::default();
    if sequence.len() < read_length {
        return Ok(spans);
    }
    let events = event_boundaries(panel, sequence, read_length, 0, sequence.len() - read_length + 1)?;
    let k = panel.syncmer_length_bp() as u64;
    for run in events.windows(2) {
        let start = run[0] as usize;
        for record in mem_records::tagged_mem_records(panel, &sequence[start..start + read_length])? {
            for i in 0..record.len() {
                for j in i..record.len() {
                    let encoded = encode_walk(&record[i..=j])?;
                    let feature = canonical(&encoded);
                    let lo = start as u64 + record[i].1;
                    let hi = start as u64 + record[j].1 + k;
                    spans.entry(feature).or_default().push((lo, hi));
                }
            }
        }
    }
    Ok(spans)
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let panel = SyngIndex::load(&args[1], SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    let routes_dir = std::path::Path::new(&args[2]);
    let graph: routes::Graph = read_json(&routes_dir.join("graph.json"))?;
    let lanes: Vec<(String, u64)> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect();
    let sources = routes::Sources::open(&graph.source_paths, lanes)?;
    println!(
        "[diag] lane 9564 = {}, lane 5231 = {}",
        graph.lanes[9564].name, graph.lanes[5231].name
    );

    // The truth route's S288C chrIII piece (the phantom's origin).
    let slice = sources.fetch(9564, 80018, 90183)?;
    println!("[diag] fetched S288C chrIII [80018,90183): {} bp", slice.len());
    let spans = oracle_segment_spans(&panel, &slice)?;
    let flood_keys: Vec<Vec<u64>> = vec![
        vec![14634400],
        vec![31524560],
        vec![4144956],
        vec![24621440],
    ];
    for key in &flood_keys {
        match spans.get(key) {
            Some(list) => println!(
                "[diag] piece spans for {:?}: {} spans, first {:?}",
                key,
                list.len(),
                list.first()
            ),
            None => println!("[diag] piece spans for {:?}: ABSENT", key),
        }
    }

    // The emitting tagged records: windows starting in [2315, 2402] place a
    // flood anchor's k-mer at piece-relative [2402,2465) iff start + rel = lo.
    let mut flood_anchor_hits = 0usize;
    for start in 2315usize..=2402usize {
        let window = &slice[start..start + READ_LENGTH];
        for record in mem_records::tagged_mem_records(&panel, window)? {
            for (anchor_index, &(node, pos)) in record.iter().enumerate() {
                let ok_token = 2u64 + 2 * (((node as i64) << 1) ^ ((node as i64) >> 63)) as u64;
                let rc_token = 2u64
                    + 2 * ((((-(node as i64))) << 1) ^ ((-(node as i64)) >> 63)) as u64;
                if flood_keys.contains(&vec![ok_token]) || flood_keys.contains(&vec![rc_token]) {
                    flood_anchor_hits += 1;
                    if flood_anchor_hits <= 4 {
                        println!(
                            "[diag] FLOOD ANCHOR in run start={}: record anchors {:?} (hit index {} pos {})",
                            start, record, anchor_index, pos
                        );
                        let pos = pos as usize;
                        println!(
                            "[diag]   read window ACGT at k-mer [{}..={}): {}",
                            start + pos,
                            start + pos + k as usize,
                            String::from_utf8_lossy(&slice[start + pos..start + pos + k as usize]).to_string().as_str()
                        );
                    }
                }
            }
        }
    }
    println!("[diag] flood-anchor tagged-record hits in the candidate windows: {}", flood_anchor_hits);

    // The panel's own view of both loci.
    let s_window = sources.fetch(9564, 82320, 82600)?;
    println!("[diag] S288C ACGT [82420,82483): {}", String::from_utf8_lossy(&s_window[100..163]).to_string().as_str());
    for m in panel.matched_syncmers_in_sequence(&s_window) {
        if (82400..82510).contains(&(82320 + m.query_pos)) {
            println!(
                "[diag] S288C matched syncmer node {} at bp {}",
                m.signed_node,
                82320 + m.query_pos
            );
        }
    }
    let b_window = sources.fetch(5231, 165996, 166185)?;
    println!("[diag] BTE ACGT [166059,166122): {}", String::from_utf8_lossy(&b_window[63..126]).to_string().as_str());
    for m in panel.matched_syncmers_in_sequence(&b_window) {
        if (166030..166140).contains(&(165996 + m.query_pos)) {
            println!(
                "[diag] BTE matched syncmer node {} at bp {}",
                m.signed_node,
                165996 + m.query_pos
            );
        }
    }
    Ok(())
}
