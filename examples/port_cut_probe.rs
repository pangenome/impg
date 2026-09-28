//! The ITEM-2 cut-gap probe (owner ruling 2026-09-28, measure-first): for
// the mosaic truth's interior junctions (a run's stage-1
/// `truth_interior_junction_census` array), query the panel's port index
// at the truth's EXACT junction cuts — the decisive measurement for the
/// cut-semantics fork: a port at the exact cut means the port-cut
// materialization can express the truth's pieces exactly (the gate can
// pass); the nearest-port distance is reported either way (the bracket's
/// port content). Diagnostic only; no production path reads it.
//
// Usage: port_cut_probe <routes-dir> <truth-junctions.json>
// (each row's left/right fields are [source, start, end, reverse]; the
// probe reports every port within +-2000bp of the cut and the at-cut
// exactness.)

use std::io;
use std::path::PathBuf;

fn main() -> io::Result<()> {
    let mut args = std::env::args().skip(1);
    let routes_dir = PathBuf::from(args.next().expect("routes dir"));
    let junctions_path = args.next().expect("truth junctions json");
    let graph: impg::genome_inference::panel_routes::Graph =
        impg::genome_inference::read_json(&routes_dir.join("graph.json"))?;
    let mut ports =
        impg::genome_inference::panel_routes::Ports::open_without_global_verification(
            &routes_dir,
            &graph,
        )?;
    let text = std::fs::read_to_string(&junctions_path)?;
    let value: serde_json::Value = {
        let start = text
            .find('{')
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "no JSON"))?;
        serde_json::from_str(&text[start..])?
    };
    let rows = value
        .get("junctions")
        .and_then(|list| list.as_array())
        .or_else(|| value.as_array())
        .cloned()
        .unwrap_or_default();
    println!("probe: {} truth junction rows", rows.len());
    // THE SEAM BRACKET (the supervisor's class-A deliverable): for each
    // truth junction, measure the CONSERVED-GAP extent — the identity
    // run between the two frames spanning the seam — by direct sequence
    // comparison around the aligned cut positions. The seam is
    // unlocalizable within the bracket (every read spanning only
    // bracket material is native-explainable on either frame); the
    // falsifiable prediction: reads longer than the bracket plus one
    // divergent anchor (k-mer) on each side WOULD attest the seam.
    let lanes: Vec<(String, u64)> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect();
    let sources =
        impg::genome_inference::panel_routes::Sources::open(&graph.source_paths, lanes)?;
    let k = graph.k;
    let window = 3000u64;
    for (index, row) in rows.iter().enumerate() {
        let field = |name: &str| -> Option<(usize, u64, u64, bool)> {
            let node = row.get(name)?;
            if let (Some(s), Some(st), Some(e), Some(r)) = (
                node.get("source")?.as_u64(),
                node.get("start")?.as_u64(),
                node.get("end")?.as_u64(),
                node.get("reverse")?.as_bool(),
            ) {
                return Some((s as usize, st, e, r));
            }
            Some((
                node.get(0)?.as_u64()? as usize,
                node.get(1)?.as_u64()?,
                node.get(2)?.as_u64()?,
                node.get(3)?.as_bool()?,
            ))
        };
        let (Some(left), Some(right)) = (field("left"), field("right")) else {
            continue;
        };
        // The truth's junction cuts: the left piece ENDS at the cut (its
        // end on the left source); the right piece STARTS at the cut (its
        // start on the right source).
        let left_cut = left.2;
        let right_cut = right.1;
        // The frames' local offset at the seam: left position x on the
        // left source corresponds to x + (right_cut - left_cut) on the
        // right source.
        let offset = right_cut as i64 - left_cut as i64;
        let fetch_left = sources.fetch(left.0, left_cut.saturating_sub(window), (left_cut + window).min(u64::MAX))?;
        let fetch_right = sources.fetch(right.0, right_cut.saturating_sub(window), (right_cut + window).min(u64::MAX))?;
        // Identity runs outward from the seam, on both sides.
        let at = |seq: &[u8], pos: u64, base_lo: u64| -> u8 {
            seq[(pos - base_lo) as usize]
        };
        let left_lo = left_cut.saturating_sub(window);
        let right_lo = right_cut.saturating_sub(window);
        let left_hi = left_lo + fetch_left.len() as u64 - 1;
        let right_hi = right_lo + fetch_right.len() as u64 - 1;
        // Walk left from the seam while both frames carry the same base.
        let mut left_run = 0u64;
        loop {
            let i = left_run;
            if left_cut < left_lo + i + 1 || right_cut < right_lo + i + 1 {
                break;
            }
            if at(&fetch_left, left_cut - i - 1, left_lo)
                != at(&fetch_right, right_cut - i - 1, right_lo)
            {
                break;
            }
            left_run += 1;
        }
        // Walk right from the seam while both frames carry the same base.
        let mut right_run = 0u64;
        loop {
            let i = right_run;
            if left_cut + i > left_hi || right_cut + i > right_hi {
                break;
            }
            if at(&fetch_left, left_cut + i, left_lo)
                != at(&fetch_right, right_cut + i, right_lo)
            {
                break;
            }
            right_run += 1;
        }
        let bracket_lo = left_cut - left_run;
        let bracket_hi = left_cut + right_run;
        let gap = left_run + right_run;
        let need = gap + 2 * k;
        println!(
            "seam-bracket {index}: left source {} cut {left_cut}, right source {} cut {right_cut}, \
             offset {offset}: conserved gap = {gap} bp (left_run {left_run} + right_run {right_run}); \
             bracket [{bracket_lo}, {bracket_hi}] on the left frame; \
             required bridging read >= {need} bp",
            left.0,
            right.0,
        );
        let mut near = |source: usize, cut: u64, label: &str, reverse: bool| {
            let lo = cut.saturating_sub(2_000);
            let hi = cut + 2_000;
            let list = ports
                .forward_ports_inside(&graph, source, lo, hi)
                .unwrap_or_default();
            let mut nearest: Option<(u64, i64)> = None;
            let mut exact = false;
            for port in &list {
                if port.reverse != reverse {
                    continue;
                }
                let port_cut = port.cut(graph.k);
                let distance = port_cut as i64 - cut as i64;
                if distance == 0 {
                    exact = true;
                }
                if nearest.is_none_or(|(_, best)| distance.abs() < best) {
                    nearest = Some((port_cut, distance));
                }
            }
            println!(
                "junction {index} {label}: source {source} cut {cut} reverse {reverse} \
                 ports_in_window {} exact_at_cut {exact} nearest {nearest:?}",
                list.len()
            );
        };
        near(left.0, left_cut, "left", left.3);
        near(right.0, right_cut, "right", right.3);
    }
    Ok(())
}
