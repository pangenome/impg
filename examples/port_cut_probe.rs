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
    for (index, row) in rows.iter().enumerate() {
        let field = |name: &str| -> Option<(usize, u64, u64, bool)> {
            let node = row.get(name)?;
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
