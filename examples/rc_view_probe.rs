//! One-off probe for the rc-symmetric territory augmentation (the
//! strand-asymmetry fix's step-26 verification): does the panel's
//! matched-syncmer extraction on the reverse complement of the S288C chrIII
//! range surface the rc frame's qualifying k-mers (the flood feature's node
//! at the mirrored position), and in which coordinate frame do the returned
//! positions sit?
//!
//! Usage: rc_view_probe <panel.syng> <routes-dir>

use impg::genome_inference::{panel_routes as routes, read_json};
use impg::syng::{SyncmerParams, SyngIndex};
use std::io;

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let panel = SyngIndex::load(&args[1], SyncmerParams::default())?;
    let k = panel.syncmer_length_bp() as u64;
    println!("[probe] k = {k}");
    let routes_dir = std::path::Path::new(&args[2]);
    let graph: routes::Graph = read_json(&routes_dir.join("graph.json"))?;
    let lanes: Vec<(String, u64)> = graph
        .lanes
        .iter()
        .map(|lane| (lane.name.clone(), lane.length))
        .collect();
    let sources = routes::Sources::open(&graph.source_paths, lanes)?;
    // lane 9564 = S288C#0#chrIII (the spelled locus path).
    let (lo, hi) = (79868u64, 90333u64);
    let seq = sources.fetch(9564, lo, hi)?;
    println!("[probe] fetched lane 9564 [{lo},{hi}): {} bp", seq.len());
    let rc_seq = impg::graph::reverse_complement(&seq);
    let rc_len = rc_seq.len() as u64;

    // The forward view's matched syncmers (the panel's own steps).
    let fwd = panel.matched_syncmers_in_sequence(&seq);
    println!("[probe] matched syncmers of the range as-is: {}", fwd.len());
    // The rc view's matched syncmers (the augmentation's input).
    let rc_view = panel.matched_syncmers_in_sequence(&rc_seq);
    println!(
        "[probe] matched syncmers of the reversed range: {}",
        rc_view.len()
    );

    // The flood feature's node (hash of the BTE k-mer Q1) and the expected
    // rc-view hit: the k-mer rc(range)[q..q+k) = Q1 at q = rc_len - k -
    // (82420 - lo).
    let expect_q = rc_len - k - (82420 - lo);
    println!("[probe] expected rc-view q for the bp-82420 k-mer: {expect_q}");
    for m in &rc_view {
        if m.query_pos + 40 >= expect_q.saturating_sub(40) && m.query_pos <= expect_q + 40 {
            println!(
                "[probe] rc-view hit near expect: node {} at q {} (bp {})",
                m.signed_node,
                m.query_pos,
                lo + seq.len() as u64 - k - m.query_pos
            );
        }
    }
    // The first 12 rc-view hits raw (frame diagnosis).
    for m in rc_view.iter().take(12) {
        println!("[probe] rc-view raw: node {} q {}", m.signed_node, m.query_pos);
    }
    // The forward view's first 12 for comparison.
    for m in fwd.iter().take(12) {
        println!("[probe] fwd raw: node {} q {}", m.signed_node, m.query_pos);
    }
    Ok(())
}
