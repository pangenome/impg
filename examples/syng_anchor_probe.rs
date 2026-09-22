//! Print raw anchor positions for a named query on a specific target path,
//! so we can compare against what the refined interval output claims.
//!
//! Usage:
//!   cargo run --release --example syng_anchor_probe -- \
//!       <syng_prefix> <query_name> <qs> <qe> <target_path>
//!
//! Example:
//!   cargo run --release --example syng_anchor_probe -- \
//!       /home/erik/scrapy/yeast235 'S288C#0#chrIV' 408000 410000 'AAA#0#chrIV'

use impg::syng::{SyncmerParams, SyngIndex};
use std::env;

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 6 {
        eprintln!(
            "usage: {} <syng_prefix> <query_name> <qs> <qe> <target_path>",
            args[0]
        );
        std::process::exit(1);
    }
    let prefix = &args[1];
    let query_name = &args[2];
    let qs: u64 = args[3].parse().expect("qs must be u64");
    let qe: u64 = args[4].parse().expect("qe must be u64");
    let target_path = &args[5];

    let idx = SyngIndex::load(prefix, SyncmerParams::default()).expect("failed to load syng index");

    eprintln!("Query: {}:{}-{}", query_name, qs, qe);
    eprintln!("Target: {}", target_path);
    eprintln!();

    let hits = idx
        .query_region_with_anchors(query_name, qs, qe, 120)
        .expect("query_region_with_anchors failed");

    let mut matched = 0usize;
    for h in &hits {
        if h.genome != *target_path {
            continue;
        }
        matched += 1;
        eprintln!(
            "Interval {}: [{}, {}) strand={} padded_width={} n_anchors={}",
            matched,
            h.start,
            h.end,
            h.strand,
            h.end - h.start,
            h.anchors.len()
        );
        let mut anchors = h.anchors.clone();
        anchors.sort_by_key(|a| a.query_pos);
        println!(
            "# {}:{}-{}  strand={}",
            target_path, h.start, h.end, h.strand
        );
        println!("# q_pos\tt_pos\tdelta_t_minus_q\tnode_id");
        for a in &anchors {
            let dt: i64 = a.target_pos as i64 - a.query_pos as i64;
            println!("{}\t{}\t{}\t{}", a.query_pos, a.target_pos, dt, a.node_id);
        }
        // First and last-anchor-derived projection for comparison
        if let (Some(first), Some(last)) = (anchors.first(), anchors.last()) {
            let q_span = last.query_pos.saturating_sub(first.query_pos);
            let t_span = (last.target_pos as i64 - first.target_pos as i64).unsigned_abs();
            eprintln!(
                "  first_anchor (q={}, t={}, dt={})",
                first.query_pos,
                first.target_pos,
                first.target_pos as i64 - first.query_pos as i64
            );
            eprintln!(
                "  last_anchor  (q={}, t={}, dt={})",
                last.query_pos,
                last.target_pos,
                last.target_pos as i64 - last.query_pos as i64
            );
            eprintln!(
                "  anchor span: q_span={}bp t_span={}bp indel_imbalance={}bp",
                q_span,
                t_span,
                (t_span as i64 - q_span as i64).abs()
            );
        }
        eprintln!();
    }
    if matched == 0 {
        eprintln!("No anchored hits found for target '{}'", target_path);
    }
}
