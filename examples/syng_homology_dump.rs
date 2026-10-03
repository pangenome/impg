//! The multi-matching census's alignment-provenance dump (the 2026-11-05
//! post-projection census: the syng's own homologous-interval machinery
//! decides whether two occurrences' source material is ALIGNED TOGETHER).
//!
//! Two modes, both loading the panel once:
//!
//! `homology <prefix> <requests.tsv>` — batch homologous-interval dump.
//!   Each request line (tab-separated):
//!     `query_path_id <window> <qs> <qe> <padding> [target_path_ids|*] [anchors]`
//!   emits one JSON line per returned homologous interval (after the
//!   optional target filter):
//!     `{"query_path":Q,"window":w,"qs":..,"qe":..,"target_path":T,"start":..,
//!        "end":..,"strand":"+|-","anchors":[[qpos,tpos,node],..]|N}`
//!   then one request-completion marker per request (the receipt
//!   self-describes what was queried — the checker matches the markers
//!   against its own re-derived request set):
//!     `{"query_path":Q,"window":w,"qs":..,"qe":..,"request_done":true,"hits":N}`
//!   (full anchor lists only when the request asks for them; otherwise the
//!   anchor count — the window-alignment cache needs only interval bounds).
//!
//! `walk <prefix> <path_id> <start> <end>` — the panel path's node steps:
//!   one `bp<TAB>signed_node` line per covered syncmer step, for the
//!   near-twin condensation measurement (shared vs parallel nodes).
//!
//! `walks <prefix> <requests.tsv>` — batch walk mode (one panel load for
//!   many walks; added for the partition-graph stage 2026-11-06). Each
//!   request line (tab-separated):
//!     `path_id <start> <end> [tag]`
//!   emits ONE JSON line per walk:
//!     `{"path":P,"start":s,"end":e,"tag":t,"steps":[[bp,signed_node],...]}`
//!   used to verify that partition-graph GFA path lines spell exactly
//!   the panel's own walks in the global interning (id continuity).
//!
//! Assessment-side only: reads the panel, writes measurements; no
//! threshold enters anything (padding is per-request, set by the caller).

use impg::syng::{SyncmerParams, SyngIndex};
use std::env;
use std::io::{BufWriter, Write};

fn usage(program: &str) {
    eprintln!(
        "usage:\n  {} homology <syng_prefix> <requests.tsv>\n  {} walk <syng_prefix> <path_id> <start> <end>\n  {} walks <syng_prefix> <requests.tsv>",
        program, program, program
    );
    std::process::exit(1);
}

fn parse_walk_request(fields: &[&str], number: usize) -> (usize, u64, u64, String) {
    if fields.len() < 3 || fields.len() > 4 {
        panic!(
            "walk request {}: need 3 or 4 fields (path, start, end [tag]), got {}",
            number,
            fields.len()
        );
    }
    let path_id: usize = fields[0].parse().expect("path id must be usize");
    let start: u64 = fields[1].parse().expect("start must be u64");
    let end: u64 = fields[2].parse().expect("end must be u64");
    let tag = fields.get(3).copied().unwrap_or("").to_string();
    (path_id, start, end, tag)
}

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() < 3 {
        usage(&args[0]);
    }
    let mode = args[1].as_str();
    let prefix = &args[2];
    let idx = SyngIndex::load(prefix, SyncmerParams::default())
        .expect("failed to load syng index");

    match mode {
        "walk" => {
            if args.len() != 6 {
                usage(&args[0]);
            }
            let path_id: usize = args[3].parse().expect("path_id must be usize");
            let start: u64 = args[4].parse().expect("start must be u64");
            let end: u64 = args[5].parse().expect("end must be u64");
            let steps = idx
                .walk_path_range(path_id, start, end)
                .expect("walk_path_range failed");
            let stdout = std::io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            for (node, bp) in steps {
                writeln!(out, "{}\t{}", bp, node).expect("write step");
            }
            out.flush().expect("flush steps");
        }
        "walks" => {
            if args.len() != 4 {
                usage(&args[0]);
            }
            let requests = std::fs::read_to_string(&args[3]).expect("read walk requests");
            let stdout = std::io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            for (number, line) in requests.lines().enumerate() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split('\t').collect();
                let (path_id, start, end, tag) =
                    parse_walk_request(&fields, number + 1);
                let steps = idx
                    .walk_path_range(path_id, start, end)
                    .expect("walk_path_range failed");
                let line = serde_json::json!({
                    "path": path_id,
                    "start": start,
                    "end": end,
                    "tag": tag,
                    "steps": steps
                        .iter()
                        .map(|(node, bp)| serde_json::json!([bp, node]))
                        .collect::<Vec<_>>(),
                });
                writeln!(out, "{}", line).expect("write walk");
            }
            out.flush().expect("flush walks");
        }
        "homology" => {
            if args.len() != 4 {
                usage(&args[0]);
            }
            let requests = std::fs::read_to_string(&args[3]).expect("read requests");
            let stdout = std::io::stdout();
            let mut out = BufWriter::new(stdout.lock());
            for (number, line) in requests.lines().enumerate() {
                let line = line.trim();
                if line.is_empty() {
                    continue;
                }
                let fields: Vec<&str> = line.split('\t').collect();
                if fields.len() < 5 {
                    panic!("request {}: need >= 5 fields (path, window, qs, qe, padding)", number);
                }
                let query_path: usize =
                    fields[0].parse().expect("query path id must be usize");
                let window: u32 = fields[1].parse().expect("window must be u32");
                let qs: u64 = fields[2].parse().expect("qs must be u64");
                let qe: u64 = fields[3].parse().expect("qe must be u64");
                let padding: u64 = fields[4].parse().expect("padding must be u64");
                let targets: Option<Vec<usize>> = if fields.len() > 5 && fields[5] != "*" {
                    Some(
                        fields[5]
                            .split(',')
                            .map(|value| value.parse().expect("target id must be usize"))
                            .collect(),
                    )
                } else {
                    None
                };
                let with_anchors =
                    fields.len() > 6 && (fields[6] == "anchors" || fields[6] == "1");
                let query_name = idx.name_map.path_to_name[query_path].clone();
                let hits = idx
                    .query_region_with_anchors(&query_name, qs, qe, padding)
                    .expect("query_region_with_anchors failed");
                let mut hits_out = 0usize;
                for hit in &hits {
                    let target_path = match idx.name_map.name_to_path.get(&hit.genome) {
                        Some(&path) => path as usize,
                        None => continue,
                    };
                    if let Some(wanted) = &targets {
                        if !wanted.contains(&target_path) {
                            continue;
                        }
                    }
                    hits_out += 1;
                    let anchors = if with_anchors {
                        serde_json::json!(hit
                            .anchors
                            .iter()
                            .map(|anchor| {
                                [
                                    anchor.query_pos,
                                    anchor.target_pos,
                                    u64::from(anchor.node_id),
                                ]
                            })
                            .collect::<Vec<_>>())
                    } else {
                        serde_json::json!(hit.anchors.len())
                    };
                    let line = serde_json::json!({
                        "query_path": query_path,
                        "window": window,
                        "qs": qs,
                        "qe": qe,
                        "target_path": target_path,
                        "start": hit.start,
                        "end": hit.end,
                        "strand": hit.strand.to_string(),
                        "anchors": anchors,
                    });
                    writeln!(out, "{}", line).expect("write hit");
                }
                let marker = serde_json::json!({
                    "query_path": query_path,
                    "window": window,
                    "qs": qs,
                    "qe": qe,
                    "request_done": true,
                    "hits": hits_out,
                });
                writeln!(out, "{}", marker).expect("write marker");
                out.flush().expect("flush hits");
            }
        }
        _ => usage(&args[0]),
    }
}
