//! NODE -> READ POSTINGS — build `<prefix>.n2r` from the evidence store.
//!
//! ```text
//! node_postings build    <prefix>    # READ <prefix>.evd{,.idx}, WRITE <prefix>.n2r
//! node_postings summary  <n2r>       # header of a built artifact
//! node_postings postings <n2r> <node>
//! ```
//!
//! ADDITIVE: the store (`.evd`, `.evd.idx`, 22.4 GB, in production use) is
//! opened read-only and its format is untouched; the artifact is a new sibling
//! file. `build` is one full pass over the store — 22.4 GB read, 168 GB
//! decompressed — because a read's node attribution lives in its packet's walk,
//! so no header carries it.
//!
//! `IMPG_NODE_POSTINGS=<prefix>` makes the restricted cut use the artifact as
//! its candidate-read supply; that is the only consumer today.

use std::io;
use std::path::Path;

use impg::genome_inference::node_postings;

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let arg = |i: usize| args.get(i).cloned().unwrap_or_default();
    match args.get(1).map(String::as_str) {
        Some("build") => {
            let prefix = arg(2);
            let summary = node_postings::build(&prefix)?;
            println!(
                "{prefix}.n2r (node -> read postings, version {})",
                node_postings::VERSION
            );
            println!(
                "  node_count      {} nodes with postings",
                summary.node_count
            );
            println!(
                "  postings_count  {} incidences ({} raw, before the multi-placement dedup)",
                summary.postings_count, summary.incidences
            );
            println!(
                "  artifact        {} ({:.2} GB)",
                summary.artifact_bytes,
                summary.artifact_bytes as f64 / 1e9
            );
            println!(
                "  store read      {} buckets, {} records, {} with walks; {:.1} GB read, {:.1} GB decompressed",
                summary.buckets,
                summary.records,
                summary.reads_with_walks,
                summary.store_bytes_read as f64 / 1e9,
                summary.store_bytes_decoded as f64 / 1e9
            );
            Ok(())
        }
        Some("summary") => {
            let postings = node_postings::NodePostings::open(Path::new(&arg(2)))?;
            println!(
                "{}: {} nodes, {} postings",
                arg(2),
                postings.node_count(),
                postings.postings_count()
            );
            Ok(())
        }
        Some("postings") => {
            let postings = node_postings::NodePostings::open(Path::new(&arg(2)))?;
            let node: i32 = arg(3).parse().map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidInput, "node must be an integer")
            })?;
            let reads = postings.postings(node);
            println!("node {node}: {} reads", reads.len());
            for (index, read) in reads.iter().enumerate() {
                if index % 16 == 0 {
                    print!("\n  ");
                }
                print!("{read} ");
            }
            println!();
            Ok(())
        }
        _ => {
            eprintln!(
                "usage: node_postings build <prefix> | summary <n2r> | postings <n2r> <node>"
            );
            Ok(())
        }
    }
}
