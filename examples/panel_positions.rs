//! PANEL POSITIONS — `node -> (path, offset)`, the key the evidence store needs.
//!
//! The syng NODE ID is a graph identifier, not a coordinate: a median
//! 3-node walk spans ~47 million node ids (measured), so bucketing
//! evidence by node id scatters every walk (63% of reads span >100
//! 10k-node buckets). Bucketing by PANEL POSITION is what makes a walk
//! land in one bucket — and `(path, offset)` is exactly the coordinate the
//! pipeline already speaks (the beds are `path_name, start, end`).
//!
//! This walks the panel's own paths through `SyngIndex::walk_path_range`
//! (the same accessor `syng_homology_dump walks` uses) and writes a sorted
//! `node -> (path, offset)` table: a PANEL-SIDE precompute, mmappable,
//! built once per panel and reused by every sample.
//!
//! ```text
//! panel_positions <syng-prefix> <out-file> [first-path] [last-path] [chunk-bp]
//! ```
//!
//! Output: a fixed-width table, sorted by node —
//! `i32 node, u32 path, u64 offset` (16 bytes), then a header with the
//! counts and the content key (path count and total bp) so a mismatch
//! against the panel fails loudly rather than silently mis-keying.

use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::time::Instant;

use impg::syng::{SyncmerParams, SyngIndex};

const MAGIC: &[u8; 8] = b"PNLPOS01";
const CHUNK_BP: u64 = 10_000_000;

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: panel_positions <syng-prefix> <out-file> [first-path] [last-path] [chunk-bp]");
        return Ok(());
    }
    let prefix = &args[1];
    let out_path = &args[2];
    let first: usize = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0);
    let last: Option<usize> = args.get(4).and_then(|s| s.parse().ok());
    let chunk: u64 = args.get(5).and_then(|s| s.parse().ok()).unwrap_or(CHUNK_BP);

    let t0 = Instant::now();
    let index = SyngIndex::load(prefix, SyncmerParams::default())?;
    println!("panel loaded in {:.1}s", t0.elapsed().as_secs_f64());

    let path_count = index.name_map.path_to_name.len();
    let end_path = last.unwrap_or(path_count).min(path_count);
    println!("paths {first}..{end_path} of {path_count}");

    // ---- walk the paths, collecting (node, path, offset) ----
    let t1 = Instant::now();
    let mut entries: Vec<(i32, u32, u64)> = Vec::new();
    let mut total_bp = 0u64;
    for path_idx in first..end_path {
        // the path's length comes from the name map's own table
        let length = index.name_map.path_to_length[path_idx];
        total_bp += length;
        let mut at = 0u64;
        while at < length {
            let stop = (at + chunk).min(length);
            let steps = index.walk_path_range(path_idx, at, stop)?;
            for (node, offset) in steps {
                entries.push((node, path_idx as u32, offset));
            }
            at = stop;
        }
    }
    let walk = t1.elapsed();

    let t2 = Instant::now();
    entries.sort_unstable();
    entries.dedup();
    let sort = t2.elapsed();

    let t3 = Instant::now();
    {
        let mut out = BufWriter::new(File::create(out_path)?);
        out.write_all(MAGIC)?;
        out.write_all(&(entries.len() as u64).to_le_bytes())?;
        out.write_all(&(path_count as u64).to_le_bytes())?;
        out.write_all(&(end_path as u64).to_le_bytes())?;
        out.write_all(&total_bp.to_le_bytes())?;
        for &(node, path, offset) in &entries {
            out.write_all(&node.to_le_bytes())?;
            out.write_all(&path.to_le_bytes())?;
            out.write_all(&offset.to_le_bytes())?;
        }
    }
    let write = t3.elapsed();
    let bytes = std::fs::metadata(out_path)?.len();

    println!(
        "  nodes {}  ({:.1} B/node incl. header)  bp {total_bp}",
        entries.len(),
        bytes as f64 / entries.len().max(1) as f64
    );
    println!(
        "  walk {:.1}s  sort {:.1}s  write {:.1}s  -> {} ({:.1} MB)",
        walk.as_secs_f64(),
        sort.as_secs_f64(),
        write.as_secs_f64(),
        out_path,
        bytes as f64 / 1e6
    );
    Ok(())
}
