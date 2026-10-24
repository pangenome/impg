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

/// Collapse an existing table onto `|node|` (one position per VERTEX, not
/// per orientation) - the fix for the reverse-strand placement hole.
fn normalise(in_path: &str, out_path: &str) -> io::Result<()> {
    use std::io::Read as _;
    let t0 = Instant::now();
    let mut bytes = Vec::new();
    File::open(in_path)?.read_to_end(&mut bytes)?;
    if &bytes[..8] != MAGIC {
        return Err(io::Error::new(io::ErrorKind::InvalidData, "bad table magic"));
    }
    let u64_at = |at: usize| u64::from_le_bytes(bytes[at..at + 8].try_into().unwrap());
    let n = u64_at(8) as usize;
    let (path_count, end_path, total_bp) = (u64_at(16), u64_at(24), u64_at(32));
    let mut entries: Vec<(i32, u32, u64)> = Vec::with_capacity(n);
    for i in 0..n {
        let at = 40 + i * 16;
        let node = i32::from_le_bytes(bytes[at..at + 4].try_into().unwrap());
        let path = u32::from_le_bytes(bytes[at + 4..at + 8].try_into().unwrap());
        let offset = u64::from_le_bytes(bytes[at + 8..at + 16].try_into().unwrap());
        entries.push((node.unsigned_abs() as i32, path, offset));
    }
    entries.sort_unstable();
    entries.dedup_by_key(|e| e.0); // one position per vertex
    let mut out = BufWriter::new(File::create(out_path)?);
    out.write_all(MAGIC)?;
    out.write_all(&(entries.len() as u64).to_le_bytes())?;
    out.write_all(&path_count.to_le_bytes())?;
    out.write_all(&end_path.to_le_bytes())?;
    out.write_all(&total_bp.to_le_bytes())?;
    for &(node, path, offset) in &entries {
        out.write_all(&node.to_le_bytes())?;
        out.write_all(&path.to_le_bytes())?;
        out.write_all(&offset.to_le_bytes())?;
    }
    out.flush()?;
    println!(
        "normalised {} -> {} entries ({:.1} MB) in {:.1}s",
        n,
        entries.len(),
        std::fs::metadata(out_path)?.len() as f64 / 1e6,
        t0.elapsed().as_secs_f64()
    );
    Ok(())
}

fn main() -> io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 3 {
        eprintln!("usage: panel_positions <syng-prefix> <out-file> [first-path] [last-path] [chunk-bp]");
        eprintln!("       panel_positions normalise <in-table> <out-table>   (collapse onto |node|)");
        return Ok(());
    }
    // the normalise pass needs NO panel: dispatch it before any syng load
    if args[1] == "normalise" {
        return normalise(&args[2], &args[3]);
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
    //
    // ONE POSITION PER NODE. The (path, node) INCIDENCES are ~58 billion
    // (1.4 Tbp of paths at ~24 bp per node) - a builder that keeps them
    // all wanted 930 GB and reached 93 GB before it was killed. We only
    // need ONE position per node, so a 34 MB bitset over the node id
    // space suppresses the repeats, and the walk STOPS once new paths
    // stop discovering nodes (the panel's vertices are shared across
    // thousands of paths, so a few dozen paths cover them).
    const NODE_SPACE: u64 = 300_000_000;
    let mut seen = vec![0u64; (NODE_SPACE / 64 + 1) as usize];
    let mut entries: Vec<(i32, u32, u64)> = Vec::new();
    let mut total_bp = 0u64;
    let mut barren = 0usize;
    // the panel's own vertex count (the GBWT load line), 99.5% of it as the
    // target; IMPG_POSITIONS_TARGET overrides for a scoped run
    let target_nodes: usize = std::env::var("IMPG_POSITIONS_TARGET")
        .ok()
        .and_then(|v| v.trim().parse().ok())
        .unwrap_or(271_000_000);
    let t1 = Instant::now();
    for path_idx in first..end_path {
        // the path's length comes from the name map's own table
        let length = index.name_map.path_to_length[path_idx];
        total_bp += length;
        let before = entries.len();
        let mut at = 0u64;
        while at < length {
            let stop = (at + chunk).min(length);
            let steps = index.walk_path_range(path_idx, at, stop)?;
            for (node, offset) in steps {
                let key = node.unsigned_abs() as u64;
                if key >= NODE_SPACE {
                    continue;
                }
                let word = (key / 64) as usize;
                let mask = 1u64 << (key % 64);
                if seen[word] & mask != 0 {
                    continue;
                }
                seen[word] |= mask;
                entries.push((node.unsigned_abs() as i32, path_idx as u32, offset));
            }
            at = stop;
        }
        // COVERAGE, not a plateau: the first 80 paths discovered 125.5M of
        // the panel's 272,213,874 vertices and the next 8 added <1% each -
        // a plateau rule stops there, but a half-covered table leaves half
        // the sample's reads unplaceable. Walk until the vertex count is
        // covered (IMPG_POSITIONS_TARGET overrides) or the paths run out.
        let discovered = entries.len() - before;
        if discovered * 100 < (length as usize / 24).max(1) {
            barren += 1;
        } else {
            barren = 0;
        }
        if path_idx % 100 == 0 {
            println!(
                "  path {path_idx}: {} nodes ({} barren in a row)",
                entries.len(),
                barren
            );
        }
        if entries.len() >= target_nodes {
            println!("  coverage target reached at path {path_idx}: {} nodes", entries.len());
            break;
        }
    }
    let walk = t1.elapsed();

    let t2 = Instant::now();
    entries.sort_unstable();
    entries.dedup_by_key(|e| e.0); // one position per vertex (orientation-free)
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
