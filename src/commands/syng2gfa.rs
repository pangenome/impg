//! Dump a loaded `SyngIndex` to GFA (1.0 with P lines by default,
//! 1.1 with W lines if `--gfa-version 1.1`).
//!
//! Lossless wrt the syng index: one S per syncmer, one L per unique
//! adjacency observed across any forward path, one P or W per indexed
//! sequence. Inter-syncmer DNA gaps (when the bp offset between two
//! consecutive syncmers exceeds the syncmer length) are not present in
//! the index — those become `0M` links instead of carrying real
//! overlap.

use std::io::{self, BufWriter, Write};

use log::{debug, info};
use rustc_hash::FxHashSet;

use crate::syng::{GbwtPathStart, SyngIndex};
use crate::syng_ffi;

/// Which GFA spec version to emit.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum GfaVersion {
    /// Strict GFA v1.0 — H/S/L/P only.
    V1_0,
    /// GFA v1.1 — H/S/L/W (paths emitted as walks with PanSN metadata).
    V1_1,
}

impl GfaVersion {
    pub fn parse(s: &str) -> io::Result<Self> {
        match s {
            "1.0" | "v1.0" | "1" => Ok(GfaVersion::V1_0),
            "1.1" | "v1.1" => Ok(GfaVersion::V1_1),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unknown GFA version '{other}' (expected 1.0 or 1.1)"),
            )),
        }
    }

    fn header_tag(self) -> &'static str {
        match self {
            GfaVersion::V1_0 => "1.0",
            GfaVersion::V1_1 => "1.1",
        }
    }
}

/// Parsed PanSN-style path name.
struct PansnParts<'a> {
    sample: &'a str,
    hap: &'a str,
    contig: &'a str,
    start: u64,
    end: u64,
}

/// Try to parse `sample#hap#contig:start-end` or `sample#contig:start-end`.
///
/// Returns `None` if the name doesn't match. The hap index defaults to "0"
/// when only two `#`-fields are present (the chm13/grch38 convention used by
/// HPRC PanSN).
fn parse_pansn(name: &str) -> Option<PansnParts<'_>> {
    let hash_parts: Vec<&str> = name.split('#').collect();
    let (sample, hap, locus) = match hash_parts.as_slice() {
        [sample, locus] => (*sample, "0", *locus),
        [sample, hap, locus] => (*sample, *hap, *locus),
        _ => return None,
    };
    let (contig, range) = locus.rsplit_once(':')?;
    let (start_s, end_s) = range.split_once('-')?;
    let start: u64 = start_s.parse().ok()?;
    let end: u64 = end_s.parse().ok()?;
    Some(PansnParts {
        sample,
        hap,
        contig,
        start,
        end,
    })
}

/// Walk every forward path in `index` and write a GFA to `writer`.
///
/// Returns `(n_segments, n_links, n_paths_emitted, n_paths_skipped)`.
/// Paths with no path-start info in the names sidecar are skipped (warned).
pub fn write_gfa<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
) -> io::Result<(usize, usize, usize, usize)> {
    let n_nodes = index.num_syncmer_nodes();
    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] writing {} segments (syncmer length {} bp), GFA {}",
        n_nodes,
        syncmer_len,
        version.header_tag()
    );

    // Silence syng's C-side debug printfs (syngBWTpathStartOld et al.)
    unsafe { syng_ffi::impg_syng_suppress_debug() };

    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    for node_id in 1..=n_nodes as i32 {
        let seq = index.syncmer_seq(node_id);
        writer.write_all(b"S\t")?;
        write!(writer, "{node_id}\t")?;
        writer.write_all(&seq)?;
        writer.write_all(b"\n")?;
    }

    // First pass: walk every forward path, accumulate per-path node sequences
    // and the global set of unique (from, from_sign, to, to_sign, overlap) edges.
    let mut all_paths: Vec<(String, Vec<(i32, u64)>)> =
        Vec::with_capacity(index.name_map.path_to_name.len());
    let mut edges: FxHashSet<(i32, i32, u32)> = FxHashSet::default();
    let mut n_skipped = 0usize;

    for (path_idx, name) in index.name_map.path_to_name.iter().enumerate() {
        let Some(start_info): Option<&GbwtPathStart> = index
            .name_map
            .path_starts
            .get(path_idx)
            .and_then(|o| o.as_ref())
        else {
            log::warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping (rebuild the syng index to include path starts)",
                name
            );
            n_skipped += 1;
            continue;
        };
        let nodes = index.walk_forward_path(start_info);
        if nodes.len() < 2 {
            // single-node or empty path: no edges to register
            all_paths.push((name.clone(), nodes));
            continue;
        }
        for w in nodes.windows(2) {
            let (a, pos_a) = w[0];
            let (b, pos_b) = w[1];
            // bp offset between consecutive syncmer starts; cap overlap at syncmer length
            let off_u64: u64 = pos_b.saturating_sub(pos_a);
            let off: u32 = if off_u64 > u32::MAX as u64 {
                u32::MAX
            } else {
                off_u64 as u32
            };
            let overlap = (syncmer_len as u32).saturating_sub(off);
            edges.insert((a, b, overlap));
        }
        all_paths.push((name.clone(), nodes));
    }

    debug!("[syng2gfa] collected {} unique edges", edges.len());

    // Emit L lines (sorted for deterministic output).
    let mut edge_vec: Vec<(i32, i32, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a, b, overlap) in &edge_vec {
        let (a_id, a_sign) = signed_to_gfa(*a);
        let (b_id, b_sign) = signed_to_gfa(*b);
        writeln!(
            writer,
            "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M"
        )?;
    }
    let n_links = edge_vec.len();

    // Emit paths.
    let n_paths_emitted = all_paths.len();
    for (name, nodes) in &all_paths {
        match version {
            GfaVersion::V1_0 => write_p_line(writer, name, nodes)?,
            GfaVersion::V1_1 => {
                if let Some(parts) = parse_pansn(name) {
                    write_w_line(writer, &parts, nodes)?;
                } else {
                    log::warn!(
                        "[syng2gfa] path '{}' is not PanSN-style; falling back to P line under GFA 1.1",
                        name
                    );
                    write_p_line(writer, name, nodes)?;
                }
            }
        }
    }

    Ok((n_nodes, n_links, n_paths_emitted, n_skipped))
}

fn signed_to_gfa(signed: i32) -> (i32, char) {
    if signed >= 0 {
        (signed, '+')
    } else {
        (-signed, '-')
    }
}

fn write_p_line<W: Write>(
    writer: &mut W,
    name: &str,
    nodes: &[(i32, u64)],
) -> io::Result<()> {
    writer.write_all(b"P\t")?;
    writer.write_all(name.as_bytes())?;
    writer.write_all(b"\t")?;
    for (i, (n, _)) in nodes.iter().enumerate() {
        if i > 0 {
            writer.write_all(b",")?;
        }
        let (id, sign) = signed_to_gfa(*n);
        write!(writer, "{id}{sign}")?;
    }
    writer.write_all(b"\t*\n")?;
    Ok(())
}

fn write_w_line<W: Write>(
    writer: &mut W,
    parts: &PansnParts<'_>,
    nodes: &[(i32, u64)],
) -> io::Result<()> {
    write!(
        writer,
        "W\t{}\t{}\t{}\t{}\t{}\t",
        parts.sample, parts.hap, parts.contig, parts.start, parts.end
    )?;
    for (n, _) in nodes {
        let (id, sign) = signed_to_gfa(*n);
        let arrow = if sign == '+' { '>' } else { '<' };
        write!(writer, "{arrow}{id}")?;
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Convenience: load a syng index from `prefix` and write the GFA to `out_path`
/// (`-` for stdout).
pub fn run(prefix: &str, out_path: &str, version: GfaVersion) -> io::Result<()> {
    info!("[syng2gfa] loading syng index from prefix '{}'", prefix);
    // SyncmerParams in load() is overridden by the .syng.meta sidecar, so any
    // placeholder works.
    let index = SyngIndex::load(prefix, crate::syng::SyncmerParams::default())?;

    if out_path == "-" {
        let stdout = io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        let (s, l, p, skipped) = write_gfa(&index, &mut w, version)?;
        w.flush()?;
        info!(
            "[syng2gfa] wrote {} S, {} L, {} P/W lines ({} paths skipped) to stdout",
            s, l, p, skipped
        );
    } else {
        let file = std::fs::File::create(out_path)?;
        let mut w = BufWriter::new(file);
        let (s, l, p, skipped) = write_gfa(&index, &mut w, version)?;
        w.flush()?;
        info!(
            "[syng2gfa] wrote {} S, {} L, {} P/W lines ({} paths skipped) to '{}'",
            s, l, p, skipped, out_path
        );
    }
    Ok(())
}
