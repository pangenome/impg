//! Dump a loaded `SyngIndex` to GFA (1.0 with P lines by default,
//! 1.1 with W lines if `--gfa-version 1.1`).
//!
//! Output: one S per syncmer (node id 1..N, canonical DNA from the
//! KmerHash) plus one S per inter-syncmer gap (id N+1, N+2, …). When
//! the bp offset between two consecutive syncmers in a path exceeds the
//! syncmer length, a gap segment of `offset - syncmer_length` bp is
//! inserted between them. Gap bases come from `--sequence-files` if
//! provided, otherwise the gap is filled with `N`s and a warning is
//! emitted. Each gap occurrence gets a fresh segment id (no
//! cross-path dedup) so that orientation/dedup never silently mixes
//! distinct DNA.

use std::io::{self, BufWriter, Write};

use log::{debug, info, warn};
use rustc_hash::FxHashSet;

use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
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

/// One step in a P/W walk: either a syncmer (signed by strand) or an
/// inserted gap segment (always on the forward strand).
#[derive(Clone)]
enum PathStep {
    Syncmer(i32),
    Gap(String),
}

/// One inserted gap segment ready to be written as an S line.
struct GapSegment {
    id: String,
    seq: Vec<u8>,
}

/// Walked path plus any inserted gap segments.
struct WalkedPath {
    name: String,
    steps: Vec<PathStep>,
}

/// Walk every forward path in `index`, gap-fill via FASTA/AGC if given,
/// and write a GFA to `writer`.
///
/// `gap_fill` selects how mid-path gaps are filled (FASTA vs. NNN).
/// Returns `(n_segments, n_links, n_paths_emitted, n_paths_skipped,
/// n_gap_segments, total_gap_bp)`.
pub fn write_gfa<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    let n_nodes = index.num_syncmer_nodes();
    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] {} syncmer segments (length {} bp), GFA {}, gap fill: {}",
        n_nodes,
        syncmer_len,
        version.header_tag(),
        if gap_fill.is_some() { "FASTA/AGC" } else { "Ns" }
    );

    // Silence syng's C-side debug printfs (syngBWTpathStartOld et al.)
    unsafe { syng_ffi::impg_syng_suppress_debug() };

    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    // Walk every forward path and collect:
    //   - per-path step list (syncmers + gap segments) for P/W emission
    //   - unique L edges with overlaps
    //   - in-order list of gap segments to emit later as S lines
    let mut walked_paths: Vec<WalkedPath> =
        Vec::with_capacity(index.name_map.path_to_name.len());
    // Edge tuple: (from_id, from_sign, to_id, to_sign, overlap_bp).
    // `String` because endpoints are integer IDs (syncmers 1..N, gaps
    // N+1, N+2, …) and we serialize them as text.
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut gap_segments: Vec<GapSegment> = Vec::new();
    let mut n_skipped = 0usize;
    let mut next_gap_id: u64 = 0;
    let mut gaps_filled_with_ns: usize = 0;
    let mut total_gap_bp: u64 = 0;

    let syncmer_len_u64 = syncmer_len as u64;

    for (path_idx, name) in index.name_map.path_to_name.iter().enumerate() {
        let Some(start_info): Option<&GbwtPathStart> = index
            .name_map
            .path_starts
            .get(path_idx)
            .and_then(|o| o.as_ref())
        else {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping (rebuild the syng index to include path starts)",
                name
            );
            n_skipped += 1;
            continue;
        };
        let nodes = index.walk_forward_path(start_info);
        let path_len_bp = index
            .name_map
            .path_to_length
            .get(path_idx)
            .copied()
            .unwrap_or(0);
        let mut steps: Vec<PathStep> = Vec::with_capacity(nodes.len());

        // Prefix: bases [0, first_syncmer_pos) before the first syncmer.
        if let Some(&(first, first_pos)) = nodes.first() {
            if first_pos > 0 {
                let gap_seq = if let Some(seq_idx) = gap_fill {
                    fetch_gap_dna(seq_idx, name, 0, first_pos)?
                } else {
                    gaps_filled_with_ns += 1;
                    vec![b'N'; first_pos as usize]
                };
                total_gap_bp += gap_seq.len() as u64;
                let gap_id = format!("{}", n_nodes as u64 + 1 + next_gap_id);
                next_gap_id += 1;
                let (first_id, first_sign) = signed_to_gfa(first);
                edges.insert((
                    gap_id.clone(),
                    '+',
                    first_id.to_string(),
                    first_sign,
                    0,
                ));
                steps.push(PathStep::Gap(gap_id.clone()));
                gap_segments.push(GapSegment {
                    id: gap_id,
                    seq: gap_seq,
                });
            }
            steps.push(PathStep::Syncmer(first));
        }

        for w in nodes.windows(2) {
            let (a, pos_a) = w[0];
            let (b, pos_b) = w[1];
            let off_u64: u64 = pos_b.saturating_sub(pos_a);
            let (a_id, a_sign) = signed_to_gfa(a);
            let (b_id, b_sign) = signed_to_gfa(b);
            if off_u64 <= syncmer_len_u64 {
                // Direct adjacency: overlap = syncmer_len - offset.
                let overlap = (syncmer_len_u64 - off_u64) as u32;
                edges.insert((
                    a_id.to_string(),
                    a_sign,
                    b_id.to_string(),
                    b_sign,
                    overlap,
                ));
                steps.push(PathStep::Syncmer(b));
            } else {
                // Inter-syncmer gap of (offset - syncmer_len) bp.
                let gap_len = off_u64 - syncmer_len_u64;
                let gap_seq = if let Some(seq_idx) = gap_fill {
                    let gap_start = pos_a + syncmer_len_u64;
                    let gap_end = pos_b;
                    fetch_gap_dna(seq_idx, name, gap_start, gap_end)?
                } else {
                    gaps_filled_with_ns += 1;
                    vec![b'N'; gap_len as usize]
                };
                total_gap_bp += gap_seq.len() as u64;
                let gap_id = format!("{}", n_nodes as u64 + 1 + next_gap_id);
                next_gap_id += 1;
                edges.insert((
                    a_id.to_string(),
                    a_sign,
                    gap_id.clone(),
                    '+',
                    0,
                ));
                edges.insert((
                    gap_id.clone(),
                    '+',
                    b_id.to_string(),
                    b_sign,
                    0,
                ));
                steps.push(PathStep::Gap(gap_id.clone()));
                steps.push(PathStep::Syncmer(b));
                gap_segments.push(GapSegment {
                    id: gap_id,
                    seq: gap_seq,
                });
            }
        }
        // Suffix: bases [last_pos + syncmer_len, path_len) after the last syncmer.
        if let Some(&(last, last_pos)) = nodes.last() {
            let suffix_start = last_pos.saturating_add(syncmer_len_u64);
            if path_len_bp > suffix_start {
                let suffix_end = path_len_bp;
                let gap_seq = if let Some(seq_idx) = gap_fill {
                    fetch_gap_dna(seq_idx, name, suffix_start, suffix_end)?
                } else {
                    gaps_filled_with_ns += 1;
                    vec![b'N'; (suffix_end - suffix_start) as usize]
                };
                total_gap_bp += gap_seq.len() as u64;
                let gap_id = format!("{}", n_nodes as u64 + 1 + next_gap_id);
                next_gap_id += 1;
                let (last_id, last_sign) = signed_to_gfa(last);
                edges.insert((
                    last_id.to_string(),
                    last_sign,
                    gap_id.clone(),
                    '+',
                    0,
                ));
                steps.push(PathStep::Gap(gap_id.clone()));
                gap_segments.push(GapSegment {
                    id: gap_id,
                    seq: gap_seq,
                });
            }
        }

        walked_paths.push(WalkedPath {
            name: name.clone(),
            steps,
        });
    }

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} inter-syncmer gap(s) with 'N' ({} bp total). \
             Concatenating these paths will NOT reconstruct the original genome. \
             Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns, total_gap_bp
        );
    } else if !gap_segments.is_empty() {
        info!(
            "[syng2gfa] spliced {} real-DNA gap segment(s), {} bp total",
            gap_segments.len(),
            total_gap_bp
        );
    }

    // S lines: syncmer nodes 1..=N first, then gap segments in walk order.
    // Syncmer DNA from the C lib is lowercase; uppercase for consistency with
    // FASTA-fetched gap sequences (already uppercased upstream).
    for node_id in 1..=n_nodes as i32 {
        let mut seq = index.syncmer_seq(node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for gap in &gap_segments {
        write_segment(writer, &gap.id, &gap.seq)?;
    }
    let n_segments = n_nodes + gap_segments.len();
    let n_gap_segments = gap_segments.len();

    debug!("[syng2gfa] collected {} unique edges", edges.len());

    // L lines, sorted for deterministic output.
    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(
            writer,
            "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M"
        )?;
    }
    let n_links = edge_vec.len();

    // P/W lines.
    let n_paths_emitted = walked_paths.len();
    for walked in &walked_paths {
        match version {
            GfaVersion::V1_0 => write_p_line(writer, &walked.name, &walked.steps)?,
            GfaVersion::V1_1 => {
                if let Some(parts) = parse_pansn(&walked.name) {
                    write_w_line(writer, &parts, &walked.steps)?;
                } else {
                    warn!(
                        "[syng2gfa] path '{}' is not PanSN-style; falling back to P line under GFA 1.1",
                        walked.name
                    );
                    write_p_line(writer, &walked.name, &walked.steps)?;
                }
            }
        }
    }

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        n_gap_segments,
        total_gap_bp,
    ))
}

/// Pull the DNA between two syncmers from a sequence index.
fn fetch_gap_dna(
    seq_idx: &UnifiedSequenceIndex,
    seq_name: &str,
    gap_start: u64,
    gap_end: u64,
) -> io::Result<Vec<u8>> {
    if gap_end <= gap_start {
        return Ok(Vec::new());
    }
    if gap_start > i32::MAX as u64 || gap_end > i32::MAX as u64 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "gap coordinates exceed i32 range on '{}': [{}, {})",
                seq_name, gap_start, gap_end
            ),
        ));
    }
    seq_idx.fetch_sequence(seq_name, gap_start as i32, gap_end as i32)
}

fn signed_to_gfa(signed: i32) -> (i32, char) {
    if signed >= 0 {
        (signed, '+')
    } else {
        (-signed, '-')
    }
}

fn write_segment<W: Write>(writer: &mut W, id: &str, seq: &[u8]) -> io::Result<()> {
    writer.write_all(b"S\t")?;
    writer.write_all(id.as_bytes())?;
    writer.write_all(b"\t")?;
    writer.write_all(seq)?;
    writer.write_all(b"\n")?;
    Ok(())
}

fn write_p_line<W: Write>(writer: &mut W, name: &str, steps: &[PathStep]) -> io::Result<()> {
    writer.write_all(b"P\t")?;
    writer.write_all(name.as_bytes())?;
    writer.write_all(b"\t")?;
    for (i, step) in steps.iter().enumerate() {
        if i > 0 {
            writer.write_all(b",")?;
        }
        match step {
            PathStep::Syncmer(signed) => {
                let (id, sign) = signed_to_gfa(*signed);
                write!(writer, "{id}{sign}")?;
            }
            PathStep::Gap(id) => {
                write!(writer, "{id}+")?;
            }
        }
    }
    writer.write_all(b"\t*\n")?;
    Ok(())
}

fn write_w_line<W: Write>(
    writer: &mut W,
    parts: &PansnParts<'_>,
    steps: &[PathStep],
) -> io::Result<()> {
    write!(
        writer,
        "W\t{}\t{}\t{}\t{}\t{}\t",
        parts.sample, parts.hap, parts.contig, parts.start, parts.end
    )?;
    for step in steps {
        match step {
            PathStep::Syncmer(signed) => {
                let (id, sign) = signed_to_gfa(*signed);
                let arrow = if sign == '+' { '>' } else { '<' };
                write!(writer, "{arrow}{id}")?;
            }
            PathStep::Gap(id) => {
                write!(writer, ">{id}")?;
            }
        }
    }
    writer.write_all(b"\n")?;
    Ok(())
}

/// Convenience: load a syng index from `prefix`, optionally build a
/// sequence index from `sequence_files`, and write the GFA to `out_path`
/// (`-` for stdout).
pub fn run(
    prefix: &str,
    out_path: &str,
    version: GfaVersion,
    sequence_files: &[String],
) -> io::Result<()> {
    info!("[syng2gfa] loading syng index from prefix '{}'", prefix);
    // SyncmerParams in load() is overridden by the .syng.meta sidecar, so any
    // placeholder works.
    let index = SyngIndex::load(prefix, crate::syng::SyncmerParams::default())?;

    let seq_idx: Option<UnifiedSequenceIndex> = if sequence_files.is_empty() {
        None
    } else {
        info!(
            "[syng2gfa] building sequence index for {} file(s) (gap fill)",
            sequence_files.len()
        );
        Some(UnifiedSequenceIndex::from_files(sequence_files)?)
    };

    let report =
        |s: usize, l: usize, p: usize, skipped: usize, gaps: usize, gap_bp: u64, where_: &str| {
            info!(
                "[syng2gfa] wrote {} S ({} syncmer + {} gap, {} gap bp), {} L, {} P/W lines ({} paths skipped) to {}",
                s,
                s - gaps,
                gaps,
                gap_bp,
                l,
                p,
                skipped,
                where_,
            );
        };

    if out_path == "-" {
        let stdout = io::stdout();
        let mut w = BufWriter::new(stdout.lock());
        let (s, l, p, skipped, gaps, gap_bp) =
            write_gfa(&index, &mut w, version, seq_idx.as_ref())?;
        w.flush()?;
        report(s, l, p, skipped, gaps, gap_bp, "stdout");
    } else {
        let file = std::fs::File::create(out_path)?;
        let mut w = BufWriter::new(file);
        let (s, l, p, skipped, gaps, gap_bp) =
            write_gfa(&index, &mut w, version, seq_idx.as_ref())?;
        w.flush()?;
        report(s, l, p, skipped, gaps, gap_bp, &format!("'{}'", out_path));
    }
    Ok(())
}
