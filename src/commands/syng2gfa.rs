//! Dump a loaded `SyngIndex` to GFA (1.0 with P lines by default,
//! 1.1 with W lines if `--gfa-version 1.1`).
//!
//! Output: one S per syncmer (node id 1..N, canonical DNA from the
//! KmerHash) plus one S per inter-syncmer gap (id N+1, N+2, …). When
//! the bp offset between two consecutive syncmers in a path exceeds the
//! syncmer length, a gap segment of `offset - syncmer_length` bp is
//! inserted between them. Gap bases come from `--sequence-files` if
//! provided, otherwise the gap is filled with `N`s and a warning is
//! emitted. Gap segments are interned by exact sequence plus local
//! signed-syncmer context, so identical terminal/inter-syncmer DNA shared
//! by the same local graph context is one node, without collapsing unrelated
//! repeated sequence elsewhere. The CLI default is blunt mode (`pangenome/bluntg`),
//! which converts native variable overlaps to 0M links; raw mode keeps
//! syng's native overlap graph.

use std::io::{self, BufWriter, Write};
use std::time::Instant;

use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::{FxHashMap, FxHashSet};

use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::syng::{GbwtPathStart, SyngIndex};
use crate::syng_ffi;

/// A source path interval to render as one GFA path.
pub struct SyngGfaPathRange {
    pub path_idx: usize,
    pub name: String,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

/// Syng GFA graph shape.
#[derive(Clone, Copy, Debug, Eq, PartialEq)]
pub enum SyngGfaMode {
    /// Native syng overlap graph, preserving nonzero link overlaps.
    Raw,
    /// Blunt graph produced by pangenome/bluntg; all link/path overlaps become 0M.
    Blunt,
}

impl SyngGfaMode {
    pub fn parse(s: &str) -> io::Result<Self> {
        match s.trim().replace('_', "-").to_ascii_lowercase().as_str() {
            "raw" | "syng:raw" | "syng-raw" | "syng-native:raw" => Ok(Self::Raw),
            "blunt" | "bluntg" | "syng" | "syng:blunt" | "syng:bluntg" | "syng-blunt"
            | "syng-bluntg" | "syng-native" | "syng-native:blunt" | "syng-native:bluntg" => {
                Ok(Self::Blunt)
            }
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("unknown syng GFA mode '{other}' (expected raw or blunt)"),
            )),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Raw => "raw",
            Self::Blunt => "blunt",
        }
    }
}

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

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
enum GapContext {
    Prefix {
        right: i32,
    },
    Between {
        left: i32,
        right: i32,
    },
    Suffix {
        left: i32,
    },
    PathOnly {
        path_idx: usize,
        start: u64,
        end: u64,
        strand: char,
    },
}

#[derive(Clone, Debug, Eq, Hash, PartialEq)]
struct GapKey {
    context: GapContext,
    seq: Vec<u8>,
}

struct GapInterner {
    n_nodes: usize,
    next_gap_id: u64,
    keys: FxHashMap<GapKey, String>,
    segments: Vec<GapSegment>,
    total_segment_bp: u64,
}

impl GapInterner {
    fn new(n_nodes: usize) -> Self {
        Self {
            n_nodes,
            next_gap_id: 0,
            keys: FxHashMap::default(),
            segments: Vec::new(),
            total_segment_bp: 0,
        }
    }

    fn intern(&mut self, context: GapContext, seq: Vec<u8>) -> String {
        let key = GapKey { context, seq };
        if let Some(id) = self.keys.get(&key) {
            return id.clone();
        }

        let id = format!("{}", self.n_nodes as u64 + 1 + self.next_gap_id);
        self.next_gap_id += 1;
        self.total_segment_bp += key.seq.len() as u64;
        self.segments.push(GapSegment {
            id: id.clone(),
            seq: key.seq.clone(),
        });
        self.keys.insert(key, id.clone());
        id
    }

    fn len(&self) -> usize {
        self.segments.len()
    }

    fn is_empty(&self) -> bool {
        self.segments.is_empty()
    }

    fn total_segment_bp(&self) -> u64 {
        self.total_segment_bp
    }
}

/// Walked path plus any inserted gap segments.
struct WalkedPath {
    name: String,
    steps: Vec<PathStep>,
}

#[derive(Clone)]
enum RawEndpoint {
    Syncmer(i32),
    Gap(usize),
}

#[derive(Clone)]
enum RawPathStep {
    Syncmer(i32),
    Gap(usize),
}

struct RawGap {
    context: GapContext,
    seq: Vec<u8>,
}

struct PathWork {
    name: String,
    gaps: Vec<RawGap>,
    steps: Vec<RawPathStep>,
    edges: Vec<(RawEndpoint, RawEndpoint, u32)>,
    gaps_filled_with_ns: usize,
    gap_occurrences: usize,
    gap_occurrence_bp: u64,
    skipped: bool,
}

fn raw_syncmer_endpoint(signed: i32) -> RawEndpoint {
    RawEndpoint::Syncmer(signed)
}

fn raw_gap_endpoint(idx: usize) -> RawEndpoint {
    RawEndpoint::Gap(idx)
}

fn raw_endpoint_to_gfa(endpoint: &RawEndpoint, gap_ids: &[String]) -> (String, char) {
    match endpoint {
        RawEndpoint::Syncmer(signed) => {
            let (id, sign) = signed_to_gfa(*signed);
            (id.to_string(), sign)
        }
        RawEndpoint::Gap(idx) => (gap_ids[*idx].clone(), '+'),
    }
}

fn collect_path_work(
    index: &SyngIndex,
    name: &str,
    start_info: Option<&GbwtPathStart>,
    path_len_bp: u64,
    syncmer_len_u64: u64,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<PathWork> {
    let Some(start_info) = start_info else {
        return Ok(PathWork {
            name: name.to_string(),
            gaps: Vec::new(),
            steps: Vec::new(),
            edges: Vec::new(),
            gaps_filled_with_ns: 0,
            gap_occurrences: 0,
            gap_occurrence_bp: 0,
            skipped: true,
        });
    };

    let nodes = index.walk_forward_path(start_info);
    let mut gaps: Vec<RawGap> = Vec::new();
    let mut steps: Vec<RawPathStep> = Vec::with_capacity(nodes.len());
    let mut edges: Vec<(RawEndpoint, RawEndpoint, u32)> = Vec::with_capacity(nodes.len() * 2);
    let mut gaps_filled_with_ns = 0usize;
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;

    let mut push_gap = |context: GapContext, seq: Vec<u8>, steps: &mut Vec<RawPathStep>| -> usize {
        gap_occurrences += 1;
        gap_occurrence_bp += seq.len() as u64;
        let gap_idx = gaps.len();
        gaps.push(RawGap { context, seq });
        steps.push(RawPathStep::Gap(gap_idx));
        gap_idx
    };

    // Prefix: bases [0, first_syncmer_pos) before the first syncmer.
    if let Some(&(first, first_pos)) = nodes.first() {
        if first_pos > 0 {
            let gap_seq = if let Some(seq_idx) = gap_fill {
                fetch_gap_dna(seq_idx, name, 0, first_pos)?
            } else {
                gaps_filled_with_ns += 1;
                vec![b'N'; first_pos as usize]
            };
            let gap_idx = push_gap(GapContext::Prefix { right: first }, gap_seq, &mut steps);
            edges.push((raw_gap_endpoint(gap_idx), raw_syncmer_endpoint(first), 0));
        }
        steps.push(RawPathStep::Syncmer(first));
    }

    for w in nodes.windows(2) {
        let (a, pos_a) = w[0];
        let (b, pos_b) = w[1];
        let off_u64: u64 = pos_b.saturating_sub(pos_a);
        if off_u64 <= syncmer_len_u64 {
            // Direct adjacency: overlap = syncmer_len - offset.
            let overlap = (syncmer_len_u64 - off_u64) as u32;
            edges.push((raw_syncmer_endpoint(a), raw_syncmer_endpoint(b), overlap));
            steps.push(RawPathStep::Syncmer(b));
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
            let gap_idx = push_gap(
                GapContext::Between { left: a, right: b },
                gap_seq,
                &mut steps,
            );
            edges.push((raw_syncmer_endpoint(a), raw_gap_endpoint(gap_idx), 0));
            edges.push((raw_gap_endpoint(gap_idx), raw_syncmer_endpoint(b), 0));
            steps.push(RawPathStep::Syncmer(b));
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
            let gap_idx = push_gap(GapContext::Suffix { left: last }, gap_seq, &mut steps);
            edges.push((raw_syncmer_endpoint(last), raw_gap_endpoint(gap_idx), 0));
        }
    }

    Ok(PathWork {
        name: name.to_string(),
        gaps,
        steps,
        edges,
        gaps_filled_with_ns,
        gap_occurrences,
        gap_occurrence_bp,
        skipped: false,
    })
}

fn fetch_oriented_gap_dna(
    gap_fill: Option<&UnifiedSequenceIndex>,
    seq_name: &str,
    gap_start: u64,
    gap_end: u64,
    reverse: bool,
) -> io::Result<Vec<u8>> {
    let mut seq = if let Some(seq_idx) = gap_fill {
        fetch_gap_dna(seq_idx, seq_name, gap_start, gap_end)?
    } else {
        vec![b'N'; gap_end.saturating_sub(gap_start) as usize]
    };
    if reverse {
        seq = crate::graph::reverse_complement(&seq);
    }
    Ok(seq)
}

fn collect_range_work(
    index: &SyngIndex,
    range: &SyngGfaPathRange,
    syncmer_len_u64: u64,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<PathWork> {
    let source_name = index
        .name_map
        .path_to_name
        .get(range.path_idx)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("path index {} is outside the syng name map", range.path_idx),
            )
        })?;
    let reverse = range.strand == '-';
    let fwd_nodes = index.walk_path_range(range.path_idx, range.start, range.end)?;
    let oriented_nodes: Vec<(i32, u64)> = if reverse {
        fwd_nodes
            .iter()
            .rev()
            .map(|(node, pos)| (-*node, *pos))
            .collect()
    } else {
        fwd_nodes
    };

    let mut gaps: Vec<RawGap> = Vec::new();
    let mut steps: Vec<RawPathStep> = Vec::with_capacity(oriented_nodes.len());
    let mut edges: Vec<(RawEndpoint, RawEndpoint, u32)> =
        Vec::with_capacity(oriented_nodes.len() * 2);
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;
    let mut gaps_filled_with_ns = 0usize;

    let mut push_gap = |context: GapContext, seq: Vec<u8>, steps: &mut Vec<RawPathStep>| -> usize {
        gap_occurrences += 1;
        gap_occurrence_bp += seq.len() as u64;
        if gap_fill.is_none() {
            gaps_filled_with_ns += 1;
        }
        let gap_idx = gaps.len();
        gaps.push(RawGap { context, seq });
        steps.push(RawPathStep::Gap(gap_idx));
        gap_idx
    };

    if oriented_nodes.is_empty() {
        if range.end > range.start {
            let seq =
                fetch_oriented_gap_dna(gap_fill, source_name, range.start, range.end, reverse)?;
            push_gap(
                GapContext::PathOnly {
                    path_idx: range.path_idx,
                    start: range.start,
                    end: range.end,
                    strand: range.strand,
                },
                seq,
                &mut steps,
            );
        }
        return Ok(PathWork {
            name: range.name.clone(),
            gaps,
            steps,
            edges,
            gaps_filled_with_ns,
            gap_occurrences,
            gap_occurrence_bp,
            skipped: false,
        });
    }

    if !reverse {
        if let Some(&(first, first_pos)) = oriented_nodes.first() {
            if first_pos > range.start {
                let gap_seq =
                    fetch_oriented_gap_dna(gap_fill, source_name, range.start, first_pos, false)?;
                let gap_idx = push_gap(GapContext::Prefix { right: first }, gap_seq, &mut steps);
                edges.push((raw_gap_endpoint(gap_idx), raw_syncmer_endpoint(first), 0));
            }
            steps.push(RawPathStep::Syncmer(first));
        }

        for w in oriented_nodes.windows(2) {
            let (a, pos_a) = w[0];
            let (b, pos_b) = w[1];
            let off_u64 = pos_b.saturating_sub(pos_a);
            if off_u64 <= syncmer_len_u64 {
                edges.push((
                    raw_syncmer_endpoint(a),
                    raw_syncmer_endpoint(b),
                    (syncmer_len_u64 - off_u64) as u32,
                ));
                steps.push(RawPathStep::Syncmer(b));
            } else {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    pos_a + syncmer_len_u64,
                    pos_b,
                    false,
                )?;
                let gap_idx = push_gap(
                    GapContext::Between { left: a, right: b },
                    gap_seq,
                    &mut steps,
                );
                edges.push((raw_syncmer_endpoint(a), raw_gap_endpoint(gap_idx), 0));
                edges.push((raw_gap_endpoint(gap_idx), raw_syncmer_endpoint(b), 0));
                steps.push(RawPathStep::Syncmer(b));
            }
        }

        if let Some(&(last, last_pos)) = oriented_nodes.last() {
            let suffix_start = last_pos.saturating_add(syncmer_len_u64);
            if range.end > suffix_start {
                let gap_seq =
                    fetch_oriented_gap_dna(gap_fill, source_name, suffix_start, range.end, false)?;
                let gap_idx = push_gap(GapContext::Suffix { left: last }, gap_seq, &mut steps);
                edges.push((raw_syncmer_endpoint(last), raw_gap_endpoint(gap_idx), 0));
            }
        }
    } else {
        if let Some(&(first_oriented, last_fwd_pos)) = oriented_nodes.first() {
            let prefix_fwd_start = last_fwd_pos.saturating_add(syncmer_len_u64);
            if range.end > prefix_fwd_start {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    prefix_fwd_start,
                    range.end,
                    true,
                )?;
                let gap_idx = push_gap(
                    GapContext::Prefix {
                        right: first_oriented,
                    },
                    gap_seq,
                    &mut steps,
                );
                edges.push((
                    raw_gap_endpoint(gap_idx),
                    raw_syncmer_endpoint(first_oriented),
                    0,
                ));
            }
            steps.push(RawPathStep::Syncmer(first_oriented));
        }

        for w in oriented_nodes.windows(2) {
            let (a_oriented, a_fwd_pos) = w[0];
            let (b_oriented, b_fwd_pos) = w[1];
            let off_u64 = a_fwd_pos.saturating_sub(b_fwd_pos);
            if off_u64 <= syncmer_len_u64 {
                edges.push((
                    raw_syncmer_endpoint(a_oriented),
                    raw_syncmer_endpoint(b_oriented),
                    (syncmer_len_u64 - off_u64) as u32,
                ));
                steps.push(RawPathStep::Syncmer(b_oriented));
            } else {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    b_fwd_pos + syncmer_len_u64,
                    a_fwd_pos,
                    true,
                )?;
                let gap_idx = push_gap(
                    GapContext::Between {
                        left: a_oriented,
                        right: b_oriented,
                    },
                    gap_seq,
                    &mut steps,
                );
                edges.push((
                    raw_syncmer_endpoint(a_oriented),
                    raw_gap_endpoint(gap_idx),
                    0,
                ));
                edges.push((
                    raw_gap_endpoint(gap_idx),
                    raw_syncmer_endpoint(b_oriented),
                    0,
                ));
                steps.push(RawPathStep::Syncmer(b_oriented));
            }
        }

        if let Some(&(last_oriented, first_fwd_pos)) = oriented_nodes.last() {
            if first_fwd_pos > range.start {
                let gap_seq = fetch_oriented_gap_dna(
                    gap_fill,
                    source_name,
                    range.start,
                    first_fwd_pos,
                    true,
                )?;
                let gap_idx = push_gap(
                    GapContext::Suffix {
                        left: last_oriented,
                    },
                    gap_seq,
                    &mut steps,
                );
                edges.push((
                    raw_syncmer_endpoint(last_oriented),
                    raw_gap_endpoint(gap_idx),
                    0,
                ));
            }
        }
    }

    Ok(PathWork {
        name: range.name.clone(),
        gaps,
        steps,
        edges,
        gaps_filled_with_ns,
        gap_occurrences,
        gap_occurrence_bp,
        skipped: false,
    })
}

/// Walk every forward path in `index`, gap-fill via FASTA/AGC if given,
/// and write a GFA to `writer`.
///
/// `gap_fill` selects how mid-path gaps are filled (FASTA vs. NNN).
/// Returns `(n_segments, n_links, n_paths_emitted, n_paths_skipped,
/// n_gap_segments, unique_gap_segment_bp)`.
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
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    // Silence syng's C-side debug printfs (syngBWTpathStartOld et al.)
    unsafe { syng_ffi::impg_syng_suppress_debug() };

    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    // Walk every forward path and collect:
    //   - per-path step list (syncmers + gap segments) for P/W emission
    //   - unique L edges with overlaps
    //   - in-order list of gap segments to emit later as S lines
    let mut walked_paths: Vec<WalkedPath> = Vec::with_capacity(index.name_map.path_to_name.len());
    // Edge tuple: (from_id, from_sign, to_id, to_sign, overlap_bp).
    // `String` because endpoints are integer IDs (syncmers 1..N, gaps
    // N+1, N+2, …) and we serialize them as text.
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut gap_interner = GapInterner::new(n_nodes);
    let mut n_skipped = 0usize;
    let mut gaps_filled_with_ns: usize = 0;
    let mut gap_occurrences: usize = 0;
    let mut gap_occurrence_bp: u64 = 0;

    let syncmer_len_u64 = syncmer_len as u64;
    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = (0..index.name_map.path_to_name.len())
        .into_par_iter()
        .map(|path_idx| {
            let name = &index.name_map.path_to_name[path_idx];
            let start_info = index
                .name_map
                .path_starts
                .get(path_idx)
                .and_then(|o| o.as_ref());
            let path_len_bp = index
                .name_map
                .path_to_length
                .get(path_idx)
                .copied()
                .unwrap_or(0);
            collect_path_work(
                index,
                name,
                start_info,
                path_len_bp,
                syncmer_len_u64,
                gap_fill,
            )
        })
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} path walks in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    let reduce_start = Instant::now();
    for work in path_work {
        if work.skipped {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping (rebuild the syng index to include path starts)",
                work.name
            );
            n_skipped += 1;
            continue;
        }

        gaps_filled_with_ns += work.gaps_filled_with_ns;
        gap_occurrences += work.gap_occurrences;
        gap_occurrence_bp += work.gap_occurrence_bp;

        let gap_ids: Vec<String> = work
            .gaps
            .into_iter()
            .map(|gap| gap_interner.intern(gap.context, gap.seq))
            .collect();

        for (from, to, overlap) in work.edges {
            let (from_id, from_sign) = raw_endpoint_to_gfa(&from, &gap_ids);
            let (to_id, to_sign) = raw_endpoint_to_gfa(&to, &gap_ids);
            edges.insert((from_id, from_sign, to_id, to_sign, overlap));
        }

        let steps = work
            .steps
            .into_iter()
            .map(|step| match step {
                RawPathStep::Syncmer(signed) => PathStep::Syncmer(signed),
                RawPathStep::Gap(idx) => PathStep::Gap(gap_ids[idx].clone()),
            })
            .collect();

        walked_paths.push(WalkedPath {
            name: work.name,
            steps,
        });
    }
    info!(
        "[syng2gfa] reduced path walks into {} unique gap segment(s) and {} edge(s) in {:.3}s",
        gap_interner.len(),
        edges.len(),
        reduce_start.elapsed().as_secs_f64()
    );

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} gap occurrence(s) with 'N' ({} bp over paths; {} unique segment bp). \
                 Concatenating these paths will NOT reconstruct the original genome. \
                 Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns,
            gap_occurrence_bp,
            gap_interner.total_segment_bp()
        );
    } else if !gap_interner.is_empty() {
        info!(
            "[syng2gfa] spliced {} unique real-DNA gap segment(s) from {} gap occurrence(s), {} unique bp ({} bp over paths)",
            gap_interner.len(),
            gap_occurrences,
            gap_interner.total_segment_bp(),
            gap_occurrence_bp
        );
    }

    // S lines: syncmer nodes 1..=N first, then gap segments in walk order.
    // Syncmer DNA from the C lib is lowercase; uppercase for consistency with
    // FASTA-fetched gap sequences (already uppercased upstream).
    let write_start = Instant::now();
    for node_id in 1..=n_nodes as i32 {
        let mut seq = index.syncmer_seq(node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for gap in &gap_interner.segments {
        write_segment(writer, &gap.id, &gap.seq)?;
    }
    let n_segments = n_nodes + gap_interner.len();
    let n_gap_segments = gap_interner.len();

    debug!("[syng2gfa] collected {} unique edges", edges.len());

    // L lines, sorted for deterministic output.
    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(writer, "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M")?;
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

    info!(
        "[syng2gfa] wrote {} S, {} L, {} path line(s) in {:.3}s",
        n_segments,
        n_links,
        n_paths_emitted,
        write_start.elapsed().as_secs_f64()
    );

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        n_gap_segments,
        gap_interner.total_segment_bp(),
    ))
}

/// Write syng GFA in either native overlap form or bluntg-processed form.
pub fn write_gfa_with_mode<W: Write>(
    index: &SyngIndex,
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    mode: SyngGfaMode,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    match mode {
        SyngGfaMode::Raw => write_gfa(index, writer, version, gap_fill),
        SyngGfaMode::Blunt => {
            let mut raw = Vec::new();
            let stats = write_gfa(index, &mut raw, version, gap_fill)?;
            let blunted = bluntify_gfa_bytes(&raw, version)?;
            info!(
                "[syng2gfa] bluntg processed {} raw GFA bytes into {} blunt GFA bytes",
                raw.len(),
                blunted.len()
            );
            writer.write_all(&blunted)?;
            Ok(stats)
        }
    }
}

pub fn write_range_gfa_with_mode<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
    mode: SyngGfaMode,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    match mode {
        SyngGfaMode::Raw => write_range_gfa(index, ranges, writer, version, gap_fill),
        SyngGfaMode::Blunt => {
            let mut raw = Vec::new();
            let stats = write_range_gfa(index, ranges, &mut raw, version, gap_fill)?;
            let blunt_start = Instant::now();
            let blunted = bluntify_gfa_bytes(&raw, version)?;
            info!(
                "[syng2gfa] bluntg processed {} raw range-GFA bytes into {} blunt GFA bytes in {:.3}s",
                raw.len(),
                blunted.len(),
                blunt_start.elapsed().as_secs_f64()
            );
            writer.write_all(&blunted)?;
            Ok(stats)
        }
    }
}

pub fn write_range_gfa<W: Write>(
    index: &SyngIndex,
    ranges: &[SyngGfaPathRange],
    writer: &mut W,
    version: GfaVersion,
    gap_fill: Option<&UnifiedSequenceIndex>,
) -> io::Result<(usize, usize, usize, usize, usize, u64)> {
    let syncmer_len = index.syncmer_length_bp();
    info!(
        "[syng2gfa] rendering {} selected path range(s), GFA {}, gap fill: {}",
        ranges.len(),
        version.header_tag(),
        if gap_fill.is_some() {
            "FASTA/AGC"
        } else {
            "Ns"
        }
    );

    unsafe { syng_ffi::impg_syng_suppress_debug() };
    writeln!(writer, "H\tVN:Z:{}", version.header_tag())?;

    let syncmer_len_u64 = syncmer_len as u64;
    let collect_start = Instant::now();
    let path_work: Vec<PathWork> = ranges
        .par_iter()
        .map(|range| collect_range_work(index, range, syncmer_len_u64, gap_fill))
        .collect::<io::Result<Vec<_>>>()?;
    info!(
        "[syng2gfa] collected {} selected path walks in {:.3}s using {} rayon threads",
        path_work.len(),
        collect_start.elapsed().as_secs_f64(),
        rayon::current_num_threads()
    );

    let reduce_start = Instant::now();
    let mut walked_paths: Vec<WalkedPath> = Vec::with_capacity(path_work.len());
    let mut edges: FxHashSet<(String, char, String, char, u32)> = FxHashSet::default();
    let mut gap_interner = GapInterner::new(index.num_syncmer_nodes());
    let mut used_syncmers: FxHashSet<i32> = FxHashSet::default();
    let mut n_skipped = 0usize;
    let mut gaps_filled_with_ns = 0usize;
    let mut gap_occurrences = 0usize;
    let mut gap_occurrence_bp = 0u64;

    for work in path_work {
        if work.skipped {
            warn!(
                "[syng2gfa] path '{}' has no GBWT start info — skipping",
                work.name
            );
            n_skipped += 1;
            continue;
        }

        gaps_filled_with_ns += work.gaps_filled_with_ns;
        gap_occurrences += work.gap_occurrences;
        gap_occurrence_bp += work.gap_occurrence_bp;

        let gap_ids: Vec<String> = work
            .gaps
            .into_iter()
            .map(|gap| gap_interner.intern(gap.context, gap.seq))
            .collect();

        for (from, to, overlap) in work.edges {
            if let RawEndpoint::Syncmer(node) = &from {
                used_syncmers.insert(node.abs());
            }
            if let RawEndpoint::Syncmer(node) = &to {
                used_syncmers.insert(node.abs());
            }
            let (from_id, from_sign) = raw_endpoint_to_gfa(&from, &gap_ids);
            let (to_id, to_sign) = raw_endpoint_to_gfa(&to, &gap_ids);
            edges.insert((from_id, from_sign, to_id, to_sign, overlap));
        }

        let steps = work
            .steps
            .into_iter()
            .map(|step| match step {
                RawPathStep::Syncmer(signed) => {
                    used_syncmers.insert(signed.abs());
                    PathStep::Syncmer(signed)
                }
                RawPathStep::Gap(idx) => PathStep::Gap(gap_ids[idx].clone()),
            })
            .collect();

        walked_paths.push(WalkedPath {
            name: work.name,
            steps,
        });
    }
    info!(
        "[syng2gfa] reduced selected path walks into {} syncmer node(s), {} unique gap segment(s), {} edge(s) in {:.3}s",
        used_syncmers.len(),
        gap_interner.len(),
        edges.len(),
        reduce_start.elapsed().as_secs_f64()
    );

    if gaps_filled_with_ns > 0 {
        warn!(
            "[syng2gfa] filled {} gap occurrence(s) with 'N' ({} bp over paths; {} unique segment bp). \
                 Concatenating these paths will NOT reconstruct the original genome. \
                 Pass --sequence-files <FASTA> to splice in real DNA.",
            gaps_filled_with_ns,
            gap_occurrence_bp,
            gap_interner.total_segment_bp()
        );
    } else if !gap_interner.is_empty() {
        info!(
            "[syng2gfa] spliced {} unique real-DNA gap segment(s) from {} gap occurrence(s), {} unique bp ({} bp over paths)",
            gap_interner.len(),
            gap_occurrences,
            gap_interner.total_segment_bp(),
            gap_occurrence_bp
        );
    }

    let write_start = Instant::now();
    let mut used_nodes: Vec<i32> = used_syncmers.into_iter().collect();
    used_nodes.sort_unstable();
    let n_syncmer_segments = used_nodes.len();
    for node_id in used_nodes {
        let mut seq = index.syncmer_seq(node_id);
        seq.make_ascii_uppercase();
        write_segment(writer, &node_id.to_string(), &seq)?;
    }
    for gap in &gap_interner.segments {
        write_segment(writer, &gap.id, &gap.seq)?;
    }
    let n_segments = n_syncmer_segments + gap_interner.len();
    let n_gap_segments = gap_interner.len();

    let mut edge_vec: Vec<(String, char, String, char, u32)> = edges.into_iter().collect();
    edge_vec.sort_unstable();
    for (a_id, a_sign, b_id, b_sign, overlap) in &edge_vec {
        writeln!(writer, "L\t{a_id}\t{a_sign}\t{b_id}\t{b_sign}\t{overlap}M")?;
    }
    let n_links = edge_vec.len();

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
    info!(
        "[syng2gfa] wrote selected range GFA: {} S, {} L, {} path line(s) in {:.3}s",
        n_segments,
        n_links,
        n_paths_emitted,
        write_start.elapsed().as_secs_f64()
    );

    Ok((
        n_segments,
        n_links,
        n_paths_emitted,
        n_skipped,
        n_gap_segments,
        gap_interner.total_segment_bp(),
    ))
}

fn bluntify_gfa_bytes(raw: &[u8], version: GfaVersion) -> io::Result<Vec<u8>> {
    if version != GfaVersion::V1_0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--gfa-mode blunt currently requires --gfa-version 1.0 because bluntg preserves P lines, not W lines",
        ));
    }

    let gfa = bluntg::parse_gfa(raw);
    let blunted = bluntg::bluntify_auto(gfa);
    let mut out = Vec::new();
    bluntg::write_gfa(&mut out, &blunted)?;
    Ok(out)
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
    mode: SyngGfaMode,
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

    let report = |s: usize,
                  l: usize,
                  p: usize,
                  skipped: usize,
                  gaps: usize,
                  gap_bp: u64,
                  where_: &str| {
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
            write_gfa_with_mode(&index, &mut w, version, seq_idx.as_ref(), mode)?;
        w.flush()?;
        report(
            s,
            l,
            p,
            skipped,
            gaps,
            gap_bp,
            &format!("stdout ({})", mode.label()),
        );
    } else {
        let file = std::fs::File::create(out_path)?;
        let mut w = BufWriter::new(file);
        let (s, l, p, skipped, gaps, gap_bp) =
            write_gfa_with_mode(&index, &mut w, version, seq_idx.as_ref(), mode)?;
        w.flush()?;
        report(
            s,
            l,
            p,
            skipped,
            gaps,
            gap_bp,
            &format!("'{}' ({})", out_path, mode.label()),
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_syng_gfa_mode_parse_defaults_syng_to_blunt() {
        assert_eq!(SyngGfaMode::parse("syng").unwrap(), SyngGfaMode::Blunt);
        assert_eq!(
            SyngGfaMode::parse("syng:blunt").unwrap(),
            SyngGfaMode::Blunt
        );
        assert_eq!(SyngGfaMode::parse("bluntg").unwrap(), SyngGfaMode::Blunt);
        assert_eq!(SyngGfaMode::parse("raw").unwrap(), SyngGfaMode::Raw);
        assert_eq!(SyngGfaMode::parse("syng:raw").unwrap(), SyngGfaMode::Raw);
    }

    #[test]
    fn test_bluntify_gfa_bytes_handles_variable_overlaps() {
        let raw = b"\
H\tVN:Z:1.0\n\
S\t1\tACGTAC\n\
S\t2\tTACGGA\n\
S\t3\tGGAACC\n\
L\t1\t+\t2\t+\t3M\n\
L\t2\t+\t3\t+\t1M\n\
P\tp\t1+,2+,3+\t3M,1M\n";
        let out = bluntify_gfa_bytes(raw, GfaVersion::V1_0).unwrap();
        let out = String::from_utf8(out).unwrap();
        assert!(out.contains("L\t1\t+\t2\t+\t0M"), "{out}");
        assert!(out.contains("L\t2\t+\t3\t+\t0M"), "{out}");
        assert!(out.contains("P\tp\t1+,2+,3+\t0M,0M"), "{out}");
        assert!(!out.contains("\t3M"), "{out}");
        assert!(!out.contains("\t1M"), "{out}");
    }

    #[test]
    fn test_gap_interner_reuses_only_same_sequence_and_context() {
        let mut gaps = GapInterner::new(10);

        let first = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        let same = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGT".to_vec());
        assert_eq!(first, same);
        assert_eq!(gaps.len(), 1);
        assert_eq!(gaps.total_segment_bp(), 4);

        let different_context = gaps.intern(GapContext::Prefix { right: 8 }, b"ACGT".to_vec());
        assert_ne!(first, different_context);
        assert_eq!(gaps.len(), 2);

        let different_sequence = gaps.intern(GapContext::Prefix { right: 7 }, b"ACGA".to_vec());
        assert_ne!(first, different_sequence);
        assert_eq!(gaps.len(), 3);
    }
}
