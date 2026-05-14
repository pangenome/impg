//! Parallel construction helpers for syng indices.
//!
//! The first implemented stage builds a deterministic global syncmer
//! dictionary in parallel. Existing serial GBWT path construction can then
//! replay sequences through that frozen dictionary.

use rayon::prelude::*;
use std::collections::BTreeMap;
use std::io;

use crate::syng::{
    packed_syncmer_word_len, GbwtPathStart, PackedSyncmer, SyngIndex, SyngNameMap,
    SyncmerParams,
};
use crate::syng_ffi;

fn base_code(base: u8) -> u8 {
    match base {
        b'a' | b'A' | 0 => 0,
        b'c' | b'C' | 1 => 1,
        b'g' | b'G' | 2 => 2,
        b't' | b'T' | 3 => 3,
        _ => 0,
    }
}

fn complement_code(code: u8) -> u8 {
    3 - code
}

fn forward_is_canonical(syncmer: &[u8]) -> bool {
    if syncmer.is_empty() {
        return true;
    }
    let mut left = 0usize;
    let mut right = syncmer.len() - 1;
    while left <= right {
        let fwd = base_code(syncmer[left]);
        let rev = complement_code(base_code(syncmer[right]));
        if fwd != rev {
            return fwd < rev;
        }
        if left == right {
            break;
        }
        left += 1;
        right -= 1;
    }
    true
}

/// Convert a sequence to syng's numeric DNA encoding.
pub fn numeric_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| base_code(b)).collect()
}

/// Pack a full syncmer in the canonical orientation expected by KmerHash.
pub fn pack_canonical_syncmer_with_orientation(syncmer: &[u8]) -> (PackedSyncmer, bool) {
    let word_len = (syncmer.len() + 31) / 32;
    let mut bytes = vec![0u8; word_len * 8];
    let use_forward = forward_is_canonical(syncmer);

    for i in 0..syncmer.len() {
        let code = if use_forward {
            base_code(syncmer[i])
        } else {
            complement_code(base_code(syncmer[syncmer.len() - 1 - i]))
        };
        bytes[i / 4] |= code << (2 * (i % 4));
    }

    let packed = bytes
        .chunks_exact(8)
        .map(|chunk| {
            let mut word = [0u8; 8];
            word.copy_from_slice(chunk);
            u64::from_le_bytes(word)
        })
        .collect::<Vec<_>>();
    (packed, use_forward)
}

/// Pack a full syncmer in the canonical orientation expected by KmerHash.
pub fn pack_canonical_syncmer(syncmer: &[u8]) -> PackedSyncmer {
    pack_canonical_syncmer_with_orientation(syncmer).0
}

/// Extract packed canonical full-syncmers from a sequence.
pub fn extract_packed_syncmers(
    params: SyncmerParams,
    seq: &[u8],
) -> io::Result<Vec<PackedSyncmer>> {
    let syncmer_len = (params.w + params.k) as usize;
    if seq.len() < syncmer_len {
        return Ok(Vec::new());
    }

    let mut seq_buf = numeric_sequence(seq);
    seq_buf.push(0);

    let mut out = Vec::new();
    unsafe {
        let seqhash = syng_ffi::impg_seqhashCreateSafe(
            params.k as i32,
            params.w as i32,
            params.seed as i32,
        );
        if seqhash.is_null() {
            return Err(io::Error::other("seqhashCreate returned null"));
        }

        let sit = syng_ffi::syncmerIterator(
            seqhash,
            seq_buf.as_mut_ptr() as *mut i8,
            seq.len() as i32,
        );
        if sit.is_null() {
            syng_ffi::impg_seqhashDestroy(seqhash);
            return Err(io::Error::other("syncmerIterator returned null"));
        }

        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(
            sit,
            std::ptr::null_mut(),
            &mut pos,
            std::ptr::null_mut(),
        ) {
            let start = pos as usize;
            if start + syncmer_len > seq.len() {
                syng_ffi::impg_seqhashIteratorDestroy(sit);
                syng_ffi::impg_seqhashDestroy(seqhash);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "syncmer iterator returned out-of-range position {} for sequence length {}",
                        start,
                        seq.len()
                    ),
                ));
            }
            out.push(pack_canonical_syncmer(
                &seq_buf[start..start + syncmer_len],
            ));
        }

        syng_ffi::impg_seqhashIteratorDestroy(sit);
        syng_ffi::impg_seqhashDestroy(seqhash);
    }

    Ok(out)
}

pub fn sort_dedup_packed_syncmers(mut syncmers: Vec<PackedSyncmer>) -> Vec<PackedSyncmer> {
    syncmers.sort_unstable();
    syncmers.dedup();
    syncmers
}

/// Build a deterministic global packed-syncmer dictionary in parallel.
pub fn build_packed_syncmer_dictionary(
    params: SyncmerParams,
    sequences: &[(String, Vec<u8>)],
) -> io::Result<Vec<PackedSyncmer>> {
    let expected_words = packed_syncmer_word_len(params);
    let syncmers = sequences
        .par_iter()
        .map(|(_, seq)| extract_packed_syncmers(params, seq))
        .try_reduce(Vec::new, |mut acc, mut part| {
            acc.append(&mut part);
            Ok(acc)
        })?;

    let dictionary = sort_dedup_packed_syncmers(syncmers);
    if let Some((i, packed)) = dictionary
        .iter()
        .enumerate()
        .find(|(_, packed)| packed.len() != expected_words)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "dictionary entry {} has {} words, expected {}",
                i,
                packed.len(),
                expected_words
            ),
        ));
    }
    Ok(dictionary)
}

#[derive(Debug, Clone)]
pub struct SequenceSyncmerPath {
    pub name: String,
    pub sequence_len: u64,
    pub syncmers: Vec<(i32, u64)>,
}

#[derive(Debug, Clone, Copy)]
pub struct ParallelGbwtBuildStats {
    pub sequences: usize,
    pub indexed_sequences: usize,
    pub gbwt_paths: usize,
    pub syncmer_occurrences: usize,
    pub reduced_nodes: usize,
    pub dictionary_nodes: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct EdgeKey {
    sync: i32,
    offset: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
struct DirectedEdge {
    from: i32,
    to: i32,
    offset: u32,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Side {
    In,
    Out,
}

#[derive(Debug, Clone)]
struct ReducedGbwtPath {
    source_path_idx: Option<usize>,
    nodes: Vec<i32>,
    positions: Vec<u64>,
}

#[derive(Debug, Clone)]
struct NodeOccurrence {
    abs_node: u32,
    path_order: usize,
    node_idx: usize,
    incoming_edge: DirectedEdge,
    incoming_rank: u32,
    incoming_side: Side,
    incoming_key: EdgeKey,
    sequence_side: Side,
    value_key: EdgeKey,
}

#[derive(Debug, Clone)]
struct SequenceItem {
    position_key: EdgeKey,
    position_rank: u32,
    value_key: EdgeKey,
}

struct NodeReduced {
    node_id: u32,
    in_edges: Vec<syng_ffi::ImpgSyngEdge>,
    in_run_sym: Vec<i64>,
    in_run_len: Vec<i64>,
    out_edges: Vec<syng_ffi::ImpgSyngEdge>,
    out_run_sym: Vec<i64>,
    out_run_len: Vec<i64>,
}

pub fn extract_syncmer_path_with_dictionary(
    params: SyncmerParams,
    seq: &[u8],
    dictionary: &[PackedSyncmer],
) -> io::Result<Vec<(i32, u64)>> {
    let syncmer_len = (params.w + params.k) as usize;
    if seq.len() < syncmer_len {
        return Ok(Vec::new());
    }

    let mut seq_buf = numeric_sequence(seq);
    seq_buf.push(0);

    let mut out = Vec::new();
    unsafe {
        let seqhash = syng_ffi::impg_seqhashCreateSafe(
            params.k as i32,
            params.w as i32,
            params.seed as i32,
        );
        if seqhash.is_null() {
            return Err(io::Error::other("seqhashCreate returned null"));
        }

        let sit = syng_ffi::syncmerIterator(
            seqhash,
            seq_buf.as_mut_ptr() as *mut i8,
            seq.len() as i32,
        );
        if sit.is_null() {
            syng_ffi::impg_seqhashDestroy(seqhash);
            return Err(io::Error::other("syncmerIterator returned null"));
        }

        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(
            sit,
            std::ptr::null_mut(),
            &mut pos,
            std::ptr::null_mut(),
        ) {
            let start = pos as usize;
            if start + syncmer_len > seq.len() {
                syng_ffi::impg_seqhashIteratorDestroy(sit);
                syng_ffi::impg_seqhashDestroy(seqhash);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "syncmer iterator returned out-of-range position {} for sequence length {}",
                        start,
                        seq.len()
                    ),
                ));
            }
            let (packed, is_forward) =
                pack_canonical_syncmer_with_orientation(&seq_buf[start..start + syncmer_len]);
            let dict_idx = dictionary.binary_search(&packed).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("syncmer at position {} is absent from the packed dictionary", start),
                )
            })?;
            let node_id = i32::try_from(dict_idx + 1).map_err(|_| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    "syncmer dictionary exceeds i32 node id range",
                )
            })?;
            let signed_node = if is_forward { node_id } else { -node_id };
            out.push((signed_node, start as u64));
        }

        syng_ffi::impg_seqhashIteratorDestroy(sit);
        syng_ffi::impg_seqhashDestroy(seqhash);
    }

    Ok(out)
}

pub fn extract_sequence_syncmer_paths(
    params: SyncmerParams,
    sequences: &[(String, Vec<u8>)],
    dictionary: &[PackedSyncmer],
) -> io::Result<Vec<SequenceSyncmerPath>> {
    sequences
        .par_iter()
        .map(|(name, seq)| {
            Ok(SequenceSyncmerPath {
                name: name.clone(),
                sequence_len: seq.len() as u64,
                syncmers: extract_syncmer_path_with_dictionary(params, seq, dictionary)?,
            })
        })
        .collect()
}

pub fn build_index_with_parallel_gbwt_reduce(
    params: SyncmerParams,
    dictionary: Vec<PackedSyncmer>,
    sequence_paths: Vec<SequenceSyncmerPath>,
    position_sample_shift: u32,
    position_sample_seed: u64,
) -> io::Result<(SyngIndex, ParallelGbwtBuildStats)> {
    let mut name_map = SyngNameMap::new();
    let mut gbwt_paths = Vec::new();
    let mut sampled_paths: Vec<Vec<(i64, i32)>> = Vec::with_capacity(sequence_paths.len());

    for path in &sequence_paths {
        let path_idx = name_map.add(path.name.clone(), path.sequence_len) as usize;
        let sampled = path
            .syncmers
            .iter()
            .map(|&(node, pos)| {
                let pos_i32 = i32::try_from(pos).map_err(|_| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("syncmer position {} exceeds i32 range", pos),
                    )
                })?;
                Ok((node as i64, pos_i32))
            })
            .collect::<io::Result<Vec<_>>>()?;
        sampled_paths.push(sampled);

        if path.syncmers.is_empty() {
            let _ = path_idx;
            continue;
        }

        let nodes: Vec<i32> = path.syncmers.iter().map(|&(node, _)| node).collect();
        let positions: Vec<u64> = path.syncmers.iter().map(|&(_, pos)| pos).collect();
        gbwt_paths.push(ReducedGbwtPath {
            source_path_idx: Some(path_idx),
            nodes: nodes.clone(),
            positions: positions.clone(),
        });

        let rc_nodes: Vec<i32> = nodes.iter().rev().map(|node| -*node).collect();
        let mut rc_positions = Vec::with_capacity(positions.len());
        let mut rc_pos = 0u64;
        rc_positions.push(rc_pos);
        for i in (1..positions.len()).rev() {
            rc_pos = rc_pos
                .checked_add(positions[i] - positions[i - 1])
                .ok_or_else(|| {
                    io::Error::new(
                        io::ErrorKind::InvalidData,
                        "reverse path syncmer coordinate overflow",
                    )
                })?;
            rc_positions.push(rc_pos);
        }
        gbwt_paths.push(ReducedGbwtPath {
            source_path_idx: None,
            nodes: rc_nodes,
            positions: rc_positions,
        });
    }

    let mut occurrences = collect_node_occurrences(&gbwt_paths)?;
    assign_incoming_ranks(&mut occurrences)?;
    install_path_starts(&gbwt_paths, &occurrences, &mut name_map)?;

    let reduced_nodes = reduce_nodes(&occurrences)?;
    let gbwt = unsafe {
        syng_ffi::syngBWTcreate(
            (params.w + params.k) as i32,
            i64::try_from(dictionary.len() + 1).map_err(|_| {
                io::Error::new(io::ErrorKind::InvalidData, "dictionary length exceeds i64 range")
            })?,
        )
    };
    if gbwt.is_null() {
        return Err(io::Error::other("syngBWTcreate returned null"));
    }

    let start_counts = start_counts_by_node(&gbwt_paths);
    unsafe {
        for (start_node, count) in start_counts {
            syng_ffi::impg_syngBWTstartCountAdd(gbwt, start_node, count);
        }
        for node in &reduced_nodes {
            syng_ffi::impg_syngBWTsetNode(
                gbwt,
                node.node_id as i32,
                node.in_edges.as_ptr(),
                node.in_edges.len() as i32,
                node.in_run_sym.as_ptr(),
                node.in_run_len.as_ptr(),
                node.in_run_sym.len() as i64,
                node.out_edges.as_ptr(),
                node.out_edges.len() as i32,
                node.out_run_sym.as_ptr(),
                node.out_run_len.as_ptr(),
                node.out_run_sym.len() as i64,
            );
        }
        syng_ffi::impg_syngBWTlocBuild(gbwt);
    }

    let mut index = SyngIndex::from_prebuilt_gbwt(params, &dictionary, gbwt, name_map)?;
    index.set_sampled_positions_from_paths(
        position_sample_shift,
        position_sample_seed,
        sampled_paths
            .iter()
            .enumerate()
            .map(|(path_idx, syncmers)| {
                (path_idx, sequence_paths[path_idx].sequence_len, syncmers.as_slice())
            }),
    )?;

    let stats = ParallelGbwtBuildStats {
        sequences: sequence_paths.len(),
        indexed_sequences: sequence_paths
            .iter()
            .filter(|path| !path.syncmers.is_empty())
            .count(),
        gbwt_paths: gbwt_paths.len(),
        syncmer_occurrences: occurrences.len(),
        reduced_nodes: reduced_nodes.len(),
        dictionary_nodes: dictionary.len(),
    };
    Ok((index, stats))
}

fn collect_node_occurrences(paths: &[ReducedGbwtPath]) -> io::Result<Vec<NodeOccurrence>> {
    let mut occurrences = Vec::new();
    for (path_order, path) in paths.iter().enumerate() {
        for node_idx in 0..path.nodes.len() {
            let current = path.nodes[node_idx];
            let abs_node = current.unsigned_abs();
            let (prev, in_offset) = if node_idx == 0 {
                (0, 0)
            } else {
                let delta = path.positions[node_idx]
                    .checked_sub(path.positions[node_idx - 1])
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "syncmer positions are not monotonic in GBWT path",
                        )
                    })?;
                (path.nodes[node_idx - 1], checked_offset(delta)?)
            };
            let (next, out_offset) = if node_idx + 1 == path.nodes.len() {
                (0, 0)
            } else {
                let delta = path.positions[node_idx + 1]
                    .checked_sub(path.positions[node_idx])
                    .ok_or_else(|| {
                        io::Error::new(
                            io::ErrorKind::InvalidData,
                            "syncmer positions are not monotonic in GBWT path",
                        )
                    })?;
                (path.nodes[node_idx + 1], checked_offset(delta)?)
            };

            let (incoming_side, incoming_key, sequence_side, value_key) = if current >= 0 {
                (
                    Side::In,
                    EdgeKey {
                        sync: prev,
                        offset: in_offset,
                    },
                    Side::Out,
                    EdgeKey {
                        sync: next,
                        offset: out_offset,
                    },
                )
            } else {
                (
                    Side::Out,
                    EdgeKey {
                        sync: -prev,
                        offset: in_offset,
                    },
                    Side::In,
                    EdgeKey {
                        sync: -next,
                        offset: out_offset,
                    },
                )
            };

            occurrences.push(NodeOccurrence {
                abs_node,
                path_order,
                node_idx,
                incoming_edge: DirectedEdge {
                    from: prev,
                    to: current,
                    offset: in_offset,
                },
                incoming_rank: 0,
                incoming_side,
                incoming_key,
                sequence_side,
                value_key,
            });
        }
    }
    Ok(occurrences)
}

fn checked_offset(delta: u64) -> io::Result<u32> {
    u32::try_from(delta).map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("syncmer offset {} exceeds u32 range", delta),
        )
    })
}

fn assign_incoming_ranks(occurrences: &mut [NodeOccurrence]) -> io::Result<()> {
    let mut order: Vec<(DirectedEdge, usize, usize, usize)> = occurrences
        .iter()
        .enumerate()
        .map(|(idx, occ)| (occ.incoming_edge, occ.path_order, occ.node_idx, idx))
        .collect();
    order.sort_unstable();

    let mut current_edge: Option<DirectedEdge> = None;
    let mut rank = 0u32;
    for (edge, _, _, idx) in order {
        if current_edge != Some(edge) {
            current_edge = Some(edge);
            rank = 0;
        }
        occurrences[idx].incoming_rank = rank;
        rank = rank.checked_add(1).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "directed edge occurrence rank exceeds u32 range",
            )
        })?;
    }
    Ok(())
}

fn install_path_starts(
    paths: &[ReducedGbwtPath],
    occurrences: &[NodeOccurrence],
    name_map: &mut SyngNameMap,
) -> io::Result<()> {
    for occ in occurrences.iter().filter(|occ| occ.node_idx == 0) {
        if let Some(path_idx) = paths[occ.path_order].source_path_idx {
            let path = &paths[occ.path_order];
            name_map.set_path_start(
                path_idx as u32,
                GbwtPathStart {
                    start_node: path.nodes[0],
                    start_count: occ.incoming_rank,
                    num_syncmers: path.nodes.len() as u32,
                    first_syncmer_pos: path.positions[0],
                },
            );
        }
    }
    Ok(())
}

fn start_counts_by_node(paths: &[ReducedGbwtPath]) -> Vec<(i32, u32)> {
    let mut counts = BTreeMap::<i32, u32>::new();
    for path in paths {
        if let Some(&start) = path.nodes.first() {
            *counts.entry(start).or_default() += 1;
        }
    }
    counts.into_iter().collect()
}

fn reduce_nodes(occurrences: &[NodeOccurrence]) -> io::Result<Vec<NodeReduced>> {
    let mut order: Vec<usize> = (0..occurrences.len()).collect();
    order.sort_unstable_by_key(|&idx| occurrences[idx].abs_node);

    let mut groups = Vec::new();
    let mut start = 0usize;
    while start < order.len() {
        let node_id = occurrences[order[start]].abs_node;
        let mut end = start + 1;
        while end < order.len() && occurrences[order[end]].abs_node == node_id {
            end += 1;
        }
        groups.push((node_id, start, end));
        start = end;
    }

    groups
        .par_iter()
        .map(|&(node_id, start, end)| reduce_one_node(node_id, &order[start..end], occurrences))
        .collect()
}

fn reduce_one_node(
    node_id: u32,
    indices: &[usize],
    occurrences: &[NodeOccurrence],
) -> io::Result<NodeReduced> {
    let mut in_counts = BTreeMap::<EdgeKey, u32>::new();
    let mut out_counts = BTreeMap::<EdgeKey, u32>::new();
    let mut in_items = Vec::new();
    let mut out_items = Vec::new();

    for &idx in indices {
        let occ = &occurrences[idx];
        match occ.incoming_side {
            Side::In => *in_counts.entry(occ.incoming_key).or_default() += 1,
            Side::Out => *out_counts.entry(occ.incoming_key).or_default() += 1,
        }

        let item = SequenceItem {
            position_key: occ.incoming_key,
            position_rank: occ.incoming_rank,
            value_key: occ.value_key,
        };
        match occ.sequence_side {
            Side::In => in_items.push(item),
            Side::Out => out_items.push(item),
        }
    }

    let (in_edges, in_run_sym, in_run_len) =
        build_side(&in_counts, &out_counts, &in_items)?;
    let (out_edges, out_run_sym, out_run_len) =
        build_side(&out_counts, &in_counts, &out_items)?;

    Ok(NodeReduced {
        node_id,
        in_edges,
        in_run_sym,
        in_run_len,
        out_edges,
        out_run_sym,
        out_run_len,
    })
}

fn build_side(
    side_counts: &BTreeMap<EdgeKey, u32>,
    order_counts: &BTreeMap<EdgeKey, u32>,
    items: &[SequenceItem],
) -> io::Result<(Vec<syng_ffi::ImpgSyngEdge>, Vec<i64>, Vec<i64>)> {
    let edges: Vec<syng_ffi::ImpgSyngEdge> = side_counts
        .iter()
        .map(|(&key, &count)| syng_ffi::ImpgSyngEdge {
            sync: key.sync,
            offset: key.offset,
            count,
        })
        .collect();
    if edges.is_empty() {
        return Ok((edges, Vec::new(), Vec::new()));
    }

    let mut value_to_index = BTreeMap::new();
    for (idx, &key) in side_counts.keys().enumerate() {
        value_to_index.insert(key, idx as i64);
    }

    let mut order_prefix = BTreeMap::new();
    let mut total_len = 0usize;
    for (&key, &count) in order_counts {
        order_prefix.insert(key, total_len);
        total_len = total_len
            .checked_add(count as usize)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "side length overflow"))?;
    }

    if total_len != items.len() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "side sequence length mismatch: order counts imply {}, got {} items",
                total_len,
                items.len()
            ),
        ));
    }

    let mut sequence = vec![None; total_len];
    for item in items {
        let prefix = *order_prefix.get(&item.position_key).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("missing order edge {:?}", item.position_key),
            )
        })?;
        let pos = prefix + item.position_rank as usize;
        let slot = sequence.get_mut(pos).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("side sequence position {} is out of range {}", pos, total_len),
            )
        })?;
        if slot.is_some() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("duplicate side sequence position {}", pos),
            ));
        }
        *slot = Some(item.value_key);
    }

    let mut run_sym = Vec::new();
    let mut run_len = Vec::new();
    let mut last_sym = None;
    for value in sequence {
        let value = value.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidData, "unfilled side sequence position")
        })?;
        let sym = *value_to_index.get(&value).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("sequence edge {:?} is absent from side directory", value),
            )
        })?;
        if last_sym == Some(sym) {
            *run_len.last_mut().expect("run length exists") += 1;
        } else {
            run_sym.push(sym);
            run_len.push(1);
            last_sym = Some(sym);
        }
    }

    Ok((edges, run_sym, run_len))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_packing_is_orientation_invariant() {
        let seq = b"ACGTACGTACGTA";
        let rc = b"TACGTACGTACGT";
        assert_eq!(
            pack_canonical_syncmer(seq),
            pack_canonical_syncmer(rc),
            "forward and reverse-complement syncmers should share packed canonical form"
        );
    }
}
