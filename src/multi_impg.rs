// src/multi_impg.rs
//! Multi-file IMPG index implementation.
//!
//! `MultiImpg` coordinates queries across multiple per-file `.impg` indices,
//! presenting a unified view while internally managing ID translation.

use crate::forest_map::ForestMap;
use crate::impg::{AdjustedInterval, CigarOp, Impg, QueryMetadata, SortedRanges};
use crate::impg_index::ImpgIndex;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::{BasicCOITree, Interval, IntervalTree};
use log::{debug, info, warn};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::VecDeque;
use std::fs::{self, File};
use std::io::{BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::path::{Path, PathBuf};
use std::sync::{Arc, RwLock};
use std::time::SystemTime;

/// Location of a tree within a specific sub-index.
#[derive(Debug, Clone)]
struct TreeLocation {
    /// Index into `sub_indices`
    index_idx: usize,
    /// Local target_id within that sub-index
    local_target_id: u32,
}

/// Serializable version of TreeLocation for cache.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct TreeLocationSer {
    index_idx: u32,
    local_target_id: u32,
}

impl From<&TreeLocation> for TreeLocationSer {
    fn from(loc: &TreeLocation) -> Self {
        TreeLocationSer {
            index_idx: loc.index_idx as u32,
            local_target_id: loc.local_target_id,
        }
    }
}

impl From<TreeLocationSer> for TreeLocation {
    fn from(ser: TreeLocationSer) -> Self {
        TreeLocation {
            index_idx: ser.index_idx as usize,
            local_target_id: ser.local_target_id,
        }
    }
}

/// File entry for staleness detection in the cache manifest.
#[derive(Debug, Clone, Serialize, Deserialize)]
struct FileEntry {
    /// Path to the index file (relative or absolute as stored in list)
    path: String,
    /// File size in bytes
    size: u64,
    /// Modification time as duration since UNIX_EPOCH
    mtime_secs: u64,
}

/// Cache for MultiImpg unified data.
///
/// This cache stores precomputed unified sequence index, forest map, and
/// local-to-unified translation tables to speed up repeated queries with
/// the same set of per-file indices.
#[derive(Debug, Serialize, Deserialize)]
pub struct MultiImpgCache {
    /// Magic bytes for identification
    magic: [u8; 10],
    /// Version for format evolution
    version: u32,
    /// Manifest of index files for staleness detection
    manifest: Vec<FileEntry>,
    /// The alignment list file path (for auto-detection)
    list_file_path: String,
    /// Unified sequence index
    unified_seq_index: SequenceIndex,
    /// Unified forest map: target_id → list of tree locations
    unified_forest_map: Vec<(u32, Vec<TreeLocationSer>)>,
    /// Local-to-unified translation tables per index
    local_to_unified: Vec<Vec<u32>>,
}

const CACHE_MAGIC: &[u8; 10] = b"MIMPGCACH1";
const CACHE_VERSION: u32 = 1;

/// Metadata loaded from a per-file index header (seq_index + forest_map only).
struct IndexHeader {
    /// Sequence index from this file
    seq_index: SequenceIndex,
    /// Forest map from this file
    forest_map: ForestMap,
    /// Whether this index is bidirectional (V2 format)
    is_bidirectional: bool,
}

/// Multi-file IMPG index.
///
/// Coordinates queries across multiple per-file indices while presenting
/// a unified interface via the `ImpgIndex` trait.
pub struct MultiImpg {
    // ============ UNIFIED VIEW (what callers see) ============
    /// Unified sequence index: name ↔ unified_id
    pub seq_index: SequenceIndex,

    /// Unified forest map: unified_target_id → Vec<TreeLocation>
    /// Multiple indices may have trees for the same sequence
    forest_map: FxHashMap<u32, Vec<TreeLocation>>,

    // ============ PER-INDEX DATA (internal only) ============
    /// Per-index file paths
    index_paths: Vec<PathBuf>,

    /// Per-index alignment file paths
    alignment_files: Vec<String>,

    /// Per-index sequence files (if any)
    sequence_files: Vec<String>,

    /// Per-index ID translation: local_id → unified_id
    /// Uses Vec for O(1) indexed access since local IDs are dense 0..n
    local_to_unified: Vec<Vec<u32>>,

    /// Lazily-loaded sub-indices (only loaded when tree data is needed)
    sub_indices: RwLock<Vec<Option<Arc<Impg>>>>,

    /// Whether all indices are bidirectional (V2 format)
    /// True only if ALL sub-indices are V2 format
    is_bidirectional: bool,
}

impl MultiImpg {
    /// Load headers from multiple per-file indices and build unified mappings.
    ///
    /// This loads ONLY the headers (seq_index + forest_map) from each file,
    /// NOT the tree data. Trees are loaded on demand.
    pub fn load_from_files(
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
    ) -> std::io::Result<Self> {
        let num_indices = index_paths.len();
        // Load headers in parallel
        let headers: Vec<IndexHeader> = index_paths
            .par_iter()
            .map(|path| {
                Self::load_header(path).map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load header from {:?}: {}", path, e),
                    )
                })
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Build unified sequence index
        let mut unified_seq_index = SequenceIndex::new();
        let mut local_to_unified: Vec<Vec<u32>> = Vec::with_capacity(num_indices);

        for header in &headers {
            // Pre-allocate Vec with capacity for all local IDs
            let mut l2u = Vec::with_capacity(header.seq_index.len());

            for local_id in 0..header.seq_index.len() as u32 {
                if let Some(name) = header.seq_index.get_name(local_id) {
                    let len = header.seq_index.get_len_from_id(local_id);
                    let unified_id = unified_seq_index.get_or_insert_id(name, len);
                    l2u.push(unified_id);
                } else {
                    // Should not happen for valid indices, but handle gracefully
                    l2u.push(u32::MAX);
                }
            }

            local_to_unified.push(l2u);
        }

        // Build unified forest map
        let mut unified_forest_map: FxHashMap<u32, Vec<TreeLocation>> = FxHashMap::default();

        for (index_idx, header) in headers.iter().enumerate() {
            let l2u = &local_to_unified[index_idx];

            for &local_target_id in header.forest_map.entries.keys() {
                let unified_id = l2u[local_target_id as usize];
                if unified_id != u32::MAX {
                    unified_forest_map
                        .entry(unified_id)
                        .or_default()
                        .push(TreeLocation {
                            index_idx,
                            local_target_id,
                        });
                }
            }
        }

        // Check if all indices are bidirectional
        let all_bidirectional = headers.iter().all(|h| h.is_bidirectional);
        if !all_bidirectional {
            let v1_count = headers.iter().filter(|h| !h.is_bidirectional).count();
            warn!(
                "{} of {} indices are V1 (unidirectional). Rebuild with default settings for full bidirectional support.",
                v1_count, num_indices
            );
        }

        info!(
            "Built unified index with {} sequences and {} targets (bidirectional: {})",
            unified_seq_index.len(),
            unified_forest_map.len(),
            all_bidirectional
        );

        Ok(Self {
            seq_index: unified_seq_index,
            forest_map: unified_forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            is_bidirectional: all_bidirectional,
        })
    }

    /// Load MultiImpg from a cache file if valid, or build from scratch.
    ///
    /// Auto-detects cache file as `{list_file}.multi_impg` and validates
    /// staleness before using. If cache is stale or missing, builds from
    /// scratch and saves a new cache.
    pub fn load_with_cache(
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
        list_file: &Path,
    ) -> std::io::Result<Self> {
        let cache_path = list_file.with_extension("multi_impg");

        // Try to load from cache
        if cache_path.exists() {
            match MultiImpgCache::load(&cache_path) {
                Ok(cache) => {
                    if cache.is_valid(index_paths, list_file)? {
                        return Self::from_cache(
                            cache,
                            index_paths,
                            alignment_files,
                            sequence_files,
                        );
                    } else {
                        info!("Cache stale, rebuilding...");
                    }
                }
                Err(e) => {
                    warn!(
                        "Failed to load cache {:?}: {}, rebuilding...",
                        cache_path, e
                    );
                }
            }
        }

        // Build from scratch
        let multi = Self::load_from_files(index_paths, alignment_files, sequence_files)?;

        // Save cache for next time
        if let Err(e) = multi.save_cache(&cache_path, list_file) {
            warn!("Failed to save cache {:?}: {}", cache_path, e);
        }

        Ok(multi)
    }

    /// Create MultiImpg from a validated cache.
    fn from_cache(
        cache: MultiImpgCache,
        index_paths: &[PathBuf],
        alignment_files: &[String],
        sequence_files: Option<&[String]>,
    ) -> std::io::Result<Self> {
        let num_indices = index_paths.len();

        // Convert serialized forest map back to FxHashMap<u32, Vec<TreeLocation>>
        let forest_map: FxHashMap<u32, Vec<TreeLocation>> = cache
            .unified_forest_map
            .into_iter()
            .map(|(target_id, locs)| {
                (
                    target_id,
                    locs.into_iter().map(TreeLocation::from).collect(),
                )
            })
            .collect();

        info!(
            "Loaded {} sequences and {} targets from cache",
            cache.unified_seq_index.len(),
            forest_map.len()
        );

        // Note: When loading from cache, we assume bidirectional=true since
        // caches are only created for new-format indices. If users have
        // V1 indices, the cache will be invalidated and rebuilt.
        Ok(Self {
            seq_index: cache.unified_seq_index,
            forest_map,
            index_paths: index_paths.to_vec(),
            alignment_files: alignment_files.to_vec(),
            sequence_files: sequence_files.map(|s| s.to_vec()).unwrap_or_default(),
            local_to_unified: cache.local_to_unified,
            sub_indices: RwLock::new(vec![None; num_indices]),
            is_bidirectional: true,
        })
    }

    /// Save the unified index data to a cache file.
    pub fn save_cache(&self, cache_path: &Path, list_file: &Path) -> std::io::Result<()> {
        // Build manifest from index files
        let manifest: Vec<FileEntry> = self
            .index_paths
            .iter()
            .map(|path| {
                let metadata = fs::metadata(path)?;
                let mtime = metadata
                    .modified()?
                    .duration_since(SystemTime::UNIX_EPOCH)
                    .map_err(std::io::Error::other)?;
                Ok(FileEntry {
                    path: path.to_string_lossy().to_string(),
                    size: metadata.len(),
                    mtime_secs: mtime.as_secs(),
                })
            })
            .collect::<std::io::Result<Vec<_>>>()?;

        // Convert forest_map to serializable format
        let unified_forest_map: Vec<(u32, Vec<TreeLocationSer>)> = self
            .forest_map
            .iter()
            .map(|(&target_id, locs)| (target_id, locs.iter().map(TreeLocationSer::from).collect()))
            .collect();

        let cache = MultiImpgCache {
            magic: *CACHE_MAGIC,
            version: CACHE_VERSION,
            manifest,
            list_file_path: list_file.to_string_lossy().to_string(),
            unified_seq_index: self.seq_index.clone(),
            unified_forest_map,
            local_to_unified: self.local_to_unified.clone(),
        };

        let file = File::create(cache_path)?;
        let mut writer = BufWriter::new(file);

        bincode::serde::encode_into_std_write(&cache, &mut writer, bincode::config::standard())
            .map_err(std::io::Error::other)?;

        writer.flush()?;
        Ok(())
    }

    /// Load only the header (seq_index + forest_map) from a single index file.
    fn load_header(path: &Path) -> std::io::Result<IndexHeader> {
        const MAGIC_V1: &[u8] = b"IMPGIDX1";
        const MAGIC_V2: &[u8] = b"IMPGIDX2";

        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        // Read and verify magic bytes (support both V1 and V2)
        let mut magic_buf = [0u8; 8];
        reader.read_exact(&mut magic_buf)?;
        if magic_buf != MAGIC_V1 && magic_buf != MAGIC_V2 {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Invalid magic bytes in {:?}", path),
            ));
        }

        // Read forest map offset
        let mut offset_buf = [0u8; 8];
        reader.read_exact(&mut offset_buf)?;
        let forest_map_offset = u64::from_le_bytes(offset_buf);

        // Read sequence index
        let seq_index: SequenceIndex =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load sequence index: {e}"),
                    )
                })?;

        // Seek to forest map and read it
        reader.seek(SeekFrom::Start(forest_map_offset))?;
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to load forest map: {e}"),
                    )
                })?;

        let is_bidirectional = magic_buf == MAGIC_V2;
        Ok(IndexHeader {
            seq_index,
            forest_map,
            is_bidirectional,
        })
    }

    /// Get or load a sub-index.
    fn get_sub_index(&self, index_idx: usize) -> std::io::Result<Arc<Impg>> {
        // Fast path: check if already loaded
        {
            let indices = self.sub_indices.read().unwrap();
            if let Some(ref impg) = indices[index_idx] {
                return Ok(Arc::clone(impg));
            }
        }

        // Slow path: load the index
        let path = &self.index_paths[index_idx];
        let alignment_files = vec![self.alignment_files[index_idx].clone()];
        let seq_files = if self.sequence_files.is_empty() {
            None
        } else {
            Some(self.sequence_files.as_slice())
        };

        let file = File::open(path)?;
        let reader = BufReader::new(file);
        let impg = Impg::load_from_file(
            reader,
            &alignment_files,
            path.to_string_lossy().to_string(),
            seq_files,
        )?;
        let impg = Arc::new(impg);

        // Store in cache
        {
            let mut indices = self.sub_indices.write().unwrap();
            indices[index_idx] = Some(Arc::clone(&impg));
        }

        Ok(impg)
    }

    /// Translate an AdjustedInterval from local IDs to unified IDs.
    fn translate_to_unified(
        &self,
        interval: AdjustedInterval,
        index_idx: usize,
    ) -> Option<AdjustedInterval> {
        let (query_interval, cigar, target_interval) = interval;
        let l2u = &self.local_to_unified[index_idx];

        // Vec indexing with bounds check
        let unified_query_id = *l2u.get(query_interval.metadata as usize)?;
        let unified_target_id = *l2u.get(target_interval.metadata as usize)?;

        // Check for invalid sentinel values
        if unified_query_id == u32::MAX || unified_target_id == u32::MAX {
            return None;
        }

        Some((
            Interval {
                first: query_interval.first,
                last: query_interval.last,
                metadata: unified_query_id,
            },
            cigar,
            Interval {
                first: target_interval.first,
                last: target_interval.last,
                metadata: unified_target_id,
            },
        ))
    }

    /// Query all sub-indices that have trees for the given target.
    fn query_all_indices(
        &self,
        unified_target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval> {
        let locations = match self.forest_map.get(&unified_target_id) {
            Some(locs) => locs,
            None => {
                return vec![self.make_self_interval(
                    unified_target_id,
                    range_start,
                    range_end,
                    store_cigar,
                )]
            }
        };

        // Query all relevant sub-indices in parallel
        let results: Vec<Vec<AdjustedInterval>> = locations
            .par_iter()
            .filter_map(|loc| {
                let impg = match self.get_sub_index(loc.index_idx) {
                    Ok(i) => i,
                    Err(e) => {
                        warn!("Failed to load sub-index {}: {}", loc.index_idx, e);
                        return None;
                    }
                };

                // Query using local target ID
                let local_results = impg.query(
                    loc.local_target_id,
                    range_start,
                    range_end,
                    store_cigar,
                    min_gap_compressed_identity,
                    sequence_index,
                    approximate_mode,
                );

                // Translate results to unified IDs
                let unified_results: Vec<AdjustedInterval> = local_results
                    .into_iter()
                    .filter_map(|r| self.translate_to_unified(r, loc.index_idx))
                    .collect();

                Some(unified_results)
            })
            .collect();

        // Merge all results, but only include ONE self-interval
        let mut final_results = Vec::new();
        let mut seen_self = false;

        for result_set in results {
            for result in result_set {
                // Check if this is a self-interval (query == target == input)
                let is_self = result.0.metadata == unified_target_id
                    && result.2.metadata == unified_target_id
                    && result.0.first == range_start
                    && result.0.last == range_end;

                if is_self {
                    if !seen_self {
                        final_results.push(result);
                        seen_self = true;
                    }
                    // Skip duplicate self-intervals
                } else {
                    final_results.push(result);
                }
            }
        }

        // If we never saw a self-interval, add one
        if !seen_self {
            final_results.insert(
                0,
                self.make_self_interval(unified_target_id, range_start, range_end, store_cigar),
            );
        }

        // Sort results for deterministic ordering (excluding the self-interval which should stay first)
        // Sort by: query_id, query_start, query_end, target_start, target_end
        if final_results.len() > 1 {
            let self_interval = final_results.remove(0);
            final_results.sort_by(|a, b| {
                let a_key = (a.0.metadata, a.0.first, a.0.last, a.2.first, a.2.last);
                let b_key = (b.0.metadata, b.0.first, b.0.last, b.2.first, b.2.last);
                a_key.cmp(&b_key)
            });
            final_results.insert(0, self_interval);
        }

        final_results
    }

    /// Create a self-referential interval for the query region.
    fn make_self_interval(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
    ) -> AdjustedInterval {
        (
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
            if store_cigar {
                vec![CigarOp::new(range_end - range_start, '=')]
            } else {
                Vec::new()
            },
            Interval {
                first: range_start,
                last: range_end,
                metadata: target_id,
            },
        )
    }
}

impl ImpgIndex for MultiImpg {
    fn seq_index(&self) -> &SequenceIndex {
        &self.seq_index
    }

    fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval> {
        self.query_all_indices(
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
        )
    }

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        _cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
        // For MultiImpg, we don't use the shared cache since each sub-index
        // has its own file offsets. Just do a normal query.
        self.query(
            target_id,
            range_start,
            range_end,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            false, // approximate_mode
        )
    }

    fn populate_cigar_cache(
        &self,
        _target_id: u32,
        _range_start: i32,
        _range_end: i32,
        _min_gap_compressed_identity: Option<f64>,
        _sequence_index: Option<&UnifiedSequenceIndex>,
        _cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
        // For MultiImpg, CIGAR caching is not implemented since each sub-index
        // has different file offsets. This is a no-op.
    }

    fn query_transitive_dfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        // Transitive query implementation using DFS with deterministic ordering
        self.transitive_query_impl(
            target_id,
            range_start,
            range_end,
            masked_regions,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
            true, // use_dfs
        )
    }

    fn query_transitive_bfs(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
    ) -> Vec<AdjustedInterval> {
        // Transitive query implementation using BFS with deterministic ordering
        self.transitive_query_impl(
            target_id,
            range_start,
            range_end,
            masked_regions,
            max_depth,
            min_transitive_len,
            min_distance_between_ranges,
            min_output_length,
            store_cigar,
            min_gap_compressed_identity,
            sequence_index,
            approximate_mode,
            subset_filter,
            false, // use_dfs
        )
    }

    fn get_or_load_tree(&self, target_id: u32) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        // For MultiImpg, we can't return a single tree since multiple indices
        // may have data for the same target. Return None and let callers use query().
        // This is only used by stats and similarity commands which may need adaptation.

        // Try to get a tree from any sub-index that has this target
        let locations = self.forest_map.get(&target_id)?;

        // Get the first location and try to load its tree
        let loc = locations.first()?;
        let local_target_id = loc.local_target_id;

        let impg = self.get_sub_index(loc.index_idx).ok()?;
        impg.get_or_load_tree(local_target_id)
    }

    fn target_ids(&self) -> Vec<u32> {
        self.forest_map.keys().copied().collect()
    }

    fn remove_cached_tree(&self, _target_id: u32) {
        // For MultiImpg, we don't cache trees at this level.
        // Sub-indices manage their own tree caches.
    }

    fn num_targets(&self) -> usize {
        self.forest_map.len()
    }

    fn sequence_files(&self) -> &[String] {
        &self.sequence_files
    }

    fn query_reverse_for_depth(&self, query_id: u32) -> Vec<(i32, i32, u32)> {
        let mut results = Vec::new();

        // Iterate through all target_ids in the forest map
        for &target_id in self.forest_map.keys() {
            // Skip if querying self
            if target_id == query_id {
                continue;
            }

            if let Some(tree) = self.get_or_load_tree(target_id) {
                for interval in tree.iter() {
                    if interval.metadata.query_id() == query_id {
                        let query_start = interval.metadata.query_start();
                        let query_end = interval.metadata.query_end();
                        results.push((query_start, query_end, target_id));
                    }
                }
            }
        }

        results
    }

    fn build_query_to_targets_map(&self) -> FxHashMap<u32, Vec<u32>> {
        let mut query_to_targets: FxHashMap<u32, Vec<u32>> = FxHashMap::default();

        for &target_id in self.forest_map.keys() {
            if let Some(tree) = self.get_or_load_tree(target_id) {
                let mut seen_queries: rustc_hash::FxHashSet<u32> = rustc_hash::FxHashSet::default();
                for interval in tree.iter() {
                    let qid = interval.metadata.query_id();
                    if qid != target_id && seen_queries.insert(qid) {
                        query_to_targets.entry(qid).or_default().push(target_id);
                    }
                }
            }
        }

        // Clear tree cache after building the map
        self.clear_tree_cache();

        query_to_targets
    }

    fn query_reverse_for_depth_with_map(
        &self,
        query_id: u32,
        query_to_targets: &FxHashMap<u32, Vec<u32>>,
    ) -> Vec<(i32, i32, u32)> {
        let mut results = Vec::new();

        if let Some(target_ids) = query_to_targets.get(&query_id) {
            for &target_id in target_ids {
                if let Some(tree) = self.get_or_load_tree(target_id) {
                    for interval in tree.iter() {
                        if interval.metadata.query_id() == query_id {
                            let query_start = interval.metadata.query_start();
                            let query_end = interval.metadata.query_end();
                            results.push((query_start, query_end, target_id));
                        }
                    }
                }
            }
        }

        results
    }

    fn clear_tree_cache(&self) {
        // For MultiImpg, clear caches of all loaded sub-indices
        let indices = self.sub_indices.read().unwrap();
        for sub_index in indices.iter().flatten() {
            sub_index.clear_tree_cache();
        }
    }

    fn is_bidirectional(&self) -> bool {
        self.is_bidirectional
    }
}

impl MultiImpg {
    /// Internal implementation of transitive queries.
    ///
    /// Matches the behavior of `Impg::query_transitive_dfs` and `Impg::query_transitive_bfs`,
    /// including min_distance_between_ranges checks and stack merging.
    fn transitive_query_impl(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        masked_regions: Option<&FxHashMap<u32, SortedRanges>>,
        max_depth: u16,
        min_transitive_len: i32,
        min_distance_between_ranges: i32,
        min_output_length: Option<i32>,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
        subset_filter: Option<&SubsetFilter>,
        use_dfs: bool,
    ) -> Vec<AdjustedInterval> {
        // Initialize visited ranges
        let mut visited_ranges: FxHashMap<u32, SortedRanges> = if let Some(m) = masked_regions {
            m.iter().map(|(&k, v)| (k, v.clone())).collect()
        } else {
            (0..self.seq_index.len() as u32)
                .into_par_iter()
                .map(|id| {
                    let len = self.seq_index.get_len_from_id(id).unwrap_or(0);
                    (id, SortedRanges::new(len as i32, 0))
                })
                .collect()
        };

        // Filter input range
        let filtered_input_range = visited_ranges
            .entry(target_id)
            .or_default()
            .insert((range_start, range_end));

        let mut results = Vec::new();
        let mut stack: VecDeque<(u32, i32, i32, u16)> = VecDeque::new();

        // Add filtered input ranges
        for (filtered_start, filtered_end) in filtered_input_range {
            results.push(self.make_self_interval(
                target_id,
                filtered_start,
                filtered_end,
                store_cigar,
            ));

            if (filtered_start - filtered_end).abs() >= min_transitive_len {
                stack.push_back((target_id, filtered_start, filtered_end, 0));
            }
        }

        while let Some((current_target_id, current_start, current_end, current_depth)) = if use_dfs
        {
            stack.pop_back() // DFS: pop from back (LIFO) - O(1)
        } else {
            stack.pop_front() // BFS: pop from front (FIFO) - O(1) with VecDeque
        } {
            if max_depth > 0 && current_depth >= max_depth {
                continue;
            }

            debug!(
                "Transitive query: {}:{}-{}, depth={}",
                self.seq_index.get_name(current_target_id).unwrap_or("?"),
                current_start,
                current_end,
                current_depth
            );

            // Query all sub-indices for this region
            let step_results = self.query_all_indices(
                current_target_id,
                current_start,
                current_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            );

            // Process results - matches original Impg behavior
            for result in step_results {
                let query_id = result.0.metadata;

                // Skip self-referential results
                if query_id == current_target_id {
                    continue;
                }

                // Apply subset filter if provided (always keep original target)
                if let Some(filter) = subset_filter {
                    if query_id != target_id {
                        if let Some(name) = self.seq_index.get_name(query_id) {
                            if !filter.matches(name) {
                                continue;
                            }
                        }
                    }
                }

                // Normalize coordinates
                let (adjusted_query_start, adjusted_query_end) = if result.0.first <= result.0.last
                {
                    (result.0.first, result.0.last)
                } else {
                    (result.0.last, result.0.first)
                };

                let length = (result.0.last - result.0.first).abs();

                // Add to results only if it passes min_output_length filter
                let should_add_to_output = if let Some(min_len) = min_output_length {
                    length >= min_len
                } else {
                    true
                };
                if should_add_to_output {
                    results.push(result);
                }

                // Only add non-overlapping portions to the stack for further exploration
                let ranges = visited_ranges.entry(query_id).or_insert_with(|| {
                    let len = self.seq_index.get_len_from_id(query_id).unwrap_or(0);
                    SortedRanges::new(len as i32, 0)
                });

                let mut should_add = true;

                // Check if the range is too close to any existing ranges
                if min_distance_between_ranges > 0 {
                    let (new_min, new_max) = (adjusted_query_start, adjusted_query_end);

                    // Find insertion point in sorted ranges
                    let idx = match ranges
                        .ranges
                        .binary_search_by_key(&new_min, |&(start, _)| start)
                    {
                        Ok(i) => i,
                        Err(i) => i,
                    };

                    // Only need to check adjacent ranges due to sorting
                    if idx > 0 {
                        // Check previous range
                        let (_, prev_end) = ranges.ranges[idx - 1];
                        if (new_min - prev_end).abs() < min_distance_between_ranges {
                            should_add = false;
                        }
                    }
                    if idx < ranges.ranges.len() {
                        // Check next range
                        let (next_start, _) = ranges.ranges[idx];
                        if (next_start - new_max).abs() < min_distance_between_ranges {
                            should_add = false;
                        }
                    }
                }

                if should_add {
                    let new_ranges = ranges.insert((adjusted_query_start, adjusted_query_end));

                    // Add non-overlapping portions to stack
                    for (new_start, new_end) in new_ranges {
                        if (new_end - new_start).abs() >= min_transitive_len {
                            stack.push_back((query_id, new_start, new_end, current_depth + 1));
                        }
                    }
                }
            }

            // Merge contiguous/overlapping ranges with same sequence_id (matches original)
            let slice = stack.make_contiguous();
            slice.sort_by_key(|(id, start, _, _)| (*id, *start));

            let mut write = 0;
            for read in 1..slice.len() {
                if slice[write].0 == slice[read].0 &&   // Same sequence_id
                    slice[write].2 >= slice[read].1
                // Overlapping or contiguous
                {
                    // Merge by extending end
                    slice[write].2 = slice[write].2.max(slice[read].2);
                } else {
                    write += 1;
                    slice.swap(write, read);
                }
            }
            if !stack.is_empty() {
                stack.truncate(write + 1);
            }
        }

        results
    }
}

impl MultiImpgCache {
    /// Load a cache from disk.
    pub fn load(path: &Path) -> std::io::Result<Self> {
        let file = File::open(path)?;
        let mut reader = BufReader::new(file);

        let cache: MultiImpgCache =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to decode cache: {e}"),
                    )
                })?;

        // Verify magic and version
        if &cache.magic != CACHE_MAGIC {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                "Invalid cache magic bytes",
            ));
        }
        if cache.version != CACHE_VERSION {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "Unsupported cache version: {} (expected {})",
                    cache.version, CACHE_VERSION
                ),
            ));
        }

        Ok(cache)
    }

    /// Check if the cache is still valid (not stale).
    ///
    /// The cache is valid if:
    /// 1. All listed index files still exist
    /// 2. No index file has been modified (mtime + size check)
    /// 3. The list of files matches exactly
    pub fn is_valid(&self, index_paths: &[PathBuf], _list_file: &Path) -> std::io::Result<bool> {
        // Check file list length matches
        if self.manifest.len() != index_paths.len() {
            debug!(
                "Cache invalid: manifest has {} files, but {} index paths provided",
                self.manifest.len(),
                index_paths.len()
            );
            return Ok(false);
        }

        // Check each file's path, mtime, and size
        for (entry, path) in self.manifest.iter().zip(index_paths) {
            let path_str = path.to_string_lossy().to_string();
            if entry.path != path_str {
                debug!(
                    "Cache invalid: path mismatch '{}' vs '{}'",
                    entry.path, path_str
                );
                return Ok(false);
            }

            let metadata = match fs::metadata(path) {
                Ok(m) => m,
                Err(_) => {
                    debug!("Cache invalid: file not found '{}'", path_str);
                    return Ok(false);
                }
            };

            if metadata.len() != entry.size {
                debug!(
                    "Cache invalid: size mismatch for '{}' ({} vs {})",
                    path_str,
                    metadata.len(),
                    entry.size
                );
                return Ok(false);
            }

            let mtime = metadata
                .modified()?
                .duration_since(SystemTime::UNIX_EPOCH)
                .map_err(std::io::Error::other)?;

            if mtime.as_secs() != entry.mtime_secs {
                debug!(
                    "Cache invalid: mtime mismatch for '{}' ({} vs {})",
                    path_str,
                    mtime.as_secs(),
                    entry.mtime_secs
                );
                return Ok(false);
            }
        }

        Ok(true)
    }
}
