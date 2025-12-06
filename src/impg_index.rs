// src/impg_index.rs
//! Trait abstraction for IMPG index operations.
//!
//! This trait allows both single-file `Impg` and multi-file `MultiImpg` to provide
//! the same interface to query commands, making the multi-index logic invisible
//! to callers.

use crate::impg::{AdjustedInterval, CigarOp, Impg, QueryMetadata, SortedRanges};
use crate::multi_impg::MultiImpg;
use crate::seqidx::SequenceIndex;
use crate::sequence_index::UnifiedSequenceIndex;
use crate::subset_filter::SubsetFilter;
use coitrees::BasicCOITree;
use rustc_hash::FxHashMap;
use std::sync::Arc;

/// Trait for IMPG index operations.
///
/// Both `Impg` (single-file) and `MultiImpg` (multi-file) implement this trait,
/// allowing query commands to work transparently with either.
pub trait ImpgIndex: Send + Sync {
    /// Get a reference to the sequence index (unified for MultiImpg).
    fn seq_index(&self) -> &SequenceIndex;

    /// Query a single region without transitivity.
    fn query(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        approximate_mode: bool,
    ) -> Vec<AdjustedInterval>;

    /// Query with pre-populated CIGAR cache for efficiency.
    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval>;

    /// Populate CIGAR cache for a region.
    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    );

    /// Transitive query using depth-first search.
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
    ) -> Vec<AdjustedInterval>;

    /// Transitive query using breadth-first search.
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
    ) -> Vec<AdjustedInterval>;

    /// Get or load an interval tree for a target sequence.
    fn get_or_load_tree(&self, target_id: u32)
        -> Option<Arc<BasicCOITree<QueryMetadata, u32>>>;

    /// Get all target IDs that have interval trees (for iteration).
    fn target_ids(&self) -> Vec<u32>;

    /// Remove a cached tree from memory (for memory management).
    fn remove_cached_tree(&self, target_id: u32);

    /// Get number of targets with trees.
    fn num_targets(&self) -> usize {
        self.target_ids().len()
    }

    /// Get the sequence files (FASTA/AGC) associated with this index.
    fn sequence_files(&self) -> &[String];
}

/// Enum wrapper that can hold either a single `Impg` or a `MultiImpg`.
///
/// This allows the CLI to work with either index type without changing
/// function signatures throughout the codebase.
pub enum ImpgWrapper {
    Single(Impg),
    Multi(MultiImpg),
}

impl ImpgWrapper {
    /// Create a wrapper from a single Impg
    pub fn from_single(impg: Impg) -> Self {
        ImpgWrapper::Single(impg)
    }

    /// Create a wrapper from a MultiImpg
    pub fn from_multi(multi: MultiImpg) -> Self {
        ImpgWrapper::Multi(multi)
    }

    /// Get a reference to the inner Impg if this is a Single variant.
    /// Returns None for Multi variant.
    pub fn as_single(&self) -> Option<&Impg> {
        match self {
            ImpgWrapper::Single(impg) => Some(impg),
            ImpgWrapper::Multi(_) => None,
        }
    }
}

impl ImpgIndex for ImpgWrapper {
    fn seq_index(&self) -> &SequenceIndex {
        match self {
            ImpgWrapper::Single(impg) => impg.seq_index(),
            ImpgWrapper::Multi(multi) => multi.seq_index(),
        }
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
        match self {
            ImpgWrapper::Single(impg) => impg.query(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            ),
            ImpgWrapper::Multi(multi) => multi.query(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                approximate_mode,
            ),
        }
    }

    fn query_with_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        store_cigar: bool,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cigar_cache: &FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) -> Vec<AdjustedInterval> {
        match self {
            ImpgWrapper::Single(impg) => impg.query_with_cache(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                cigar_cache,
            ),
            ImpgWrapper::Multi(multi) => multi.query_with_cache(
                target_id,
                range_start,
                range_end,
                store_cigar,
                min_gap_compressed_identity,
                sequence_index,
                cigar_cache,
            ),
        }
    }

    fn populate_cigar_cache(
        &self,
        target_id: u32,
        range_start: i32,
        range_end: i32,
        min_gap_compressed_identity: Option<f64>,
        sequence_index: Option<&UnifiedSequenceIndex>,
        cache: &mut FxHashMap<(u32, u64), Vec<CigarOp>>,
    ) {
        match self {
            ImpgWrapper::Single(impg) => impg.populate_cigar_cache(
                target_id,
                range_start,
                range_end,
                min_gap_compressed_identity,
                sequence_index,
                cache,
            ),
            ImpgWrapper::Multi(multi) => multi.populate_cigar_cache(
                target_id,
                range_start,
                range_end,
                min_gap_compressed_identity,
                sequence_index,
                cache,
            ),
        }
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
        match self {
            ImpgWrapper::Single(impg) => impg.query_transitive_dfs(
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
            ),
            ImpgWrapper::Multi(multi) => multi.query_transitive_dfs(
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
            ),
        }
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
        match self {
            ImpgWrapper::Single(impg) => impg.query_transitive_bfs(
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
            ),
            ImpgWrapper::Multi(multi) => multi.query_transitive_bfs(
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
            ),
        }
    }

    fn get_or_load_tree(
        &self,
        target_id: u32,
    ) -> Option<Arc<BasicCOITree<QueryMetadata, u32>>> {
        match self {
            ImpgWrapper::Single(impg) => impg.get_or_load_tree(target_id),
            ImpgWrapper::Multi(multi) => multi.get_or_load_tree(target_id),
        }
    }

    fn target_ids(&self) -> Vec<u32> {
        match self {
            ImpgWrapper::Single(impg) => impg.target_ids(),
            ImpgWrapper::Multi(multi) => multi.target_ids(),
        }
    }

    fn remove_cached_tree(&self, target_id: u32) {
        match self {
            ImpgWrapper::Single(impg) => impg.remove_cached_tree(target_id),
            ImpgWrapper::Multi(multi) => multi.remove_cached_tree(target_id),
        }
    }

    fn num_targets(&self) -> usize {
        match self {
            ImpgWrapper::Single(impg) => impg.num_targets(),
            ImpgWrapper::Multi(multi) => multi.num_targets(),
        }
    }

    fn sequence_files(&self) -> &[String] {
        match self {
            ImpgWrapper::Single(impg) => impg.sequence_files(),
            ImpgWrapper::Multi(multi) => multi.sequence_files(),
        }
    }
}
