use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};

/// Forest map entry: target_id -> tree_offset in the serialized IMPG index
#[derive(Serialize, Deserialize)]
struct ForestMapEntry {
    target_id: u32,
    tree_offset: u64,
}

/// Forest map structure containing mappings from target_id to tree offsets
#[derive(Serialize, Deserialize)]
pub struct ForestMap {
    pub entries: FxHashMap<u32, u64>,
}

impl Default for ForestMap {
    fn default() -> Self {
        Self::new()
    }
}

impl ForestMap {
    /// Create a new empty forest map
    pub fn new() -> Self {
        Self {
            entries: FxHashMap::default(),
        }
    }

    /// Add a target_id -> tree_offset mapping
    pub fn add_entry(&mut self, target_id: u32, tree_offset: u64) {
        self.entries.insert(target_id, tree_offset);
    }

    /// Get the tree offset for a target_id
    pub fn get_tree_offset(&self, target_id: u32) -> Option<u64> {
        self.entries.get(&target_id).copied()
    }
}
