use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, BufWriter};

/// Forest map entry: target_id -> tree_offset in the serialized IMPG index
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForestMapEntry {
    pub target_id: u32,
    pub tree_offset: u64,
}

/// Forest map structure containing mappings from target_id to tree offsets
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ForestMap {
    pub entries: FxHashMap<u32, u64>,
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

    /// Check if the forest map contains a target_id
    pub fn contains_target(&self, target_id: u32) -> bool {
        self.entries.contains_key(&target_id)
    }

    /// Get the number of entries in the forest map
    pub fn len(&self) -> usize {
        self.entries.len()
    }

    /// Check if the forest map is empty
    pub fn is_empty(&self) -> bool {
        self.entries.is_empty()
    }

    /// Save the forest map to a file
    pub fn save_to_file(&self, file_path: &str) -> std::io::Result<()> {
        let file = File::create(file_path)?;
        let mut writer = BufWriter::new(file);
        bincode::serde::encode_into_std_write(self, &mut writer, bincode::config::standard())
            .map_err(|e| std::io::Error::other(format!("Failed to serialize forest map: {:?}", e)))?;
        Ok(())
    }

    /// Load the forest map from a file
    pub fn load_from_file(file_path: &str) -> std::io::Result<Self> {
        let file = File::open(file_path)?;
        let mut reader = BufReader::new(file);
        let forest_map: ForestMap =
            bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
                .map_err(|e| {
                    std::io::Error::new(
                        std::io::ErrorKind::InvalidData,
                        format!("Failed to deserialize forest map: {:?}", e),
                    )
                })?;
        Ok(forest_map)
    }

    /// Generate forest map filename based on the index filename
    pub fn get_forest_map_filename(index_filename: &str) -> String {
        format!("{}.map", index_filename)
    }
}

impl Default for ForestMap {
    fn default() -> Self {
        Self::new()
    }
}