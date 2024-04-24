use std::collections::HashMap;
use serde::{Serialize, Deserialize};

#[derive(Debug, Clone, Serialize, Deserialize, Default)]
pub struct SequenceIndex {
    name_to_id: HashMap<String, u32>,
    id_to_name: HashMap<u32, String>,
    id_to_len: HashMap<u32, usize>,
    next_id: u32,
}

impl SequenceIndex {
    pub fn new() -> Self {
        SequenceIndex {
            name_to_id: HashMap::new(),
            id_to_name: HashMap::new(),
            id_to_len: HashMap::new(),
            next_id: 0,
        }
    }

    pub fn get_or_insert_id(&mut self, name: &str, length: Option<usize>) -> u32 {
        let id = *self.name_to_id.entry(name.to_owned()).or_insert_with(|| {
            let id = self.next_id;
            self.id_to_name.insert(id, name.to_owned());
            self.next_id += 1;
            id
        });

        if let Some(len) = length {
            self.id_to_len.entry(id).or_insert(len);
        }

        id
    }

    pub fn get_id(&self, name: &str) -> Option<u32> {
        self.name_to_id.get(name).copied()
    }

    pub fn get_name(&self, id: u32) -> Option<&str> {
        self.id_to_name.get(&id).map(|s| s.as_str())
    }

    pub fn get_len_from_id(&self, id: u32) -> Option<usize> {
        self.id_to_len.get(&id).copied()
    }

    pub fn is_empty(&self) -> bool {
        self.name_to_id.is_empty()
    }
    
    pub fn len(&self) -> usize {
        self.name_to_id.len()
    }
}
