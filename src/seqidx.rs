use std::collections::HashMap;

pub struct SequenceIndex {
    name_to_id: HashMap<String, u32>,
    next_id: u32,
}

impl SequenceIndex {
    pub fn new() -> Self {
        SequenceIndex {
            name_to_id: HashMap::new(),
            next_id: 0,
        }
    }

    pub fn get_or_insert_id(&mut self, name: &str) -> u32 {
        *self.name_to_id.entry(name.to_owned()).or_insert_with(|| {
            let id = self.next_id;
            self.next_id += 1;
            id
        })
    }

    pub fn get_id(&self, name: &str) -> Option<u32> {
        self.name_to_id.get(name).copied()
    }
}
