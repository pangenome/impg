use serde::{Deserialize, Serialize};
use std::io;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct SourceSequenceId(pub u32);

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PanSn {
    pub sample: String,
    pub haplotype: String,
    pub contig: String,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct PathIdentity {
    pub full_name: String,
    pub pansn: Option<PanSn>,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceSequenceRecord {
    pub id: SourceSequenceId,
    pub name: String,
    pub length: u64,
    pub identity: PathIdentity,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct SourceInterval {
    pub source_sequence_id: SourceSequenceId,
    pub start: u64,
    pub end: u64,
    pub strand: char,
}

#[derive(Debug, Clone, Default, PartialEq, Eq, Serialize, Deserialize)]
pub struct SequenceNamespace {
    pub sequences: Vec<SourceSequenceRecord>,
}

impl PanSn {
    pub fn parse(name: &str) -> Option<Self> {
        let mut parts = name.split('#');
        let first = parts.next()?;
        let second = parts.next()?;
        let third = parts.next()?;
        if parts.next().is_some() || first.is_empty() || second.is_empty() || third.is_empty() {
            return None;
        }
        Some(Self {
            sample: first.to_string(),
            haplotype: second.to_string(),
            contig: third.to_string(),
        })
    }
}

impl PathIdentity {
    pub fn from_name(name: impl Into<String>) -> Self {
        let full_name = name.into();
        let pansn = PanSn::parse(&full_name);
        Self { full_name, pansn }
    }
}

impl SourceInterval {
    pub fn new(
        source_sequence_id: SourceSequenceId,
        start: u64,
        end: u64,
        strand: char,
    ) -> io::Result<Self> {
        if end < start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("source interval end {end} is before start {start}"),
            ));
        }
        if !matches!(strand, '+' | '-') {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("source interval strand must be '+' or '-', got '{strand}'"),
            ));
        }
        Ok(Self {
            source_sequence_id,
            start,
            end,
            strand,
        })
    }
}

impl SequenceNamespace {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn add_sequence(&mut self, name: impl Into<String>, length: u64) -> SourceSequenceId {
        let name = name.into();
        if let Some(existing) = self.sequences.iter().find(|record| record.name == name) {
            return existing.id;
        }
        let id = SourceSequenceId(self.sequences.len() as u32);
        self.sequences.push(SourceSequenceRecord {
            id,
            identity: PathIdentity::from_name(name.clone()),
            name,
            length,
        });
        id
    }

    pub fn get(&self, id: SourceSequenceId) -> Option<&SourceSequenceRecord> {
        self.sequences.get(id.0 as usize)
    }

    pub fn get_by_name(&self, name: &str) -> Option<&SourceSequenceRecord> {
        self.sequences.iter().find(|record| record.name == name)
    }

    pub fn from_named_lengths<I, N>(records: I) -> Self
    where
        I: IntoIterator<Item = (N, u64)>,
        N: Into<String>,
    {
        let mut namespace = Self::new();
        for (name, length) in records {
            namespace.add_sequence(name, length);
        }
        namespace
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn pansn_parses_three_part_names() {
        let pansn = PanSn::parse("HG002#1#chr6").unwrap();
        assert_eq!(pansn.sample, "HG002");
        assert_eq!(pansn.haplotype, "1");
        assert_eq!(pansn.contig, "chr6");
    }

    #[test]
    fn two_part_names_are_not_silently_parsed_as_pansn() {
        assert!(PanSn::parse("GRCh38#chr6").is_none());
    }

    #[test]
    fn non_pansn_names_are_still_source_sequences() {
        let mut namespace = SequenceNamespace::new();
        let id = namespace.add_sequence("fragment_001", 123);
        let record = namespace.get(id).unwrap();
        assert_eq!(record.name, "fragment_001");
        assert!(record.identity.pansn.is_none());
    }

    #[test]
    fn source_interval_validates_coordinates_and_strand() {
        let id = SourceSequenceId(0);
        assert!(SourceInterval::new(id, 0, 10, '+').is_ok());
        assert!(SourceInterval::new(id, 10, 0, '+').is_err());
        assert!(SourceInterval::new(id, 0, 10, '?').is_err());
    }
}
