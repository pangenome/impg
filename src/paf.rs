use crate::seqidx::SequenceIndex;
use std::io::{BufRead, Error as IoError};
use std::num::ParseIntError;

pub struct PartialPafRecord {
    pub query_id: u32,
    pub query_start: usize,
    pub query_end: usize,
    pub target_id: u32,
    pub target_start: usize,
    pub target_end: usize,
    pub strand_and_cigar_offset: u64, // Track strand and cigar offset
    pub cigar_bytes: usize,
}

#[derive(Default, PartialEq, Clone, Copy)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

impl PartialPafRecord {
    const STRAND_BIT: u64 = 0x8000000000000000; // Most significant bit for u64

    pub fn strand(&self) -> Strand {
        if (self.strand_and_cigar_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
    pub fn set_strand(&mut self, strand: Strand) {
        match strand {
            Strand::Forward => self.strand_and_cigar_offset &= !Self::STRAND_BIT,
            Strand::Reverse => self.strand_and_cigar_offset |= Self::STRAND_BIT,
        }
    }

    pub fn parse(
        line: &str,
        file_pos: u64,
        seq_index: &mut SequenceIndex,
    ) -> Result<Self, ParseErr> {
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(ParseErr::NotEnoughFields);
        }

        let query_name = fields[0].to_string();
        let query_length = fields[1].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let query_start = fields[2].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let query_end = fields[3].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let target_name = fields[5].to_string();
        let target_length = fields[6].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let target_start = fields[7].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let target_end = fields[8].parse::<usize>().map_err(ParseErr::InvalidField)?;
        let strand_char = fields[4]
            .chars()
            .next()
            .ok_or_else(|| ParseErr::InvalidFormat("Expected '+' or '-' for strand".to_string()))?;
        let strand = match strand_char {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => return Err(ParseErr::InvalidStrand),
        };

        // Convert names to IDs using the SequenceIndex
        let query_id = seq_index.get_or_insert_id(&query_name, Some(query_length));
        let target_id = seq_index.get_or_insert_id(&target_name, Some(target_length));

        let mut cigar_offset: u64 = file_pos;
        let mut cigar_bytes: usize = 0;

        for tag_str in fields.iter() {
            if tag_str.starts_with("cg:Z:") {
                cigar_offset += 5;
                cigar_bytes = tag_str.len() - 5;
                break;
            } else {
                cigar_offset += (tag_str.len() + 1) as u64;
            }
        }

        // Create the record and set strand
        let mut record = Self {
            query_id,
            query_start,
            query_end,
            target_id,
            target_start,
            target_end,
            strand_and_cigar_offset: cigar_offset,
            cigar_bytes,
        };
        record.set_strand(strand);

        Ok(record)
    }
}

#[derive(Debug)]
pub enum ParseErr {
    NotEnoughFields,
    IoError(IoError),
    InvalidField(ParseIntError),
    InvalidStrand,
    InvalidCigarFormat,
    UnsupportedCigarOperation,
    InvalidFormat(String),
}

pub fn parse_paf<R: BufRead>(
    reader: R,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<PartialPafRecord>, ParseErr> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(ParseErr::IoError)?;
        let record = PartialPafRecord::parse(&line, bytes_read, seq_index)?;
        records.push(record);

        // Size of line plus newline
        bytes_read += (line.len() + 1) as u64;
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_paf_valid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255";
        let mut seq_index = SequenceIndex::new();
        let record = PartialPafRecord::parse(line, 0, &mut seq_index).unwrap();

        // IDs should be 0 and 1 as they're the first entries in the SequenceIndex
        let query_id = seq_index.get_id("seq1").unwrap();
        let target_id = seq_index.get_id("seq2").unwrap();

        assert_eq!(
            record,
            PartialPafRecord {
                query_id,
                query_start: 0,
                query_end: 100,
                target_id,
                target_start: 0,
                target_end: 100,
                // If no cigar, then the offset is just the length of the line and cigar_bytes=0
                // Should we use Option<> instead?
                strand_and_cigar_offset: (line.len() + 1) as u64,
                cigar_bytes: 0,
            }
        );
    }

    #[test]
    fn test_parse_paf_valid_2() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10=";
        let mut seq_index = SequenceIndex::new();
        assert!(PartialPafRecord::parse(line, 0, &mut seq_index).is_ok());
    }

    #[test]
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        let mut seq_index = SequenceIndex::new();
        assert!(PartialPafRecord::parse(line, 0, &mut seq_index).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10Q";
        let mut seq_index = SequenceIndex::new();
        assert!(PartialPafRecord::parse(line, 0, &mut seq_index).is_err());
    }
}
