use crate::seqidx::SequenceIndex;
use std::io::{BufRead, Error as IoError};
use std::num::ParseIntError;

#[derive(Debug, PartialEq)]
pub struct PartialPafRecord {
    pub query_id: u32,
    pub query_start: usize,
    pub query_end: usize,
    pub target_id: u32,
    pub target_start: usize,
    pub target_end: usize,
    pub cigar_offset: u64, // Track cigar offset
    pub cigar_bytes: usize,
}

#[derive(Default, PartialEq, Clone, Copy, Debug)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

impl PartialPafRecord {
    pub fn strand(&self) -> Strand {
        if self.query_start <= self.query_end {
            Strand::Forward
        } else {
            Strand::Reverse
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
        // Adjust query coordinates based on strand to encode strand information
        let (adjusted_query_start, adjusted_query_end) = match strand_char {
            '+' => (query_start, query_end),
            '-' => (query_end, query_start), // Swap to encode reverse strand
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

        // Create the record
        let record = Self {
            query_id,
            query_start: adjusted_query_start,
            query_end: adjusted_query_end,
            target_id,
            target_start,
            target_end,
            cigar_offset,
            cigar_bytes,
        };

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
                cigar_offset: (line.len() + 1) as u64,
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

    #[test]
    fn test_strand_encoding() {
        let mut seq_index = SequenceIndex::new();
        
        // Test forward strand - query coordinates should remain in order
        let forward_line = "seq1\t100\t10\t20\t+\tseq2\t100\t30\t40\t10\t20\t255";
        let forward_record = PartialPafRecord::parse(forward_line, 0, &mut seq_index).unwrap();
        assert_eq!(forward_record.strand(), Strand::Forward);
        assert_eq!(forward_record.query_start, 10);
        assert_eq!(forward_record.query_end, 20);
        assert!(forward_record.query_start <= forward_record.query_end);
        
        // Test reverse strand - query coordinates should be swapped  
        let reverse_line = "seq3\t100\t10\t20\t-\tseq4\t100\t30\t40\t10\t20\t255";
        let reverse_record = PartialPafRecord::parse(reverse_line, 0, &mut seq_index).unwrap();
        assert_eq!(reverse_record.strand(), Strand::Reverse);
        assert_eq!(reverse_record.query_start, 20);
        assert_eq!(reverse_record.query_end, 10);
        assert!(reverse_record.query_start > reverse_record.query_end);
    }
}
