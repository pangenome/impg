use serde::{Deserialize, Serialize};
use std::io::{BufRead, Error as IoError};
use std::num::ParseIntError;

#[derive(Debug, PartialEq, Clone)]
pub struct PafRecord {
    pub query_name: String,
    pub query_length: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub target_name: String,
    pub target_length: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub strand: Strand,
    pub cigar_offset: u64,
    pub cigar_bytes: usize,
}

#[derive(Default, Debug, PartialEq, Clone, Copy, Serialize, Deserialize)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

impl PafRecord {
    pub fn parse(line: &str, file_pos: u64) -> Result<Self, ParseErr> {
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

        Ok(Self {
            query_name,
            query_length,
            query_start,
            query_end,
            target_name,
            target_length,
            target_start,
            target_end,
            strand,
            cigar_offset,
            cigar_bytes,
        })
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

pub fn parse_paf<R: BufRead>(reader: R) -> Result<Vec<PafRecord>, ParseErr> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(ParseErr::IoError)?;
        let record = PafRecord::parse(&line, bytes_read)?;
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
        let record = PafRecord::parse(line, 0).unwrap();
        assert_eq!(
            record,
            PafRecord {
                query_name: "seq1".to_string(),
                query_length: 100,
                query_start: 0,
                query_end: 100,
                target_name: "seq2".to_string(),
                target_length: 100,
                target_start: 0,
                target_end: 100,
                strand: Strand::Forward,
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
        assert!(PafRecord::parse(line, 0).is_ok());
    }

    #[test]
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        assert!(PafRecord::parse(line, 0).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10Q";
        assert!(PafRecord::parse(line, 0).is_err());
    }
}
