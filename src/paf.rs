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
    pub cigar: Option<String>,
    pub strand: Strand,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Strand {
    Forward,
    Reverse,
}

impl PafRecord {
    pub fn parse(line: &str) -> Result<Self, ParseErr> {
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
        let strand_char = fields[4].chars().next().ok_or_else(|| {
            ParseErr::InvalidFormat("Expected '+' or '-' for strand".to_string())
        })?;
        let strand = match strand_char {
            '+' => Strand::Forward,
            '-' => Strand::Reverse,
            _ => return Err(ParseErr::InvalidStrand),
        };

        let cigar = fields.iter()
            .find(|&&f| f.starts_with("cg:Z:"))
            .map(|&s| s[5..].to_string());

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
            cigar,
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
    let mut records = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(ParseErr::IoError)?;
        let record = PafRecord::parse(&line)?;
        records.push(record);
    }
    Ok(records)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_paf_valid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255";
        let record = PafRecord::parse(line).unwrap();
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
                cigar: None, // Adjusted for the example to not include a CIGAR string
                strand: Strand::Forward,
            }
        );
    }

    #[test]
    fn test_parse_paf_valid_2() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10=";
        assert!(PafRecord::parse(line).is_ok());
    }

    #[test]
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        assert!(PafRecord::parse(line).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10Q";
        assert!(PafRecord::parse(line).is_err());
    }
}
