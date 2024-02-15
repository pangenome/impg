use std::io::BufRead;
use std::num::ParseIntError;
use std::str::FromStr;

#[derive(Debug, PartialEq)]
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
        let query_length = fields[1].parse::<usize>()?;
        let query_start = fields[2].parse::<usize>()?;
        let query_end = fields[3].parse::<usize>()?;
        let target_name = fields[5].to_string();
        let target_length = fields[6].parse::<usize>()?;
        let target_start = fields[7].parse::<usize>()?;
        let target_end = fields[8].parse::<usize>()?;
        let strand_char = fields[4].chars().next().ok_or(ParseErr::InvalidField)?;
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
    IoError(std::io::Error),
    InvalidField(std::num::ParseIntError),
    InvalidStrand,
}

impl From<std::num::ParseIntError> for ParseErr {
    fn from(err: std::num::ParseIntError) -> Self {
        ParseErr::InvalidField(err)
    }
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
                cigar: Some("60M".to_string()),
            }
        );
    }

    #[test]
    fn test_parse_paf_invalid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10M";
        assert!(PafRecord::parse(line).is_err());
    }
}
