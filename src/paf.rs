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
    pub strand_and_format_and_offset: u64, // Track strand, format (cigar/tp), and offset
    pub cg_or_tp_bytes: usize,
}

#[derive(Default, PartialEq, Clone, Copy)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum AlignmentFormat {
    Cigar,
    Tracepoints,
}

impl PartialPafRecord {
    // Constants for bit manipulation
    const STRAND_BIT: u64 = 0x8000000000000000;     // Most significant bit for u64
    const FORMAT_BIT: u64 = 0x4000000000000000;     // Second most significant bit for format (0=CIGAR, 1=Tracepoints)
    const OFFSET_MASK: u64 = 0x3FFFFFFFFFFFFFFF;    // Remaining 62 bits for offset

    pub fn strand(&self) -> Strand {
        if (self.strand_and_format_and_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }
    pub fn set_strand(&mut self, strand: Strand) {
        match strand {
            Strand::Forward => self.strand_and_format_and_offset &= !Self::STRAND_BIT,
            Strand::Reverse => self.strand_and_format_and_offset |= Self::STRAND_BIT,
        }
    }

    pub fn is_tracepoints(&self) -> bool {
        (self.strand_and_format_and_offset & Self::FORMAT_BIT) != 0
    }

    pub fn set_format_tracepoints(&mut self, is_tracepoints: bool) {
        if is_tracepoints {
            self.strand_and_format_and_offset |= Self::FORMAT_BIT;
        } else {
            self.strand_and_format_and_offset &= !Self::FORMAT_BIT;
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

        let mut cg_or_tp_offset: u64 = file_pos;
        let mut cg_or_tp_bytes: usize = 0;
        let mut is_tracepoints = false;

        // Look for both cg:Z: and tp:Z: tags (the former has precedence)
        for tag_str in fields.iter() {
            if tag_str.starts_with("cg:Z:") {
                cg_or_tp_offset += 5;
                cg_or_tp_bytes = tag_str.len() - 5;
                is_tracepoints = false;
                break;
            } else if tag_str.starts_with("tp:Z:") {
                cg_or_tp_offset += 5;
                cg_or_tp_bytes = tag_str.len() - 5;
                is_tracepoints = true;
                break;
            } else {
                cg_or_tp_offset += (tag_str.len() + 1) as u64;
            }
        }

        // Create the record
        let mut record = Self {
            query_id,
            query_start,
            query_end,
            target_id,
            target_start,
            target_end,
            strand_and_format_and_offset: cg_or_tp_offset & Self::OFFSET_MASK,
            cg_or_tp_bytes,
        };
        
        // Set strand and format
        record.set_strand(strand);
        record.set_format_tracepoints(is_tracepoints);

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
                strand_and_format_and_offset: (line.len() + 1) as u64,
                cg_or_tp_bytes: 0,
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
