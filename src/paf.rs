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

/// Parse PAF from a BGZF-compressed file, storing virtual positions for seeking.
/// If a GZI index is provided, uses it for faster multithreaded decompression and converts
/// uncompressed offsets to virtual positions. Otherwise, reads with single-threaded BGZF reader.
pub fn parse_paf_bgzf<R: std::io::Read + std::io::Seek>(
    mut reader: noodles::bgzf::io::Reader<R>,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<PartialPafRecord>, ParseErr> {
    use std::io::{BufRead, Read};

    let mut records = Vec::new();
    let mut line_bytes = Vec::new();

    loop {
        // Get virtual position BEFORE reading
        let line_start_vpos = reader.virtual_position();
        line_bytes.clear();

        let bytes_read = reader
            .read_until(b'\n', &mut line_bytes)
            .map_err(ParseErr::IoError)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Capture position after reading the line (for next iteration)
        let line_end_vpos = reader.virtual_position();

        // Convert to string for parsing (excluding newline)
        let line_len = if line_bytes.ends_with(&[b'\n']) {
            line_bytes.len() - 1
        } else {
            line_bytes.len()
        };
        let line = std::str::from_utf8(&line_bytes[..line_len])
            .map_err(|_| ParseErr::InvalidFormat("Invalid UTF-8".to_string()))?;

        if line.is_empty() {
            continue;
        }

        // Parse to get byte offset to CIGAR within the line
        let mut record = PartialPafRecord::parse(line, 0, seq_index)?;
        let cigar_byte_offset = record.strand_and_cigar_offset & !PartialPafRecord::STRAND_BIT;

        // Compute CIGAR virtual position by seeking back and advancing by the offset
        // This correctly handles BGZF block boundaries
        if cigar_byte_offset > 0 {
            reader.seek(line_start_vpos).map_err(ParseErr::IoError)?;
            std::io::copy(
                &mut reader.by_ref().take(cigar_byte_offset),
                &mut std::io::sink(),
            )
            .map_err(ParseErr::IoError)?;
        } else {
            reader.seek(line_start_vpos).map_err(ParseErr::IoError)?;
        }

        let cigar_vpos = reader.virtual_position();

        // Update record with correct BGZF virtual position
        let strand_bit = record.strand_and_cigar_offset & PartialPafRecord::STRAND_BIT;
        record.strand_and_cigar_offset = u64::from(cigar_vpos) | strand_bit;

        records.push(record);

        // Restore reader position to end of line for next iteration
        reader.seek(line_end_vpos).map_err(ParseErr::IoError)?;
    }

    Ok(records)
}

/// Parse PAF from a BGZF-compressed file using a GZI index for faster multithreaded decompression.
/// After parsing with uncompressed offsets, converts them to virtual positions for seeking.
pub fn parse_paf_bgzf_with_gzi<R: std::io::Read>(
    reader: R,
    gzi_index: noodles::bgzf::gzi::Index,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<PartialPafRecord>, ParseErr> {
    // First pass: parse with uncompressed byte offsets
    let reader = std::io::BufReader::new(reader);
    let mut records = parse_paf(reader, seq_index)?;

    // Second pass: convert uncompressed offsets to virtual positions using GZI
    for record in &mut records {
        // Extract the uncompressed offset (ignoring the strand bit)
        let uncompressed_offset = record.strand_and_cigar_offset & !PartialPafRecord::STRAND_BIT;

        // Convert to virtual position using GZI query
        let virtual_pos = gzi_index.query(uncompressed_offset).map_err(|e| {
            ParseErr::InvalidFormat(format!(
                "Failed to find virtual position for offset {}: {:?}",
                uncompressed_offset, e
            ))
        })?;

        // Update the record with virtual position, preserving strand bit
        let strand_bit = record.strand_and_cigar_offset & PartialPafRecord::STRAND_BIT;
        record.strand_and_cigar_offset = u64::from(virtual_pos) | strand_bit;
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
