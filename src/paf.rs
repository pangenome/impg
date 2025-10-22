//! PAF (Pairwise Alignment Format) parsing
//!
//! This module provides functions for parsing PAF format alignment files.
//! Supports both uncompressed and BGZF-compressed files with optional GZI indices.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::seqidx::SequenceIndex;
use log::debug;
use noodles::bgzf;
use std::cell::RefCell;
use std::collections::{hash_map::Entry, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, Read, Seek, SeekFrom};
use std::num::{NonZeroUsize, ParseIntError};

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

enum PafHandle {
    Plain(File),
    Compressed(bgzf::io::Reader<File>),
}

thread_local! {
    static PAF_FILE_CACHE: RefCell<HashMap<String, PafHandle>> = RefCell::new(HashMap::new());
}

pub fn read_cigar_data(
    alignment_file: &str,
    offset: u64,
    buffer: &mut [u8],
) -> Result<(), String> {
    let is_compressed = [".gz", ".bgz"]
        .iter()
        .any(|extension| alignment_file.ends_with(extension));

    with_paf_file_handle(alignment_file, is_compressed, |handle| match handle {
        PafHandle::Compressed(reader) => {
            let virtual_position = bgzf::VirtualPosition::from(offset);
            reader.seek(virtual_position).map_err(|e| {
                format!(
                    "Failed to seek in compressed file '{}': {}",
                    alignment_file, e
                )
            })?;
            reader.read_exact(buffer).map_err(|e| {
                format!(
                    "Failed to read data from compressed file '{}': {}",
                    alignment_file, e
                )
            })
        }
        PafHandle::Plain(file) => {
            file.seek(SeekFrom::Start(offset)).map_err(|e| {
                format!("Failed to seek in file '{}': {}", alignment_file, e)
            })?;
            file.read_exact(buffer).map_err(|e| {
                format!("Failed to read data from file '{}': {}", alignment_file, e)
            })
        }
    })
}

fn with_paf_file_handle<R, F>(
    alignment_file: &str,
    is_compressed: bool,
    f: F,
) -> Result<R, String>
where
    F: FnOnce(&mut PafHandle) -> Result<R, String>,
{
    PAF_FILE_CACHE.with(|cache_cell| -> Result<R, String> {
        let mut cache = cache_cell.borrow_mut();
        let entry = cache.entry(alignment_file.to_string());
        let handle = match entry {
            Entry::Occupied(entry) => entry.into_mut(),
            Entry::Vacant(entry) => {
                let handle = if is_compressed {
                    let file = File::open(alignment_file)
                        .map_err(|e| format!("Failed to open compressed file '{}': {}", alignment_file, e))?;
                    PafHandle::Compressed(bgzf::io::Reader::new(file))
                } else {
                    let file = File::open(alignment_file)
                        .map_err(|e| format!("Failed to open file '{}': {}", alignment_file, e))?;
                    PafHandle::Plain(file)
                };
                entry.insert(handle)
            }
        };
        f(handle)
    })
}

/// Parse a single PAF line into an AlignmentRecord
/// This is PAF-format specific parsing logic
fn parse_paf_line(
    line: &str,
    file_pos: u64,
    seq_index: &mut SequenceIndex,
) -> Result<AlignmentRecord, ParseErr> {
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
    let mut record = AlignmentRecord {
        query_id,
        query_start,
        query_end,
        target_id,
        target_start,
        target_end,
        strand_and_data_offset: cigar_offset,
        data_bytes: cigar_bytes,
    };
    record.set_strand(strand);

    Ok(record)
}

pub fn parse_paf<R: BufRead>(
    reader: R,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    let mut bytes_read: u64 = 0;
    let mut records = Vec::new();
    for line_result in reader.lines() {
        let line = line_result.map_err(ParseErr::IoError)?;
        let record = parse_paf_line(&line, bytes_read, seq_index)?;
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
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    use std::io::BufRead;

    let mut records = Vec::new();
    let mut line_buf = String::new();

    loop {
        // Get virtual position BEFORE reading the line
        let virtual_pos = reader.virtual_position();
        line_buf.clear();

        let bytes_read = reader.read_line(&mut line_buf).map_err(ParseErr::IoError)?;
        if bytes_read == 0 {
            break; // EOF
        }

        // Remove trailing newline
        let line = line_buf.trim_end();
        if line.is_empty() {
            continue;
        }

        // Parse the record using the virtual position
        let record = parse_paf_line(line, virtual_pos.into(), seq_index)?;
        records.push(record);
    }

    Ok(records)
}

/// Parse PAF from a BGZF-compressed file using a GZI index for faster multithreaded decompression.
/// After parsing with uncompressed offsets, converts them to virtual positions for seeking.
pub fn parse_paf_bgzf_with_gzi<R: std::io::Read>(
    reader: R,
    gzi_index: noodles::bgzf::gzi::Index,
    seq_index: &mut SequenceIndex,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    // First pass: parse with uncompressed byte offsets
    let reader = std::io::BufReader::new(reader);
    let mut records = parse_paf(reader, seq_index)?;

    // Second pass: convert uncompressed offsets to virtual positions using GZI
    for record in &mut records {
        // Extract the uncompressed offset (ignoring the strand bit)
        let uncompressed_offset = record.strand_and_data_offset & !AlignmentRecord::STRAND_BIT;

        // Convert to virtual position using GZI query
        let virtual_pos = gzi_index.query(uncompressed_offset).map_err(|e| {
            ParseErr::InvalidFormat(format!(
                "Failed to find virtual position for offset {}: {:?}",
                uncompressed_offset, e
            ))
        })?;

        // Update the record with virtual position, preserving strand bit
        let strand_bit = record.strand_and_data_offset & AlignmentRecord::STRAND_BIT;
        record.strand_and_data_offset = u64::from(virtual_pos) | strand_bit;
    }

    Ok(records)
}

/// Parse a PAF file with automatic format detection (compressed or uncompressed)
/// Handles .gz/.bgz compression with optional GZI index for multithreaded decompression
pub fn parse_paf_file(
    paf_file: &str,
    file: File,
    threads: NonZeroUsize,
    seq_index: &mut SequenceIndex,
) -> std::io::Result<Vec<AlignmentRecord>> {
    if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        let gzi_path = format!("{}.gzi", paf_file);
        if std::path::Path::new(&gzi_path).exists() {
            debug!(
                "Found GZI index for {}, using multithreaded decompression",
                paf_file
            );
            let gzi_index = noodles::bgzf::gzi::fs::read(&gzi_path).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to read GZI index {}: {}", gzi_path, e),
                )
            })?;
            let mt_reader =
                noodles::bgzf::io::MultithreadedReader::with_worker_count(threads, file);
            parse_paf_bgzf_with_gzi(mt_reader, gzi_index, seq_index).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse PAF from {}: {:?}", paf_file, e),
                )
            })
        } else {
            debug!("No GZI index for {}, using BGZF reader", paf_file);
            let bgzf_reader = noodles::bgzf::io::Reader::new(file);
            parse_paf_bgzf(bgzf_reader, seq_index).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to parse PAF from {}: {:?}", paf_file, e),
                )
            })
        }
    } else {
        let reader = BufReader::new(file);
        parse_paf(reader, seq_index).map_err(|e| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("Failed to parse PAF from {}: {:?}", paf_file, e),
            )
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_parse_paf_valid() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255";
        let mut seq_index = SequenceIndex::new();
        let record = parse_paf_line(line, 0, &mut seq_index).unwrap();

        // IDs should be 0 and 1 as they're the first entries in the SequenceIndex
        let query_id = seq_index.get_id("seq1").unwrap();
        let target_id = seq_index.get_id("seq2").unwrap();

        assert_eq!(
            record,
            AlignmentRecord {
                query_id,
                query_start: 0,
                query_end: 100,
                target_id,
                target_start: 0,
                target_end: 100,
                // If no cigar, offset is line length; data_bytes=0
                strand_and_data_offset: (line.len() + 1) as u64,
                data_bytes: 0,
            }
        );
    }

    #[test]
    fn test_parse_paf_valid_2() {
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10=";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_ok());
    }

    #[test]
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10Q";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index).is_err());
    }
}
