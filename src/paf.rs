//! PAF (Pairwise Alignment Format) parsing
//!
//! This module provides functions for parsing PAF format alignment files.
//! Supports both uncompressed and BGZF-compressed files with optional GZI indices.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::seqidx::SequenceIndex;
use log::{debug, warn};
use rayon::prelude::*;
use noodles::bgzf;
use std::fs::File;
use std::io::{BufRead, BufReader, Error as IoError, Read, Seek, SeekFrom};
use std::num::{NonZeroUsize, ParseIntError};

#[derive(Debug)]
pub enum ParseErr {
    NotEnoughFields,
    IoError(IoError),
    InvalidField(ParseIntError),
    InvalidStrand,
    UnsupportedCigarOperation,
    InvalidFormat(String),
}

impl std::fmt::Display for ParseErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseErr::NotEnoughFields => write!(f, "Not enough fields in PAF record"),
            ParseErr::IoError(e) => write!(f, "IO error: {}", e),
            ParseErr::InvalidField(e) => write!(f, "Invalid field: {}", e),
            ParseErr::InvalidStrand => write!(f, "Invalid strand"),
            ParseErr::UnsupportedCigarOperation => write!(f, "Unsupported CIGAR operation"),
            ParseErr::InvalidFormat(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for ParseErr {}

enum PafHandle {
    Plain(File),
    Compressed(bgzf::io::Reader<File>),
}

const BGZF_HEADER_SIZE: usize = 18;

/// Check whether a file starts with a valid BGZF header.
/// Returns `Ok(false)` for regular gzip, too-small files, or plain text.
fn is_bgzf<R: Read + Seek>(reader: &mut R) -> std::io::Result<bool> {
    let mut header = [0u8; BGZF_HEADER_SIZE];
    let result = match reader.read_exact(&mut header) {
        Ok(()) => {
            Ok(header[0..2] == [0x1f, 0x8b]      // gzip magic
                && header[2] == 0x08              // DEFLATE
                && header[3] == 0x04              // FEXTRA
                && header[10..12] == [0x06, 0x00] // XLEN=6
                && header[12..14] == [b'B', b'C'] // BC subfield
                && header[14..16] == [0x02, 0x00]) // SLEN=2
        }
        Err(e) if e.kind() == std::io::ErrorKind::UnexpectedEof => Ok(false),
        Err(e) => Err(e),
    };
    reader.seek(SeekFrom::Start(0))?;
    result
}

pub fn read_cigar_data(alignment_file: &str, offset: u64, buffer: &mut [u8]) -> Result<(), String> {
    let is_compressed = [".gz", ".bgz"]
        .iter()
        .any(|extension| alignment_file.ends_with(extension));

    let handle = if is_compressed {
        let mut file = File::open(alignment_file)
            .map_err(|e| format!("Failed to open compressed file '{}': {}", alignment_file, e))?;
        if !is_bgzf(&mut file)
            .map_err(|e| format!("Failed to read header of '{}': {}", alignment_file, e))?
        {
            return Err(format!(
                "'{}' is regular gzip, not BGZF. Convert with: zcat '{}' | bgzip > output.paf.gz",
                alignment_file, alignment_file
            ));
        }
        PafHandle::Compressed(bgzf::io::Reader::new(file))
    } else {
        let file = File::open(alignment_file)
            .map_err(|e| format!("Failed to open file '{}': {}", alignment_file, e))?;
        PafHandle::Plain(file)
    };

    match handle {
        PafHandle::Compressed(mut reader) => {
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
        PafHandle::Plain(mut file) => {
            file.seek(SeekFrom::Start(offset))
                .map_err(|e| format!("Failed to seek in file '{}': {}", alignment_file, e))?;
            file.read_exact(buffer)
                .map_err(|e| format!("Failed to read data from file '{}': {}", alignment_file, e))
        }
    }
}

/// Validate the payload of a `cg:Z` tag and return the query and target bases
/// it spans and whether it uses `M`, in one pass over its bytes.
///
/// impg re-reads the CIGAR out of the file on every query. Without this check a
/// malformed CIGAR was only discovered inside a rayon worker, where
/// `CigarOp::new` panics on an unknown operation and the run dies without
/// naming the offending file or record. The spans let the parser refuse a
/// CIGAR that does not cover the record coordinates; a separate pass for them
/// cost 12.6% of the instructions of `impg index` on 400,000 records.
fn cigar_tag_spans(cigar: &str) -> Result<(u64, u64, bool), String> {
    // `CigarOp` packs the length into 29 bits
    const MAX_OP_LEN: u64 = 1 << 29;
    let (mut q, mut t, mut len, mut has_len, mut has_m) = (0u64, 0u64, 0u64, false, false);
    let mut ok = !cigar.is_empty();
    for &c in cigar.as_bytes() {
        let d = c.wrapping_sub(b'0');
        if d < 10 {
            len = len.saturating_mul(10).saturating_add(u64::from(d));
            has_len = true;
            continue;
        }
        // A length out of range would overflow the span sums
        if !has_len || len >= MAX_OP_LEN {
            return Err(describe_cigar_error(cigar.as_bytes()));
        }
        match c {
            b'=' | b'X' => (q, t) = (q + len, t + len),
            b'M' => (q, t, has_m) = (q + len, t + len, true),
            b'I' => q += len,
            b'D' => t += len,
            _ => ok = false,
        }
        (len, has_len) = (0, false);
    }
    if ok & !has_len {
        Ok((q, t, has_m))
    } else {
        Err(describe_cigar_error(cigar.as_bytes()))
    }
}

/// Reason a CIGAR fails, or an empty string when it is valid.
fn describe_cigar_error(bytes: &[u8]) -> String {
    // `CigarOp` packs the length into 29 bits
    const MAX_OP_LEN: u64 = 1 << 29;
    if bytes.is_empty() {
        return "it is empty".to_string();
    }
    let mut run_start = 0;
    for (i, &byte) in bytes.iter().enumerate() {
        if byte.is_ascii_digit() {
            continue;
        }
        let digits = &bytes[run_start..i];
        if digits.is_empty() {
            return format!("operation '{}' has no length", byte as char);
        }
        if !matches!(byte, b'=' | b'X' | b'I' | b'D' | b'M') {
            return format!(
                "operation '{}' is not one of =XIDM (S/H clips, N and P are unsupported)",
                byte as char
            );
        }
        let len = digits.iter().fold(0u64, |n, d| n.saturating_mul(10).saturating_add(u64::from(d - b'0')));
        if len >= MAX_OP_LEN {
            return format!("an operation length reaches {MAX_OP_LEN}, the packed maximum");
        }
        run_start = i + 1;
    }
    if run_start < bytes.len() {
        return "it ends with a length and no operation".to_string();
    }
    String::new()
}

/// Parse a single PAF line into an AlignmentRecord
/// This is PAF-format specific parsing logic. `m_records` counts the records
/// whose CIGAR uses `M`.
fn parse_paf_line(
    line: &str,
    file_pos: u64,
    seq_index: &mut SequenceIndex,
    m_records: &mut u64,
) -> Result<AlignmentRecord, ParseErr> {
    let mut fields = line.split('\t');
    let query_name = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let query_len_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let query_start_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let query_end_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let strand_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let target_name = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let target_len_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let target_start_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let target_end_str = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let f9 = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let f10 = fields.next().ok_or(ParseErr::NotEnoughFields)?;
    let f11 = fields.next().ok_or(ParseErr::NotEnoughFields)?;

    let query_length = query_len_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let query_start = query_start_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let query_end = query_end_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_length = target_len_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_start = target_start_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let target_end = target_end_str.parse::<usize>().map_err(ParseErr::InvalidField)?;
    let strand_char = strand_str
        .chars()
        .next()
        .ok_or_else(|| ParseErr::InvalidFormat("Expected '+' or '-' for strand".to_string()))?;
    let strand = match strand_char {
        '+' => Strand::Forward,
        '-' => Strand::Reverse,
        _ => return Err(ParseErr::InvalidStrand),
    };

    // An end past the sequence length gave a PAF with qend > qlen downstream
    if query_start > query_end
        || query_end > query_length
        || target_start > target_end
        || target_end > target_length
    {
        return Err(ParseErr::InvalidFormat(format!(
            "record {query_name}:{query_start}-{query_end} (length {query_length}) -> \
             {target_name}:{target_start}-{target_end} (length {target_length}): \
             every start must be <= its end, and every end <= its sequence length"
        )));
    }

    // Convert names to IDs using the SequenceIndex
    let query_id = seq_index
        .get_or_insert_id_checked(query_name, query_length)
        .map_err(ParseErr::InvalidFormat)?;
    let target_id = seq_index
        .get_or_insert_id_checked(target_name, target_length)
        .map_err(ParseErr::InvalidFormat)?;

    // Byte offset up to the start of optional tags (column 12)
    let prefix_len: u64 = [
        query_name,
        query_len_str,
        query_start_str,
        query_end_str,
        strand_str,
        target_name,
        target_len_str,
        target_start_str,
        target_end_str,
        f9,
        f10,
        f11,
    ]
    .iter()
    .map(|s| (s.len() + 1) as u64)
    .sum();
    let mut cigar_offset: u64 = file_pos + prefix_len;
    let mut cigar_bytes: usize = 0;

    for tag_str in fields {
        if let Some(cigar) = tag_str.strip_prefix("cg:Z:") {
            cigar_offset += 5;
            let (q_span, t_span, has_m) = cigar_tag_spans(cigar).map_err(|reason| {
                ParseErr::InvalidFormat(format!(
                    "invalid cg:Z CIGAR in record {query_name}:{query_start}-{query_end} -> \
                     {target_name}:{target_start}-{target_end}: {reason}"
                ))
            })?;
            cigar_bytes = cigar.len();
            // A CIGAR shorter than the coordinates lost the projection of the last bases
            let (q_len, t_len) = (query_end - query_start, target_end - target_start);
            if q_span != q_len as u64 || t_span != t_len as u64 {
                return Err(ParseErr::InvalidFormat(format!(
                    "record {query_name}:{query_start}-{query_end} -> \
                     {target_name}:{target_start}-{target_end}: the CIGAR spans {q_span} query \
                     and {t_span} target bases, the coordinates {q_len} and {t_len}"
                )));
            }
            *m_records += u64::from(has_m);
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
    read_paf(reader, seq_index, &mut 0, rayon::current_num_threads())
}

/// Bytes of PAF lines one parallel task parses.
const PAF_BLOCK_BYTES: usize = 1 << 20;

/// Parse PAF with `threads` blocks of 1 MiB per round. The caller gives each of
/// the files it parses at once its share of the threads, so the text held is
/// about 1 MiB per thread in total: 16 MiB per file grew to 336 MB against
/// 155 MB for 16 files at -t 16.
fn read_paf<R: Read>(
    reader: R,
    seq_index: &mut SequenceIndex,
    m_records: &mut u64,
    threads: usize,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    read_paf_blocks(reader, seq_index, m_records, PAF_BLOCK_BYTES, threads.max(1))
}

/// Parse PAF lines in blocks of about `block_bytes`, one rayon task per block.
///
/// The CIGAR span check costs 2 ns per CIGAR byte against 0.4 ns for the
/// syntax check, so a parse on one thread took +57% index time on a real
/// minimap2 PAF. Each block interns its names locally and the blocks merge in
/// file order: ids, offsets and the first error are those of a sequential read.
fn read_paf_blocks<R: Read>(
    mut reader: R,
    seq_index: &mut SequenceIndex,
    m_records: &mut u64,
    block_bytes: usize,
    blocks_per_round: usize,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    let batch = block_bytes * blocks_per_round;
    let mut records = Vec::new();
    let mut buf: Vec<u8> = Vec::new();
    let mut base = 0u64; // file offset of buf[0]
    loop {
        // Read a batch, then on until a line ends: a line can be longer than a batch
        let eof = loop {
            let old_len = buf.len();
            let want = batch.saturating_sub(old_len).max(1 << 20);
            let n = (&mut reader)
                .take(want as u64)
                .read_to_end(&mut buf)
                .map_err(ParseErr::IoError)?;
            if n == 0 {
                break true;
            }
            if buf.len() >= batch && buf[old_len..].contains(&b'\n') {
                break false;
            }
        };
        let cut = if eof {
            buf.len()
        } else {
            buf.iter().rposition(|&b| b == b'\n').map_or(0, |i| i + 1)
        };
        // Blocks of whole lines
        let mut blocks = Vec::new();
        let mut start = 0;
        while start < cut {
            let end = match buf[(start + block_bytes).min(cut)..cut].iter().position(|&b| b == b'\n') {
                Some(i) if start + block_bytes < cut => start + block_bytes + i + 1,
                _ => cut,
            };
            blocks.push((start, end));
            start = end;
        }
        let parsed: Vec<_> = blocks
            .par_iter()
            .map(|&(start, end)| -> Result<(Vec<AlignmentRecord>, SequenceIndex, u64), ParseErr> {
                let (mut local, mut m, mut out) = (SequenceIndex::new(), 0u64, Vec::new());
                let mut pos = base + start as u64;
                for line in buf[start..end].split_inclusive(|&b| b == b'\n') {
                    let text = std::str::from_utf8(line).map_err(|_| {
                        ParseErr::InvalidFormat(format!("the line at byte {pos} is not valid UTF-8"))
                    })?;
                    // Strip the terminator, but advance by the raw byte count: the
                    // CIGAR offsets are absolute file positions, also with CRLF
                    let trimmed = text.trim_end_matches(['\n', '\r']);
                    if trimmed.is_empty() {
                        return Err(ParseErr::InvalidFormat(format!("the line at byte {pos} is blank")));
                    }
                    out.push(parse_paf_line(trimmed, pos, &mut local, &mut m)?);
                    pos += line.len() as u64;
                }
                Ok((out, local, m))
            })
            .collect();
        for block in parsed {
            let (mut block_records, local, m) = block?;
            // Local ids are dense and in order of first appearance
            let to_global = (0..local.len() as u32)
                .map(|id| {
                    let name = local.get_name(id).expect("dense local sequence ids");
                    let len = local.get_len_from_id(id).expect("every PAF name has a length");
                    seq_index.get_or_insert_id_checked(name, len).map_err(ParseErr::InvalidFormat)
                })
                .collect::<Result<Vec<u32>, _>>()?;
            for record in &mut block_records {
                record.query_id = to_global[record.query_id as usize];
                record.target_id = to_global[record.target_id as usize];
            }
            records.extend(block_records);
            *m_records += m;
        }
        base += cut as u64;
        buf.drain(..cut);
        if eof {
            break;
        }
    }
    Ok(records)
}

/// Build the GZI index of a BGZF file from its block headers and footers,
/// without inflating any block. Without a `.gzi`, the BGZF reader seeked back
/// to the start of every line to find its CIGAR: 3.66 s against 0.11 s with a
/// `.gzi` on 100,000 records.
fn scan_gzi<R: Read + Seek>(file: &mut R, paf_file: &str) -> std::io::Result<bgzf::gzi::Index> {
    let mut reader = BufReader::with_capacity(1 << 20, &mut *file);
    let (mut cpos, mut upos) = (0u64, 0u64);
    let mut entries = Vec::new();
    let (mut header, mut footer, mut extra) = ([0u8; 12], [0u8; 8], Vec::new());
    while !reader.fill_buf()?.is_empty() {
        let invalid = |what: &str| {
            std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!("{paf_file}: BGZF block at byte {cpos}: {what}"),
            )
        };
        reader.read_exact(&mut header).map_err(|_| invalid("truncated header"))?;
        if header[..4] != [0x1f, 0x8b, 0x08, 0x04] {
            return Err(invalid("not a BGZF block header"));
        }
        let xlen = usize::from(u16::from_le_bytes([header[10], header[11]]));
        extra.resize(xlen, 0);
        reader.read_exact(&mut extra).map_err(|_| invalid("truncated header"))?;
        // The BC subfield holds the block size minus 1
        let (mut i, mut bsize) = (0, None);
        while i + 4 <= xlen {
            let slen = usize::from(u16::from_le_bytes([extra[i + 2], extra[i + 3]]));
            if extra[i..i + 2] == *b"BC" && slen == 2 && i + 6 <= xlen {
                bsize = Some(u64::from(u16::from_le_bytes([extra[i + 4], extra[i + 5]])) + 1);
            }
            i += 4 + slen;
        }
        let block_size = bsize.ok_or_else(|| invalid("no BC subfield"))?;
        let cdata = block_size
            .checked_sub(12 + xlen as u64 + 8)
            .ok_or_else(|| invalid("block size smaller than its header"))?;
        reader.seek_relative(cdata as i64)?;
        reader.read_exact(&mut footer).map_err(|_| invalid("truncated block"))?;
        if cpos > 0 {
            entries.push((cpos, upos));
        }
        cpos += block_size;
        upos += u64::from(u32::from_le_bytes([footer[4], footer[5], footer[6], footer[7]]));
    }
    drop(reader);
    file.seek(SeekFrom::Start(0))?;
    debug!("{paf_file}: scanned {} BGZF blocks", entries.len() + 1);
    Ok(bgzf::gzi::Index::from(entries))
}

/// Parse PAF from a BGZF-compressed file using a GZI index for faster multithreaded decompression.
/// After parsing with uncompressed offsets, converts them to virtual positions for seeking.
pub fn parse_paf_bgzf_with_gzi<R: std::io::Read>(
    reader: R,
    gzi_index: noodles::bgzf::gzi::Index,
    seq_index: &mut SequenceIndex,
    m_records: &mut u64,
    threads: usize,
) -> Result<Vec<AlignmentRecord>, ParseErr> {
    // First pass: parse with uncompressed byte offsets
    let reader = std::io::BufReader::new(reader);
    let mut records = read_paf(reader, seq_index, m_records, threads)?;

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
    let mut m_records = 0u64;
    let records = if [".gz", ".bgz"].iter().any(|e| paf_file.ends_with(e)) {
        let mut file = file;
        if !is_bgzf(&mut file)? {
            return Err(std::io::Error::new(
                std::io::ErrorKind::InvalidData,
                format!(
                    "'{}' is regular gzip, not BGZF. Convert with: zcat '{}' | bgzip > output.paf.gz",
                    paf_file, paf_file
                ),
            ));
        }
        let gzi_path = format!("{}.gzi", paf_file);
        let gzi_index = if std::path::Path::new(&gzi_path).exists() {
            debug!("Found GZI index for {}", paf_file);
            noodles::bgzf::gzi::fs::read(&gzi_path).map_err(|e| {
                std::io::Error::new(
                    std::io::ErrorKind::InvalidData,
                    format!("Failed to read GZI index {}: {}", gzi_path, e),
                )
            })?
        } else {
            debug!("No GZI index for {}, scanning its BGZF blocks", paf_file);
            scan_gzi(&mut file, paf_file)?
        };
        let mt_reader = noodles::bgzf::io::MultithreadedReader::with_worker_count(threads, file);
        parse_paf_bgzf_with_gzi(mt_reader, gzi_index, seq_index, &mut m_records, threads.get())
    } else {
        read_paf(BufReader::new(file), seq_index, &mut m_records, threads.get())
    }
    .map_err(|e| {
        std::io::Error::new(
            std::io::ErrorKind::InvalidData,
            format!("Failed to parse PAF from {}: {}", paf_file, e),
        )
    })?;
    if m_records > 0 {
        warn!(
            "{paf_file}: {m_records} of {} records use 'M' CIGAR operations. impg counts every M base \
             as a match, so their identity (gi:f, --min-result-identity) is overestimated. \
             Use =/X CIGARs (wfmash, or minimap2 --eqx).",
            records.len()
        );
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
        let record = parse_paf_line(line, 0, &mut seq_index, &mut 0).unwrap();

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
    fn test_parse_paf_invalid() {
        // it's got a character 'z' in the length field
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\tz\t100\t60\t100\t255\tcg:Z:10M";
        let mut seq_index = SequenceIndex::new();
        assert!(parse_paf_line(line, 0, &mut seq_index, &mut 0).is_err());
    }

    #[test]
    fn test_parse_paf_cigar_invalid() {
        // it's got Q in the CIGAR string (all other fields are valid)
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10Q";
        let mut seq_index = SequenceIndex::new();
        let err = parse_paf_line(line, 0, &mut seq_index, &mut 0).unwrap_err().to_string();
        assert!(err.contains("'Q'") && err.contains("seq1:0-100 -> seq2:0-100"), "{err}");
    }

    /// Write `paf` as BGZF and parse it as `impg index` does, without a `.gzi`
    fn parse_bgzf(paf: &str) -> (tempfile::TempDir, String, Vec<AlignmentRecord>) {
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("x.paf.gz").to_string_lossy().into_owned();
        let mut writer = noodles::bgzf::io::Writer::new(File::create(&path).unwrap());
        std::io::Write::write_all(&mut writer, paf.as_bytes()).unwrap();
        writer.finish().unwrap();
        let one = NonZeroUsize::new(1).unwrap();
        let records =
            parse_paf_file(&path, File::open(&path).unwrap(), one, &mut SequenceIndex::new()).unwrap();
        (dir, path, records)
    }

    #[test]
    fn test_parse_paf_bgzf_crlf() {
        // The BGZF reader kept the `\r`, so the CIGAR ended in `\r` and was rejected
        let (_dir, _, records) = parse_bgzf("seq1\t100\t0\t50\t+\tseq2\t100\t0\t50\t50\t50\t60\tcg:Z:50=\r\n");
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].data_bytes, 3, "the CIGAR is `50=` without the `\\r`");
        // A blank line is an error here too (the reader without .gzi skipped it)
        let dir = tempfile::TempDir::new().unwrap();
        let path = dir.path().join("blank.paf.gz").to_string_lossy().into_owned();
        let mut writer = noodles::bgzf::io::Writer::new(File::create(&path).unwrap());
        std::io::Write::write_all(&mut writer, b"seq1\t100\t0\t50\t+\tseq2\t100\t0\t50\t50\t50\t60\tcg:Z:50=\n\n").unwrap();
        writer.finish().unwrap();
        let one = NonZeroUsize::new(1).unwrap();
        let err = parse_paf_file(&path, File::open(&path).unwrap(), one, &mut SequenceIndex::new()).unwrap_err();
        assert!(err.to_string().contains("is blank"), "{err}");
    }

    #[test]
    fn test_parse_paf_bgzf_without_gzi_finds_every_cigar() {
        // 3000 lines over several 64 KiB blocks, a last line without newline
        // and without cg:Z (the old reader underflowed on it)
        let mut paf = String::new();
        let mut want = Vec::new();
        for i in 0..3000usize {
            let (a, b) = (10 + i % 50, 39 - i % 30);
            let span = a + 1 + b;
            let cg = format!("{a}={}X{b}=", 1);
            paf += &format!("q\t100000\t{i}\t{}\t+\tt\t100000\t{i}\t{}\t{span}\t{span}\t60\tcg:Z:{cg}\n", i + span, i + span);
            want.push(cg);
        }
        paf += "q\t100000\t0\t50\t+\tt\t100000\t0\t50\t50\t50\t60\ttp:A:P";
        let (_dir, path, records) = parse_bgzf(&paf);
        assert_eq!(records.len(), 3001);
        assert_eq!(records[3000].data_bytes, 0);
        for (record, cg) in records.iter().zip(&want) {
            let mut buf = vec![0u8; record.data_bytes];
            read_cigar_data(&path, record.data_offset(), &mut buf).unwrap();
            assert_eq!(std::str::from_utf8(&buf).unwrap(), cg);
        }
    }

    #[test]
    fn test_parse_paf_cigar_with_clip_is_rejected() {
        // impg cannot project soft/hard clips (`CigarOp::new` panics); reject
        // them at index time instead of inside a query worker.
        let line = "seq1\t100\t0\t100\t+\tseq2\t100\t0\t100\t60\t100\t255\tcg:Z:10S90=";
        let mut seq_index = SequenceIndex::new();
        let err = parse_paf_line(line, 0, &mut seq_index, &mut 0).unwrap_err().to_string();
        assert!(err.contains("'S'") && err.contains("clips"), "{err}");
    }

    #[test]
    fn test_cigar_tag_spans_edges() {
        for (ok, spans) in [
            ("10=20X3I4D5M", (38, 39, true)),
            ("536870911=", (536870911, 536870911, false)),
            ("0000000001=", (1, 1, false)),
            ("1=", (1, 1, false)),
        ] {
            assert_eq!(cigar_tag_spans(ok), Ok(spans), "{ok}");
        }
        for bad in ["", "=", "100", "10=5", "536870912=", "99999999999X", "10=X", "5N", "3P", "18446744073709551615=1="] {
            assert!(cigar_tag_spans(bad).is_err(), "{bad:?} must be rejected");
        }
    }

    #[test]
    fn test_cigar_tag_spans_agrees_with_the_exact_check() {
        // The one-pass check and `describe_cigar_error` (its error messages)
        // must agree. 200,000 random strings over digits and ops, lengths
        // 0-14, many near the 29-bit limit
        let alphabet = b"0123456789=XIDMSN";
        let mut x: u64 = 0x2545F4914F6CDD1D;
        let mut next = || {
            x ^= x << 13;
            x ^= x >> 7;
            x ^= x << 17;
            x
        };
        let mut rejected = 0;
        for _ in 0..200_000 {
            let len = (next() % 15) as usize;
            let s: String = (0..len)
                .map(|_| alphabet[(next() % alphabet.len() as u64) as usize] as char)
                .collect();
            let exact_ok = describe_cigar_error(s.as_bytes()).is_empty();
            assert_eq!(cigar_tag_spans(&s).is_ok(), exact_ok, "{s:?}");
            rejected += usize::from(!exact_ok);
        }
        assert!(rejected > 0 && rejected < 200_000, "{rejected} rejected");
    }

    #[test]
    fn test_parse_paf_refuses_inconsistent_records() {
        let cases = [
            // 8= on a declared span of 10 lost the projection of the last 2 bp
            ("q\t100\t0\t10\t+\tt\t100\t0\t10\t8\t10\t60\tcg:Z:8=", "the CIGAR spans 8 query"),
            ("q\t100\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:10=2D", "and 12 target"),
            ("q\t100\t20\t10\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:10=", "start must be <= its end"),
            ("q\t100\t95\t105\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:10=", "end <= its sequence length"),
            ("q\t100\t0\t10\t+\tt\t5\t0\t10\t10\t10\t60\tcg:Z:10=", "end <= its sequence length"),
        ];
        for (line, want) in cases {
            let err = parse_paf_line(line, 0, &mut SequenceIndex::new(), &mut 0)
                .unwrap_err()
                .to_string();
            assert!(err.contains(want), "{line:?}: {err}");
        }
        // The same contig with two lengths gave a PAF with qend > qlen
        let paf = "q\t100\t0\t10\t+\tt\t100\t0\t10\t10\t10\t60\tcg:Z:10=\n\
q\t50\t0\t10\t+\tt\t100\t10\t20\t10\t10\t60\tcg:Z:10=\n";
        let err = parse_paf(BufReader::new(paf.as_bytes()), &mut SequenceIndex::new()).unwrap_err();
        assert!(err.to_string().contains("'q' has length 100 in one record and 50"), "{err}");
    }

    #[test]
    fn test_parse_paf_counts_m_records() {
        let paf = "q\t100\t0\t10\t+\tt\t100\t0\t10\t8\t10\t60\tcg:Z:10M\n\
q\t100\t10\t20\t+\tt\t100\t10\t20\t10\t10\t60\tcg:Z:10=\n\
q\t100\t20\t30\t+\tt\t100\t20\t31\t10\t11\t60\tcg:Z:4M1D6=\n";
        let mut m_records = 0;
        let records =
            read_paf(BufReader::new(paf.as_bytes()), &mut SequenceIndex::new(), &mut m_records, 1)
                .unwrap();
        assert_eq!((records.len(), m_records), (3, 2));
    }

    #[test]
    fn test_parse_paf_blocks_match_a_sequential_read() {
        // 300 CRLF and LF lines over 20 names, cut in blocks of 64 bytes
        let mut paf = String::new();
        for i in 0..300usize {
            let (q, t, eol) = (i % 7, i % 13 + 7, if i % 3 == 0 { "\r\n" } else { "\n" });
            paf += &format!("q{q}\t1000\t{i}\t{}\t+\tt{t}\t2000\t{i}\t{}\t20\t20\t60\tcg:Z:{}={}X{}={eol}", i + 20, i + 20, 1 + i % 9, 1, 18 - i % 9);
        }
        let (mut one, mut many) = (SequenceIndex::new(), SequenceIndex::new());
        let (mut m_one, mut m_many) = (0, 0);
        let whole = read_paf_blocks(paf.as_bytes(), &mut one, &mut m_one, 1 << 30, 1).unwrap();
        let blocks = read_paf_blocks(paf.as_bytes(), &mut many, &mut m_many, 64, 4).unwrap();
        assert_eq!(whole, blocks);
        let names = |idx: &SequenceIndex| (0..idx.len() as u32).map(|id| (idx.get_name(id).unwrap().to_string(), idx.get_len_from_id(id))).collect::<Vec<_>>();
        assert_eq!(names(&one), names(&many));
        assert_eq!(whole[299].data_offset(), paf.rfind("cg:Z:").unwrap() as u64 + 5);
        // The first error in file order wins, as in a sequential read
        let bad = format!("{paf}q0\t1000\t0\t10\t+\tt7\t2000\t0\t10\t8\t10\t60\tcg:Z:8=\nq0\t5\t0\t10\t+\tt7\t2000\t0\t10\t10\t10\t60\tcg:Z:10=\n");
        let err = read_paf_blocks(bad.as_bytes(), &mut SequenceIndex::new(), &mut 0, 64, 4).unwrap_err();
        assert!(err.to_string().contains("the CIGAR spans 8 query"), "{err}");
        // A blank line, LF or CRLF, is an error
        for blank in ["\n", "\r\n"] {
            let text = format!("{}\n{blank}", paf.lines().next().unwrap());
            let err = read_paf_blocks(text.as_bytes(), &mut SequenceIndex::new(), &mut 0, 64, 4).unwrap_err();
            assert!(err.to_string().contains("is blank"), "{err}");
        }
    }

    #[test]
    fn test_parse_paf_crlf_offsets_stay_absolute() {
        // Regression: `BufRead::lines()` strips `\r\n`, so advancing by
        // `len + 1` drifted every following CIGAR offset by one byte per line
        // and read the wrong bytes back out of the file.
        let paf = "seq1\t100\t0\t50\t+\tseq2\t100\t0\t50\t50\t50\t60\tcg:Z:50=\r\n\
seq1\t100\t50\t100\t+\tseq2\t100\t50\t100\t50\t50\t60\tcg:Z:50=\r\n";
        let mut seq_index = SequenceIndex::new();
        let records = parse_paf(std::io::BufReader::new(paf.as_bytes()), &mut seq_index).unwrap();
        assert_eq!(records.len(), 2);

        let second_offset = records[1].data_offset();
        let expected = paf.match_indices("cg:Z:").nth(1).unwrap().0 + 5;
        assert_eq!(
            second_offset, expected as u64,
            "second CIGAR offset must point at its own `cg:Z:` field"
        );
    }
}
