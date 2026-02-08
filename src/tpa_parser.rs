//! TPA format parsing
//!
//! Handles parsing of TracePoint Alignment format (.tpa files) using the tpa crate.

use crate::alignment_record::{AlignmentRecord, Strand};
use crate::onealn::OneAlnAlignment;
use crate::seqidx::SequenceIndex;
use log::debug;
use tpa::{TracepointData, TpaReader};

/// TPA file parser
pub struct TpaParser {
    file_path: String,
}

/// Error type for TPA parsing
#[derive(Debug)]
pub enum ParseErr {
    InvalidFormat(String),
}

impl std::fmt::Display for ParseErr {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseErr::InvalidFormat(msg) => write!(f, "{}", msg),
        }
    }
}

impl std::error::Error for ParseErr {}

impl TpaParser {
    pub fn new(file_path: String) -> Result<Self, ParseErr> {
        Ok(TpaParser { file_path })
    }

    /// Read trace spacing from TPA header (max_complexity)
    pub fn read_trace_spacing(file_path: &str) -> Result<i64, ParseErr> {
        let reader = TpaReader::new(file_path).map_err(|e| {
            ParseErr::InvalidFormat(format!("Failed to open TPA file '{}': {}", file_path, e))
        })?;
        Ok(reader.header().max_complexity() as i64)
    }

    /// Parse all alignments from TPA file into AlignmentRecords.
    /// `threads` controls BGZF decompression parallelism (1 = single-threaded).
    pub fn parse_alignments(
        &self,
        seq_index: &mut SequenceIndex,
        threads: usize,
    ) -> Result<Vec<AlignmentRecord>, ParseErr> {
        let mut reader = TpaReader::new(&self.file_path).map_err(|e| {
            ParseErr::InvalidFormat(format!("Failed to open TPA file: {}", e))
        })?;

        reader.load_string_table().map_err(|e| {
            ParseErr::InvalidFormat(format!("Failed to load string table: {}", e))
        })?;

        let num_records = reader.len();
        let mut records = Vec::with_capacity(num_records);

        // Use sequential iteration (much faster than per-record random access)
        let string_table = reader.string_table().map_err(|e| {
            ParseErr::InvalidFormat(format!("Failed to get string table: {}", e))
        })?;

        // Pre-register all sequences from the string table to avoid lookups in the loop
        let mut name_id_to_seq_id: Vec<Option<(u32, usize)>> = vec![None; string_table.len()];
        for name_id in 0..string_table.len() as u64 {
            if let Some((name, len)) = string_table.get_name_and_len(name_id) {
                let seq_id = seq_index.get_or_insert_id(name, Some(len as usize));
                name_id_to_seq_id[name_id as usize] = Some((seq_id, len as usize));
            }
        }

        let mut record_id: u64 = 0;
        let metadata_iter = reader.iter_record_metadata(threads).map_err(|e| {
            ParseErr::InvalidFormat(format!("Failed to create metadata iterator: {}", e))
        })?;
        for result in metadata_iter {
            let compact = result.map_err(|e| {
                ParseErr::InvalidFormat(format!(
                    "Failed to read record metadata {}: {}",
                    record_id, e
                ))
            })?;

            let (query_id, _) =
                name_id_to_seq_id[compact.query_name_id as usize].ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Query name ID {} not found in string table",
                        compact.query_name_id
                    ))
                })?;
            let (target_id, _) =
                name_id_to_seq_id[compact.target_name_id as usize].ok_or_else(|| {
                    ParseErr::InvalidFormat(format!(
                        "Target name ID {} not found in string table",
                        compact.target_name_id
                    ))
                })?;

            let strand = if compact.strand == '-' {
                Strand::Reverse
            } else {
                Strand::Forward
            };

            // For TPA, coordinates are already in scaffold space (like PAF)
            let query_start = compact.query_start as usize;
            let query_end = compact.query_end as usize;
            let target_start = compact.target_start as usize;
            let target_end = compact.target_end as usize;

            let mut record = AlignmentRecord {
                query_id,
                query_start,
                query_end,
                target_id,
                target_start,
                target_end,
                strand_and_data_offset: record_id, // Store record ID for O(1) seeking
                data_bytes: 0, // Tracepoints not loaded during indexing
            };
            record.set_strand(strand);

            records.push(record);
            record_id += 1;
        }

        debug!(
            "Parsed {} alignment records from TPA file: {}",
            records.len(),
            self.file_path
        );

        Ok(records)
    }

    /// Fetch a single alignment from a TPA reader, converting to OneAlnAlignment
    pub fn fetch_alignment(
        reader: &mut TpaReader,
        record_id: u64,
    ) -> Result<OneAlnAlignment, ParseErr> {
        let header = reader.header();
        let tp_type = header.tp_type();
        let complexity_metric = header.complexity_metric();
        let max_complexity = header.max_complexity();

        let compact = reader.get_compact_record(record_id).map_err(|e| {
            ParseErr::InvalidFormat(format!(
                "Failed to read compact record {}: {}",
                record_id, e
            ))
        })?;

        let strand = compact.strand;

        // TPA coordinates are scaffold coordinates (like PAF), so contig_offset = 0
        let query_start = compact.query_start as i64;
        let query_end = compact.query_end as i64;
        let target_start = compact.target_start as i64;
        let target_end = compact.target_end as i64;

        let mut tracepoints = Vec::new();
        let mut trace_diffs = Vec::new();
        let mut query_deltas = Vec::new();
        let trace_spacing = max_complexity as i64;

        match compact.tracepoints {
            TracepointData::Fastga(pairs) => {
                // Fastga: (num_diffs, target_delta)
                for (diffs, target_delta) in pairs {
                    tracepoints.push(target_delta as i64);
                    trace_diffs.push(diffs as i64);
                }
                // trace_spacing = max_complexity (fixed query spacing)
            }
            TracepointData::Standard(pairs) => {
                // Standard: (query_delta, target_delta)
                for (qd, td) in pairs {
                    tracepoints.push(td as i64);
                    query_deltas.push(qd as i64);
                }
                // trace_spacing = max_complexity (used as max_value for banded alignment)
            }
            TracepointData::Variable(pairs) => {
                // Variable: (first_val, Option<second_val>)
                for (first, second) in pairs {
                    tracepoints.push(first as i64);
                    if let Some(s) = second {
                        query_deltas.push(s as i64);
                    }
                }
            }
            TracepointData::Mixed(items) => {
                for item in items {
                    match item {
                        tpa::MixedRepresentation::Tracepoint(a, b) => {
                            tracepoints.push(b as i64);
                            query_deltas.push(a as i64);
                        }
                        tpa::MixedRepresentation::CigarOp(_, _) => {
                            // Skip CIGAR ops in mixed mode
                        }
                    }
                }
            }
        }

        let is_fastga = matches!(tp_type, tpa::TracepointType::Fastga);

        Ok(OneAlnAlignment {
            query_name: String::new(),
            query_length: 0,
            query_contig_start: query_start,
            query_contig_end: query_end,
            query_contig_offset: 0, // PAF-style coords = scaffold coords
            target_name: String::new(),
            target_length: 0,
            target_contig_start: target_start,
            target_contig_end: target_end,
            target_contig_offset: 0,
            target_contig_len: 0,
            strand,
            differences: 0,
            tracepoints,
            trace_diffs,
            trace_spacing,
            query_deltas,
            complexity_metric: if !is_fastga {
                Some(complexity_metric)
            } else {
                None
            },
            max_complexity: if !is_fastga {
                Some(max_complexity)
            } else {
                None
            },
        })
    }
}
