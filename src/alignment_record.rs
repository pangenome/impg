/// A generic alignment record that can represent alignments from multiple formats
///
/// ## Field interpretation by format:
/// - **PAF format**:
///   - `strand_and_data_offset`: MSB=strand; remaining bits=file offset to CIGAR string
///   - `data_bytes`: length of CIGAR string in bytes
///
/// - **1aln format**:
///   - `strand_and_data_offset`: MSB=strand; remaining bits=alignment index for O(1) seeking
///   - `data_bytes`: number of tracepoints
#[derive(Debug, PartialEq)]
pub struct AlignmentRecord {
    pub query_id: u32,
    pub query_start: usize,
    pub query_end: usize,
    pub query_contig_offset: i64,
    pub query_contig_len: i64,
    pub target_id: u32,
    pub target_start: usize,
    pub target_end: usize,
    pub target_contig_offset: i64,
    pub target_contig_len: i64,
    pub strand_and_data_offset: u64,
    pub data_bytes: usize,
}

/// Strand orientation for alignments
#[derive(Default, PartialEq, Clone, Copy, Debug)]
#[repr(u8)]
pub enum Strand {
    #[default]
    Forward,
    Reverse,
}

impl AlignmentRecord {
    /// Bit flag for encoding strand in the MSB of strand_and_data_offset
    pub const STRAND_BIT: u64 = 0x8000000000000000;

    /// Get the strand orientation (format-agnostic)
    pub fn strand(&self) -> Strand {
        if (self.strand_and_data_offset & Self::STRAND_BIT) != 0 {
            Strand::Reverse
        } else {
            Strand::Forward
        }
    }

    /// Set the strand orientation (format-agnostic)
    pub fn set_strand(&mut self, strand: Strand) {
        match strand {
            Strand::Forward => self.strand_and_data_offset &= !Self::STRAND_BIT,
            Strand::Reverse => self.strand_and_data_offset |= Self::STRAND_BIT,
        }
    }

    /// Get the data offset without the strand bit
    pub fn data_offset(&self) -> u64 {
        self.strand_and_data_offset & !Self::STRAND_BIT
    }
}

/// Alignment file format types
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum AlignmentFormat {
    Paf,
    OneAln,
}

impl AlignmentFormat {
    /// Detect format from file extension
    pub fn from_path(path: &str) -> Option<Self> {
        if path.ends_with(".paf") || path.ends_with(".paf.gz") || path.ends_with(".paf.bgz") {
            Some(AlignmentFormat::Paf)
        } else if path.ends_with(".1aln") {
            Some(AlignmentFormat::OneAln)
        } else {
            None
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_format_detection() {
        assert_eq!(
            AlignmentFormat::from_path("test.paf"),
            Some(AlignmentFormat::Paf)
        );
        assert_eq!(
            AlignmentFormat::from_path("test.paf.gz"),
            Some(AlignmentFormat::Paf)
        );
        assert_eq!(
            AlignmentFormat::from_path("test.paf.bgz"),
            Some(AlignmentFormat::Paf)
        );
        assert_eq!(
            AlignmentFormat::from_path("test.1aln"),
            Some(AlignmentFormat::OneAln)
        );
        assert_eq!(AlignmentFormat::from_path("test.txt"), None);
    }

    #[test]
    fn test_strand_operations() {
        let mut record = AlignmentRecord {
            query_id: 0,
            query_start: 0,
            query_end: 100,
            query_contig_offset: 0,
            query_contig_len: 100,
            target_id: 1,
            target_start: 0,
            target_end: 100,
            target_contig_offset: 0,
            target_contig_len: 100,
            strand_and_data_offset: 12345,
            data_bytes: 10,
        };

        // Default should be forward
        assert_eq!(record.strand(), Strand::Forward);
        assert_eq!(record.data_offset(), 12345);

        // Set to reverse
        record.set_strand(Strand::Reverse);
        assert_eq!(record.strand(), Strand::Reverse);
        assert_eq!(record.data_offset(), 12345); // Data offset unchanged

        // Set back to forward
        record.set_strand(Strand::Forward);
        assert_eq!(record.strand(), Strand::Forward);
        assert_eq!(record.data_offset(), 12345);
    }
}
