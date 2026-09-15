use super::{atomic_write, checksum, PanelIdentity, COUNT_POLICY, FORMAT_VERSION};
use crate::sample_mem_bwt::{canonical, encode_walk, invalid, WeightedBwt};
use crate::syng::{SyngIndex, SyngWalkStep};
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader};
use std::path::{Path, PathBuf};

#[derive(Debug, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleStats {
    pub reads: u64,
    pub bases: u64,
    pub reads_with_mems: u64,
    pub mem_records: u64,
    pub distinct_mems: usize,
    pub collection_symbols: usize,
    pub read_lengths: BTreeMap<usize, u64>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SampleIndex {
    pub version: u32,
    pub panel: PanelIdentity,
    pub count_policy: String,
    pub stats: SampleStats,
    pub counts: WeightedBwt,
}

/// The existing matcher chooses its higher-anchor orientation (input wins ties).
/// Invoke it on both input orientations without changing sketch policy. Normalize
/// returned MEMs to input coordinates, deduplicate exact coordinate+node vectors,
/// then remove only exact contiguous coordinate+node subwalks of longer records.
/// Envelope containment alone never removes differing content. Overlaps survive.
/// Return the complete canonical maximal-MEM records retained for one read.
/// Coordinates and read identity are intentionally absent, matching `WeightedBwt`
/// construction exactly. Experimental variable-length scorers should derive
/// node-to-node subwalk keys from these records rather than `observed_pairs`.
pub fn canonical_mem_records(panel: &SyngIndex, sequence: &[u8]) -> io::Result<Vec<Vec<u64>>> {
    collect_tagged_read(panel, sequence)?
        .iter()
        .map(|r| encode_walk(r).map(|t| canonical(&t)))
        .collect()
}

pub(crate) type TaggedWalk = Vec<(i32, u64)>;

pub(crate) fn collect_tagged_read(
    panel: &SyngIndex,
    sequence: &[u8],
) -> io::Result<Vec<TaggedWalk>> {
    let rc = crate::graph::reverse_complement(sequence);
    let mut records = BTreeSet::new();
    for (reverse, seq) in [(false, sequence), (true, rc.as_slice())] {
        let mut matched = panel.matched_syncmers_in_sequence(seq);
        matched.sort_unstable_by_key(|m| m.query_pos);
        let walk: Vec<_> = matched
            .iter()
            .map(|m| SyngWalkStep {
                signed_node: m.signed_node,
                bp_pos: m.query_pos,
            })
            .collect();
        for mem in panel.gbwt_mems_for_walk(&walk)? {
            let mut record: Vec<_> = walk[mem.step_start..mem.step_end]
                .iter()
                .map(|m| (m.signed_node, m.bp_pos))
                .collect();
            if reverse {
                record.reverse();
                for (node, pos) in &mut record {
                    *node = -*node;
                    *pos = sequence.len() as u64 - panel.syncmer_length_bp() as u64 - *pos;
                }
            }
            records.insert(record);
        }
    }
    Ok(maximal_content_records(records))
}

/// Both raw views use increasing input-forward coordinates. Selection reproduces
/// the public matcher for each input orientation (input wins ties). Coordinates
/// remain transient; reverse MEM traversal is queried separately, never interleaved.
pub(crate) fn collect_raw_views(
    panel: &SyngIndex,
    forward: &[(i32, u64)],
    reverse: &[(i32, u64)],
    length: u64,
) -> io::Result<Vec<TaggedWalk>> {
    let mut records = BTreeSet::new();
    for rev in [false, true] {
        let selected = if forward.len() > reverse.len() || (forward.len() == reverse.len() && !rev)
        {
            forward
        } else {
            reverse
        };
        let mut selected = selected.to_vec();
        if rev {
            selected = normalize_reverse(&selected, length, panel.syncmer_length_bp() as u64);
        }
        let walk: Vec<_> = selected
            .iter()
            .map(|&(signed_node, bp_pos)| SyngWalkStep {
                signed_node,
                bp_pos,
            })
            .collect();
        for mem in panel.gbwt_mems_for_walk(&walk)? {
            let mut record = selected[mem.step_start..mem.step_end].to_vec();
            if rev {
                record = normalize_reverse(&record, length, panel.syncmer_length_bp() as u64);
            }
            records.insert(record);
        }
    }
    Ok(maximal_content_records(records))
}

pub(crate) fn normalize_reverse(walk: &[(i32, u64)], length: u64, k: u64) -> TaggedWalk {
    walk.iter()
        .rev()
        .map(|&(n, p)| (-n, length - k - p))
        .collect()
}

fn maximal_content_records(records: BTreeSet<Vec<(i32, u64)>>) -> Vec<Vec<(i32, u64)>> {
    records
        .iter()
        .filter(|record| {
            !records.iter().any(|other| {
                other.len() > record.len()
                    && other.windows(record.len()).any(|w| w == record.as_slice())
            })
        })
        .cloned()
        .collect()
}

pub fn build(
    panel: &SyngIndex,
    identity: PanelIdentity,
    reads: &[PathBuf],
) -> io::Result<SampleIndex> {
    if reads.is_empty() {
        return Err(invalid("at least one reads file is required"));
    }
    let mut stats = SampleStats::default();
    let mut records: BTreeMap<Vec<u64>, u64> = BTreeMap::new();
    for path in reads {
        stream_reads(path, |sequence| {
            stats.reads += 1;
            stats.bases = stats
                .bases
                .checked_add(sequence.len() as u64)
                .ok_or_else(|| invalid("base count overflow"))?;
            *stats.read_lengths.entry(sequence.len()).or_default() += 1;
            let mems = canonical_mem_records(panel, sequence)?;
            stats.reads_with_mems += u64::from(!mems.is_empty());
            for mem in mems {
                stats.mem_records += 1;
                let count = records.entry(mem).or_default();
                *count = count
                    .checked_add(1)
                    .ok_or_else(|| invalid("MEM count overflow"))?;
            }
            Ok(())
        })?;
    }
    let counts = WeightedBwt::build(&records)?;
    stats.distinct_mems = records.len();
    stats.collection_symbols = counts.symbols();
    Ok(SampleIndex {
        version: FORMAT_VERSION,
        panel: identity,
        count_policy: COUNT_POLICY.into(),
        stats,
        counts,
    })
}

impl SampleIndex {
    pub fn save(&self, path: &Path) -> io::Result<()> {
        let payload = bincode::serde::encode_to_vec(self, bincode::config::standard())
            .map_err(io::Error::other)?;
        let mut bytes = b"IMPGMEM1".to_vec();
        bytes.extend((payload.len() as u64).to_le_bytes());
        bytes.extend(checksum(&payload).to_le_bytes());
        bytes.extend(payload);
        atomic_write(path, &bytes)
    }
    pub fn load(path: &Path, panel: &PanelIdentity) -> io::Result<Self> {
        Self::load_with_checksum(path, panel).map(|(sample, _)| sample)
    }
    pub fn load_with_checksum(path: &Path, panel: &PanelIdentity) -> io::Result<(Self, String)> {
        let bytes = fs::read(path)?;
        if bytes.len() < 24 || &bytes[..8] != b"IMPGMEM1" {
            return Err(invalid("sample index header mismatch"));
        }
        let length = u64::from_le_bytes(bytes[8..16].try_into().unwrap());
        let digest = u64::from_le_bytes(bytes[16..24].try_into().unwrap());
        if length != (bytes.len() - 24) as u64 || digest != checksum(&bytes[24..]) {
            return Err(invalid("sample index length/checksum mismatch"));
        }
        let (mut sample, used): (Self, usize) = bincode::serde::decode_from_slice(
            &bytes[24..],
            bincode::config::standard().with_limit::<8589934592>(),
        )
        .map_err(io::Error::other)?;
        if used != bytes.len() - 24
            || sample.version != FORMAT_VERSION
            || sample.count_policy != COUNT_POLICY
            || &sample.panel != panel
        {
            return Err(invalid(
                "incompatible sample version, count policy or panel dictionary",
            ));
        }
        sample.counts.rebuild()?;
        if sample.stats.collection_symbols != sample.counts.symbols() {
            return Err(invalid("sample symbol count mismatch"));
        }
        Ok((sample, format!("{digest:016x}")))
    }
}

/// Streaming FASTA or strict four-line FASTQ, optionally compressed by niffler.
/// Names and quality strings are validated transiently and never retained.
fn stream_reads(path: &Path, mut visit: impl FnMut(&[u8]) -> io::Result<()>) -> io::Result<()> {
    let (reader, _) = niffler::get_reader(Box::new(File::open(path)?)).map_err(io::Error::other)?;
    let mut lines = BufReader::new(reader).lines();
    let first = lines
        .next()
        .transpose()?
        .ok_or_else(|| invalid("empty read file"))?;
    let valid_seq = |seq: &[u8]| -> io::Result<()> {
        if seq.is_empty()
            || seq.len() > i32::MAX as usize
            || seq
                .iter()
                .any(|c| !b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(c))
        {
            return Err(invalid("empty, oversized or invalid DNA read"));
        }
        Ok(())
    };
    let header = |s: &str, prefix: char| -> io::Result<()> {
        if !s.starts_with(prefix) || s[1..].trim().is_empty() {
            Err(invalid("invalid read header"))
        } else {
            Ok(())
        }
    };
    if first.starts_with('>') {
        header(&first, '>')?;
        let mut seq = Vec::new();
        for line in lines {
            let line = line?;
            if line.starts_with('>') {
                valid_seq(&seq)?;
                visit(&seq)?;
                seq.clear();
                header(&line, '>')?;
            } else {
                seq.extend(line.as_bytes());
            }
        }
        valid_seq(&seq)?;
        visit(&seq)?;
    } else {
        let mut current = first;
        loop {
            header(&current, '@')?;
            let seq = lines
                .next()
                .transpose()?
                .ok_or_else(|| invalid("truncated FASTQ sequence"))?;
            let plus = lines
                .next()
                .transpose()?
                .ok_or_else(|| invalid("truncated FASTQ separator"))?;
            let qual = lines
                .next()
                .transpose()?
                .ok_or_else(|| invalid("truncated FASTQ quality"))?;
            valid_seq(seq.as_bytes())?;
            if !plus.starts_with('+')
                || seq.len() != qual.len()
                || qual.bytes().any(|c| !(33..=126).contains(&c))
            {
                return Err(invalid("invalid FASTQ separator or quality"));
            }
            visit(seq.as_bytes())?;
            match lines.next().transpose()? {
                Some(line) => current = line,
                None => break,
            }
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn containment_requires_coordinate_and_content_not_envelope() {
        let long = vec![(1, 0), (2, 5), (3, 10)];
        let distinct = vec![(9, 5)];
        let overlapping = vec![(3, 10), (4, 15)];
        let records = BTreeSet::from([
            long.clone(),
            distinct.clone(),
            vec![(2, 5)],
            overlapping.clone(),
        ]);
        assert_eq!(
            maximal_content_records(records),
            vec![long, overlapping, distinct]
        );
    }
    #[test]
    fn invalid_fastq_is_rejected() {
        let temp = tempfile::tempdir().unwrap();
        for text in [
            "@r\nACGT\n+\n!\n",
            "@r\nACGT\n+\n",
            ">r\nAC?T\n",
            ">r\n>empty\n",
        ] {
            let path = temp.path().join("reads");
            fs::write(&path, text).unwrap();
            assert!(stream_reads(&path, |_| Ok(())).is_err());
        }
    }
}
