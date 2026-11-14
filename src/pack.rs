use rustc_hash::FxHashMap;
use std::fs::File;
use std::io::{self, BufRead, Read, Seek, SeekFrom, Write};

pub const BINARY_MAGIC: &[u8; 8] = b"IMPGPKB1";
pub const BINARY_VERSION: u32 = 1;
pub const BINARY_HEADER_LEN: u32 = 96;
pub const DEFAULT_BINARY_BLOCK_SIZE: usize = 1 << 20;

#[derive(Debug, Clone, Copy)]
pub struct WriteStats {
    pub retained_records: u64,
    pub syncmer_anchors: u64,
}

#[derive(Debug)]
pub struct Coverage {
    pub counts: FxHashMap<u32, u64>,
    pub universe_nodes: Option<u64>,
    pub nonzero_nodes: u64,
    pub retained_records: Option<u64>,
    pub syncmer_anchors: Option<u64>,
}

pub fn is_binary_format(output_format: &str) -> bool {
    matches!(output_format, "pack" | "packbin" | "pack-bin" | "bpack")
}

pub fn is_tsv_format(output_format: &str) -> bool {
    matches!(output_format, "pack-tsv" | "pack-text" | "packtsv")
}

fn write_u32_le<W: Write>(out: &mut W, value: u32) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}

fn write_i32_le<W: Write>(out: &mut W, value: i32) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}

fn write_u64_le<W: Write>(out: &mut W, value: u64) -> io::Result<()> {
    out.write_all(&value.to_le_bytes())
}

fn read_u32_le<R: Read>(input: &mut R) -> io::Result<u32> {
    let mut bytes = [0u8; 4];
    input.read_exact(&mut bytes)?;
    Ok(u32::from_le_bytes(bytes))
}

fn read_i32_le<R: Read>(input: &mut R) -> io::Result<i32> {
    let mut bytes = [0u8; 4];
    input.read_exact(&mut bytes)?;
    Ok(i32::from_le_bytes(bytes))
}

fn read_u64_le<R: Read>(input: &mut R) -> io::Result<u64> {
    let mut bytes = [0u8; 8];
    input.read_exact(&mut bytes)?;
    Ok(u64::from_le_bytes(bytes))
}

pub fn write_tsv<W: Write>(out: &mut W, counts: FxHashMap<u32, u64>) -> io::Result<usize> {
    writeln!(out, "#node_id\tcount")?;
    let mut rows: Vec<(u32, u64)> = counts.into_iter().filter(|&(_, c)| c > 0).collect();
    rows.sort_unstable_by_key(|&(node_id, _)| node_id);
    for (node_id, count) in &rows {
        writeln!(out, "{node_id}\t{count}")?;
    }
    Ok(rows.len())
}

pub fn write_binary<W: Write>(
    out: &mut W,
    counts: &FxHashMap<u32, u64>,
    stats: WriteStats,
    universe_nodes: usize,
    compression_level: i32,
    block_size: usize,
) -> io::Result<usize> {
    if !(1..=22).contains(&compression_level) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("--pack-compression-level must be in 1..=22, got {compression_level}"),
        ));
    }
    if block_size == 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--pack-block-size must be greater than 0",
        ));
    }
    if block_size > u32::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--pack-block-size must fit in u32",
        ));
    }
    if universe_nodes > u32::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "packbin currently supports at most 2^32 graph features",
        ));
    }

    let mut dense = vec![0u8; universe_nodes];
    let mut overflow = Vec::new();
    let mut nonzero_nodes = 0u64;
    for (&node_id, &count) in counts {
        if node_id == 0 || node_id as usize > universe_nodes {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("pack count node {node_id} is outside 1..={universe_nodes}"),
            ));
        }
        if count == 0 {
            continue;
        }
        nonzero_nodes += 1;
        dense[(node_id - 1) as usize] = count.min(255) as u8;
        if count > 255 {
            overflow.push((node_id, count));
        }
    }
    overflow.sort_unstable_by_key(|&(node_id, _)| node_id);

    let block_count = dense.len().div_ceil(block_size);
    let mut block_offsets = Vec::with_capacity(block_count + 1);
    let mut compressed_blocks = Vec::with_capacity(block_count);
    let mut compressed_offset = 0u64;
    block_offsets.push(compressed_offset);
    for block in dense.chunks(block_size) {
        let compressed = zstd::bulk::compress(block, compression_level)?;
        compressed_offset = compressed_offset
            .checked_add(compressed.len() as u64)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin offset overflow"))?;
        block_offsets.push(compressed_offset);
        compressed_blocks.push(compressed);
    }

    let header_len = u64::from(BINARY_HEADER_LEN);
    let block_index_offset = header_len;
    let block_index_bytes = (block_offsets.len() as u64)
        .checked_mul(8)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin index overflow"))?;
    let overflow_offset = block_index_offset
        .checked_add(block_index_bytes)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin offset overflow"))?;
    let overflow_bytes = (overflow.len() as u64)
        .checked_mul(12)
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin overflow overflow"))?;
    let data_offset = overflow_offset.checked_add(overflow_bytes).ok_or_else(|| {
        io::Error::new(io::ErrorKind::InvalidData, "packbin data offset overflow")
    })?;

    out.write_all(BINARY_MAGIC)?;
    write_u32_le(out, BINARY_VERSION)?;
    write_u32_le(out, BINARY_HEADER_LEN)?;
    write_u64_le(out, universe_nodes as u64)?;
    write_u64_le(out, nonzero_nodes)?;
    write_u64_le(out, stats.retained_records)?;
    write_u64_le(out, stats.syncmer_anchors)?;
    write_u32_le(out, block_size as u32)?;
    write_i32_le(out, compression_level)?;
    write_u64_le(out, block_count as u64)?;
    write_u64_le(out, overflow.len() as u64)?;
    write_u64_le(out, block_index_offset)?;
    write_u64_le(out, overflow_offset)?;
    write_u64_le(out, data_offset)?;

    for offset in block_offsets {
        write_u64_le(out, offset)?;
    }
    for (node_id, count) in overflow {
        write_u32_le(out, node_id)?;
        write_u64_le(out, count)?;
    }
    for block in compressed_blocks {
        out.write_all(&block)?;
    }

    Ok(nonzero_nodes as usize)
}

pub fn read(path: &str) -> io::Result<Coverage> {
    let mut file = File::open(path)?;
    let mut magic = [0u8; 8];
    let n = file.read(&mut magic)?;
    if n == BINARY_MAGIC.len() && &magic == BINARY_MAGIC {
        read_binary(file)
    } else {
        read_tsv(path)
    }
}

fn read_tsv(path: &str) -> io::Result<Coverage> {
    let reader = auto_reader(path)?;
    let mut counts = FxHashMap::default();
    for (line_no, line) in reader.lines().enumerate() {
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let mut fields = line.split('\t');
        let node_id = fields
            .next()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("pack TSV line {} is missing node_id", line_no + 1),
                )
            })?
            .parse::<u32>()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid pack TSV node_id on line {}: {}", line_no + 1, e),
                )
            })?;
        let count = fields
            .next()
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("pack TSV line {} is missing count", line_no + 1),
                )
            })?
            .parse::<u64>()
            .map_err(|e| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("invalid pack TSV count on line {}: {}", line_no + 1, e),
                )
            })?;
        if count > 0 {
            counts.insert(node_id, count);
        }
    }
    Ok(Coverage {
        nonzero_nodes: counts.len() as u64,
        counts,
        universe_nodes: None,
        retained_records: None,
        syncmer_anchors: None,
    })
}

fn read_binary(mut file: File) -> io::Result<Coverage> {
    let version = read_u32_le(&mut file)?;
    if version != BINARY_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported packbin version {}", version),
        ));
    }
    let header_len = read_u32_le(&mut file)?;
    if header_len != BINARY_HEADER_LEN {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "unsupported packbin header length {}; expected {}",
                header_len, BINARY_HEADER_LEN
            ),
        ));
    }
    let universe_nodes = read_u64_le(&mut file)?;
    if universe_nodes > u64::from(u32::MAX) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "packbin universe_nodes exceeds u32::MAX",
        ));
    }
    let universe_usize = universe_nodes as usize;
    let nonzero_nodes = read_u64_le(&mut file)?;
    let retained_records = read_u64_le(&mut file)?;
    let syncmer_anchors = read_u64_le(&mut file)?;
    let block_size = read_u32_le(&mut file)? as usize;
    let _zstd_level = read_i32_le(&mut file)?;
    let block_count_u64 = read_u64_le(&mut file)?;
    let overflow_count_u64 = read_u64_le(&mut file)?;
    let block_index_offset = read_u64_le(&mut file)?;
    let overflow_offset = read_u64_le(&mut file)?;
    let data_offset = read_u64_le(&mut file)?;

    if block_size == 0 && universe_nodes > 0 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "packbin block size is zero",
        ));
    }
    let expected_blocks = if universe_usize == 0 {
        0
    } else {
        universe_usize.div_ceil(block_size)
    };
    let block_count = usize::try_from(block_count_u64)
        .ok()
        .filter(|&c| c == expected_blocks)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "packbin block count {} does not match expected {}",
                    block_count_u64, expected_blocks
                ),
            )
        })?;
    let overflow_count = usize::try_from(overflow_count_u64)
        .ok()
        .filter(|&c| (c as u64) <= universe_nodes)
        .ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "packbin overflow count exceeds universe_nodes",
            )
        })?;

    let block_index_end = (block_count as u64)
        .checked_add(1)
        .and_then(|n| n.checked_mul(8))
        .and_then(|bytes| block_index_offset.checked_add(bytes))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin index overflow"))?;
    let overflow_end = (overflow_count as u64)
        .checked_mul(12)
        .and_then(|bytes| overflow_offset.checked_add(bytes))
        .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin overflow table overflow"))?;
    let file_len = file.metadata()?.len();
    if block_index_end > file_len || overflow_end > file_len {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "packbin table extends past end of file",
        ));
    }

    file.seek(SeekFrom::Start(block_index_offset))?;
    let mut block_offsets = Vec::with_capacity(block_count + 1);
    for _ in 0..=block_count {
        block_offsets.push(read_u64_le(&mut file)?);
    }

    file.seek(SeekFrom::Start(overflow_offset))?;
    let mut overflow = Vec::with_capacity(overflow_count);
    for _ in 0..overflow_count {
        let node_id = read_u32_le(&mut file)?;
        let count = read_u64_le(&mut file)?;
        if node_id == 0 || u64::from(node_id) > universe_nodes || count <= 255 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "invalid packbin overflow entry (node_id={node_id}, count={count}, universe={universe_nodes})"
                ),
            ));
        }
        overflow.push((node_id, count));
    }

    let initial_cap = usize::try_from(nonzero_nodes)
        .unwrap_or(usize::MAX)
        .min(universe_usize)
        .min(1 << 20);
    let mut counts = FxHashMap::with_capacity_and_hasher(initial_cap, Default::default());
    let mut compressed_buf = Vec::new();
    let mut decompressed_buf = Vec::new();
    for block_idx in 0..block_count {
        let start = block_offsets[block_idx];
        let end = block_offsets[block_idx + 1];
        if end < start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packbin block offsets are not monotonic",
            ));
        }
        let compressed_len = usize::try_from(end - start)
            .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "packbin block too large"))?;
        let seek_pos = data_offset
            .checked_add(start)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin offset overflow"))?;
        if seek_pos
            .checked_add(compressed_len as u64)
            .is_none_or(|block_end| block_end > file_len)
        {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packbin block extends past end of file",
            ));
        }
        file.seek(SeekFrom::Start(seek_pos))?;
        compressed_buf.resize(compressed_len, 0);
        file.read_exact(&mut compressed_buf)?;
        let base_idx = block_idx
            .checked_mul(block_size)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "packbin block overflow"))?;
        let remaining = universe_usize.saturating_sub(base_idx);
        let expected_len = remaining.min(block_size);
        decompressed_buf.resize(expected_len, 0);
        let decompressed_len =
            zstd::bulk::decompress_to_buffer(&compressed_buf, &mut decompressed_buf)?;
        if decompressed_len != expected_len {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                "packbin block decompressed to unexpected length",
            ));
        }
        for (offset, &count) in decompressed_buf[..expected_len].iter().enumerate() {
            if count == 0 {
                continue;
            }
            let node_id = (base_idx + offset + 1) as u32;
            counts.insert(node_id, u64::from(count));
        }
    }
    for (node_id, count) in overflow {
        if counts.insert(node_id, count) != Some(255) {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "packbin overflow node {node_id} is duplicate or lacks saturated count 255"
                ),
            ));
        }
    }

    Ok(Coverage {
        counts,
        universe_nodes: Some(universe_nodes),
        nonzero_nodes,
        retained_records: Some(retained_records),
        syncmer_anchors: Some(syncmer_anchors),
    })
}

fn auto_reader(path: &str) -> io::Result<Box<dyn BufRead>> {
    Ok(Box::new(crate::fasta::open_reader(std::path::Path::new(
        path,
    ))?))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn binary_pack_round_trip_handles_overflow_boundary_and_zero_counts() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("sample.pack");
        let mut counts = FxHashMap::default();
        counts.insert(1, 0);
        counts.insert(2, 255);
        counts.insert(3, 256);
        counts.insert(7, u64::MAX);

        let written = {
            let mut f = File::create(&path).unwrap();
            write_binary(
                &mut f,
                &counts,
                WriteStats {
                    retained_records: 12,
                    syncmer_anchors: 34,
                },
                8,
                3,
                3,
            )
            .unwrap()
        };
        assert_eq!(written, 3);

        let loaded = read(path.to_str().unwrap()).unwrap();
        assert_eq!(loaded.nonzero_nodes, 3);
        assert_eq!(loaded.universe_nodes, Some(8));
        assert_eq!(loaded.retained_records, Some(12));
        assert_eq!(loaded.syncmer_anchors, Some(34));
        assert_eq!(loaded.counts.len(), 3);
        assert_eq!(loaded.counts.get(&1), None);
        assert_eq!(loaded.counts.get(&2), Some(&255));
        assert_eq!(loaded.counts.get(&3), Some(&256));
        assert_eq!(loaded.counts.get(&7), Some(&u64::MAX));
    }

    #[test]
    fn binary_pack_rejects_corrupt_overflow_and_block_count() {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("bad.pack");
        let mut counts = FxHashMap::default();
        counts.insert(1, 300);
        {
            let mut f = File::create(&path).unwrap();
            write_binary(
                &mut f,
                &counts,
                WriteStats {
                    retained_records: 1,
                    syncmer_anchors: 1,
                },
                4,
                3,
                4,
            )
            .unwrap();
        }
        // 1. Corrupt overflow entry's count (offset 116) to 100 (<= 255).
        let orig = std::fs::read(&path).unwrap();
        let mut bytes = orig.clone();
        bytes[116..124].copy_from_slice(&100u64.to_le_bytes());
        std::fs::write(&path, &bytes).unwrap();
        assert!(read(path.to_str().unwrap()).is_err());

        // 2. Corrupt overflow entry's node_id (offset 112) to 2 (dense count 0 != 255).
        let mut bytes = orig.clone();
        bytes[112..116].copy_from_slice(&2u32.to_le_bytes());
        std::fs::write(&path, &bytes).unwrap();
        assert!(read(path.to_str().unwrap()).is_err());

        // 3. Corrupt overflow_count (offset 64) to 99 (> universe_nodes).
        let mut bytes = orig.clone();
        bytes[64..72].copy_from_slice(&99u64.to_le_bytes());
        std::fs::write(&path, &bytes).unwrap();
        assert!(read(path.to_str().unwrap()).is_err());

        // 4. Corrupt block_count (offset 56) to 99 (!= expected_blocks).
        let mut bytes = orig.clone();
        bytes[56..64].copy_from_slice(&99u64.to_le_bytes());
        std::fs::write(&path, &bytes).unwrap();
        assert!(read(path.to_str().unwrap()).is_err());

        // 5. Truncate before end of overflow table to verify file-bounds guard.
        std::fs::write(&path, &orig[..118]).unwrap();
        assert!(read(path.to_str().unwrap()).is_err());
    }
}
