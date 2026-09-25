//! Write-ahead journal that lets `impg partition` survive being killed.
//!
//! Why this exists: partition is a serial greedy loop that keeps every
//! partition in memory and writes the BED only at the very end. On a cluster
//! whose large-memory queues cap at 48 h, a run longer than that is killed
//! with nothing on disk -- and there is no resume.
//!
//! What is journaled: for every window, the overlap set that goes INTO
//! `mask_and_update_regions` (after merge and boundary extension). That
//! function, the final re-merge, the stall guard and the next window selection
//! are all deterministic given this input and the current masks, so replaying
//! the journal through the unchanged loop rebuilds exactly the state -- masks,
//! missing regions, partition numbers, collected partitions -- that the killed
//! run had. Only the expensive part, the transitive query, is skipped.
//!
//! What would break it, and what is done about it:
//! - Resuming with different parameters or a different alignment set would
//!   replay one run's masks into another's. The header carries every
//!   parameter that shapes the partition plus a hash of the sequence index
//!   (names and lengths, in id order, since records store ids) and a caller
//!   tag; any difference is a hard error, never a warning.
//! - Window selection that is not reproducible would silently diverge. Every
//!   record carries its window, and replay demands the loop reselect that exact
//!   window; a mismatch is a hard error.
//! - A kill mid-write leaves a torn tail record. Each record ends in a
//!   checksum; replay stops at the first short or corrupt record, the file is
//!   truncated back to the last good record, and those windows are simply
//!   queried again.
//!
//! Enabled only by IMPG_PARTITION_JOURNAL=<path>. Without it, partition runs
//! the upstream code path unchanged.

use coitrees::Interval;
use log::info;
use std::fs::{File, OpenOptions};
use std::io::{self, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use std::time::{Duration, Instant};

const FILE_MAGIC: &[u8; 8] = b"IMPGPJ01";
const REC_MAGIC: u32 = 0x5057_4A31; // "PWJ1"
const SYNC_EVERY: Duration = Duration::from_secs(60);

fn fnv1a(bytes: &[u8], mut h: u64) -> u64 {
    for b in bytes {
        h ^= *b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}
const FNV_INIT: u64 = 0xcbf2_9ce4_8422_2325;

/// Hash of the sequence index in id order: records store numeric ids, so an
/// index whose ids map to other sequences must not be replayed into.
pub fn seq_index_fingerprint(names_lens: impl Iterator<Item = (String, usize)>) -> u64 {
    let mut h = FNV_INIT;
    for (name, len) in names_lens {
        h = fnv1a(name.as_bytes(), h);
        h = fnv1a(&[0], h);
        h = fnv1a(&(len as u64).to_le_bytes(), h);
    }
    h
}

pub type Replayed = Vec<Interval<u32>>;

pub struct Journal {
    path: String,
    reader: Option<BufReader<File>>,
    good_offset: u64,
    writer: Option<BufWriter<File>>,
    last_sync: Instant,
    pub replayed_windows: u64,
    pub appended_windows: u64,
}

impl Journal {
    /// Open (or create) the journal named by IMPG_PARTITION_JOURNAL.
    /// Returns Ok(None) when the variable is unset.
    pub fn from_env(header: &str) -> io::Result<Option<Journal>> {
        let path = match std::env::var("IMPG_PARTITION_JOURNAL") {
            Ok(p) if !p.is_empty() => p,
            _ => return Ok(None),
        };
        let exists = std::fs::metadata(&path).map(|m| m.len() > 0).unwrap_or(false);
        if !exists {
            let mut f = File::create(&path)?;
            f.write_all(FILE_MAGIC)?;
            f.write_all(&(header.len() as u32).to_le_bytes())?;
            f.write_all(header.as_bytes())?;
            f.sync_all()?;
            info!("partition journal: NEW {path}");
            info!("partition journal header: {header}");
            let good = (8 + 4 + header.len()) as u64;
            let writer = OpenOptions::new().append(true).open(&path)?;
            return Ok(Some(Journal {
                path,
                reader: None,
                good_offset: good,
                writer: Some(BufWriter::with_capacity(8 << 20, writer)),
                last_sync: Instant::now(),
                replayed_windows: 0,
                appended_windows: 0,
            }));
        }
        let mut r = BufReader::with_capacity(8 << 20, File::open(&path)?);
        let mut magic = [0u8; 8];
        r.read_exact(&mut magic)?;
        if &magic != FILE_MAGIC {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("{path}: not a partition journal (magic {magic:?})"),
            ));
        }
        let mut len = [0u8; 4];
        r.read_exact(&mut len)?;
        let mut old = vec![0u8; u32::from_le_bytes(len) as usize];
        r.read_exact(&mut old)?;
        let old = String::from_utf8_lossy(&old).into_owned();
        if old != header {
            return Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "{path}: journal was written by a run with different parameters or inputs; \
                     refusing to replay it.\n  journal: {old}\n  this run: {header}"
                ),
            ));
        }
        info!("partition journal: RESUME {path}");
        info!("partition journal header: {header}");
        let good = (8 + 4 + old.len()) as u64;
        Ok(Some(Journal {
            path,
            reader: Some(r),
            good_offset: good,
            writer: None,
            last_sync: Instant::now(),
            replayed_windows: 0,
            appended_windows: 0,
        }))
    }

    /// Next journaled window, if replay is still in progress. The loop's own
    /// choice of window is passed in and must match the journal's.
    pub fn take_replay(&mut self, seq_id: u32, start: i32, end: i32) -> io::Result<Option<Replayed>> {
        let Some(r) = self.reader.as_mut() else {
            return Ok(None);
        };
        match read_record(r)? {
            Some((w, ivs, nbytes)) => {
                if w != (seq_id, start, end) {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!(
                            "{}: replay diverged at record {}: journal has window {:?}, the loop \
                             selected {:?}. Window selection is not reproducible; the journal \
                             cannot be trusted and must not be resumed.",
                            self.path,
                            self.replayed_windows,
                            w,
                            (seq_id, start, end)
                        ),
                    ));
                }
                self.good_offset += nbytes;
                self.replayed_windows += 1;
                Ok(Some(ivs))
            }
            None => {
                self.switch_to_append()?;
                Ok(None)
            }
        }
    }

    fn switch_to_append(&mut self) -> io::Result<()> {
        self.reader = None;
        let f = OpenOptions::new().write(true).open(&self.path)?;
        let len = f.metadata()?.len();
        if len > self.good_offset {
            info!(
                "partition journal: dropping {} bytes of torn tail after record {}",
                len - self.good_offset,
                self.replayed_windows
            );
        }
        f.set_len(self.good_offset)?;
        let mut f = f;
        f.seek(SeekFrom::Start(self.good_offset))?;
        info!(
            "partition journal: replayed {} windows, now querying",
            self.replayed_windows
        );
        self.writer = Some(BufWriter::with_capacity(8 << 20, f));
        Ok(())
    }

    pub fn append(&mut self, seq_id: u32, start: i32, end: i32, ivs: &[Interval<u32>]) -> io::Result<()> {
        if self.reader.is_some() {
            // Replay must be exhausted (take_replay returned None) first.
            return Err(io::Error::new(io::ErrorKind::Other, "journal append during replay"));
        }
        let mut buf = Vec::with_capacity(20 + ivs.len() * 12 + 8);
        buf.extend_from_slice(&REC_MAGIC.to_le_bytes());
        buf.extend_from_slice(&seq_id.to_le_bytes());
        buf.extend_from_slice(&start.to_le_bytes());
        buf.extend_from_slice(&end.to_le_bytes());
        buf.extend_from_slice(&(ivs.len() as u32).to_le_bytes());
        for iv in ivs {
            buf.extend_from_slice(&iv.metadata.to_le_bytes());
            buf.extend_from_slice(&iv.first.to_le_bytes());
            buf.extend_from_slice(&iv.last.to_le_bytes());
        }
        let h = fnv1a(&buf, FNV_INIT);
        buf.extend_from_slice(&h.to_le_bytes());
        let w = self.writer.as_mut().expect("journal writer");
        w.write_all(&buf)?;
        self.appended_windows += 1;
        if self.last_sync.elapsed() >= SYNC_EVERY {
            self.sync()?;
        }
        Ok(())
    }

    pub fn sync(&mut self) -> io::Result<()> {
        if let Some(w) = self.writer.as_mut() {
            w.flush()?;
            w.get_ref().sync_data()?;
        }
        self.last_sync = Instant::now();
        Ok(())
    }

    /// Replay must have consumed the whole journal before the loop can end;
    /// otherwise the journal holds windows the loop never reached, which means
    /// replay diverged in a way the per-window check could not see.
    pub fn finish(&mut self) -> io::Result<()> {
        if let Some(r) = self.reader.as_mut() {
            if read_record(r)?.is_some() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "{}: partitioning finished with journal records left unreplayed",
                        self.path
                    ),
                ));
            }
            self.switch_to_append()?;
        }
        self.sync()?;
        info!(
            "partition journal: {} windows replayed, {} queried this run",
            self.replayed_windows, self.appended_windows
        );
        Ok(())
    }
}

/// One record, or None at EOF / a torn or corrupt record (both end replay).
fn read_record(r: &mut BufReader<File>) -> io::Result<Option<((u32, i32, i32), Replayed, u64)>> {
    let mut head = [0u8; 20];
    if !read_full(r, &mut head)? {
        return Ok(None);
    }
    let u32_at = |b: &[u8], o: usize| u32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    let i32_at = |b: &[u8], o: usize| i32::from_le_bytes(b[o..o + 4].try_into().unwrap());
    if u32_at(&head, 0) != REC_MAGIC {
        return Ok(None);
    }
    let n = u32_at(&head, 16) as usize;
    let mut body = vec![0u8; n * 12];
    if !read_full(r, &mut body)? {
        return Ok(None);
    }
    let mut tail = [0u8; 8];
    if !read_full(r, &mut tail)? {
        return Ok(None);
    }
    let h = fnv1a(&body, fnv1a(&head, FNV_INIT));
    if h != u64::from_le_bytes(tail) {
        return Ok(None);
    }
    let ivs = (0..n)
        .map(|k| Interval {
            metadata: u32_at(&body, k * 12),
            first: i32_at(&body, k * 12 + 4),
            last: i32_at(&body, k * 12 + 8),
        })
        .collect();
    Ok(Some((
        (u32_at(&head, 4), i32_at(&head, 8), i32_at(&head, 12)),
        ivs,
        (20 + n * 12 + 8) as u64,
    )))
}

/// read_exact that reports a clean or torn EOF as `false` instead of an error.
fn read_full(r: &mut BufReader<File>, buf: &mut [u8]) -> io::Result<bool> {
    let mut got = 0;
    while got < buf.len() {
        match r.read(&mut buf[got..]) {
            Ok(0) => return Ok(false),
            Ok(k) => got += k,
            Err(e) if e.kind() == io::ErrorKind::Interrupted => {}
            Err(e) => return Err(e),
        }
    }
    Ok(true)
}
