use super::*;
use std::io::BufRead;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Seal {
    pub path: String,
    pub bytes: u64,
    pub hash: String,
}
impl Seal {
    pub fn create(root: &Path, path: &str) -> io::Result<Self> {
        let mut reader = BufReader::new(File::open(root.join(path))?);
        let mut hash = 0xcbf29ce484222325;
        let mut bytes = 0;
        let mut buf = vec![0; 1048576];
        loop {
            let n = reader.read(&mut buf)?;
            if n == 0 {
                break;
            }
            super::super::hash_update(&mut hash, &buf[..n]);
            add(&mut bytes, n as u64)?;
        }
        Ok(Self {
            path: path.into(),
            bytes,
            hash: format!("{hash:016x}"),
        })
    }
    pub fn validate(&self) -> io::Result<()> {
        require(
            !self.path.is_empty()
                && Path::new(&self.path)
                    .components()
                    .all(|c| matches!(c, std::path::Component::Normal(_)))
                && self.hash.len() == 16
                && self.hash.bytes().all(|b| b.is_ascii_hexdigit()),
            "invalid route artifact seal/path",
        )
    }
    pub fn verify(&self, root: &Path) -> io::Result<()> {
        self.validate()?;
        let actual = Self::create(root, &self.path)?;
        require(
            self.bytes == actual.bytes && self.hash == actual.hash,
            "route artifact content mismatch",
        )
    }
}
pub fn line(w: &mut impl Write, value: &impl Serialize) -> io::Result<()> {
    serde_json::to_writer(&mut *w, value).map_err(io::Error::other)?;
    w.write_all(b"\n")
}
pub fn next<T: serde::de::DeserializeOwned>(r: &mut impl io::BufRead) -> io::Result<Option<T>> {
    let mut line = String::new();
    let n = r.take(67108865).read_line(&mut line)?;
    require(n <= 67108864, "route row exceeds 64MiB resource cap")?;
    if n == 0 {
        Ok(None)
    } else {
        serde_json::from_str(&line)
            .map(Some)
            .map_err(io::Error::other)
    }
}
/// A disk hub member; equal words form a bucket in the globally sorted fixed-width
/// file. Both orientations reference ONE source resource, never extra copies.
#[derive(Clone, Debug, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Port {
    pub word: Vec<u8>,
    pub source: usize,
    pub anchor: u64,
    pub reverse: bool,
}
impl Port {
    pub fn cut(&self, k: u64) -> u64 {
        self.anchor + if self.reverse { k - k / 2 } else { k / 2 }
    }
    pub fn write(&self, w: &mut impl Write, k: usize) -> io::Result<()> {
        require(self.word.len() == k, "port word width mismatch")?;
        w.write_all(&self.word)?;
        w.write_all(&(self.source as u64).to_le_bytes())?;
        w.write_all(&self.anchor.to_le_bytes())?;
        w.write_all(&[u8::from(self.reverse)])
    }
    pub fn read(r: &mut impl Read, k: usize) -> io::Result<Option<Self>> {
        let mut word = vec![0; k];
        let n = r.read(&mut word[..1])?;
        if n == 0 {
            return Ok(None);
        }
        r.read_exact(&mut word[1..])?;
        let mut src = [0; 8];
        let mut anchor = [0; 8];
        let mut rev = [0];
        r.read_exact(&mut src)?;
        r.read_exact(&mut anchor)?;
        r.read_exact(&mut rev)?;
        require(
            rev[0] <= 1 && word.iter().all(|b| b"ACGT".contains(b)),
            "invalid DNA port record",
        )?;
        Ok(Some(Self {
            word,
            source: usize::try_from(u64::from_le_bytes(src)).map_err(io::Error::other)?,
            anchor: u64::from_le_bytes(anchor),
            reverse: rev[0] == 1,
        }))
    }
}
pub fn port_at(file: &mut File, k: usize, index: u64) -> io::Result<Port> {
    file.seek(SeekFrom::Start(
        index
            .checked_mul(k as u64 + 17)
            .ok_or_else(|| invalid("port offset overflow"))?,
    ))?;
    Port::read(file, k)?.ok_or_else(|| invalid("incomplete port index"))
}
pub fn bucket(file: &mut File, k: usize, count: u64, word: &[u8]) -> io::Result<(u64, u64)> {
    let bound = |file: &mut File, upper: bool| -> io::Result<u64> {
        let (mut lo, mut hi) = (0, count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let p = port_at(file, k, mid)?;
            if p.word.as_slice() < word || (upper && p.word == word) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        Ok(lo)
    };
    Ok((bound(file, false)?, bound(file, true)?))
}
/// External merge with bounded fan-in. Never expand hub membership into edges.
pub fn merge_ports(root: &Path, mut runs: Vec<String>, k: usize) -> io::Result<String> {
    let mut round = 0;
    if runs.is_empty() {
        let path = "ports-empty.bin";
        File::create(root.join(path))?;
        return Ok(path.into());
    }
    while runs.len() > 1 {
        let mut next = Vec::new();
        for (group, chunk) in runs.chunks(32).enumerate() {
            let name = format!("port-merge-{round}-{group}.bin");
            let mut out = BufWriter::new(File::create(root.join(&name))?);
            let mut readers = chunk
                .iter()
                .map(|p| File::open(root.join(p)).map(BufReader::new))
                .collect::<io::Result<Vec<_>>>()?;
            let mut heap = std::collections::BinaryHeap::new();
            for (i, r) in readers.iter_mut().enumerate() {
                if let Some(p) = Port::read(r, k)? {
                    heap.push(std::cmp::Reverse((p, i)));
                }
            }
            while let Some(std::cmp::Reverse((p, i))) = heap.pop() {
                p.write(&mut out, k)?;
                if let Some(p) = Port::read(&mut readers[i], k)? {
                    heap.push(std::cmp::Reverse((p, i)));
                }
            }
            out.flush()?;
            for old in chunk {
                fs::remove_file(root.join(old))?;
            }
            next.push(name);
        }
        runs = next;
        round += 1;
    }
    Ok(runs.pop().unwrap())
}
