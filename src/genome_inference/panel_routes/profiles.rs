//! Indexed exact native source-start runs; query domains, not source incidences.
use super::*;
use storage::Seal;
pub type Counts = BTreeMap<[u64; 3], u64>;
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Event {
    pub start: u64,
    pub end: u64,
    pub counts: Vec<([u64; 3], u64)>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct NativeIndex {
    pub length: u64,
    pub starts: u64,
    pub runs: u64,
    pub events: Seal,
    pub index: Seal,
    pub totals: Seal,
}
impl NativeIndex {
    pub fn validate(&self, length: u64, source_length: u64) -> io::Result<()> {
        let starts = if length <= source_length {
            source_length - length + 1
        } else {
            0
        };
        require(
            self.length == length
                && self.starts == starts
                && self.index.bytes
                    == self
                        .runs
                        .checked_mul(16)
                        .ok_or_else(|| invalid("native index size overflow"))?
                && ((starts == 0 && self.runs == 0)
                    || (starts > 0 && self.runs > 0 && self.runs <= starts)),
            "incomplete native start profile",
        )?;
        for s in [&self.events, &self.index, &self.totals] {
            s.validate()?;
        }
        Ok(())
    }
}
pub fn accumulate(
    counts: &mut Counts,
    t: [u64; 3],
    n: u64,
    subtract: bool,
    max_terms: usize,
) -> io::Result<()> {
    if n == 0 {
        return Ok(());
    }
    if subtract {
        let q = counts
            .get_mut(&t)
            .ok_or_else(|| invalid("negative native context correction"))?;
        *q = q
            .checked_sub(n)
            .ok_or_else(|| invalid("negative native context correction"))?;
        if *q == 0 {
            counts.remove(&t);
        }
    } else {
        if !counts.contains_key(&t) {
            require(
                counts.len() < max_terms,
                "native feature resource cap exhausted (no truncation)",
            )?;
        }
        add(counts.entry(t).or_default(), n)?;
    }
    Ok(())
}
/// Both raw views are extracted once per bounded core+halo; their unchanged
/// native union is RC invariant. Each event counts all surviving MEM pairs.
pub fn events(
    panel: &SyngIndex,
    molecule: u64,
    length: u64,
    lo: u64,
    hi: u64,
    core: u64,
    mut fetch: impl FnMut(u64, u64) -> io::Result<Vec<u8>>,
    mut emit: impl FnMut(Event) -> io::Result<()>,
) -> io::Result<()> {
    let starts = if length <= molecule {
        molecule - length + 1
    } else {
        0
    };
    require(
        length > 0 && core > 0 && lo <= hi && hi <= starts,
        "invalid admitted event domain",
    )?;
    let k = panel.syncmer_length_bp() as u64;
    let mut base = lo;
    while base < hi {
        let end = base.saturating_add(core).min(hi);
        let dna = fetch(base, end - 1 + length)?;
        require(
            dna.len() as u64 == end - base - 1 + length,
            "event halo mismatch",
        )?;
        let views = super::super::observations::profile::raw_views(panel, &dna)?;
        let mut cuts = vec![0, end - base];
        for view in &views {
            for &(_, p) in view {
                for e in [(p + k).saturating_sub(length), p + 1] {
                    if e > 0 && e < end - base {
                        cuts.push(e);
                    }
                }
            }
        }
        cuts.sort_unstable();
        cuts.dedup();
        for w in cuts.windows(2) {
            let f = super::super::observations::profile::restrict(&views[0], w[0], length, k);
            let r = super::super::observations::profile::restrict(&views[1], w[0], length, k);
            let mut counts = Counts::new();
            for record in super::super::sample::collect_raw_views(panel, &f, &r, length)? {
                for pair in record.windows(2) {
                    let t = canonical(&encode_walk(pair)?).try_into().unwrap();
                    add(counts.entry(t).or_default(), 1)?;
                }
            }
            emit(Event {
                start: base + w[0],
                end: base + w[1],
                counts: counts.into_iter().collect(),
            })?;
        }
        base = end;
    }
    Ok(())
}
pub fn compile(
    root: &Path,
    source: usize,
    n: u64,
    length: u64,
    panel: &SyngIndex,
    sources: &graph::Sources,
    core: u64,
    max_terms: usize,
) -> io::Result<NativeIndex> {
    let stem = format!("native/{source}-{length}");
    let epath = format!("{stem}.events.jsonl");
    let ipath = format!("{stem}.index");
    let tpath = format!("{stem}.totals.jsonl");
    let mut out = BufWriter::new(File::create(root.join(&epath))?);
    let mut index = BufWriter::new(File::create(root.join(&ipath))?);
    let mut total = Counts::new();
    let mut offset = 0u64;
    let mut runs = 0;
    let mut next_start = 0;
    let mut flush = |e: Event| -> io::Result<()> {
        require(
            e.start == next_start && e.end > e.start,
            "nonpartitioning source events",
        )?;
        next_start = e.end;
        for &(t, q) in &e.counts {
            accumulate(
                &mut total,
                t,
                q.checked_mul(e.end - e.start)
                    .ok_or_else(|| invalid("native total overflow"))?,
                false,
                max_terms,
            )?;
        }
        index.write_all(&e.start.to_le_bytes())?;
        index.write_all(&offset.to_le_bytes())?;
        let bytes = serde_json::to_vec(&e).map_err(io::Error::other)?;
        out.write_all(&bytes)?;
        out.write_all(b"\n")?;
        add(&mut offset, bytes.len() as u64 + 1)?;
        add(&mut runs, 1)
    };
    let starts = if length <= n { n - length + 1 } else { 0 };
    let mut previous: Option<Event> = None;
    events(
        panel,
        n,
        length,
        0,
        starts,
        core,
        |a, b| sources.fetch(source, a, b),
        |e| {
            if let Some(old) = previous.as_mut() {
                if old.end == e.start && old.counts == e.counts {
                    old.end = e.end;
                    return Ok(());
                }
            }
            if let Some(old) = previous.replace(e) {
                flush(old)?;
            }
            Ok(())
        },
    )?;
    if let Some(last) = previous {
        flush(last)?;
    }
    out.flush()?;
    index.flush()?;
    require(next_start == starts, "incomplete compiled start domain")?;
    let mut writer = BufWriter::new(File::create(root.join(&tpath))?);
    for row in total {
        storage::line(&mut writer, &row)?;
    }
    writer.flush()?;
    Ok(NativeIndex {
        length,
        starts,
        runs,
        events: Seal::create(root, &epath)?,
        index: Seal::create(root, &ipath)?,
        totals: Seal::create(root, &tpath)?,
    })
}
pub struct Store {
    root: PathBuf,
    verified: BTreeSet<(usize, usize)>,
    pub queried_starts: u64,
}
impl Store {
    pub fn new(root: &Path) -> Self {
        Self {
            root: root.into(),
            verified: BTreeSet::new(),
            queried_starts: 0,
        }
    }
    /// The sorted start+offset index seeks to the first intersecting event. Full
    /// domains use persisted exact totals. No read-location matrix or prefix map.
    pub fn integrate(
        &mut self,
        g: &Graph,
        source: usize,
        li: usize,
        lo: u64,
        hi: u64,
        counts: &mut Counts,
        subtract: bool,
        max_terms: usize,
    ) -> io::Result<()> {
        let p = &g.lanes[source].profiles[li];
        require(lo <= hi && hi <= p.starts, "native start-range mismatch")?;
        if lo == hi {
            return Ok(());
        }
        add(&mut self.queried_starts, hi - lo)?;
        if !self.verified.contains(&(source, li)) {
            for s in [&p.events, &p.index, &p.totals] {
                s.verify(&self.root)?;
            }
            self.verified.insert((source, li));
        }
        if lo == 0 && hi == p.starts {
            let mut r = BufReader::new(File::open(self.root.join(&p.totals.path))?);
            let mut last = None;
            while let Some((t, q)) = storage::next::<([u64; 3], u64)>(&mut r)? {
                require(
                    q > 0 && last.is_none_or(|old| old < t),
                    "invalid native totals order",
                )?;
                last = Some(t);
                crate::sample_mem_bwt::validate_pattern(&t)?;
                require(canonical(&t) == t, "noncanonical native total")?;
                accumulate(counts, t, q, subtract, max_terms)?;
            }
            return Ok(());
        }
        let mut index = File::open(self.root.join(&p.index.path))?;
        let row = |f: &mut File, i: u64| -> io::Result<(u64, u64)> {
            f.seek(SeekFrom::Start(i * 16))?;
            let mut a = [0; 8];
            let mut b = [0; 8];
            f.read_exact(&mut a)?;
            f.read_exact(&mut b)?;
            Ok((u64::from_le_bytes(a), u64::from_le_bytes(b)))
        };
        let (mut left, mut right) = (0, p.runs);
        while left < right {
            let mid = left + (right - left) / 2;
            if row(&mut index, mid)?.0 <= lo {
                left = mid + 1
            } else {
                right = mid
            }
        }
        require(left > 0, "missing native range index")?;
        let (_, offset) = row(&mut index, left - 1)?;
        let mut reader = BufReader::new(File::open(self.root.join(&p.events.path))?);
        reader.seek(SeekFrom::Start(offset))?;
        let mut next = lo;
        while next < hi {
            let e: Event = storage::next(&mut reader)?
                .ok_or_else(|| invalid("incomplete native event stream"))?;
            require(
                e.start <= next && e.end > next && e.end <= p.starts,
                "nonpartitioning native range",
            )?;
            let end = e.end.min(hi);
            let mut last = None;
            for (t, q) in e.counts {
                crate::sample_mem_bwt::validate_pattern(&t)?;
                require(
                    q > 0 && canonical(&t) == t && last.is_none_or(|old| old < t),
                    "invalid native event counts",
                )?;
                last = Some(t);
                accumulate(
                    counts,
                    t,
                    q.checked_mul(end - next)
                        .ok_or_else(|| invalid("restricted profile overflow"))?,
                    subtract,
                    max_terms,
                )?;
            }
            next = end;
        }
        Ok(())
    }
}
