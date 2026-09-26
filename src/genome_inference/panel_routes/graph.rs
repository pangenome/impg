use super::*;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use rustc_hash::FxHashMap;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::{Arc, Mutex, OnceLock};
use storage::{Port, Seal};
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Lane {
    pub id: usize,
    pub name: String,
    pub length: u64,
    pub sequence_fnv1a64: String,
    pub family: usize,
    pub ports: Seal,
    pub port_count: u64,
    pub profiles: Vec<profiles::NativeIndex>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Family {
    pub identity: String,
    pub paths: Vec<usize>,
}
/// Auxiliary files consumed by source lookup, distinct from sequence data inputs.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SourceAccess {
    pub source_file: usize,
    pub path: String,
    pub fingerprint: serde_json::Value,
}

pub fn source_access_paths(paths: &[String]) -> io::Result<Vec<(usize, String)>> {
    let mut result = Vec::new();
    for (i, path) in paths.iter().enumerate() {
        if Path::new(path).extension().is_some_and(|x| x == "agc") {
            continue;
        }
        require(
            [".fa", ".fasta", ".fna", ".fa.gz", ".fasta.gz", ".fna.gz"]
                .iter()
                .any(|s| path.ends_with(s)),
            "unsupported source access format",
        )?;
        result.push((i, format!("{path}.fai")));
        let mut magic = [0; 2];
        File::open(path)?.read_exact(&mut magic)?;
        if magic == [0x1f, 0x8b] {
            result.push((i, format!("{path}.gzi")));
        }
    }
    Ok(result)
}

/// Build-only initialization: HTSlib may create FAI/GZI here, never during reuse.
pub fn bind_source_access(paths: &[String]) -> io::Result<Vec<SourceAccess>> {
    for path in paths {
        if !Path::new(path).extension().is_some_and(|x| x == "agc") {
            rust_htslib::faidx::Reader::from_path(path).map_err(io::Error::other)?;
        }
    }
    source_access_paths(paths)?
        .into_iter()
        .map(|(source_file, path)| {
            let fingerprint = super::super::reconstruction::fingerprint(Path::new(&path))?;
            Ok(SourceAccess {
                source_file,
                path,
                fingerprint,
            })
        })
        .collect()
}

pub fn verify_source_access(paths: &[String], bindings: &[SourceAccess]) -> io::Result<()> {
    let expected = source_access_paths(paths)?;
    require(
        expected.len() == bindings.len(),
        "incomplete source-access sidecar bindings",
    )?;
    for ((i, path), binding) in expected.iter().zip(bindings) {
        require(
            *i == binding.source_file && *path == binding.path,
            "source-access sidecar identity mismatch",
        )?;
        require(
            binding.fingerprint == super::super::reconstruction::fingerprint(Path::new(path))?,
            "source-access sidecar content mismatch",
        )?;
    }
    Ok(())
}

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Graph {
    pub version: u32,
    pub model: String,
    pub panel: PanelIdentity,
    pub compiler_identity: String,
    pub count_policy: String,
    pub generation_rule: String,
    pub generation_complete: bool,
    pub endpoint_semantics: String,
    pub k: u64,
    pub cut_offset: u64,
    pub core_bp: u64,
    pub read_lengths: Vec<u64>,
    pub lanes: Vec<Lane>,
    pub families: Vec<Family>,
    pub ports: Seal,
    pub port_count: u64,
    pub ownership: Seal,
    pub original_ownership: Seal,
    pub registry: Seal,
    pub registry_count: u64,
    pub catalog_provenance: serde_json::Value,
    pub source_files: Vec<serde_json::Value>,
    pub source_paths: Vec<String>,
    pub source_access: Vec<SourceAccess>,
    pub total_source_bp: u64,
}
impl Graph {
    pub fn digest(&self) -> io::Result<String> {
        digest(self)
    }
    pub fn validate(&self) -> io::Result<()> {
        require(self.version==VERSION && self.model==MODEL && self.count_policy==COUNT_POLICY && self.compiler_identity==compiler_identity() && self.generation_rule==RULE && self.generation_complete && self.endpoint_semantics=="panel-identity-native-linear-assembly-path-end-pairings-not-biological-termini","incompatible/incomplete route graph")?;
        require(
            self.k >= 2
                && self.k <= 1048576
                && self.cut_offset == self.k / 2
                && self.core_bp > 0
                && self.core_bp <= 1048576
                && !self.read_lengths.is_empty()
                && self.read_lengths.len() <= 32
                && self.read_lengths.iter().all(|&l| l > 0 && l <= 1048576)
                && self.read_lengths.windows(2).all(|w| w[0] < w[1]),
            "invalid native route parameters",
        )?;
        require(
            !self.lanes.is_empty()
                && !self.families.is_empty()
                && self.source_paths.len() == self.source_files.len()
                && !self.source_paths.is_empty(),
            "empty route graph/sources",
        )?;
        let mut names = BTreeSet::new();
        let mut bp = 0;
        let mut ports = 0;
        for (id, l) in self.lanes.iter().enumerate() {
            require(
                l.id == id
                    && names.insert(&l.name)
                    && l.length > 0
                    && l.length <= i32::MAX as u64
                    && l.sequence_fnv1a64.len() == 16
                    && l.sequence_fnv1a64.bytes().all(|b| b.is_ascii_hexdigit())
                    && l.family < self.families.len()
                    && l.profiles.len() == self.read_lengths.len(),
                "invalid native lane inventory",
            )?;
            require(
                super::super::genotype::source_identity(&l.name)?
                    == self.families[l.family].identity,
                "lane identity mismatch",
            )?;
            l.ports.validate()?;
            require(
                l.ports.bytes
                    == l.port_count
                        .checked_mul(self.k + 17)
                        .ok_or_else(|| invalid("port size overflow"))?
                    && l.port_count % 2 == 0,
                "invalid physical port count",
            )?;
            for (&length, p) in self.read_lengths.iter().zip(&l.profiles) {
                p.validate(length, l.length)?;
            }
            add(&mut bp, l.length)?;
            add(&mut ports, l.port_count)?;
        }
        let mut membership = BTreeSet::new();
        let mut identities = BTreeSet::new();
        for (f, family) in self.families.iter().enumerate() {
            require(
                !family.paths.is_empty() && identities.insert(&family.identity),
                "duplicate/empty topology family",
            )?;
            for &s in &family.paths {
                require(
                    s < self.lanes.len() && self.lanes[s].family == f && membership.insert(s),
                    "invalid/duplicate topology path inventory",
                )?;
            }
        }
        require(
            membership.len() == self.lanes.len()
                && bp == self.total_source_bp
                && ports == self.port_count
                && self.ports.bytes
                    == ports
                        .checked_mul(self.k + 17)
                        .ok_or_else(|| invalid("port size overflow"))?,
            "incomplete native path/port inventory",
        )?;
        for s in [
            &self.ports,
            &self.ownership,
            &self.original_ownership,
            &self.registry,
        ] {
            s.validate()?;
        }
        Ok(())
    }
    pub fn load(root: &Path, identity: &PanelIdentity) -> io::Result<Self> {
        let manifest: serde_json::Value = read_json(&root.join("manifest.json"))?;
        require(
            manifest["status"] == "succeeded" && manifest["model"] == MODEL,
            "route graph is incomplete/failed/incompatible",
        )?;
        let seal: Seal =
            serde_json::from_value(manifest["stats"]["graph"].clone()).map_err(io::Error::other)?;
        seal.verify(root)?;
        let graph: Self = read_json(&root.join("graph.json"))?;
        graph.validate()?;
        require(graph.panel == *identity, "route panel identity mismatch")?;
        for s in [
            &graph.ports,
            &graph.ownership,
            &graph.original_ownership,
            &graph.registry,
        ] {
            s.verify(root)?;
        }
        Ok(graph)
    }
    pub fn native_assignment(&self, family: usize) -> io::Result<Assignment> {
        require(family < self.families.len(), "unknown topology family")?;
        Ok(Assignment {
            version: VERSION,
            model: MODEL.into(),
            graph_checksum: self.digest()?,
            family,
            routes: self.families[family]
                .paths
                .iter()
                .map(|&s| Route {
                    segments: vec![Segment {
                        source: s,
                        start: 0,
                        end: self.lanes[s].length,
                        reverse: false,
                    }],
                })
                .collect(),
        })
    }
}
/// Probe counters for the `Sources::fetch` memo (IMPG_DB only): total
/// unmemoized fetches and memo hits, printed by callers that quote the
/// boundary-endpoint stage.
static PROBE_FETCH_CALLS: AtomicU64 = AtomicU64::new(0);
static PROBE_FETCH_MEMO_HITS: AtomicU64 = AtomicU64::new(0);
static PROBE_FETCH_ON: OnceLock<bool> = OnceLock::new();
fn probe_sources_fetch() -> bool {
    *PROBE_FETCH_ON.get_or_init(|| std::env::var("IMPG_DB").is_ok())
}
/// (unmemoized fetches, memo hits) since process start; IMPG_DB probes only.
pub fn sources_fetch_probe_stats() -> (u64, u64) {
    (
        PROBE_FETCH_CALLS.load(Ordering::Relaxed),
        PROBE_FETCH_MEMO_HITS.load(Ordering::Relaxed),
    )
}

pub struct Sources {
    indexes: Vec<UnifiedSequenceIndex>,
    chosen: Vec<usize>,
    pub lanes: Vec<(String, u64)>,
    /// Memoized answers to small-crop `fetch` queries, keyed by the public
    /// query (source, start, end). The boundary-endpoint derivation issues
    /// ~100k short flank fetches through this path with heavily repeated
    /// coordinates (shared segment heads across alleles of a partition);
    /// source files are immutable during a run, so a memoized answer equals
    /// the recomputed one bit for bit (the same argument as the port
    /// `at_cut` memo). Large crops bypass the memo to bound its memory.
    fetch_memo: Mutex<FxHashMap<(usize, u64, u64), Arc<Vec<u8>>>>,
}
impl Sources {
    pub fn open(paths: &[String], lanes: Vec<(String, u64)>) -> io::Result<Self> {
        let indexes = paths
            .iter()
            .map(|p| UnifiedSequenceIndex::from_files(&[p.clone()]))
            .collect::<io::Result<Vec<_>>>()?;
        let mut chosen = Vec::new();
        for (name, length) in &lanes {
            require(
                *length > 0 && *length <= i32::MAX as u64,
                "unsupported source length",
            )?;
            let matches: Vec<_> = indexes
                .iter()
                .enumerate()
                .filter_map(|(i, x)| x.get_sequence_length(name).ok().map(|l| (i, l as u64)))
                .collect();
            require(
                matches.len() == 1 && matches[0].1 == *length,
                "missing/ambiguous source name or length binding",
            )?;
            chosen.push(matches[0].0);
        }
        Ok(Self {
            indexes,
            chosen,
            lanes,
            fetch_memo: Mutex::new(FxHashMap::default()),
        })
    }
    pub fn fetch(&self, source: usize, start: u64, end: u64) -> io::Result<Vec<u8>> {
        require(
            source < self.lanes.len() && start < end && end <= self.lanes[source].1,
            "invalid source crop",
        )?;
        // Only small crops are memoized: the boundary flanks are ~150 bp and
        // repeat heavily, while full-sequence spells are large and mostly
        // distinct (memoizing them would waste RSS for no hits).
        const FETCH_MEMO_MAX_CROP: u64 = 1024;
        let memoizable = end - start <= FETCH_MEMO_MAX_CROP;
        let key = (source, start, end);
        if memoizable {
            if let Ok(memo) = self.fetch_memo.lock() {
                if let Some(dna) = memo.get(&key) {
                    if probe_sources_fetch() {
                        PROBE_FETCH_MEMO_HITS.fetch_add(1, Ordering::Relaxed);
                    }
                    return Ok(dna.as_ref().clone());
                }
            }
        }
        if probe_sources_fetch() {
            PROBE_FETCH_CALLS.fetch_add(1, Ordering::Relaxed);
        }
        let mut dna = self.indexes[self.chosen[source]].fetch_sequence(
            &self.lanes[source].0,
            start as i32,
            end as i32,
        )?;
        require(
            dna.len() as u64 == end - start
                && dna
                    .iter()
                    .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b)),
            "invalid source DNA",
        )?;
        dna.make_ascii_uppercase();
        if memoizable {
            if let Ok(mut memo) = self.fetch_memo.lock() {
                memo.insert(key, Arc::new(dna.clone()));
            }
        }
        Ok(dna)
    }
    pub fn spell(&self, route: &Route, lo: u64, hi: u64) -> io::Result<Vec<u8>> {
        require(lo < hi && hi <= route.length()?, "invalid route crop")?;
        let mut out = Vec::new();
        let mut offset = 0;
        for p in &route.segments {
            let end = offset + p.end - p.start;
            let a = lo.max(offset);
            let b = hi.min(end);
            if a < b {
                let (start, stop) = if p.reverse {
                    (p.end - (b - offset), p.end - (a - offset))
                } else {
                    (p.start + a - offset, p.start + b - offset)
                };
                let mut dna = self.fetch(p.source, start, stop)?;
                if p.reverse {
                    dna = crate::graph::reverse_complement(&dna);
                }
                out.extend(dna);
            }
            offset = end;
            if offset >= hi {
                break;
            }
        }
        require(out.len() as u64 == hi - lo, "incomplete route spelling")?;
        Ok(out)
    }
}
/// Verified native port indexes, lazily opened and checked, with implicit on-disk
/// binary-search hubs. No all-pairs or all-members RAM structure is constructed.
pub struct Ports {
    root: PathBuf,
    verified: BTreeSet<usize>,
    /// Open per-source port file handles, reused across lookups: re-opening
    /// the file on every `at_cut`/`source_file` call dominated the seam
    /// derivation stage (~100k opens on a tract slice). The port files are
    /// immutable during a run, so handle reuse is behavior-identical.
    files: BTreeMap<usize, File>,
    /// Memoized `at_cut` answers keyed by the public query
    /// (source, cut, reverse). Port files are immutable during a run, so the
    /// memoized answer equals the recomputed one bit for bit.
    at_cut_memo: BTreeMap<(usize, u64, bool), Option<Port>>,
    pub global: File,
}
impl Ports {
    pub fn open(root: &Path, g: &Graph) -> io::Result<Self> {
        g.ports.verify(root)?;
        Ok(Self {
            root: root.into(),
            verified: BTreeSet::new(),
            files: BTreeMap::new(),
            at_cut_memo: BTreeMap::new(),
            global: File::open(root.join(&g.ports.path))?,
        })
    }
    /// Component-local open for subrange smoke runs: skips the whole-artifact
    /// global ports verification (an O(35 GB) FNV re-read) while keeping the
    /// per-source lazy verification that actually guards seam lookups. The
    /// global file is used only by the route search engine, never by the
    /// genome smoke's per-cut `at_cut`/`forward_ports_inside` lookups.
    pub fn open_without_global_verification(root: &Path, g: &Graph) -> io::Result<Self> {
        Ok(Self {
            root: root.into(),
            verified: BTreeSet::new(),
            files: BTreeMap::new(),
            at_cut_memo: BTreeMap::new(),
            global: File::open(root.join(&g.ports.path))?,
        })
    }
    /// Cached per-source port file handle: verifies the source's port index
    /// once, then keeps the handle open for reuse.
    fn source_handle(&mut self, g: &Graph, source: usize) -> io::Result<&mut File> {
        if !self.verified.contains(&source) {
            g.lanes[source].ports.verify(&self.root)?;
            self.verified.insert(source);
        }
        let path = self.root.join(&g.lanes[source].ports.path);
        if let std::collections::btree_map::Entry::Vacant(entry) = self.files.entry(source) {
            entry.insert(File::open(&path)?);
        }
        Ok(self
            .files
            .get_mut(&source)
            .expect("per-source port handle just inserted"))
    }
    /// A duplicate of the cached per-source port file handle, positioned at
    /// offset 0 so sequential readers see a fresh file exactly as the old
    /// open-per-call behavior provided (dup vs open — same bytes; every
    /// binary-search reader re-seeks explicitly anyway).
    pub fn source_file(&mut self, g: &Graph, source: usize) -> io::Result<File> {
        let mut file = self.source_handle(g, source)?.try_clone()?;
        file.seek(SeekFrom::Start(0))?;
        Ok(file)
    }
    pub fn forward_ports_inside(
        &mut self,
        g: &Graph,
        source: usize,
        start: u64,
        end: u64,
    ) -> io::Result<Vec<Port>> {
        require(
            source < g.lanes.len() && start < end && end <= g.lanes[source].length,
            "invalid source port interval",
        )?;
        let mut file = self.source_file(g, source)?;
        let first_anchor = start.saturating_add(1).saturating_sub(g.cut_offset);
        let (mut lo, mut hi) = (0, g.lanes[source].port_count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let p = storage::port_at(&mut file, g.k as usize, mid)?;
            if (p.anchor, p.reverse) < (first_anchor, false) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        let mut result = Vec::new();
        while lo < g.lanes[source].port_count {
            let p = storage::port_at(&mut file, g.k as usize, lo)?;
            let cut = p.anchor.saturating_add(g.cut_offset);
            if cut >= end {
                break;
            }
            if !p.reverse && start < cut {
                result.push(p);
            }
            lo += 1;
        }
        Ok(result)
    }

    pub fn at_cut(
        &mut self,
        g: &Graph,
        source: usize,
        cut: u64,
        reverse: bool,
    ) -> io::Result<Option<Port>> {
        let key = (source, cut, reverse);
        if let Some(answer) = self.at_cut_memo.get(&key) {
            return Ok(answer.clone());
        }
        let shift = if reverse {
            g.k - g.cut_offset
        } else {
            g.cut_offset
        };
        let Some(anchor) = cut.checked_sub(shift) else {
            self.at_cut_memo.insert(key, None);
            return Ok(None);
        };
        let file = self.source_handle(g, source)?;
        let (mut lo, mut hi) = (0, g.lanes[source].port_count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let p = storage::port_at(file, g.k as usize, mid)?;
            if (p.anchor, p.reverse) < (anchor, reverse) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        let answer = if lo == g.lanes[source].port_count {
            None
        } else {
            let p = storage::port_at(file, g.k as usize, lo)?;
            ((p.anchor, p.reverse) == (anchor, reverse)).then_some(p)
        };
        self.at_cut_memo.insert(key, answer.clone());
        Ok(answer)
    }
}
