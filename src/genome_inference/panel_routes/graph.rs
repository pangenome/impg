use super::*;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
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
pub struct Sources {
    indexes: Vec<UnifiedSequenceIndex>,
    chosen: Vec<usize>,
    pub lanes: Vec<(String, u64)>,
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
        })
    }
    pub fn fetch(&self, source: usize, start: u64, end: u64) -> io::Result<Vec<u8>> {
        require(
            source < self.lanes.len() && start < end && end <= self.lanes[source].1,
            "invalid source crop",
        )?;
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
    pub global: File,
}
impl Ports {
    pub fn open(root: &Path, g: &Graph) -> io::Result<Self> {
        g.ports.verify(root)?;
        Ok(Self {
            root: root.into(),
            verified: BTreeSet::new(),
            global: File::open(root.join(&g.ports.path))?,
        })
    }
    pub fn source_file(&mut self, g: &Graph, source: usize) -> io::Result<File> {
        let lane = &g.lanes[source];
        if !self.verified.contains(&source) {
            lane.ports.verify(&self.root)?;
            self.verified.insert(source);
        }
        File::open(self.root.join(&lane.ports.path))
    }
    pub fn at_cut(
        &mut self,
        g: &Graph,
        source: usize,
        cut: u64,
        reverse: bool,
    ) -> io::Result<Option<Port>> {
        let shift = if reverse {
            g.k - g.cut_offset
        } else {
            g.cut_offset
        };
        let Some(anchor) = cut.checked_sub(shift) else {
            return Ok(None);
        };
        let mut file = self.source_file(g, source)?;
        let (mut lo, mut hi) = (0, g.lanes[source].port_count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let p = storage::port_at(&mut file, g.k as usize, mid)?;
            if (p.anchor, p.reverse) < (anchor, reverse) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        if lo == g.lanes[source].port_count {
            return Ok(None);
        }
        let p = storage::port_at(&mut file, g.k as usize, lo)?;
        Ok(((p.anchor, p.reverse) == (anchor, reverse)).then_some(p))
    }
}
