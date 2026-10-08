//! The legacy JSON envelope is traversed once; large location/multiplicity values
//! are consumed as IgnoredAny, never a serde_json::Value or legacy Catalog.
use super::*;
use serde::de::{IgnoredAny, SeqAccess, Visitor};
use std::fmt;

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Definition {
    pub id: usize,
    pub tokens: [u64; 3],
    pub original_owner: Option<usize>,
}
#[derive(Deserialize)]
struct SlimFeature {
    tokens: [u64; 3],
    owning_group: Option<usize>,
    locations: IgnoredAny,
    exclusion_reasons: IgnoredAny,
}
fn definitions<'de, D: serde::Deserializer<'de>>(d: D) -> Result<Vec<Definition>, D::Error> {
    struct Features;
    impl<'de> Visitor<'de> for Features {
        type Value = Vec<Definition>;
        fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
            write!(f, "legacy features")
        }
        fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<Self::Value, A::Error> {
            let mut out = Vec::new();
            while let Some(f) = a.next_element::<SlimFeature>()? {
                let _ = (f.locations, f.exclusion_reasons);
                out.push(Definition {
                    id: out.len(),
                    tokens: f.tokens,
                    original_owner: f.owning_group,
                });
            }
            Ok(out)
        }
    }
    d.deserialize_seq(Features)
}
#[derive(Deserialize)]
struct SlimOccurrence {
    id: usize,
    group: usize,
    source: usize,
    interval: catalog::Interval,
    fully_contained_anchors: usize,
    anchor_status: String,
    feature_multiplicities: IgnoredAny,
}
#[derive(Deserialize)]
struct Payload {
    version: u32,
    panel: PanelIdentity,
    catalog_accepted: bool,
    validation: String,
    sources: Vec<catalog::Source>,
    groups: Vec<catalog::Group>,
    occurrences: Vec<SlimOccurrence>,
    #[serde(deserialize_with = "definitions")]
    features: Vec<Definition>,
    links: IgnoredAny,
}
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct LegacyEnvelope {
    version: u32,
    checksum: String,
    payload: Payload,
}
struct HashReader<R> {
    reader: R,
    hash: u64,
    bytes: u64,
}
impl<R: Read> Read for HashReader<R> {
    fn read(&mut self, bytes: &mut [u8]) -> io::Result<usize> {
        let n = self.reader.read(bytes)?;
        super::super::hash_update(&mut self.hash, &bytes[..n]);
        self.bytes += n as u64;
        Ok(n)
    }
}

pub fn load(
    path: &Path,
    identity: &PanelIdentity,
    scope: &[String],
) -> io::Result<(catalog::Catalog, Vec<Definition>, Value)> {
    let mut reader = HashReader {
        reader: BufReader::with_capacity(4 * 1024 * 1024, File::open(path)?),
        hash: 0xcbf29ce484222325,
        bytes: 0,
    };
    let mut d = serde_json::Deserializer::from_reader(&mut reader);
    let envelope = LegacyEnvelope::deserialize(&mut d).map_err(io::Error::other)?;
    d.end().map_err(io::Error::other)?;
    let p = envelope.payload;
    require(
        envelope.version == FORMAT_VERSION
            && p.version == FORMAT_VERSION
            && !p.catalog_accepted
            && p.panel == *identity,
        "incompatible legacy catalog/panel",
    )?;
    require(
        envelope.checksum.len() == 16 && envelope.checksum.bytes().all(|b| b.is_ascii_hexdigit()),
        "invalid legacy checksum header",
    )?;
    let _ = (p.links, p.validation);
    let mut catalog = catalog::Catalog {
        version: FORMAT_VERSION,
        panel: p.panel,
        catalog_accepted: false,
        validation: "fixed physical ownership; slim observation metadata".into(),
        sources: p.sources,
        groups: p.groups,
        occurrences: Vec::new(),
        features: Vec::new(),
        links: Vec::new(),
    };
    for o in p.occurrences {
        let _ = o.feature_multiplicities;
        catalog.occurrences.push(catalog::Occurrence {
            id: o.id,
            group: o.group,
            source: o.source,
            interval: o.interval,
            fully_contained_anchors: o.fully_contained_anchors,
            anchor_status: o.anchor_status,
            feature_multiplicities: BTreeMap::new(),
        });
    }
    catalog.validate()?;
    let mut names = BTreeSet::new();
    for (id, source) in catalog.sources.iter().enumerate() {
        require(
            source.id == id && names.insert(&source.path),
            "duplicate/non-dense source namespace",
        )?;
        genotype::source_identity(&source.path)?;
    }
    let mut names = BTreeSet::new();
    let mut occurrences = BTreeSet::new();
    for group in &catalog.groups {
        require(
            !group.id.trim().is_empty()
                && !group.id.chars().any(char::is_control)
                && names.insert(&group.id),
            "invalid/duplicate group ID",
        )?;
        for &id in &group.occurrences {
            require(occurrences.insert(id), "duplicate original occurrence ID")?;
        }
    }
    require(
        occurrences.len() == catalog.occurrences.len(),
        "incomplete original group membership",
    )?;
    for o in &catalog.occurrences {
        require(
            o.interval
                .strand
                .as_deref()
                .is_none_or(|s| s == "+" || s == "-"),
            "unknown source strand spelling",
        )?;
    }
    let wanted: BTreeSet<_> = scope
        .iter()
        .map(|name| {
            catalog
                .groups
                .iter()
                .position(|g| &g.id == name)
                .ok_or_else(|| invalid(format!("unknown feature group {name}")))
        })
        .collect::<io::Result<_>>()?;
    let total = p.features.len();
    let mut features = p.features;
    for f in &features {
        require(
            f.original_owner.is_none_or(|g| g < catalog.groups.len()),
            "foreign original feature owner",
        )?;
        let t = f.tokens;
        require(
            t[0] >= 4
                && t[0] <= 4 * i32::MAX as u64 + 2
                && t[0] % 2 == 0
                && t[2] >= 4
                && t[2] <= 4 * i32::MAX as u64 + 2
                && t[2] % 2 == 0
                && t[1] >= 3
                && t[1] % 2 == 1
                && t[1] / 2 <= u32::MAX as u64,
            "invalid two-node feature tokens",
        )?;
        require(canonical(&t) == t, "noncanonical feature definition")?;
    }
    if !scope.is_empty() {
        features.retain(|f| f.original_owner.is_some_and(|g| wanted.contains(&g)));
    }
    require(!features.is_empty(), "empty selected feature universe")?;
    let provenance = json!({"path":path,"raw_file_bytes":reader.bytes,"raw_file_fnv1a64":format!("{:016x}",reader.hash),
        "legacy_payload_checksum_header_unverified":envelope.checksum,"original_feature_count":total,
        "selected_feature_count":features.len(),"feature_groups":scope,"full_feature_scope":scope.is_empty()});
    Ok((catalog, features, provenance))
}

/// Normalize same-owner descriptions into disjoint maximal pieces. Preserve all
/// original occurrence IDs separately. Different-owner overlaps and uncovered
/// source bases are fatal, including paths with no requested features.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Core {
    pub source: usize,
    pub group: usize,
    pub start: u64,
    pub end: u64,
}
pub fn ownership(catalog: &catalog::Catalog) -> io::Result<Vec<Vec<Core>>> {
    let mut rows = vec![Vec::new(); catalog.sources.len()];
    for o in &catalog.occurrences {
        rows[o.source].push(o);
    }
    let mut result = Vec::new();
    for (source, mut rows) in rows.into_iter().enumerate() {
        rows.sort_by_key(|o| (o.interval.start, o.interval.end, o.id));
        let mut cores: Vec<Core> = Vec::new();
        for o in rows {
            let i = &o.interval;
            if let Some(last) = cores.last_mut() {
                require(i.start <= last.end, "gap in fixed partition ownership")?;
                if i.start < last.end {
                    require(last.group == o.group, "conflicting partition ownership")?;
                }
                if last.group == o.group {
                    last.end = last.end.max(i.end);
                    continue;
                }
            } else {
                require(i.start == 0, "unowned source prefix")?;
            }
            cores.push(Core {
                source,
                group: o.group,
                start: i.start,
                end: i.end,
            });
        }
        require(
            cores
                .last()
                .is_some_and(|c| c.end == catalog.sources[source].length),
            "unowned source/suffix",
        )?;
        result.push(cores);
    }
    Ok(result)
}

pub struct Sources {
    indexes: Vec<UnifiedSequenceIndex>,
    chosen: Vec<usize>,
}
impl Sources {
    pub fn open(paths: &[String], catalog: &catalog::Catalog) -> io::Result<Self> {
        require(!paths.is_empty(), "source files required")?;
        let indexes = paths
            .iter()
            .map(|s| UnifiedSequenceIndex::from_files(&[s.clone()]))
            .collect::<io::Result<Vec<_>>>()?;
        let mut chosen = Vec::new();
        for source in &catalog.sources {
            require(
                source.length <= i32::MAX as u64,
                "source exceeds native coordinate cap",
            )?;
            let found: Vec<_> = indexes
                .iter()
                .enumerate()
                .filter_map(|(id, index)| {
                    index
                        .get_sequence_length(&source.path)
                        .ok()
                        .map(|l| (id, l))
                })
                .collect();
            require(
                found.len() == 1 && found[0].1 as u64 == source.length,
                "missing/duplicate source binding or length mismatch",
            )?;
            chosen.push(found[0].0);
        }
        Ok(Self { indexes, chosen })
    }
    pub fn fetch(
        &self,
        catalog: &catalog::Catalog,
        source: usize,
        start: u64,
        end: u64,
    ) -> io::Result<Vec<u8>> {
        require(
            start < end && end <= catalog.sources[source].length,
            "invalid source fetch",
        )?;
        let dna = self.indexes[self.chosen[source]].fetch_sequence(
            &catalog.sources[source].path,
            start as i32,
            end as i32,
        )?;
        require(
            dna.len() as u64 == end - start
                && dna
                    .iter()
                    .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b)),
            "invalid source DNA/crop length",
        )?;
        Ok(dna)
    }
}
