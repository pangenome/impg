//! Stream legacy feature definitions while discarding large location/multiplicity
//! values. Only compact source/group/occurrence metadata is retained for ownership.
use super::super::catalog;
use super::*;
use serde::de::{DeserializeSeed, IgnoredAny, MapAccess, SeqAccess, Visitor};
use std::fmt;
#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
struct SlimFeature {
    tokens: Vec<u64>,
    owning_group: Option<usize>,
    locations: IgnoredAny,
    exclusion_reasons: IgnoredAny,
}
struct Features<'a>(&'a mut BufWriter<File>);
impl<'de> DeserializeSeed<'de> for Features<'_> {
    type Value = u64;
    fn deserialize<D: serde::Deserializer<'de>>(self, d: D) -> Result<u64, D::Error> {
        struct FeatureVisitor<'a>(&'a mut BufWriter<File>);
        impl<'de> Visitor<'de> for FeatureVisitor<'_> {
            type Value = u64;
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "complete legacy feature array")
            }
            fn visit_seq<A: SeqAccess<'de>>(self, mut a: A) -> Result<u64, A::Error> {
                let mut count = 0;
                while let Some(f) = a.next_element::<SlimFeature>()? {
                    let _ = (f.locations, f.exclusion_reasons);
                    crate::sample_mem_bwt::validate_pattern(&f.tokens)
                        .map_err(serde::de::Error::custom)?;
                    if f.tokens.len() != 3 || canonical(&f.tokens) != f.tokens {
                        return Err(serde::de::Error::custom("invalid registry pair"));
                    }
                    storage::line(self.0,&serde_json::json!({"id":count,"tokens":f.tokens,"original_owner":f.owning_group})).map_err(serde::de::Error::custom)?;
                    count += 1;
                }
                Ok(count)
            }
        }
        d.deserialize_seq(FeatureVisitor(self.0))
    }
}
struct Payload<'a>(&'a mut BufWriter<File>);
impl<'de> DeserializeSeed<'de> for Payload<'_> {
    type Value = (catalog::Catalog, u64);
    fn deserialize<D: serde::Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        struct PayloadVisitor<'a>(&'a mut BufWriter<File>);
        impl<'de> Visitor<'de> for PayloadVisitor<'_> {
            type Value = (catalog::Catalog, u64);
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "legacy catalog payload")
            }
            fn visit_map<A: MapAccess<'de>>(self, mut m: A) -> Result<Self::Value, A::Error> {
                let (
                    mut version,
                    mut panel,
                    mut accepted,
                    mut validation,
                    mut sources,
                    mut groups,
                    mut occurrences,
                    mut count,
                ) = (None, None, None, None, None, None, None, None);
                let mut seen = BTreeSet::new();
                while let Some(key) = m.next_key::<String>()? {
                    if !seen.insert(key.clone()) {
                        return Err(serde::de::Error::custom("duplicate catalog field"));
                    }
                    match key.as_str() {
                        "version" => version = Some(m.next_value()?),
                        "panel" => panel = Some(m.next_value()?),
                        "catalog_accepted" => accepted = Some(m.next_value()?),
                        "validation" => validation = Some(m.next_value()?),
                        "sources" => sources = Some(m.next_value()?),
                        "groups" => groups = Some(m.next_value()?),
                        "occurrences" => occurrences = Some(m.next_value::<Vec<SlimOccurrence>>()?),
                        "features" => count = Some(m.next_value_seed(Features(self.0))?),
                        "links" => {
                            m.next_value::<IgnoredAny>()?;
                        }
                        _ => return Err(serde::de::Error::custom("unknown catalog field")),
                    }
                }
                let missing = || serde::de::Error::custom("missing catalog field");
                Ok((
                    catalog::Catalog {
                        version: version.ok_or_else(missing)?,
                        panel: panel.ok_or_else(missing)?,
                        catalog_accepted: accepted.ok_or_else(missing)?,
                        validation: validation.ok_or_else(missing)?,
                        sources: sources.ok_or_else(missing)?,
                        groups: groups.ok_or_else(missing)?,
                        occurrences: occurrences
                            .ok_or_else(missing)?
                            .into_iter()
                            .map(|o| {
                                let _ = o.feature_multiplicities;
                                catalog::Occurrence {
                                    id: o.id,
                                    group: o.group,
                                    source: o.source,
                                    interval: o.interval,
                                    fully_contained_anchors: o.fully_contained_anchors,
                                    anchor_status: o.anchor_status,
                                    feature_multiplicities: BTreeMap::new(),
                                }
                            })
                            .collect(),
                        features: Vec::new(),
                        links: Vec::new(),
                    },
                    count.ok_or_else(missing)?,
                ))
            }
        }
        d.deserialize_map(PayloadVisitor(self.0))
    }
}
struct Envelope<'a>(&'a mut BufWriter<File>);
impl<'de> DeserializeSeed<'de> for Envelope<'_> {
    type Value = (catalog::Catalog, u64, String);
    fn deserialize<D: serde::Deserializer<'de>>(self, d: D) -> Result<Self::Value, D::Error> {
        struct EnvelopeVisitor<'a>(&'a mut BufWriter<File>);
        impl<'de> Visitor<'de> for EnvelopeVisitor<'_> {
            type Value = (catalog::Catalog, u64, String);
            fn expecting(&self, f: &mut fmt::Formatter) -> fmt::Result {
                write!(f, "legacy envelope")
            }
            fn visit_map<A: MapAccess<'de>>(self, mut m: A) -> Result<Self::Value, A::Error> {
                let (mut version, mut checksum, mut payload) = (None, None, None);
                let mut seen = BTreeSet::new();
                while let Some(key) = m.next_key::<String>()? {
                    if !seen.insert(key.clone()) {
                        return Err(serde::de::Error::custom("duplicate envelope field"));
                    }
                    match key.as_str() {
                        "version" => version = Some(m.next_value::<u32>()?),
                        "checksum" => checksum = Some(m.next_value::<String>()?),
                        "payload" => payload = Some(m.next_value_seed(Payload(self.0))?),
                        _ => return Err(serde::de::Error::custom("unknown envelope field")),
                    }
                }
                if version != Some(1) {
                    return Err(serde::de::Error::custom("invalid envelope version"));
                }
                let checksum =
                    checksum.ok_or_else(|| serde::de::Error::custom("missing checksum"))?;
                if checksum.len() != 16 || !checksum.bytes().all(|b| b.is_ascii_hexdigit()) {
                    return Err(serde::de::Error::custom("invalid old checksum header"));
                }
                let (c, n) = payload.ok_or_else(|| serde::de::Error::custom("missing payload"))?;
                Ok((c, n, checksum))
            }
        }
        d.deserialize_map(EnvelopeVisitor(self.0))
    }
}
pub fn load(
    path: &Path,
    identity: &PanelIdentity,
    out: &Path,
) -> io::Result<(catalog::Catalog, u64, serde_json::Value)> {
    let before = super::super::reconstruction::fingerprint(path)?;
    let mut writer = BufWriter::new(File::create(out.join("registry.jsonl"))?);
    let mut d = serde_json::Deserializer::from_reader(BufReader::new(File::open(path)?));
    let (catalog, count, header) = Envelope(&mut writer)
        .deserialize(&mut d)
        .map_err(io::Error::other)?;
    d.end().map_err(io::Error::other)?;
    writer.flush()?;
    require(
        catalog.version == 1 && catalog.panel == *identity && !catalog.catalog_accepted,
        "incompatible route ownership catalog",
    )?;
    catalog.validate()?;
    let mut groups = BTreeSet::new();
    let mut membership = BTreeSet::new();
    for group in &catalog.groups {
        require(
            !group.id.trim().is_empty()
                && !group.id.chars().any(char::is_control)
                && groups.insert(&group.id),
            "invalid ownership group namespace",
        )?;
        for &id in &group.occurrences {
            require(
                membership.insert(id),
                "duplicate original ownership membership",
            )?;
        }
    }
    require(
        membership.len() == catalog.occurrences.len(),
        "incomplete original ownership membership",
    )?;
    for o in &catalog.occurrences {
        require(
            o.interval
                .strand
                .as_deref()
                .is_none_or(|s| s == "+" || s == "-"),
            "unsupported ownership strand spelling",
        )?;
    }
    #[derive(Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Definition {
        id: u64,
        tokens: [u64; 3],
        original_owner: Option<usize>,
    }
    let mut definitions = BufReader::new(File::open(out.join("registry.jsonl"))?);
    let mut expected = 0;
    while let Some(d) = storage::next::<Definition>(&mut definitions)? {
        require(
            d.id == expected
                && d.original_owner.is_none_or(|g| g < catalog.groups.len())
                && canonical(&d.tokens) == d.tokens,
            "invalid original registry provenance",
        )?;
        add(&mut expected, 1)?;
    }
    require(expected == count, "incomplete registry definition import")?;
    require(
        before == super::super::reconstruction::fingerprint(path)?,
        "catalog changed during import",
    )?;
    Ok((
        catalog,
        count,
        serde_json::json!({"raw_input":before,"old_checksum_header_unverified":header,"full_feature_scope":true}),
    ))
}
