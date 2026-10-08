use super::{PanelIdentity, FORMAT_VERSION};
use crate::sample_mem_bwt::{canonical, encode_walk, invalid};
use crate::syng::SyngIndex;
use serde::{Deserialize, Serialize};
use std::collections::{BTreeMap, BTreeSet};
use std::fs;
use std::io;
use std::path::Path;

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scaffold {
    pub component: String,
    pub start: u64,
    pub end: u64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Interval {
    pub path: String,
    pub start: u64,
    pub end: u64,
    /// None for BED3: ownership intervals do not establish homolog orientation.
    pub strand: Option<String>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct GroupInput {
    pub id: String,
    pub scaffold: Option<Scaffold>,
    pub occurrences: Vec<Interval>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CatalogInput {
    pub version: u32,
    pub groups: Vec<GroupInput>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub id: usize,
    pub path: String,
    pub length: u64,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Group {
    pub id: String,
    pub scaffold: Option<Scaffold>,
    pub occurrences: Vec<usize>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Occurrence {
    pub id: usize,
    pub group: usize,
    pub source: usize,
    pub interval: Interval,
    pub fully_contained_anchors: usize,
    pub anchor_status: String,
    pub feature_multiplicities: BTreeMap<usize, u64>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FeatureLocation {
    pub source: usize,
    pub start: u64,
    pub end: u64,
    pub signed_nodes: [i32; 2],
    pub containing_occurrences: Vec<usize>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Feature {
    pub tokens: Vec<u64>,
    pub locations: Vec<FeatureLocation>,
    pub owning_group: Option<usize>,
    pub exclusion_reasons: Vec<String>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Link {
    pub from: usize,
    pub to: usize,
    pub source_gap_bp: i64,
    pub relation: String,
    pub oriented_continuation: bool,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Catalog {
    pub version: u32,
    pub panel: PanelIdentity,
    pub catalog_accepted: bool,
    pub validation: String,
    pub sources: Vec<Source>,
    pub groups: Vec<Group>,
    pub occurrences: Vec<Occurrence>,
    pub features: Vec<Feature>,
    pub links: Vec<Link>,
}

pub fn import_beds(directory: &Path) -> io::Result<CatalogInput> {
    let mut files = fs::read_dir(directory)?
        .map(|r| r.map(|e| e.path()))
        .collect::<io::Result<Vec<_>>>()?;
    files.retain(|p| p.extension().is_some_and(|e| e == "bed"));
    files.sort();
    let mut groups = Vec::new();
    for path in files {
        let id = path
            .file_stem()
            .and_then(|s| s.to_str())
            .ok_or_else(|| invalid("non-UTF8 BED filename"))?
            .to_string();
        let mut occurrences = Vec::new();
        for (i, line) in fs::read_to_string(&path)?.lines().enumerate() {
            if line.is_empty() {
                continue;
            }
            let fields: Vec<_> = line.split('\t').collect();
            if fields.len() != 3 {
                return Err(invalid(format!(
                    "{}:{}: expected BED3",
                    path.display(),
                    i + 1
                )));
            }
            occurrences.push(Interval {
                path: fields[0].into(),
                start: fields[1].parse().map_err(io::Error::other)?,
                end: fields[2].parse().map_err(io::Error::other)?,
                strand: None,
            });
        }
        groups.push(GroupInput {
            id,
            scaffold: None,
            occurrences,
        });
    }
    Ok(CatalogInput {
        version: FORMAT_VERSION,
        groups,
    })
}
fn valid_name(name: &str) -> bool {
    !name.trim().is_empty() && !name.chars().any(char::is_control)
}

pub fn build(
    panel: &SyngIndex,
    identity: PanelIdentity,
    input: CatalogInput,
) -> io::Result<Catalog> {
    if input.version != FORMAT_VERSION || input.groups.is_empty() {
        return Err(invalid("empty/incompatible catalog input"));
    }
    let mut catalog = Catalog { version: FORMAT_VERSION, panel: identity, catalog_accepted: false,
        validation: "source-bounds-and-stored-anchors-only;homology-spanning-and-copy-structure-unresolved;no-source-sequence-validation".into(),
        sources: panel.name_map.path_to_name.iter().enumerate().map(|(id, path)| Source {
            id, path: path.clone(), length: panel.name_map.path_to_length[id] }).collect(),
        groups: Vec::new(), occurrences: Vec::new(), features: Vec::new(), links: Vec::new() };
    let mut names = BTreeSet::new();
    for source in &catalog.sources {
        if !valid_name(&source.path)
            || !names.insert(&source.path)
            || source.length > i64::MAX as u64
        {
            return Err(invalid(
                "invalid/duplicate panel source name or oversized source",
            ));
        }
    }
    let mut group_names = BTreeSet::new();
    let mut by_source: Vec<Vec<usize>> = vec![Vec::new(); catalog.sources.len()];
    let mut feature_ids = BTreeMap::new();
    for group in input.groups {
        if !valid_name(&group.id)
            || !group_names.insert(group.id.clone())
            || group.occurrences.is_empty()
        {
            return Err(invalid(
                "invalid/duplicate group ID or empty ownership group",
            ));
        }
        if let Some(s) = &group.scaffold {
            if !valid_name(&s.component) || s.start >= s.end {
                return Err(invalid("invalid scaffold interval"));
            }
        }
        let group_id = catalog.groups.len();
        let mut ids = Vec::new();
        for interval in group.occurrences {
            let source = *panel
                .name_map
                .name_to_path
                .get(&interval.path)
                .ok_or_else(|| invalid(format!("unknown full source path: {}", interval.path)))?
                as usize;
            if interval.start >= interval.end
                || interval.end > catalog.sources[source].length
                || interval
                    .strand
                    .as_deref()
                    .is_some_and(|s| s != "+" && s != "-")
            {
                return Err(invalid("invalid occurrence bounds or orientation"));
            }
            let anchors: Vec<_> = panel
                .walk_path_range(source, interval.start, interval.end)?
                .into_iter()
                .filter(|&(_, p)| {
                    p >= interval.start && p + panel.syncmer_length_bp() as u64 <= interval.end
                })
                .collect();
            let mut multiplicities = BTreeMap::new();
            for pair in anchors.windows(2) {
                let tokens = canonical(&encode_walk(pair)?);
                let next = feature_ids.len();
                let f = *feature_ids.entry(tokens).or_insert(next);
                *multiplicities.entry(f).or_default() += 1;
            }
            let id = catalog.occurrences.len();
            ids.push(id);
            by_source[source].push(id);
            catalog.occurrences.push(Occurrence {
                id,
                group: group_id,
                source,
                interval,
                fully_contained_anchors: anchors.len(),
                anchor_status: if anchors.len() >= 2 {
                    "stored-context-present"
                } else {
                    "unsupported-no-two-node-context"
                }
                .into(),
                feature_multiplicities: multiplicities,
            });
        }
        catalog.groups.push(Group {
            id: group.id,
            scaffold: group.scaffold,
            occurrences: ids,
        });
    }
    // IDs reflect deterministic source occurrence input order, not tree key order.
    let mut ordered: Vec<_> = feature_ids.iter().collect();
    ordered.sort_by_key(|(_, id)| **id);
    catalog.features = ordered
        .into_iter()
        .map(|(tokens, _)| Feature {
            tokens: tokens.clone(),
            locations: Vec::new(),
            owning_group: None,
            exclusion_reasons: Vec::new(),
        })
        .collect();
    let k = panel.syncmer_length_bp() as u64;
    for (source, ids) in by_source.iter_mut().enumerate() {
        ids.sort_by_key(|&id| {
            let i = &catalog.occurrences[id].interval;
            (i.start, i.end, id)
        });
        // Scan complete stored paths, not cropped candidate walks. A feature's
        // complete bp span must fit an owning occurrence; crossing ranges survive
        // as explicit excluded global occurrences even if their starts lie inside.
        let walk = panel.walk_path_range(source, 0, catalog.sources[source].length)?;
        let mut cursor = 0;
        let mut active: Vec<usize> = Vec::new();
        for pair in walk.windows(2) {
            let start = pair[0].1;
            let end = pair[1].1 + k;
            let Some(&f) = feature_ids.get(&canonical(&encode_walk(pair)?)) else {
                continue;
            };
            // Include all intervals overlapping the span, including intervals
            // beginning between the two anchors. Starts progress monotonically.
            while cursor < ids.len() && catalog.occurrences[ids[cursor]].interval.start < end {
                active.push(ids[cursor]);
                cursor += 1;
            }
            active.retain(|&id| catalog.occurrences[id].interval.end > start);
            let overlaps: Vec<_> = active
                .iter()
                .copied()
                .filter(|&id| catalog.occurrences[id].interval.start < end)
                .collect();
            let containing: Vec<_> = overlaps
                .iter()
                .copied()
                .filter(|&id| {
                    let i = &catalog.occurrences[id].interval;
                    i.start <= start && end <= i.end
                })
                .collect();
            let feature = &mut catalog.features[f];
            if overlaps.is_empty() {
                feature.exclusion_reasons.push("outside-catalog".into());
            }
            if overlaps.len() != containing.len() {
                feature.exclusion_reasons.push("boundary-crossing".into());
            }
            feature.locations.push(FeatureLocation {
                source,
                start,
                end,
                signed_nodes: [pair[0].0, pair[1].0],
                containing_occurrences: containing,
            });
        }
        // Source-neighbor relations only, never joins between different paths.
        // Duplicate intervals retain ALL neighbor links, not just the arbitrarily
        // last member of an equal-coordinate bucket. No biological joins inferred.
        let mut buckets: Vec<Vec<usize>> = Vec::new();
        for &id in ids.iter() {
            let interval = &catalog.occurrences[id].interval;
            if buckets.last().is_some_and(|b| {
                let previous = &catalog.occurrences[b[0]].interval;
                (previous.start, previous.end) == (interval.start, interval.end)
            }) {
                buckets.last_mut().unwrap().push(id);
            } else {
                buckets.push(vec![id]);
            }
        }
        for pair in buckets.windows(2) {
            for &a_id in &pair[0] {
                for &b_id in &pair[1] {
                    let (a, b) = (&catalog.occurrences[a_id], &catalog.occurrences[b_id]);
                    let gap = b.interval.start as i64 - a.interval.end as i64;
                    let reverse = a.interval.strand.as_deref() == Some("-")
                        && b.interval.strand.as_deref() == Some("-");
                    let forward = a.interval.strand.as_deref() == Some("+")
                        && b.interval.strand.as_deref() == Some("+");
                    catalog.links.push(Link {
                        from: if reverse { b.id } else { a.id },
                        to: if reverse { a.id } else { b.id },
                        source_gap_bp: gap,
                        relation: if gap < 0 {
                            "overlap"
                        } else if gap > 0 {
                            "gap"
                        } else {
                            "adjacent"
                        }
                        .into(),
                        oriented_continuation: (reverse || forward)
                            && b.interval.start > a.interval.start
                            && b.interval.end > a.interval.end,
                    });
                }
            }
        }
    }
    for feature in &mut catalog.features {
        let owners: BTreeSet<_> = feature
            .locations
            .iter()
            .flat_map(|l| {
                l.containing_occurrences
                    .iter()
                    .map(|&id| catalog.occurrences[id].group)
            })
            .collect();
        if owners.len() > 1 {
            feature
                .exclusion_reasons
                .push("shared-between-groups".into());
        }
        if feature.locations.is_empty() {
            return Err(invalid(
                "candidate feature absent from full stored-path scan",
            ));
        }
        feature.exclusion_reasons.sort();
        feature.exclusion_reasons.dedup();
        if feature.exclusion_reasons.is_empty() && owners.len() == 1 {
            feature.owning_group = owners.first().copied();
        }
    }
    catalog.validate()?;
    Ok(catalog)
}
impl Catalog {
    pub fn validate(&self) -> io::Result<()> {
        if self.version != FORMAT_VERSION || self.catalog_accepted || self.groups.is_empty() {
            return Err(invalid(
                "incompatible/empty catalog or unsupported certification claim",
            ));
        }
        for (id, o) in self.occurrences.iter().enumerate() {
            if o.id != id
                || o.source >= self.sources.len()
                || o.group >= self.groups.len()
                || o.interval.path != self.sources[o.source].path
                || o.interval.start >= o.interval.end
                || o.interval.end > self.sources[o.source].length
                || o.feature_multiplicities
                    .keys()
                    .any(|&f| f >= self.features.len())
            {
                return Err(invalid("malformed catalog occurrence"));
            }
        }
        for (g, group) in self.groups.iter().enumerate() {
            if group.occurrences.is_empty()
                || group
                    .occurrences
                    .iter()
                    .any(|&id| id >= self.occurrences.len() || self.occurrences[id].group != g)
            {
                return Err(invalid("malformed catalog group"));
            }
        }
        for feature in &self.features {
            if feature.tokens.len() != 3
                || feature.owning_group.is_some_and(|g| g >= self.groups.len())
                || feature.locations.iter().any(|l| {
                    l.source >= self.sources.len()
                        || l.start >= l.end
                        || l.end > self.sources[l.source].length
                        || l.containing_occurrences
                            .iter()
                            .any(|&id| id >= self.occurrences.len())
                })
            {
                return Err(invalid("malformed catalog feature"));
            }
        }
        Ok(())
    }
}
