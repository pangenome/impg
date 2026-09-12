//! Explicit hypotheses, not routes inferred from ownership BEDs.
use super::*;
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub id: usize,
    pub name: String,
    pub length: u64,
}
/// One physical resource. Different IDs assert distinct copies, even if their
/// source intervals coincide. Repeated descriptions must use the same ID.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Instance {
    pub id: String,
    pub source: usize,
    pub start: u64,
    pub end: u64,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Piece {
    pub instance: String,
    pub start: u64,
    pub end: u64,
    pub strand: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Adjacency {
    pub from: usize,
    pub to: usize,
    pub kind: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Walk {
    pub id: String,
    pub topology: String,
    pub left_endpoint: String,
    pub right_endpoint: String,
    pub pieces: Vec<Piece>,
    pub adjacencies: Vec<Adjacency>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Slot {
    pub id: String,
    pub alternatives: Vec<Walk>,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Layout {
    pub version: u32,
    pub model: String,
    pub panel: PanelIdentity,
    pub sources: Vec<Source>,
    pub instances: Vec<Instance>,
    pub slots: Vec<Slot>,
}
fn name(s: &str) -> bool {
    !s.trim().is_empty() && !s.chars().any(char::is_control)
}
impl Layout {
    pub fn validate(&self) -> io::Result<()> {
        require(
            self.version == VERSION && self.model == LAYOUT_MODEL,
            "incompatible joint layout",
        )?;
        require(
            !self.sources.is_empty() && !self.slots.is_empty(),
            "empty joint layout",
        )?;
        let mut names = BTreeSet::new();
        for (id, s) in self.sources.iter().enumerate() {
            require(
                s.id == id
                    && name(&s.name)
                    && names.insert(&s.name)
                    && s.length > 0
                    && s.length <= i32::MAX as u64,
                "invalid source namespace",
            )?;
        }
        let mut instances = BTreeMap::new();
        for i in &self.instances {
            require(
                name(&i.id)
                    && instances.insert(&i.id, i).is_none()
                    && i.source < self.sources.len()
                    && i.start < i.end
                    && i.end <= self.sources[i.source].length,
                "duplicate/invalid physical instance",
            )?;
        }
        let mut slots = BTreeSet::new();
        for s in &self.slots {
            require(
                name(&s.id) && slots.insert(&s.id) && !s.alternatives.is_empty(),
                "invalid/duplicate slot",
            )?;
            let mut alternatives = BTreeSet::new();
            for w in &s.alternatives {
                require(
                    name(&w.id) && alternatives.insert(&w.id),
                    "invalid/duplicate walk ID",
                )?;
                require(
                    w.topology == "linear"
                        && w.left_endpoint == "asserted-molecule-terminus"
                        && w.right_endpoint == "asserted-molecule-terminus",
                    "unsupported topology or unknown endpoints",
                )?;
                require(
                    !w.pieces.is_empty() && w.adjacencies.len() == w.pieces.len() - 1,
                    "linear walk requires every adjacency explicitly",
                )?;
                for (n, a) in w.adjacencies.iter().enumerate() {
                    require(
                        a.from == n && a.to == n + 1 && a.kind == "abut",
                        "unsupported/undeclared adjacency",
                    )?;
                }
                let mut occupied: BTreeMap<&str, Vec<(u64, u64)>> = BTreeMap::new();
                for p in &w.pieces {
                    let i = instances
                        .get(&p.instance)
                        .ok_or_else(|| invalid("unknown physical instance"))?;
                    require(
                        p.start < p.end
                            && p.start >= i.start
                            && p.end <= i.end
                            && (p.strand == "+" || p.strand == "-"),
                        "invalid piece binding/strand",
                    )?;
                    let spans = occupied.entry(&p.instance).or_default();
                    require(
                        spans.iter().all(|&(a, b)| p.end <= a || p.start >= b),
                        "duplicate/overlapping description of one physical instance in walk",
                    )?;
                    spans.push((p.start, p.end));
                }
                w.length()?;
            }
        }
        Ok(())
    }
    pub fn resources(&self, slot: usize, alternative: usize) -> BTreeSet<String> {
        self.slots[slot].alternatives[alternative]
            .pieces
            .iter()
            .map(|p| p.instance.clone())
            .collect()
    }
    pub fn feasible(&self, choices: &[usize]) -> io::Result<bool> {
        require(
            choices.len() == self.slots.len(),
            "assignment must choose exactly one alternative per slot",
        )?;
        let mut used = BTreeSet::new();
        for (s, &a) in choices.iter().enumerate() {
            require(
                a < self.slots[s].alternatives.len(),
                "foreign assignment alternative",
            )?;
            for id in self.resources(s, a) {
                if !used.insert(id) {
                    return Ok(false);
                }
            }
        }
        Ok(true)
    }
}
impl Walk {
    pub fn length(&self) -> io::Result<u64> {
        self.pieces.iter().try_fold(0u64, |n, p| {
            n.checked_add(p.end - p.start)
                .ok_or_else(|| invalid("walk length overflow"))
        })
    }
}

pub struct Sources {
    indexes: Vec<UnifiedSequenceIndex>,
    chosen: Vec<usize>,
}
impl Sources {
    pub fn open(paths: &[String], layout: &Layout, panel: &SyngIndex) -> io::Result<Self> {
        let indexes = paths
            .iter()
            .map(|p| UnifiedSequenceIndex::from_files(&[p.clone()]))
            .collect::<io::Result<Vec<_>>>()?;
        let mut chosen = Vec::new();
        for source in &layout.sources {
            require(
                panel.name_map.path_to_name.get(source.id) == Some(&source.name)
                    && panel.name_map.path_to_length.get(source.id) == Some(&source.length),
                "layout source differs from exact panel namespace",
            )?;
            let matches: Vec<_> = indexes
                .iter()
                .enumerate()
                .filter_map(|(i, index)| {
                    index
                        .get_sequence_length(&source.name)
                        .ok()
                        .map(|n| (i, n as u64))
                })
                .collect();
            require(
                matches.len() == 1 && matches[0].1 == source.length,
                "missing/ambiguous source file binding",
            )?;
            chosen.push(matches[0].0);
        }
        Ok(Self { indexes, chosen })
    }
    /// Fetch only the requested walk crop. A crop may cross arbitrarily many joins.
    pub fn spell(&self, layout: &Layout, walk: &Walk, lo: u64, hi: u64) -> io::Result<Vec<u8>> {
        require(lo < hi && hi <= walk.length()?, "invalid walk crop")?;
        let mut dna = Vec::new();
        let mut offset = 0;
        for p in &walk.pieces {
            let end = offset + p.end - p.start;
            let a = lo.max(offset);
            let b = hi.min(end);
            if a < b {
                let instance = layout
                    .instances
                    .iter()
                    .find(|i| i.id == p.instance)
                    .unwrap();
                let (start, stop) = if p.strand == "+" {
                    (p.start + a - offset, p.start + b - offset)
                } else {
                    (p.end - (b - offset), p.end - (a - offset))
                };
                let mut part = self.indexes[self.chosen[instance.source]].fetch_sequence(
                    &layout.sources[instance.source].name,
                    start as i32,
                    stop as i32,
                )?;
                require(
                    part.len() as u64 == b - a
                        && part
                            .iter()
                            .all(|b| b"ACGTRYSWKMBDHVNacgtryswkmbdhvn".contains(b)),
                    "invalid source DNA/crop",
                )?;
                if p.strand == "-" {
                    part = crate::graph::reverse_complement(&part);
                }
                dna.extend(part);
            }
            offset = end;
            if offset >= hi {
                break;
            }
        }
        require(dna.len() as u64 == hi - lo, "incomplete walk spelling")?;
        Ok(dna)
    }
}
