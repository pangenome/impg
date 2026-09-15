//! Finite explicit-layout joint evidence. This is not a whole-panel layout builder
//! or a scalable chromosome-mosaic solver; see docs/syng-gmem-bwt/joint-walks.md.
pub mod layout;
mod replay;
mod solver;
use super::{checksum, observations, read_json, sample, write_json, PanelIdentity, COUNT_POLICY};
use crate::{
    sample_mem_bwt::{canonical, encode_walk, invalid},
    syng::SyngIndex,
};
pub use replay::{Contribution, Profile};
use serde::{Deserialize, Serialize};
use serde_json::Value;
pub use solver::{evaluate, solve, Assignment, Evaluation, Factor, Problem, SearchResult};
use std::{
    collections::{BTreeMap, BTreeSet},
    io,
    path::Path,
};

pub const VERSION: u32 = 1;
pub const MODEL: &str = "joint-walk-poisson-v1";
pub const LAYOUT_MODEL: &str = "explicit-physical-linear-layout-v1";
const CORE_BP: u64 = 65536;
const MAX_LENGTH: u64 = 1048576;
fn require(ok: bool, message: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(message))
    }
}
pub fn compiler_identity() -> String {
    let mut h = checksum(observations::compiler_identity().as_bytes());
    for text in [
        include_str!("mod.rs"),
        include_str!("layout.rs"),
        include_str!("replay.rs"),
        include_str!("solver.rs"),
    ] {
        super::hash_update(&mut h, text.as_bytes());
    }
    format!("joint-fresh-native-v1-{h:016x}")
}
fn valid_tokens(t: &[u64; 3]) -> io::Result<()> {
    crate::sample_mem_bwt::validate_pattern(t)?;
    require(canonical(t) == *t, "noncanonical joint feature")
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Definition {
    pub tokens: [u64; 3],
    pub original_ids: Vec<usize>,
}
#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Compiled {
    pub version: u32,
    pub model: String,
    pub compiler_identity: String,
    pub count_policy: String,
    pub layout: layout::Layout,
    pub read_lengths: Vec<u64>,
    pub definitions: Vec<Definition>,
    pub profiles: Vec<Vec<Profile>>,
    pub source_files: Vec<Value>,
    pub registry: Option<Value>,
    pub context_complete: bool,
    pub candidate_universe: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Envelope {
    version: u32,
    model: String,
    checksum: String,
    payload: Compiled,
}
impl Compiled {
    pub fn validate(&self) -> io::Result<()> {
        require(
            self.version == VERSION
                && self.model == MODEL
                && self.count_policy == COUNT_POLICY
                && self.compiler_identity == compiler_identity()
                && self.context_complete
                && self.candidate_universe == "finite-explicit-layout-only",
            "incompatible/incomplete joint profiles",
        )?;
        self.layout.validate()?;
        require(
            !self.read_lengths.is_empty()
                && self.read_lengths.len() <= 32
                && self.read_lengths.iter().all(|&l| l > 0 && l <= MAX_LENGTH)
                && self.read_lengths.windows(2).all(|w| w[0] < w[1]),
            "invalid joint read lengths",
        )?;
        require(
            !self.source_files.is_empty()
                && self.source_files.iter().all(|v| {
                    v["bytes"].as_u64().is_some()
                        && v["fnv1a64"].as_str().is_some_and(|s| s.len() == 16)
                }),
            "missing source provenance",
        )?;
        let mut features = BTreeSet::new();
        let mut ids = BTreeSet::new();
        for d in &self.definitions {
            valid_tokens(&d.tokens)?;
            require(
                features.insert(d.tokens) && d.original_ids.iter().all(|id| ids.insert(*id)),
                "duplicate joint feature/registry ID",
            )?;
        }
        require(
            self.profiles.len() == self.layout.slots.len(),
            "incomplete slot profiles",
        )?;
        for (s, profiles) in self.profiles.iter().enumerate() {
            require(
                profiles.len() == self.layout.slots[s].alternatives.len(),
                "incomplete alternative profiles",
            )?;
            for (a, p) in profiles.iter().enumerate() {
                let m = self.layout.slots[s].alternatives[a].length()?;
                let expected: Vec<_> = self
                    .read_lengths
                    .iter()
                    .map(|&l| if l <= m { m - l + 1 } else { 0 })
                    .collect();
                require(
                    p.admitted_starts == expected
                        && p.max_crop_bp <= CORE_BP + MAX_LENGTH - 1
                        && (expected.iter().all(|&n| n == 0) || p.event_runs > 0),
                    "invalid/incomplete start-domain replay",
                )?;
                let mut seen = BTreeSet::new();
                for c in &p.contributions {
                    require(
                        features.contains(&c.tokens)
                            && seen.insert(c.tokens)
                            && c.totals.len() == self.read_lengths.len()
                            && c.totals
                                .iter()
                                .enumerate()
                                .all(|(i, q)| q[0] == q[1] && (expected[i] > 0 || q[0] == 0))
                            && c.totals.iter().any(|q| q[0] > 0),
                        "invalid sparse physical contribution",
                    )?;
                }
            }
        }
        Ok(())
    }
    pub fn digest(&self) -> io::Result<String> {
        Ok(format!(
            "{:016x}",
            checksum(&serde_json::to_vec(self).map_err(io::Error::other)?)
        ))
    }
    pub fn save(&self, path: &Path) -> io::Result<()> {
        self.validate()?;
        let digest = self.digest()?;
        #[derive(Serialize)]
        struct Borrowed<'a> {
            version: u32,
            model: &'a str,
            checksum: String,
            payload: &'a Compiled,
        }
        write_json(
            path,
            &Borrowed {
                version: VERSION,
                model: MODEL,
                checksum: digest,
                payload: self,
            },
        )
    }
    pub fn load(path: &Path, panel: &PanelIdentity) -> io::Result<(Self, String)> {
        let e: Envelope = read_json(path)?;
        require(
            e.version == VERSION
                && e.model == MODEL
                && e.payload.layout.panel == *panel
                && e.checksum
                    == format!(
                        "{:016x}",
                        checksum(&serde_json::to_vec(&e.payload).map_err(io::Error::other)?)
                    ),
            "joint artifact checksum/model/panel mismatch",
        )?;
        e.payload.validate()?;
        Ok((e.payload, e.checksum))
    }
}

pub fn compile(
    panel: &SyngIndex,
    identity: PanelIdentity,
    layout: layout::Layout,
    source_paths: &[String],
    mut read_lengths: Vec<u64>,
    registry_catalog: Option<&Path>,
) -> io::Result<Compiled> {
    layout.validate()?;
    require(layout.panel == identity, "layout/panel identity mismatch")?;
    read_lengths.sort_unstable();
    read_lengths.dedup();
    require(
        !read_lengths.is_empty()
            && read_lengths.len() <= 32
            && read_lengths.iter().all(|&l| l > 0 && l <= MAX_LENGTH),
        "invalid joint read lengths",
    )?;
    let sources = layout::Sources::open(source_paths, &layout, panel)?;
    let source_files = source_paths
        .iter()
        .map(|p| super::reconstruction::fingerprint(Path::new(p)))
        .collect::<io::Result<Vec<_>>>()?;
    let mut definitions: BTreeMap<[u64; 3], Vec<usize>> = BTreeMap::new();
    let registry = if let Some(path) = registry_catalog {
        // Full slim stream, no location matrix. Retain all definitions and original
        // IDs, never import old conditional exposures or spoof v2 compatibility.
        let (_, rows, provenance) = observations::input::load(path, &identity, &[])?;
        for d in rows {
            definitions.entry(d.tokens).or_default().push(d.id);
        }
        Some(provenance)
    } else {
        None
    };
    let mut profiles = Vec::new();
    for slot in &layout.slots {
        let mut alternatives = Vec::new();
        for walk in &slot.alternatives {
            let p =
                replay::compile_walk(panel, walk.length()?, &read_lengths, CORE_BP, |lo, hi| {
                    sources.spell(&layout, walk, lo, hi)
                })?;
            for c in &p.contributions {
                definitions.entry(c.tokens).or_default();
            }
            alternatives.push(p);
        }
        profiles.push(alternatives);
    }
    // Detect source mutation during replay rather than blessing mixed provenance.
    for (path, before) in source_paths.iter().zip(&source_files) {
        let after = super::reconstruction::fingerprint(Path::new(path))?;
        require(before == &after, "source changed during compilation")?;
    }
    let compiled = Compiled {
        version: VERSION,
        model: MODEL.into(),
        compiler_identity: compiler_identity(),
        count_policy: COUNT_POLICY.into(),
        layout,
        read_lengths,
        definitions: definitions
            .into_iter()
            .map(|(tokens, original_ids)| Definition {
                tokens,
                original_ids,
            })
            .collect(),
        profiles,
        source_files,
        registry,
        context_complete: true,
        candidate_universe: "finite-explicit-layout-only".into(),
    };
    compiled.validate()?;
    Ok(compiled)
}

#[cfg(test)]
mod tests;
