//! Partition-major observation-compatible, occurrence-weighted error-free replay.
//! Source ownership is fixed. Raw sketch views are internal observers, not copies.
mod consumer;
mod orientation;
pub use orientation::{prepare_orientations, write_derived_calls};
pub(crate) mod input;
pub(crate) mod profile;
use super::{catalog, genotype, sample, threading, PanelIdentity, COUNT_POLICY, FORMAT_VERSION};
use crate::{
    sample_mem_bwt::{canonical, encode_walk, invalid},
    sequence_index::{SequenceIndex, UnifiedSequenceIndex},
    syng::SyngIndex,
};
pub use consumer::{call, classify, eligible, thread, Provenance};
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet, BinaryHeap},
    fs::{self, File},
    io::{self, BufRead, BufReader, BufWriter, Read, Write},
    path::{Path, PathBuf},
};

pub const MODEL: &str = "partition-observation-poisson-v2";
pub const THREAD_MODEL: &str = "partition-observation-source-min-sum-v2";
pub const CHECKSUM_ALGORITHM: &str =
    "fnv1a64-content-v1; sample=bincode-payload; observations=metadata-content";
const CORE_BP: u64 = 1024 * 1024;
const MERGE_FAN_IN: usize = 32;
const MAX_LENGTH: u64 = 1024 * 1024;
pub(super) fn require(ok: bool, message: &str) -> io::Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(message))
    }
}

#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Incidence {
    pub feature: usize,
    pub source: usize,
    pub start: u64,
    pub end: u64,
    /// 1 canonical signed pair, 2 reverse pair, 3 conflicting; 0 palindrome.
    pub signed_mask: u8,
    pub group: Option<usize>,
    pub context_nonlocal: bool,
    pub source_terminal_context: bool,
    /// Exact surviving substrings by supported length, then input orientation.
    pub contributions: Vec<[u64; 2]>,
}
impl Incidence {
    fn key(&self) -> (usize, usize, u64, u64) {
        (self.feature, self.source, self.start, self.end)
    }
}
#[derive(Serialize, Deserialize)]
pub struct Metadata {
    pub version: u32,
    pub model: String,
    pub count_policy: String,
    pub compiler_identity: String,
    pub capability: String,
    pub catalog: catalog::Catalog,
    pub legacy_input: Value,
    pub sources: Vec<Value>,
    pub syncmer_length_bp: u64,
    pub read_lengths: Vec<u64>,
    pub feature_groups: Vec<String>,
    pub full_feature_scope: bool,
    pub feature_count: usize,
    pub files: BTreeMap<String, Value>,
    pub stats: Value,
}

fn json_line(writer: &mut impl Write, value: &impl Serialize) -> io::Result<()> {
    serde_json::to_writer(&mut *writer, value).map_err(io::Error::other)?;
    writer.write_all(b"\n")
}
fn finish(writer: &mut BufWriter<File>, temporary: &Path, path: &Path) -> io::Result<()> {
    writer.flush()?;
    writer.get_ref().sync_all()?;
    fs::rename(temporary, path)
}
fn write_lines<T: Serialize>(path: &Path, rows: impl IntoIterator<Item = T>) -> io::Result<()> {
    let temporary = path.with_extension("incomplete");
    let mut writer = BufWriter::new(File::create(&temporary)?);
    for row in rows {
        json_line(&mut writer, &row)?;
    }
    finish(&mut writer, &temporary, path)
}
pub(super) fn next<T: serde::de::DeserializeOwned>(
    reader: &mut impl BufRead,
) -> io::Result<Option<T>> {
    let mut line = String::new();
    // One incidence/definition has fixed metadata plus <=32 length entries.
    let n = reader.take(65537).read_line(&mut line)?;
    require(n <= 65536, "observation row exceeds 64KiB cap")?;
    if n == 0 {
        Ok(None)
    } else {
        serde_json::from_str(&line)
            .map(Some)
            .map_err(io::Error::other)
    }
}
fn fingerprint(path: &Path) -> io::Result<Value> {
    super::reconstruction::fingerprint(path)
}
fn verify(path: &Path, expected: &Value) -> io::Result<()> {
    let actual = fingerprint(path)?;
    require(
        actual["bytes"] == expected["bytes"] && actual["fnv1a64"] == expected["fnv1a64"],
        "observation content fingerprint mismatch",
    )
}
pub(crate) fn compiler_identity() -> String {
    // Source-bound algorithm identity, independent of build path or timestamps.
    let mut hash = 0xcbf29ce484222325;
    for source in [
        include_str!("mod.rs"),
        include_str!("input.rs"),
        include_str!("profile.rs"),
        include_str!("consumer.rs"),
        include_str!("orientation.rs"),
        include_str!("../sample.rs"),
        include_str!("../../syng.rs"),
        include_str!("../../sample_mem_bwt.rs"),
        include_str!("../../../vendor/syng/seqhash.c"),
        include_str!("../../../vendor/syng/kmerhash.c"),
        include_str!("../../../vendor/syng/syngbwt3.c"),
        include_str!("../../../vendor/syng/impg_syng_helpers.c"),
    ] {
        super::hash_update(&mut hash, source.as_bytes());
    }
    format!("partition-replay-v2-{hash:016x}")
}
#[derive(Clone)]
struct Endpoint {
    source: usize,
    start: u64,
    view: usize,
    first: (i32, u64),
    last: (i32, u64),
}

fn candidate(
    pair: &[(i32, u64)],
    source: usize,
    k: u64,
    definitions: &[input::Definition],
    lookup: &[usize],
) -> io::Result<Option<(usize, Incidence)>> {
    let tokens = encode_walk(pair)?;
    let canonical_tokens = canonical(&tokens);
    let Ok(p) =
        lookup.binary_search_by(|&i| definitions[i].tokens.as_slice().cmp(&canonical_tokens))
    else {
        return Ok(None);
    };
    let f = lookup[p];
    let palindrome = tokens == crate::sample_mem_bwt::reverse_complement(&tokens);
    Ok(Some((
        f,
        Incidence {
            feature: definitions[f].id,
            source,
            start: pair[0].1,
            end: pair[1].1 + k,
            signed_mask: if palindrome {
                0
            } else if tokens == canonical_tokens {
                1
            } else {
                2
            },
            group: None,
            context_nonlocal: false,
            source_terminal_context: false,
            contributions: Vec::new(),
        },
    )))
}
fn shard(
    path: &Path,
    mut hits: Vec<(usize, Incidence)>,
    panel: &SyngIndex,
    sources: &input::Sources,
    catalog: &catalog::Catalog,
    cores: &[Vec<input::Core>],
    lengths: &[u64],
    definitions: &[input::Definition],
    events: &mut u64,
    context: Option<&profile::SourceContext>,
    profile_stats: &mut profile::CompilationStats,
) -> io::Result<usize> {
    hits.sort_by_key(|(_, h)| h.key());
    let mut unique: Vec<(usize, Incidence)> = Vec::new();
    for (f, hit) in hits {
        if let Some((_, last)) = unique
            .last_mut()
            .filter(|(_, last)| last.key() == hit.key())
        {
            last.signed_mask |= hit.signed_mask;
        } else {
            unique.push((f, hit));
        }
    }
    for (f, hit) in &mut unique {
        profile::add(
            events,
            profile::compile_with_context(
                panel,
                sources,
                catalog,
                cores,
                lengths,
                &definitions[*f],
                hit,
                context,
                profile_stats,
            )?,
        )?;
    }
    let n = unique.len();
    write_lines(path, unique.into_iter().map(|(_, h)| h))?;
    Ok(n)
}

/// Merge only bounded fan-in sorted streams. Duplicate raw descriptions cannot
/// add contributions: their source-context profiles must agree exactly.
fn merge(paths: &[PathBuf], out: &Path) -> io::Result<()> {
    let mut readers = paths
        .iter()
        .map(|p| File::open(p).map(BufReader::new))
        .collect::<io::Result<Vec<_>>>()?;
    let mut heads = Vec::new();
    let mut heap = BinaryHeap::new();
    for (i, r) in readers.iter_mut().enumerate() {
        let head: Option<Incidence> = next(r)?;
        if let Some(h) = &head {
            heap.push(std::cmp::Reverse((h.key(), i)));
        }
        heads.push(head);
    }
    let tmp = out.with_extension("incomplete");
    let mut writer = BufWriter::new(File::create(&tmp)?);
    let mut previous: Option<Incidence> = None;
    while let Some(std::cmp::Reverse((key, i))) = heap.pop() {
        let h = heads[i].take().unwrap();
        if let Some(old) = previous.as_mut().filter(|old| old.key() == key) {
            require(
                old.contributions == h.contributions
                    && old.group == h.group
                    && old.context_nonlocal == h.context_nonlocal
                    && old.source_terminal_context == h.source_terminal_context,
                "duplicate physical profile disagreement",
            )?;
            old.signed_mask |= h.signed_mask;
        } else {
            if let Some(old) = previous.replace(h) {
                json_line(&mut writer, &old)?;
            }
        }
        heads[i] = next(&mut readers[i])?;
        if let Some(h) = &heads[i] {
            require(h.key() > key, "unsorted or duplicate shard")?;
            heap.push(std::cmp::Reverse((h.key(), i)));
        }
    }
    if let Some(h) = previous {
        json_line(&mut writer, &h)?;
    }
    finish(&mut writer, &tmp, out)
}

/// Run inside with_output_model: pre-existing directories are always refused;
/// failed/running manifests preserve diagnostic shards, never reusable success.
pub fn build(
    panel: &SyngIndex,
    identity: PanelIdentity,
    legacy: &Path,
    source_paths: &[String],
    mut lengths: Vec<u64>,
    scope: &[String],
    out: &Path,
) -> io::Result<Value> {
    lengths.sort_unstable();
    lengths.dedup();
    require(
        !lengths.is_empty()
            && lengths.len() <= 32
            && lengths.iter().all(|&l| l > 0 && l <= MAX_LENGTH),
        "supported lengths require 1..=32 values in 1..=1048576",
    )?;
    let (catalog, mut definitions, legacy_input) = input::load(legacy, &identity, scope)?;
    definitions.shrink_to_fit();
    for f in &definitions {
        for t in [f.tokens[0], f.tokens[2]] {
            let z = (t - 2) / 2;
            let n = ((z >> 1) as i64) ^ -((z & 1) as i64);
            require(
                n != 0 && n.unsigned_abs() <= panel.num_syncmer_nodes() as u64,
                "feature contains foreign dictionary node",
            )?;
        }
    }
    require(
        catalog.sources.len() == panel.name_map.path_to_name.len()
            && catalog.sources.iter().enumerate().all(|(i, s)| {
                s.id == i
                    && s.path == panel.name_map.path_to_name[i]
                    && s.length == panel.name_map.path_to_length[i]
            }),
        "catalog/panel source namespace mismatch",
    )?;
    let cores = input::ownership(&catalog)?;
    let sources = input::Sources::open(source_paths, &catalog)?;
    let source_fingerprints = source_paths
        .iter()
        .map(|p| fingerprint(Path::new(p)))
        .collect::<io::Result<Vec<_>>>()?;
    let mut lookup: Vec<_> = (0..definitions.len()).collect();
    lookup.sort_unstable_by_key(|&i| definitions[i].tokens);
    require(
        lookup
            .windows(2)
            .all(|p| definitions[p[0]].tokens != definitions[p[1]].tokens),
        "duplicate feature tokens",
    )?;
    write_lines(&out.join("definitions.jsonl"), &definitions)?;
    super::write_json(&out.join("ownership.json"), &cores)?;
    let mut grouped = vec![Vec::new(); catalog.groups.len()];
    for core in cores.iter().flatten() {
        grouped[core.group].push(core);
    }
    #[derive(Serialize)]
    struct InputBinding<'a> {
        model: &'a str,
        compiler_identity: String,
        catalog: &'a catalog::Catalog,
        legacy_input: &'a Value,
        source_files: &'a [Value],
        read_lengths: &'a [u64],
        feature_groups: &'a [String],
        definitions: Value,
        ownership: Value,
    }
    super::write_json(
        &out.join("inputs.json"),
        &InputBinding {
            model: MODEL,
            compiler_identity: compiler_identity(),
            catalog: &catalog,
            legacy_input: &legacy_input,
            source_files: &source_fingerprints,
            read_lengths: &lengths,
            feature_groups: scope,
            definitions: fingerprint(&out.join("definitions.jsonl"))?,
            ownership: fingerprint(&out.join("ownership.json"))?,
        },
    )?;
    fs::create_dir(out.join("shards"))?;
    let mut shards = Vec::new();
    let mut endpoints = Vec::new();
    let mut events = 0;
    let mut profile_stats = profile::CompilationStats::default();
    let mut context_bp_fetched = 0u64;
    let mut physical = 0usize;
    let k = panel.syncmer_length_bp() as u64;
    let mut core_count = 0;
    let halo = (*lengths.last().unwrap()).max(k - 1);
    for (g, group) in grouped.into_iter().enumerate() {
        let mut checkpoint = Vec::new();
        for core in group {
            let mut start = core.start;
            while start < core.end {
                let end = (start + CORE_BP).min(core.end);
                let crop_start = start.saturating_sub(halo);
                let crop_end = end
                    .checked_add(halo)
                    .ok_or_else(|| invalid("core context overflow"))?
                    .min(catalog.sources[core.source].length);
                let dna = sources.fetch(&catalog, core.source, crop_start, crop_end)?;
                profile::add(&mut context_bp_fetched, crop_end - crop_start)?;
                // Validate every stored anchor's actual signed DNA, without
                // assuming that stored forward walks enumerate both raw views.
                for (node, p) in panel.walk_path_range(core.source, crop_start, crop_end)? {
                    if p >= crop_start && p + k <= crop_end {
                        require(
                            panel.syncmer_seq(node).eq_ignore_ascii_case(
                                &dna[(p - crop_start) as usize..(p + k - crop_start) as usize],
                            ),
                            "source DNA contradicts panel dictionary anchor",
                        )?;
                    }
                }
                let context = profile::SourceContext {
                    source: core.source,
                    crop_start,
                    crop_end,
                    views: profile::raw_views(panel, &dna)?,
                };
                let mut hits = Vec::new();
                for (view, walk) in context.views.iter().enumerate() {
                    // Halos are read-only context, never additional ownership or
                    // discovery responsibility. Endpoints contain core anchors only.
                    let walk: Vec<_> = walk
                        .iter()
                        .copied()
                        .map(|(n, p)| (n, p + crop_start))
                        .filter(|&(_, p)| p >= start && p < end)
                        .collect();
                    if let (Some(&first), Some(&last)) = (walk.first(), walk.last()) {
                        endpoints.push(Endpoint {
                            source: core.source,
                            start,
                            view,
                            first,
                            last,
                        });
                    }
                    for pair in walk.windows(2) {
                        if let Some(h) = candidate(pair, core.source, k, &definitions, &lookup)? {
                            hits.push(h);
                        }
                    }
                }
                let name = format!("shards/partition-{g:08}-core-{core_count:08}.jsonl");
                let path = out.join(&name);
                if !hits.is_empty() {
                    physical += shard(
                        &path,
                        hits,
                        panel,
                        &sources,
                        &catalog,
                        &cores,
                        &lengths,
                        &definitions,
                        &mut events,
                        Some(&context),
                        &mut profile_stats,
                    )?;
                    checkpoint.push(fingerprint(&path)?);
                    shards.push(path);
                }
                core_count += 1;
                start = end;
            }
        }
        super::write_json(
            &out.join(format!("shards/partition-{g:08}.json")),
            &json!({"group":g,"id":catalog.groups[g].id,"completed":true,"shards":checkpoint}),
        )?;
    }
    // The first/last anchor of each nonempty core in each actual raw orientation
    // links consecutive cores even across long anchor deserts. No chromosome DNA
    // or chromosome anchor vector is materialized, and no anchor streams mix.
    endpoints.sort_by_key(|e| (e.source, e.view, e.start));
    let mut hits = Vec::new();
    let mut boundary_shard = 0;
    for pair in endpoints.windows(2) {
        if pair[0].source == pair[1].source && pair[0].view == pair[1].view {
            if let Some(h) = candidate(
                &[pair[0].last, pair[1].first],
                pair[0].source,
                k,
                &definitions,
                &lookup,
            )? {
                hits.push(h);
            }
        }
        if hits.len() >= 16384 {
            let path = out.join(format!("shards/endpoints-{boundary_shard:08}.jsonl"));
            physical += shard(
                &path,
                std::mem::take(&mut hits),
                panel,
                &sources,
                &catalog,
                &cores,
                &lengths,
                &definitions,
                &mut events,
                None,
                &mut profile_stats,
            )?;
            shards.push(path);
            boundary_shard += 1;
        }
    }
    if !hits.is_empty() {
        let path = out.join(format!("shards/endpoints-{boundary_shard:08}.jsonl"));
        physical += shard(
            &path,
            hits,
            panel,
            &sources,
            &catalog,
            &cores,
            &lengths,
            &definitions,
            &mut events,
            None,
            &mut profile_stats,
        )?;
        shards.push(path);
    }
    let shard_count = shards.len();
    let mut pass = 0;
    while shards.len() > MERGE_FAN_IN {
        let mut merged = Vec::new();
        for (i, chunk) in shards.chunks(MERGE_FAN_IN).enumerate() {
            let path = out.join(format!("shards/merge-{pass:04}-{i:08}.jsonl"));
            merge(chunk, &path)?;
            merged.push(path);
        }
        // Preserve original checkpoint shards, remove only previous scratch runs.
        if pass > 0 {
            for path in &shards {
                fs::remove_file(path)?;
            }
        }
        shards = merged;
        pass += 1;
    }
    merge(&shards, &out.join("incidences.jsonl"))?;
    if pass > 0 {
        for path in &shards {
            fs::remove_file(path)?;
        }
    }
    for (path, expected) in source_paths.iter().zip(&source_fingerprints) {
        verify(Path::new(path), expected)?;
    }
    let mut files = BTreeMap::new();
    for name in [
        "definitions.jsonl",
        "incidences.jsonl",
        "ownership.json",
        "inputs.json",
    ] {
        files.insert(name.into(), fingerprint(&out.join(name))?);
    }
    let shard_bytes: u64 = fs::read_dir(out.join("shards"))?
        .map(|e| e.and_then(|e| e.metadata()).map(|m| m.len()))
        .collect::<io::Result<Vec<_>>>()?
        .into_iter()
        .sum();
    let stats = json!({"checkpoint_shard_bytes":shard_bytes,"merged_incidence_bytes":fs::metadata(out.join("incidences.jsonl"))?.len(),"partitions_scanned":catalog.groups.len(),"sources_scanned":catalog.sources.len(),"source_bp":catalog.sources.iter().map(|s| s.length).sum::<u64>(),
        "selected_features":definitions.len(),"full_feature_scope":scope.is_empty(),"core_chunks":core_count,"core_context_fetches":core_count,"core_context_bp_fetched":context_bp_fetched,
        "profile_execution":profile_stats,"raw_endpoint_records":endpoints.len(),
        "physical_shard_records_before_merge":physical,"exact_profile_events":events,"shards":shard_count,"core_bp_cap":CORE_BP,"merge_fan_in":MERGE_FAN_IN,
        "resume":"unavailable; fresh-only; failed/running output is not reusable"});
    let metadata = Metadata { version:FORMAT_VERSION,model:MODEL.into(),count_policy:COUNT_POLICY.into(),compiler_identity:compiler_identity(),
        capability:"error-free-source-replay; exact-local-raw-window restriction; full-source conditional profiles; all admissible spanning contexts required for local validity; no noisy-read completeness or novel-junction solver".into(),
        catalog,legacy_input,sources:source_fingerprints,syncmer_length_bp:k,read_lengths:lengths,feature_groups:scope.to_vec(),full_feature_scope:scope.is_empty(),feature_count:definitions.len(),files,stats:stats.clone() };
    super::write_json(&out.join("metadata.json"), &metadata)?;
    Ok(json!({"model":MODEL,"metadata":fingerprint(&out.join("metadata.json"))?,"compiler":stats}))
}

pub fn load(out: &Path, identity: &PanelIdentity) -> io::Result<(Metadata, String)> {
    let manifest: Value = super::read_json(&out.join("manifest.json"))?;
    require(
        manifest["status"] == "succeeded" && manifest["model"] == MODEL,
        "observations incomplete/failed/incompatible",
    )?;
    verify(&out.join("metadata.json"), &manifest["stats"]["metadata"])?;
    let metadata: Metadata = super::read_json(&out.join("metadata.json"))?;
    require(
        metadata.version == FORMAT_VERSION
            && metadata.model == MODEL
            && metadata.count_policy == COUNT_POLICY
            && metadata.catalog.panel == *identity
            && metadata.compiler_identity == compiler_identity(),
        "incompatible observation model/panel/compiler",
    )?;
    require(
        metadata.full_feature_scope == metadata.feature_groups.is_empty(),
        "inconsistent feature scope",
    )?;
    require(
        metadata.syncmer_length_bp > 0
            && metadata.syncmer_length_bp <= i32::MAX as u64
            && !metadata.read_lengths.is_empty()
            && metadata.read_lengths.len() <= 32
            && metadata
                .read_lengths
                .iter()
                .all(|&l| l > 0 && l <= MAX_LENGTH)
            && metadata.read_lengths.windows(2).all(|w| w[0] < w[1]),
        "invalid compiled lengths",
    )?;
    for name in [
        "definitions.jsonl",
        "incidences.jsonl",
        "ownership.json",
        "inputs.json",
    ] {
        verify(
            &out.join(name),
            metadata
                .files
                .get(name)
                .ok_or_else(|| invalid("missing artifact fingerprint"))?,
        )?;
    }
    metadata.catalog.validate()?;
    let checksum = manifest["stats"]["metadata"]["fnv1a64"]
        .as_str()
        .ok_or_else(|| invalid("missing metadata checksum"))?
        .to_string();
    Ok((metadata, checksum))
}

#[cfg(test)]
mod tests;
