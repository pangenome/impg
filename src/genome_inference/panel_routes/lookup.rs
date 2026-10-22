//! The call-path route inventory derived by LOOKUP (the re-anchor
//! principle applied to the routing layer).
//!
//! Survey: the calling machinery (`call-diplotypes`, `chain-molecules`,
//! and the derive-cache builders) consumes from a routes directory ONLY
//! `graph.json`'s lane inventory (panel path id/name/length) plus the
//! source paths it opens `Sources` with — the record->locus routing the
//! call consumes lives in the census's own `partitions` field (already
//! re-anchored by `reanchor-census`). `build-panel-routes` spends its
//! whole wall deriving native ports and per-read-length start profiles
//! — the route-search machinery's artifacts, which the calling path
//! never reads.
//!
//! `derive-routes` therefore emits the call-consumable inventory
//! directly from the target panel's own name map (BY NAME, the same
//! authority the re-anchor maps occurrence paths by), seals honest
//! empty native artifacts (zero ports, zero profile runs — stated in
//! the receipt and in `generation_rule`), and cross-checks the
//! re-anchored census: every occurrence's path id must resolve in the
//! target panel's name map (the routing bind the call would otherwise
//! fail late). The artifact FAILS CLOSED under `Graph::load` (the
//! search machinery's validate pins `generation_rule == RULE`); the
//! calling path's own contract — `graph.panel == identity &&
//! graph.generation_complete` — is satisfied by construction.
use super::*;
use rayon::prelude::*;
use std::io::BufRead;
use storage::Seal;

/// The honest generation rule of the lookup-derived artifact (distinct
/// from the search machinery's RULE by design; `Graph::validate` fails
/// closed on it).
pub const LOOKUP_RULE: &str = "call-path-lane-inventory-lookup;panel-name-map-by-name;census-routing-bind-verified;zero-native-ports;zero-native-profiles;consumed-by-call-chain-derive-builders-only";
pub const LOOKUP_ENDPOINTS: &str =
    "lookup-derived-call-path-lane-inventory;no-endpoint-ports-derived";

#[derive(Deserialize)]
struct CensusOccurrenceLine {
    path: usize,
    #[serde(default)]
    partitions: Vec<u32>,
}

#[derive(Deserialize)]
struct CensusRecordLine {
    #[serde(default)]
    occurrences: Vec<CensusOccurrenceLine>,
}

pub fn lookup(
    panel: &SyngIndex,
    identity: PanelIdentity,
    census: &Path,
    source_paths: &[String],
    mut lengths: Vec<u64>,
    out: &Path,
    core_bp: u64,
) -> io::Result<serde_json::Value> {
    require(
        core_bp > 0 && core_bp <= 1048576,
        "invalid route compiler resource bounds",
    )?;
    lengths.sort_unstable();
    lengths.dedup();
    require(
        !lengths.is_empty()
            && lengths.len() <= 32
            && lengths.iter().all(|&l| l > 0 && l <= 1048576),
        "invalid route read lengths",
    )?;
    let absolute_sources = source_paths
        .iter()
        .map(|p| fs::canonicalize(p).map(|p| p.to_string_lossy().into_owned()))
        .collect::<io::Result<Vec<_>>>()?;
    let source_paths = absolute_sources.as_slice();
    let k = panel.syncmer_length_bp() as u64;
    require(k >= 2 && k <= 1048576, "unsupported switch DNA word width")?;
    fs::create_dir(out.join("native"))?;
    fs::create_dir(out.join("source-ports"))?;
    let source_access = graph::bind_source_access(source_paths)?;
    let fingerprints = source_paths
        .iter()
        .map(|p| super::super::reconstruction::fingerprint(Path::new(p)))
        .collect::<io::Result<Vec<_>>>()?;
    // THE LANE INVENTORY BY NAME: the target panel's own name map is the
    // authority (the same rule build asserts against the catalog: every
    // panel source path, in path-id order, with the panel's own lengths).
    let lanes: Vec<(usize, String, u64)> = panel
        .name_map
        .path_to_name
        .iter()
        .enumerate()
        .map(|(id, name)| {
            let length = panel.name_map.path_to_length[id];
            require(
                length > 0 && length <= i32::MAX as u64,
                "unsupported source length",
            )?;
            Ok((id, name.clone(), length))
        })
        .collect::<io::Result<Vec<_>>>()?;
    require(!lanes.is_empty(), "empty panel name map")?;
    let sources = graph::Sources::open(
        source_paths,
        lanes
            .iter()
            .map(|&(_, ref name, length)| (name.clone(), length))
            .collect(),
    )?;
    // The topology families (the same identity rule build uses).
    let mut family_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for &(id, ref name, _) in &lanes {
        family_map
            .entry(super::super::genotype::source_identity(name)?)
            .or_default()
            .push(id);
    }
    let families: Vec<Family> = family_map
        .into_iter()
        .map(|(identity, paths)| Family { identity, paths })
        .collect();
    let mut source_family = vec![0usize; lanes.len()];
    for (f, family) in families.iter().enumerate() {
        for &s in &family.paths {
            source_family[s] = f;
        }
    }
    // THE CENSUS ROUTING CROSS-CHECK: every re-anchored occurrence's
    // path id must resolve in the target panel's name map (the routing
    // lookup the call performs through the lane inventory).
    let mut census_records = 0u64;
    let mut census_occurrences = 0u64;
    let mut census_partitions: BTreeSet<u32> = BTreeSet::new();
    let mut census_paths: BTreeSet<usize> = BTreeSet::new();
    {
        let (reader, _) =
            niffler::get_reader(Box::new(File::open(census)?)).map_err(io::Error::other)?;
        for line in BufReader::new(reader).lines() {
            let record: CensusRecordLine = serde_json::from_str(&line?)?;
            census_records += 1;
            for occurrence in &record.occurrences {
                require(
                    occurrence.path < lanes.len(),
                    "census occurrence path has no panel name-map entry",
                )?;
                census_occurrences += 1;
                census_paths.insert(occurrence.path);
                census_partitions.extend(occurrence.partitions.iter().copied());
            }
        }
    }
    require(census_records > 0, "census carries no records")?;
    // THE LANES: the sequence hash recomputed over the same source
    // bytes build hashes (one linear read per lane; no syncmer scan,
    // no MEM event enumeration), honest empty port files sealed.
    // (PARALLEL ACROSS LANES: each lane's hash is a pure function of
    // its own bytes — the per-lane crop loop is unchanged, the
    // per-lane files are distinct, and the fetch path is thread-safe
    // by construction (Sources::fetch is &self; the AGC keeps
    // per-rayon-thread decompressor pools; core_bp crops exceed the
    // small-crop fetch memo) — so the lanes hash concurrently with
    // byte-identical results, and `collect` preserves lane order.
    // The serial form measured ~55 Mbp/s over the whole 1.403 Tbp /
    // 38,790-lane committed panel — a ~7.5h panel-side wall; the
    // parallel form emits the SAME artifact at the fleet's thread
    // count, and the serial run's graph.json remains the cmp
    // receipt. Serves the whole-panel routes inventory (lever 1)
    // and the shared-index cut lanes (lever 2) alike.)
    let mut total_source_bp = 0u64;
    for &(_id, ref _name, length) in &lanes {
        add(&mut total_source_bp, length)?;
    }
    let graph_lanes: Vec<Lane> = lanes
        .par_iter()
        .map(|&(id, ref name, length)| -> io::Result<Lane> {
            let path = format!("source-ports/{id}.bin");
            let mut ports = BufWriter::new(File::create(out.join(&path))?);
            ports.flush()?;
            let mut sequence_hash = 0xcbf29ce484222325u64;
            let mut lo = 0u64;
            while lo < length {
                let end = lo.saturating_add(core_bp).min(length);
                let crop = end.saturating_add(k - 1).min(length);
                let dna = sources.fetch(id, lo, crop)?;
                super::super::hash_update(&mut sequence_hash, &dna[..(end - lo) as usize]);
                lo = end;
            }
            let mut profiles = Vec::with_capacity(lengths.len());
            for &read_length in &lengths {
                let stem = format!("native/{id}-{read_length}");
                for suffix in ["events.jsonl", "index", "totals.jsonl"] {
                    let file = format!("{stem}.{suffix}");
                    BufWriter::new(File::create(out.join(&file))?).flush()?;
                }
                profiles.push(profiles::NativeIndex {
                    length: read_length,
                    starts: if read_length <= length {
                        length - read_length + 1
                    } else {
                        0
                    },
                    runs: 0,
                    events: Seal::create(out, &format!("{stem}.events.jsonl"))?,
                    index: Seal::create(out, &format!("{stem}.index"))?,
                    totals: Seal::create(out, &format!("{stem}.totals.jsonl"))?,
                });
            }
            Ok(Lane {
                id,
                name: name.clone(),
                length,
                sequence_fnv1a64: format!("{sequence_hash:016x}"),
                family: source_family[id],
                ports: Seal::create(out, &path)?,
                port_count: 0,
                profiles,
            })
        })
        .collect::<io::Result<Vec<_>>>()?;
    // Honest empty native artifacts, sealed and counted.
    BufWriter::new(File::create(out.join("ports.bin"))?).flush()?;
    BufWriter::new(File::create(out.join("ownership.jsonl"))?).flush()?;
    BufWriter::new(File::create(out.join("ownership-input.json"))?).flush()?;
    BufWriter::new(File::create(out.join("registry.jsonl"))?).flush()?;
    for (p, before) in source_paths.iter().zip(&fingerprints) {
        require(
            *before == super::super::reconstruction::fingerprint(Path::new(p))?,
            "source mutation during lookup derivation",
        )?;
    }
    graph::verify_source_access(source_paths, &source_access)?;
    let graph = Graph {
        version: VERSION,
        model: MODEL.into(),
        panel: identity,
        compiler_identity: compiler_identity(),
        count_policy: COUNT_POLICY.into(),
        generation_rule: LOOKUP_RULE.into(),
        generation_complete: true,
        endpoint_semantics: LOOKUP_ENDPOINTS.into(),
        k,
        cut_offset: k / 2,
        core_bp,
        read_lengths: lengths.clone(),
        lanes: graph_lanes,
        families,
        ports: Seal::create(out, "ports.bin")?,
        port_count: 0,
        ownership: Seal::create(out, "ownership.jsonl")?,
        original_ownership: Seal::create(out, "ownership-input.json")?,
        registry: Seal::create(out, "registry.jsonl")?,
        registry_count: 0,
        catalog_provenance: serde_json::json!({
            "derivation": "derive-routes lookup",
            "native_ports": false,
            "native_profiles": false,
            "routing": "the census partitions field (re-anchored by reanchor-census)"
        }),
        source_files: fingerprints,
        source_paths: source_paths.to_vec(),
        source_access,
        total_source_bp,
    };
    write_json(&out.join("graph.json"), &graph)?;
    Ok(serde_json::json!({
        "graph": Seal::create(out, "graph.json")?,
        "derivation": "lookup",
        "generation_rule": LOOKUP_RULE,
        "sources": graph.lanes.len(),
        "source_bp": graph.total_source_bp,
        "families": graph.families.len(),
        "census_records": census_records,
        "census_occurrences": census_occurrences,
        "census_paths": census_paths.len(),
        "census_partitions": census_partitions.len(),
        "ports_derived": 0,
        "native_profiles_derived": 0,
        "routing_bind": "every census occurrence path id resolves in the target panel name map",
        "consumers": "call-diplotypes, chain-molecules, derive-cache builders (graph.json lanes + source paths only)",
        "search_machinery": "rejected by design (Graph::load pins the built generation rule)"
    }))
}
