use super::*;
use storage::{Port, Seal};
/// Sample-independent builder. Resource caps fail explicitly; none are frequency
/// masks or changes to the declared switch/endpoint domain.
pub fn build(
    panel: &SyngIndex,
    identity: PanelIdentity,
    catalog_path: &Path,
    source_paths: &[String],
    mut lengths: Vec<u64>,
    out: &Path,
    core_bp: u64,
    max_terms: usize,
) -> io::Result<serde_json::Value> {
    require(
        core_bp > 0 && core_bp <= 1048576 && max_terms > 0,
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
    let port_buffer_records = 65536usize.min(67108864usize / (k as usize + 64)).max(1);
    require(k >= 2 && k <= 1048576, "unsupported switch DNA word width")?;
    fs::create_dir(out.join("native"))?;
    fs::create_dir(out.join("source-ports"))?;
    let (catalog, registry_count, catalog_provenance) = import::load(catalog_path, &identity, out)?;
    require(
        catalog.sources.len() == panel.name_map.path_to_name.len(),
        "ownership omits panel sources",
    )?;
    for (id, s) in catalog.sources.iter().enumerate() {
        require(
            s.id == id
                && s.path == panel.name_map.path_to_name[id]
                && s.length == panel.name_map.path_to_length[id],
            "catalog/panel native source namespace mismatch",
        )?;
    }
    let ownership = super::super::observations::input::ownership(&catalog)?;
    write_json(&out.join("ownership-input.json"), &catalog)?;
    let mut writer = BufWriter::new(File::create(out.join("ownership.jsonl"))?);
    for cores in &ownership {
        for c in cores {
            storage::line(&mut writer, c)?;
        }
    }
    writer.flush()?;
    drop(ownership);
    let sources = graph::Sources::open(
        source_paths,
        catalog
            .sources
            .iter()
            .map(|s| (s.path.clone(), s.length))
            .collect(),
    )?;
    let source_access = graph::bind_source_access(source_paths)?;
    let fingerprints = source_paths
        .iter()
        .map(|p| super::super::reconstruction::fingerprint(Path::new(p)))
        .collect::<io::Result<Vec<_>>>()?;
    let mut family_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for s in &catalog.sources {
        family_map
            .entry(super::super::genotype::source_identity(&s.path)?)
            .or_default()
            .push(s.id);
    }
    let families: Vec<_> = family_map
        .into_iter()
        .map(|(identity, paths)| Family { identity, paths })
        .collect();
    let mut source_family = vec![0; catalog.sources.len()];
    for (f, family) in families.iter().enumerate() {
        for &s in &family.paths {
            source_family[s] = f;
        }
    }
    let mut runs = Vec::new();
    let mut buffer = Vec::new();
    let flush = |buffer: &mut Vec<Port>, runs: &mut Vec<String>| -> io::Result<()> {
        if buffer.is_empty() {
            return Ok(());
        }
        buffer.sort_unstable();
        let name = format!("port-run-{}.bin", runs.len());
        let mut w = BufWriter::new(File::create(out.join(&name))?);
        for p in buffer.iter() {
            p.write(&mut w, k as usize)?;
        }
        w.flush()?;
        buffer.clear();
        runs.push(name);
        Ok(())
    };
    let mut lanes = Vec::new();
    let mut total_source_bp = 0;
    let mut port_count = 0;
    for s in &catalog.sources {
        add(&mut total_source_bp, s.length)?;
        let path = format!("source-ports/{}.bin", s.id);
        let mut ports = BufWriter::new(File::create(out.join(&path))?);
        let mut count = 0;
        let mut sequence_hash = 0xcbf29ce484222325;
        let mut lo = 0;
        while lo < s.length {
            let end = lo.saturating_add(core_bp).min(s.length);
            let crop = end.saturating_add(k - 1).min(s.length);
            let dna = sources.fetch(s.id, lo, crop)?;
            super::super::hash_update(&mut sequence_hash, &dna[..(end - lo) as usize]);
            let views = super::super::observations::profile::raw_views(panel, &dna)?;
            let positions: BTreeSet<_> = views
                .iter()
                .flat_map(|v| v.iter().map(|&(_, p)| p))
                .filter(|&p| p < end - lo && p + k <= dna.len() as u64)
                .collect();
            for p in positions {
                let word = dna[p as usize..(p + k) as usize].to_vec();
                require(
                    word.iter().all(|b| b"ACGT".contains(b)),
                    "raw anchor does not have a complete unambiguous DNA witness",
                )?;
                for reverse in [false, true] {
                    let port = Port {
                        word: if reverse {
                            crate::graph::reverse_complement(&word)
                        } else {
                            word.clone()
                        },
                        source: s.id,
                        anchor: lo + p,
                        reverse,
                    };
                    port.write(&mut ports, k as usize)?;
                    buffer.push(port);
                    add(&mut count, 1)?;
                    if buffer.len() >= port_buffer_records {
                        flush(&mut buffer, &mut runs)?;
                    }
                }
            }
            lo = end;
        }
        ports.flush()?;
        add(&mut port_count, count)?;
        let mut native = Vec::new();
        for &length in &lengths {
            native.push(profiles::compile(
                out, s.id, s.length, length, panel, &sources, core_bp, max_terms,
            )?);
        }
        lanes.push(Lane {
            id: s.id,
            name: s.path.clone(),
            length: s.length,
            sequence_fnv1a64: format!("{sequence_hash:016x}"),
            family: source_family[s.id],
            ports: Seal::create(out, &path)?,
            port_count: count,
            profiles: native,
        });
        // Durable per-lane checkpoint remains inspectable if later compilation fails.
        write_json(
            &out.join(format!("native/{}-checkpoint.json", s.id)),
            lanes.last().unwrap(),
        )?;
    }
    flush(&mut buffer, &mut runs)?;
    let sorted = storage::merge_ports(out, runs, k as usize)?;
    fs::rename(out.join(sorted), out.join("ports.bin"))?;
    // Validate full membership/order/DNA-derived coordinate invariants without an
    // anchor-pair catalog. Source files and their opposite views were emitted once.
    let mut r = BufReader::new(File::open(out.join("ports.bin"))?);
    let mut previous = None;
    let mut seen = 0;
    while let Some(p) = Port::read(&mut r, k as usize)? {
        require(
            p.source < lanes.len()
                && p.anchor + k <= lanes[p.source].length
                && previous.as_ref().is_none_or(|old| old < &p),
            "duplicate or malformed hub member",
        )?;
        previous = Some(p);
        add(&mut seen, 1)?;
    }
    require(seen == port_count, "incomplete switch hub merge")?;
    for (p, before) in source_paths.iter().zip(&fingerprints) {
        require(
            *before == super::super::reconstruction::fingerprint(Path::new(p))?,
            "source mutation during native compilation",
        )?;
    }
    graph::verify_source_access(source_paths, &source_access)?;
    let graph = Graph {
        version: VERSION,
        model: MODEL.into(),
        panel: identity,
        compiler_identity: compiler_identity(),
        count_policy: COUNT_POLICY.into(),
        generation_rule: RULE.into(),
        generation_complete: true,
        endpoint_semantics:
            "panel-identity-native-linear-assembly-path-end-pairings-not-biological-termini".into(),
        k,
        cut_offset: k / 2,
        core_bp,
        read_lengths: lengths,
        lanes,
        families,
        ports: Seal::create(out, "ports.bin")?,
        port_count,
        ownership: Seal::create(out, "ownership.jsonl")?,
        original_ownership: Seal::create(out, "ownership-input.json")?,
        registry: Seal::create(out, "registry.jsonl")?,
        registry_count,
        catalog_provenance,
        source_files: fingerprints,
        source_paths: source_paths.to_vec(),
        source_access,
        total_source_bp,
    };
    graph.validate()?;
    write_json(&out.join("graph.json"), &graph)?;
    Ok(
        serde_json::json!({"graph":Seal::create(out,"graph.json")?,"sources":graph.lanes.len(),"source_bp":graph.total_source_bp,"families":graph.families.len(),"oriented_ports":graph.port_count,"physical_anchors":graph.port_count/2,"zero_anchor_paths":graph.lanes.iter().filter(|l|l.port_count==0).count(),"generation_complete":true,"all_pairs_edges_materialized":0}),
    )
}
