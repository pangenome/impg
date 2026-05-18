use crate::commands::{graph, partition, syng2gfa};
use crate::graph::reverse_complement;
use crate::render_bundle::{
    self, RenderBundlePaths, RenderManifest, RenderTranslationTables, RenderedPathRecord,
    StepTranslationRecord,
};
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::sequence_namespace::{SequenceNamespace, SourceInterval};
use crate::syng::{self, SyngIndex};
use log::{info, warn};
use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::Path;

#[derive(Debug, Clone)]
pub struct RenderConfig<'a> {
    pub syng_prefix: &'a str,
    pub target_range: &'a str,
    pub output: &'a str,
    pub sequence_files: &'a [String],
    pub engine: &'a str,
    pub syng_gfa_mode: syng2gfa::SyngGfaMode,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub emit_gfa: bool,
    pub keep_existing: bool,
    pub threads: usize,
}

#[derive(Debug, Clone)]
struct RenderInterval {
    source_name: String,
    start: u64,
    end: u64,
    strand: char,
}

struct RenderInputs {
    namespace: SequenceNamespace,
    rendered_paths: Vec<RenderedPathRecord>,
    fetched: Vec<(String, Vec<u8>)>,
}

pub fn run(config: &RenderConfig<'_>) -> io::Result<()> {
    match normalize_render_engine(config.engine).as_str() {
        "syng-native" => render_syng_native(config),
        "poa" | "seqwish" | "pggb" => render_local_graph(config),
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "unsupported render engine '{other}'; expected one of: syng-native, poa, seqwish, pggb"
            ),
        )),
    }
}

fn render_syng_native(config: &RenderConfig<'_>) -> io::Result<()> {
    if config.sequence_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "impg render --engine syng-native currently requires --sequence-files or --sequence-list",
        ));
    }

    let paths = RenderBundlePaths::new(config.output);
    prepare_bundle_dir(&paths.root, config.keep_existing)?;

    info!("Loading syng index from prefix: {}", config.syng_prefix);
    let index = SyngIndex::load(config.syng_prefix, syng::SyncmerParams::default())?;

    let mut inputs = collect_render_inputs(&index, config)?;
    for (path_idx, rendered_path) in inputs.rendered_paths.iter_mut().enumerate() {
        rendered_path.gbwt_path_id = Some(path_idx as u32);
    }
    write_rendered_fasta(&paths.rendered_fasta, &inputs.fetched)?;

    let fetched_refs: Vec<(String, &[u8])> = inputs
        .fetched
        .iter()
        .map(|(name, seq)| (name.clone(), seq.as_slice()))
        .collect();
    let render_index = SyngIndex::new(index.params);
    let syng_prefix = paths.syng_prefix.to_string_lossy().to_string();
    render_index.build_region_gbwt(&fetched_refs, &syng_prefix)?;

    let region_index = SyngIndex::load(&syng_prefix, syng::SyncmerParams::default())?;
    let step_samples = collect_region_step_samples(&region_index, &inputs.rendered_paths)?;
    let tables = RenderTranslationTables {
        namespace: inputs.namespace,
        rendered_paths: inputs.rendered_paths,
        step_samples,
    };

    render_bundle::write_namespace(&paths.namespace, &tables.namespace)?;
    render_bundle::write_translation_binary(&paths.translation, &tables)?;
    render_bundle::write_translation_tsv(&paths.translation_tsv, &tables)?;

    if config.emit_gfa {
        syng2gfa::run(
            &syng_prefix,
            paths.graph_gfa.to_str().ok_or_else(|| {
                io::Error::new(io::ErrorKind::InvalidInput, "graph GFA path is not UTF-8")
            })?,
            syng2gfa::GfaVersion::V1_0,
            &[],
            config.syng_gfa_mode,
        )?;
    }

    let mut manifest = RenderManifest::new_syng_native(
        config.syng_prefix.to_string(),
        config.target_range.to_string(),
        &paths,
        config.emit_gfa,
        tables.namespace.sequences.len(),
        tables.rendered_paths.len(),
        tables.step_samples.len(),
    );
    manifest.engine = format!("syng:{}", config.syng_gfa_mode.label());
    render_bundle::write_manifest(&paths.manifest, &manifest)?;

    info!(
        "Wrote syng-native render bundle to {}: {} source sequences, {} rendered paths, {} step samples",
        paths.root.display(),
        manifest.source_sequences,
        manifest.rendered_paths,
        manifest.step_samples
    );
    Ok(())
}

fn render_local_graph(config: &RenderConfig<'_>) -> io::Result<()> {
    if config.sequence_files.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "impg render local graph engines require --sequence-files or --sequence-list",
        ));
    }
    if !config.emit_gfa {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            "--no-gfa is not supported with local graph render engines",
        ));
    }

    let engine = normalize_render_engine(config.engine);
    let paths = RenderBundlePaths::new(config.output);
    prepare_bundle_dir(&paths.root, config.keep_existing)?;

    info!("Loading syng index from prefix: {}", config.syng_prefix);
    let index = SyngIndex::load(config.syng_prefix, syng::SyncmerParams::default())?;
    let inputs = collect_render_inputs(&index, config)?;
    write_rendered_fasta(&paths.rendered_fasta, &inputs.fetched)?;
    build_local_graph(&engine, &paths.rendered_fasta, &paths.graph_gfa, config.threads)?;
    let step_samples = collect_gfa_step_samples(&paths.graph_gfa, &inputs.rendered_paths)?;

    let tables = RenderTranslationTables {
        namespace: inputs.namespace,
        rendered_paths: inputs.rendered_paths,
        step_samples,
    };
    render_bundle::write_namespace(&paths.namespace, &tables.namespace)?;
    render_bundle::write_translation_binary(&paths.translation, &tables)?;
    render_bundle::write_translation_tsv(&paths.translation_tsv, &tables)?;

    let manifest = RenderManifest::new_local_graph(
        engine,
        config.syng_prefix.to_string(),
        config.target_range.to_string(),
        &paths,
        tables.namespace.sequences.len(),
        tables.rendered_paths.len(),
        tables.step_samples.len(),
    );
    render_bundle::write_manifest(&paths.manifest, &manifest)?;

    info!(
        "Wrote local graph render bundle to {}: {} source sequences, {} rendered paths, {} step samples",
        paths.root.display(),
        manifest.source_sequences,
        manifest.rendered_paths,
        manifest.step_samples
    );
    Ok(())
}

fn prepare_bundle_dir(root: &Path, keep_existing: bool) -> io::Result<()> {
    if root.exists() {
        if !keep_existing {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!(
                    "render output '{}' already exists; pass --keep-existing to write into it",
                    root.display()
                ),
            ));
        }
        if !root.is_dir() {
            return Err(io::Error::new(
                io::ErrorKind::AlreadyExists,
                format!(
                    "render output '{}' exists and is not a directory",
                    root.display()
                ),
            ));
        }
        warn!(
            "Writing render bundle into existing directory {}",
            root.display()
        );
    } else {
        fs::create_dir_all(root)?;
    }
    Ok(())
}

fn normalize_render_engine(engine: &str) -> String {
    match engine.trim().replace('_', "-").as_str() {
        "syng" | "syng-native" => "syng-native".to_string(),
        other => other.to_string(),
    }
}

fn collect_render_inputs(index: &SyngIndex, config: &RenderConfig<'_>) -> io::Result<RenderInputs> {
    info!(
        "Building sequence index for {} source file(s)",
        config.sequence_files.len()
    );
    let sequence_index = UnifiedSequenceIndex::from_files(config.sequence_files)?;

    let intervals = collect_render_intervals(
        index,
        config.target_range,
        config.syng_padding,
        config.syng_extension,
    )?;
    if intervals.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "no source intervals found for render target",
        ));
    }

    let mut namespace = SequenceNamespace::new();
    for (name, &length) in index
        .name_map
        .path_to_name
        .iter()
        .zip(index.name_map.path_to_length.iter())
    {
        namespace.add_sequence(name.clone(), length);
    }

    let mut fetched = Vec::with_capacity(intervals.len());
    let mut rendered_paths = Vec::with_capacity(intervals.len());
    for interval in &intervals {
        let source = namespace
            .get_by_name(&interval.source_name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("source '{}' is missing from namespace", interval.source_name),
                )
            })?;
        let mut seq = fetch_source_interval(
            &sequence_index,
            &interval.source_name,
            interval.start,
            interval.end,
        )?;
        if interval.strand == '-' {
            seq = reverse_complement(&seq);
        }
        let rendered_path_id = rendered_paths.len() as u32;
        let rendered_name = format!(
            "{}:{}-{}({})",
            interval.source_name, interval.start, interval.end, interval.strand
        );
        rendered_paths.push(RenderedPathRecord {
            rendered_path_id,
            rendered_name: rendered_name.clone(),
            source_interval: SourceInterval::new(
                source.id,
                interval.start,
                interval.end,
                interval.strand,
            )?,
            gbwt_path_id: None,
        });
        fetched.push((rendered_name, seq));
    }

    Ok(RenderInputs {
        namespace,
        rendered_paths,
        fetched,
    })
}

fn write_rendered_fasta(path: &Path, sequences: &[(String, Vec<u8>)]) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    for (name, seq) in sequences {
        writeln!(writer, ">{name}")?;
        for chunk in seq.chunks(80) {
            writer.write_all(chunk)?;
            writer.write_all(b"\n")?;
        }
    }
    writer.flush()
}

fn build_local_graph(engine: &str, fasta: &Path, gfa: &Path, threads: usize) -> io::Result<()> {
    let fasta_files = vec![fasta.to_string_lossy().to_string()];
    let mut config = graph::GraphBuildConfig::default();
    config.num_threads = threads.max(1);
    config.show_progress = false;

    match engine {
        "poa" => {
            let mut writer = BufWriter::new(File::create(gfa)?);
            graph::run_graph_build_poa(
                fasta_files,
                &mut writer,
                (5, 4, 6, 2, 24, 1),
                &config,
            )
        }
        "seqwish" => {
            let gfa_path = gfa.to_string_lossy().to_string();
            graph::run_graph_build(fasta_files, &gfa_path, config)
        }
        "pggb" => {
            let mut writer = BufWriter::new(File::create(gfa)?);
            graph::run_graph_build_pggb(
                fasta_files,
                &mut writer,
                &config,
                vec![700, 1100],
                100,
                0.001,
            )
        }
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported local graph render engine '{other}'"),
        )),
    }
}

fn collect_render_intervals(
    index: &SyngIndex,
    target_range: &str,
    padding: u64,
    extension: u64,
) -> io::Result<Vec<RenderInterval>> {
    let (target_name, (start, end), _) = partition::parse_target_range(target_range)?;
    if start < 0 || end <= start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid render target range '{target_range}'"),
        ));
    }
    let hits = index.query_region_with_anchors_ext(
        &target_name,
        start as u64,
        end as u64,
        padding,
        extension,
    )?;
    let mut seen = BTreeSet::new();
    let mut intervals = Vec::new();
    for hit in hits {
        if hit.end <= hit.start {
            continue;
        }
        let key = (hit.genome.clone(), hit.start, hit.end, hit.strand);
        if seen.insert(key.clone()) {
            intervals.push(RenderInterval {
                source_name: key.0,
                start: key.1,
                end: key.2,
                strand: key.3,
            });
        }
    }
    intervals.sort_by(|a, b| {
        a.source_name
            .cmp(&b.source_name)
            .then(a.start.cmp(&b.start))
            .then(a.end.cmp(&b.end))
            .then(a.strand.cmp(&b.strand))
    });
    Ok(intervals)
}

fn fetch_source_interval(
    sequence_index: &UnifiedSequenceIndex,
    source_name: &str,
    start: u64,
    end: u64,
) -> io::Result<Vec<u8>> {
    if start > i32::MAX as u64 || end > i32::MAX as u64 {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("render interval exceeds i32 coordinate range: {source_name}:{start}-{end}"),
        ));
    }
    sequence_index.fetch_sequence(source_name, start as i32, end as i32)
}

fn collect_region_step_samples(
    index: &SyngIndex,
    rendered_paths: &[RenderedPathRecord],
) -> io::Result<Vec<StepTranslationRecord>> {
    let mut records = Vec::new();
    let syncmer_len = index.syncmer_length_bp() as u64;
    for (path_idx, &length) in index.name_map.path_to_length.iter().enumerate() {
        let rendered_path = rendered_paths.get(path_idx).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("region syng path {path_idx} has no rendered path metadata"),
            )
        })?;
        let walk = index.walk_path_range(path_idx, 0, length)?;
        for (step_idx, (signed_node, source_bp)) in walk.into_iter().enumerate() {
            let source_bp = match rendered_path.source_interval.strand {
                '+' => rendered_path
                    .source_interval
                    .start
                    .saturating_add(source_bp),
                '-' => rendered_path
                    .source_interval
                    .end
                    .saturating_sub(source_bp.saturating_add(syncmer_len)),
                _ => source_bp,
            };
            records.push(StepTranslationRecord {
                rendered_path_id: path_idx as u32,
                rendered_step: step_idx as u32,
                source_bp,
                feature_id: signed_node.unsigned_abs(),
                orientation: if signed_node >= 0 { '+' } else { '-' },
            });
        }
    }
    Ok(records)
}

fn collect_gfa_step_samples(
    gfa: &Path,
    rendered_paths: &[RenderedPathRecord],
) -> io::Result<Vec<StepTranslationRecord>> {
    let mut segment_lengths = HashMap::new();
    let mut path_walks: BTreeMap<String, Vec<(u32, char)>> = BTreeMap::new();
    let reader = BufReader::new(File::open(gfa)?);
    for line in reader.lines() {
        let line = line?;
        if line.is_empty() {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                let segment_id = parse_gfa_segment_id(fields[1])?;
                segment_lengths.insert(segment_id, fields[2].len() as u64);
            }
            Some("P") if fields.len() >= 3 => {
                let walk = parse_gfa_p_walk(fields[2])?;
                path_walks.insert(fields[1].to_string(), walk);
            }
            Some("W") if fields.len() >= 7 => {
                let path_name = format!(
                    "{}#{}#{}:{}-{}(+)",
                    fields[1], fields[2], fields[3], fields[4], fields[5]
                );
                let walk = parse_gfa_w_walk(fields[6])?;
                path_walks.insert(path_name, walk);
            }
            _ => {}
        }
    }

    let mut records = Vec::new();
    for rendered_path in rendered_paths {
        let walk = find_gfa_path_walk(&path_walks, &rendered_path.rendered_name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GFA graph '{}' has no path for rendered interval '{}'",
                        gfa.display(),
                        rendered_path.rendered_name
                    ),
                )
            })?;
        let mut rendered_offset = 0u64;
        for (rendered_step, (feature_id, orientation)) in walk.iter().copied().enumerate() {
            let segment_len = *segment_lengths.get(&feature_id).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GFA path references missing segment {feature_id}"),
                )
            })?;
            let source_bp = match rendered_path.source_interval.strand {
                '+' => rendered_path
                    .source_interval
                    .start
                    .saturating_add(rendered_offset),
                '-' => rendered_path
                    .source_interval
                    .end
                    .saturating_sub(rendered_offset.saturating_add(segment_len)),
                _ => rendered_offset,
            };
            records.push(StepTranslationRecord {
                rendered_path_id: rendered_path.rendered_path_id,
                rendered_step: rendered_step as u32,
                source_bp,
                feature_id,
                orientation,
            });
            rendered_offset = rendered_offset.saturating_add(segment_len);
        }
    }
    Ok(records)
}

fn find_gfa_path_walk<'a>(
    path_walks: &'a BTreeMap<String, Vec<(u32, char)>>,
    rendered_name: &str,
) -> Option<&'a Vec<(u32, char)>> {
    if let Some(walk) = path_walks.get(rendered_name) {
        return Some(walk);
    }

    let metadata_prefix = format!("{rendered_name}:");
    let mut matches = path_walks
        .range(metadata_prefix.clone()..)
        .take_while(|(name, _)| name.starts_with(&metadata_prefix))
        .map(|(_, walk)| walk);
    let first = matches.next()?;
    if matches.next().is_none() {
        Some(first)
    } else {
        None
    }
}

fn parse_gfa_p_walk(walk: &str) -> io::Result<Vec<(u32, char)>> {
    if walk == "*" || walk.is_empty() {
        return Ok(Vec::new());
    }
    walk.split(',')
        .map(|step| {
            let (id, orientation) = step.split_at(step.len().saturating_sub(1));
            let orientation = match orientation {
                "+" => '+',
                "-" => '-',
                _ => {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("invalid GFA P-line step orientation in '{step}'"),
                    ));
                }
            };
            Ok((parse_gfa_segment_id(id)?, orientation))
        })
        .collect()
}

fn parse_gfa_w_walk(walk: &str) -> io::Result<Vec<(u32, char)>> {
    let mut steps = Vec::new();
    let mut orientation = None;
    let mut start = 0usize;
    for (idx, ch) in walk.char_indices() {
        if ch == '>' || ch == '<' {
            if let Some(prev_orientation) = orientation {
                let id = &walk[start..idx];
                if !id.is_empty() {
                    steps.push((parse_gfa_segment_id(id)?, prev_orientation));
                }
            }
            orientation = Some(if ch == '>' { '+' } else { '-' });
            start = idx + ch.len_utf8();
        }
    }
    if let Some(prev_orientation) = orientation {
        let id = &walk[start..];
        if !id.is_empty() {
            steps.push((parse_gfa_segment_id(id)?, prev_orientation));
        }
    }
    Ok(steps)
}

fn parse_gfa_segment_id(id: &str) -> io::Result<u32> {
    id.parse::<u32>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "GFA segment id '{id}' is not a u32; render translation currently requires numeric segment ids"
            ),
        )
    })
}
