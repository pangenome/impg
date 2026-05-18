use crate::commands::{partition, syng2gfa};
use crate::graph::reverse_complement;
use crate::render_bundle::{
    self, RenderBundlePaths, RenderManifest, RenderTranslationTables, RenderedPathRecord,
    StepTranslationRecord,
};
use crate::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use crate::sequence_namespace::{SequenceNamespace, SourceInterval};
use crate::syng::{self, SyngIndex};
use log::{info, warn};
use std::collections::BTreeSet;
use std::fs;
use std::io;
use std::path::Path;

#[derive(Debug, Clone)]
pub struct RenderConfig<'a> {
    pub syng_prefix: &'a str,
    pub target_range: &'a str,
    pub output: &'a str,
    pub sequence_files: &'a [String],
    pub engine: &'a str,
    pub syng_padding: u64,
    pub syng_extension: u64,
    pub emit_gfa: bool,
    pub keep_existing: bool,
}

#[derive(Debug, Clone)]
struct RenderInterval {
    source_name: String,
    start: u64,
    end: u64,
    strand: char,
}

pub fn run(config: &RenderConfig<'_>) -> io::Result<()> {
    match config.engine {
        "syng-native" => render_syng_native(config),
        other => Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("unsupported render engine '{other}'; currently expected 'syng-native'"),
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

    info!(
        "Building sequence index for {} source file(s)",
        config.sequence_files.len()
    );
    let sequence_index = UnifiedSequenceIndex::from_files(config.sequence_files)?;

    let intervals = collect_render_intervals(
        &index,
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
                    format!(
                        "source '{}' is missing from namespace",
                        interval.source_name
                    ),
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
            gbwt_path_id: Some(rendered_path_id),
        });
        fetched.push((rendered_name, seq));
    }

    let fetched_refs: Vec<(String, &[u8])> = fetched
        .iter()
        .map(|(name, seq)| (name.clone(), seq.as_slice()))
        .collect();
    let render_index = SyngIndex::new(index.params);
    let syng_prefix = paths.syng_prefix.to_string_lossy().to_string();
    render_index.build_region_gbwt(&fetched_refs, &syng_prefix)?;

    let region_index = SyngIndex::load(&syng_prefix, syng::SyncmerParams::default())?;
    let step_samples = collect_region_step_samples(&region_index, &rendered_paths)?;
    let tables = RenderTranslationTables {
        namespace,
        rendered_paths,
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
        )?;
    }

    let manifest = RenderManifest::new_syng_native(
        config.syng_prefix.to_string(),
        config.target_range.to_string(),
        &paths,
        config.emit_gfa,
        tables.namespace.sequences.len(),
        tables.rendered_paths.len(),
        tables.step_samples.len(),
    );
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
