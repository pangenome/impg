use crate::sequence_namespace::{SequenceNamespace, SourceInterval};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufWriter, Read, Write};
use std::path::{Path, PathBuf};

pub const RENDER_FORMAT: &str = "impg-render-bundle";
pub const RENDER_VERSION: u32 = 1;
pub const TRANSLATION_MAGIC: &[u8; 8] = b"IMPGTRN1";

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RenderManifest {
    pub format: String,
    pub version: u32,
    pub engine: String,
    pub graph_kind: String,
    pub source_index: String,
    pub target_range: String,
    pub feature_space: String,
    pub namespace: String,
    pub translation: String,
    pub translation_tsv: Option<String>,
    pub rendered_fasta: Option<String>,
    pub graph_gfa: Option<String>,
    pub syng_prefix: Option<String>,
    pub source_sequences: usize,
    pub rendered_paths: usize,
    pub step_samples: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RenderedPathRecord {
    pub rendered_path_id: u32,
    pub rendered_name: String,
    pub source_interval: SourceInterval,
    pub gbwt_path_id: Option<u32>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct StepTranslationRecord {
    pub rendered_path_id: u32,
    pub rendered_step: u32,
    pub source_bp: u64,
    pub feature_id: u32,
    pub orientation: char,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RenderTranslationTables {
    pub namespace: SequenceNamespace,
    pub rendered_paths: Vec<RenderedPathRecord>,
    pub step_samples: Vec<StepTranslationRecord>,
}

#[derive(Debug, Clone)]
pub struct RenderBundlePaths {
    pub root: PathBuf,
    pub manifest: PathBuf,
    pub namespace: PathBuf,
    pub translation: PathBuf,
    pub translation_tsv: PathBuf,
    pub rendered_fasta: PathBuf,
    pub graph_gfa: PathBuf,
    pub syng_prefix: PathBuf,
}

#[derive(Debug, Clone)]
pub struct LoadedRenderBundle {
    pub root: PathBuf,
    pub manifest: RenderManifest,
    pub tables: RenderTranslationTables,
}

impl RenderManifest {
    pub fn new_syng_native(
        source_index: String,
        target_range: String,
        paths: &RenderBundlePaths,
        include_gfa: bool,
        source_sequences: usize,
        rendered_paths: usize,
        step_samples: usize,
    ) -> Self {
        Self {
            format: RENDER_FORMAT.to_string(),
            version: RENDER_VERSION,
            engine: "syng-native".to_string(),
            graph_kind: "syng-syncmer-graph".to_string(),
            source_index,
            target_range,
            feature_space: "syng-syncmer-node".to_string(),
            namespace: relative_file(&paths.root, &paths.namespace),
            translation: relative_file(&paths.root, &paths.translation),
            translation_tsv: Some(relative_file(&paths.root, &paths.translation_tsv)),
            rendered_fasta: Some(relative_file(&paths.root, &paths.rendered_fasta)),
            graph_gfa: include_gfa.then(|| relative_file(&paths.root, &paths.graph_gfa)),
            syng_prefix: Some(relative_file(&paths.root, &paths.syng_prefix)),
            source_sequences,
            rendered_paths,
            step_samples,
        }
    }

    pub fn new_local_graph(
        engine: String,
        source_index: String,
        target_range: String,
        paths: &RenderBundlePaths,
        source_sequences: usize,
        rendered_paths: usize,
        step_samples: usize,
    ) -> Self {
        Self {
            format: RENDER_FORMAT.to_string(),
            version: RENDER_VERSION,
            engine,
            graph_kind: "local-sequence-graph".to_string(),
            source_index,
            target_range,
            feature_space: "gfa-segment".to_string(),
            namespace: relative_file(&paths.root, &paths.namespace),
            translation: relative_file(&paths.root, &paths.translation),
            translation_tsv: Some(relative_file(&paths.root, &paths.translation_tsv)),
            rendered_fasta: Some(relative_file(&paths.root, &paths.rendered_fasta)),
            graph_gfa: Some(relative_file(&paths.root, &paths.graph_gfa)),
            syng_prefix: None,
            source_sequences,
            rendered_paths,
            step_samples,
        }
    }
}

impl RenderBundlePaths {
    pub fn new(root: impl AsRef<Path>) -> Self {
        let root = root.as_ref().to_path_buf();
        Self {
            manifest: root.join("manifest.json"),
            namespace: root.join("namespace.json"),
            translation: root.join("translation.bin"),
            translation_tsv: root.join("translation.tsv"),
            rendered_fasta: root.join("rendered.fa"),
            graph_gfa: root.join("graph.gfa"),
            syng_prefix: root.join("paths"),
            root,
        }
    }
}

impl LoadedRenderBundle {
    pub fn syng_prefix_path(&self) -> io::Result<PathBuf> {
        let syng_prefix = self.manifest.syng_prefix.as_deref().ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                "render bundle does not contain a syng prefix",
            )
        })?;
        Ok(resolve_bundle_path(&self.root, syng_prefix))
    }

    pub fn rendered_target_range(&self, source_range: Option<&str>) -> io::Result<String> {
        let source_range = source_range.unwrap_or(&self.manifest.target_range);
        let (source_name, source_start, source_end) = parse_source_range(source_range)?;
        let source = self
            .tables
            .namespace
            .get_by_name(&source_name)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("source range references unknown sequence '{source_name}'"),
                )
            })?;

        let mut best = None;
        for path in &self.tables.rendered_paths {
            let interval = &path.source_interval;
            if interval.source_sequence_id != source.id {
                continue;
            }
            if interval.start > source_start || interval.end < source_end {
                continue;
            }
            let span = interval.end.saturating_sub(interval.start);
            if best
                .as_ref()
                .map(|(_, best_span)| span < *best_span)
                .unwrap_or(true)
            {
                best = Some((path, span));
            }
        }

        let (path, _) = best.ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidInput,
                format!("render bundle has no path containing source range '{source_range}'"),
            )
        })?;
        let (local_start, local_end) = match path.source_interval.strand {
            '+' => (
                source_start - path.source_interval.start,
                source_end - path.source_interval.start,
            ),
            '-' => (
                path.source_interval.end - source_end,
                path.source_interval.end - source_start,
            ),
            strand => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("rendered path has invalid strand '{strand}'"),
                ));
            }
        };
        Ok(format!(
            "{}:{}-{}",
            path.rendered_name, local_start, local_end
        ))
    }
}

fn parse_source_range(range: &str) -> io::Result<(String, u64, u64)> {
    let (name, coords) = range.rsplit_once(':').ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid source range '{range}'; expected NAME:START-END"),
        )
    })?;
    let (start, end) = coords.split_once('-').ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid source range '{range}'; expected NAME:START-END"),
        )
    })?;
    let start = start.parse::<u64>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid start coordinate in source range '{range}'"),
        )
    })?;
    let end = end.parse::<u64>().map_err(|_| {
        io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid end coordinate in source range '{range}'"),
        )
    })?;
    if end <= start {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!("invalid source range '{range}': end must be greater than start"),
        ));
    }
    Ok((name.to_string(), start, end))
}

pub fn write_manifest(path: &Path, manifest: &RenderManifest) -> io::Result<()> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, manifest).map_err(io::Error::other)
}

pub fn read_manifest(path: &Path) -> io::Result<RenderManifest> {
    let file = File::open(path)?;
    let manifest: RenderManifest = serde_json::from_reader(file).map_err(io::Error::other)?;
    validate_manifest(&manifest)?;
    Ok(manifest)
}

pub fn load_bundle(root: impl AsRef<Path>) -> io::Result<LoadedRenderBundle> {
    let root = root.as_ref().to_path_buf();
    let manifest = read_manifest(&root.join("manifest.json"))?;
    let translation = resolve_bundle_path(&root, &manifest.translation);
    let tables = read_translation_binary(&translation)?;
    Ok(LoadedRenderBundle {
        root,
        manifest,
        tables,
    })
}

pub fn resolve_bundle_path(root: &Path, relative_or_absolute: &str) -> PathBuf {
    let path = Path::new(relative_or_absolute);
    if path.is_absolute() {
        path.to_path_buf()
    } else {
        root.join(path)
    }
}

fn validate_manifest(manifest: &RenderManifest) -> io::Result<()> {
    if manifest.format != RENDER_FORMAT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "render manifest has format '{}', expected '{}'",
                manifest.format, RENDER_FORMAT
            ),
        ));
    }
    if manifest.version != RENDER_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "render manifest has version {}, expected {}",
                manifest.version, RENDER_VERSION
            ),
        ));
    }
    if manifest.namespace.is_empty() || manifest.translation.is_empty() {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            "render manifest is missing namespace or translation path",
        ));
    }
    Ok(())
}

pub fn write_namespace(path: &Path, namespace: &SequenceNamespace) -> io::Result<()> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, namespace).map_err(io::Error::other)
}

pub fn write_translation_binary(path: &Path, tables: &RenderTranslationTables) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writer.write_all(TRANSLATION_MAGIC)?;
    bincode::serde::encode_into_std_write(tables, &mut writer, bincode::config::standard())
        .map_err(io::Error::other)?;
    writer.flush()
}

pub fn read_translation_binary(path: &Path) -> io::Result<RenderTranslationTables> {
    let mut reader = std::io::BufReader::new(File::open(path)?);
    let mut magic = [0u8; 8];
    reader.read_exact(&mut magic)?;
    if &magic != TRANSLATION_MAGIC {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "render translation '{}' has invalid magic bytes",
                path.display()
            ),
        ));
    }
    bincode::serde::decode_from_std_read(&mut reader, bincode::config::standard())
        .map_err(io::Error::other)
}

pub fn write_translation_tsv(path: &Path, tables: &RenderTranslationTables) -> io::Result<()> {
    let mut writer = BufWriter::new(File::create(path)?);
    writeln!(
        writer,
        "#record_type\trendered_path_id\trendered_name\tsource_sequence_id\tsource_name\tstart\tend\tstrand\tstep\tfeature_id\torientation"
    )?;
    for rendered in &tables.rendered_paths {
        let source = tables
            .namespace
            .get(rendered.source_interval.source_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "rendered path {} references missing source sequence {}",
                        rendered.rendered_path_id, rendered.source_interval.source_sequence_id.0
                    ),
                )
            })?;
        writeln!(
            writer,
            "path\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t.\t.\t.",
            rendered.rendered_path_id,
            rendered.rendered_name,
            rendered.source_interval.source_sequence_id.0,
            source.name,
            rendered.source_interval.start,
            rendered.source_interval.end,
            rendered.source_interval.strand
        )?;
    }
    for step in &tables.step_samples {
        let rendered = tables
            .rendered_paths
            .get(step.rendered_path_id as usize)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "step references missing rendered path {}",
                        step.rendered_path_id
                    ),
                )
            })?;
        let source = tables
            .namespace
            .get(rendered.source_interval.source_sequence_id)
            .ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "step references missing source sequence {}",
                        rendered.source_interval.source_sequence_id.0
                    ),
                )
            })?;
        writeln!(
            writer,
            "step\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            step.rendered_path_id,
            rendered.rendered_name,
            rendered.source_interval.source_sequence_id.0,
            source.name,
            step.source_bp,
            step.source_bp,
            rendered.source_interval.strand,
            step.rendered_step,
            step.feature_id,
            step.orientation
        )?;
    }
    writer.flush()
}

fn relative_file(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::sequence_namespace::{SequenceNamespace, SourceInterval};

    #[test]
    fn translation_binary_roundtrips() {
        let dir = std::env::temp_dir().join("impg_test_render_translation_roundtrip");
        std::fs::create_dir_all(&dir).unwrap();
        let path = dir.join("translation.bin");
        let mut namespace = SequenceNamespace::new();
        let source_id = namespace.add_sequence("HG002#1#chr6", 1000);
        let tables = RenderTranslationTables {
            namespace,
            rendered_paths: vec![RenderedPathRecord {
                rendered_path_id: 0,
                rendered_name: "HG002#1#chr6:10-100".to_string(),
                source_interval: SourceInterval::new(source_id, 10, 100, '+').unwrap(),
                gbwt_path_id: Some(0),
            }],
            step_samples: vec![StepTranslationRecord {
                rendered_path_id: 0,
                rendered_step: 0,
                source_bp: 10,
                feature_id: 42,
                orientation: '+',
            }],
        };

        write_translation_binary(&path, &tables).unwrap();
        let decoded = read_translation_binary(&path).unwrap();

        assert_eq!(
            decoded.rendered_paths[0].rendered_name,
            "HG002#1#chr6:10-100"
        );
        assert_eq!(decoded.step_samples[0].feature_id, 42);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn bundle_loader_validates_manifest_and_translation() {
        let dir = std::env::temp_dir().join("impg_test_render_bundle_loader");
        std::fs::remove_dir_all(&dir).ok();
        std::fs::create_dir_all(&dir).unwrap();
        let paths = RenderBundlePaths::new(&dir);

        let mut namespace = SequenceNamespace::new();
        let source_id = namespace.add_sequence("HG002#1#chr6", 1000);
        let tables = RenderTranslationTables {
            namespace,
            rendered_paths: vec![RenderedPathRecord {
                rendered_path_id: 0,
                rendered_name: "HG002#1#chr6:10-100".to_string(),
                source_interval: SourceInterval::new(source_id, 10, 100, '+').unwrap(),
                gbwt_path_id: Some(0),
            }],
            step_samples: vec![StepTranslationRecord {
                rendered_path_id: 0,
                rendered_step: 0,
                source_bp: 10,
                feature_id: 42,
                orientation: '+',
            }],
        };
        let manifest = RenderManifest::new_local_graph(
            "poa".to_string(),
            "panel.syng".to_string(),
            "HG002#1#chr6:10-100".to_string(),
            &paths,
            tables.namespace.sequences.len(),
            tables.rendered_paths.len(),
            tables.step_samples.len(),
        );

        write_manifest(&paths.manifest, &manifest).unwrap();
        write_translation_binary(&paths.translation, &tables).unwrap();
        let loaded = load_bundle(&dir).unwrap();

        assert_eq!(loaded.manifest.engine, "poa");
        assert_eq!(
            loaded.tables.rendered_paths[0].rendered_name,
            "HG002#1#chr6:10-100"
        );
        assert_eq!(loaded.tables.step_samples[0].source_bp, 10);
        std::fs::remove_dir_all(&dir).ok();
    }

    #[test]
    fn rendered_target_range_projects_source_coordinates() {
        let mut namespace = SequenceNamespace::new();
        let source_id = namespace.add_sequence("HG002#1#chr6", 1000);
        let root = PathBuf::from("/tmp/render-test");
        let paths = RenderBundlePaths::new(&root);
        let tables = RenderTranslationTables {
            namespace,
            rendered_paths: vec![
                RenderedPathRecord {
                    rendered_path_id: 0,
                    rendered_name: "HG002#1#chr6:100-500(+)".to_string(),
                    source_interval: SourceInterval::new(source_id, 100, 500, '+').unwrap(),
                    gbwt_path_id: Some(0),
                },
                RenderedPathRecord {
                    rendered_path_id: 1,
                    rendered_name: "HG002#1#chr6:100-500(-)".to_string(),
                    source_interval: SourceInterval::new(source_id, 100, 500, '-').unwrap(),
                    gbwt_path_id: Some(1),
                },
            ],
            step_samples: Vec::new(),
        };
        let bundle = LoadedRenderBundle {
            root,
            manifest: RenderManifest::new_syng_native(
                "panel.syng".to_string(),
                "HG002#1#chr6:150-250".to_string(),
                &paths,
                true,
                tables.namespace.sequences.len(),
                tables.rendered_paths.len(),
                tables.step_samples.len(),
            ),
            tables,
        };

        assert_eq!(
            bundle.rendered_target_range(None).unwrap(),
            "HG002#1#chr6:100-500(+):50-150"
        );
        assert_eq!(
            bundle
                .rendered_target_range(Some("HG002#1#chr6:200-300"))
                .unwrap(),
            "HG002#1#chr6:100-500(+):100-200"
        );
    }
}
