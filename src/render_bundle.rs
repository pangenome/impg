use crate::sequence_namespace::{SequenceNamespace, SourceInterval};
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{self, BufWriter, Write};
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
    pub graph_gfa: PathBuf,
    pub syng_prefix: PathBuf,
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
            graph_gfa: include_gfa.then(|| relative_file(&paths.root, &paths.graph_gfa)),
            syng_prefix: Some(relative_file(&paths.root, &paths.syng_prefix)),
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
            graph_gfa: root.join("graph.gfa"),
            syng_prefix: root.join("paths"),
            root,
        }
    }
}

pub fn write_manifest(path: &Path, manifest: &RenderManifest) -> io::Result<()> {
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, manifest).map_err(io::Error::other)
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
    std::io::Read::read_exact(&mut reader, &mut magic)?;
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
}
