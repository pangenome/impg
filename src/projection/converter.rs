use rustc_hash::FxHashMap;
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};

use crate::commands::genotype::{parse_normalized_gfa, GraphContributionModel, GraphFeatureIdMode};
use crate::genotyping::FeatureSpace;

use super::{PROJECTION_FORMAT, PROJECTION_VERSION};

pub const GFA_PROJECTION_METHOD: &str = "gaf-to-gfa";
pub const DEFAULT_PACK_NAME: &str = "sample.pack.tsv";
pub const DEFAULT_GAF_NAME: &str = "alignments.gaf";
pub const DEFAULT_READ_CONTRIBUTIONS_NAME: &str = "read-contributions.tsv";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GfaProjectionOutputFormat {
    ProjectionBundle,
    PackTsv,
}

impl GfaProjectionOutputFormat {
    pub fn parse(raw: &str) -> io::Result<Self> {
        match raw {
            "proj" | "projection" | "projection-bundle" | "bundle" => Ok(Self::ProjectionBundle),
            "pack" | "pack-tsv" | "pack-text" | "packtsv" => Ok(Self::PackTsv),
            other => Err(io::Error::new(
                io::ErrorKind::InvalidInput,
                format!(
                    "unsupported project output format '{}'; expected 'proj' or 'pack-tsv'",
                    other
                ),
            )),
        }
    }

    pub fn as_str(self) -> &'static str {
        match self {
            Self::ProjectionBundle => "proj",
            Self::PackTsv => "pack-tsv",
        }
    }
}

pub struct GfaProjectionConfig<'a> {
    pub gfa_path: &'a Path,
    pub gaf_path: &'a Path,
    pub output_path: &'a Path,
    pub output_format: GfaProjectionOutputFormat,
    pub feature_id_mode: GraphFeatureIdMode,
    pub contribution_model: GraphContributionModel,
    pub read_contributions_path: Option<&'a Path>,
}

#[derive(Debug, Clone)]
pub struct GfaProjectionSummary {
    pub output_path: PathBuf,
    pub pack_path: PathBuf,
    pub read_contributions_path: Option<PathBuf>,
    pub feature_space: String,
    pub graph_id: String,
    pub feature_id_mode: String,
    pub contribution_model: String,
    pub total_records: u64,
    pub retained_records: u64,
    pub contributed_steps: u64,
    pub nonzero_features: usize,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct GfaProjectionManifest {
    pub format: String,
    pub version: u32,
    pub projection_method: String,
    pub graph: String,
    pub gaf: String,
    pub pack: String,
    pub read_contributions: Option<String>,
    pub feature_space: String,
    pub graph_id: String,
    pub feature_id_mode: String,
    pub contribution_model: String,
    pub read_space: String,
}

#[derive(Debug, Clone)]
pub struct GfaProjectionPaths {
    pub root: PathBuf,
    pub graph_path: PathBuf,
    pub gaf_path: PathBuf,
    pub pack_path: PathBuf,
    pub read_contributions_path: Option<PathBuf>,
    pub feature_space: String,
    pub graph_id: String,
    pub feature_id_mode: String,
    pub contribution_model: String,
}

#[derive(Debug, Clone)]
struct SegmentInfo {
    feature_id: u32,
    length: u64,
}

#[derive(Debug, Clone)]
struct GafWalkStep {
    segment_name: String,
    orientation: char,
}

#[derive(Debug, Clone)]
struct ReadContribution {
    read_name: String,
    read_ordinal: u64,
    step_index: usize,
    segment_name: String,
    orientation: char,
    feature_id: u32,
    segment_visit_in_read: u64,
    count_delta: u64,
    explanation: String,
}

#[derive(Debug, Default)]
struct ProjectionCounts {
    counts: FxHashMap<u32, u64>,
    contributions: Vec<ReadContribution>,
    total_records: u64,
    retained_records: u64,
    contributed_steps: u64,
}

fn create_parent_dir(path: &Path) -> io::Result<()> {
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent)?;
        }
    }
    Ok(())
}

fn resolve_relative(root: &Path, path: &str) -> PathBuf {
    let p = Path::new(path);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        root.join(p)
    }
}

fn path_string(path: &Path) -> String {
    path.to_string_lossy().to_string()
}

fn parse_u64_field(value: &str, field_name: &str, line_no: usize) -> io::Result<u64> {
    value.parse::<u64>().map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid GAF {} on line {}: {}", field_name, line_no, e),
        )
    })
}

fn parse_p_token(token: &str, line_no: usize) -> io::Result<GafWalkStep> {
    if let Some(name) = token.strip_suffix('+') {
        if !name.is_empty() {
            return Ok(GafWalkStep {
                segment_name: name.to_string(),
                orientation: '+',
            });
        }
    }
    if let Some(name) = token.strip_suffix('-') {
        if !name.is_empty() {
            return Ok(GafWalkStep {
                segment_name: name.to_string(),
                orientation: '-',
            });
        }
    }
    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "GAF line {} has invalid comma-walk token '{}'",
            line_no, token
        ),
    ))
}

fn parse_gaf_walk(path_field: &str, line_no: usize) -> io::Result<Vec<GafWalkStep>> {
    if path_field.is_empty() || path_field == "*" {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("GAF line {} does not contain a graph walk path", line_no),
        ));
    }

    if path_field.starts_with('>') || path_field.starts_with('<') {
        let mut steps = Vec::new();
        let mut orientation: Option<char> = None;
        let mut name = String::new();
        for ch in path_field.chars() {
            if ch == '>' || ch == '<' {
                if let Some(prev) = orientation {
                    if name.is_empty() {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!("GAF line {} has an empty graph-walk step", line_no),
                        ));
                    }
                    steps.push(GafWalkStep {
                        segment_name: std::mem::take(&mut name),
                        orientation: prev,
                    });
                }
                orientation = Some(if ch == '>' { '+' } else { '-' });
            } else {
                if orientation.is_none() {
                    return Err(io::Error::new(
                        io::ErrorKind::InvalidData,
                        format!("GAF line {} graph walk must start with '>' or '<'", line_no),
                    ));
                }
                name.push(ch);
            }
        }
        if let Some(prev) = orientation {
            if name.is_empty() {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GAF line {} has an empty trailing graph-walk step", line_no),
                ));
            }
            steps.push(GafWalkStep {
                segment_name: name,
                orientation: prev,
            });
        }
        if steps.is_empty() {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("GAF line {} has an empty graph walk", line_no),
            ));
        }
        return Ok(steps);
    }

    if path_field.contains(',') {
        let steps = path_field
            .split(',')
            .map(|token| parse_p_token(token, line_no))
            .collect::<io::Result<Vec<_>>>()?;
        if !steps.is_empty() {
            return Ok(steps);
        }
    }

    Err(io::Error::new(
        io::ErrorKind::InvalidData,
        format!(
            "GAF line {} path field '{}' is not an oriented graph walk over GFA segment names",
            line_no, path_field
        ),
    ))
}

fn project_gaf_records(
    gaf_path: &Path,
    segment_by_name: &FxHashMap<String, SegmentInfo>,
) -> io::Result<ProjectionCounts> {
    let file = File::open(gaf_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("failed to open GAF '{}': {}", gaf_path.display(), e),
        )
    })?;
    let mut projection = ProjectionCounts::default();

    for (line_idx, line) in BufReader::new(file).lines().enumerate() {
        let line_no = line_idx + 1;
        let line = line?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() < 12 {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GAF line {} has {} field(s); expected at least 12",
                    line_no,
                    fields.len()
                ),
            ));
        }

        projection.total_records += 1;
        let read_ordinal = projection.total_records;
        let read_name = fields[0];
        let path_field = fields[5];
        let path_start = parse_u64_field(fields[7], "path_start", line_no)?;
        let path_end = parse_u64_field(fields[8], "path_end", line_no)?;
        if path_end <= path_start {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GAF line {} has empty path interval {}-{}",
                    line_no, path_start, path_end
                ),
            ));
        }

        let walk = parse_gaf_walk(path_field, line_no)?;
        let mut cursor = 0u64;
        let mut retained = false;
        let mut segment_visits: FxHashMap<String, u64> = FxHashMap::default();
        for (step_offset, step) in walk.iter().enumerate() {
            let segment = segment_by_name.get(&step.segment_name).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "GAF line {} references unknown GFA segment '{}'",
                        line_no, step.segment_name
                    ),
                )
            })?;
            let step_start = cursor;
            let step_end = step_start.checked_add(segment.length).ok_or_else(|| {
                io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("GAF line {} walk length overflows u64", line_no),
                )
            })?;
            cursor = step_end;
            let overlap_start = step_start.max(path_start);
            let overlap_end = step_end.min(path_end);
            if overlap_end <= overlap_start {
                continue;
            }

            retained = true;
            projection.contributed_steps += 1;
            *projection.counts.entry(segment.feature_id).or_insert(0) += 1;
            let visit = segment_visits.entry(step.segment_name.clone()).or_insert(0);
            *visit += 1;
            let explanation = if *visit == 1 {
                "first visit to segment in read; counted".to_string()
            } else {
                format!(
                    "repeated visit {} to segment in read; counted again",
                    *visit
                )
            };
            projection.contributions.push(ReadContribution {
                read_name: read_name.to_string(),
                read_ordinal,
                step_index: step_offset + 1,
                segment_name: step.segment_name.clone(),
                orientation: step.orientation,
                feature_id: segment.feature_id,
                segment_visit_in_read: *visit,
                count_delta: 1,
                explanation,
            });
        }

        if path_end > cursor {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!(
                    "GAF line {} path_end {} exceeds graph-walk length {}",
                    line_no, path_end, cursor
                ),
            ));
        }
        if retained {
            projection.retained_records += 1;
        }
    }

    Ok(projection)
}

fn sorted_counts(counts: &FxHashMap<u32, u64>) -> BTreeMap<u32, u64> {
    counts
        .iter()
        .filter_map(|(&feature_id, &count)| (count > 0).then_some((feature_id, count)))
        .collect()
}

fn write_typed_pack_tsv(
    path: &Path,
    counts: &FxHashMap<u32, u64>,
    graph_id: &str,
    feature_id_mode: &str,
    contribution_model: GraphContributionModel,
) -> io::Result<usize> {
    create_parent_dir(path)?;
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(out, "#feature_space\t{}", FeatureSpace::GfaSegment.as_str())?;
    writeln!(out, "#graph_id\t{}", graph_id)?;
    writeln!(out, "#feature_id_mode\t{}", feature_id_mode)?;
    writeln!(
        out,
        "#graph_contribution_model\t{}",
        contribution_model.as_str()
    )?;
    writeln!(out, "#projection_method\t{}", GFA_PROJECTION_METHOD)?;
    writeln!(out, "#node_id\tcount")?;
    let sorted = sorted_counts(counts);
    for (feature_id, count) in &sorted {
        writeln!(out, "{feature_id}\t{count}")?;
    }
    out.flush()?;
    Ok(sorted.len())
}

fn write_read_contributions(path: &Path, rows: &[ReadContribution]) -> io::Result<()> {
    create_parent_dir(path)?;
    let mut out = BufWriter::new(File::create(path)?);
    writeln!(
        out,
        "read_name\tread_ordinal\tstep_index\tsegment_name\torientation\tfeature_id\tsegment_visit_in_read\tcount_delta\texplanation"
    )?;
    for row in rows {
        writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            row.read_name,
            row.read_ordinal,
            row.step_index,
            row.segment_name,
            row.orientation,
            row.feature_id,
            row.segment_visit_in_read,
            row.count_delta,
            row.explanation
        )?;
    }
    out.flush()
}

fn write_manifest(path: &Path, manifest: &GfaProjectionManifest) -> io::Result<()> {
    create_parent_dir(path)?;
    let file = File::create(path)?;
    serde_json::to_writer_pretty(file, manifest).map_err(io::Error::other)
}

pub fn load_gfa_projection_bundle(path: impl AsRef<Path>) -> io::Result<GfaProjectionPaths> {
    let root = path.as_ref();
    let manifest_path = root.join("manifest.json");
    let file = File::open(&manifest_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "failed to open projection manifest '{}': {}",
                manifest_path.display(),
                e
            ),
        )
    })?;
    let manifest: GfaProjectionManifest = serde_json::from_reader(file).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse GFA projection manifest '{}': {}",
                manifest_path.display(),
                e
            ),
        )
    })?;
    if manifest.format != PROJECTION_FORMAT {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "projection manifest has format '{}', expected '{}'",
                manifest.format, PROJECTION_FORMAT
            ),
        ));
    }
    if manifest.version != PROJECTION_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "unsupported projection version {}; expected {}",
                manifest.version, PROJECTION_VERSION
            ),
        ));
    }
    if manifest.projection_method != GFA_PROJECTION_METHOD {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "projection method '{}' is not a GFA GAF projection",
                manifest.projection_method
            ),
        ));
    }
    if !matches!(
        manifest.feature_space.as_str(),
        "gfa-segment" | "variation-graph-node"
    ) {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "projection feature_space '{}' is not a graph-node feature space",
                manifest.feature_space
            ),
        ));
    }

    let parent = manifest_path.parent().unwrap_or(root);
    Ok(GfaProjectionPaths {
        root: root.to_path_buf(),
        graph_path: resolve_relative(parent, &manifest.graph),
        gaf_path: resolve_relative(parent, &manifest.gaf),
        pack_path: resolve_relative(parent, &manifest.pack),
        read_contributions_path: manifest
            .read_contributions
            .as_deref()
            .map(|p| resolve_relative(parent, p)),
        feature_space: manifest.feature_space,
        graph_id: manifest.graph_id,
        feature_id_mode: manifest.feature_id_mode,
        contribution_model: manifest.contribution_model,
    })
}

pub fn project_gaf_to_gfa(config: &GfaProjectionConfig<'_>) -> io::Result<GfaProjectionSummary> {
    let gfa_text = fs::read_to_string(config.gfa_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!("failed to read GFA '{}': {}", config.gfa_path.display(), e),
        )
    })?;
    let graph = parse_normalized_gfa(
        &gfa_text,
        "gfa",
        Some(config.gfa_path.to_path_buf()),
        None,
        FeatureSpace::GfaSegment.as_str(),
        config.feature_id_mode,
    )?;

    let mut segment_by_name = FxHashMap::default();
    for segment in &graph.segments {
        segment_by_name.insert(
            segment.name.clone(),
            SegmentInfo {
                feature_id: segment.feature_id,
                length: segment.length,
            },
        );
    }

    let projection = project_gaf_records(config.gaf_path, &segment_by_name)?;
    let feature_id_mode = graph.effective_feature_id_mode.as_str().to_string();
    let contribution_model = config.contribution_model.as_str().to_string();

    match config.output_format {
        GfaProjectionOutputFormat::ProjectionBundle => {
            fs::create_dir_all(config.output_path)?;
            let pack_path = config.output_path.join(DEFAULT_PACK_NAME);
            let debug_path = config.output_path.join(DEFAULT_READ_CONTRIBUTIONS_NAME);
            let gaf_copy_path = config.output_path.join(DEFAULT_GAF_NAME);
            let nonzero_features = write_typed_pack_tsv(
                &pack_path,
                &projection.counts,
                &graph.graph_id,
                &feature_id_mode,
                config.contribution_model,
            )?;
            write_read_contributions(&debug_path, &projection.contributions)?;
            fs::copy(config.gaf_path, &gaf_copy_path).map_err(|e| {
                io::Error::new(
                    e.kind(),
                    format!(
                        "failed to copy GAF '{}' to '{}': {}",
                        config.gaf_path.display(),
                        gaf_copy_path.display(),
                        e
                    ),
                )
            })?;
            let manifest = GfaProjectionManifest {
                format: PROJECTION_FORMAT.to_string(),
                version: PROJECTION_VERSION,
                projection_method: GFA_PROJECTION_METHOD.to_string(),
                graph: path_string(config.gfa_path),
                gaf: DEFAULT_GAF_NAME.to_string(),
                pack: DEFAULT_PACK_NAME.to_string(),
                read_contributions: Some(DEFAULT_READ_CONTRIBUTIONS_NAME.to_string()),
                feature_space: FeatureSpace::GfaSegment.as_str().to_string(),
                graph_id: graph.graph_id.clone(),
                feature_id_mode: feature_id_mode.clone(),
                contribution_model: contribution_model.clone(),
                read_space: "gaf-graph-walk".to_string(),
            };
            write_manifest(&config.output_path.join("manifest.json"), &manifest)?;
            Ok(GfaProjectionSummary {
                output_path: config.output_path.to_path_buf(),
                pack_path,
                read_contributions_path: Some(debug_path),
                feature_space: FeatureSpace::GfaSegment.as_str().to_string(),
                graph_id: graph.graph_id,
                feature_id_mode,
                contribution_model,
                total_records: projection.total_records,
                retained_records: projection.retained_records,
                contributed_steps: projection.contributed_steps,
                nonzero_features,
            })
        }
        GfaProjectionOutputFormat::PackTsv => {
            let nonzero_features = write_typed_pack_tsv(
                config.output_path,
                &projection.counts,
                &graph.graph_id,
                &feature_id_mode,
                config.contribution_model,
            )?;
            let debug_path = if let Some(path) = config.read_contributions_path {
                write_read_contributions(path, &projection.contributions)?;
                Some(path.to_path_buf())
            } else {
                None
            };
            Ok(GfaProjectionSummary {
                output_path: config.output_path.to_path_buf(),
                pack_path: config.output_path.to_path_buf(),
                read_contributions_path: debug_path,
                feature_space: FeatureSpace::GfaSegment.as_str().to_string(),
                graph_id: graph.graph_id,
                feature_id_mode,
                contribution_model,
                total_records: projection.total_records,
                retained_records: projection.retained_records,
                contributed_steps: projection.contributed_steps,
                nonzero_features,
            })
        }
    }
}
