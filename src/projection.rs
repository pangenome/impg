use serde::{Deserialize, Serialize};
use std::fs;
use std::io;
use std::path::{Path, PathBuf};

pub mod converter;

pub const PROJECTION_FORMAT: &str = "impg-projection";
pub const PROJECTION_VERSION: u32 = 1;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProjectionManifest {
    pub format: String,
    pub version: u32,
    pub syng_prefix: String,
    pub pack: String,
    pub gaf: Option<String>,
    pub feature_space: String,
    pub read_space: Option<String>,
}

#[derive(Debug, Clone)]
pub struct ProjectionPaths {
    pub root: PathBuf,
    pub syng_prefix: String,
    pub pack_path: PathBuf,
    pub gaf_path: Option<PathBuf>,
}

fn resolve_relative(root: &Path, path: &str) -> PathBuf {
    let p = Path::new(path);
    if p.is_absolute() {
        p.to_path_buf()
    } else {
        root.join(p)
    }
}

pub fn write_manifest(
    root: &Path,
    syng_prefix: &str,
    pack: &str,
    gaf: Option<&str>,
) -> io::Result<()> {
    fs::create_dir_all(root)?;
    let manifest = ProjectionManifest {
        format: PROJECTION_FORMAT.to_string(),
        version: PROJECTION_VERSION,
        syng_prefix: syng_prefix.to_string(),
        pack: pack.to_string(),
        gaf: gaf.map(str::to_string),
        feature_space: "syng-syncmer-node".to_string(),
        read_space: gaf.map(|_| "syng-gaf-walk".to_string()),
    };
    let manifest_path = root.join("manifest.json");
    let file = fs::File::create(manifest_path)?;
    serde_json::to_writer_pretty(file, &manifest).map_err(io::Error::other)?;
    Ok(())
}

pub fn load(path: impl AsRef<Path>) -> io::Result<ProjectionPaths> {
    let root = path.as_ref();
    let manifest_path = root.join("manifest.json");
    let file = fs::File::open(&manifest_path).map_err(|e| {
        io::Error::new(
            e.kind(),
            format!(
                "failed to open projection manifest '{}': {}",
                manifest_path.display(),
                e
            ),
        )
    })?;
    let manifest: ProjectionManifest = serde_json::from_reader(file).map_err(|e| {
        io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "failed to parse projection manifest '{}': {}",
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
    let parent = manifest_path.parent().unwrap_or(root);
    let pack_path = resolve_relative(parent, &manifest.pack);
    let gaf_path = manifest
        .gaf
        .as_deref()
        .map(|gaf| resolve_relative(parent, gaf));
    Ok(ProjectionPaths {
        root: root.to_path_buf(),
        syng_prefix: manifest.syng_prefix,
        pack_path,
        gaf_path,
    })
}
