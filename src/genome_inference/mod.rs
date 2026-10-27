//! Experimental ownership diagnostics, quantitative haploid calls and source-path threads.
pub mod calling;
pub mod catalog;
pub mod diplotype_chain;
pub mod diplotype_dosage;
pub mod evidence_store;
pub mod genotype;
pub mod mem_records;
pub mod joint;
pub mod observations;
pub mod panel_routes;
pub mod realign_score;
pub mod reanchor;
pub mod reconstruction;
pub mod sample;
pub mod sample_projection;
pub mod sequence_evaluation;
pub mod territory_frame;
pub mod threading;

use crate::sample_mem_bwt::invalid;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::fs::{self, File, OpenOptions};
use std::io::{self, Read, Write};
use std::path::Path;

pub const FORMAT_VERSION: u32 = 1;
pub const MODEL: &str = "presence-compatibility";
pub const COUNT_POLICY: &str = "alternating-node-gap-v1;best-query-orientation-both-inputs;content-subwalk-maximal;canonical-records;rc-orbit-occurrences";

/// Stable FNV-1a content checksum. Detects accidental corruption/identity changes;
/// not an adversarial authenticity or cryptographic digest contract.
pub fn checksum(bytes: &[u8]) -> u64 {
    let mut hash = 0xcbf29ce484222325;
    hash_update(&mut hash, bytes);
    hash
}
fn hash_update(hash: &mut u64, bytes: &[u8]) {
    for b in bytes {
        *hash = (*hash ^ *b as u64).wrapping_mul(0x100000001b3);
    }
}

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub struct PanelIdentity {
    pub checksum_algorithm: String,
    pub sidecars: Vec<(String, u64, String)>,
}
impl PanelIdentity {
    /// Includes native graph and dictionary, names, parameters, and positional
    /// sidecars. Same parameters alone cannot make two panels compatible.
    pub fn read(prefix: &str) -> io::Result<Self> {
        let mut sidecars = Vec::new();
        for suffix in ["1gbwt", "1khash", "names", "meta", "spos", "pstep"] {
            let primary = if suffix.starts_with('1') || prefix.ends_with(".syng") {
                format!("{prefix}.{suffix}")
            } else {
                format!("{prefix}.syng.{suffix}")
            };
            // Match SyngIndex's primary/legacy resolution exactly.
            let path = if Path::new(&primary).exists() {
                primary
            } else {
                format!("{prefix}.syng.{suffix}")
            };
            let mut file = File::open(&path)?;
            let mut hash = 0xcbf29ce484222325;
            let mut size = 0;
            let mut buffer = vec![0; 1024 * 1024];
            loop {
                let n = file.read(&mut buffer)?;
                if n == 0 {
                    break;
                }
                size += n as u64;
                hash_update(&mut hash, &buffer[..n]);
            }
            sidecars.push((suffix.to_string(), size, format!("{hash:016x}")));
        }
        // Legacy ambiguity-derived dictionaries are not acceptable for new samples.
        let primary_meta = crate::syng::syng_meta_path(prefix);
        let meta_path = if Path::new(&primary_meta).exists() {
            primary_meta
        } else {
            format!("{prefix}.syng.meta")
        };
        let meta_text = fs::read_to_string(&meta_path)?;
        let metadata_version: u32 = meta_text
            .lines()
            .find_map(|l| l.strip_prefix("version\t").and_then(|v| v.parse().ok()))
            .unwrap_or(0);
        if metadata_version != 2 {
            // THE LEGACY-PANEL OPT-IN (the global-panel call mode's
            // committed-evidence continuity, fail-closed by default):
            // the v2 gate exists so legacy ambiguity-derived
            // dictionaries are never used to derive NEW samples'
            // evidence. The ONE committed global panel (the May
            // projection build every committed cut/call receipt is
            // keyed to) is metadata v1 — and byte-identity against
            // those receipts is unattainable against any OTHER
            // dictionary (a fresh rebuild renumbers every node id,
            // and the projection's records ARE this panel's ids). The
            // opt-in states that continuity EXPLICITLY; it never
            // becomes the default, and every use receipts itself
            // loudly. The honest retirement path (a metadata
            // re-stamp of the SAME dictionary, no id churn) is
            // named in the design doc.
            if std::env::var("IMPG_PANEL_LEGACY_V1_OK").ok().as_deref() != Some("1") {
                return Err(invalid(
                    "genome-infer requires freshly built metadata-v2 panel \
                     (this one is metadata v1; set IMPG_PANEL_LEGACY_V1_OK=1 ONLY for \
                     committed-evidence continuity against this panel's own receipts)",
                ));
            }
            let built = fs::metadata(&meta_path)
                .and_then(|m| m.modified())
                .map(|t| format!("{}", t.duration_since(std::time::UNIX_EPOCH).map(|d| d.as_secs()).unwrap_or(0)))
                .unwrap_or_else(|_| "<unrecorded>".into());
            eprintln!(
                "[panel] LEGACY-PANEL OPT-IN: {prefix} is metadata v{metadata_version} \
                 (meta dated unix {built}); accepted ONLY for identity-gate continuity against the \
                 committed evidence base keyed to this panel's own node ids"
            );
        }
        Ok(Self {
            checksum_algorithm: "fnv1a64-content-v1".into(),
            sidecars,
        })
    }
}

/// The panel's own metadata version statement (the syng `.meta`
/// `version` line; 0 when the meta carries none) — a provenance
/// receipt helper, NOT a gate (the gate lives in
/// `PanelIdentity::read`, fail-closed by default).
pub fn panel_metadata_version(prefix: &str) -> io::Result<u32> {
    let primary_meta = crate::syng::syng_meta_path(prefix);
    let meta_path = if Path::new(&primary_meta).exists() {
        primary_meta
    } else {
        format!("{prefix}.syng.meta")
    };
    Ok(fs::read_to_string(meta_path)?
        .lines()
        .find_map(|l| l.strip_prefix("version\t").and_then(|v| v.parse().ok()))
        .unwrap_or(0))
}

pub fn write_json(path: &Path, value: &impl Serialize) -> io::Result<()> {
    atomic_write(
        path,
        &serde_json::to_vec_pretty(value).map_err(io::Error::other)?,
    )
}
pub fn read_json<T: DeserializeOwned>(path: &Path) -> io::Result<T> {
    serde_json::from_reader(io::BufReader::new(File::open(path)?)).map_err(io::Error::other)
}
pub fn atomic_write(path: &Path, bytes: &[u8]) -> io::Result<()> {
    let tmp = path.with_extension("incomplete");
    let mut f = OpenOptions::new().write(true).create_new(true).open(&tmp)?;
    f.write_all(bytes)?;
    f.sync_all()?;
    fs::rename(tmp, path)
}

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Envelope<T> {
    version: u32,
    checksum: String,
    payload: T,
}
pub fn save_catalog(path: &Path, catalog: &catalog::Catalog) -> io::Result<()> {
    let bytes = serde_json::to_vec(catalog).map_err(io::Error::other)?;
    write_json(
        path,
        &Envelope {
            version: FORMAT_VERSION,
            checksum: format!("{:016x}", checksum(&bytes)),
            payload: catalog,
        },
    )
}
pub fn load_catalog(path: &Path) -> io::Result<catalog::Catalog> {
    load_catalog_with_checksum(path).map(|(catalog, _)| catalog)
}
pub fn load_catalog_with_checksum(path: &Path) -> io::Result<(catalog::Catalog, String)> {
    let envelope: Envelope<catalog::Catalog> = read_json(path)?;
    let bytes = serde_json::to_vec(&envelope.payload).map_err(io::Error::other)?;
    if envelope.version != FORMAT_VERSION
        || envelope.checksum != format!("{:016x}", checksum(&bytes))
    {
        return Err(invalid("catalog version/checksum mismatch"));
    }
    envelope.payload.validate()?;
    Ok((envelope.payload, envelope.checksum))
}

/// Reserve a NEW directory atomically; even empty preexisting directories are
/// refused. Status is authoritative: absent/running/failed never means success.
pub fn with_output(
    out: &Path,
    operation: impl FnOnce() -> io::Result<serde_json::Value>,
) -> io::Result<()> {
    with_output_model(out, MODEL, operation)
}
pub fn with_output_model(
    out: &Path,
    model: &str,
    operation: impl FnOnce() -> io::Result<serde_json::Value>,
) -> io::Result<()> {
    fs::create_dir(out)?;
    write_json(
        &out.join("manifest.json"),
        &serde_json::json!({"version": FORMAT_VERSION, "status": "running", "experimental": true, "model": model}),
    )?;
    let started = std::time::Instant::now();
    match operation() {
        Ok(stats) => {
            let result = write_json(
                &out.join("manifest.json"),
                &serde_json::json!({
                "version": FORMAT_VERSION, "status": "succeeded", "experimental": true,
                "model": model, "elapsed_seconds": started.elapsed().as_secs_f64(), "stats": stats,
                "argv": std::env::args().collect::<Vec<_>>() }),
            );
            if result.is_err() {
                let _ = fs::remove_file(out.join("manifest.json"));
            }
            result
        }
        Err(error) => {
            // Intermediate reusable objects remain diagnostic only under failed status.
            let _ = fs::remove_file(out.join("calls.json"));
            let _ = fs::remove_file(out.join("evaluation.json"));
            let _ = fs::remove_file(out.join("threads.json"));
            for name in [
                "reconstruction.fa",
                "provenance.json",
                "unresolved.json",
                "unresolved.bed",
                "evaluation.json",
            ] {
                let _ = fs::remove_file(out.join(name));
                let _ = fs::remove_file(out.join(name).with_extension("incomplete"));
            }
            let _ = write_json(
                &out.join("manifest.json"),
                &serde_json::json!({"version": FORMAT_VERSION,
                "status": "failed", "experimental": true, "model": model, "error": error.to_string()}),
            );
            Err(error)
        }
    }
}
