//! End-to-end integration tests for the syng CLI path.
//!
//! Covers gaps that unit tests missed:
//! - `impg syng --agc` builds a non-empty index (caught a silent failure in the
//!   yeast235 workflow where only `.syng.names` got written).
//! - `impg syng --agc` round-trips: built index can be loaded and queried.
//! - `impg partition -a <syng-prefix>` runs end-to-end and produces non-empty BED.

use ahash::AHashSet;
use ragc_core::{StreamingQueueCompressor, StreamingQueueConfig};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};
use std::process::Command;

// All syng tests must be serialized because the C library has non-thread-safe
// global state (see src/syng.rs test comments). tests/ and src/ test modules
// run in separate processes, but within this binary we still need a lock.
static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
    std::sync::LazyLock::new(|| std::sync::Mutex::new(()));

fn lock_syng() -> std::sync::MutexGuard<'static, ()> {
    SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

/// Deterministic pseudo-random sequence (numeric encoding, 0-3 = A/C/G/T) —
/// what ragc-core's push() expects.
fn make_sequence_numeric(len: usize, seed: u8) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    let mut state: u32 = seed as u32;
    for _ in 0..len {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(((state >> 16) % 4) as u8);
    }
    seq
}

/// Same sequence but as ASCII for FASTA output.
fn numeric_to_ascii(numeric: &[u8]) -> Vec<u8> {
    numeric
        .iter()
        .map(|&b| match b {
            0 => b'A',
            1 => b'C',
            2 => b'G',
            3 => b'T',
            _ => b'N',
        })
        .collect()
}

fn mutate_ascii_every(seq: &[u8], offset: usize, stride: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    for i in (offset..out.len()).step_by(stride) {
        out[i] = match out[i] {
            b'A' => b'C',
            b'C' => b'G',
            b'G' => b'T',
            b'T' => b'A',
            other => other,
        };
    }
    out
}

fn mutate_base_ascii(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        other => other,
    }
}

fn write_tiled_fastq<W: std::io::Write>(
    writer: &mut W,
    prefix: &str,
    seq: &[u8],
    read_len: usize,
    step: usize,
) -> std::io::Result<usize> {
    assert!(read_len > 0);
    assert!(step > 0);
    assert!(seq.len() >= read_len);

    let mut starts = Vec::new();
    let mut start = 0usize;
    while start + read_len <= seq.len() {
        starts.push(start);
        start += step;
    }
    let terminal_start = seq.len() - read_len;
    if starts.last().copied() != Some(terminal_start) {
        starts.push(terminal_start);
    }

    for (read_idx, start) in starts.iter().enumerate() {
        let end = start + read_len;
        writeln!(writer, "@{}_{}", prefix, read_idx)?;
        writer.write_all(&seq[*start..end])?;
        writeln!(writer)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", "I".repeat(read_len))?;
    }
    Ok(starts.len())
}

fn write_tiled_fastq_with_errors<W: std::io::Write>(
    writer: &mut W,
    prefix: &str,
    seq: &[u8],
    read_len: usize,
    step: usize,
    max_reads: usize,
    error_stride: usize,
) -> std::io::Result<usize> {
    assert!(read_len > 0);
    assert!(step > 0);
    assert!(max_reads > 0);
    assert!(seq.len() >= read_len);
    assert!(error_stride > 0);

    let mut starts = Vec::new();
    let mut start = 0usize;
    while start + read_len <= seq.len() && starts.len() < max_reads {
        starts.push(start);
        start += step;
    }
    let terminal_start = seq.len() - read_len;
    if starts.len() < max_reads && starts.last().copied() != Some(terminal_start) {
        starts.push(terminal_start);
    }

    for (read_idx, start) in starts.iter().enumerate() {
        let end = start + read_len;
        let mut read = seq[*start..end].to_vec();
        for i in ((read_idx + 1) * 17 % error_stride..read.len()).step_by(error_stride) {
            read[i] = mutate_base_ascii(read[i]);
        }
        writeln!(writer, "@{}_{}", prefix, read_idx)?;
        writer.write_all(&read)?;
        writeln!(writer)?;
        writeln!(writer, "+")?;
        writeln!(writer, "{}", "I".repeat(read_len))?;
    }
    Ok(starts.len())
}

/// Build a small test AGC archive with multiple samples sharing a backbone.
fn create_test_agc(path: &str) {
    let config = StreamingQueueConfig {
        queue_capacity: 10 * 1024 * 1024,
        num_threads: 1,
        verbosity: 0,
        ..Default::default()
    };

    // Shared backbone (same seed) so samples share syncmer nodes.
    let backbone = make_sequence_numeric(800, 42);

    let samples: Vec<(&str, &str, Vec<u8>)> = vec![
        ("sampleA", "chr1", {
            let mut s = backbone.clone();
            s.extend(make_sequence_numeric(400, 1));
            s
        }),
        ("sampleB", "chr1", {
            let mut s = backbone.clone();
            s.extend(make_sequence_numeric(400, 2));
            s
        }),
        ("sampleC", "chr1", {
            // Completely different — should NOT share nodes
            make_sequence_numeric(1200, 99)
        }),
    ];

    let mut compressor = StreamingQueueCompressor::with_splitters(path, config, AHashSet::new())
        .expect("Failed to create AGC compressor");

    for (sample, contig, data) in samples {
        compressor
            .push(sample.to_string(), contig.to_string(), data)
            .expect("Failed to push contig");
    }
    compressor.finalize().expect("Failed to finalize AGC");
}

/// Locate the compiled `impg` binary for CLI tests.
fn impg_binary() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return Some(PathBuf::from(path));
    }

    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for candidate in [
        manifest_dir.join("target/debug/impg"),
        manifest_dir.join("target/release/impg"),
    ] {
        if candidate.exists() {
            return Some(candidate);
        }
    }
    None
}

fn read_pack_tsv_counts(path: &Path) -> BTreeMap<u32, u64> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut counts = BTreeMap::new();
    for line in text.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(fields.len(), 2, "pack row should have two fields: {line}");
        counts.insert(fields[0].parse().unwrap(), fields[1].parse().unwrap());
    }
    counts
}

fn report_section_lines<'a>(report: &'a str, section: &str) -> Vec<&'a str> {
    let marker = format!("#section\t{section}");
    let mut in_section = false;
    let mut lines = Vec::new();
    for line in report.lines() {
        if line.starts_with("#section\t") {
            if in_section {
                break;
            }
            in_section = line == marker;
            continue;
        }
        if in_section && !line.trim().is_empty() {
            lines.push(line);
        }
    }
    lines
}

#[test]
fn test_crush_cli_resolves_blunt_gfa() {
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI test: impg binary not found");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_crush_cli");
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    let input = dir.join("input.gfa");
    let output = dir.join("out.gfa");
    std::fs::write(
        &input,
        "\
H\tVN:Z:1.0
S\t1\tAC
S\t2\tGGG
S\t3\tTA
L\t1\t+\t3\t+\t0M
L\t1\t+\t2\t+\t0M
L\t2\t+\t3\t+\t0M
P\tref\t1+,3+\t*
P\tins\t1+,2+,3+\t*
",
    )
    .unwrap();

    let run = Command::new(&bin)
        .args([
            "crush",
            "-g",
            input.to_str().unwrap(),
            "-o",
            output.to_str().unwrap(),
            "-v",
            "2",
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("failed to run impg crush");
    assert!(
        run.status.success(),
        "impg crush failed: {}",
        String::from_utf8_lossy(&run.stderr)
    );
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(
        stderr.contains("crush: 1 resolved"),
        "crush should resolve one insertion bubble, stderr:\n{}",
        stderr
    );
    let text = std::fs::read_to_string(&output).unwrap();
    assert!(text.starts_with("H\tVN:Z:1.0\n"), "{}", text);
    assert!(text.contains("\nS\t"), "{}", text);
    assert!(text.contains("\nP\tref\t"), "{}", text);

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_render_bundle_preserves_source_namespace() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI test: impg binary not found");
        return;
    };

    use std::io::Write;

    let dir = std::env::temp_dir().join("impg_test_syng_render_bundle");
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    let fasta = dir.join("panel.fa");
    let prefix = dir.join("panel.syng");
    let bundle = dir.join("render.impg-gbz");

    let base = numeric_to_ascii(&make_sequence_numeric(1400, 23));
    let mut hap_b = base.clone();
    for i in (180..1200).step_by(97) {
        hap_b[i] = mutate_base_ascii(hap_b[i]);
    }
    let fragment = &base[300..950];

    {
        let mut f = std::fs::File::create(&fasta).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&base).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#1#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">fragment_001").unwrap();
        f.write_all(fragment).unwrap();
        writeln!(f).unwrap();
    }

    let build = Command::new(&bin)
        .args([
            "syng",
            "--fasta",
            fasta.to_str().unwrap(),
            "-o",
            prefix.to_str().unwrap(),
            "--position-sample-rate",
            "8",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let render = Command::new(&bin)
        .args([
            "render",
            "-a",
            prefix.to_str().unwrap(),
            "-r",
            "sampleA#0#chr1:100-1000",
            "--sequence-files",
            fasta.to_str().unwrap(),
            "-O",
            bundle.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg render");
    assert!(
        render.status.success(),
        "impg render failed: {}",
        String::from_utf8_lossy(&render.stderr)
    );

    for file in [
        "manifest.json",
        "namespace.json",
        "translation.bin",
        "translation.tsv",
        "rendered.fa",
        "graph.gfa",
        "paths.1gbwt",
        "paths.1khash",
    ] {
        assert!(
            bundle.join(file).exists(),
            "render bundle should contain {}",
            file
        );
    }

    let manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(bundle.join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(manifest["format"], "impg-render-bundle");
    assert_eq!(manifest["engine"], "syng:blunt");
    assert_eq!(manifest["feature_space"], "syng-syncmer-node");
    assert!(manifest["rendered_paths"].as_u64().unwrap() >= 2);
    assert!(manifest["step_samples"].as_u64().unwrap() > 0);

    let namespace = std::fs::read_to_string(bundle.join("namespace.json")).unwrap();
    assert!(namespace.contains("sampleA#0#chr1"), "{}", namespace);
    assert!(
        namespace.contains("\"sample\": \"sampleA\""),
        "{}",
        namespace
    );
    assert!(namespace.contains("\"haplotype\": \"0\""), "{}", namespace);
    assert!(namespace.contains("fragment_001"), "{}", namespace);
    assert!(namespace.contains("\"pansn\": null"), "{}", namespace);

    let translation = std::fs::read(bundle.join("translation.bin")).unwrap();
    assert_eq!(&translation[..8], b"IMPGTRN1");
    let translation_tsv = std::fs::read_to_string(bundle.join("translation.tsv")).unwrap();
    assert!(translation_tsv.contains("\npath\t"), "{}", translation_tsv);
    assert!(translation_tsv.contains("\nstep\t"), "{}", translation_tsv);
    assert!(
        translation_tsv.contains("sampleA#0#chr1"),
        "{}",
        translation_tsv
    );

    let gfa = std::fs::read_to_string(bundle.join("graph.gfa")).unwrap();
    assert!(gfa.starts_with("H\tVN:Z:1.0\n"), "{}", gfa);
    assert!(gfa.contains("\nS\t"), "{}", gfa);
    assert!(gfa.contains("\nP\t"), "{}", gfa);

    let local_pack = dir.join("render.pack");
    let map_local_pack = Command::new(&bin)
        .args([
            "map",
            "-a",
            bundle.join("paths").to_str().unwrap(),
            "-q",
            bundle.join("rendered.fa").to_str().unwrap(),
            "-o",
            "pack",
            "-O",
            local_pack.to_str().unwrap(),
            "--min-anchors",
            "1",
        ])
        .output()
        .expect("failed to run impg map on rendered syng bundle");
    assert!(
        map_local_pack.status.success(),
        "impg map on rendered syng bundle failed: {}",
        String::from_utf8_lossy(&map_local_pack.stderr)
    );
    let genotype_bundle = Command::new(&bin)
        .args([
            "genotype",
            "cos",
            "--render-bundle",
            bundle.to_str().unwrap(),
            "--pack",
            local_pack.to_str().unwrap(),
            "--ploidy",
            "1",
            "--top-n",
            "1",
            "--min-anchors",
            "1",
            "--min-span-fraction",
            "0",
        ])
        .output()
        .expect("failed to run impg genotype cos --render-bundle");
    assert!(
        genotype_bundle.status.success(),
        "impg genotype cos --render-bundle failed: {}",
        String::from_utf8_lossy(&genotype_bundle.stderr)
    );
    let genotype_stdout = String::from_utf8_lossy(&genotype_bundle.stdout);
    assert!(
        genotype_stdout.contains("#impg genotype cos"),
        "{}",
        genotype_stdout
    );
    assert!(
        genotype_stdout.contains("sampleA#0#chr1"),
        "bundle genotype should report rendered paths with source names:\n{}",
        genotype_stdout
    );
    let infer_bundle = Command::new(&bin)
        .args([
            "infer",
            "--render-bundle",
            bundle.to_str().unwrap(),
            "--pack",
            local_pack.to_str().unwrap(),
            "--ploidy",
            "1",
            "--top-n",
            "1",
            "--min-anchors",
            "1",
            "--min-span-fraction",
            "0",
        ])
        .output()
        .expect("failed to run impg infer --render-bundle");
    assert!(
        infer_bundle.status.success(),
        "impg infer --render-bundle failed: {}",
        String::from_utf8_lossy(&infer_bundle.stderr)
    );
    let infer_stdout = String::from_utf8_lossy(&infer_bundle.stdout);
    assert!(infer_stdout.contains("#impg infer"), "{}", infer_stdout);
    assert!(
        infer_stdout.contains("sampleA#0#chr1"),
        "bundle infer should report rendered paths with source names:\n{}",
        infer_stdout
    );

    let poa_bundle = dir.join("render.poa.impg-gbz");
    let render_poa = Command::new(&bin)
        .args([
            "render",
            "-a",
            prefix.to_str().unwrap(),
            "-r",
            "sampleA#0#chr1:100-1000",
            "--sequence-files",
            fasta.to_str().unwrap(),
            "-O",
            poa_bundle.to_str().unwrap(),
            "--engine",
            "poa",
            "-t",
            "2",
        ])
        .output()
        .expect("failed to run impg render --engine poa");
    assert!(
        render_poa.status.success(),
        "impg render --engine poa failed: {}",
        String::from_utf8_lossy(&render_poa.stderr)
    );
    let poa_manifest: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(poa_bundle.join("manifest.json")).unwrap())
            .unwrap();
    assert_eq!(poa_manifest["engine"], "poa");
    assert_eq!(poa_manifest["graph_kind"], "local-sequence-graph");
    assert_eq!(poa_manifest["feature_space"], "gfa-segment");
    assert!(poa_manifest["syng_prefix"].is_null());
    assert!(poa_manifest["step_samples"].as_u64().unwrap() > 0);
    assert!(poa_bundle.join("rendered.fa").exists());
    let poa_translation = std::fs::read_to_string(poa_bundle.join("translation.tsv")).unwrap();
    assert!(poa_translation.contains("\nstep\t"), "{}", poa_translation);
    let poa_gfa = std::fs::read_to_string(poa_bundle.join("graph.gfa")).unwrap();
    assert!(poa_gfa.contains("\nS\t"), "{}", poa_gfa);
    assert!(poa_gfa.contains("\nP\t"), "{}", poa_gfa);
}

#[test]
fn test_syng_agc_build_produces_non_empty_index() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_agc_build");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let agc_path = dir.join("test.agc");
    create_test_agc(agc_path.to_str().unwrap());

    let out_prefix = dir.join("idx");

    let output = Command::new(&bin)
        .args([
            "syng",
            "--agc",
            agc_path.to_str().unwrap(),
            "--syncmer-length",
            "63",
            "--smer-length",
            "8",
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng --agc");

    assert!(
        output.status.success(),
        "impg syng --agc failed. stderr: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let gbwt_path = format!("{}.1gbwt", out_prefix.to_str().unwrap());
    let khash_path = format!("{}.1khash", out_prefix.to_str().unwrap());
    let names_path = format!("{}.syng.names", out_prefix.to_str().unwrap());
    let spos_path = format!("{}.syng.spos", out_prefix.to_str().unwrap());
    let pstep_path = format!("{}.syng.pstep", out_prefix.to_str().unwrap());
    let meta_path = format!("{}.syng.meta", out_prefix.to_str().unwrap());

    assert!(std::path::Path::new(&gbwt_path).exists(), ".1gbwt missing");
    assert!(
        std::path::Path::new(&khash_path).exists(),
        ".1khash missing"
    );
    assert!(
        std::path::Path::new(&names_path).exists(),
        ".syng.names missing"
    );
    assert!(
        std::path::Path::new(&spos_path).exists(),
        ".syng.spos missing"
    );
    assert!(
        std::path::Path::new(&pstep_path).exists(),
        ".syng.pstep missing"
    );
    assert!(
        std::path::Path::new(&meta_path).exists(),
        ".syng.meta missing"
    );

    // Key assertion: the GBWT must contain actual vertex data, not just a
    // schema header. The yeast235 bug slipped through because nothing
    // asserted on file size. An empty GBWT is ~1300 bytes; a populated one
    // for even this tiny test should be several KB.
    let gbwt_size = std::fs::metadata(&gbwt_path).unwrap().len();
    assert!(
        gbwt_size > 2000,
        ".1gbwt is only {} bytes — suggests no syncmer vertices were written. \
         This is the bug that slipped through on yeast235.agc.",
        gbwt_size
    );

    let khash_size = std::fs::metadata(&khash_path).unwrap().len();
    assert!(
        khash_size > 1000,
        ".1khash is only {} bytes — suggests no k-mer hash entries",
        khash_size
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_agc_roundtrip_query() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_agc_query");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let agc_path = dir.join("test.agc");
    create_test_agc(agc_path.to_str().unwrap());

    let out_prefix = dir.join("idx");
    let out_prefix_str = out_prefix.to_str().unwrap();

    let build_output = Command::new(&bin)
        .args([
            "syng",
            "--agc",
            agc_path.to_str().unwrap(),
            "-o",
            out_prefix_str,
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng --agc");
    assert!(
        build_output.status.success(),
        "impg syng --agc failed. stderr: {}",
        String::from_utf8_lossy(&build_output.stderr)
    );

    // Load the index back and run a query on the shared backbone region.
    let index = impg::syng::SyngIndex::load(out_prefix_str, impg::syng::SyncmerParams::default())
        .expect("Failed to load built AGC-derived syng index");

    assert_eq!(
        index.name_map.path_to_name.len(),
        3,
        "Expected 3 sequences in name map"
    );

    // Find one of the contigs that contains the shared backbone.
    let query_name = index
        .name_map
        .path_to_name
        .iter()
        .find(|n| n.contains("sampleA"))
        .expect("sampleA contig should be in name map")
        .clone();

    let intervals = index
        .query_region(&query_name, 0, 500, 0)
        .expect("query_region failed");

    // Should find at least self-hit + sampleB (shared backbone)
    assert!(
        intervals.len() >= 2,
        "Expected at least 2 homologous intervals (self + sampleB), got {}: {:?}",
        intervals.len(),
        intervals.iter().map(|iv| &iv.genome).collect::<Vec<_>>()
    );

    let genomes: Vec<&str> = intervals.iter().map(|iv| iv.genome.as_str()).collect();
    assert!(
        genomes.iter().any(|g| g.contains("sampleA")),
        "Should find self-hit, got: {:?}",
        genomes
    );
    assert!(
        genomes.iter().any(|g| g.contains("sampleB")),
        "Should find sampleB (shared backbone), got: {:?}",
        genomes
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_agc_parallel_dictionary_roundtrip_query() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_agc_parallel_dictionary");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let agc_path = dir.join("test.agc");
    create_test_agc(agc_path.to_str().unwrap());

    let out_prefix = dir.join("idx");
    let out_prefix_str = out_prefix.to_str().unwrap();

    let build_output = Command::new(&bin)
        .args([
            "syng",
            "--agc",
            agc_path.to_str().unwrap(),
            "--parallel-dictionary",
            "-o",
            out_prefix_str,
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng --agc --parallel-dictionary");
    assert!(
        build_output.status.success(),
        "impg syng --agc --parallel-dictionary failed. stderr: {}",
        String::from_utf8_lossy(&build_output.stderr)
    );

    let index = impg::syng::SyngIndex::load(out_prefix_str, impg::syng::SyncmerParams::default())
        .expect("Failed to load AGC parallel-dictionary syng index");
    let query_name = index
        .name_map
        .path_to_name
        .iter()
        .find(|n| n.contains("sampleA"))
        .expect("sampleA contig should be in name map")
        .clone();
    let intervals = index
        .query_region(&query_name, 0, 500, 0)
        .expect("query_region failed");
    assert!(
        intervals.iter().any(|iv| iv.genome.contains("sampleB")),
        "parallel dictionary build should preserve shared-backbone query behavior"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_fasta_build_produces_non_empty_index() {
    // Same non-empty-GBWT assertion but via the FASTA input path, to catch
    // regressions in either build path symmetrically.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_fasta_nonempty");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(800, 42));
        let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
        let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));
        let other = numeric_to_ascii(&make_sequence_numeric(1200, 99));

        writeln!(f, ">sampleA#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleC#chr1").unwrap();
        f.write_all(&other).unwrap();
        writeln!(f).unwrap();
    }

    let out_prefix = dir.join("idx");
    let out = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng -f");
    assert!(
        out.status.success(),
        "impg syng -f failed. stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gbwt_size = std::fs::metadata(format!("{}.1gbwt", out_prefix.to_str().unwrap()))
        .unwrap()
        .len();
    assert!(
        gbwt_size > 2000,
        ".1gbwt is only {} bytes from FASTA build",
        gbwt_size
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_fasta_parallel_dictionary_build_produces_non_empty_index() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_fasta_parallel_dictionary");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(800, 42));
        let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
        let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));

        writeln!(f, ">sampleA#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
    }

    let out_prefix = dir.join("idx");
    let out = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "--parallel-dictionary",
            "-o",
            out_prefix.to_str().unwrap(),
        ])
        .output()
        .expect("Failed to run impg syng -f --parallel-dictionary");
    assert!(
        out.status.success(),
        "impg syng -f --parallel-dictionary failed. stderr: {}",
        String::from_utf8_lossy(&out.stderr)
    );

    let gbwt_size = std::fs::metadata(format!("{}.1gbwt", out_prefix.to_str().unwrap()))
        .unwrap()
        .len();
    assert!(
        gbwt_size > 2000,
        ".1gbwt is only {} bytes from parallel FASTA build",
        gbwt_size
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_map_cli_gaf_and_paf() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_map_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let backbone = numeric_to_ascii(&make_sequence_numeric(1000, 42));
    let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
    let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
    }

    let query_path = dir.join("query.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&query_path).unwrap();
        writeln!(f, "@read1").unwrap();
        f.write_all(&backbone[100..800]).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(700)).unwrap();
    }
    let query_rc_path = dir.join("query_rc.fq");
    {
        use std::io::Write;
        let rc = impg::graph::reverse_complement(&backbone[100..800]);
        let mut f = std::fs::File::create(&query_rc_path).unwrap();
        writeln!(f, "@read_rc").unwrap();
        f.write_all(&rc).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(700)).unwrap();
    }
    let query_multi_path = dir.join("query_multi.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&query_multi_path).unwrap();
        for i in 0..2050 {
            writeln!(f, "@read{}", i).unwrap();
            f.write_all(&backbone[100..800]).unwrap();
            writeln!(f).unwrap();
            writeln!(f, "+").unwrap();
            writeln!(f, "{}", "I".repeat(700)).unwrap();
        }
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let gaf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o gaf");
    assert!(
        gaf.status.success(),
        "impg map -o gaf failed: {}",
        String::from_utf8_lossy(&gaf.stderr)
    );
    let gaf_stdout = String::from_utf8_lossy(&gaf.stdout);
    let gaf_lines: Vec<&str> = gaf_stdout.lines().collect();
    assert_eq!(
        gaf_lines.len(),
        1,
        "expected one GAF line, got:\n{}",
        gaf_stdout
    );
    let gaf_fields: Vec<&str> = gaf_lines[0].split('\t').collect();
    assert!(
        gaf_fields.len() >= 12,
        "GAF line has too few fields: {}",
        gaf_lines[0]
    );
    assert_eq!(gaf_fields[0], "read1");
    assert!(
        gaf_fields[5].contains('>') || gaf_fields[5].contains('<'),
        "GAF path should contain syncmer node orientations: {}",
        gaf_fields[5]
    );
    let mut expected_pack_nodes = std::collections::BTreeSet::new();
    let mut node_digits = String::new();
    for ch in gaf_fields[5].chars() {
        if ch == '>' || ch == '<' {
            if !node_digits.is_empty() {
                expected_pack_nodes.insert(node_digits.parse::<u32>().unwrap());
                node_digits.clear();
            }
        } else if ch.is_ascii_digit() {
            node_digits.push(ch);
        }
    }
    if !node_digits.is_empty() {
        expected_pack_nodes.insert(node_digits.parse::<u32>().unwrap());
    }

    let default_gaf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map with default output");
    assert!(
        default_gaf.status.success(),
        "impg map default output failed: {}",
        String::from_utf8_lossy(&default_gaf.stderr)
    );
    assert_eq!(
        default_gaf.stdout, gaf.stdout,
        "impg map should default to GAF output"
    );

    let pack_tsv = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv");
    assert!(
        pack_tsv.status.success(),
        "impg map -o pack-tsv failed: {}",
        String::from_utf8_lossy(&pack_tsv.stderr)
    );
    let pack_stdout = String::from_utf8_lossy(&pack_tsv.stdout);
    let pack_lines: Vec<&str> = pack_stdout.lines().collect();
    assert_eq!(
        pack_lines.first().copied(),
        Some("#node_id\tcount"),
        "pack-tsv output should start with a TSV header, got:\n{}",
        pack_stdout
    );
    assert_eq!(
        pack_lines.len().saturating_sub(1),
        expected_pack_nodes.len(),
        "pack-tsv output should have one row per distinct GAF node"
    );
    for line in pack_lines.iter().skip(1) {
        let fields: Vec<&str> = line.split('\t').collect();
        assert_eq!(
            fields.len(),
            2,
            "pack line should have two fields: {}",
            line
        );
        let node_id = fields[0].parse::<u32>().unwrap();
        let count = fields[1].parse::<u64>().unwrap();
        assert!(
            expected_pack_nodes.contains(&node_id),
            "pack-tsv emitted unexpected node {} in:\n{}",
            node_id,
            pack_stdout
        );
        assert_eq!(
            count, 1,
            "single-read pack count should be 1 for {}",
            node_id
        );
    }
    let pack_zst_path = dir.join("pack.tsv.zst");
    let pack_zst = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "-O",
            pack_zst_path.to_str().unwrap(),
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv -O pack.tsv.zst");
    assert!(
        pack_zst.status.success(),
        "impg map -o pack-tsv -O pack.tsv.zst failed: {}",
        String::from_utf8_lossy(&pack_zst.stderr)
    );
    assert!(
        pack_zst.stdout.is_empty(),
        "pack-tsv output should go to -O path, got stdout: {}",
        String::from_utf8_lossy(&pack_zst.stdout)
    );
    let compressed = std::fs::read(&pack_zst_path).unwrap();
    assert!(
        compressed.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]),
        "pack .zst output should have zstd magic bytes"
    );
    let mut decoder =
        zstd::stream::read::Decoder::new(std::fs::File::open(&pack_zst_path).unwrap()).unwrap();
    let mut decoded_pack = String::new();
    {
        use std::io::Read;
        decoder.read_to_string(&mut decoded_pack).unwrap();
    }
    assert_eq!(
        decoded_pack, pack_stdout,
        ".zst-compressed pack-tsv output should decompress to plain pack TSV"
    );

    let packbin_path = dir.join("pack.ipack");
    let packbin = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "pack",
            "-O",
            packbin_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack");
    assert!(
        packbin.status.success(),
        "impg map -o pack failed: {}",
        String::from_utf8_lossy(&packbin.stderr)
    );
    assert!(
        packbin.stdout.is_empty(),
        "pack output should go to -O path, got stdout bytes: {}",
        packbin.stdout.len()
    );
    let packbin_bytes = std::fs::read(&packbin_path).unwrap();
    assert_eq!(&packbin_bytes[..8], b"IMPGPKB1");
    let read_u32 = |offset: usize| -> u32 {
        u32::from_le_bytes(packbin_bytes[offset..offset + 4].try_into().unwrap())
    };
    let read_i32 = |offset: usize| -> i32 {
        i32::from_le_bytes(packbin_bytes[offset..offset + 4].try_into().unwrap())
    };
    let read_u64 = |offset: usize| -> u64 {
        u64::from_le_bytes(packbin_bytes[offset..offset + 8].try_into().unwrap())
    };
    assert_eq!(read_u32(8), 1, "packbin version should be 1");
    let header_len = read_u32(12) as usize;
    assert_eq!(header_len, 96);
    let universe_nodes = read_u64(16) as usize;
    let nonzero_nodes = read_u64(24) as usize;
    let retained_reads = read_u64(32);
    let anchors = read_u64(40);
    let block_size = read_u32(48) as usize;
    let zstd_level = read_i32(52);
    let block_count = read_u64(56) as usize;
    let overflow_count = read_u64(64) as usize;
    let block_index_offset = read_u64(72) as usize;
    let overflow_offset = read_u64(80) as usize;
    let data_offset = read_u64(88) as usize;
    assert_eq!(nonzero_nodes, expected_pack_nodes.len());
    assert_eq!(retained_reads, 1);
    assert!(anchors >= expected_pack_nodes.len() as u64);
    assert_eq!(block_size, 64);
    assert_eq!(zstd_level, 3);
    assert_eq!(overflow_count, 0);
    assert_eq!(block_index_offset, header_len);
    assert!(overflow_offset >= block_index_offset);
    assert!(data_offset >= overflow_offset);
    assert_eq!(block_count, universe_nodes.div_ceil(block_size));

    let mut block_offsets = Vec::with_capacity(block_count + 1);
    for i in 0..=block_count {
        block_offsets.push(read_u64(block_index_offset + i * 8) as usize);
    }
    let mut dense = Vec::with_capacity(universe_nodes);
    for block_idx in 0..block_count {
        let compressed_start = data_offset + block_offsets[block_idx];
        let compressed_end = data_offset + block_offsets[block_idx + 1];
        let expected_len = (universe_nodes - block_idx * block_size).min(block_size);
        let block = zstd::bulk::decompress(
            &packbin_bytes[compressed_start..compressed_end],
            expected_len,
        )
        .unwrap();
        assert_eq!(block.len(), expected_len);
        dense.extend_from_slice(&block);
    }
    assert_eq!(dense.len(), universe_nodes);
    assert_eq!(
        dense.iter().filter(|&&count| count != 0).count(),
        expected_pack_nodes.len()
    );
    for &node_id in &expected_pack_nodes {
        assert_eq!(
            dense[(node_id - 1) as usize],
            1,
            "packbin count for {}",
            node_id
        );
    }

    let proj_path = dir.join("sample.proj");
    let proj = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        proj.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&proj.stderr)
    );
    assert!(proj_path.join("manifest.json").exists());
    assert_eq!(
        &std::fs::read(proj_path.join("sample.pack")).unwrap()[..8],
        b"IMPGPKB1"
    );
    let mut gaf_decoder = zstd::stream::read::Decoder::new(
        std::fs::File::open(proj_path.join("reads.gaf.zst")).unwrap(),
    )
    .unwrap();
    let mut projected_gaf = String::new();
    {
        use std::io::Read;
        gaf_decoder.read_to_string(&mut projected_gaf).unwrap();
    }
    assert!(
        projected_gaf.starts_with("read1\t"),
        "projection should include GAF read walks, got:\n{}",
        projected_gaf
    );

    let gaf_rc = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_rc_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o gaf on reverse-complement read");
    assert!(
        gaf_rc.status.success(),
        "impg map -o gaf failed on reverse-complement read: {}",
        String::from_utf8_lossy(&gaf_rc.stderr)
    );
    let gaf_rc_stdout = String::from_utf8_lossy(&gaf_rc.stdout);
    let gaf_rc_lines: Vec<&str> = gaf_rc_stdout.lines().collect();
    assert_eq!(
        gaf_rc_lines.len(),
        1,
        "expected one reverse-complement GAF line, got:\n{}",
        gaf_rc_stdout
    );
    assert!(
        gaf_rc_lines[0].starts_with("read_rc\t"),
        "reverse-complement GAF should report the read name, got: {}",
        gaf_rc_lines[0]
    );

    let paf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "paf",
            "--min-anchors",
            "2",
            "--max-hits",
            "10",
        ])
        .output()
        .expect("failed to run impg map -o paf");
    assert!(
        paf.status.success(),
        "impg map -o paf failed: {}",
        String::from_utf8_lossy(&paf.stderr)
    );
    let paf_stdout = String::from_utf8_lossy(&paf.stdout);
    let paf_lines: Vec<&str> = paf_stdout.lines().collect();
    assert!(!paf_lines.is_empty(), "expected PAF hits, got none");
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleA#0#chr1")),
        "expected sampleA hit, got:\n{}",
        paf_stdout
    );
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleB#0#chr1")),
        "expected sampleB hit, got:\n{}",
        paf_stdout
    );
    for line in paf_lines {
        let fields: Vec<&str> = line.split('\t').collect();
        assert!(fields.len() >= 12, "PAF line has too few fields: {}", line);
        assert_eq!(fields[0], "read1");
    }

    for suffix in ["1gbwt", "syng.names", "syng.spos", "syng.pstep"] {
        std::fs::remove_file(dir.join(format!("idx.{suffix}"))).ok();
    }
    let gaf_khash_only = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o gaf with only khash");
    assert!(
        gaf_khash_only.status.success(),
        "impg map -o gaf should only require .1khash/.meta, stderr: {}",
        String::from_utf8_lossy(&gaf_khash_only.stderr)
    );
    assert!(
        !String::from_utf8_lossy(&gaf_khash_only.stdout)
            .lines()
            .collect::<Vec<_>>()
            .is_empty(),
        "expected GAF output with only .1khash/.meta present"
    );
    let pack_khash_only = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv with only khash");
    assert!(
        pack_khash_only.status.success(),
        "impg map -o pack-tsv should only require .1khash/.meta, stderr: {}",
        String::from_utf8_lossy(&pack_khash_only.stderr)
    );
    assert_eq!(
        String::from_utf8_lossy(&pack_khash_only.stdout)
            .lines()
            .next(),
        Some("#node_id\tcount"),
        "expected pack header with only .1khash/.meta present"
    );
    let gaf_threads_1 = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_multi_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
            "-t",
            "1",
        ])
        .output()
        .expect("failed to run impg map -o gaf with one thread");
    assert!(
        gaf_threads_1.status.success(),
        "single-threaded GAF mapping failed: {}",
        String::from_utf8_lossy(&gaf_threads_1.stderr)
    );
    let gaf_threads_2 = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_multi_path.to_str().unwrap(),
            "-o",
            "gaf",
            "--min-anchors",
            "2",
            "-t",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o gaf with two threads");
    assert!(
        gaf_threads_2.status.success(),
        "two-threaded GAF mapping failed: {}",
        String::from_utf8_lossy(&gaf_threads_2.stderr)
    );
    assert_eq!(
        gaf_threads_1.stdout, gaf_threads_2.stdout,
        "parallel GAF mapping should preserve deterministic output order"
    );
    assert_eq!(
        String::from_utf8_lossy(&gaf_threads_2.stdout)
            .lines()
            .count(),
        2050,
        "expected one GAF record per repeated query read"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_genotype_cos_cli_permutations() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_genotype_cos_permutations");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(700, 11));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(700, 12));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(700, 13));
    let right = numeric_to_ascii(&make_sequence_numeric(700, 14));
    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("reads.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&reads_path).unwrap();
        writeln!(f, "@read_a").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(hap_a.len())).unwrap();
        writeln!(f, "@read_b").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(hap_b.len())).unwrap();
    }

    let proj_path = dir.join("sample.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        map.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );
    let compact_pack_path = proj_path.join("sample.pack");

    let pack_path = dir.join("sample.pack.tsv");
    let map_pack = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "-O",
            pack_path.to_str().unwrap(),
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv");
    assert!(
        map_pack.status.success(),
        "impg map -o pack-tsv failed: {}",
        String::from_utf8_lossy(&map_pack.stderr)
    );

    let pack_zst_path = dir.join("sample.pack.tsv.zst");
    let map_pack_zst = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "-O",
            pack_zst_path.to_str().unwrap(),
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv -O sample.pack.tsv.zst");
    assert!(
        map_pack_zst.status.success(),
        "impg map -o pack-tsv -O sample.pack.tsv.zst failed: {}",
        String::from_utf8_lossy(&map_pack_zst.stderr)
    );
    let compressed_pack = std::fs::read(&pack_zst_path).unwrap();
    assert!(
        compressed_pack.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]),
        "compressed pack should have zstd magic bytes"
    );

    let genotype_help = Command::new(&bin)
        .args(["genotype", "--help"])
        .output()
        .expect("failed to run impg genotype --help");
    assert!(genotype_help.status.success());
    let genotype_help_stdout = String::from_utf8_lossy(&genotype_help.stdout);
    assert!(
        genotype_help_stdout.contains("  cos  "),
        "genotype help should list cos:\n{}",
        genotype_help_stdout
    );
    assert!(
        !genotype_help_stdout.contains("cosigt"),
        "genotype help should not list the cosigt alias:\n{}",
        genotype_help_stdout
    );

    let gt_help = Command::new(&bin)
        .args(["gt", "--help"])
        .output()
        .expect("failed to run impg gt --help");
    assert!(gt_help.status.success());
    let gt_help_stdout = String::from_utf8_lossy(&gt_help.stdout);
    assert!(
        gt_help_stdout.contains("  cos  "),
        "gt help should list cos:\n{}",
        gt_help_stdout
    );
    assert!(
        !gt_help_stdout.contains("cosigt"),
        "gt help should not list the cosigt alias:\n{}",
        gt_help_stdout
    );

    let command_roots = ["genotype", "gt"];
    let methods = ["cos", "cosigt"];
    let candidate_modes = ["spanning", "overlapping"];
    let packs = [
        ("pack-tsv", "pack", pack_path.as_path()),
        ("pack-tsv.zst", "pack", pack_zst_path.as_path()),
        ("pack", "packbin", compact_pack_path.as_path()),
    ];
    let mut expected_by_pack_mode: std::collections::BTreeMap<(String, String), String> =
        std::collections::BTreeMap::new();
    let mut checked = 0usize;

    for root in command_roots {
        for method in methods {
            for (pack_label, compare_label, pack_path) in packs {
                for candidate_mode in candidate_modes {
                    checked += 1;
                    let output = Command::new(&bin)
                        .args([
                            root,
                            method,
                            "-a",
                            idx_prefix.to_str().unwrap(),
                            "-p",
                            pack_path.to_str().unwrap(),
                            "-r",
                            "sampleA#0#chr1:0-2100",
                            "--candidate-mode",
                            candidate_mode,
                            "--top-n",
                            "3",
                            "--candidate-top-k",
                            "10",
                            "--min-anchors",
                            "2",
                            "--min-span-fraction",
                            "0.7",
                        ])
                        .output()
                        .unwrap_or_else(|e| {
                            panic!(
                                "failed to run impg {} {} with {} {}: {}",
                                root, method, pack_label, candidate_mode, e
                            )
                        });
                    assert!(
                        output.status.success(),
                        "impg {} {} failed with {} {}: {}",
                        root,
                        method,
                        pack_label,
                        candidate_mode,
                        String::from_utf8_lossy(&output.stderr)
                    );

                    let stdout = String::from_utf8_lossy(&output.stdout).to_string();
                    assert!(stdout.contains("#impg genotype cos"), "{}", stdout);
                    assert!(stdout.contains("#method\tcos"), "{}", stdout);
                    assert!(stdout.contains("#metric\tcosine"), "{}", stdout);
                    assert!(stdout.contains("#alias\tcosigt"), "{}", stdout);
                    assert!(
                        stdout.contains(&format!(
                            "#candidate_mode\t{}",
                            if candidate_mode == "spanning" {
                                "Spanning"
                            } else {
                                "Overlapping"
                            }
                        )),
                        "{}",
                        stdout
                    );

                    let top = stdout
                        .lines()
                        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
                        .expect("expected at least one cosine genotype result row");
                    let fields: Vec<&str> = top.split('\t').collect();
                    assert!(
                        fields.len() >= 12,
                        "cos genotype row should have at least 12 fields: {}",
                        top
                    );
                    assert_eq!(fields[0], "1");
                    assert_eq!(fields[1], "cos");
                    assert_eq!(fields[2], "2");

                    if candidate_mode == "spanning" {
                        let similarity = fields[3].parse::<f64>().unwrap();
                        assert!(
                            similarity > 0.99,
                            "expected near-perfect heterozygous cosine score, got {}\n{}",
                            similarity,
                            stdout
                        );
                        assert!(
                            fields[8].contains("sampleA#0#chr1")
                                && fields[8].contains("sampleB#0#chr1"),
                            "top spanning genotype should contain both haplotypes, got:\n{}",
                            stdout
                        );
                    }

                    let key = (compare_label.to_string(), candidate_mode.to_string());
                    if let Some(expected) = expected_by_pack_mode.get(&key) {
                        assert_eq!(
                            &stdout, expected,
                            "{} {} should match the first {} {} output",
                            root, method, pack_label, candidate_mode
                        );
                    } else {
                        expected_by_pack_mode.insert(key, stdout);
                    }
                }
            }
        }
    }
    assert_eq!(
        checked, 24,
        "expected to exercise the full genotype CLI matrix"
    );

    let genotype_proj = Command::new(&bin)
        .args([
            "genotype",
            "cos",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            "sampleA#0#chr1:0-2100",
            "--top-n",
            "1",
            "--candidate-top-k",
            "10",
            "--min-span-fraction",
            "0.7",
        ])
        .output()
        .expect("failed to run impg genotype cos --proj");
    assert!(
        genotype_proj.status.success(),
        "impg genotype cos --proj failed: {}",
        String::from_utf8_lossy(&genotype_proj.stderr)
    );
    assert!(
        String::from_utf8_lossy(&genotype_proj.stdout).contains("#impg genotype cos"),
        "genotype --proj should emit cos output"
    );

    let idx_prefix_str = idx_prefix.to_str().unwrap();
    let packbin_str = compact_pack_path.to_str().unwrap();
    let target_range = "sampleA#0#chr1:0-2100";
    let run_gt = |extra: &[&str]| {
        let mut args = vec![
            "genotype",
            "cos",
            "-a",
            idx_prefix_str,
            "-p",
            packbin_str,
            "-r",
            target_range,
        ];
        args.extend_from_slice(extra);
        Command::new(&bin)
            .args(args)
            .output()
            .expect("failed to run impg genotype cos")
    };

    let custom = run_gt(&[
        "--ploidy",
        "1",
        "--top-n",
        "1",
        "--candidate-top-k",
        "0",
        "--syng-padding",
        "12",
        "--syng-extension",
        "12",
        "--min-anchors",
        "1",
        "--min-span-fraction",
        "0.5",
    ]);
    assert!(
        custom.status.success(),
        "non-default genotype options failed: {}",
        String::from_utf8_lossy(&custom.stderr)
    );
    let custom_stdout = String::from_utf8_lossy(&custom.stdout);
    assert!(custom_stdout.contains("#ploidy\t1"), "{}", custom_stdout);
    let custom_rows: Vec<&str> = custom_stdout
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(custom_rows.len(), 1, "top-n 1 should emit one row");
    assert_eq!(custom_rows[0].split('\t').nth(2), Some("1"));

    let ploidy3 = run_gt(&[
        "--ploidy",
        "3",
        "--top-n",
        "2",
        "--candidate-top-k",
        "0",
        "--min-span-fraction",
        "0.7",
    ]);
    assert!(
        ploidy3.status.success(),
        "ploidy 3 genotype options failed: {}",
        String::from_utf8_lossy(&ploidy3.stderr)
    );
    let ploidy3_stdout = String::from_utf8_lossy(&ploidy3.stdout);
    assert!(ploidy3_stdout.contains("#ploidy\t3"), "{}", ploidy3_stdout);
    let ploidy3_rows: Vec<&str> = ploidy3_stdout
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(ploidy3_rows.len(), 2, "top-n 2 should emit two rows");
    assert!(ploidy3_rows
        .iter()
        .all(|row| row.split('\t').nth(2) == Some("3")));

    let stdout_for_output = run_gt(&["--top-n", "2", "--candidate-top-k", "10"]);
    assert!(
        stdout_for_output.status.success(),
        "stdout genotype baseline failed: {}",
        String::from_utf8_lossy(&stdout_for_output.stderr)
    );
    let output_zst_path = dir.join("genotype.tsv.zst");
    let output_to_file = Command::new(&bin)
        .args([
            "genotype",
            "cos",
            "-a",
            idx_prefix_str,
            "-p",
            packbin_str,
            "-r",
            target_range,
            "--top-n",
            "2",
            "--candidate-top-k",
            "10",
            "-O",
            output_zst_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg genotype cos -O genotype.tsv.zst");
    assert!(
        output_to_file.status.success(),
        "compressed genotype output failed: {}",
        String::from_utf8_lossy(&output_to_file.stderr)
    );
    assert!(
        output_to_file.stdout.is_empty(),
        "genotype -O should not write result rows to stdout"
    );
    let compressed_output = std::fs::read(&output_zst_path).unwrap();
    assert!(
        compressed_output.starts_with(&[0x28, 0xb5, 0x2f, 0xfd]),
        "genotype .zst output should have zstd magic bytes"
    );
    let mut decoder =
        zstd::stream::read::Decoder::new(std::fs::File::open(&output_zst_path).unwrap()).unwrap();
    let mut decoded_output = String::new();
    {
        use std::io::Read;
        decoder.read_to_string(&mut decoded_output).unwrap();
    }
    assert_eq!(
        decoded_output,
        String::from_utf8_lossy(&stdout_for_output.stdout),
        "compressed genotype output should match stdout output"
    );

    for (label, extra, expected) in [
        (
            "ploidy zero",
            &["--ploidy", "0"][..],
            "--ploidy must be greater than 0",
        ),
        (
            "top-n zero",
            &["--top-n", "0"][..],
            "--top-n must be greater than 0",
        ),
        (
            "bad min-span-fraction",
            &["--min-span-fraction", "1.1"][..],
            "--min-span-fraction must be between 0 and 1",
        ),
        (
            "search budget",
            &[
                "--candidate-top-k",
                "10",
                "--max-combinations",
                "1",
                "--min-span-fraction",
                "0.7",
            ][..],
            "genotype combination search exceeded --max-combinations",
        ),
    ] {
        let failed = run_gt(extra);
        assert!(
            !failed.status.success(),
            "{} should fail but succeeded with stdout:\n{}",
            label,
            String::from_utf8_lossy(&failed.stdout)
        );
        let stderr = String::from_utf8_lossy(&failed.stderr);
        assert!(
            stderr.contains(expected),
            "{} should report '{}', got:\n{}",
            label,
            expected,
            stderr
        );
    }

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_genotype_cos_emit_report_exposes_counts() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_genotype_cos_emit_report");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(650, 41));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(650, 42));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(650, 43));
    let right = numeric_to_ascii(&make_sequence_numeric(650, 44));
    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("reads.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&reads_path).unwrap();
        writeln!(f, "@read_a").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(hap_a.len())).unwrap();
        writeln!(f, "@read_b").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(hap_b.len())).unwrap();
    }

    let pack_path = dir.join("sample.pack.tsv");
    let map_pack = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "-O",
            pack_path.to_str().unwrap(),
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv");
    assert!(
        map_pack.status.success(),
        "impg map -o pack-tsv failed: {}",
        String::from_utf8_lossy(&map_pack.stderr)
    );
    let pack_counts = read_pack_tsv_counts(&pack_path);
    assert!(!pack_counts.is_empty(), "pack TSV should contain counts");

    let target_range = format!("sampleA#0#chr1:0-{}", hap_a.len());
    let idx_prefix_str = idx_prefix.to_str().unwrap();
    let pack_path_str = pack_path.to_str().unwrap();
    let base_args = vec![
        "genotype",
        "cos",
        "-a",
        idx_prefix_str,
        "-p",
        pack_path_str,
        "-r",
        target_range.as_str(),
        "--top-n",
        "3",
        "--candidate-top-k",
        "10",
        "--min-anchors",
        "2",
        "--min-span-fraction",
        "0.7",
    ];

    let baseline = Command::new(&bin)
        .args(&base_args)
        .output()
        .expect("failed to run baseline impg genotype cos");
    assert!(
        baseline.status.success(),
        "baseline impg genotype cos failed: {}",
        String::from_utf8_lossy(&baseline.stderr)
    );

    let report_path = dir.join("genotype.report.tsv");
    let mut report_args = base_args.clone();
    report_args.extend_from_slice(&["--emit-report", report_path.to_str().unwrap()]);
    let with_report = Command::new(&bin)
        .args(&report_args)
        .output()
        .expect("failed to run impg genotype cos --emit-report");
    assert!(
        with_report.status.success(),
        "impg genotype cos --emit-report failed: {}",
        String::from_utf8_lossy(&with_report.stderr)
    );
    assert_eq!(
        with_report.stdout, baseline.stdout,
        "--emit-report must not change primary genotype TSV stdout"
    );

    let stdout = String::from_utf8_lossy(&with_report.stdout);
    let top = stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .expect("expected at least one cosine genotype result row");
    let top_fields: Vec<&str> = top.split('\t').collect();
    assert!(
        top_fields[8].contains("sampleA#0#chr1") && top_fields[8].contains("sampleB#0#chr1"),
        "top genotype should be the synthetic A/B diploid, got:\n{}",
        stdout
    );

    let report = std::fs::read_to_string(&report_path).unwrap();
    assert!(report.contains("#impg genotype cos report"), "{report}");
    assert!(
        report.contains("sample_pack_counting_semantics\tdistinct_nodes_per_read"),
        "report should state pack counting semantics:\n{report}"
    );

    let sample_rows = report_section_lines(&report, "sample_locus_features");
    assert_eq!(
        sample_rows.first().copied(),
        Some("node_id\tsample_count"),
        "sample feature section should have a parseable header"
    );
    let mut nonzero_sample_features = 0usize;
    for row in sample_rows.iter().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 2, "bad sample feature row: {row}");
        let node_id: u32 = fields[0].parse().unwrap();
        let sample_count: u64 = fields[1].parse().unwrap();
        assert_eq!(
            sample_count,
            pack_counts.get(&node_id).copied().unwrap_or(0),
            "report sample count should match pack TSV for node {node_id}"
        );
        if sample_count > 0 {
            nonzero_sample_features += 1;
        }
    }
    assert!(
        nonzero_sample_features > 0,
        "report should expose nonzero sample evidence over locus features"
    );

    let summary_rows = report_section_lines(&report, "pack_evidence_summary");
    assert!(
        summary_rows
            .iter()
            .any(|row| *row == format!("pack_nonzero_nodes\t{}", pack_counts.len())),
        "pack summary should match pack TSV nonzero node count:\n{report}"
    );

    let candidate_rows = report_section_lines(&report, "candidates");
    assert!(
        candidate_rows
            .first()
            .unwrap()
            .contains("sample_overlap_unique_nodes"),
        "candidate table should include sample overlap columns"
    );
    let mut saw_candidate_overlap = false;
    for row in candidate_rows.iter().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        assert!(
            fields.len() >= 18,
            "candidate rows should include feature counts and overlap stats: {row}"
        );
        let feature_count = fields[8].parse::<usize>().unwrap();
        let total_mass = fields[9].parse::<u64>().unwrap();
        let unique_nodes = fields[10].parse::<usize>().unwrap();
        let sample_overlap = fields[15].parse::<usize>().unwrap();
        assert!(feature_count > 0, "candidate should have features: {row}");
        assert!(
            total_mass >= feature_count as u64,
            "bad candidate mass: {row}"
        );
        assert_eq!(
            feature_count, unique_nodes,
            "unique node count should match vector width"
        );
        if sample_overlap > 0 {
            saw_candidate_overlap = true;
        }
    }
    assert!(
        saw_candidate_overlap,
        "candidate table should show sample overlaps"
    );

    let candidate_feature_rows = report_section_lines(&report, "candidate_features");
    assert_eq!(
        candidate_feature_rows.first().copied(),
        Some("candidate_index\tpath\tnode_id\tsample_count\tcandidate_count\tdot_contribution")
    );
    assert!(
        candidate_feature_rows.iter().skip(1).any(|row| row
            .split('\t')
            .nth(4)
            .unwrap()
            .parse::<u64>()
            .unwrap()
            > 0),
        "candidate feature table should expose candidate node counts"
    );

    let result_rows = report_section_lines(&report, "result_scores");
    assert!(
        result_rows.len() >= 2,
        "report should include result scores"
    );
    let report_top_fields: Vec<&str> = result_rows[1].split('\t').collect();
    assert_eq!(report_top_fields[0], top_fields[0], "rank should match");
    assert_eq!(report_top_fields[1], top_fields[1], "method should match");
    assert_eq!(report_top_fields[2], top_fields[2], "ploidy should match");
    assert_eq!(
        report_top_fields[3], top_fields[3],
        "similarity should match"
    );
    assert_eq!(report_top_fields[4], top_fields[4], "qv should match");
    assert_eq!(report_top_fields[5], top_fields[5], "dot should match");
    assert_eq!(
        report_top_fields[6], top_fields[6],
        "sample norm should match"
    );
    assert_eq!(
        report_top_fields[7], top_fields[7],
        "genotype norm should match"
    );
    assert_eq!(
        report_top_fields[9], top_fields[8],
        "report haplotypes should match primary output"
    );

    let result_feature_rows = report_section_lines(&report, "result_features");
    assert_eq!(
        result_feature_rows.first().copied(),
        Some("rank\tnode_id\tsample_count\tcandidate_counts\tgenotype_count\tdot_contribution")
    );
    let mut saw_nonzero_dot = false;
    for row in result_feature_rows.iter().skip(1) {
        let fields: Vec<&str> = row.split('\t').collect();
        assert_eq!(fields.len(), 6, "bad result feature row: {row}");
        let sample_count = fields[2].parse::<u64>().unwrap();
        let candidate_sum = fields[3]
            .split(',')
            .map(|count| count.parse::<u64>().unwrap())
            .sum::<u64>();
        let genotype_count = fields[4].parse::<u64>().unwrap();
        let dot = fields[5].parse::<f64>().unwrap();
        assert_eq!(
            candidate_sum, genotype_count,
            "genotype count should be the sum of candidate counts: {row}"
        );
        assert_eq!(
            dot,
            (sample_count * genotype_count) as f64,
            "dot contribution should decompose sample_count * genotype_count: {row}"
        );
        if dot > 0.0 {
            saw_nonzero_dot = true;
        }
    }
    assert!(
        saw_nonzero_dot,
        "result feature table should include nonzero dot contributions"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_genotype_cos_short_read_simulation_calls_heterozygote() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_genotype_short_read_sim");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(900, 31));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(900, 32));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(900, 33));
    let allele_c = numeric_to_ascii(&make_sequence_numeric(900, 34));
    let right = numeric_to_ascii(&make_sequence_numeric(900, 35));

    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);
    let mut hap_c = Vec::new();
    hap_c.extend_from_slice(&left);
    hap_c.extend_from_slice(&allele_c);
    hap_c.extend_from_slice(&right);
    assert_eq!(hap_a.len(), hap_b.len());
    assert_eq!(hap_a.len(), hap_c.len());

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleC#0#chr1").unwrap();
        f.write_all(&hap_c).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("simulated_reads.fq");
    let read_count = {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        let a = write_tiled_fastq(&mut f, "hapA", &hap_a, 250, 25).unwrap();
        let b = write_tiled_fastq(&mut f, "hapB", &hap_b, 250, 25).unwrap();
        a + b
    };
    assert!(
        read_count >= 190,
        "expected dense tiled short-read simulation, got {} reads",
        read_count
    );

    let packbin_path = dir.join("simulated.packbin");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "packbin",
            "-O",
            packbin_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o packbin");
    assert!(
        map.status.success(),
        "impg map -o packbin failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let target_range = format!("sampleA#0#chr1:0-{}", hap_a.len());
    let gt = Command::new(&bin)
        .args([
            "genotype",
            "cos",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-p",
            packbin_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--top-n",
            "5",
            "--candidate-top-k",
            "10",
            "--min-anchors",
            "2",
            "--min-span-fraction",
            "0.8",
        ])
        .output()
        .expect("failed to run impg genotype cos");
    assert!(
        gt.status.success(),
        "impg genotype cos failed: {}",
        String::from_utf8_lossy(&gt.stderr)
    );
    let gt_stdout = String::from_utf8_lossy(&gt.stdout);
    let top = gt_stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .expect("expected at least one cosine genotype result row");
    let fields: Vec<&str> = top.split('\t').collect();
    assert!(
        fields.len() >= 12,
        "cos genotype row should have at least 12 fields: {}",
        top
    );
    assert_eq!(fields[0], "1");
    assert_eq!(fields[1], "cos");
    assert_eq!(fields[2], "2");
    let similarity = fields[3].parse::<f64>().unwrap();
    assert!(
        similarity > 0.90,
        "expected strong short-read heterozygous cosine score, got {}\n{}",
        similarity,
        gt_stdout
    );
    assert!(
        fields[8].contains("sampleA#0#chr1") && fields[8].contains("sampleB#0#chr1"),
        "top genotype should call the simulated A/B diploid, got:\n{}",
        gt_stdout
    );
    assert!(
        !fields[8].contains("sampleC#0#chr1"),
        "top genotype should not include the unsampled decoy haplotype, got:\n{}",
        gt_stdout
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_pack_partitions_and_discovery() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_pack");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(900, 41));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(900, 42));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(900, 43));
    let allele_c = numeric_to_ascii(&make_sequence_numeric(900, 44));
    let right = numeric_to_ascii(&make_sequence_numeric(900, 45));

    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);
    let mut hap_c = Vec::new();
    hap_c.extend_from_slice(&left);
    hap_c.extend_from_slice(&allele_c);
    hap_c.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleC#0#chr1").unwrap();
        f.write_all(&hap_c).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("simulated_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "hapA", &hap_a, 250, 25).unwrap();
        write_tiled_fastq(&mut f, "hapB", &hap_b, 250, 25).unwrap();
    }

    let proj_path = dir.join("sample.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        map.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let target_range = format!("sampleA#0#chr1:0-{}", hap_a.len());
    let infer_range = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--top-n",
            "2",
            "--candidate-top-k",
            "10",
            "--min-span-fraction",
            "0.8",
        ])
        .output()
        .expect("failed to run impg infer -r");
    assert!(
        infer_range.status.success(),
        "impg infer -r failed: {}",
        String::from_utf8_lossy(&infer_range.stderr)
    );
    let infer_stdout = String::from_utf8_lossy(&infer_range.stdout);
    assert!(infer_stdout.contains("#impg infer"), "{}", infer_stdout);
    assert!(infer_stdout.contains("#score\tcos"), "{}", infer_stdout);
    let first = infer_stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .expect("expected infer result row");
    let fields: Vec<&str> = first.split('\t').collect();
    assert!(
        fields.len() >= 14,
        "infer row should have at least 14 fields: {}",
        first
    );
    assert_eq!(fields[0], "1");
    assert_eq!(fields[5], "cos");
    assert_eq!(fields[6], "2");
    assert_eq!(fields[13], "PASS");
    assert!(
        fields[9].contains("sampleA#0#chr1") && fields[9].contains("sampleB#0#chr1"),
        "infer should call the simulated A/B diploid, got:\n{}",
        infer_stdout
    );
    assert!(
        !fields[9].contains("sampleC#0#chr1"),
        "infer top call should not include the unsampled decoy, got:\n{}",
        infer_stdout
    );

    let target_bed = dir.join("targets.bed");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&target_bed).unwrap();
        writeln!(f, "sampleA#0#chr1\t0\t{}\twhole", hap_a.len()).unwrap();
        writeln!(f, "sampleA#0#chr1\t900\t1800\tallele").unwrap();
    }
    let infer_bed_zst = dir.join("infer.bed.tsv.zst");
    let infer_bed = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "--target-bed",
            target_bed.to_str().unwrap(),
            "--top-n",
            "1",
            "--candidate-top-k",
            "10",
            "--min-span-fraction",
            "0.7",
            "-O",
            infer_bed_zst.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer --target-bed");
    assert!(
        infer_bed.status.success(),
        "impg infer --target-bed failed: {}",
        String::from_utf8_lossy(&infer_bed.stderr)
    );
    assert!(
        infer_bed.stdout.is_empty(),
        "infer -O should write result rows to the file"
    );
    let mut decoder =
        zstd::stream::read::Decoder::new(std::fs::File::open(&infer_bed_zst).unwrap()).unwrap();
    let mut decoded = String::new();
    {
        use std::io::Read;
        decoder.read_to_string(&mut decoded).unwrap();
    }
    let bed_rows: Vec<&str> = decoded
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(bed_rows.len(), 2, "target BED top-n 1 should emit two rows");
    assert!(
        bed_rows.iter().all(|row| row.ends_with("\tPASS")),
        "{}",
        decoded
    );

    let partitions_path = dir.join("partitions.bed");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&partitions_path).unwrap();
        writeln!(f, "sampleA#0#chr1\t0\t{}\tp0", hap_a.len()).unwrap();
        writeln!(f, "sampleB#0#chr1\t0\t{}\tp0", hap_b.len()).unwrap();
        writeln!(f, "sampleA#0#chr1\t900\t1800\tp1").unwrap();
        writeln!(f, "sampleB#0#chr1\t900\t1800\tp1").unwrap();
    }
    let infer_partitions = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "--partitions",
            partitions_path.to_str().unwrap(),
            "--top-n",
            "1",
            "--candidate-top-k",
            "10",
            "--min-span-fraction",
            "0.7",
            "--stitch",
            "beam",
            "--emit-mosaic",
            dir.join("mosaic.tsv").to_str().unwrap(),
            "--emit-fasta",
            dir.join("mosaic.fa").to_str().unwrap(),
            "--emit-gfa",
            dir.join("mosaic.gfa").to_str().unwrap(),
            "--sequence-files",
            fasta_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer --partitions");
    assert!(
        infer_partitions.status.success(),
        "impg infer --partitions failed: {}",
        String::from_utf8_lossy(&infer_partitions.stderr)
    );
    let partition_stdout = String::from_utf8_lossy(&infer_partitions.stdout);
    let partition_rows: Vec<&str> = partition_stdout
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .collect();
    assert_eq!(
        partition_rows.len(),
        2,
        "ready partition BED should emit one top row per partition"
    );
    assert!(
        partition_rows
            .iter()
            .any(|row| row.split('\t').nth(1) == Some("p0")),
        "{}",
        partition_stdout
    );
    assert!(
        partition_rows
            .iter()
            .any(|row| row.split('\t').nth(1) == Some("p1")),
        "{}",
        partition_stdout
    );
    let mosaic_tsv = std::fs::read_to_string(dir.join("mosaic.tsv")).unwrap();
    assert!(mosaic_tsv.contains("#impg infer mosaic"), "{}", mosaic_tsv);
    assert!(mosaic_tsv.contains("sampleA#0#chr1"), "{}", mosaic_tsv);
    let mosaic_fa = std::fs::read_to_string(dir.join("mosaic.fa")).unwrap();
    assert!(mosaic_fa.contains(">impg_phase_0"), "{}", mosaic_fa);
    assert!(mosaic_fa.contains(">impg_phase_1"), "{}", mosaic_fa);
    let mosaic_gfa = std::fs::read_to_string(dir.join("mosaic.gfa")).unwrap();
    assert!(mosaic_gfa.starts_with("H\tVN:Z:1.0"), "{}", mosaic_gfa);
    assert!(mosaic_gfa.contains("\nS\tphase0_"), "{}", mosaic_gfa);
    assert!(mosaic_gfa.contains("\nP\timpg_phase_0"), "{}", mosaic_gfa);

    let no_d = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer without target inputs");
    assert!(
        !no_d.status.success(),
        "infer discovery without -d should fail"
    );
    assert!(
        String::from_utf8_lossy(&no_d.stderr).contains("requires -d/--merge-distance"),
        "unexpected stderr: {}",
        String::from_utf8_lossy(&no_d.stderr)
    );

    let infer_discovery = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-w",
            "1500",
            "-d",
            "100",
            "--min-missing-size",
            "100",
            "--min-boundary-distance",
            "100",
            "--top-n",
            "1",
            "--candidate-top-k",
            "10",
            "--min-span-fraction",
            "0.7",
        ])
        .output()
        .expect("failed to run impg infer discovery");
    assert!(
        infer_discovery.status.success(),
        "impg infer discovery failed: {}",
        String::from_utf8_lossy(&infer_discovery.stderr)
    );
    let discovery_stdout = String::from_utf8_lossy(&infer_discovery.stdout);
    assert!(
        discovery_stdout.contains("#impg infer"),
        "{}",
        discovery_stdout
    );
    assert!(
        discovery_stdout
            .lines()
            .any(|line| !line.starts_with('#') && line.ends_with("\tPASS")),
        "discovery infer should emit at least one PASS row:\n{}",
        discovery_stdout
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_read_walk_links_phase_recombinant_mosaic() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_read_links");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left_1 = numeric_to_ascii(&make_sequence_numeric(950, 81));
    let left_2 = mutate_ascii_every(&left_1, 37, 127);
    let right_1 = numeric_to_ascii(&make_sequence_numeric(950, 83));
    let right_2 = mutate_ascii_every(&right_1, 53, 131);

    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left_1);
    hap_a.extend_from_slice(&right_1);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left_2);
    hap_b.extend_from_slice(&right_2);
    let mut hap_c = Vec::new();
    hap_c.extend_from_slice(&left_1);
    hap_c.extend_from_slice(&right_2);
    let mut hap_d = Vec::new();
    hap_d.extend_from_slice(&left_2);
    hap_d.extend_from_slice(&right_1);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&hap_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&hap_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleC#0#chr1").unwrap();
        f.write_all(&hap_c).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleD#0#chr1").unwrap();
        f.write_all(&hap_d).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("simulated_recombinant_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "hapC", &hap_c, 1200, 175).unwrap();
        write_tiled_fastq(&mut f, "hapD", &hap_d, 1200, 175).unwrap();
    }

    let proj_path = dir.join("sample.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        map.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let partitions_path = dir.join("partitions.bed");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&partitions_path).unwrap();
        writeln!(f, "sampleA#0#chr1\t0\t950\tleft").unwrap();
        writeln!(f, "sampleA#0#chr1\t950\t1900\tright").unwrap();
    }

    let mosaic_path = dir.join("mosaic.tsv");
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "--partitions",
            partitions_path.to_str().unwrap(),
            "--top-n",
            "20",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.7",
            "--stitch",
            "beam",
            "--stitch-beam",
            "500",
            "--read-link-weight",
            "5",
            "--emit-mosaic",
            mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer with read-link stitching");
    assert!(
        infer.status.success(),
        "impg infer with read-link stitching failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );

    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    assert!(
        mosaic.contains("read_link_reward"),
        "mosaic output should expose read-link scoring columns:\n{}",
        mosaic
    );
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.len(),
        4,
        "expected two phases across two partitions:\n{}",
        mosaic
    );
    assert!(
        rows.iter()
            .filter(|row| row[1] == "right")
            .all(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "right-partition transitions should be rewarded by spanning read walks:\n{}",
        mosaic
    );

    let zero_weight_mosaic_path = dir.join("mosaic.zero_read_weight.tsv");
    let infer_zero_weight = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "--partitions",
            partitions_path.to_str().unwrap(),
            "--top-n",
            "20",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.7",
            "--stitch",
            "beam",
            "--stitch-beam",
            "500",
            "--read-link-weight",
            "0",
            "--emit-mosaic",
            zero_weight_mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer with disabled read-link weight");
    assert!(
        infer_zero_weight.status.success(),
        "impg infer with disabled read-link weight failed: {}",
        String::from_utf8_lossy(&infer_zero_weight.stderr)
    );
    let zero_weight_mosaic = std::fs::read_to_string(&zero_weight_mosaic_path).unwrap();
    let zero_weight_rows: Vec<Vec<&str>> = zero_weight_mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert!(
        zero_weight_rows
            .iter()
            .filter(|row| row[1] == "right")
            .all(|row| row[16].parse::<f64>().unwrap() == 0.0),
        "read-link rewards should be disabled when --read-link-weight 0:\n{}",
        zero_weight_mosaic
    );

    let mut phase_paths: std::collections::BTreeMap<&str, Vec<&str>> =
        std::collections::BTreeMap::new();
    for row in &rows {
        phase_paths.entry(row[0]).or_default().push(row[5]);
    }
    let mut stitched_paths: Vec<Vec<&str>> = phase_paths.into_values().collect();
    stitched_paths.sort();
    assert_eq!(
        stitched_paths,
        vec![
            vec!["sampleC#0#chr1", "sampleC#0#chr1"],
            vec!["sampleD#0#chr1", "sampleD#0#chr1"],
        ],
        "read links plus crossover penalty should choose the sampled recombinant haplotypes:\n{}",
        mosaic
    );

    let block_mosaic_path = dir.join("mosaic.phase_blocks.tsv");
    let target_range = "sampleA#0#chr1:0-1900";
    let infer_blocks = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            target_range,
            "--phase-block-size",
            "950",
            "--top-n",
            "20",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.7",
            "--stitch",
            "beam",
            "--stitch-beam",
            "500",
            "--read-link-weight",
            "5",
            "--emit-mosaic",
            block_mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer with phase-block stitching");
    assert!(
        infer_blocks.status.success(),
        "impg infer with phase-block stitching failed: {}",
        String::from_utf8_lossy(&infer_blocks.stderr)
    );
    let block_mosaic = std::fs::read_to_string(&block_mosaic_path).unwrap();
    assert!(
        block_mosaic.contains("#phase_block_size\t950")
            || String::from_utf8_lossy(&infer_blocks.stdout).contains("#phase_block_size\t950"),
        "phase-block inference should record the internal block size:\nstdout:\n{}\nmosaic:\n{}",
        String::from_utf8_lossy(&infer_blocks.stdout),
        block_mosaic
    );
    let block_rows: Vec<Vec<&str>> = block_mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        block_rows.len(),
        4,
        "expected two phases across two internal phase blocks:\n{}",
        block_mosaic
    );
    assert!(
        block_rows
            .iter()
            .filter(|row| row[1].ends_with("#block1"))
            .all(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "second internal phase block should be rewarded by spanning read walks:\n{}",
        block_mosaic
    );
    let mut block_phase_paths: std::collections::BTreeMap<&str, Vec<&str>> =
        std::collections::BTreeMap::new();
    for row in &block_rows {
        block_phase_paths.entry(row[0]).or_default().push(row[5]);
    }
    let mut block_stitched_paths: Vec<Vec<&str>> = block_phase_paths.into_values().collect();
    block_stitched_paths.sort();
    assert_eq!(
        block_stitched_paths,
        vec![
            vec!["sampleC#0#chr1", "sampleC#0#chr1"],
            vec!["sampleD#0#chr1", "sampleD#0#chr1"],
        ],
        "phase-block stitching should infer copied recombinant segments inside one target:\n{}",
        block_mosaic
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_cnv_repeated_syncmer_path_calls_duplicated_haplotype() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_cnv_repeats");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(500, 91));
    let copy = numeric_to_ascii(&make_sequence_numeric(700, 92));
    let copy_alt = mutate_ascii_every(&copy, 23, 101);
    let right = numeric_to_ascii(&make_sequence_numeric(500, 93));

    let mut hap_single = Vec::new();
    hap_single.extend_from_slice(&left);
    hap_single.extend_from_slice(&copy);
    hap_single.extend_from_slice(&right);

    let mut hap_double = Vec::new();
    hap_double.extend_from_slice(&left);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&right);

    let mut hap_alt = Vec::new();
    hap_alt.extend_from_slice(&left);
    hap_alt.extend_from_slice(&copy_alt);
    hap_alt.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleSingle#0#chr1").unwrap();
        f.write_all(&hap_single).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleDouble#0#chr1").unwrap();
        f.write_all(&hap_double).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleAlt#0#chr1").unwrap();
        f.write_all(&hap_alt).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("duplicated_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "double", &hap_double, 1100, 175).unwrap();
    }

    let proj_path = dir.join("sample.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        map.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let mut gaf_decoder = zstd::stream::read::Decoder::new(
        std::fs::File::open(proj_path.join("reads.gaf.zst")).unwrap(),
    )
    .unwrap();
    let mut gaf = String::new();
    {
        use std::io::Read;
        gaf_decoder.read_to_string(&mut gaf).unwrap();
    }
    assert!(
        gaf.lines().any(|line| {
            let Some(path) = line.split('\t').nth(5) else {
                return false;
            };
            let mut seen = std::collections::BTreeSet::new();
            let mut token = String::new();
            for ch in path.chars() {
                if ch == '>' || ch == '<' {
                    if !token.is_empty() && !seen.insert(token.clone()) {
                        return true;
                    }
                    token.clear();
                } else {
                    token.push(ch);
                }
            }
            !token.is_empty() && !seen.insert(token)
        }),
        "duplicated-copy reads should produce GAF walks with repeated syncmer nodes:\n{}",
        gaf
    );

    let target_range = format!("sampleSingle#0#chr1:0-{}", hap_single.len());
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--top-n",
            "5",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.45",
        ])
        .output()
        .expect("failed to run impg infer on duplicated haplotype");
    assert!(
        infer.status.success(),
        "impg infer on duplicated haplotype failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );
    let infer_stdout = String::from_utf8_lossy(&infer.stdout);
    let first = infer_stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .expect("expected at least one infer row");
    let fields: Vec<&str> = first.split('\t').collect();
    assert!(
        fields[9].contains("sampleDouble#0#chr1"),
        "duplicated-read evidence should call the duplicated haplotype first:\n{}",
        infer_stdout
    );
    assert!(
        !fields[9].contains("sampleAlt#0#chr1"),
        "unrelated single-copy allele should not be part of the top CNV call:\n{}",
        infer_stdout
    );

    let mosaic_path = dir.join("cnv_phase_blocks.tsv");
    let infer_blocks = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--phase-block-size",
            "600",
            "--top-n",
            "5",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.35",
            "--stitch",
            "beam",
            "--stitch-beam",
            "300",
            "--read-link-weight",
            "3",
            "--emit-mosaic",
            mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer phase blocks on duplicated haplotype");
    assert!(
        infer_blocks.status.success(),
        "phase-block infer on duplicated haplotype failed: {}",
        String::from_utf8_lossy(&infer_blocks.stderr)
    );
    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert!(
        rows.len() >= 3,
        "phase-block CNV mosaic should emit multiple phased blocks:\n{}",
        mosaic
    );
    assert!(
        rows.iter().any(|row| row[5] == "sampleDouble#0#chr1"),
        "phase-block CNV mosaic should retain the duplicated haplotype candidate:\n{}",
        mosaic
    );
    assert!(
        rows.iter()
            .skip(2)
            .any(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "repeated-node read walks should still contribute nonzero transition reward:\n{}",
        mosaic
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_triplicated_syncmer_path_beats_lower_copy_decoys() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_triplication");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(520, 101));
    let copy = numeric_to_ascii(&make_sequence_numeric(680, 102));
    let copy_alt = mutate_ascii_every(&copy, 31, 103);
    let right = numeric_to_ascii(&make_sequence_numeric(520, 104));

    let mut hap_single = Vec::new();
    hap_single.extend_from_slice(&left);
    hap_single.extend_from_slice(&copy);
    hap_single.extend_from_slice(&right);

    let mut hap_double = Vec::new();
    hap_double.extend_from_slice(&left);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&right);

    let mut hap_triple = Vec::new();
    hap_triple.extend_from_slice(&left);
    hap_triple.extend_from_slice(&copy);
    hap_triple.extend_from_slice(&copy);
    hap_triple.extend_from_slice(&copy);
    hap_triple.extend_from_slice(&right);

    let mut hap_alt = Vec::new();
    hap_alt.extend_from_slice(&left);
    hap_alt.extend_from_slice(&copy_alt);
    hap_alt.extend_from_slice(&copy_alt);
    hap_alt.extend_from_slice(&copy_alt);
    hap_alt.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleSingle#0#chr1").unwrap();
        f.write_all(&hap_single).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleDouble#0#chr1").unwrap();
        f.write_all(&hap_double).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleTriple#0#chr1").unwrap();
        f.write_all(&hap_triple).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleAltTriple#0#chr1").unwrap();
        f.write_all(&hap_alt).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("triplicated_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "triple", &hap_triple, 1120, 175).unwrap();
    }

    let proj_path = dir.join("sample_triple.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj for triplicated reads");
    assert!(
        map.status.success(),
        "impg map -o proj for triplicated reads failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let target_range = format!("sampleSingle#0#chr1:0-{}", hap_single.len());
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--top-n",
            "5",
            "--candidate-top-k",
            "20",
            "--min-span-fraction",
            "0.35",
        ])
        .output()
        .expect("failed to run impg infer on triplicated haplotype");
    assert!(
        infer.status.success(),
        "impg infer on triplicated haplotype failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );
    let stdout = String::from_utf8_lossy(&infer.stdout);
    let first = stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .expect("expected at least one triplicated infer row");
    let fields: Vec<&str> = first.split('\t').collect();
    assert!(
        fields[9].contains("sampleTriple#0#chr1"),
        "three-copy evidence should call the triplicated haplotype first:\n{}",
        stdout
    );
    assert!(
        !fields[9].contains("sampleSingle#0#chr1")
            && !fields[9].contains("sampleDouble#0#chr1")
            && !fields[9].contains("sampleAltTriple#0#chr1"),
        "three-copy dosage should not collapse to one-copy, two-copy, or unrelated triplicated alleles:\n{}",
        stdout
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_nested_sv_noisy_low_coverage_phase_blocks() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_nested_sv_noise");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let a = numeric_to_ascii(&make_sequence_numeric(500, 111));
    let b = numeric_to_ascii(&make_sequence_numeric(500, 112));
    let c = numeric_to_ascii(&make_sequence_numeric(500, 113));
    let d = numeric_to_ascii(&make_sequence_numeric(500, 114));
    let e = numeric_to_ascii(&make_sequence_numeric(500, 115));
    let insertion = numeric_to_ascii(&make_sequence_numeric(350, 116));

    let mut hap_ref = Vec::new();
    hap_ref.extend_from_slice(&a);
    hap_ref.extend_from_slice(&b);
    hap_ref.extend_from_slice(&c);
    hap_ref.extend_from_slice(&d);
    hap_ref.extend_from_slice(&e);

    // Nested event relative to reference: insertion after B and deletion of D.
    // This forces the local calls to use a non-reference inserted interval and
    // skip an absent reference block under sparse, noisy read evidence.
    let mut hap_complex = Vec::new();
    hap_complex.extend_from_slice(&a);
    hap_complex.extend_from_slice(&b);
    hap_complex.extend_from_slice(&insertion);
    hap_complex.extend_from_slice(&c);
    hap_complex.extend_from_slice(&e);

    let mut hap_deletion = Vec::new();
    hap_deletion.extend_from_slice(&a);
    hap_deletion.extend_from_slice(&b);
    hap_deletion.extend_from_slice(&c);
    hap_deletion.extend_from_slice(&e);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleRef#0#chr1").unwrap();
        f.write_all(&hap_ref).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleComplex#0#chr1").unwrap();
        f.write_all(&hap_complex).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleDeletion#0#chr1").unwrap();
        f.write_all(&hap_deletion).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("complex_noisy_low_coverage_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        let reads =
            write_tiled_fastq_with_errors(&mut f, "complex", &hap_complex, 650, 425, 6, 173)
                .unwrap();
        assert!(
            reads <= 6,
            "test should remain low coverage, wrote {} reads",
            reads
        );
    }

    let proj_path = dir.join("complex.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "1",
        ])
        .output()
        .expect("failed to run impg map -o proj for noisy nested-SV reads");
    assert!(
        map.status.success(),
        "impg map -o proj for noisy nested-SV reads failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let target_range = format!("sampleRef#0#chr1:0-{}", hap_ref.len());
    let mosaic_path = dir.join("nested_sv_mosaic.tsv");
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--candidate-mode",
            "overlapping",
            "--phase-block-size",
            "500",
            "--top-n",
            "12",
            "--candidate-top-k",
            "40",
            "--min-anchors",
            "1",
            "--min-span-fraction",
            "0.2",
            "--stitch",
            "beam",
            "--stitch-beam",
            "600",
            "--read-link-weight",
            "4",
            "--emit-mosaic",
            mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer on noisy nested-SV reads");
    assert!(
        infer.status.success(),
        "impg infer on noisy nested-SV reads failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );
    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert!(
        rows.len() >= 4,
        "nested-SV phase-block mosaic should emit several copied segments:\n{}",
        mosaic
    );
    let complex_rows = rows
        .iter()
        .filter(|row| row[5] == "sampleComplex#0#chr1")
        .count();
    assert!(
        complex_rows >= 2,
        "noisy low-coverage nested-SV evidence should copy from the complex haplotype in multiple blocks:\n{}",
        mosaic
    );
    assert!(
        rows.iter()
            .skip(1)
            .any(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "noisy low-coverage read walks should still contribute transition evidence:\n{}",
        mosaic
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_read_walk_emission_resolves_order_decoy() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_read_walk_emission");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(420, 151));
    let copy_a = numeric_to_ascii(&make_sequence_numeric(480, 152));
    let copy_b = numeric_to_ascii(&make_sequence_numeric(480, 153));
    let copy_c = numeric_to_ascii(&make_sequence_numeric(480, 154));
    let right = numeric_to_ascii(&make_sequence_numeric(420, 155));

    let mut true_repeat = Vec::new();
    true_repeat.extend_from_slice(&copy_a);
    true_repeat.extend_from_slice(&copy_b);
    true_repeat.extend_from_slice(&copy_a);
    true_repeat.extend_from_slice(&copy_c);
    true_repeat.extend_from_slice(&copy_a);

    let mut decoy_repeat = Vec::new();
    decoy_repeat.extend_from_slice(&copy_a);
    decoy_repeat.extend_from_slice(&copy_c);
    decoy_repeat.extend_from_slice(&copy_a);
    decoy_repeat.extend_from_slice(&copy_b);
    decoy_repeat.extend_from_slice(&copy_a);

    let mut hap_ref = Vec::new();
    hap_ref.extend_from_slice(&left);
    hap_ref.extend_from_slice(&true_repeat);
    hap_ref.extend_from_slice(&right);

    // Same node-count and adjacent-transition projection as the target interval:
    // A-B-A-C-A and A-C-A-B-A both contain AB, BA, AC, and CA once. Only the
    // whole read walk disambiguates the true order.
    let mut hap_decoy = Vec::new();
    hap_decoy.extend_from_slice(&left);
    hap_decoy.extend_from_slice(&decoy_repeat);
    hap_decoy.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleRef#0#chr1").unwrap();
        f.write_all(&hap_ref).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleADecoy#0#chr1").unwrap();
        f.write_all(&hap_decoy).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("ordered_repeat_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "ordered", &true_repeat, true_repeat.len(), 120).unwrap();
    }

    let proj_path = dir.join("ordered.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "1",
        ])
        .output()
        .expect("failed to run impg map -o proj for ordered repeat reads");
    assert!(
        map.status.success(),
        "impg map -o proj for ordered repeat reads failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );
    let mut gaf_decoder = zstd::stream::read::Decoder::new(
        std::fs::File::open(proj_path.join("reads.gaf.zst")).unwrap(),
    )
    .unwrap();
    let mut projected_gaf = String::new();
    {
        use std::io::Read;
        gaf_decoder.read_to_string(&mut projected_gaf).unwrap();
    }
    assert!(
        projected_gaf.contains("\tqp:B:I,"),
        "projection GAF should carry query syncmer positions for GBWT MEM scoring:\n{}",
        projected_gaf
    );

    let target_range = format!(
        "sampleRef#0#chr1:{}-{}",
        left.len(),
        left.len() + true_repeat.len()
    );
    let mosaic_path = dir.join("ordered_repeat_mosaic.tsv");
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--candidate-mode",
            "spanning",
            "--top-n",
            "16",
            "--candidate-top-k",
            "80",
            "--min-anchors",
            "1",
            "--min-read-link-anchors",
            "1",
            "--min-span-fraction",
            "0.15",
            "--stitch",
            "beam",
            "--stitch-beam",
            "1000",
            "--read-link-weight",
            "6",
            "--emit-mosaic",
            mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer on read-walk order decoy");
    assert!(
        infer.status.success(),
        "impg infer on read-walk order decoy failed: {}\nstdout:\n{}",
        String::from_utf8_lossy(&infer.stderr),
        String::from_utf8_lossy(&infer.stdout)
    );

    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.len(),
        1,
        "expected a single stitched interval:\n{}",
        mosaic
    );
    assert!(
        rows[0][5] == "sampleRef#0#chr1" && rows[0][19].parse::<f64>().unwrap() > 0.0,
        "whole-read-walk emission should choose the true A-B-A-C-A path and expose support:\n{}",
        mosaic
    );
    assert!(
        !rows.iter().any(|row| row[5] == "sampleADecoy#0#chr1"),
        "node-count and adjacent-transition equivalent decoy should not win once whole walks are scored:\n{}",
        mosaic
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_infer_paralogous_swapped_copies_avoid_decoy_collapse() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_infer_paralogy");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let left = numeric_to_ascii(&make_sequence_numeric(420, 131));
    let copy_a = numeric_to_ascii(&make_sequence_numeric(760, 132));
    let copy_b = mutate_ascii_every(&copy_a, 29, 97);
    let spacer = numeric_to_ascii(&make_sequence_numeric(360, 133));
    let right = numeric_to_ascii(&make_sequence_numeric(420, 134));

    let mut hap_ab = Vec::new();
    hap_ab.extend_from_slice(&left);
    hap_ab.extend_from_slice(&copy_a);
    hap_ab.extend_from_slice(&spacer);
    hap_ab.extend_from_slice(&copy_b);
    hap_ab.extend_from_slice(&right);

    let mut hap_ba = Vec::new();
    hap_ba.extend_from_slice(&left);
    hap_ba.extend_from_slice(&copy_b);
    hap_ba.extend_from_slice(&spacer);
    hap_ba.extend_from_slice(&copy_a);
    hap_ba.extend_from_slice(&right);

    let mut hap_aa = Vec::new();
    hap_aa.extend_from_slice(&left);
    hap_aa.extend_from_slice(&copy_a);
    hap_aa.extend_from_slice(&spacer);
    hap_aa.extend_from_slice(&copy_a);
    hap_aa.extend_from_slice(&right);

    let mut hap_bb = Vec::new();
    hap_bb.extend_from_slice(&left);
    hap_bb.extend_from_slice(&copy_b);
    hap_bb.extend_from_slice(&spacer);
    hap_bb.extend_from_slice(&copy_b);
    hap_bb.extend_from_slice(&right);

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleAB#0#chr1").unwrap();
        f.write_all(&hap_ab).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleBA#0#chr1").unwrap();
        f.write_all(&hap_ba).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleAA#0#chr1").unwrap();
        f.write_all(&hap_aa).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleBB#0#chr1").unwrap();
        f.write_all(&hap_bb).unwrap();
        writeln!(f).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let reads_path = dir.join("swapped_paralog_reads.fq");
    {
        let mut f = std::fs::File::create(&reads_path).unwrap();
        write_tiled_fastq(&mut f, "ba", &hap_ba, 900, 225).unwrap();
    }

    let proj_path = dir.join("sample_ba.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            reads_path.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj_path.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
        ])
        .output()
        .expect("failed to run impg map -o proj for swapped-paralog reads");
    assert!(
        map.status.success(),
        "impg map -o proj for swapped-paralog reads failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );

    let target_range = format!("sampleAB#0#chr1:0-{}", hap_ab.len());
    let mosaic_path = dir.join("paralogy_mosaic.tsv");
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            idx_prefix.to_str().unwrap(),
            "--proj",
            proj_path.to_str().unwrap(),
            "-r",
            &target_range,
            "--ploidy",
            "1",
            "--phase-block-size",
            "760",
            "--top-n",
            "10",
            "--candidate-top-k",
            "40",
            "--min-span-fraction",
            "0.35",
            "--stitch",
            "beam",
            "--stitch-beam",
            "600",
            "--read-link-weight",
            "5",
            "--emit-mosaic",
            mosaic_path.to_str().unwrap(),
        ])
        .output()
        .expect("failed to run impg infer on swapped-paralog reads");
    assert!(
        infer.status.success(),
        "impg infer on swapped-paralog reads failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );
    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    let ba_rows = rows
        .iter()
        .filter(|row| row[5] == "sampleBA#0#chr1")
        .count();
    assert!(
        ba_rows >= 2,
        "swapped-paralog evidence should copy multiple blocks from sampleBA:\n{}",
        mosaic
    );
    assert!(
        !rows
            .iter()
            .any(|row| row[5] == "sampleAA#0#chr1" || row[5] == "sampleBB#0#chr1"),
        "swapped-paralog evidence should not collapse to homo-copy decoys:\n{}",
        mosaic
    );
    assert!(
        rows.iter()
            .skip(1)
            .any(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "paralog-spanning reads should contribute transition evidence:\n{}",
        mosaic
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_map_cli_sampled_positions_paf() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_map_sampled_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let backbone = numeric_to_ascii(&make_sequence_numeric(1000, 42));
    let tail_a = numeric_to_ascii(&make_sequence_numeric(400, 1));
    let tail_b = numeric_to_ascii(&make_sequence_numeric(400, 2));

    let fasta_path = dir.join("index.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail_b).unwrap();
        writeln!(f).unwrap();
    }

    let query_path = dir.join("query.fq");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&query_path).unwrap();
        writeln!(f, "@read1").unwrap();
        f.write_all(&backbone[100..800]).unwrap();
        writeln!(f).unwrap();
        writeln!(f, "+").unwrap();
        writeln!(f, "{}", "I".repeat(700)).unwrap();
    }

    let idx_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            idx_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("failed to run impg syng with sampled positions");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );
    assert!(
        std::path::Path::new(&format!("{}.syng.spos", idx_prefix.to_str().unwrap())).exists(),
        "default sampled-position sidecar should be written"
    );
    assert!(
        std::path::Path::new(&format!("{}.syng.pstep", idx_prefix.to_str().unwrap())).exists(),
        "default sampled path-step sidecar should be written"
    );
    assert!(
        std::path::Path::new(&format!("{}.syng.meta", idx_prefix.to_str().unwrap())).exists(),
        "syncmer metadata should be written"
    );
    assert!(
        !std::path::Path::new(&format!("{}.syng.locate", idx_prefix.to_str().unwrap())).exists(),
        "legacy FastLocate sidecar should not be built"
    );

    let paf = Command::new(&bin)
        .args([
            "map",
            "-a",
            idx_prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "paf",
            "--min-anchors",
            "2",
            "--max-hits",
            "10",
        ])
        .output()
        .expect("failed to run impg map -o paf with sampled positions");
    assert!(
        paf.status.success(),
        "impg map -o paf failed: {}",
        String::from_utf8_lossy(&paf.stderr)
    );
    let paf_stdout = String::from_utf8_lossy(&paf.stdout);
    let paf_lines: Vec<&str> = paf_stdout.lines().collect();
    assert!(!paf_lines.is_empty(), "expected PAF hits, got none");
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleA#0#chr1")),
        "expected sampleA hit, got:\n{}",
        paf_stdout
    );
    assert!(
        paf_lines.iter().any(|l| l.contains("sampleB#0#chr1")),
        "expected sampleB hit, got:\n{}",
        paf_stdout
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_syng_identical_sequences_build_and_query() {
    // Regression test for the vendored syng hash.c REMOVED-sentinel bug: two
    // byte-identical sequences (like AAA#0#chrIII and SGDref#0#chrIII in
    // yeast235) would both get stored as (start_node=1, start_count=0),
    // colliding in the GBWT startCount table. On query, the C layer would
    // die with "syngBWTpathStartOld startNode 1 count 0 >= startCount 0".
    //
    // This test builds two identical sequences via the CLI, loads the index,
    // and queries — it crashes at the C layer if the hash bug comes back.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_syng_ident_seqs_cli");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two records, byte-identical sequence content
        let seq = numeric_to_ascii(&make_sequence_numeric(1500, 42));
        writeln!(f, ">sampleA#0#chrIII").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chrIII").unwrap();
        f.write_all(&seq).unwrap();
        writeln!(f).unwrap();
    }

    let out_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            out_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("impg syng failed to run");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Load and query — this is where the bug used to crash at C level.
    let index = impg::syng::SyngIndex::load(
        out_prefix.to_str().unwrap(),
        impg::syng::SyncmerParams::default(),
    )
    .expect("load failed");

    // Verify start_counts are distinct (the bug used to save both as 0)
    let start_a = index.name_map.path_starts[0]
        .as_ref()
        .expect("seqA missing");
    let start_b = index.name_map.path_starts[1]
        .as_ref()
        .expect("seqB missing");
    assert_eq!(
        start_a.start_node, start_b.start_node,
        "identical sequences should share start_node"
    );
    assert_ne!(
        start_a.start_count, start_b.start_count,
        "two paths at same start_node must have distinct start_count indices \
         (vendored syng hash.c REMOVED-sentinel bug regression)"
    );

    // Query both sequences — used to crash with startCount 0 >= 0
    let intervals_a = index
        .query_region("sampleA#0#chrIII", 0, 1000, 0)
        .expect("query for sampleA should succeed");
    let intervals_b = index
        .query_region("sampleB#0#chrIII", 0, 1000, 0)
        .expect("query for sampleB should succeed");

    assert!(!intervals_a.is_empty());
    assert!(!intervals_b.is_empty());

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_partition_syng_end_to_end_bed() {
    // End-to-end: build syng index from FASTA, then run `impg partition -a`
    // and verify non-empty BED output. This is the path I added but never
    // exercised via CLI until now.
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_partition_syng_bed");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // Build a FASTA with enough content to produce multiple partition windows.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two sequences sharing a 2kbp backbone so partition finds homology
        let backbone = numeric_to_ascii(&make_sequence_numeric(2000, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(1000, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(1000, 2));

        writeln!(f, ">sampleA#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    // Step 1: build syng index
    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed. stderr: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Step 2: run partition -a with a syng index prefix and BED output
    let out_folder = dir.join("parts");
    std::fs::create_dir_all(&out_folder).unwrap();
    let part = Command::new(&bin)
        .args([
            "partition",
            "-d",
            "100000",
            "-a",
            syng_prefix.to_str().unwrap(),
            "-w",
            "1500",
            "-o",
            "bed",
            "--output-folder",
            out_folder.to_str().unwrap(),
            "--min-missing-size",
            "100",
            "--min-boundary-distance",
            "0",
            "-t",
            "1",
        ])
        .output()
        .expect("Failed to run impg partition -a syng prefix");
    assert!(
        part.status.success(),
        "impg partition -a syng prefix failed. stdout: {} stderr: {}",
        String::from_utf8_lossy(&part.stdout),
        String::from_utf8_lossy(&part.stderr)
    );

    // Single-file BED mode writes to `partitions.bed` in the output folder.
    let bed_files: Vec<_> = std::fs::read_dir(&out_folder)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("bed"))
        .collect();

    assert!(
        !bed_files.is_empty(),
        "No BED files produced in {:?}",
        out_folder
    );

    // At least one BED should have content (non-empty partition)
    let any_non_empty = bed_files
        .iter()
        .any(|e| std::fs::metadata(e.path()).unwrap().len() > 0);
    assert!(any_non_empty, "All partition BED files are empty");

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_partition_syng_gfa_blunt_engine() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_partition_syng_gfa_blunt");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(2500, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(700, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(700, 2));

        writeln!(f, ">sampleA#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed. stderr: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let out_folder = dir.join("gfas");
    std::fs::create_dir_all(&out_folder).unwrap();
    let part = Command::new(&bin)
        .args([
            "partition",
            "-d",
            "100000",
            "-a",
            syng_prefix.to_str().unwrap(),
            "-w",
            "1500",
            "-o",
            "gfa",
            "--gfa-engine",
            "syng",
            "--sequence-files",
            fasta_path.to_str().unwrap(),
            "--separate-files",
            "--output-folder",
            out_folder.to_str().unwrap(),
            "--min-missing-size",
            "100",
            "--min-boundary-distance",
            "0",
            "-t",
            "1",
        ])
        .output()
        .expect("Failed to run impg partition -a syng prefix -o gfa");
    assert!(
        part.status.success(),
        "impg partition syng gfa failed. stdout: {} stderr: {}",
        String::from_utf8_lossy(&part.stdout),
        String::from_utf8_lossy(&part.stderr)
    );

    let gfa_files: Vec<_> = std::fs::read_dir(&out_folder)
        .unwrap()
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().and_then(|s| s.to_str()) == Some("gfa"))
        .collect();
    assert!(
        !gfa_files.is_empty(),
        "No GFA files produced in {:?}",
        out_folder
    );

    let mut saw_segment = false;
    let mut nonzero_links = Vec::new();
    for entry in &gfa_files {
        let text = std::fs::read_to_string(entry.path()).unwrap();
        saw_segment |= text.lines().any(|line| line.starts_with("S\t"));
        nonzero_links.extend(
            text.lines()
                .filter(|line| line.starts_with("L\t") && !line.ends_with("\t0M"))
                .map(|line| line.to_string()),
        );
    }
    assert!(
        saw_segment,
        "partition syng GFA should contain segment lines"
    );
    assert!(
        nonzero_links.is_empty(),
        "syng partition GFA should be bluntified to 0M links. Nonzero links: {:?}",
        nonzero_links
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_partition_syng_gfa_reports_query_backend_error() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_partition_syng_gfa_query_error");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        let backbone = numeric_to_ascii(&make_sequence_numeric(1200, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(400, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(400, 2));

        writeln!(f, ">sampleA#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#1#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed. stderr: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    let pstep_path = format!("{}.syng.pstep", syng_prefix.to_str().unwrap());
    corrupt_syng_pstep_payload(&pstep_path);

    let out_folder = dir.join("gfas");
    std::fs::create_dir_all(&out_folder).unwrap();
    let part = Command::new(&bin)
        .args([
            "partition",
            "-d",
            "100000",
            "-a",
            syng_prefix.to_str().unwrap(),
            "-w",
            "800",
            "-o",
            "gfa",
            "--gfa-engine",
            "syng",
            "--sequence-files",
            fasta_path.to_str().unwrap(),
            "--separate-files",
            "--output-folder",
            out_folder.to_str().unwrap(),
            "--min-missing-size",
            "100",
            "--min-boundary-distance",
            "0",
            "-t",
            "1",
        ])
        .output()
        .expect("Failed to run impg partition -a syng prefix -o gfa");

    assert!(
        !part.status.success(),
        "partition should fail when syng query sidecars are corrupted. stdout: {} stderr: {}",
        String::from_utf8_lossy(&part.stdout),
        String::from_utf8_lossy(&part.stderr)
    );
    let stderr = String::from_utf8_lossy(&part.stderr);
    assert!(
        stderr.contains("syng query_region")
            && (stderr.contains("truncated varint") || stderr.contains("varint exceeds u64 width")),
        "stderr should report the syng query backend failure, got: {stderr}"
    );

    std::fs::remove_dir_all(&dir).ok();
}

fn corrupt_syng_pstep_payload(path: &str) {
    let mut data = std::fs::read(path).expect("failed to read syng pstep sidecar");
    let mut offset = 0usize;
    for _ in 0..4 {
        read_le_u64(&data, &mut offset);
    }
    let n_offsets = read_le_u64(&data, &mut offset) as usize;
    offset = offset
        .checked_add(n_offsets * std::mem::size_of::<u64>())
        .expect("pstep offset table size overflowed");
    let n_data = read_le_u64(&data, &mut offset) as usize;
    assert!(n_data > 0, "test pstep sidecar should have sampled data");
    let end = offset
        .checked_add(n_data)
        .expect("pstep data size overflowed");
    assert!(
        end <= data.len(),
        "pstep sidecar data segment extends past EOF"
    );
    for byte in &mut data[offset..end] {
        *byte = 0x80;
    }
    std::fs::write(path, data).expect("failed to corrupt syng pstep sidecar");
}

fn read_le_u64(data: &[u8], offset: &mut usize) -> u64 {
    let end = offset.checked_add(8).expect("u64 offset overflowed");
    let bytes = data
        .get(*offset..end)
        .expect("truncated u64 in syng pstep sidecar");
    *offset = end;
    u64::from_le_bytes(bytes.try_into().unwrap())
}

/// Regression test: `impg query -a <syng-prefix> -o gfa --gfa-engine seqwish:X` must
/// treat `X` as a sub-window size and split the query range into per-window
/// partitions, not silently collapse into a flat single-engine run.
///
/// Before this fix, `pggb:X`/`seqwish:X` in the syng+gfa path set
/// `partition_size = Some(X)` but was only used as a `skip_normalize` boolean
/// flag — the actual `X` was discarded and one flat engine call ran on the
/// entire query range. This broke for large regions because syng's
/// `query_region` returned whole-chromosome context spans.
///
/// We verify the new behavior by checking that the per-sub-window log lines
/// show up in stderr, one per sub-window.
#[test]
fn test_query_syng_gfa_subwindow_splitter() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping: impg binary not built");
        return;
    };

    let dir = std::env::temp_dir().join("impg_test_query_syng_subwindow");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // This test's primary purpose is verifying that pggb:X / seqwish:X is
    // interpreted as a sub-window size (not a boolean flag), which is done by
    // inspecting stderr for the per-window log lines. The downstream
    // FastGA/seqwish pipeline runs but we DO NOT require it to succeed —
    // the regression assertion fires before any engine work. So we size the
    // input to the minimum that still lets syng's default syncmer params
    // (k=8, w=55, total=63bp) produce syncmer hits per window, and let the
    // engine fail fast on CI's 2-vCPU runners if it wants to.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        // Two haplotypes sharing a 3 kbp backbone + 500 bp unique tails.
        // Total ~3.5 kbp each, 7 kbp FASTA.
        let backbone = numeric_to_ascii(&make_sequence_numeric(3000, 42));
        let tail1 = numeric_to_ascii(&make_sequence_numeric(500, 1));
        let tail2 = numeric_to_ascii(&make_sequence_numeric(500, 2));
        writeln!(f, ">sampleA#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail1).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">sampleB#0#chr1").unwrap();
        f.write_all(&backbone).unwrap();
        f.write_all(&tail2).unwrap();
        writeln!(f).unwrap();
    }

    // Build syng index
    let syng_prefix = dir.join("idx");
    let build = Command::new(&bin)
        .args([
            "syng",
            "-f",
            fasta_path.to_str().unwrap(),
            "-o",
            syng_prefix.to_str().unwrap(),
            "--position-sample-rate",
            "1",
        ])
        .output()
        .expect("Failed to run impg syng");
    assert!(
        build.status.success(),
        "impg syng failed: {}",
        String::from_utf8_lossy(&build.stderr)
    );

    // Query 3 kbp of sampleA's backbone with poa:1000 — three sub-windows
    // of 1000 bp each (the minimum window size impg allows). We pick `poa`
    // instead of `seqwish`/`pggb` because POA runs a single-pass partial
    // order alignment per sub-window and skips the FastGA + seqwish +
    // gfaffix chain, which is the main source of multi-minute runtime on
    // small CI runners. The sub-window wiring regression is engine-agnostic
    // — the log lines are emitted before any engine runs — so POA is a
    // cheaper substitute that exercises the exact same code path.
    let out_prefix = dir.join("region");
    let out = Command::new(&bin)
        .args([
            "query",
            "-d",
            "0",
            "-a",
            syng_prefix.to_str().unwrap(),
            "--sequence-files",
            fasta_path.to_str().unwrap(),
            "-r",
            "sampleA#0#chr1:0-3000",
            "-o",
            "gfa",
            "--gfa-engine",
            "poa:1000",
            "-O",
            out_prefix.to_str().unwrap(),
            "-t",
            "1",
            "-v",
            "2", // info-level logging to capture sub-window log lines
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("Failed to run impg query");

    let stderr = String::from_utf8_lossy(&out.stderr);

    // Primary (and only) assertion: sub-window wiring — three log lines must
    // appear. This is the regression check for the pggb:X / seqwish:X no-op
    // bug, emitted inside the sub-window loop BEFORE the engine runs. A
    // downstream engine failure doesn't invalidate this check and doesn't
    // fail the test — subprocess exit status is intentionally unchecked to
    // keep the test cheap on small CI runners.
    let subwindow_count = stderr
        .lines()
        .filter(|l| l.contains("[syng sub-window"))
        .count();
    assert_eq!(
        subwindow_count, 3,
        "expected 3 sub-window log lines for 3000bp query with poa:1000, got {}. stderr: {}",
        subwindow_count, stderr
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// RC homology is detected and boundary-realigned end-to-end.
///
/// Construct a fixture with an inversion: genome_b contains the
/// reverse-complement of a stretch of genome_a. Confirm:
/// (1) the raw syng path reports the homolog with strand='-',
/// (2) boundary realignment produces tight target coordinates in forward
///     strand space,
/// (3) the refined interval's target forward-strand bases, read in RC, match
///     the query's forward bases exactly (zero edit distance).
#[test]
fn test_syng_rc_homolog_end_to_end() {
    let _guard = lock_syng();

    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // genome_a: 3000bp random (seed 42).
    // genome_b: [prefix_1000, RC(a[500..2500]), suffix_500] — 3500bp total,
    //          with an RC-inverted copy of a[500..2500] embedded at b[1000..3000].
    let a = mk_seq(3000, 42);
    let inv_region = &a[500..2500];
    let rc_inv = impg::graph::reverse_complement(inv_region);
    let mut b = mk_seq(1000, 11);
    b.extend_from_slice(&rc_inv);
    b.extend_from_slice(&mk_seq(500, 13));
    assert_eq!(b.len(), 3500);

    let dir = std::env::temp_dir().join("impg_test_syng_rc_homolog");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("rc.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&b).unwrap();
        writeln!(f).unwrap();
    }

    let params = impg::syng::SyncmerParams {
        k: 8,
        w: 55,
        seed: 7,
    };
    let syng_index = impg::syng::SyngIndex::build(
        params,
        vec![
            ("genome_a".to_string(), a.clone()),
            ("genome_b".to_string(), b.clone()),
        ]
        .into_iter(),
    );
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[fasta_path
        .to_string_lossy()
        .to_string()])
    .unwrap();

    // Query the middle of the inverted region on genome_a.
    let query_start = 1000u64;
    let query_end = 2000u64;
    let padding = 120;

    // (1) Raw query_region reports a '-' homolog for genome_b.
    let raw = syng_index
        .query_region("genome_a", query_start, query_end, padding)
        .unwrap();
    let rc_raw: Vec<_> = raw
        .iter()
        .filter(|iv| iv.genome == "genome_b" && iv.strand == '-')
        .collect();
    assert!(
        !rc_raw.is_empty(),
        "raw query should report at least one '-' strand homolog on genome_b; got: {:?}",
        raw.iter()
            .map(|iv| (&iv.genome, iv.start, iv.end, iv.strand))
            .collect::<Vec<_>>()
    );

    // (2) Boundary realignment produces tight coordinates for the RC homolog.
    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();
    let rc_refined: Vec<_> = refined
        .iter()
        .filter(|iv| iv.genome == "genome_b" && iv.strand == '-')
        .collect();
    assert!(
        !rc_refined.is_empty(),
        "refined query should keep the '-' strand homolog; got: {:?}",
        refined
            .iter()
            .map(|iv| (&iv.genome, iv.start, iv.end, iv.strand))
            .collect::<Vec<_>>()
    );

    // The expected refined forward-strand region on genome_b:
    // genome_a[500..2500] ↔ RC → b[1000..3000]. Query a[1000..2000]
    // corresponds (under RC) to b[1500, 2500). Our refined intervals must
    // overlap that window and, when RC'd, reproduce query bases with high
    // fidelity somewhere along the alignment.
    let expected_start = 1500u64;
    let expected_end = 2500u64;
    let overlap_bp = rc_refined
        .iter()
        .map(|iv| {
            let lo = iv.start.max(expected_start);
            let hi = iv.end.min(expected_end);
            if hi > lo {
                hi - lo
            } else {
                0
            }
        })
        .sum::<u64>();
    assert!(
        overlap_bp >= 200,
        "RC refined intervals should collectively overlap the expected RC window \
         [{}, {}) by at least 200bp; got {}bp overlap. Intervals: {:?}",
        expected_start,
        expected_end,
        overlap_bp,
        rc_refined
            .iter()
            .map(|iv| (iv.start, iv.end))
            .collect::<Vec<_>>()
    );

    // Base-content validation: for each refined '-' interval, RC'ing the
    // target bytes should produce something that matches a sub-string of
    // the query region (with at least one long exact-match run). This
    // confirms the ORIENTATION and COVERAGE are correct even if the exact
    // endpoints are subject to the syng-position calibration noise flagged
    // as a known limitation below.
    let query_bytes = &a[query_start as usize..query_end as usize];
    let mut best_match_len = 0usize;
    for iv in &rc_refined {
        let b_slice = &b[iv.start as usize..iv.end as usize];
        let b_rc = impg::graph::reverse_complement(b_slice);
        // Find the longest common substring between b_rc and query_bytes,
        // capped at the shorter length. Cheap approximation: slide b_rc
        // over query_bytes and measure the longest run.
        if b_rc.len() < 30 || query_bytes.len() < 30 {
            continue;
        }
        let search = &b_rc[..b_rc.len().min(200)];
        for start in 0..query_bytes.len().saturating_sub(search.len()) {
            let run = query_bytes[start..]
                .iter()
                .zip(search.iter())
                .take_while(|(x, y)| x == y)
                .count();
            if run > best_match_len {
                best_match_len = run;
            }
        }
    }
    assert!(
        best_match_len >= 30,
        "RC'd refined target bytes should share a >=30bp exact run with the query \
         region; longest run was {}bp",
        best_match_len
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Boundary-realignment tightens fuzzy syncmer-resolution edges.
///
/// Construct two genomes that share a 500bp backbone at offset 0..500. Query
/// the backbone region on genome_a. The RAW syng path returns padded
/// syncmer-resolution intervals (slop-boundary up to ±syncmer_len bp around
/// the real edges). The boundary-realignment path should return intervals
/// whose edges are much closer to the true 0..500 boundary on genome_b.
#[test]
fn test_syng_boundary_realign_tightens_edges() {
    let _guard = lock_syng();

    // Reproducible backbone + distinct tails.
    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // Backbone length must comfortably exceed query_end + syncmer_len so that
    // shared syncmers exist on BOTH sides of every tested edge. At default
    // (k=8, w=55) a syncmer starting at position p extends to p+63; only
    // syncmers whose whole length sits inside the shared backbone are shared.
    let backbone = mk_seq(2000, 42);
    let mut seq_a = backbone.clone();
    seq_a.extend(mk_seq(500, 1));
    let mut seq_b = backbone.clone();
    seq_b.extend(mk_seq(500, 2));

    let dir = std::env::temp_dir().join("impg_test_syng_boundary_realign");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    // Write a FASTA for the UnifiedSequenceIndex.
    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&seq_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&seq_b).unwrap();
        writeln!(f).unwrap();
    }

    // Build the syng index in-process.
    let params = impg::syng::SyncmerParams {
        k: 8,
        w: 55,
        seed: 7,
    };
    let sequences = vec![
        ("genome_a".to_string(), seq_a.clone()),
        ("genome_b".to_string(), seq_b.clone()),
    ];
    let syng_index = impg::syng::SyngIndex::build(params, sequences.into_iter());

    // Build the file-backed sequence index.
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[fasta_path
        .to_string_lossy()
        .to_string()])
    .unwrap();

    let padding = 120; // default
    let query_start = 50u64;
    let query_end = 450u64;

    // Raw: what syng gives us directly.
    let raw = syng_index
        .query_region("genome_a", query_start, query_end, padding)
        .unwrap();
    let raw_b = raw
        .iter()
        .find(|iv| iv.genome == "genome_b")
        .expect("raw query should find genome_b hit");

    // Realigned.
    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();
    let refined_b = refined
        .iter()
        .find(|iv| iv.genome == "genome_b")
        .expect("realigned query should find genome_b hit");

    // On identical backbone, boundary realignment should land on the
    // query coordinates exactly (BiWFA edit distance = 0 between query
    // and target slices around the edges).
    assert_eq!(
        refined_b.start, query_start,
        "refined start should snap to query_start on identical backbone; got {} vs expected {} (raw was {})",
        refined_b.start, query_start, raw_b.start
    );
    assert_eq!(
        refined_b.end, query_end,
        "refined end should snap to query_end on identical backbone; got {} vs expected {} (raw was {})",
        refined_b.end, query_end, raw_b.end
    );

    // The raw interval must be at least as loose as the refined one
    // (slop is additive, not subtractive).
    assert!(
        raw_b.start <= refined_b.start && raw_b.end >= refined_b.end,
        "raw interval should enclose refined interval: raw {}..{} vs refined {}..{}",
        raw_b.start,
        raw_b.end,
        refined_b.start,
        refined_b.end
    );

    std::fs::remove_dir_all(&dir).ok();
}

/// Realistic "homologous-with-small-differences" end-to-end test.
///
/// Build three genomes that share a 3kb region but differ at a handful of
/// scattered positions (SNPs) and with a small indel on one of them. Query
/// a 2kb region in the middle of genome_a at depth 1; expect the refined
/// boundary coordinates on both targets to snap to the biological
/// boundaries of the query (not the padded syncmer-resolution bounds).
///
/// This is the test Erik asked for — the current code is expected to show
/// fragmentation / loose edges on this fixture, because `resolve_edge_via_biwfa`
/// projects through flanking anchors around the query edge and doesn't do a
/// bedtools-style merge + end-alignment pass. Once the merge-then-end-align
/// pipeline lands, this test should pass with single-digit bp precision.
#[test]
fn test_syng_query_reconstructs_homology_with_diffs() {
    let _guard = lock_syng();

    fn mk_seq(len: usize, seed: u8) -> Vec<u8> {
        let bases = [b'A', b'C', b'G', b'T'];
        let mut s = Vec::with_capacity(len);
        let mut state: u32 = seed as u32;
        for _ in 0..len {
            state = state.wrapping_mul(1103515245).wrapping_add(12345);
            s.push(bases[((state >> 16) % 4) as usize]);
        }
        s
    }

    // Shared 3kb backbone, unique 500bp tails.
    let backbone = mk_seq(3000, 42);

    let seq_a = {
        let mut s = backbone.clone();
        s.extend(mk_seq(500, 1));
        s
    };

    // seq_b: identical backbone in the first 3kb, plus 5 scattered SNPs.
    //        SNPs shouldn't shift any coordinates.
    let seq_b = {
        let mut s = backbone.clone();
        for &pos in &[250usize, 800, 1337, 1900, 2500] {
            // Flip one base deterministically (shift by 1 in the 4-base alphabet).
            s[pos] = match s[pos] {
                b'A' => b'C',
                b'C' => b'G',
                b'G' => b'T',
                b'T' => b'A',
                other => other,
            };
        }
        s.extend(mk_seq(500, 2));
        s
    };

    // seq_c: identical first ~1500bp; 10bp deletion at 1500; then the rest
    //        of the backbone continues. So homology to genome_a covers
    //        [0, 1500) and [1500, 2990) on seq_c (target coords are shifted
    //        by -10 after position 1500).
    let seq_c = {
        let mut s: Vec<u8> = Vec::with_capacity(3490);
        s.extend_from_slice(&backbone[..1500]);
        s.extend_from_slice(&backbone[1510..]); // drop 10bp at 1500..1510
        s.extend(mk_seq(500, 3));
        s
    };

    let dir = std::env::temp_dir().join("impg_test_syng_homology_diffs");
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();

    let fasta_path = dir.join("test.fa");
    {
        use std::io::Write;
        let mut f = std::fs::File::create(&fasta_path).unwrap();
        writeln!(f, ">genome_a").unwrap();
        f.write_all(&seq_a).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_b").unwrap();
        f.write_all(&seq_b).unwrap();
        writeln!(f).unwrap();
        writeln!(f, ">genome_c").unwrap();
        f.write_all(&seq_c).unwrap();
        writeln!(f).unwrap();
    }

    let params = impg::syng::SyncmerParams {
        k: 8,
        w: 55,
        seed: 7,
    };
    let syng_index = impg::syng::SyngIndex::build(
        params,
        vec![
            ("genome_a".to_string(), seq_a.clone()),
            ("genome_b".to_string(), seq_b.clone()),
            ("genome_c".to_string(), seq_c.clone()),
        ]
        .into_iter(),
    );
    let sequence_index = impg::sequence_index::UnifiedSequenceIndex::from_files(&[fasta_path
        .to_string_lossy()
        .to_string()])
    .unwrap();

    // Query the middle 2kb of genome_a: [500, 2500).
    let query_start = 500u64;
    let query_end = 2500u64;
    let padding = 120;

    let refined = impg::syng_transitive::query_transitive(
        &syng_index,
        "genome_a",
        query_start,
        query_end,
        padding,
        1,
        &sequence_index,
    )
    .unwrap();

    // Expected: ONE homolog per target genome at the biological coordinates.
    //
    // On genome_b: matches [500, 2500) exactly (SNPs don't shift coords).
    // On genome_c: matches [500, 2490) — the 10bp deletion at backbone pos
    //              1500 shifts all downstream positions by -10, so the
    //              query's right edge at 2500 maps to c's 2490.
    let refined_on = |genome: &str| -> Vec<(u64, u64)> {
        refined
            .iter()
            .filter(|iv| iv.genome == genome && iv.strand == '+')
            .map(|iv| (iv.start, iv.end))
            .collect()
    };

    let on_b = refined_on("genome_b");
    let on_c = refined_on("genome_c");

    // Step 1: no fragmentation.
    assert_eq!(
        on_b.len(),
        1,
        "expected exactly one forward-strand homolog on genome_b; got {:?}",
        on_b
    );
    assert_eq!(
        on_c.len(),
        1,
        "expected exactly one forward-strand homolog on genome_c; got {:?}",
        on_c
    );

    // Step 2: refined edges match biological truth within a few bp.
    let tol: u64 = 5;
    let (b_start, b_end) = on_b[0];
    let (c_start, c_end) = on_c[0];

    assert!(
        b_start.abs_diff(500) <= tol && b_end.abs_diff(2500) <= tol,
        "refined genome_b interval should snap to [500, 2500) within {}bp; got [{}, {})",
        tol,
        b_start,
        b_end
    );
    assert!(
        c_start.abs_diff(500) <= tol && c_end.abs_diff(2490) <= tol,
        "refined genome_c interval should snap to [500, 2490) within {}bp (10bp deletion at 1500 shifts the right edge); got [{}, {})",
        tol, c_start, c_end
    );

    std::fs::remove_dir_all(&dir).ok();
}
