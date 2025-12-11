//! Integration test for the full impg pipeline: index -> partition -> graph -> lace
//! Uses yeast chrV test data (7 strains)
//! Requires wfmash and samtools to be installed

use flate2::read::GzDecoder;
use std::fs::{self, File};
use std::io::{Read, Write};
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

fn get_impg_binary() -> PathBuf {
    // CARGO_BIN_EXE_impg is set by cargo test for the binary crate
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return PathBuf::from(path);
    }

    // Get manifest dir and look for binary relative to it
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));

    let candidates = [
        manifest_dir.join("target/release/impg"),
        manifest_dir.join("target/debug/impg"),
    ];

    for path in &candidates {
        if path.exists() {
            return path.clone();
        }
    }

    // Fall back to PATH
    PathBuf::from("impg")
}

fn decompress_gz(gz_path: &PathBuf, out_path: &PathBuf) -> std::io::Result<()> {
    let gz_file = File::open(gz_path)?;
    let mut decoder = GzDecoder::new(gz_file);
    let mut contents = Vec::new();
    decoder.read_to_end(&mut contents)?;
    let mut out_file = File::create(out_path)?;
    out_file.write_all(&contents)?;
    Ok(())
}

fn run_impg(work_dir: &PathBuf, args: &[&str]) -> std::io::Result<std::process::Output> {
    let impg = get_impg_binary();
    Command::new(&impg).current_dir(work_dir).args(args).output()
}

#[test]
fn test_full_pipeline() -> std::io::Result<()> {
    // Create temp directory
    let temp_dir = TempDir::new()?;
    let work_dir = temp_dir.path().to_path_buf();

    // Test data paths (7-strain yeast chrV)
    let test_data_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("tests")
        .join("test_data");

    let fasta_gz = test_data_dir.join("yeast.chrV.fa.gz");

    // Check test data exists
    assert!(fasta_gz.exists(), "Test FASTA not found: {:?}", fasta_gz);

    // Decompress FASTA
    let fasta_path = work_dir.join("seqs.fa");
    decompress_gz(&fasta_gz, &fasta_path)?;

    // Create FASTA index using samtools
    let samtools_result = Command::new("samtools")
        .args(["faidx", fasta_path.to_str().unwrap()])
        .output();

    if samtools_result.is_err() || !samtools_result.as_ref().unwrap().status.success() {
        eprintln!("Warning: samtools faidx failed, test may fail");
    }

    // Generate PAF alignments using wfmash
    let paf_path = work_dir.join("alignments.paf");
    let wfmash_output = Command::new("wfmash")
        .args([
            fasta_path.to_str().unwrap(),
            "-p",
            "90",
            "-n",
            "2",
            "-t",
            "2",
        ])
        .stdout(std::process::Stdio::piped())
        .stderr(std::process::Stdio::null())
        .output()?;

    if !wfmash_output.status.success() {
        panic!(
            "wfmash failed - is it installed? This test requires wfmash to generate alignments."
        );
    }

    // Write PAF output
    let mut paf_file = File::create(&paf_path)?;
    paf_file.write_all(&wfmash_output.stdout)?;

    // Step 1: Create IMPG index
    let output = run_impg(
        &work_dir,
        &["index", "-a", "alignments.paf", "-i", "test.impg"],
    )?;
    assert!(
        output.status.success(),
        "Index failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let index_path = work_dir.join("test.impg");
    assert!(index_path.exists(), "Index file not created");

    // Step 2: Partition sequences
    let partitions_dir = work_dir.join("partitions");
    fs::create_dir_all(&partitions_dir)?;

    let output = run_impg(
        &work_dir,
        &[
            "partition",
            "-a",
            "alignments.paf",
            "-i",
            "test.impg",
            "-w",
            "200000",
            "--sequence-files",
            "seqs.fa",
            "-o",
            "fasta",
            "--output-folder",
            "partitions",
            "--separate-files",
            "-t",
            "2",
        ],
    )?;
    assert!(
        output.status.success(),
        "Partition failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    // Count partitions
    let partition_count = fs::read_dir(&partitions_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "fasta"))
        .count();
    assert!(partition_count >= 1, "No partitions created");

    // Step 3: Build graphs for each partition
    let gfas_dir = work_dir.join("gfas");
    fs::create_dir_all(&gfas_dir)?;

    for entry in fs::read_dir(&partitions_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(false, |ext| ext == "fasta") {
            let stem = path.file_stem().unwrap().to_str().unwrap();
            let gfa_name = format!("{}.gfa", stem);

            let output = run_impg(
                &work_dir,
                &[
                    "graph",
                    "--fasta-files",
                    &format!("partitions/{}.fasta", stem),
                    "-g",
                    &format!("gfas/{}", gfa_name),
                    "-t",
                    "2",
                ],
            )?;
            assert!(
                output.status.success(),
                "Graph failed for {}: {}",
                stem,
                String::from_utf8_lossy(&output.stderr)
            );
        }
    }

    // Count GFAs
    let gfa_count = fs::read_dir(&gfas_dir)?
        .filter_map(|e| e.ok())
        .filter(|e| e.path().extension().map_or(false, |ext| ext == "gfa"))
        .count();
    assert_eq!(
        gfa_count, partition_count,
        "Expected {} GFAs, got {}",
        partition_count, gfa_count
    );

    // Step 4: Create GFA list file
    let gfa_list_path = work_dir.join("gfa_list.txt");
    let mut gfa_list = File::create(&gfa_list_path)?;
    for entry in fs::read_dir(&gfas_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.extension().map_or(false, |ext| ext == "gfa") {
            writeln!(gfa_list, "{}", path.display())?;
        }
    }
    drop(gfa_list);

    // Step 5: Lace GFAs into pangenome
    let output = run_impg(
        &work_dir,
        &[
            "lace",
            "--file-list",
            "gfa_list.txt",
            "--sequence-files",
            "seqs.fa",
            "-o",
            "pangenome.gfa",
            "-t",
            "2",
        ],
    )?;
    assert!(
        output.status.success(),
        "Lace failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );

    let pangenome_path = work_dir.join("pangenome.gfa");
    assert!(pangenome_path.exists(), "Pangenome GFA not created");

    // Validate: check that the GFA has content
    let pangenome_size = fs::metadata(&pangenome_path)?.len();
    assert!(
        pangenome_size > 1000,
        "Pangenome too small: {} bytes",
        pangenome_size
    );

    // Parse GFA header to count nodes
    let gfa_content = fs::read_to_string(&pangenome_path)?;
    let node_count = gfa_content.lines().filter(|l| l.starts_with("S\t")).count();
    let path_count = gfa_content.lines().filter(|l| l.starts_with("P\t")).count();

    assert!(node_count >= 10, "Too few nodes: {}", node_count);
    assert_eq!(path_count, 7, "Expected 7 paths, got {}", path_count);

    eprintln!(
        "Pipeline test passed: {} nodes, {} paths, {} bytes",
        node_count, path_count, pangenome_size
    );

    Ok(())
}
