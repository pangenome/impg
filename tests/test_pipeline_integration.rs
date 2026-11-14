//! Integration test for the full impg pipeline: index -> partition -> graph -> lace
//! Uses yeast chrV test data (7 strains)
//! Uses the wfmash that the impg build puts beside the impg binary

use std::fs::{self, File};
use std::io::Write;
use std::path::PathBuf;
use std::process::Command;
use tempfile::TempDir;

fn get_impg_binary() -> PathBuf {
    // The binary cargo built for this test run
    PathBuf::from(env!("CARGO_BIN_EXE_impg"))
}

fn decompress_gz(gz_path: &PathBuf, out_path: &PathBuf) -> std::io::Result<()> {
    // gzip -dc reads standard gzip and BGZIP
    let output = Command::new("gzip").arg("-dc").arg(gz_path).output()?;
    if !output.status.success() {
        return Err(std::io::Error::new(
            std::io::ErrorKind::Other,
            "gzip -dc failed",
        ));
    }
    let mut out_file = File::create(out_path)?;
    out_file.write_all(&output.stdout)?;
    Ok(())
}

/// Write the first `bases` of the first `records` FASTA records.
fn write_subset(
    full: &PathBuf,
    subset: &PathBuf,
    records: usize,
    bases: usize,
) -> std::io::Result<Vec<(String, String)>> {
    let text = fs::read_to_string(full)?;
    let mut kept: Vec<(String, String)> = Vec::new();
    for line in text.lines() {
        if let Some(name) = line.strip_prefix('>') {
            if kept.len() == records {
                break;
            }
            kept.push((name.to_string(), String::new()));
        } else if let Some((_, seq)) = kept.last_mut() {
            let room = bases - seq.len();
            seq.push_str(&line[..line.len().min(room)].to_ascii_uppercase());
        }
    }
    let mut out = File::create(subset)?;
    for (name, seq) in &kept {
        writeln!(out, ">{name}\n{seq}")?;
    }
    Ok(kept)
}

fn run_impg(work_dir: &PathBuf, args: &[&str]) -> std::io::Result<std::process::Output> {
    let impg = get_impg_binary();
    Command::new(&impg)
        .current_dir(work_dir)
        .args(args)
        .output()
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

    // Decompress FASTA, keep 4 strains x 60 kb
    let full_path = work_dir.join("full.fa");
    decompress_gz(&fasta_gz, &full_path)?;
    let fasta_path = work_dir.join("seqs.fa");
    let inputs = write_subset(&full_path, &fasta_path, 4, 60_000)?;
    assert_eq!(inputs.len(), 4, "subset records: {}", inputs.len());

    // Generate PAF alignments using wfmash
    let paf_path = work_dir.join("alignments.paf");
    let wfmash = get_impg_binary().with_file_name("wfmash");
    assert!(wfmash.exists(), "bundled wfmash not found: {:?}", wfmash);
    let wfmash_output = Command::new(&wfmash)
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
        panic!("bundled wfmash failed with {}", wfmash_output.status);
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
            "-d",
            "10000",
            "-a",
            "alignments.paf",
            "-i",
            "test.impg",
            "-w",
            "20000",
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
                    "--sequence-files",
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
    assert_eq!(path_count, 4, "Expected 4 paths, got {}", path_count);

    // Every laced path spells its input sequence
    let segments: std::collections::HashMap<&str, &str> = gfa_content
        .lines()
        .filter(|l| l.starts_with("S\t"))
        .map(|l| {
            let mut f = l.split('\t');
            (f.nth(1).unwrap(), f.next().unwrap())
        })
        .collect();
    let complement = |b: u8| match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        other => other,
    };
    for (name, expected) in &inputs {
        let line = gfa_content
            .lines()
            .find(|l| l.starts_with("P\t") && l.split('\t').nth(1) == Some(name.as_str()))
            .unwrap_or_else(|| panic!("no path named {name}"));
        let mut spelled = String::with_capacity(expected.len());
        for step in line.split('\t').nth(2).unwrap().split(',') {
            let (id, orient) = step.split_at(step.len() - 1);
            let seq = segments[id].to_ascii_uppercase();
            if orient == "+" {
                spelled.push_str(&seq);
            } else {
                spelled.extend(seq.bytes().rev().map(|b| complement(b) as char));
            }
        }
        assert!(spelled == *expected, "path {name} does not spell its input");
    }

    eprintln!(
        "Pipeline test passed: {} nodes, {} paths, {} bytes",
        node_count, path_count, pangenome_size
    );

    Ok(())
}
