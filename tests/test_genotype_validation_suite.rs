//! Focused truth-known validation for syng pack/genotype/infer evidence.
//!
//! These cases are intentionally small and deterministic. They complement the
//! broader syng integration smoke tests by asserting exact evidence vectors and
//! simple genotype order expectations on synthetic panels.

use rustc_hash::FxHashMap;
use std::collections::{BTreeMap, BTreeSet};
use std::io::{self, Read, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

// Syng's C layer has process-global state. Keep all tests in this binary
// serialized even when cargo runs this integration test with multiple threads.
static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
    std::sync::LazyLock::new(|| std::sync::Mutex::new(()));

fn lock_syng() -> std::sync::MutexGuard<'static, ()> {
    SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner())
}

fn impg_binary() -> Option<PathBuf> {
    if let Ok(path) = std::env::var("CARGO_BIN_EXE_impg") {
        return Some(PathBuf::from(path));
    }

    let manifest_dir = Path::new(env!("CARGO_MANIFEST_DIR"));
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

fn make_sequence_numeric(len: usize, seed: u8) -> Vec<u8> {
    let mut seq = Vec::with_capacity(len);
    let mut state = seed as u32;
    for _ in 0..len {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(((state >> 16) % 4) as u8);
    }
    seq
}

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

fn mutate_base_ascii(base: u8) -> u8 {
    match base {
        b'A' => b'C',
        b'C' => b'G',
        b'G' => b'T',
        b'T' => b'A',
        other => other,
    }
}

fn mutate_ascii_every(seq: &[u8], offset: usize, stride: usize) -> Vec<u8> {
    let mut out = seq.to_vec();
    for i in (offset..out.len()).step_by(stride) {
        out[i] = mutate_base_ascii(out[i]);
    }
    out
}

fn temp_case_dir(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("impg_{name}_{}", std::process::id()));
    std::fs::remove_dir_all(&dir).ok();
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn build_syng_index(prefix: &Path, sequences: &[(&str, Vec<u8>)]) {
    let mut index = impg::syng::SyngIndex::build(
        impg::syng::SyncmerParams::default(),
        sequences
            .iter()
            .map(|(name, seq)| ((*name).to_string(), seq.clone())),
    );
    index.save(prefix.to_str().unwrap()).unwrap();
}

fn write_fastq_record<W: Write>(writer: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writeln!(writer, "@{name}")?;
    writer.write_all(seq)?;
    writeln!(writer)?;
    writeln!(writer, "+")?;
    writeln!(writer, "{}", "I".repeat(seq.len()))
}

fn write_fasta_record<W: Write>(writer: &mut W, name: &str, seq: &[u8]) -> io::Result<()> {
    writeln!(writer, ">{name}")?;
    writer.write_all(seq)?;
    writeln!(writer)
}

fn write_fastq_reads(path: &Path, reads: &[(&str, &[u8])]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (name, seq) in reads {
        write_fastq_record(&mut f, name, seq).unwrap();
    }
}

fn write_fasta_reads(path: &Path, reads: &[(&str, &[u8])]) {
    let mut f = std::fs::File::create(path).unwrap();
    for (name, seq) in reads {
        write_fasta_record(&mut f, name, seq).unwrap();
    }
}

fn write_tiled_fastq<W: Write>(
    writer: &mut W,
    prefix: &str,
    seq: &[u8],
    read_len: usize,
    step: usize,
) -> io::Result<usize> {
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
        write_fastq_record(
            writer,
            &format!("{prefix}_{read_idx}"),
            &seq[*start..start + read_len],
        )?;
    }
    Ok(starts.len())
}

fn parse_query_records(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).unwrap();
    let mut records = Vec::new();
    let first = text
        .lines()
        .find(|line| !line.trim().is_empty())
        .expect("query should not be empty")
        .as_bytes()[0];
    if first == b'>' {
        let mut name = String::new();
        let mut seq = Vec::new();
        for line in text.lines() {
            if let Some(header) = line.strip_prefix('>') {
                if !name.is_empty() {
                    records.push((std::mem::take(&mut name), std::mem::take(&mut seq)));
                }
                name = header.split_whitespace().next().unwrap_or("").to_string();
            } else {
                seq.extend_from_slice(line.trim().as_bytes());
            }
        }
        if !name.is_empty() {
            records.push((name, seq));
        }
    } else {
        let mut lines = text.lines();
        while let Some(header) = lines.next() {
            if header.trim().is_empty() {
                continue;
            }
            assert!(header.starts_with('@'), "bad FASTQ header: {header}");
            let seq = lines.next().unwrap().as_bytes().to_vec();
            let plus = lines.next().unwrap();
            assert!(plus.starts_with('+'), "bad FASTQ plus line: {plus}");
            let _qual = lines.next().unwrap();
            records.push((
                header[1..]
                    .split_whitespace()
                    .next()
                    .unwrap_or("")
                    .to_string(),
                seq,
            ));
        }
    }
    records
}

fn expected_pack_counts_from_matcher(
    prefix: &Path,
    query_path: &Path,
    min_anchors: usize,
) -> BTreeMap<u32, u64> {
    let matcher = impg::syng::SyngMatcher::load(
        prefix.to_str().unwrap(),
        impg::syng::SyncmerParams::default(),
    )
    .unwrap();
    let mut counts = BTreeMap::new();
    for (_, seq) in parse_query_records(query_path) {
        let syncmers = matcher.matched_syncmers_in_sequence(&seq);
        let distinct_nodes: BTreeSet<u32> =
            syncmers.iter().map(|syncmer| syncmer.node_id).collect();
        if distinct_nodes.len() < min_anchors {
            continue;
        }
        for node_id in distinct_nodes {
            *counts.entry(node_id).or_insert(0) += 1;
        }
    }
    counts
}

fn occurrence_counts_from_matcher(prefix: &Path, query_path: &Path) -> BTreeMap<u32, u64> {
    let matcher = impg::syng::SyngMatcher::load(
        prefix.to_str().unwrap(),
        impg::syng::SyncmerParams::default(),
    )
    .unwrap();
    let mut counts = BTreeMap::new();
    for (_, seq) in parse_query_records(query_path) {
        for syncmer in matcher.matched_syncmers_in_sequence(&seq) {
            *counts.entry(syncmer.node_id).or_insert(0) += 1;
        }
    }
    counts
}

fn read_pack_tsv_counts_from_text(text: &str) -> BTreeMap<u32, u64> {
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

fn run_pack_tsv(
    bin: &Path,
    prefix: &Path,
    query_path: &Path,
    min_anchors: usize,
) -> BTreeMap<u32, u64> {
    let output = Command::new(bin)
        .args([
            "map",
            "-a",
            prefix.to_str().unwrap(),
            "-q",
            query_path.to_str().unwrap(),
            "-o",
            "pack-tsv",
            "--min-anchors",
            &min_anchors.to_string(),
            "-t",
            "1",
        ])
        .output()
        .expect("failed to run impg map -o pack-tsv");
    assert!(
        output.status.success(),
        "impg map -o pack-tsv failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    read_pack_tsv_counts_from_text(&String::from_utf8_lossy(&output.stdout))
}

fn top_result_fields(stdout: &[u8]) -> Vec<String> {
    let stdout = String::from_utf8_lossy(stdout);
    let top = stdout
        .lines()
        .find(|line| !line.starts_with('#') && !line.trim().is_empty())
        .unwrap_or_else(|| panic!("expected a result row in:\n{stdout}"));
    top.split('\t').map(str::to_string).collect()
}

fn run_genotype(
    bin: &Path,
    prefix: &Path,
    pack_path: &Path,
    target_range: &str,
    extra: &[&str],
) -> Vec<String> {
    let mut args = vec![
        "genotype",
        "cos",
        "-a",
        prefix.to_str().unwrap(),
        "-p",
        pack_path.to_str().unwrap(),
        "-r",
        target_range,
        "--top-n",
        "5",
        "--candidate-top-k",
        "20",
        "--min-anchors",
        "2",
        "--min-span-fraction",
        "0.6",
    ];
    args.extend_from_slice(extra);
    let output = Command::new(bin)
        .args(args)
        .output()
        .expect("failed to run impg genotype cos");
    assert!(
        output.status.success(),
        "impg genotype cos failed: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    top_result_fields(&output.stdout)
}

fn write_pack_tsv(path: &Path, counts: &BTreeMap<u32, u64>) {
    let mut f = std::fs::File::create(path).unwrap();
    writeln!(f, "#node_id\tcount").unwrap();
    for (node_id, count) in counts {
        writeln!(f, "{node_id}\t{count}").unwrap();
    }
}

#[test]
fn test_pack_tsv_matches_independent_syncmer_vector_for_fasta_and_fastq() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI validation: impg binary not found");
        return;
    };

    let dir = temp_case_dir("genotype_validation_pack_vector");
    let left = numeric_to_ascii(&make_sequence_numeric(720, 1));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(520, 2));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(520, 3));
    let right = numeric_to_ascii(&make_sequence_numeric(720, 4));

    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);

    let prefix = dir.join("panel.syng");
    build_syng_index(
        &prefix,
        &[
            ("sampleA#0#chr1", hap_a.clone()),
            ("sampleB#0#chr1", hap_b.clone()),
        ],
    );

    let short = &hap_a[0..40];
    let reads: Vec<(&str, &[u8])> = vec![
        ("a_left", &hap_a[100..540]),
        ("a_allele", &hap_a[760..1250]),
        ("b_allele", &hap_b[760..1250]),
        ("too_short_for_syncmers", short),
    ];
    let fastq = dir.join("reads.fq");
    let fasta = dir.join("reads.fa");
    write_fastq_reads(&fastq, &reads);
    write_fasta_reads(&fasta, &reads);

    let expected_fastq = expected_pack_counts_from_matcher(&prefix, &fastq, 2);
    let expected_fasta = expected_pack_counts_from_matcher(&prefix, &fasta, 2);
    assert_eq!(
        expected_fastq, expected_fasta,
        "FASTA and FASTQ parsing should feed the same matched-syncmer oracle"
    );
    assert!(
        !expected_fastq.is_empty(),
        "synthetic reads should produce nonempty expected pack counts"
    );

    let cli_fastq = run_pack_tsv(&bin, &prefix, &fastq, 2);
    let cli_fasta = run_pack_tsv(&bin, &prefix, &fasta, 2);
    assert_eq!(
        cli_fastq, expected_fastq,
        "pack-tsv must equal independently accumulated distinct matched syncmers"
    );
    assert_eq!(
        cli_fasta, expected_fasta,
        "FASTA query input should produce the same pack vector as FASTQ"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_genotype_truth_known_homozygote_heterozygote_and_decoy() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI validation: impg binary not found");
        return;
    };

    let dir = temp_case_dir("genotype_validation_truth_known");
    let left = numeric_to_ascii(&make_sequence_numeric(850, 11));
    let allele_a = numeric_to_ascii(&make_sequence_numeric(760, 12));
    let allele_b = numeric_to_ascii(&make_sequence_numeric(760, 13));
    let allele_decoy = numeric_to_ascii(&make_sequence_numeric(760, 14));
    let right = numeric_to_ascii(&make_sequence_numeric(850, 15));

    let mut hap_a = Vec::new();
    hap_a.extend_from_slice(&left);
    hap_a.extend_from_slice(&allele_a);
    hap_a.extend_from_slice(&right);
    let mut hap_b = Vec::new();
    hap_b.extend_from_slice(&left);
    hap_b.extend_from_slice(&allele_b);
    hap_b.extend_from_slice(&right);
    let mut hap_decoy = Vec::new();
    hap_decoy.extend_from_slice(&left);
    hap_decoy.extend_from_slice(&allele_decoy);
    hap_decoy.extend_from_slice(&right);

    let prefix = dir.join("panel.syng");
    build_syng_index(
        &prefix,
        &[
            ("sampleA#0#chr1", hap_a.clone()),
            ("sampleB#0#chr1", hap_b.clone()),
            ("sampleDecoy#0#chr1", hap_decoy),
        ],
    );

    let target = format!("sampleA#0#chr1:0-{}", hap_a.len());

    let hom_reads = dir.join("hom_a.fq");
    write_fastq_reads(&hom_reads, &[("a_copy_0", &hap_a), ("a_copy_1", &hap_a)]);
    let hom_pack = dir.join("hom_a.pack.tsv");
    write_pack_tsv(
        &hom_pack,
        &expected_pack_counts_from_matcher(&prefix, &hom_reads, 2),
    );
    let hom_top = run_genotype(&bin, &prefix, &hom_pack, &target, &[]);
    assert_eq!(hom_top[0], "1");
    assert_eq!(hom_top[2], "2");
    assert!(
        hom_top[8] == "sampleA#0#chr1,sampleA#0#chr1",
        "easy homozygote should rank A/A first, got {:?}",
        hom_top
    );
    assert!(
        !hom_top[8].contains("sampleDecoy"),
        "decoy haplotype must not appear in the homozygote top call"
    );

    let het_reads = dir.join("het_ab.fq");
    write_fastq_reads(&het_reads, &[("a_copy", &hap_a), ("b_copy", &hap_b)]);
    let het_pack = dir.join("het_ab.pack.tsv");
    write_pack_tsv(
        &het_pack,
        &expected_pack_counts_from_matcher(&prefix, &het_reads, 2),
    );
    let het_top = run_genotype(&bin, &prefix, &het_pack, &target, &[]);
    assert_eq!(het_top[0], "1");
    assert_eq!(het_top[2], "2");
    assert!(
        het_top[8].contains("sampleA#0#chr1") && het_top[8].contains("sampleB#0#chr1"),
        "heterozygote should rank A/B first, got {:?}",
        het_top
    );
    assert!(
        !het_top[8].contains("sampleDecoy"),
        "decoy haplotype must not appear in the heterozygote top call"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_repeated_node_pack_dedup_and_cnv_counting_counterfactual() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI validation: impg binary not found");
        return;
    };

    let dir = temp_case_dir("genotype_validation_cnv_semantics");
    let left = numeric_to_ascii(&make_sequence_numeric(540, 21));
    let copy = numeric_to_ascii(&make_sequence_numeric(720, 22));
    let right = numeric_to_ascii(&make_sequence_numeric(540, 23));

    let mut hap_single = Vec::new();
    hap_single.extend_from_slice(&left);
    hap_single.extend_from_slice(&copy);
    hap_single.extend_from_slice(&right);
    let mut hap_double = Vec::new();
    hap_double.extend_from_slice(&left);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&copy);
    hap_double.extend_from_slice(&right);

    let prefix = dir.join("panel.syng");
    build_syng_index(
        &prefix,
        &[
            ("sampleSingle#0#chr1", hap_single),
            ("sampleDouble#0#chr1", hap_double.clone()),
        ],
    );

    let reads = dir.join("one_double_read.fq");
    write_fastq_reads(&reads, &[("double_full", &hap_double)]);

    let expected_dedup = expected_pack_counts_from_matcher(&prefix, &reads, 2);
    let occurrence = occurrence_counts_from_matcher(&prefix, &reads);
    let cli_pack = run_pack_tsv(&bin, &prefix, &reads, 2);
    assert_eq!(
        cli_pack, expected_dedup,
        "current pack semantics are distinct syng nodes per retained read"
    );
    let repeated_nodes: Vec<u32> = occurrence
        .iter()
        .filter_map(|(&node_id, &occ_count)| {
            let dedup_count = expected_dedup.get(&node_id).copied().unwrap_or(0);
            (occ_count > dedup_count).then_some(node_id)
        })
        .collect();
    assert!(
        !repeated_nodes.is_empty(),
        "CNV read should contain repeated syncmer-node occurrences"
    );

    let repeat = repeated_nodes[0];
    let flank_nodes: Vec<u32> = expected_dedup
        .keys()
        .copied()
        .filter(|node| *node != repeat)
        .take(2)
        .collect();
    assert_eq!(
        flank_nodes.len(),
        2,
        "CNV fixture should have at least two non-repeated flank nodes"
    );
    let single_copy = vec![(flank_nodes[0], 1), (repeat, 1), (flank_nodes[1], 1)];
    let double_copy = vec![(flank_nodes[0], 1), (repeat, 2), (flank_nodes[1], 1)];
    let candidate_features: Vec<&[(u32, u64)]> = vec![&single_copy, &double_copy];

    let mut dedup_sample = FxHashMap::default();
    dedup_sample.insert(flank_nodes[0], 1);
    dedup_sample.insert(repeat, 1);
    dedup_sample.insert(flank_nodes[1], 1);
    let features = impg::genotyping::feature_universe(candidate_features.iter().copied());
    let dedup_norm = impg::genotyping::sample_norm_sq_for_features(&dedup_sample, &features);
    let dedup_rank = impg::genotyping::run_cosine_combination_search(
        &candidate_features,
        &dedup_sample,
        dedup_norm,
        1,
        10,
    )
    .unwrap();
    assert_eq!(
        dedup_rank[0].combination,
        vec![0],
        "per-read dedup sample vector ranks the single-copy model first"
    );

    let mut occurrence_sample = FxHashMap::default();
    occurrence_sample.insert(flank_nodes[0], 1);
    occurrence_sample.insert(repeat, 2);
    occurrence_sample.insert(flank_nodes[1], 1);
    let occurrence_norm =
        impg::genotyping::sample_norm_sq_for_features(&occurrence_sample, &features);
    let occurrence_rank = impg::genotyping::run_cosine_combination_search(
        &candidate_features,
        &occurrence_sample,
        occurrence_norm,
        1,
        10,
    )
    .unwrap();
    assert_eq!(
        occurrence_rank[0].combination,
        vec![1],
        "per-occurrence sample vector ranks the duplicated model first"
    );

    std::fs::remove_dir_all(&dir).ok();
}

#[test]
fn test_projection_bundle_gaf_read_walk_evidence_rewards_recombinant_stitching() {
    let _guard = lock_syng();
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI validation: impg binary not found");
        return;
    };

    let dir = temp_case_dir("genotype_validation_projection_gaf");
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

    let prefix = dir.join("panel.syng");
    build_syng_index(
        &prefix,
        &[
            ("sampleA#0#chr1", hap_a),
            ("sampleB#0#chr1", hap_b),
            ("sampleC#0#chr1", hap_c.clone()),
            ("sampleD#0#chr1", hap_d.clone()),
        ],
    );

    let reads = dir.join("recombinant_reads.fq");
    {
        let mut f = std::fs::File::create(&reads).unwrap();
        write_tiled_fastq(&mut f, "hapC", &hap_c, 1200, 175).unwrap();
        write_tiled_fastq(&mut f, "hapD", &hap_d, 1200, 175).unwrap();
    }

    let proj = dir.join("sample.proj");
    let map = Command::new(&bin)
        .args([
            "map",
            "-a",
            prefix.to_str().unwrap(),
            "-q",
            reads.to_str().unwrap(),
            "-o",
            "proj",
            "-O",
            proj.to_str().unwrap(),
            "--pack-compression-level",
            "3",
            "--pack-block-size",
            "64",
            "--min-anchors",
            "2",
            "-t",
            "1",
        ])
        .output()
        .expect("failed to run impg map -o proj");
    assert!(
        map.status.success(),
        "impg map -o proj failed: {}",
        String::from_utf8_lossy(&map.stderr)
    );
    let manifest = std::fs::read_to_string(proj.join("manifest.json")).unwrap();
    assert!(
        manifest.contains("\"gaf\"") && manifest.contains("\"reads.gaf.zst\""),
        "projection manifest should record the read-walk GAF path:\n{manifest}"
    );
    let mut gaf_decoder =
        zstd::stream::read::Decoder::new(std::fs::File::open(proj.join("reads.gaf.zst")).unwrap())
            .unwrap();
    let mut gaf = String::new();
    gaf_decoder.read_to_string(&mut gaf).unwrap();
    assert!(
        gaf.lines().any(|line| line
            .split('\t')
            .nth(5)
            .is_some_and(|path| { path.contains('>') || path.contains('<') })),
        "projection bundle should contain oriented syncmer read walks:\n{gaf}"
    );

    let partitions = dir.join("partitions.bed");
    {
        let mut f = std::fs::File::create(&partitions).unwrap();
        writeln!(f, "sampleA#0#chr1\t0\t950\tleft").unwrap();
        writeln!(f, "sampleA#0#chr1\t950\t1900\tright").unwrap();
    }

    let mosaic_path = dir.join("mosaic.tsv");
    let infer = Command::new(&bin)
        .args([
            "infer",
            "-a",
            prefix.to_str().unwrap(),
            "--proj",
            proj.to_str().unwrap(),
            "--partitions",
            partitions.to_str().unwrap(),
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
        "impg infer with projection GAF evidence failed: {}",
        String::from_utf8_lossy(&infer.stderr)
    );
    let mosaic = std::fs::read_to_string(&mosaic_path).unwrap();
    assert!(
        mosaic.contains("read_link_reward"),
        "mosaic output should expose read-link scoring columns:\n{mosaic}"
    );
    let rows: Vec<Vec<&str>> = mosaic
        .lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| line.split('\t').collect())
        .collect();
    assert_eq!(
        rows.len(),
        4,
        "expected two phases across two partitions:\n{mosaic}"
    );
    assert!(
        rows.iter()
            .filter(|row| row[1] == "right")
            .all(|row| row[16].parse::<f64>().unwrap() > 0.0),
        "right-partition transitions should be rewarded by spanning GAF read walks:\n{mosaic}"
    );

    let mut phase_paths: BTreeMap<&str, Vec<&str>> = BTreeMap::new();
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
        "read-walk links should choose the sampled recombinant haplotypes:\n{mosaic}"
    );

    std::fs::remove_dir_all(&dir).ok();
}
