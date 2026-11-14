//! Regression tests for the crush engine on real C4 data.
//!
//! The checks for the open failures listed in docs/crush-audit.md needed the
//! full C4A graph, which is not in the repository, and are not here.
//!
//! ## Real data used:
//!
//! - `tests/test_data/crush/c4_slice_1500_3000.gfa` — 2942 segments / 64k bp
//!   extracted from the C4A blunt pangenome GFA (chr6:31891045-32123783, 465
//!   haplotypes, seqwish-k=311). Real chromosomal data, NOT synthetic.

use impg::resolution::{
    path_sequences, resolve_gfa_bubbles, ReplacementMinMatchLenPolicy, ResolutionConfig,
    ResolutionMethod,
};
use std::collections::{HashMap, HashSet};
use std::path::Path;
use std::process::Command;
use std::sync::{Mutex, MutexGuard, Once};

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

static C4_FRAGMENT_TEST_THREADS: Once = Once::new();
static C4_FRAGMENT_TEST_LOCK: Mutex<()> = Mutex::new(());

fn c4_fragment_test_guard() -> MutexGuard<'static, ()> {
    C4_FRAGMENT_TEST_THREADS.call_once(|| {
        let _ = rayon::ThreadPoolBuilder::new()
            .num_threads(32)
            .build_global();
    });
    C4_FRAGMENT_TEST_LOCK
        .lock()
        .unwrap_or_else(|poisoned| poisoned.into_inner())
}

/// Locate the compiled `impg` binary.
fn impg_binary() -> Option<std::path::PathBuf> {
    // The binary cargo built for this test run
    Some(std::path::PathBuf::from(env!("CARGO_BIN_EXE_impg")))
}

/// Count S-lines and sum their sequence lengths in a GFA text.
fn gfa_segment_stats(gfa: &str) -> (usize, usize) {
    let mut count = 0usize;
    let mut bp = 0usize;
    for line in gfa.lines() {
        if line.starts_with("S\t") {
            let mut parts = line.splitn(4, '\t');
            parts.next(); // "S"
            parts.next(); // id
            if let Some(seq) = parts.next() {
                bp += seq.len();
                count += 1;
            }
        }
    }
    (count, bp)
}

fn seq_map(gfa: &str) -> HashMap<String, String> {
    path_sequences(gfa)
        .expect("path_sequences failed")
        .into_iter()
        .collect()
}

fn fasta_records(path: &Path) -> Vec<(String, Vec<u8>)> {
    let text = std::fs::read_to_string(path).expect("failed to read FASTA");
    let mut records = Vec::new();
    let mut name: Option<String> = None;
    let mut seq = Vec::new();
    for line in text.lines() {
        if let Some(header) = line.strip_prefix('>') {
            if let Some(prev) = name.replace(header.to_string()) {
                records.push((prev, std::mem::take(&mut seq)));
            }
        } else {
            seq.extend_from_slice(line.trim().as_bytes());
        }
    }
    if let Some(name) = name {
        records.push((name, seq));
    }
    records
}

fn gfa_shared_path_segment_count(gfa: &str) -> usize {
    let mut support: HashMap<String, HashSet<String>> = HashMap::new();
    for line in gfa.lines().filter(|line| line.starts_with("P\t")) {
        let mut fields = line.split('\t');
        fields.next();
        let Some(path_name) = fields.next() else {
            continue;
        };
        let Some(steps) = fields.next() else {
            continue;
        };
        for step in steps.split(',').filter(|step| !step.is_empty()) {
            let node = step
                .strip_suffix('+')
                .or_else(|| step.strip_suffix('-'))
                .unwrap_or(step);
            support
                .entry(node.to_string())
                .or_default()
                .insert(path_name.to_string());
        }
    }
    support.values().filter(|paths| paths.len() > 1).count()
}

fn gfa_max_path_segment_depth(gfa: &str) -> usize {
    let mut support: HashMap<String, HashSet<String>> = HashMap::new();
    for line in gfa.lines().filter(|line| line.starts_with("P\t")) {
        let mut fields = line.split('\t');
        fields.next();
        let Some(path_name) = fields.next() else {
            continue;
        };
        let Some(steps) = fields.next() else {
            continue;
        };
        for step in steps.split(',').filter(|step| !step.is_empty()) {
            let node = step
                .strip_suffix('+')
                .or_else(|| step.strip_suffix('-'))
                .unwrap_or(step);
            support
                .entry(node.to_string())
                .or_default()
                .insert(path_name.to_string());
        }
    }
    support.values().map(HashSet::len).max().unwrap_or_default()
}

fn gfa_path_count(gfa: &str) -> usize {
    gfa.lines().filter(|line| line.starts_with("P\t")).count()
}

fn gfa_duplicate_segment_seqs(gfa: &str) -> usize {
    let mut seqs: HashSet<String> = HashSet::new();
    let mut duplicates = 0usize;
    for line in gfa.lines().filter(|line| line.starts_with("S\t")) {
        let mut fields = line.splitn(4, '\t');
        fields.next();
        fields.next();
        if let Some(seq) = fields.next() {
            if !seqs.insert(seq.to_string()) {
                duplicates += 1;
            }
        }
    }
    duplicates
}

#[derive(Default, Debug)]
struct PafFixtureStats {
    records: usize,
    aligned_bp: usize,
    name_mismatches: usize,
    length_mismatches: usize,
}

fn paf_fixture_stats(paf: &str, seqs: &[(String, Vec<u8>)]) -> PafFixtureStats {
    let lengths: HashMap<&str, usize> = seqs
        .iter()
        .map(|(name, seq)| (name.as_str(), seq.len()))
        .collect();
    let mut stats = PafFixtureStats::default();
    for line in paf.lines().filter(|line| !line.trim().is_empty()) {
        let fields: Vec<_> = line.split('\t').collect();
        if fields.len() < 12 {
            continue;
        }
        stats.records += 1;
        stats.aligned_bp += fields[10].parse::<usize>().unwrap_or_default();
        let q_len = lengths.get(fields[0]).copied();
        let t_len = lengths.get(fields[5]).copied();
        if q_len.is_none() || t_len.is_none() {
            stats.name_mismatches += 1;
            continue;
        }
        if fields[1].parse::<usize>().ok() != q_len || fields[6].parse::<usize>().ok() != t_len {
            stats.length_mismatches += 1;
        }
    }
    stats
}

fn assert_gfa_paths_spell_fasta(gfa: &str, seqs: &[(String, Vec<u8>)], context: &str) {
    let expected: HashMap<String, String> = seqs
        .iter()
        .map(|(name, seq)| {
            (
                name.clone(),
                String::from_utf8(seq.clone()).expect("fixture sequence is not UTF-8"),
            )
        })
        .collect();
    let observed: HashMap<String, String> = path_sequences(gfa)
        .expect("failed to spell GFA paths")
        .into_iter()
        .collect();
    assert_eq!(
        observed.len(),
        expected.len(),
        "{context}: path count changed"
    );
    for (name, expected_seq) in &expected {
        let observed_seq = observed
            .get(name)
            .unwrap_or_else(|| panic!("{context}: path {name} disappeared"));
        assert_eq!(
            observed_seq, expected_seq,
            "{context}: path {name} changed spelling"
        );
    }
}

fn common_prefix_len(seqs: &[Vec<u8>]) -> usize {
    let min_len = seqs.iter().map(Vec::len).min().unwrap_or_default();
    (0..min_len)
        .take_while(|&idx| seqs.iter().all(|seq| seq[idx] == seqs[0][idx]))
        .count()
}

fn common_suffix_len(seqs: &[Vec<u8>], max_len: usize) -> usize {
    let min_len = seqs.iter().map(Vec::len).min().unwrap_or_default();
    (0..min_len.min(max_len))
        .take_while(|&idx| {
            seqs.iter()
                .all(|seq| seq[seq.len() - 1 - idx] == seqs[0][seqs[0].len() - 1 - idx])
        })
        .count()
}

fn one_bubble_gfa_from_fasta_records(seqs: &[(String, Vec<u8>)]) -> String {
    let raw_seqs: Vec<Vec<u8>> = seqs.iter().map(|(_, seq)| seq.clone()).collect();
    let left_len = common_prefix_len(&raw_seqs).min(20);
    let right_len = common_suffix_len(
        &raw_seqs,
        raw_seqs.iter().map(Vec::len).min().unwrap() - left_len,
    )
    .min(20);
    assert!(
        left_len > 0 && right_len > 0,
        "C4 fragment fixture must have shared flanks to build a bubble"
    );
    let exit_id = seqs.len() + 2;
    let mut gfa = String::from("H\tVN:Z:1.0\n");
    gfa.push_str(&format!(
        "S\t1\t{}\n",
        String::from_utf8(raw_seqs[0][..left_len].to_vec()).unwrap()
    ));
    for (idx, (_, seq)) in seqs.iter().enumerate() {
        assert!(
            left_len + right_len < seq.len(),
            "C4 fragment fixture traversal must contain non-empty bubble interior"
        );
        gfa.push_str(&format!(
            "S\t{}\t{}\n",
            idx + 2,
            String::from_utf8(seq[left_len..seq.len() - right_len].to_vec()).unwrap()
        ));
    }
    gfa.push_str(&format!(
        "S\t{exit_id}\t{}\n",
        String::from_utf8(raw_seqs[0][raw_seqs[0].len() - right_len..].to_vec()).unwrap()
    ));
    for idx in 0..seqs.len() {
        let seg_id = idx + 2;
        gfa.push_str(&format!("L\t1\t+\t{seg_id}\t+\t0M\n"));
        gfa.push_str(&format!("L\t{seg_id}\t+\t{exit_id}\t+\t0M\n"));
    }
    for (idx, (name, _)) in seqs.iter().enumerate() {
        let seg_id = idx + 2;
        gfa.push_str(&format!("P\t{name}\t1+,{seg_id}+,{exit_id}+\t*\n"));
    }
    gfa
}

struct C4FragmentFixture {
    name: &'static str,
    fasta_rel: &'static str,
    paf_rel: &'static str,
    protects: &'static str,
    min_paf_records: usize,
    min_shared_segments: usize,
    min_shared_depth: usize,
    max_segments: usize,
    max_segment_bp: usize,
    max_duplicate_segment_seqs: usize,
}

const C4_FRAGMENT_FIXTURES: &[C4FragmentFixture] = &[
    C4FragmentFixture {
        name: "easy_shared_flank",
        fasta_rel: "tests/test_data/crush/c4_fragments/easy_shared_flank.fa",
        paf_rel: "tests/test_data/crush/c4_fragments/easy_shared_flank.paf",
        protects: "simple top-flubble block with obvious shared flanks and two allelic interiors",
        min_paf_records: 12,
        min_shared_segments: 1,
        min_shared_depth: 3,
        max_segments: 4,
        max_segment_bp: 700,
        max_duplicate_segment_seqs: 0,
    },
    C4FragmentFixture {
        name: "bounded_multi_bubble",
        fasta_rel: "tests/test_data/crush/c4_fragments/bounded_multi_bubble.fa",
        paf_rel: "tests/test_data/crush/c4_fragments/bounded_multi_bubble.paf",
        protects: "locally complex but bounded C4 top-flubble block with three compact allelic classes",
        min_paf_records: 12,
        min_shared_segments: 1,
        min_shared_depth: 4,
        max_segments: 12,
        max_segment_bp: 2_000,
        max_duplicate_segment_seqs: 6,
    },
    C4FragmentFixture {
        name: "unfolded_minrun",
        fasta_rel: "tests/test_data/crush/top_flubble_seqwish_minrun.fa",
        paf_rel: "tests/test_data/crush/top_flubble_seqwish_minrun.paf",
        protects: "previous broken output class where nonempty C4 PAF became path-specific unfolded segments",
        min_paf_records: 16,
        min_shared_segments: 1,
        min_shared_depth: 2,
        max_segments: 24,
        max_segment_bp: 1800,
        max_duplicate_segment_seqs: 6,
    },
    C4FragmentFixture {
        name: "short_floor",
        fasta_rel: "tests/test_data/crush/c4_fragments/short_floor.fa",
        paf_rel: "tests/test_data/crush/c4_fragments/short_floor.paf",
        protects: "short C4 block at the 159-300 bp length floor where filtering and seqwish-k must not erase alignments",
        min_paf_records: 20,
        min_shared_segments: 2,
        min_shared_depth: 4,
        max_segments: 10,
        max_segment_bp: 900,
        max_duplicate_segment_seqs: 1,
    },
    C4FragmentFixture {
        name: "duplicated_repeat",
        fasta_rel: "tests/test_data/crush/c4_fragments/duplicated_repeat.fa",
        paf_rel: "tests/test_data/crush/c4_fragments/duplicated_repeat.paf",
        protects: "repeated/duplicated C4 sequence where alignments must survive filtering and be consumed by seqwish",
        min_paf_records: 40,
        min_shared_segments: 2,
        min_shared_depth: 4,
        max_segments: 10,
        max_segment_bp: 900,
        max_duplicate_segment_seqs: 2,
    },
];

#[test]
fn c4_fragment_seqwish_regressions_induce_shared_graphs() {
    let _guard = c4_fragment_test_guard();
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for fixture in C4_FRAGMENT_FIXTURES {
        let fasta_path = manifest_dir.join(fixture.fasta_rel);
        let paf_path = manifest_dir.join(fixture.paf_rel);
        let seqs = fasta_records(&fasta_path);
        let paf = std::fs::read_to_string(&paf_path)
            .unwrap_or_else(|err| panic!("failed to read {}: {err}", paf_path.display()));
        let input_bp: usize = seqs.iter().map(|(_, seq)| seq.len()).sum();
        let paf_stats = paf_fixture_stats(&paf, &seqs);
        assert!(
            paf_stats.records >= fixture.min_paf_records,
            "{}: expected at least {} real PAF records, got {}; protects: {}",
            fixture.name,
            fixture.min_paf_records,
            paf_stats.records,
            fixture.protects
        );
        assert_eq!(
            paf_stats.name_mismatches, 0,
            "{}: PAF query/target names must match FASTA IDs",
            fixture.name
        );
        assert_eq!(
            paf_stats.length_mismatches, 0,
            "{}: PAF query/target lengths must match FASTA lengths",
            fixture.name
        );
        assert!(
            paf_stats.aligned_bp > 0,
            "{}: homologous C4 fragment should have nonempty aligned bp",
            fixture.name
        );

        let config = impg::commands::graph::GraphBuildConfig {
            num_threads: 1,
            min_match_len: 1,
            adaptive_min_match_len: false,
            min_map_length: 0,
            no_filter: true,
            num_mappings: "many:many".to_string(),
            scaffold_filter: "many:many".to_string(),
            scaffold_mass: 0,
            show_progress: false,
            ..impg::commands::graph::GraphBuildConfig::default()
        };
        let gfa = impg::syng_graph::build_gfa_from_paf_and_sequences(&seqs, &paf, &config)
            .unwrap_or_else(|err| panic!("{}: seqwish induction failed: {err}", fixture.name));
        assert_gfa_paths_spell_fasta(&gfa, &seqs, fixture.name);

        let (segments, segment_bp) = gfa_segment_stats(&gfa);
        let shared_segments = gfa_shared_path_segment_count(&gfa);
        let shared_depth = gfa_max_path_segment_depth(&gfa);
        let duplicate_segment_seqs = gfa_duplicate_segment_seqs(&gfa);
        eprintln!(
            "{} seqwish: paths={}, segments={}, bp={}, shared={}, depth={}, duplicate_seqs={}",
            fixture.name,
            gfa_path_count(&gfa),
            segments,
            segment_bp,
            shared_segments,
            shared_depth,
            duplicate_segment_seqs
        );
        assert_eq!(
            gfa_path_count(&gfa),
            seqs.len(),
            "{}: seqwish path count changed",
            fixture.name
        );
        assert!(
            shared_segments >= fixture.min_shared_segments,
            "{}: induced graph has only {} shared segment(s), expected >= {}; protects: {}",
            fixture.name,
            shared_segments,
            fixture.min_shared_segments,
            fixture.protects
        );
        assert!(
            shared_depth >= fixture.min_shared_depth,
            "{}: max shared-node depth {} below expected {}; GFA was:\n{}",
            fixture.name,
            shared_depth,
            fixture.min_shared_depth,
            gfa
        );
        assert!(
            segments <= fixture.max_segments,
            "{}: induced graph has {} segment(s), expected <= {}",
            fixture.name,
            segments,
            fixture.max_segments
        );
        assert!(
            segment_bp <= fixture.max_segment_bp,
            "{}: induced graph has {} segment bp, expected <= {}",
            fixture.name,
            segment_bp,
            fixture.max_segment_bp
        );
        assert!(
            segment_bp < input_bp,
            "{}: induced graph looks unfolded/concatenated (segment bp {} >= input bp {})",
            fixture.name,
            segment_bp,
            input_bp
        );
        assert!(
            duplicate_segment_seqs <= fixture.max_duplicate_segment_seqs,
            "{}: duplicate segment sequence count {} exceeds bound {}; duplicate/unfolded signature not eliminated",
            fixture.name,
            duplicate_segment_seqs,
            fixture.max_duplicate_segment_seqs
        );
    }
}

#[test]
fn c4_fragment_lacing_uses_pairwise_induced_graphs() {
    let _guard = c4_fragment_test_guard();
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    for fixture in C4_FRAGMENT_FIXTURES {
        let fasta_path = manifest_dir.join(fixture.fasta_rel);
        let seqs = fasta_records(&fasta_path);
        let gfa = one_bubble_gfa_from_fasta_records(&seqs);
        let before: HashMap<_, _> = path_sequences(&gfa).unwrap().into_iter().collect();
        let (before_segments, before_bp) = gfa_segment_stats(&gfa);
        let config = ResolutionConfig {
            max_iterations: 1,
            method: ResolutionMethod::TopFlubbleSweepga,
            replacement_seqwish_min_match_len: 311,
            replacement_min_match_len_policy: ReplacementMinMatchLenPolicy::Adaptive,
            replacement_min_map_length: 311,
            replacement_num_mappings: "many:many".to_string(),
            replacement_scaffold_filter: "many:many".to_string(),
            replacement_scaffold_mass: 0,
            sweepga_aligner: "fastga".to_string(),
            sweepga_min_aln_length: 1,
            max_pair_alignments: 0,
            max_replacement_paf_bytes: 0,
            polish_iterations: 0,
            polish_max_traversal_len: 0,
            polish_max_median_traversal_len: 0,
            polish_max_total_sequence: 0,
            polish_max_traversals: 0,
            ..ResolutionConfig::default()
        };
        let resolved = resolve_gfa_bubbles(&gfa, &config)
            .unwrap_or_else(|err| panic!("{}: top-flubble crush failed: {err}", fixture.name));
        assert_eq!(
            resolved.stats.bailed, 0,
            "{}: no quality-guard fallback or bailout should hide bad behavior; stats={:?}",
            fixture.name, resolved.stats
        );
        assert_eq!(
            resolved.stats.resolved, 1,
            "{}: expected exactly one top-flubble replacement to be applied; stats={:?}",
            fixture.name, resolved.stats
        );
        let after: HashMap<_, _> = path_sequences(&resolved.gfa).unwrap().into_iter().collect();
        assert_eq!(
            after, before,
            "{}: replacement/lacing changed path spellings",
            fixture.name
        );
        let (after_segments, after_bp) = gfa_segment_stats(&resolved.gfa);
        let after_shared = gfa_shared_path_segment_count(&resolved.gfa);
        assert!(
            after_shared > 0,
            "{}: laced graph has no shared path-supported nodes",
            fixture.name
        );
        assert!(
            after_bp < before_bp,
            "{}: laced graph still looks like concatenated traversal sequence (bp {} -> {})",
            fixture.name,
            before_bp,
            after_bp
        );
        assert!(
            after_segments < before_segments,
            "{}: replacement should reduce segment count on this bounded fixture ({} -> {})",
            fixture.name,
            before_segments,
            after_segments
        );
    }
}

#[test]
fn c4_motif_zero_paf_fastga_seqwish_variant_is_rejected_without_replacement() {
    let gfa = "\
H\tVN:Z:1.0
S\tleft\tAA
S\tshared\tC
S\tright\tTT
S\tf00\tG
S\tf01\tG
S\tf02\tG
S\tf03\tG
S\tf04\tG
S\tf05\tG
S\tf06\tG
S\tf07\tG
S\tf08\tG
S\tf09\tG
S\tf10\tG
S\tf11\tG
S\tf12\tG
S\tf13\tG
S\tf14\tG
S\tf15\tG
S\tprivate\tG
L\tleft\t+\tshared\t+\t0M
L\tshared\t+\tright\t+\t0M
L\tleft\t+\tprivate\t+\t0M
L\tprivate\t+\tright\t+\t0M
P\tHG00001#1#chr6:100-105\tleft+,shared+,right+\t*
P\tHG00002#1#chr6:200-205\tleft+,shared+,right+\t*
P\tHG00003#1#chr6:300-305\tleft+,shared+,right+\t*
P\tHG00004#1#chr6:400-405\tleft+,shared+,right+\t*
P\tHG00005#1#chr6:500-505\tleft+,shared+,right+\t*
P\tHG00006#1#chr6:600-605\tleft+,private+,right+\t*
";
    let before_paths = seq_map(gfa);
    let resolved = resolve_gfa_bubbles(
        gfa,
        &ResolutionConfig {
            max_iterations: 1,
            method: ResolutionMethod::MotifLocal,
            auto_spoa_max_traversal_len: 0,
            auto_poasta_max_traversal_len: 0,
            motif_max_sparse_paths: 1,
            motif_min_flank_paths: 4,
            motif_min_order_jump: 8,
            motif_max_window_bp: 100,
            sweepga_aligner: "fastga".to_string(),
            sweepga_min_aln_length: 1_000,
            polish_iterations: 0,
            ..ResolutionConfig::default()
        },
    )
    .expect("motif-local zero-PAF candidate should be rejected, not fatal");

    assert_eq!(resolved.stats.resolved, 0, "{:?}", resolved.stats);
    assert!(
        resolved.stats.candidates_seen > 0,
        "synthetic motif should be admitted before FastGA rejects it: {:?}",
        resolved.stats
    );
    assert!(
        resolved.stats.bailed > 0,
        "zero-PAF FastGA/seqwish variant should count as a rejected build: {:?}",
        resolved.stats
    );
    assert_eq!(
        resolved.gfa, gfa,
        "zero-PAF motif-local variant must not replace the graph"
    );
    assert_eq!(seq_map(&resolved.gfa), before_paths);
}

// ---------------------------------------------------------------------------
// Companion test: small C4 slice with auto method should not degrade quality
//
// This test uses the committed test data slice and PASSES on HEAD because
// auto-SPOA on small bubbles (<500 bp) actually improves quality. It confirms
// the crush API is callable and path sequences are preserved.
// ---------------------------------------------------------------------------

#[test]
fn c4_slice_auto_crush_preserves_path_sequences() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read C4 slice GFA");
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();

    let resolved = resolve_gfa_bubbles(&gfa, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on C4 slice");

    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();

    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs.get(name).expect("path disappeared after crush");
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after crush",
            name
        );
    }

    let (before_segs, before_bp) = gfa_segment_stats(&gfa);
    let (after_segs, after_bp) = gfa_segment_stats(&resolved.gfa);
    eprintln!(
        "slice crush: {} paths preserved, {} resolved, {} bailed",
        before_seqs.len(),
        resolved.stats.resolved,
        resolved.stats.bailed
    );
    eprintln!(
        "slice: segs {} → {}, bp {} → {}",
        before_segs, after_segs, before_bp, after_bp
    );
}

#[test]
#[ignore = "needs IMPG_ABPOA_BIN"]
fn c4_slice_abpoa_crush_preserves_path_sequences_when_configured() {
    let abpoa_bin = std::env::var("IMPG_ABPOA_BIN")
        .expect("IMPG_ABPOA_BIN is not set; abPOA C4 slice smoke needs it");
    if let Err(err) = Command::new(&abpoa_bin).arg("-v").output() {
        panic!("configured IMPG_ABPOA_BIN is not runnable: {abpoa_bin}: {err}");
    }

    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read C4 slice GFA");
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();

    let resolved = resolve_gfa_bubbles(
        &gfa,
        &ResolutionConfig {
            method: ResolutionMethod::Abpoa,
            abpoa_bin,
            max_iterations: 1,
            max_traversal_len: 10_000,
            max_median_traversal_len: 10_000,
            max_total_sequence: 1_000_000,
            ..ResolutionConfig::default()
        },
    )
    .expect("resolve_gfa_bubbles failed on C4 slice with abPOA");

    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();
    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs
            .get(name)
            .unwrap_or_else(|| panic!("path {name} disappeared after abPOA crush"));
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after abPOA crush",
            name
        );
    }

    assert!(
        resolved.stats.resolved > 0,
        "abPOA C4 slice smoke reached candidates but accepted no replacements: {:?}",
        resolved.stats
    );
    eprintln!(
        "abPOA slice crush: {} paths preserved, {} resolved, {} bailed",
        before_seqs.len(),
        resolved.stats.resolved,
        resolved.stats.bailed
    );
}

// ---------------------------------------------------------------------------
// Wider-context bubble resolution — see docs/crush-wider-context-bubbles.md.
// Verifies that enabling `replacement_flank_bp` on the canonical C4 slice still
// preserves every path sequence exactly. The integration code clips the
// flanking aligned portion before substitution; if that clipping is wrong, the
// resolved path sequence will not equal the input path sequence.
// ---------------------------------------------------------------------------

#[test]
fn c4_slice_auto_crush_with_flank_preserves_path_sequences() {
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }

    let gfa = std::fs::read_to_string(&gfa_path).expect("failed to read C4 slice GFA");
    let before_seqs: std::collections::HashMap<_, _> =
        path_sequences(&gfa).unwrap().into_iter().collect();

    let mut config = ResolutionConfig::default();
    config.replacement_flank_bp = 500;
    let resolved = resolve_gfa_bubbles(&gfa, &config)
        .expect("resolve_gfa_bubbles with flank=500 failed on C4 slice");

    let after_seqs: std::collections::HashMap<_, _> =
        path_sequences(&resolved.gfa).unwrap().into_iter().collect();

    for (name, before_seq) in &before_seqs {
        let after_seq = after_seqs.get(name).expect("path disappeared after crush");
        assert_eq!(
            before_seq, after_seq,
            "path {} changed sequence after crush with flank=500",
            name
        );
    }
    eprintln!(
        "slice crush (flank=500): {} paths preserved, {} resolved, {} bailed",
        before_seqs.len(),
        resolved.stats.resolved,
        resolved.stats.bailed
    );
}

#[test]
fn c4_slice_outward_admission_only_cli_reports_without_building() {
    let Some(bin) = impg_binary() else {
        eprintln!("Skipping CLI admission-only test: impg binary not found");
        return;
    };
    let manifest_dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR"));
    let gfa_path = manifest_dir.join("tests/test_data/crush/c4_slice_1500_3000.gfa");
    if !gfa_path.exists() {
        eprintln!("SKIP: test data slice not found at {:?}", gfa_path);
        return;
    }
    let out_path = std::env::temp_dir().join(format!(
        "impg_c4_outward_admission_only_{}.gfa",
        std::process::id()
    ));
    std::fs::remove_file(&out_path).ok();

    let run = Command::new(&bin)
        .args([
            "crush",
            "--gfa",
            gfa_path.to_str().unwrap(),
            "--output",
            out_path.to_str().unwrap(),
            "--method",
            "iterative-multi-level",
            "--window-mode",
            "outward",
            "--window-target-bp",
            "200k",
            "--max-window-sites",
            "8",
            "--candidate-limit",
            "16",
            "--admission-only",
            "--max-iterations",
            "1",
            "-v",
            "2",
        ])
        .env("RUST_LOG", "info")
        .output()
        .expect("failed to run impg crush admission-only CLI");
    let stderr = String::from_utf8_lossy(&run.stderr);
    assert!(
        run.status.success(),
        "admission-only crush failed:\n{}",
        stderr
    );
    assert!(
        stderr.contains("admission-only dry-run admitted"),
        "missing admission-only summary:\n{}",
        stderr
    );
    assert!(
        stderr.contains("source=outward-residual-window"),
        "missing outward residual admission detail:\n{}",
        stderr
    );
    assert!(
        !stderr.contains("replacement build concurrency"),
        "dry-run must not schedule replacement builds:\n{}",
        stderr
    );
    assert!(
        !stderr.contains("building candidate"),
        "dry-run must not enter replacement builders:\n{}",
        stderr
    );

    let input = std::fs::read_to_string(&gfa_path).expect("failed to reread C4 slice fixture");
    let output = std::fs::read_to_string(&out_path).expect("failed to read dry-run output");
    assert_eq!(
        input, output,
        "admission-only dry-run should not rewrite GFA"
    );
    std::fs::remove_file(&out_path).ok();
}

#[test]
fn crush_public_api_preserves_paths_p_line_gfa() {
    let gfa = include_str!("test_data/crush/small_insertion.gfa");
    let before = seq_map(gfa);
    assert!(!before.is_empty(), "fixture must have at least one path");

    let result = resolve_gfa_bubbles(gfa, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on P-line GFA");

    assert!(
        result.stats.resolved >= 1,
        "expected at least 1 resolved bubble, got {:?}",
        result.stats,
    );
    assert_eq!(before, seq_map(&result.gfa));
}

#[test]
fn crush_public_api_preserves_paths_w_line_gfa() {
    let gfa = include_str!("test_data/crush/small_insertion_walks.gfa");
    let before = seq_map(gfa);
    assert!(!before.is_empty(), "fixture must have at least one walk");

    let result = resolve_gfa_bubbles(gfa, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on W-line GFA");

    assert!(
        result.stats.resolved >= 1,
        "expected at least 1 resolved bubble in W-line GFA, got {:?}",
        result.stats,
    );
    assert_eq!(before, seq_map(&result.gfa));
}

#[test]
fn crush_public_api_p_and_w_line_produce_same_sequences() {
    let gfa_p = include_str!("test_data/crush/small_insertion.gfa");
    let gfa_w = include_str!("test_data/crush/small_insertion_walks.gfa");

    let result_p = resolve_gfa_bubbles(gfa_p, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on P-line GFA");
    let result_w = resolve_gfa_bubbles(gfa_w, &ResolutionConfig::default())
        .expect("resolve_gfa_bubbles failed on W-line GFA");

    assert_eq!(
        result_p.stats.resolved, result_w.stats.resolved,
        "both GFA variants should resolve the same number of bubbles",
    );
    assert_eq!(seq_map(&result_p.gfa), seq_map(&result_w.gfa));
}
