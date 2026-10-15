//! THE SYNG NODE-ID PROJECTION ONTO THE PGGB SUBSTRATE (stage 3 of the
//! owner's approved four, the build-substrate ruling).
//!
//! The ruling's key property: the syng's syncmer nodes are interned BY
//! EXACT SUBSEQUENCE (the k-mer content of the window), so they can be
//! PROJECTED onto ANY graph built from the same panel subsequences.
//! The pggb graphs (the impg graph machinery's pggb engine: FastGA/
//! SweepGA all-pairs alignment + seqwish induction + smoothxg-style
//! block smoothing + gfaffix normalization, over the partition's
//! member-row sequences) carry their OWN node interning; this tool
//! projects the global syng ids back onto them and PROVES the
//! projection end-to-end.
//!
//! The proof is SPELL EQUALITY, the id-continuity discipline: every
//! member row's sequence reconstructed from the pggb graph's own P
//! line (segments concatenated in path order, orientation honored)
//! must be BYTE-IDENTICAL to the panel AGC row it claims to be, and
//! the row's syng walk — the raw forward-frame matched-syncmer
//! extraction of the reconstructed sequence (the same interned
//! windows the syng stores for that path range) — must spell EXACTLY
//! the signed global syncmer ids of the export graph's P line for
//! the same row (the raw syng range GFA whose S names ARE the global
//! interning). If every row spells equal, the pggb substrate's
//! coordinates ARE the syng coordinates under the projection — the
//! instrument, the anchor projection and the truth index work
//! unchanged over the projected graphs.
//!
//! Assessment-side only: this tool reads the syng, the AGC and the
//! two graph families and writes receipts; no product file, no
//! instrument change, no thresholds.

use std::collections::{BTreeMap, HashMap, HashSet};
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, BufWriter, Write};

use clap::Parser;
use impg::sequence_index::{SequenceIndex, UnifiedSequenceIndex};
use impg::syng::{SyncmerParams, SyngIndex};

#[derive(Parser)]
struct Options {
    /// The panel syng prefix.
    #[arg(long)]
    panel: String,
    /// The panel AGC archive (the sequence source of truth).
    #[arg(long)]
    agc: String,
    /// The export (raw syng range) partition-graph directory.
    #[arg(long)]
    export_graphs: std::path::PathBuf,
    /// The pggb partition-graph build root (partition<N>/graph.gfa).
    #[arg(long)]
    pggb_dir: std::path::PathBuf,
    /// The output directory for the projected receipts.
    #[arg(long)]
    out_dir: std::path::PathBuf,
    /// The partition ids (comma-separated).
    #[arg(long)]
    partitions: String,
}

fn parse_range_name(name: &str) -> Option<(String, u64, u64)> {
    let colon = name.rfind(':')?;
    let dash = name[colon + 1..].rfind('-')?;
    let start = name[colon + 1..colon + 1 + dash].parse::<u64>().ok()?;
    let end = name[colon + 2 + dash..].parse::<u64>().ok()?;
    Some((name[..colon].to_string(), start, end))
}

/// One parsed GFA: segment sequences, links and paths.
struct Gfa {
    segments: BTreeMap<String, Vec<u8>>,
    /// (from_id, from_orient, to_id, to_orient)
    links: Vec<(String, bool, String, bool)>,
    /// path name -> signed steps (id, orientation)
    paths: BTreeMap<String, Vec<(String, bool)>>,
}

fn read_gfa(path: &std::path::Path) -> io::Result<Gfa> {
    let file = File::open(path)?;
    let reader = BufReader::with_capacity(8 * 1024 * 1024, file);
    let mut gfa = Gfa {
        segments: BTreeMap::new(),
        links: Vec::new(),
        paths: BTreeMap::new(),
    };
    for line in reader.lines() {
        let line = line?;
        if let Some(rest) = line.strip_prefix("S\t") {
            let mut fields = rest.splitn(3, '\t');
            let id = fields.next().unwrap_or_default();
            let seq = fields.next().unwrap_or_default();
            let seq = seq.as_bytes().to_vec();
            if gfa.segments.insert(id.to_string(), seq).is_some() {
                return Err(io::Error::other(format!("duplicate S id {id} in {}", path.display())));
            }
        } else if let Some(rest) = line.strip_prefix("L\t") {
            let fields: Vec<&str> = rest.split('\t').collect();
            if fields.len() < 4 {
                return Err(io::Error::other("short L line"));
            }
            gfa.links.push((
                fields[0].to_string(),
                fields[1] == "+",
                fields[2].to_string(),
                fields[3] == "+",
            ));
        } else if let Some(rest) = line.strip_prefix("P\t") {
            let mut fields = rest.splitn(3, '\t');
            if let (Some(name), Some(steps)) = (fields.next(), fields.next()) {
                let mut parsed = Vec::new();
                for step in steps.split(',') {
                    let step = step.trim_end_matches('\r');
                    if step.is_empty() {
                        continue;
                    }
                    let (id, orient) = step.split_at(step.len() - 1);
                    parsed.push((id.to_string(), orient == "+"));
                }
                if gfa.paths.insert(name.to_string(), parsed).is_some() {
                    return Err(io::Error::other(format!("duplicate P name {name}")));
                }
            }
        }
    }
    Ok(gfa)
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'N' | b'n' => b'N',
            other => other,
        })
        .collect()
}

fn ensure_slice(got: &[u8], want: &[u8], msg: &str) {
    if got != want {
        panic!("{msg} ({} vs {} bytes)", got.len(), want.len());
    }
}

fn ensure_slice_walk(got: &[(i32, u64)], want: &[(i32, u64)], msg: &str) {
    if got != want {
        panic!(
            "{msg} ({} vs {} steps; first diff {:?})",
            got.len(),
            want.len(),
            got.iter().zip(want.iter()).find(|(a, b)| a != b)
        );
    }
}

/// Union-find over segment ids by name.
struct UnionFind {
    parent: HashMap<String, String>,
}

impl UnionFind {
    fn new() -> Self {
        UnionFind {
            parent: HashMap::new(),
        }
    }
    fn add(&mut self, key: &str) {
        self.parent
            .entry(key.to_string())
            .or_insert_with(|| key.to_string());
    }
    fn find(&mut self, key: &str) -> String {
        self.add(key);
        let mut root = key.to_string();
        loop {
            let next = self.parent[&root].clone();
            if next == root {
                break;
            }
            root = next;
        }
        // path compression
        let mut walked = key.to_string();
        while walked != root {
            let next = self.parent[&walked].clone();
            self.parent.insert(walked.clone(), root.clone());
            walked = next;
        }
        root
    }
    fn union(&mut self, a: &str, b: &str) {
        let ra = self.find(a);
        let rb = self.find(b);
        if ra != rb {
            self.parent.insert(ra, rb);
        }
    }
}

fn main() {
    let options = Options::parse();
    let partitions: Vec<u32> = options
        .partitions
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| s.parse().expect("partition id"))
        .collect();
    assert!(!partitions.is_empty(), "no partitions given");

    eprintln!("[pggb-projection] loading syng index from {}", options.panel);
    let panel = SyngIndex::load(&options.panel, SyncmerParams::default())
        .expect("failed to load syng index");
    let n_syncmer_nodes = panel.num_syncmer_nodes() as u64;
    eprintln!(
        "[pggb-projection] index loaded: {} syncmer nodes, syncmer span {} bp",
        n_syncmer_nodes,
        panel.syncmer_length_bp()
    );
    let sequences = vec![options.agc.clone()];
    let sequence_index =
        UnifiedSequenceIndex::from_files(&sequences).expect("failed to load AGC");

    fs::create_dir_all(&options.out_dir).expect("create output dir");
    let tables_path = options.out_dir.join("pggb-projection.tables.txt");
    let mut tables = BufWriter::new(File::create(&tables_path).expect("create tables"));

    writeln!(
        tables,
        "# per-partition structural comparison: the pggb substrate vs the export (raw syng range) substrate"
    )
    .unwrap();
    writeln!(
        tables,
        "# columns: partition members | export S syncmer_S gap_S L P | pggb S L P segment_bp revisiting_paths | pggb shared_segments shared_bp unique_bp | export shared_syncmer_nodes shared_occ | pggb bubbles_in_out | export bubbles_in_out | pggb components paths_in_largest"
    )
    .unwrap();

    let mut totals = (
        0u64, 0u64, 0u64, // members, spell_equal, spell_diff
    );

    for &partition in &partitions {
        let export_path = options
            .export_graphs
            .join(format!("partition{partition}.gfa"));
        let pggb_path = options
            .pggb_dir
            .join(format!("partition{partition}/graph.gfa"));
        let export = read_gfa(&export_path).expect("read export gfa");
        let pggb = read_gfa(&pggb_path).expect("read pggb gfa");

        let out_prefix = options.out_dir.join(format!("partition{partition}"));
        let projection_path = out_prefix.with_extension("pggb.projection.jsonl");
        let map_path = out_prefix.with_extension("pggb.map.json");

        let mut projection =
            BufWriter::new(File::create(&projection_path).expect("create projection"));
        let mut spell_equal = 0u64;
        let mut spell_diff = 0u64;
        let mut sequence_mismatch = 0u64;
        let mut rows_total = 0u64;

        // ---- the row-level proof: reconstruct, verify, intern, compare
        for (name, steps) in &pggb.paths {
            let Some((path_name, start, end)) = parse_range_name(name) else {
                panic!("pggb path name {name} is not a range name");
            };
            rows_total += 1;
            let mut seq: Vec<u8> = Vec::with_capacity(16 * 1024);
            let mut revisits = false;
            let mut seen = HashSet::new();
            for (id, forward) in steps {
                let segment = pggb
                    .segments
                    .get(id)
                    .unwrap_or_else(|| panic!("pggb path {name} steps absent segment {id}"));
                if !seen.insert(id.clone()) {
                    revisits = true;
                }
                if *forward {
                    seq.extend_from_slice(segment);
                } else {
                    seq.extend_from_slice(&revcomp(segment));
                }
            }
            for b in seq.iter_mut() {
                b.make_ascii_uppercase();
            }
            let fetched = sequence_index
                .fetch_sequence(&path_name, start as i32, end as i32)
                .unwrap_or_else(|e| panic!("AGC fetch {path_name}:{start}-{end}: {e}"));
            let mut fetched = fetched;
            for b in fetched.iter_mut() {
                b.make_ascii_uppercase();
            }
            let seq_equal = fetched == seq;
            if !seq_equal {
                sequence_mismatch += 1;
                panic!(
                    "pggb path {name}: reconstructed sequence differs from the AGC row \
                     ({} vs {} bp)",
                    seq.len(),
                    fetched.len()
                );
            }

            // the projection: the row's syng walk. The export graphs'
            // P lines spell the stored walk over [start,end) under the
            // writer's own OVERLAP convention (walk_path_range keeps a
            // window iff bp < end && bp + (params.k + params.w) > start,
            // the full scan span), so
            // the projected walk is derived over the row's CONTEXT:
            // the pggb-reconstructed interior (verified byte-identical
            // to the AGC row above) plus the panel's own flanking
            // material from the AGC — the same source the export
            // writer's stored walk spells. The interior windows are
            // cross-checked to be exactly the raw extraction of the
            // pggb graph's own reconstructed sequence.
            let overlap_len = (panel.params.k + panel.params.w) as u64;
            let path_len = sequence_index
                .get_sequence_length(&path_name)
                .unwrap_or_else(|e| panic!("AGC length {path_name}: {e}")) as u64;
            let ctx_lo = start.saturating_sub(overlap_len - 1);
            let ctx_hi = (end + overlap_len - 1).min(path_len);
            let mut ctx = sequence_index
                .fetch_sequence(&path_name, ctx_lo as i32, ctx_hi as i32)
                .unwrap_or_else(|e| panic!("AGC fetch {path_name}:{ctx_lo}-{ctx_hi}: {e}"));
            for b in ctx.iter_mut() {
                b.make_ascii_uppercase();
            }
            // the pggb-reconstructed interior must be the context's own
            // interior slice (the row's material is the panel's row)
            let interior_offset = (start - ctx_lo) as usize;
            ensure_slice(
                &ctx[interior_offset..interior_offset + seq.len()],
                &seq,
                &format!("pggb path {name}: interior mismatch vs AGC context"),
            );
            let ctx_walk: Vec<(i32, u64)> =
                impg::genome_inference::mem_records::raw_matched_syncmers(&panel, &ctx)
                    .expect("raw matched syncmers");
            let mut walk: Vec<(i32, u64)> = ctx_walk
                .into_iter()
                .map(|(signed, pos)| (signed, ctx_lo + pos))
                .filter(|&(_signed, bp)| bp < end && bp + overlap_len > start)
                .collect();
            walk.sort_unstable_by_key(|&(_signed, bp)| bp);
            let spelled: Vec<i32> = walk.iter().map(|&(signed, _)| signed).collect();
            // the interior windows must be EXACTLY the raw extraction of
            // the pggb graph's own reconstructed sequence (the walk goes
            // THROUGH the pggb graph, not around it)
            let interior_walk: Vec<(i32, u64)> =
                impg::genome_inference::mem_records::raw_matched_syncmers(&panel, &seq)
                    .expect("raw matched syncmers");
            let interior_of_walk: Vec<(i32, u64)> = walk
                .iter()
                .copied()
                .filter(|&(_signed, bp)| bp >= start && bp + overlap_len <= end)
                .collect();
            let interior_rebased: Vec<(i32, u64)> = interior_walk
                .into_iter()
                .map(|(signed, pos)| (signed, start + pos))
                .collect();
            ensure_slice_walk(
                &interior_rebased,
                &interior_of_walk,
                &format!("pggb path {name}: interior walk mismatch"),
            );

            // the export substrate's spelling of the same row: signed
            // global syncmer ids (gap segments > n_syncmer_nodes are
            // render-local splices, excluded).
            let export_steps: Vec<i32> = export
                .paths
                .get(name)
                .unwrap_or_else(|| panic!("export gfa lacks path {name}"))
                .iter()
                .filter_map(|(id, forward)| {
                    let abs: u64 = id.parse().expect("numeric export S id");
                    if abs > n_syncmer_nodes {
                        None
                    } else {
                        Some(if *forward { abs as i32 } else { -(abs as i32) })
                    }
                })
                .collect();
            let equal = spelled == export_steps;
            if equal {
                spell_equal += 1;
            } else {
                spell_diff += 1;
                eprintln!(
                    "[pggb-projection] partition {partition} path {name}: SPELL DIFF \
                     (pggb-projected {} steps, export {} steps)",
                    spelled.len(),
                    export_steps.len()
                );
            }

            serde_json_ish::write_path(
                &mut projection,
                name,
                &path_name,
                start,
                end,
                seq.len(),
                seq_equal,
                equal,
                &walk,
                revisits,
            )
            .expect("write path projection");
        }
        projection.flush().expect("flush projection");

        totals.0 += rows_total;
        totals.1 += spell_equal;
        totals.2 += spell_diff;

        // ---- the segment-level projection map
        let mut segment_lines = 0u64;
        let mut segment_windows = 0u64;
        {
            let mut seg_projection =
                BufWriter::new(File::create(out_prefix.with_extension("pggb.segments.jsonl")).expect("create segments"));
            for (id, seq) in &pggb.segments {
                let walk: Vec<(i32, u64)> =
                    impg::genome_inference::mem_records::raw_matched_syncmers(&panel, seq)
                        .expect("raw matched syncmers");
                segment_windows += walk.len() as u64;
                if !walk.is_empty() {
                    segment_lines += 1;
                }
                serde_json_ish::write_segment(&mut seg_projection, partition, id, seq.len(), &walk)
                    .expect("write segment projection");
            }
            seg_projection.flush().expect("flush segments");
        }

        // ---- the structural comparison (receipt-side tables)
        // export counts
        let mut export_gap_segments = 0u64;
        let mut export_syncmer_segments = 0u64;
        for id in export.segments.keys() {
            let abs: u64 = id.parse().expect("numeric export S id");
            if abs > n_syncmer_nodes {
                export_gap_segments += 1;
            } else {
                export_syncmer_segments += 1;
            }
        }
        // sharing: how many distinct paths visit each node/segment
        let mut pggb_path_count: HashMap<&String, u64> = HashMap::new();
        for (_name, steps) in &pggb.paths {
            let mut visited: HashSet<&String> = HashSet::new();
            for (id, _) in steps {
                visited.insert(id);
            }
            for id in visited {
                *pggb_path_count.entry(id).or_insert(0) += 1;
            }
        }
        let mut shared_segments = 0u64;
        let mut shared_bp = 0u64;
        let mut unique_bp = 0u64;
        for (id, k_paths) in &pggb_path_count {
            let len = pggb.segments[id.as_str()].len() as u64;
            if *k_paths >= 2 {
                shared_segments += 1;
                shared_bp += len;
            } else {
                unique_bp += len;
            }
        }
        let mut export_path_count: HashMap<u64, u64> = HashMap::new();
        for (_name, steps) in &export.paths {
            let mut visited: HashSet<u64> = HashSet::new();
            for (id, _) in steps {
                let abs: u64 = id.parse().expect("numeric export S id");
                if abs <= n_syncmer_nodes {
                    visited.insert(abs);
                }
            }
            for abs in visited {
                *export_path_count.entry(abs).or_insert(0) += 1;
            }
        }
        let mut export_shared_nodes = 0u64;
        let mut export_shared_occ = 0u64;
        for (&abs, &k_paths) in &export_path_count {
            if k_paths >= 2 {
                export_shared_nodes += 1;
                export_shared_occ += k_paths;
            }
        }
        // bubbles: segments with in-degree>1 or out-degree>1
        let mut pggb_out: HashMap<&str, HashSet<&str>> = HashMap::new();
        let mut pggb_in: HashMap<&str, HashSet<&str>> = HashMap::new();
        for (a, _ao, b, _bo) in &pggb.links {
            pggb_out.entry(a.as_str()).or_default().insert(b.as_str());
            pggb_in.entry(b.as_str()).or_default().insert(a.as_str());
        }
        let pggb_bubbles = pggb_out
            .values()
            .chain(pggb_in.values())
            .filter(|s| s.len() > 1)
            .count() as u64;
        let mut export_out: HashMap<&str, HashSet<&str>> = HashMap::new();
        let mut export_in: HashMap<&str, HashSet<&str>> = HashMap::new();
        for (a, _ao, b, _bo) in &export.links {
            export_out.entry(a.as_str()).or_default().insert(b.as_str());
            export_in.entry(b.as_str()).or_default().insert(a.as_str());
        }
        let export_bubbles = export_out
            .values()
            .chain(export_in.values())
            .filter(|s| s.len() > 1)
            .count() as u64;
        // components (union-find over pggb links)
        let mut uf = UnionFind::new();
        for id in pggb.segments.keys() {
            uf.add(id);
        }
        for (a, _ao, b, _bo) in &pggb.links {
            uf.union(a, b);
        }
        let mut comp_paths: HashMap<String, u64> = HashMap::new();
        for (name, steps) in &pggb.paths {
            let Some((_, _, _)) = parse_range_name(name) else {
                unreachable!()
            };
            let root = uf.find(steps.first().map(|s| s.0.as_str()).unwrap_or(""));
            *comp_paths.entry(root).or_insert(0) += 1;
        }
        let n_components = comp_paths.keys().count() as u64;
        let paths_in_largest = comp_paths.values().copied().max().unwrap_or(0);
        let revisiting_paths = pggb
            .paths
            .values()
            .filter(|steps| {
                let mut seen = HashSet::new();
                steps.iter().any(|(id, _)| !seen.insert(id))
            })
            .count() as u64;
        let segment_bp: u64 = pggb.segments.values().map(|s| s.len() as u64).sum();

        writeln!(
            tables,
            "partition{partition}\t{}\tS={}\tsyncmerS={}\tgapS={}\tL={}\tP={}\tS={}\tL={}\tP={}\tbp={}\trevisit={}\tsharedS={}\tsharedbp={}\tuniquebp={}\texportSharedNodes={}\texportSharedOcc={}\tpggbBubbles={}\texportBubbles={}\tcomponents={}\tlargestCompPaths={}",
            pggb.paths.len(),
            export.segments.len(),
            export_syncmer_segments,
            export_gap_segments,
            export.links.len(),
            export.paths.len(),
            pggb.segments.len(),
            pggb.links.len(),
            pggb.paths.len(),
            segment_bp,
            revisiting_paths,
            shared_segments,
            shared_bp,
            unique_bp,
            export_shared_nodes,
            export_shared_occ,
            pggb_bubbles,
            export_bubbles,
            n_components,
            paths_in_largest,
        )
        .unwrap();

        // ---- the projected map.json (the scorer interface, same schema)
        let bed_path = options.pggb_dir.join(format!("partition{partition}/members.bed"));
        let bed_text = fs::read_to_string(&bed_path).expect("read members.bed");
        let members: Vec<serde_json_ish::Member> = bed_text
            .lines()
            .filter(|l| !l.trim().is_empty())
            .map(|l| {
                let f: Vec<&str> = l.split('\t').collect();
                serde_json_ish::Member {
                    path_name: f[0].to_string(),
                    start: f[1].parse().expect("bed start"),
                    end: f[2].parse().expect("bed end"),
                }
            })
            .collect();
        // THE COMPLETENESS CHECK (the substrate must carry every member
        // row of the locality: the scorer's candidate domain is the map's
        // member rows, so a pggb build that dropped or renamed a row
        // would silently shrink the domain — fail loudly instead).
        let expected_names: HashSet<String> = members
            .iter()
            .map(|m| format!("{}:{}-{}", m.path_name, m.start, m.end))
            .collect();
        let observed_names: HashSet<String> = pggb.paths.keys().cloned().collect();
        for missing in expected_names.difference(&observed_names) {
            panic!(
                "partition {partition}: members.bed row {missing} has no pggb path"
            );
        }
        for extra in observed_names.difference(&expected_names) {
            panic!(
                "partition {partition}: pggb path {extra} is not a members.bed row"
            );
        }
        let map = serde_json_ish::Map {
            partition,
            construction_material:
                "impg partition membership BED (alignment-induced; syng chain criteria, commit 295bca9)",
            construction_graph:
                "impg graph --gfa-engine pggb (FastGA/SweepGA + seqwish + smooth + gfaffix) over the members' AGC-restricted sequences",
            n_syncmer_nodes,
            members,
            counts_rows: rows_total,
            spell_equal,
            spell_diff,
            sequence_mismatch,
            segment_lines,
            segment_windows,
        };
        map.write(&map_path).expect("write map");

        eprintln!(
            "[pggb-projection] partition {partition}: {} rows, spell_equal {} spell_diff {}, \
             pggb S={} L={} P={}, export S={} L={}",
            rows_total,
            spell_equal,
            spell_diff,
            pggb.segments.len(),
            pggb.links.len(),
            pggb.paths.len(),
            export.segments.len(),
            export.links.len(),
        );
    }

    tables.flush().expect("flush tables");
    eprintln!(
        "[pggb-projection] TOTALS: {} rows, spell_equal {}, spell_diff {}",
        totals.0, totals.1, totals.2
    );
    if totals.2 > 0 {
        panic!(
            "SPELL EQUALITY FAILED at {} of {} rows — the projection is not proven",
            totals.2, totals.0
        );
    }
}

// Small local JSON emitters (no serde dependency games in examples).
mod serde_json_ish {
    use std::io::{self, BufWriter, Write};

    pub struct Member {
        pub path_name: String,
        pub start: u64,
        pub end: u64,
    }

    pub struct Map<'a> {
        pub partition: u32,
        pub construction_material: &'a str,
        pub construction_graph: &'a str,
        pub n_syncmer_nodes: u64,
        pub members: Vec<Member>,
        pub counts_rows: u64,
        pub spell_equal: u64,
        pub spell_diff: u64,
        pub sequence_mismatch: u64,
        pub segment_lines: u64,
        pub segment_windows: u64,
    }

    impl Map<'_> {
        pub fn write(&self, path: &std::path::Path) -> io::Result<()> {
            let mut out = BufWriter::new(std::fs::File::create(path)?);
            writeln!(out, "{{")?;
            writeln!(
                out,
                "  \"partition\": {},",
                self.partition
            )?;
            writeln!(
                out,
                "  \"construction\": {{\n    \"material\": \"{}\",\n    \"graph\": \"{}\",\n    \"projection\": \"syng node ids projected by exact-subsequence interning (raw matched syncmers of each path's reconstructed sequence); spell-equality proven vs the export graphs\"\n  }},",
                self.construction_material, self.construction_graph
            )?;
            writeln!(out, "  \"n_syncmer_nodes\": {},", self.n_syncmer_nodes)?;
            writeln!(
                out,
                "  \"node_id_mapping\": {{\n    \"gfa_nodes\": \"pggb segments carry the build's OWN interning; the syng projection lives in the .pggb.projection.jsonl sidecar\",\n    \"path_steps\": \"member walks through the pggb graph spell the same signed global syncmer ids as the export graphs (the spell-equality proof)\"\n  }},"
            )?;
            writeln!(
                out,
                "  \"counts\": {{\n    \"members\": {},\n    \"rows\": {},\n    \"spell_equal\": {},\n    \"spell_diff\": {},\n    \"sequence_mismatch\": {},\n    \"segments_with_projection\": {},\n    \"segment_windows\": {}\n  }},",
                self.members.len(),
                self.counts_rows,
                self.spell_equal,
                self.spell_diff,
                self.sequence_mismatch,
                self.segment_lines,
                self.segment_windows
            )?;
            writeln!(out, "  \"members\": [")?;
            for (i, m) in self.members.iter().enumerate() {
                writeln!(
                    out,
                    "    {{\"path_name\": \"{}\", \"start\": {}, \"end\": {}}}{}",
                    m.path_name,
                    m.start,
                    m.end,
                    if i + 1 == self.members.len() { "" } else { "," }
                )?;
            }
            writeln!(out, "  ]")?;
            writeln!(out, "}}")?;
            out.flush()
        }
    }

    pub fn write_path<W: Write>(
        out: &mut W,
        name: &str,
        path_name: &str,
        start: u64,
        end: u64,
        length: usize,
        seq_equal: bool,
        spell_equal: bool,
        walk: &[(i32, u64)],
        revisits: bool,
    ) -> io::Result<()> {
        write!(
            out,
            "{{\"p\":\"{}\",\"n\":\"{}\",\"s\":{},\"e\":{},\"len\":{},\"seq\":{},\"spell\":{},\"rev\":{},\"w\":[",
            name, path_name, start, end, length, seq_equal, spell_equal, revisits
        )?;
        for (i, (signed, pos)) in walk.iter().enumerate() {
            if i > 0 {
                write!(out, ",")?;
            }
            write!(out, "[{},{}]", signed, pos)?;
        }
        writeln!(out, "]}}")
    }

    pub fn write_segment<W: Write>(
        out: &mut W,
        partition: u32,
        id: &str,
        length: usize,
        walk: &[(i32, u64)],
    ) -> io::Result<()> {
        write!(
            out,
            "{{\"part\":{},\"seg\":\"{}\",\"len\":{},\"w\":[",
            partition, id, length
        )?;
        for (i, (signed, pos)) in walk.iter().enumerate() {
            if i > 0 {
                write!(out, ",")?;
            }
            write!(out, "[{},{}]", signed, pos)?;
        }
        writeln!(out, "]}}")
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

/// Range names parse to (path, start, end) and reject malformed forms.
#[test]
fn range_names_parse_and_reject() {
    assert_eq!(
        parse_range_name("S288C#0#chrI:23809-29028"),
        Some(("S288C#0#chrI".to_string(), 23809, 29028))
    );
    assert_eq!(
        parse_range_name("BTE#3#block28_contig1:40000-50000"),
        Some(("BTE#3#block28_contig1".to_string(), 40000, 50000))
    );
    assert_eq!(parse_range_name("no-range"), None);
    assert_eq!(parse_range_name("x:a-b"), None);
}

/// The GFA reader parses S/L/P with orientation and tolerates extra
/// columns; duplicate ids are rejected loudly.
#[test]
fn gfa_reader_parses_lines() {
    let dir = std::env::temp_dir().join("pggb_projection_test_gfa");
    std::fs::create_dir_all(&dir).unwrap();
    let path = dir.join("t.gfa");
    std::fs::write(
        &path,
        "H\tVN:Z:1.0\nS\t1\tACGTACGT\nS\t2\tTTTTAAAA\tLN:i:8\nL\t1\t+\t2\t-\t0M\nP\tp:0-16\t1+,2-\t*\n",
    )
    .unwrap();
    let gfa = read_gfa(&path).unwrap();
    assert_eq!(gfa.segments["1"], b"ACGTACGT".to_vec());
    assert_eq!(gfa.segments["2"].len(), 8);
    assert_eq!(gfa.links.len(), 1);
    assert_eq!(
        gfa.paths["p:0-16"],
        vec![("1".to_string(), true), ("2".to_string(), false)]
    );
    std::fs::write(&path, "S\t1\tA\nS\t1\tC\n").unwrap();
    assert!(read_gfa(&path).is_err());
}

/// The reverse complement used for '-'-oriented path steps.
#[test]
fn revcomp_handles_orientation() {
    assert_eq!(revcomp(b"ACGTNN"), b"NNACGT".to_vec());
    assert_eq!(revcomp(b"acgt"), b"ACGT".to_vec());
}

/// Union-find connects linked segments and reports one root per
/// component (the alignment-induced connectivity measurement).
#[test]
fn union_find_connects_components() {
    let mut uf = UnionFind::new();
    uf.union("a", "b");
    uf.union("b", "c");
    uf.union("d", "d");
    assert_eq!(uf.find("a"), uf.find("c"));
    assert_ne!(uf.find("a"), uf.find("d"));
    uf.find("e");
    assert_eq!(uf.find("e"), "e".to_string());
}
