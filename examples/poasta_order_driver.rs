use impg::graph_report::{describe_gfa, GraphReportOptions};
use impg::resolution::{abpoa_sequences_to_gfa_in_order, poasta_sequences_to_gfa_in_order};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::env;
use std::fs;
use std::io;
use std::path::{Path, PathBuf};
use std::process::Command;

const SCORING: (u8, u8, u8, u8, u8, u8) = (1, 4, 6, 2, 26, 1);
const KMER_SIZE: usize = 31;

#[derive(Clone, Debug)]
struct Args {
    input: PathBuf,
    output_dir: PathBuf,
    abpoa_bin: Option<String>,
    render_svg: bool,
    random_orders: usize,
}

#[derive(Clone, Debug)]
struct Link {
    from_segment: String,
    from_orient: char,
    to_segment: String,
    to_orient: char,
}

#[derive(Clone, Debug)]
struct GfaData {
    segments: BTreeMap<String, Vec<u8>>,
    links: Vec<Link>,
    paths: Vec<(String, Vec<(String, char)>)>,
}

#[derive(Clone, Debug)]
struct Experiment {
    method: Method,
    label: String,
    order: Vec<usize>,
}

#[derive(Clone, Copy, Debug)]
enum Method {
    Poasta,
    Abpoa,
}

#[derive(Clone, Debug)]
struct CycleSummary {
    directed_sccs_gt1: usize,
    directed_self_loop_edges: usize,
    segment_self_loop_links: usize,
    largest_directed_scc: usize,
    cyclic_oriented_nodes: usize,
}

#[derive(Clone, Debug)]
struct PreservationSummary {
    missing_paths: usize,
    extra_paths: usize,
    changed_paths: usize,
}

fn main() -> io::Result<()> {
    let args = parse_args()?;
    fs::create_dir_all(&args.output_dir)?;
    let gfa_dir = args.output_dir.join("gfas");
    let render_dir = args.output_dir.join("renders");
    fs::create_dir_all(&gfa_dir)?;
    if args.render_svg {
        fs::create_dir_all(&render_dir)?;
    }

    let input_text = fs::read_to_string(&args.input)?;
    let input_gfa = parse_gfa(&input_text)?;
    let input_sequences = path_sequences(&input_gfa)?;
    let headers = input_sequences
        .iter()
        .map(|(name, _)| name.clone())
        .collect::<Vec<_>>();
    let sequences = input_sequences
        .iter()
        .map(|(_, sequence)| sequence.clone())
        .collect::<Vec<_>>();
    write_fasta(
        &args.output_dir.join("input_sequences.fa"),
        &input_sequences,
    )?;

    let baseline = longest_then_name_order(&headers, &sequences);
    let medoid = medoid_index(&sequences, &headers);
    let guide = nearest_neighbor_order(&sequences, &headers, medoid);
    let mut experiments = Vec::new();
    experiments.push(Experiment {
        method: Method::Poasta,
        label: "longest_then_name".to_string(),
        order: baseline.clone(),
    });
    experiments.push(Experiment {
        method: Method::Poasta,
        label: "reverse_longest_then_name".to_string(),
        order: baseline.iter().rev().copied().collect(),
    });
    experiments.push(Experiment {
        method: Method::Poasta,
        label: "medoid_first_then_longest".to_string(),
        order: medoid_first_order(medoid, &baseline),
    });
    experiments.push(Experiment {
        method: Method::Poasta,
        label: "nearest_neighbor_guidetree".to_string(),
        order: guide.clone(),
    });

    let mut skipped = Vec::new();
    if let Some(reference_order) = reference_first_order(&headers, &baseline) {
        experiments.push(Experiment {
            method: Method::Poasta,
            label: "reference_chm13_first".to_string(),
            order: reference_order,
        });
    } else {
        skipped.push((
            "reference_chm13_first",
            "no CHM13 or GRCh38 path was present in the 47 focused paths",
        ));
    }

    for seed in 0..args.random_orders {
        experiments.push(Experiment {
            method: Method::Poasta,
            label: format!("random_seed_{seed:02}"),
            order: deterministic_shuffle(headers.len(), seed as u64),
        });
    }

    if args.abpoa_bin.is_some() {
        experiments.push(Experiment {
            method: Method::Abpoa,
            label: "longest_then_name".to_string(),
            order: baseline.clone(),
        });
        experiments.push(Experiment {
            method: Method::Abpoa,
            label: "nearest_neighbor_guidetree".to_string(),
            order: guide.clone(),
        });
    } else {
        skipped.push(("abpoa", "no --abpoa-bin supplied"));
    }

    let mut order_lines = vec!["method\tlabel\trank\tpath\tlength".to_string()];
    let mut summary_lines = vec![summary_header()];

    for experiment in &experiments {
        for (rank, &idx) in experiment.order.iter().enumerate() {
            order_lines.push(format!(
                "{}\t{}\t{}\t{}\t{}",
                method_name(experiment.method),
                experiment.label,
                rank + 1,
                headers[idx],
                sequences[idx].len()
            ));
        }

        let gfa = match experiment.method {
            Method::Poasta => {
                poasta_sequences_to_gfa_in_order(&headers, &sequences, &experiment.order, SCORING)?
            }
            Method::Abpoa => abpoa_sequences_to_gfa_in_order(
                &headers,
                &sequences,
                &experiment.order,
                SCORING,
                args.abpoa_bin
                    .as_deref()
                    .expect("abPOA experiment needs binary"),
            )?,
        };
        let file_stem = format!(
            "{}_{}",
            method_name(experiment.method).to_ascii_lowercase(),
            sanitize_label(&experiment.label)
        );
        let gfa_path = gfa_dir.join(format!("{file_stem}.gfa"));
        fs::write(&gfa_path, &gfa)?;
        let svg_path = if args.render_svg {
            let path = render_dir.join(format!("{file_stem}.svg"));
            run_gfalook(&gfa_path, &path)?;
            path.display().to_string()
        } else {
            String::new()
        };

        let output_gfa = parse_gfa(&gfa)?;
        let preservation = compare_path_sequences(&input_sequences, &path_sequences(&output_gfa)?);
        let cycles = summarize_cycles(&output_gfa);
        let report = describe_gfa(&experiment.label, &gfa, &GraphReportOptions::default())?;
        let m = report.metrics;
        summary_lines.push(format!(
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            method_name(experiment.method),
            experiment.label,
            gfa_path.display(),
            svg_path,
            preservation.missing_paths,
            preservation.extra_paths,
            preservation.changed_paths,
            m.segments,
            m.links,
            m.paths,
            m.path_steps,
            m.total_segment_bp,
            cycles.directed_sccs_gt1,
            cycles.directed_self_loop_edges,
            cycles.segment_self_loop_links,
            cycles.largest_directed_scc,
            cycles.cyclic_oriented_nodes,
            m.node_coverage_bp_weighted_mean,
            m.segment_white_space_bp_fraction,
            m.segment_white_space_bp_total,
            m.path_white_space_bp_p95,
            m.path_white_space_bp_p99,
            m.path_white_space_bp_max,
            m.path_white_space_bridges_ge_threshold,
            m.path_depth_median,
            m.path_depth_p95,
            m.path_depth_max,
            headers[experiment.order[0]],
            sequences[experiment.order[0]].len(),
        ));
    }

    fs::write(
        args.output_dir.join("orders.tsv"),
        format!("{}\n", order_lines.join("\n")),
    )?;
    fs::write(
        args.output_dir.join("summary.tsv"),
        format!("{}\n", summary_lines.join("\n")),
    )?;
    let skipped_lines = std::iter::once("label\treason".to_string())
        .chain(
            skipped
                .into_iter()
                .map(|(label, reason)| format!("{label}\t{reason}")),
        )
        .collect::<Vec<_>>();
    fs::write(
        args.output_dir.join("skipped.tsv"),
        format!("{}\n", skipped_lines.join("\n")),
    )?;
    Ok(())
}

fn parse_args() -> io::Result<Args> {
    let mut input = None;
    let mut output_dir = None;
    let mut abpoa_bin = None;
    let mut render_svg = false;
    let mut random_orders = 5usize;
    let mut iter = env::args().skip(1);
    while let Some(arg) = iter.next() {
        match arg.as_str() {
            "--input" => input = iter.next().map(PathBuf::from),
            "--output-dir" => output_dir = iter.next().map(PathBuf::from),
            "--abpoa-bin" => abpoa_bin = iter.next(),
            "--render-svg" => render_svg = true,
            "--random-orders" => {
                let value = iter.next().ok_or_else(|| {
                    io::Error::new(io::ErrorKind::InvalidInput, "--random-orders needs a value")
                })?;
                random_orders = value.parse().map_err(|err| {
                    io::Error::new(
                        io::ErrorKind::InvalidInput,
                        format!("invalid --random-orders value {value:?}: {err}"),
                    )
                })?;
            }
            "-h" | "--help" => {
                print_usage();
                std::process::exit(0);
            }
            other => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidInput,
                    format!("unknown argument {other:?}"),
                ));
            }
        }
    }
    Ok(Args {
        input: input
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "--input is required"))?,
        output_dir: output_dir.ok_or_else(|| {
            io::Error::new(io::ErrorKind::InvalidInput, "--output-dir is required")
        })?,
        abpoa_bin,
        render_svg,
        random_orders,
    })
}

fn print_usage() {
    eprintln!(
        "usage: cargo run --release --example poasta_order_driver -- --input IN.gfa --output-dir OUT [--abpoa-bin PATH] [--render-svg] [--random-orders N]"
    );
}

fn parse_gfa(text: &str) -> io::Result<GfaData> {
    let mut segments = BTreeMap::new();
    let mut links = Vec::new();
    let mut paths = Vec::new();
    for (line_number, line) in text.lines().enumerate() {
        if line.is_empty() {
            continue;
        }
        let fields = line.split('\t').collect::<Vec<_>>();
        match fields.first().copied() {
            Some("S") if fields.len() >= 3 => {
                segments.insert(fields[1].to_string(), fields[2].as_bytes().to_vec());
            }
            Some("L") if fields.len() >= 5 => {
                links.push(Link {
                    from_segment: fields[1].to_string(),
                    from_orient: parse_orient(fields[2], line_number + 1)?,
                    to_segment: fields[3].to_string(),
                    to_orient: parse_orient(fields[4], line_number + 1)?,
                });
            }
            Some("P") if fields.len() >= 3 => {
                paths.push((fields[1].to_string(), parse_path_steps(fields[2])?));
            }
            Some("H") => {}
            Some(record) => {
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "unsupported or malformed GFA record {record:?} at line {}",
                        line_number + 1
                    ),
                ));
            }
            None => {}
        }
    }
    Ok(GfaData {
        segments,
        links,
        paths,
    })
}

fn parse_orient(text: &str, line_number: usize) -> io::Result<char> {
    match text {
        "+" => Ok('+'),
        "-" => Ok('-'),
        _ => Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("invalid orientation {text:?} at line {line_number}"),
        )),
    }
}

fn parse_path_steps(text: &str) -> io::Result<Vec<(String, char)>> {
    if text.is_empty() || text == "*" {
        return Ok(Vec::new());
    }
    text.split(',')
        .map(|token| {
            let (name, orient) = token.split_at(token.len().saturating_sub(1));
            Ok((name.to_string(), parse_orient(orient, 0)?))
        })
        .collect()
}

fn path_sequences(gfa: &GfaData) -> io::Result<Vec<(String, Vec<u8>)>> {
    gfa.paths
        .iter()
        .map(|(name, steps)| Ok((name.clone(), path_sequence(gfa, steps)?)))
        .collect()
}

fn path_sequence(gfa: &GfaData, steps: &[(String, char)]) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    for (segment, orient) in steps {
        let sequence = gfa.segments.get(segment).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("path references missing segment {segment:?}"),
            )
        })?;
        if *orient == '+' {
            out.extend_from_slice(sequence);
        } else {
            out.extend(reverse_complement(sequence));
        }
    }
    Ok(out)
}

fn reverse_complement(sequence: &[u8]) -> Vec<u8> {
    sequence
        .iter()
        .rev()
        .map(|base| match base {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            b'N' | b'n' => b'N',
            other => *other,
        })
        .collect()
}

fn write_fasta(path: &Path, records: &[(String, Vec<u8>)]) -> io::Result<()> {
    let mut out = String::new();
    for (name, sequence) in records {
        out.push('>');
        out.push_str(name);
        out.push('\n');
        out.push_str(&String::from_utf8_lossy(sequence));
        out.push('\n');
    }
    fs::write(path, out)
}

fn longest_then_name_order(headers: &[String], sequences: &[Vec<u8>]) -> Vec<usize> {
    let mut order = (0..headers.len()).collect::<Vec<_>>();
    order.sort_by(|&a, &b| {
        sequences[b]
            .len()
            .cmp(&sequences[a].len())
            .then_with(|| headers[a].cmp(&headers[b]))
    });
    order
}

fn medoid_first_order(medoid: usize, baseline: &[usize]) -> Vec<usize> {
    std::iter::once(medoid)
        .chain(baseline.iter().copied().filter(|&idx| idx != medoid))
        .collect()
}

fn reference_first_order(headers: &[String], baseline: &[usize]) -> Option<Vec<usize>> {
    let reference_idx = baseline.iter().copied().find(|&idx| {
        let name = &headers[idx];
        name.contains("CHM13") || name.contains("GRCh38")
    })?;
    Some(
        std::iter::once(reference_idx)
            .chain(baseline.iter().copied().filter(|&idx| idx != reference_idx))
            .collect(),
    )
}

fn medoid_index(sequences: &[Vec<u8>], headers: &[String]) -> usize {
    let sketches = sequences
        .iter()
        .map(|sequence| kmer_set(sequence, KMER_SIZE))
        .collect::<Vec<_>>();
    (0..sequences.len())
        .min_by(|&a, &b| {
            let da = total_distance(a, &sketches);
            let db = total_distance(b, &sketches);
            da.total_cmp(&db).then_with(|| headers[a].cmp(&headers[b]))
        })
        .unwrap_or(0)
}

fn nearest_neighbor_order(sequences: &[Vec<u8>], headers: &[String], start: usize) -> Vec<usize> {
    if sequences.is_empty() {
        return Vec::new();
    }
    let sketches = sequences
        .iter()
        .map(|sequence| kmer_set(sequence, KMER_SIZE))
        .collect::<Vec<_>>();
    let mut used = vec![false; sequences.len()];
    used[start] = true;
    let mut order = vec![start];
    while order.len() < sequences.len() {
        let next = (0..sequences.len())
            .filter(|&idx| !used[idx])
            .min_by(|&a, &b| {
                let da = order
                    .iter()
                    .map(|&placed| jaccard_distance(&sketches[a], &sketches[placed]))
                    .fold(f64::INFINITY, f64::min);
                let db = order
                    .iter()
                    .map(|&placed| jaccard_distance(&sketches[b], &sketches[placed]))
                    .fold(f64::INFINITY, f64::min);
                da.total_cmp(&db)
                    .then_with(|| sequences[b].len().cmp(&sequences[a].len()))
                    .then_with(|| headers[a].cmp(&headers[b]))
            })
            .expect("unused sequence exists");
        used[next] = true;
        order.push(next);
    }
    order
}

fn total_distance(idx: usize, sketches: &[HashSet<Vec<u8>>]) -> f64 {
    sketches
        .iter()
        .enumerate()
        .filter(|(other, _)| *other != idx)
        .map(|(_, sketch)| jaccard_distance(&sketches[idx], sketch))
        .sum()
}

fn kmer_set(sequence: &[u8], k: usize) -> HashSet<Vec<u8>> {
    if sequence.len() <= k {
        return std::iter::once(sequence.to_vec()).collect();
    }
    sequence.windows(k).map(|window| window.to_vec()).collect()
}

fn jaccard_distance(a: &HashSet<Vec<u8>>, b: &HashSet<Vec<u8>>) -> f64 {
    if a.is_empty() && b.is_empty() {
        return 0.0;
    }
    let intersection = a.intersection(b).count();
    let union = a.len() + b.len() - intersection;
    1.0 - (intersection as f64 / union.max(1) as f64)
}

fn deterministic_shuffle(len: usize, seed: u64) -> Vec<usize> {
    let mut order = (0..len).collect::<Vec<_>>();
    let mut state = seed ^ 0x9e3779b97f4a7c15;
    for i in (1..len).rev() {
        let value = splitmix64(&mut state);
        let j = (value as usize) % (i + 1);
        order.swap(i, j);
    }
    order
}

fn splitmix64(state: &mut u64) -> u64 {
    *state = state.wrapping_add(0x9e3779b97f4a7c15);
    let mut z = *state;
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
    z ^ (z >> 31)
}

fn compare_path_sequences(
    expected: &[(String, Vec<u8>)],
    observed: &[(String, Vec<u8>)],
) -> PreservationSummary {
    let expected_map = expected.iter().cloned().collect::<BTreeMap<_, _>>();
    let observed_map = observed.iter().cloned().collect::<BTreeMap<_, _>>();
    let missing_paths = expected_map
        .keys()
        .filter(|name| !observed_map.contains_key(*name))
        .count();
    let extra_paths = observed_map
        .keys()
        .filter(|name| !expected_map.contains_key(*name))
        .count();
    let changed_paths = expected_map
        .iter()
        .filter(|(name, expected_sequence)| {
            observed_map
                .get(*name)
                .is_some_and(|observed_sequence| observed_sequence != *expected_sequence)
        })
        .count();
    PreservationSummary {
        missing_paths,
        extra_paths,
        changed_paths,
    }
}

fn summarize_cycles(gfa: &GfaData) -> CycleSummary {
    let mut nodes = BTreeSet::new();
    for segment in gfa.segments.keys() {
        nodes.insert((segment.clone(), '+'));
        nodes.insert((segment.clone(), '-'));
    }
    let mut adjacency: HashMap<(String, char), Vec<(String, char)>> = HashMap::new();
    let mut directed_self_loop_edges = 0usize;
    let mut segment_self_loop_links = 0usize;
    for (source, target, link) in directed_edges(gfa) {
        nodes.insert(source.clone());
        nodes.insert(target.clone());
        adjacency
            .entry(source.clone())
            .or_default()
            .push(target.clone());
        if source == target {
            directed_self_loop_edges += 1;
        }
        if link.from_segment == link.to_segment {
            segment_self_loop_links += 1;
        }
    }
    let components = tarjan_scc(nodes.iter().cloned(), &adjacency);
    let sccs_gt1 = components
        .iter()
        .filter(|component| component.len() > 1)
        .collect::<Vec<_>>();
    let mut cyclic_nodes = BTreeSet::new();
    for (source, target, _) in directed_edges(gfa) {
        if source == target {
            cyclic_nodes.insert(source);
        }
    }
    for component in &sccs_gt1 {
        for node in component.iter() {
            cyclic_nodes.insert((*node).clone());
        }
    }
    CycleSummary {
        directed_sccs_gt1: sccs_gt1.len(),
        directed_self_loop_edges,
        segment_self_loop_links,
        largest_directed_scc: sccs_gt1
            .iter()
            .map(|component| component.len())
            .max()
            .unwrap_or(0),
        cyclic_oriented_nodes: cyclic_nodes.len(),
    }
}

fn directed_edges(gfa: &GfaData) -> Vec<((String, char), (String, char), Link)> {
    let mut edges = Vec::new();
    for link in &gfa.links {
        edges.push((
            (link.from_segment.clone(), link.from_orient),
            (link.to_segment.clone(), link.to_orient),
            link.clone(),
        ));
        edges.push((
            (link.to_segment.clone(), reverse_orient(link.to_orient)),
            (link.from_segment.clone(), reverse_orient(link.from_orient)),
            link.clone(),
        ));
    }
    edges
}

fn reverse_orient(orient: char) -> char {
    if orient == '+' {
        '-'
    } else {
        '+'
    }
}

fn tarjan_scc<I>(
    nodes: I,
    adjacency: &HashMap<(String, char), Vec<(String, char)>>,
) -> Vec<Vec<(String, char)>>
where
    I: IntoIterator<Item = (String, char)>,
{
    struct Tarjan<'a> {
        index: usize,
        stack: Vec<(String, char)>,
        on_stack: HashSet<(String, char)>,
        indices: HashMap<(String, char), usize>,
        lowlinks: HashMap<(String, char), usize>,
        components: Vec<Vec<(String, char)>>,
        adjacency: &'a HashMap<(String, char), Vec<(String, char)>>,
    }

    impl Tarjan<'_> {
        fn visit(&mut self, node: (String, char)) {
            self.indices.insert(node.clone(), self.index);
            self.lowlinks.insert(node.clone(), self.index);
            self.index += 1;
            self.stack.push(node.clone());
            self.on_stack.insert(node.clone());

            for target in self.adjacency.get(&node).into_iter().flatten() {
                if !self.indices.contains_key(target) {
                    self.visit(target.clone());
                    let target_low = self.lowlinks[target];
                    let node_low = self.lowlinks[&node].min(target_low);
                    self.lowlinks.insert(node.clone(), node_low);
                } else if self.on_stack.contains(target) {
                    let target_index = self.indices[target];
                    let node_low = self.lowlinks[&node].min(target_index);
                    self.lowlinks.insert(node.clone(), node_low);
                }
            }

            if self.lowlinks[&node] == self.indices[&node] {
                let mut component = Vec::new();
                while let Some(item) = self.stack.pop() {
                    self.on_stack.remove(&item);
                    let done = item == node;
                    component.push(item);
                    if done {
                        break;
                    }
                }
                self.components.push(component);
            }
        }
    }

    let mut tarjan = Tarjan {
        index: 0,
        stack: Vec::new(),
        on_stack: HashSet::new(),
        indices: HashMap::new(),
        lowlinks: HashMap::new(),
        components: Vec::new(),
        adjacency,
    };
    for node in nodes {
        if !tarjan.indices.contains_key(&node) {
            tarjan.visit(node);
        }
    }
    tarjan.components
}

fn run_gfalook(gfa_path: &Path, svg_path: &Path) -> io::Result<()> {
    let status = Command::new("gfalook")
        .arg("-i")
        .arg(gfa_path)
        .arg("-o")
        .arg(svg_path)
        .arg("-x")
        .arg("2400")
        .arg("-y")
        .arg("1400")
        .arg("-a")
        .arg("8")
        .arg("-t")
        .arg("4")
        .arg("-v")
        .arg("1")
        .status()
        .map_err(|err| io::Error::other(format!("failed to run gfalook: {err}")))?;
    if !status.success() {
        return Err(io::Error::other(format!(
            "gfalook failed for {} with status {status}",
            gfa_path.display()
        )));
    }
    Ok(())
}

fn method_name(method: Method) -> &'static str {
    match method {
        Method::Poasta => "POASTA",
        Method::Abpoa => "abPOA",
    }
}

fn sanitize_label(label: &str) -> String {
    label
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect()
}

fn summary_header() -> String {
    [
        "method",
        "label",
        "gfa_path",
        "svg_path",
        "missing_paths",
        "extra_paths",
        "changed_paths",
        "segments",
        "links",
        "paths",
        "path_steps",
        "segment_bp",
        "directed_sccs_gt1",
        "directed_self_loop_edges",
        "segment_self_loop_links",
        "largest_directed_scc",
        "cyclic_oriented_nodes",
        "bp_weighted_path_depth",
        "segment_white_space_bp_fraction",
        "segment_white_space_bp_total",
        "path_white_space_bp_p95",
        "path_white_space_bp_p99",
        "path_white_space_bp_max",
        "path_white_space_bridges_ge_1000bp",
        "path_depth_median",
        "path_depth_p95",
        "path_depth_max",
        "first_path",
        "first_path_bp",
    ]
    .join("\t")
}
