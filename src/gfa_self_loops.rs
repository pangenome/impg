use serde::Serialize;
use std::cmp::Reverse;
use std::collections::{BTreeMap, HashMap, HashSet};
use std::io;

#[derive(Clone, Debug, Default, Serialize)]
pub struct SelfLoopReport {
    pub direct_self_loop_edges: usize,
    pub direct_self_loop_nodes: usize,
    pub adjacent_same_node_path_steps: usize,
    pub adjacent_same_step_path_steps: usize,
    pub repeated_path_runs: usize,
    pub max_repeat_run_len: usize,
    pub length_buckets: Vec<SelfLoopLengthBucket>,
    pub top_nodes: Vec<SelfLoopNodeReport>,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SelfLoopLengthBucket {
    pub sequence_len: usize,
    pub nodes: usize,
    pub direct_self_loop_edges: usize,
}

#[derive(Clone, Debug, Default, Serialize)]
pub struct SelfLoopNodeReport {
    pub node: String,
    pub order: usize,
    pub sequence_len: usize,
    pub sequence_preview: String,
    pub direct_self_loop_edges: usize,
    pub link_orientations: Vec<String>,
    pub total_path_visits: usize,
    pub distinct_paths: usize,
    pub adjacent_same_node_path_steps: usize,
    pub adjacent_same_step_path_steps: usize,
    pub repeated_path_runs: usize,
    pub max_repeat_run_len: usize,
    pub max_repeat_run_path: String,
}

#[derive(Clone, Debug, Eq, PartialEq, Serialize)]
pub struct NormalizeSelfLoopsStats {
    pub input_direct_self_loop_edges: usize,
    pub output_direct_self_loop_edges: usize,
    pub input_adjacent_same_node_path_steps: usize,
    pub output_adjacent_same_node_path_steps: usize,
    pub input_adjacent_same_step_path_steps: usize,
    pub output_adjacent_same_step_path_steps: usize,
    pub normalized_nodes: usize,
    pub collapsed_runs: usize,
    pub created_segments: usize,
    pub added_links: usize,
    pub removed_self_loop_links: usize,
    pub paths_changed: usize,
    pub path_steps_before: usize,
    pub path_steps_after: usize,
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NormalizeSelfLoopsConfig {
    /// Maximum sequence length eligible for run-node normalization. A value of
    /// 0 means no length cap; validation is still exact path preservation.
    pub max_unit_len: usize,
}

impl Default for NormalizeSelfLoopsConfig {
    fn default() -> Self {
        Self { max_unit_len: 0 }
    }
}

#[derive(Clone, Debug, Eq, PartialEq)]
pub struct NormalizeSelfLoopsResult {
    pub gfa: String,
    pub stats: NormalizeSelfLoopsStats,
}

#[derive(Clone, Debug)]
struct SegmentRecord {
    id: String,
    raw: String,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
struct Step {
    node: String,
    orientation: char,
}

#[derive(Clone, Debug, Eq, PartialEq, Hash, Ord, PartialOrd)]
struct LinkKey {
    from: String,
    from_orientation: char,
    to: String,
    to_orientation: char,
}

#[derive(Clone, Debug)]
struct LinkRecord {
    key: LinkKey,
    raw: String,
}

#[derive(Clone, Debug)]
enum PathRecord {
    P {
        name: String,
        steps: Vec<Step>,
        overlaps: String,
        tags: Vec<String>,
    },
    W {
        fields: Vec<String>,
        steps: Vec<Step>,
    },
}

#[derive(Clone, Debug)]
enum RawRecord {
    Header(String),
    Comment(String),
    Blank(String),
    Other(String),
    Segment(usize),
    Link(usize),
    Path(usize),
}

#[derive(Clone, Debug, Default)]
struct ParsedGfa {
    records: Vec<RawRecord>,
    segments: Vec<SegmentRecord>,
    segment_order: HashMap<String, usize>,
    segment_sequences: HashMap<String, String>,
    links: Vec<LinkRecord>,
    paths: Vec<PathRecord>,
}

#[derive(Clone, Debug, Default)]
struct NodePathAcc {
    total_visits: usize,
    distinct_paths: usize,
    adjacent_same_node_path_steps: usize,
    adjacent_same_step_path_steps: usize,
    repeated_path_runs: usize,
    max_repeat_run_len: usize,
    max_repeat_run_path: String,
}

pub fn diagnose_gfa(gfa_text: &str) -> io::Result<SelfLoopReport> {
    diagnose_gfa_with_top_n(gfa_text, 10)
}

pub fn diagnose_gfa_with_top_n(gfa_text: &str, top_n: usize) -> io::Result<SelfLoopReport> {
    let parsed = parse_gfa(gfa_text)?;
    Ok(diagnose_parsed_gfa(&parsed, top_n))
}

pub fn normalize_repeat_self_loops(
    gfa_text: &str,
    config: &NormalizeSelfLoopsConfig,
) -> io::Result<NormalizeSelfLoopsResult> {
    let before_paths = path_sequence_map(gfa_text)?;
    let parsed = parse_gfa(gfa_text)?;
    let before_report = diagnose_parsed_gfa(&parsed, 10);

    let eligible_nodes = eligible_self_loop_nodes(&parsed, config.max_unit_len);
    let targeted_nodes = targeted_run_nodes(&parsed, &eligible_nodes);

    let mut existing_segment_ids = parsed
        .segments
        .iter()
        .map(|segment| segment.id.clone())
        .collect::<HashSet<_>>();
    let mut next_numeric_segment_id = next_numeric_segment_id(&existing_segment_ids);
    let mut run_segment_ids = BTreeMap::<(String, usize), String>::new();
    let mut new_segments = Vec::<SegmentRecord>::new();
    let mut transformed_paths = Vec::<PathRecord>::with_capacity(parsed.paths.len());
    let mut collapsed_runs = 0usize;
    let mut paths_changed = 0usize;
    let path_steps_before = parsed
        .paths
        .iter()
        .map(|path| path.steps().len())
        .sum::<usize>();

    for path in &parsed.paths {
        let (new_steps, path_collapsed_runs) = transform_path_steps(
            path.steps(),
            &targeted_nodes,
            &parsed.segment_sequences,
            &mut existing_segment_ids,
            &mut next_numeric_segment_id,
            &mut run_segment_ids,
            &mut new_segments,
        )?;
        collapsed_runs = collapsed_runs.saturating_add(path_collapsed_runs);
        if new_steps != path.steps() {
            paths_changed += 1;
        }
        transformed_paths.push(path.with_steps(new_steps));
    }

    let mut retained_links = Vec::new();
    let mut retained_link_keys = HashSet::<LinkKey>::new();
    let mut removed_self_loop_links = 0usize;
    for link in &parsed.links {
        if link.key.from == link.key.to && targeted_nodes.contains(&link.key.from) {
            removed_self_loop_links += 1;
            continue;
        }
        retained_link_keys.insert(link.key.clone());
        retained_links.push(link.clone());
    }

    let mut added_links = Vec::<LinkKey>::new();
    for path in &transformed_paths {
        for pair in path.steps().windows(2) {
            let key = LinkKey {
                from: pair[0].node.clone(),
                from_orientation: pair[0].orientation,
                to: pair[1].node.clone(),
                to_orientation: pair[1].orientation,
            };
            if retained_link_keys.insert(key.clone()) {
                added_links.push(key);
            }
        }
    }

    let normalized = emit_normalized_gfa(
        &parsed,
        &new_segments,
        &retained_links,
        &added_links,
        &transformed_paths,
    );
    let after_paths = path_sequence_map(&normalized)?;
    if before_paths != after_paths {
        return Err(io::Error::other(
            "self-loop normalization changed one or more GFA path spellings",
        ));
    }

    let after_report = diagnose_gfa(&normalized)?;
    let path_steps_after = transformed_paths
        .iter()
        .map(|path| path.steps().len())
        .sum::<usize>();
    let stats = NormalizeSelfLoopsStats {
        input_direct_self_loop_edges: before_report.direct_self_loop_edges,
        output_direct_self_loop_edges: after_report.direct_self_loop_edges,
        input_adjacent_same_node_path_steps: before_report.adjacent_same_node_path_steps,
        output_adjacent_same_node_path_steps: after_report.adjacent_same_node_path_steps,
        input_adjacent_same_step_path_steps: before_report.adjacent_same_step_path_steps,
        output_adjacent_same_step_path_steps: after_report.adjacent_same_step_path_steps,
        normalized_nodes: targeted_nodes.len(),
        collapsed_runs,
        created_segments: new_segments.len(),
        added_links: added_links.len(),
        removed_self_loop_links,
        paths_changed,
        path_steps_before,
        path_steps_after,
    };

    Ok(NormalizeSelfLoopsResult {
        gfa: normalized,
        stats,
    })
}

fn path_sequence_map(gfa_text: &str) -> io::Result<BTreeMap<String, String>> {
    crate::resolution::path_sequences(gfa_text)
        .map(|paths| paths.into_iter().collect::<BTreeMap<_, _>>())
}

fn parse_gfa(gfa_text: &str) -> io::Result<ParsedGfa> {
    let mut parsed = ParsedGfa::default();
    for (line_no, line) in gfa_text.lines().enumerate() {
        if line.is_empty() {
            parsed.records.push(RawRecord::Blank(line.to_string()));
            continue;
        }
        if line.starts_with('#') {
            parsed.records.push(RawRecord::Comment(line.to_string()));
            continue;
        }
        let fields = line.split('\t').map(str::to_string).collect::<Vec<_>>();
        let Some(record_type) = fields.first().map(String::as_str) else {
            parsed.records.push(RawRecord::Other(line.to_string()));
            continue;
        };
        match record_type {
            "H" => parsed.records.push(RawRecord::Header(line.to_string())),
            "S" => {
                if fields.len() < 3 {
                    return Err(invalid_gfa(line_no, "S line has fewer than 3 fields"));
                }
                let id = fields[1].clone();
                if parsed.segment_order.contains_key(&id) {
                    return Err(invalid_gfa(line_no, "duplicate segment ID"));
                }
                let idx = parsed.segments.len();
                parsed.segment_order.insert(id.clone(), idx);
                parsed
                    .segment_sequences
                    .insert(id.clone(), fields[2].clone());
                parsed.segments.push(SegmentRecord {
                    id,
                    raw: line.to_string(),
                });
                parsed.records.push(RawRecord::Segment(idx));
            }
            "L" => {
                if fields.len() < 6 {
                    return Err(invalid_gfa(line_no, "L line has fewer than 6 fields"));
                }
                let key = LinkKey {
                    from: fields[1].clone(),
                    from_orientation: parse_orientation(line_no, &fields[2], "L from")?,
                    to: fields[3].clone(),
                    to_orientation: parse_orientation(line_no, &fields[4], "L to")?,
                };
                let idx = parsed.links.len();
                parsed.links.push(LinkRecord {
                    key,
                    raw: line.to_string(),
                });
                parsed.records.push(RawRecord::Link(idx));
            }
            "P" => {
                if fields.len() < 3 {
                    return Err(invalid_gfa(line_no, "P line has fewer than 3 fields"));
                }
                let steps = parse_p_steps(line_no, &fields[2])?;
                let overlaps = fields.get(3).cloned().unwrap_or_else(|| "*".to_string());
                let tags = fields.get(4..).unwrap_or(&[]).to_vec();
                let idx = parsed.paths.len();
                parsed.paths.push(PathRecord::P {
                    name: fields[1].clone(),
                    steps,
                    overlaps,
                    tags,
                });
                parsed.records.push(RawRecord::Path(idx));
            }
            "W" => {
                if fields.len() < 7 {
                    return Err(invalid_gfa(line_no, "W line has fewer than 7 fields"));
                }
                let steps = parse_w_steps(line_no, &fields[6])?;
                let idx = parsed.paths.len();
                parsed.paths.push(PathRecord::W { fields, steps });
                parsed.records.push(RawRecord::Path(idx));
            }
            _ => parsed.records.push(RawRecord::Other(line.to_string())),
        }
    }
    Ok(parsed)
}

fn invalid_gfa(line_no: usize, msg: &str) -> io::Error {
    io::Error::new(
        io::ErrorKind::InvalidData,
        format!("invalid GFA at line {}: {}", line_no + 1, msg),
    )
}

fn parse_orientation(line_no: usize, raw: &str, field: &str) -> io::Result<char> {
    let mut chars = raw.chars();
    let Some(orientation) = chars.next() else {
        return Err(invalid_gfa(
            line_no,
            &format!("{field} orientation is empty"),
        ));
    };
    if chars.next().is_some() || !matches!(orientation, '+' | '-') {
        return Err(invalid_gfa(
            line_no,
            &format!("{field} orientation must be '+' or '-'"),
        ));
    }
    Ok(orientation)
}

fn parse_p_steps(line_no: usize, raw: &str) -> io::Result<Vec<Step>> {
    if raw == "*" || raw.is_empty() {
        return Ok(Vec::new());
    }
    raw.split(',')
        .filter(|token| !token.is_empty())
        .map(|token| {
            let split = token.len().saturating_sub(1);
            if !token.is_char_boundary(split) {
                return Err(invalid_gfa(line_no, "P line has malformed path step"));
            }
            let (node, orientation) = token.split_at(split);
            let orientation = orientation
                .chars()
                .next()
                .ok_or_else(|| invalid_gfa(line_no, "P line has malformed path step"))?;
            if node.is_empty() || !matches!(orientation, '+' | '-') {
                return Err(invalid_gfa(line_no, "P line has malformed path step"));
            }
            Ok(Step {
                node: node.to_string(),
                orientation,
            })
        })
        .collect()
}

fn parse_w_steps(line_no: usize, raw: &str) -> io::Result<Vec<Step>> {
    if raw == "*" || raw.is_empty() {
        return Ok(Vec::new());
    }
    let mut steps = Vec::new();
    let mut current_orientation = None;
    let mut current_start = 0usize;
    for (idx, ch) in raw.char_indices() {
        if !matches!(ch, '>' | '<') {
            continue;
        }
        if let Some(orientation) = current_orientation {
            if current_start == idx {
                return Err(invalid_gfa(line_no, "W line has empty path step"));
            }
            steps.push(Step {
                node: raw[current_start..idx].to_string(),
                orientation,
            });
        }
        current_orientation = Some(if ch == '>' { '+' } else { '-' });
        current_start = idx + ch.len_utf8();
    }
    let Some(orientation) = current_orientation else {
        return Err(invalid_gfa(line_no, "W line has no oriented path steps"));
    };
    if current_start == raw.len() {
        return Err(invalid_gfa(line_no, "W line has empty final path step"));
    }
    steps.push(Step {
        node: raw[current_start..].to_string(),
        orientation,
    });
    Ok(steps)
}

fn diagnose_parsed_gfa(parsed: &ParsedGfa, top_n: usize) -> SelfLoopReport {
    let mut self_loop_edges_by_node = HashMap::<String, usize>::new();
    let mut orientations_by_node = HashMap::<String, Vec<String>>::new();
    for link in &parsed.links {
        if link.key.from != link.key.to {
            continue;
        }
        *self_loop_edges_by_node
            .entry(link.key.from.clone())
            .or_default() += 1;
        orientations_by_node
            .entry(link.key.from.clone())
            .or_default()
            .push(format!(
                "{}{}",
                link.key.from_orientation, link.key.to_orientation
            ));
    }

    let mut acc_by_node = HashMap::<String, NodePathAcc>::new();
    let mut adjacent_same_node_path_steps = 0usize;
    let mut adjacent_same_step_path_steps = 0usize;
    let mut repeated_path_runs = 0usize;
    let mut max_repeat_run_len = 0usize;

    for path in &parsed.paths {
        let mut seen_in_path = HashSet::<&str>::new();
        for step in path.steps() {
            let acc = acc_by_node.entry(step.node.clone()).or_default();
            acc.total_visits += 1;
            if seen_in_path.insert(step.node.as_str()) {
                acc.distinct_paths += 1;
            }
        }
        for pair in path.steps().windows(2) {
            if pair[0].node == pair[1].node {
                adjacent_same_node_path_steps += 1;
                acc_by_node
                    .entry(pair[0].node.clone())
                    .or_default()
                    .adjacent_same_node_path_steps += 1;
                if pair[0].orientation == pair[1].orientation {
                    adjacent_same_step_path_steps += 1;
                    acc_by_node
                        .entry(pair[0].node.clone())
                        .or_default()
                        .adjacent_same_step_path_steps += 1;
                }
            }
        }
        let steps = path.steps();
        let mut i = 0usize;
        while i < steps.len() {
            let mut j = i + 1;
            while j < steps.len()
                && steps[j].node == steps[i].node
                && steps[j].orientation == steps[i].orientation
            {
                j += 1;
            }
            let run_len = j - i;
            if run_len > 1 {
                repeated_path_runs += 1;
                max_repeat_run_len = max_repeat_run_len.max(run_len);
                let acc = acc_by_node.entry(steps[i].node.clone()).or_default();
                acc.repeated_path_runs += 1;
                if run_len > acc.max_repeat_run_len {
                    acc.max_repeat_run_len = run_len;
                    acc.max_repeat_run_path = path.name().to_string();
                }
            }
            i = j;
        }
    }

    let mut length_counts = BTreeMap::<usize, (usize, usize)>::new();
    for (node, edges) in &self_loop_edges_by_node {
        let len = parsed
            .segment_sequences
            .get(node)
            .map(|seq| sequence_len(seq))
            .unwrap_or(0);
        let entry = length_counts.entry(len).or_default();
        entry.0 += 1;
        entry.1 += *edges;
    }
    let length_buckets = length_counts
        .into_iter()
        .map(
            |(sequence_len, (nodes, direct_self_loop_edges))| SelfLoopLengthBucket {
                sequence_len,
                nodes,
                direct_self_loop_edges,
            },
        )
        .collect::<Vec<_>>();

    let mut top_nodes = self_loop_edges_by_node
        .iter()
        .map(|(node, edges)| {
            let acc = acc_by_node.get(node).cloned().unwrap_or_default();
            let sequence = parsed
                .segment_sequences
                .get(node)
                .map(String::as_str)
                .unwrap_or("");
            let mut link_orientations = orientations_by_node.get(node).cloned().unwrap_or_default();
            link_orientations.sort();
            link_orientations.dedup();
            SelfLoopNodeReport {
                node: node.clone(),
                order: parsed
                    .segment_order
                    .get(node)
                    .copied()
                    .unwrap_or(usize::MAX),
                sequence_len: sequence_len(sequence),
                sequence_preview: sequence_preview(sequence, 48),
                direct_self_loop_edges: *edges,
                link_orientations,
                total_path_visits: acc.total_visits,
                distinct_paths: acc.distinct_paths,
                adjacent_same_node_path_steps: acc.adjacent_same_node_path_steps,
                adjacent_same_step_path_steps: acc.adjacent_same_step_path_steps,
                repeated_path_runs: acc.repeated_path_runs,
                max_repeat_run_len: acc.max_repeat_run_len,
                max_repeat_run_path: acc.max_repeat_run_path,
            }
        })
        .collect::<Vec<_>>();
    top_nodes.sort_unstable_by_key(|node| {
        Reverse((
            node.adjacent_same_step_path_steps,
            node.total_path_visits,
            node.direct_self_loop_edges,
            usize::MAX.saturating_sub(node.order),
        ))
    });
    top_nodes.truncate(top_n);

    SelfLoopReport {
        direct_self_loop_edges: self_loop_edges_by_node.values().sum(),
        direct_self_loop_nodes: self_loop_edges_by_node.len(),
        adjacent_same_node_path_steps,
        adjacent_same_step_path_steps,
        repeated_path_runs,
        max_repeat_run_len,
        length_buckets,
        top_nodes,
    }
}

fn eligible_self_loop_nodes(parsed: &ParsedGfa, max_unit_len: usize) -> HashSet<String> {
    parsed
        .links
        .iter()
        .filter(|link| link.key.from == link.key.to)
        .filter_map(|link| {
            let sequence = parsed.segment_sequences.get(&link.key.from)?;
            let len = sequence_len(sequence);
            (sequence != "*" && len > 0 && (max_unit_len == 0 || len <= max_unit_len))
                .then(|| link.key.from.clone())
        })
        .collect()
}

fn targeted_run_nodes(parsed: &ParsedGfa, eligible_nodes: &HashSet<String>) -> HashSet<String> {
    let mut targeted = HashSet::new();
    for path in &parsed.paths {
        let steps = path.steps();
        let mut i = 0usize;
        while i < steps.len() {
            let mut j = i + 1;
            while j < steps.len()
                && steps[j].node == steps[i].node
                && steps[j].orientation == steps[i].orientation
            {
                j += 1;
            }
            if j - i > 1 && eligible_nodes.contains(&steps[i].node) {
                targeted.insert(steps[i].node.clone());
            }
            i = j;
        }
    }
    targeted
}

fn transform_path_steps(
    steps: &[Step],
    targeted_nodes: &HashSet<String>,
    segment_sequences: &HashMap<String, String>,
    existing_segment_ids: &mut HashSet<String>,
    next_numeric_segment_id: &mut Option<usize>,
    run_segment_ids: &mut BTreeMap<(String, usize), String>,
    new_segments: &mut Vec<SegmentRecord>,
) -> io::Result<(Vec<Step>, usize)> {
    let mut out = Vec::with_capacity(steps.len());
    let mut collapsed_runs = 0usize;
    let mut i = 0usize;
    while i < steps.len() {
        let mut j = i + 1;
        while j < steps.len()
            && steps[j].node == steps[i].node
            && steps[j].orientation == steps[i].orientation
        {
            j += 1;
        }
        let run_len = j - i;
        if run_len > 1 && targeted_nodes.contains(&steps[i].node) {
            let run_node = run_segment_id(
                &steps[i].node,
                run_len,
                segment_sequences,
                existing_segment_ids,
                next_numeric_segment_id,
                run_segment_ids,
                new_segments,
            )?;
            out.push(Step {
                node: run_node,
                orientation: steps[i].orientation,
            });
            collapsed_runs += 1;
        } else {
            out.extend_from_slice(&steps[i..j]);
        }
        i = j;
    }
    Ok((out, collapsed_runs))
}

fn run_segment_id(
    node: &str,
    run_len: usize,
    segment_sequences: &HashMap<String, String>,
    existing_segment_ids: &mut HashSet<String>,
    next_numeric_segment_id: &mut Option<usize>,
    run_segment_ids: &mut BTreeMap<(String, usize), String>,
    new_segments: &mut Vec<SegmentRecord>,
) -> io::Result<String> {
    let key = (node.to_string(), run_len);
    if let Some(id) = run_segment_ids.get(&key) {
        return Ok(id.clone());
    }
    let Some(unit) = segment_sequences.get(node) else {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("path references unknown self-loop segment '{node}'"),
        ));
    };
    let id = unique_run_segment_id(node, run_len, existing_segment_ids, next_numeric_segment_id);
    let sequence = unit.repeat(run_len);
    new_segments.push(SegmentRecord {
        id: id.clone(),
        raw: format!("S\t{id}\t{sequence}"),
    });
    run_segment_ids.insert(key, id.clone());
    Ok(id)
}

fn next_numeric_segment_id(existing: &HashSet<String>) -> Option<usize> {
    let mut max_id = 0usize;
    for id in existing {
        let parsed = id.parse::<usize>().ok()?;
        max_id = max_id.max(parsed);
    }
    Some(max_id.saturating_add(1))
}

fn unique_run_segment_id(
    node: &str,
    run_len: usize,
    existing: &mut HashSet<String>,
    next_numeric: &mut Option<usize>,
) -> String {
    if let Some(mut candidate_num) = *next_numeric {
        while candidate_num < usize::MAX {
            let candidate = candidate_num.to_string();
            candidate_num = candidate_num.saturating_add(1);
            *next_numeric = Some(candidate_num);
            if existing.insert(candidate.clone()) {
                return candidate;
            }
        }
        *next_numeric = None;
    }

    let base = format!("{}_selfloop_run{}", sanitize_segment_id(node), run_len);
    if existing.insert(base.clone()) {
        return base;
    }
    for suffix in 1usize.. {
        let candidate = format!("{base}_{suffix}");
        if existing.insert(candidate.clone()) {
            return candidate;
        }
    }
    unreachable!("unbounded suffix search should always find a segment ID")
}

fn sanitize_segment_id(value: &str) -> String {
    let mut out = String::with_capacity(value.len().max(1));
    for ch in value.chars() {
        if ch.is_ascii_alphanumeric() || matches!(ch, '_' | '.' | ':' | '#') {
            out.push(ch);
        } else {
            out.push('_');
        }
    }
    while out.ends_with(['+', '-']) {
        out.pop();
    }
    if out.is_empty() {
        "node".to_string()
    } else {
        out
    }
}

fn emit_normalized_gfa(
    parsed: &ParsedGfa,
    new_segments: &[SegmentRecord],
    retained_links: &[LinkRecord],
    added_links: &[LinkKey],
    transformed_paths: &[PathRecord],
) -> String {
    let mut out = String::with_capacity(
        parsed
            .records
            .iter()
            .map(|record| match record {
                RawRecord::Header(raw)
                | RawRecord::Comment(raw)
                | RawRecord::Blank(raw)
                | RawRecord::Other(raw) => raw.len() + 1,
                RawRecord::Segment(idx) => parsed.segments[*idx].raw.len() + 1,
                RawRecord::Link(idx) => parsed.links[*idx].raw.len() + 1,
                RawRecord::Path(idx) => parsed.paths[*idx].steps().len() * 16 + 32,
            })
            .sum::<usize>()
            + new_segments
                .iter()
                .map(|segment| segment.raw.len() + 1)
                .sum::<usize>()
            + added_links.len() * 32,
    );

    for record in &parsed.records {
        match record {
            RawRecord::Header(raw) | RawRecord::Comment(raw) | RawRecord::Blank(raw) => {
                push_line(&mut out, raw);
            }
            _ => {}
        }
    }
    for segment in &parsed.segments {
        push_line(&mut out, &segment.raw);
    }
    for segment in new_segments {
        push_line(&mut out, &segment.raw);
    }
    for record in &parsed.records {
        if let RawRecord::Other(raw) = record {
            push_line(&mut out, raw);
        }
    }
    for link in retained_links {
        push_line(&mut out, &link.raw);
    }
    for link in added_links {
        push_line(
            &mut out,
            &format!(
                "L\t{}\t{}\t{}\t{}\t0M",
                link.from, link.from_orientation, link.to, link.to_orientation
            ),
        );
    }
    for path in transformed_paths {
        push_line(&mut out, &path.to_gfa_line());
    }
    out
}

fn push_line(out: &mut String, line: &str) {
    out.push_str(line);
    out.push('\n');
}

impl PathRecord {
    fn name(&self) -> &str {
        match self {
            PathRecord::P { name, .. } => name,
            PathRecord::W { fields, .. } => fields.get(3).map(String::as_str).unwrap_or("*"),
        }
    }

    fn steps(&self) -> &[Step] {
        match self {
            PathRecord::P { steps, .. } | PathRecord::W { steps, .. } => steps,
        }
    }

    fn with_steps(&self, steps: Vec<Step>) -> Self {
        match self {
            PathRecord::P {
                name,
                overlaps,
                tags,
                ..
            } => PathRecord::P {
                name: name.clone(),
                steps,
                overlaps: overlaps.clone(),
                tags: tags.clone(),
            },
            PathRecord::W { fields, .. } => PathRecord::W {
                fields: fields.clone(),
                steps,
            },
        }
    }

    fn to_gfa_line(&self) -> String {
        match self {
            PathRecord::P {
                name,
                steps,
                overlaps,
                tags,
            } => {
                let mut fields = vec![
                    "P".to_string(),
                    name.clone(),
                    steps_to_p_walk(steps),
                    overlaps.clone(),
                ];
                fields.extend(tags.iter().cloned());
                fields.join("\t")
            }
            PathRecord::W { fields, steps } => {
                let mut out_fields = fields.clone();
                if out_fields.len() < 7 {
                    out_fields.resize(7, "*".to_string());
                }
                out_fields[6] = steps_to_w_walk(steps);
                out_fields.join("\t")
            }
        }
    }
}

fn steps_to_p_walk(steps: &[Step]) -> String {
    if steps.is_empty() {
        "*".to_string()
    } else {
        steps
            .iter()
            .map(step_to_p_token)
            .collect::<Vec<_>>()
            .join(",")
    }
}

fn steps_to_w_walk(steps: &[Step]) -> String {
    if steps.is_empty() {
        "*".to_string()
    } else {
        steps.iter().map(step_to_w_token).collect::<String>()
    }
}

fn step_to_p_token(step: &Step) -> String {
    format!("{}{}", step.node, step.orientation)
}

fn step_to_w_token(step: &Step) -> String {
    let orientation = if step.orientation == '+' { '>' } else { '<' };
    format!("{orientation}{}", step.node)
}

fn sequence_len(sequence: &str) -> usize {
    if sequence == "*" {
        0
    } else {
        sequence.len()
    }
}

fn sequence_preview(sequence: &str, max_len: usize) -> String {
    if sequence.len() <= max_len {
        sequence.to_string()
    } else {
        format!("{}...", &sequence[..max_len])
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn path_map(gfa: &str) -> BTreeMap<String, String> {
        crate::resolution::path_sequences(gfa)
            .unwrap()
            .into_iter()
            .collect()
    }

    #[test]
    fn diagnoses_self_loop_edges_and_path_repeats() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tC
L\t1\t+\t1\t+\t0M
L\t1\t+\t2\t+\t0M
P\tp1\t1+,1+,1+,2+\t*
P\tp2\t1+,2+\t*
";
        let report = diagnose_gfa(gfa).unwrap();
        assert_eq!(report.direct_self_loop_edges, 1);
        assert_eq!(report.direct_self_loop_nodes, 1);
        assert_eq!(report.adjacent_same_node_path_steps, 2);
        assert_eq!(report.adjacent_same_step_path_steps, 2);
        assert_eq!(report.repeated_path_runs, 1);
        assert_eq!(report.max_repeat_run_len, 3);
        assert_eq!(report.length_buckets[0].sequence_len, 1);
        assert_eq!(report.top_nodes[0].node, "1");
        assert_eq!(report.top_nodes[0].sequence_preview, "A");
        assert_eq!(report.top_nodes[0].total_path_visits, 4);
    }

    #[test]
    fn normalizes_one_base_self_loop_runs_and_preserves_path_spelling() {
        let gfa = "\
H\tVN:Z:1.0
S\tstart\tGG
S\ta\tA
S\tend\tTT
L\tstart\t+\ta\t+\t0M
L\ta\t+\ta\t+\t0M
L\ta\t+\tend\t+\t0M
P\ttriple\tstart+,a+,a+,a+,end+\t*
P\tsingle\tstart+,a+,end+\t*
";
        let before = path_map(gfa);
        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig::default()).unwrap();
        let after = path_map(&result.gfa);
        assert_eq!(before, after);
        assert_eq!(after["triple"], "GGAAATT");
        assert!(result.gfa.contains("\nS\ta_selfloop_run3\tAAA\n"));
        assert!(!result.gfa.contains("\nL\ta\t+\ta\t+\t0M\n"));
        assert_eq!(result.stats.input_direct_self_loop_edges, 1);
        assert_eq!(result.stats.output_direct_self_loop_edges, 0);
        assert_eq!(result.stats.input_adjacent_same_step_path_steps, 2);
        assert_eq!(result.stats.output_adjacent_same_step_path_steps, 0);
        assert_eq!(result.stats.collapsed_runs, 1);
    }

    #[test]
    fn preserves_p_line_star_overlap_and_tags_when_rewriting_paths() {
        let gfa = "\
H\tVN:Z:1.0
S\ta\tA
S\tb\tT
L\ta\t+\ta\t+\t0M
L\ta\t+\tb\t+\t0M
P\tp\ta+,a+,b+\t*\tTP:Z:kept
";
        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig::default()).unwrap();
        assert_eq!(path_map(gfa), path_map(&result.gfa));
        assert!(
            result
                .gfa
                .contains("\nP\tp\ta_selfloop_run2+,b+\t*\tTP:Z:kept\n"),
            "{}",
            result.gfa
        );
    }

    #[test]
    fn normalizes_two_base_self_loop_runs_and_preserves_reverse_walks() {
        let gfa = "\
H\tVN:Z:1.0
S\tleft\tG
S\tunit\tAT
S\tright\tC
L\tleft\t+\tunit\t+\t0M
L\tunit\t+\tunit\t+\t0M
L\tunit\t-\tunit\t-\t0M
L\tunit\t+\tright\t+\t0M
W\ts\t0\tchr1\t0\t5\t>left>unit>unit>right
W\tt\t0\tchr2\t0\t5\t>left<unit<unit>right
";
        let before = path_map(gfa);
        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig::default()).unwrap();
        let after = path_map(&result.gfa);
        assert_eq!(before, after);
        assert_eq!(after["chr1"], "GATATC");
        assert_eq!(after["chr2"], "GATATC");
        assert!(result.gfa.contains("\nS\tunit_selfloop_run2\tATAT\n"));
        assert_eq!(result.stats.output_adjacent_same_step_path_steps, 0);
        assert_eq!(result.stats.created_segments, 1);
        assert_eq!(result.stats.collapsed_runs, 2);
    }

    #[test]
    fn max_unit_len_limits_normalization_scope() {
        let gfa = "\
H\tVN:Z:1.0
S\tunit\tAT
L\tunit\t+\tunit\t+\t0M
P\tp\tunit+,unit+\t*
";
        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig { max_unit_len: 1 })
                .unwrap();
        assert_eq!(result.stats.normalized_nodes, 0);
        assert_eq!(result.stats.output_direct_self_loop_edges, 1);
        assert_eq!(path_map(gfa), path_map(&result.gfa));
    }

    #[test]
    fn numeric_segment_graphs_get_numeric_run_node_ids() {
        let gfa = "\
H\tVN:Z:1.0
S\t1\tA
S\t2\tT
L\t1\t+\t1\t+\t0M
L\t1\t+\t2\t+\t0M
P\tp\t1+,1+,2+\t*
";
        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig::default()).unwrap();
        assert_eq!(path_map(gfa), path_map(&result.gfa));
        assert!(result.gfa.contains("\nS\t3\tAA\n"), "{}", result.gfa);
        assert!(
            !result.gfa.contains("selfloop_run"),
            "numeric C4-style graphs should remain numeric-ID compatible:\n{}",
            result.gfa
        );
    }

    #[test]
    fn c4_like_direct_loop_repeat_metrics_drop_to_zero_and_paths_compare() {
        let gfa = "\
H\tVN:Z:1.0
S\t7067\tAC
S\t7068\tG
S\t7069\tTT
L\t7067\t+\t7068\t+\t0M
L\t7068\t+\t7068\t+\t0M
L\t7068\t+\t7069\t+\t0M
P\tHG01358#2#CM089068.1:31846088-32072428(+)\t7067+,7068+,7068+,7068+,7068+,7069+\t*
P\tHG02055#1#CM088166.1:31846088-32072428(+)\t7067+,7068+,7068+,7069+\t*
P\tGRCh38#0#chr6:31846088-32072428(+)\t7067+,7068+,7069+\t*
";
        let before_report = diagnose_gfa(gfa).unwrap();
        assert_eq!(before_report.direct_self_loop_edges, 1);
        assert_eq!(before_report.repeated_path_runs, 2);
        assert_eq!(before_report.adjacent_same_step_path_steps, 4);

        let result =
            normalize_repeat_self_loops(gfa, &NormalizeSelfLoopsConfig::default()).unwrap();
        let after_report = diagnose_gfa(&result.gfa).unwrap();

        assert_eq!(path_map(gfa), path_map(&result.gfa));
        assert_eq!(result.stats.input_direct_self_loop_edges, 1);
        assert_eq!(result.stats.output_direct_self_loop_edges, 0);
        assert_eq!(after_report.direct_self_loop_edges, 0);
        assert_eq!(after_report.repeated_path_runs, 0);
        assert_eq!(after_report.adjacent_same_step_path_steps, 0);
        assert!(
            path_map(gfa) == path_map(&result.gfa),
            "compare_gfa_paths-style path spelling check must pass"
        );
    }
}
