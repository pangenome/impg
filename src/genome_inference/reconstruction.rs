//! Conservative source-supported blocks, never a chromosome assembly or truth-guided splice.
use super::{catalog::*, genotype::*, threading::*, *};
use crate::{
    sample_mem_bwt::invalid,
    sequence_index::{SequenceIndex, UnifiedSequenceIndex},
};
use serde_json::{json, Value};
use std::{
    collections::{BTreeMap, BTreeSet},
    fs::File,
    io::{self, Read},
    path::{Path, PathBuf},
};

pub const MODEL: &str = "source-supported-block-reconstruction-v1";

/// Streaming content fingerprint; FNV is provenance/corruption detection, not authentication.
pub fn fingerprint(path: &Path) -> io::Result<Value> {
    let mut file = File::open(path)?;
    let mut hash = 0xcbf29ce484222325;
    let mut size = 0u64;
    let mut buffer = vec![0; 1024 * 1024];
    loop {
        let n = file.read(&mut buffer)?;
        if n == 0 {
            break;
        }
        super::hash_update(&mut hash, &buffer[..n]);
        size += n as u64;
    }
    Ok(json!({"path":path,"bytes":size,"fnv1a64":format!("{hash:016x}")}))
}
fn same_interval(a: &Interval, b: &Interval) -> bool {
    a.path == b.path && a.start == b.start && a.end == b.end
}
fn digest(s: &Option<String>) -> bool {
    s.as_ref()
        .is_some_and(|s| s.len() == 16 && s.bytes().all(|b| b.is_ascii_hexdigit()))
}

/// Recreate only small source/group metadata, NOT the feature catalog. Re-run the
/// original DP and block construction to check all selected marginals and blocks.
/// Signed-context orientations are frozen evidence from threading (not re-inferred).
pub fn validate(
    calls: &Genotypes,
    threads: &Threads,
    lengths: &BTreeMap<String, u64>,
) -> io::Result<()> {
    calls.parameters.validate()?;
    let new_model = calls.model == super::observations::MODEL;
    let linked = if new_model {
        calls.observations.as_ref().is_some_and(|p| {
            p.full_feature_scope
                && p.feature_groups.is_empty()
                && p.orientation_support
                    .windows(2)
                    .all(|w| (w[0].0, w[0].1) < (w[1].0, w[1].1))
                && digest(&Some(p.metadata_checksum.clone()))
                && threads.observation_metadata_checksum.as_ref() == Some(&p.metadata_checksum)
        }) && calls.catalog_payload_checksum.is_none()
            && calls.artifact_checksum_algorithm == super::observations::CHECKSUM_ALGORITHM
    } else {
        calls.observations.is_none()
            && threads.observation_metadata_checksum.is_none()
            && digest(&calls.catalog_payload_checksum)
            && calls.artifact_checksum_algorithm
                == "fnv1a64-payload-v1; sample=bincode-payload; catalog=compact-json-payload"
    };
    if calls.version != FORMAT_VERSION
        || (!new_model && calls.model != super::genotype::MODEL)
        || calls.ploidy != 1
        || calls.count_policy != COUNT_POLICY
        || !linked
        || !calls.experimental
        || calls.catalog_accepted
        || !digest(&calls.sample_payload_checksum)
        || threads.panel != calls.panel
        || threads.sample_payload_checksum != calls.sample_payload_checksum
        || threads.catalog_payload_checksum != calls.catalog_payload_checksum
        || threads.count_policy != calls.count_policy
    {
        return Err(invalid(
            "incompatible calls/threads version, model, panel or payload linkage",
        ));
    }
    let mut catalog = Catalog {
        version: FORMAT_VERSION,
        panel: calls.panel.clone(),
        catalog_accepted: false,
        validation: "reconstruction metadata only".into(),
        sources: Vec::new(),
        occurrences: Vec::new(),
        groups: Vec::new(),
        links: Vec::new(),
        features: Vec::new(),
    };
    let mut groups = BTreeMap::new();
    let mut paths = BTreeMap::new();
    for c in &calls.calls {
        if c.group.is_empty()
            || groups
                .insert(c.group.clone(), catalog.groups.len())
                .is_some()
            || c.bundles.is_empty()
        {
            return Err(invalid("duplicate/empty call group or bundles"));
        }
        catalog.groups.push(Group {
            id: c.group.clone(),
            scaffold: None,
            occurrences: Vec::new(),
        });
    }
    for (id, o) in calls.source_occurrences.iter().enumerate() {
        let &g = groups
            .get(&o.group)
            .ok_or_else(|| invalid("foreign occurrence group"))?;
        let &length = lengths
            .get(&o.interval.path)
            .ok_or_else(|| invalid("source absent from panel names"))?;
        if o.id != id
            || o.identity != source_identity(&o.interval.path)?
            || o.interval.start >= o.interval.end
            || o.interval.end > length
            || o.interval.end > i32::MAX as u64
            || o.interval
                .strand
                .as_deref()
                .is_some_and(|s| s != "+" && s != "-")
        {
            return Err(invalid(
                "invalid source ID, identity, coordinates or strand",
            ));
        }
        let source = *paths.entry(o.interval.path.clone()).or_insert_with(|| {
            let id = catalog.sources.len();
            catalog.sources.push(Source {
                id,
                path: o.interval.path.clone(),
                length,
            });
            id
        });
        catalog.groups[g].occurrences.push(id);
        catalog.occurrences.push(Occurrence {
            id,
            group: g,
            source,
            interval: o.interval.clone(),
            fully_contained_anchors: 0,
            anchor_status: "frozen-thread-orientations".into(),
            feature_multiplicities: BTreeMap::new(),
        });
    }
    for (g, c) in calls.calls.iter().enumerate() {
        let mut members = BTreeSet::new();
        let mut identities = BTreeSet::new();
        for b in &c.bundles {
            if !identities.insert(&b.identity)
                || b.source_occurrences.is_empty()
                || !b.score.is_finite()
                || !b.mean_deviance.is_finite()
                || b.mean_deviance < 0.0
            {
                return Err(invalid("invalid genotype bundle"));
            }
            for &id in &b.source_occurrences {
                let o = calls
                    .source_occurrences
                    .get(id)
                    .ok_or_else(|| invalid("unknown bundle occurrence"))?;
                if o.group != c.group || o.identity != b.identity || !members.insert(id) {
                    return Err(invalid("bundle membership mismatch"));
                }
            }
        }
        if members != catalog.groups[g].occurrences.iter().copied().collect() {
            return Err(invalid("incomplete bundle membership"));
        }
        let best = c
            .bundles
            .iter()
            .map(|b| b.score)
            .fold(f64::INFINITY, f64::min);
        let expected: Vec<_> = c
            .bundles
            .iter()
            .enumerate()
            .filter(|(_, b)| b.score - best <= TIE_EPSILON)
            .map(|(i, _)| i)
            .collect();
        let mut reported = c.best_bundles.clone();
        reported.sort_unstable();
        let status = if c.retained_features == 0 {
            "no-call-no-local-features"
        } else if c.observed_total == 0.0 {
            "no-call-no-positive-features"
        } else if c
            .bundles
            .iter()
            .min_by(|a, b| a.score.total_cmp(&b.score))
            .unwrap()
            .mean_deviance
            > calls.parameters.max_mean_deviance
        {
            "poor-fit"
        } else if expected.len() > 1 {
            "tied"
        } else {
            "informative"
        };
        let (status, mut expected) = if new_model {
            let (status, expected, _) = super::observations::classify(c, &calls.parameters);
            (status, expected)
        } else {
            (status.to_string(), expected)
        };
        // Classifiers may order epsilon ties by score rather than bundle index.
        // Normalize order only: duplicate reported indices must still fail.
        expected.sort_unstable();
        if c.status != status
            || (if status.starts_with("no-call") {
                !reported.is_empty()
            } else {
                reported != expected
            })
        {
            return Err(invalid("inconsistent local call status/best bundles"));
        }
    }
    let mut axis = Axis {
        version: FORMAT_VERSION,
        coordinate_system: threads.coordinate_system.clone(),
        intervals: Vec::new(),
    };
    for row in &threads.intervals {
        let mut a = row.axis.clone();
        for s in &row.states {
            if s.source_occurrences.is_empty() {
                return Err(invalid("empty physical state"));
            }
            for &id in &s.source_occurrences {
                let o = calls
                    .source_occurrences
                    .get(id)
                    .ok_or_else(|| invalid("foreign state occurrence"))?;
                if !same_interval(&o.interval, &s.interval)
                    || (id == s.source_occurrences[0] && o.interval.strand != s.interval.strand)
                    || o.group != a.group
                    || o.identity != s.identity
                {
                    return Err(invalid("state/source provenance mismatch"));
                }
                let reference = calls
                    .source_occurrences
                    .get(a.reference_occurrence)
                    .ok_or_else(|| invalid("unknown reference occurrence"))?;
                let physical_reference = o.interval.path == reference.interval.path
                    && o.interval.start == reference.interval.start
                    && o.interval.end == reference.interval.end;
                let established = if physical_reference {
                    Some((a.reference_strand.clone(), "reference-identity"))
                } else if let Some(strand) = row.axis.orientations.get(&id) {
                    Some((strand.clone(), "explicit-axis-input"))
                } else if let (Some(strand), Some(r)) =
                    (&o.interval.strand, &reference.interval.strand)
                {
                    Some((
                        if strand == r {
                            a.reference_strand.clone()
                        } else if a.reference_strand == "+" {
                            "-".into()
                        } else {
                            "+".into()
                        },
                        "explicit-catalog-relative-to-reference",
                    ))
                } else {
                    None
                };
                if established.is_none() {
                    if let Some(p) = &calls.observations {
                        let mask = p
                            .orientation_support
                            .binary_search_by_key(&(a.reference_occurrence, id), |&(r, m, _)| {
                                (r, m)
                            })
                            .ok()
                            .map(|i| p.orientation_support[i].2)
                            .unwrap_or(0);
                        let (strand, evidence) = match mask {
                            1 => (
                                Some(a.reference_strand.clone()),
                                "consistent-signed-contexts",
                            ),
                            2 => (
                                Some(if a.reference_strand == "+" { "-" } else { "+" }.into()),
                                "consistent-signed-contexts",
                            ),
                            3 => (None, "conflicting-signed-contexts"),
                            _ => (None, "unsupported-or-palindromic-contexts"),
                        };
                        if s.strand != strand
                            || (id == s.source_occurrences[0] && s.orientation_evidence != evidence)
                        {
                            return Err(invalid(
                                "state contradicts frozen observation orientation evidence",
                            ));
                        }
                    }
                }
                if let Some((strand, evidence)) = established {
                    if s.strand.as_ref() != Some(&strand)
                        || (id == s.source_occurrences[0] && s.orientation_evidence != evidence)
                    {
                        return Err(invalid(
                            "state contradicts established orientation evidence",
                        ));
                    }
                } else if id == s.source_occurrences[0]
                    && ((s.strand.is_some()
                        && s.orientation_evidence != "consistent-signed-contexts")
                        || (s.strand.is_none()
                            && !matches!(
                                s.orientation_evidence.as_str(),
                                "conflicting-signed-contexts"
                                    | "unsupported-or-palindromic-contexts"
                            )))
                {
                    return Err(invalid("state orientation/evidence mismatch"));
                }
                // Do not let an override contradict independently checkable orientation evidence.
                if let Some(strand) = &s.strand {
                    if strand != "+" && strand != "-" {
                        return Err(invalid("unknown state strand spelling"));
                    }
                    if a.orientations.get(&id).is_some_and(|old| old != strand) {
                        return Err(invalid("contradictory axis orientation"));
                    }
                    a.orientations.insert(id, strand.clone());
                }
            }
        }
        axis.intervals.push(a);
    }
    let mut expected = super::threading::thread(&catalog, calls, axis, threads.switch_penalty)?;
    if expected.intervals.len() != threads.intervals.len() {
        return Err(invalid("thread row mismatch"));
    }
    for (e, r) in expected.intervals.iter_mut().zip(&threads.intervals) {
        if e.states.len() != r.states.len() {
            return Err(invalid("thread candidate set mismatch"));
        }
        e.axis = r.axis.clone();
        for (a, b) in e.states.iter_mut().zip(&r.states) {
            if !matches!(
                b.orientation_evidence.as_str(),
                "reference-identity"
                    | "explicit-axis-input"
                    | "explicit-catalog-relative-to-reference"
                    | "consistent-signed-contexts"
                    | "conflicting-signed-contexts"
                    | "unsupported-or-palindromic-contexts"
            ) {
                return Err(invalid("unknown orientation evidence"));
            }
            a.orientation_evidence = b.orientation_evidence.clone();
        }
    }
    // JSON round trips can change a last bit before the summed DP score is
    // recomputed. Only the descriptive segment score gets numeric tolerance;
    // selected marginal sets, statuses, source states and blocks stay exact.
    for (a, b) in expected.segments.iter_mut().zip(&threads.segments) {
        if (a.score - b.score).abs() <= TIE_EPSILON + a.score.abs() * 1e-12 {
            a.score = b.score;
        }
    }
    if serde_json::to_value(expected).map_err(io::Error::other)?
        != serde_json::to_value(threads).map_err(io::Error::other)?
    {
        return Err(invalid("threads disagree with calls or replayed DP/blocks"));
    }
    Ok(())
}

pub struct Reconstruction {
    pub fasta: Vec<u8>,
    pub provenance: Value,
    pub unresolved: Value,
    pub bed: Vec<u8>,
}
fn extract(index: &impl SequenceIndex, interval: &Interval, strand: &str) -> io::Result<Vec<u8>> {
    if interval.start >= interval.end
        || interval.end > i32::MAX as u64
        || !matches!(strand, "+" | "-")
    {
        return Err(invalid("invalid extraction bounds/strand"));
    }
    let mut seq =
        index.fetch_sequence(&interval.path, interval.start as i32, interval.end as i32)?;
    seq.make_ascii_uppercase();
    if seq.len() as u64 != interval.end - interval.start
        || seq.iter().any(|b| !b"ACGTRYSWKMBDHVN".contains(b))
    {
        return Err(invalid("source extraction length/alphabet mismatch"));
    }
    if strand == "-" {
        seq = super::sequence_evaluation::reverse_complement(&seq);
    }
    Ok(seq)
}
fn union(mut spans: Vec<(u64, u64)>) -> u64 {
    spans.sort_unstable();
    let mut end = 0;
    let mut n = 0;
    for (s, e) in spans {
        n += e.saturating_sub(s.max(end));
        end = end.max(e);
    }
    n
}

/// Caller must validate artifacts and source lengths first. Pure sequence spelling,
/// intentionally has no truth argument or access to an alignment score.
pub fn reconstruct(
    calls: &Genotypes,
    threads: &Threads,
    index: &impl SequenceIndex,
    copy_source: bool,
) -> io::Result<Reconstruction> {
    let mut records = Vec::new();
    let mut fasta = Vec::new();
    let mut emitted = BTreeSet::new();
    let mut emitted_bp = 0u64;
    let mut supported_bp = 0u64;
    let mut bridged_bp = 0u64;
    let mut suppressed_bp = 0u64;
    let mut add = |rows: &[usize],
                   alternatives: Vec<&State>,
                   equivalent: bool,
                   dp_block: Option<usize>,
                   previous: Option<&State>|
     -> io::Result<()> {
        let first = alternatives[0];
        let spans: Vec<_> = rows
            .iter()
            .map(|&i| {
                let r = &threads.intervals[i];
                let s = &r.states[r.optimal_states[0]];
                (s.interval.start, s.interval.end)
            })
            .collect();
        let mut interval = if equivalent {
            first.interval.clone()
        } else {
            Interval {
                path: first.interval.path.clone(),
                start: spans.iter().map(|x| x.0).min().unwrap(),
                end: spans.iter().map(|x| x.1).max().unwrap(),
                strand: first.strand.clone(),
            }
        };
        let original_interval = interval.clone();
        if let Some(previous) = previous {
            if first.strand.as_deref() == Some("+") {
                interval.start = interval.start.max(previous.interval.end);
            } else {
                interval.end = interval.end.min(previous.interval.start);
            }
        }
        let suppressed =
            (original_interval.end - original_interval.start) - (interval.end - interval.start);
        let sequence = extract(index, &interval, first.strand.as_deref().unwrap())?;
        let directly_supported = if equivalent {
            sequence.len() as u64
        } else {
            union(
                spans
                    .into_iter()
                    .map(|(s, e)| (s.max(interval.start), e.min(interval.end)))
                    .filter(|(s, e)| s < e)
                    .collect(),
            )
        };
        let bridge = sequence.len() as u64 - directly_supported;
        let id = format!("block{:06}", records.len() + 1);
        fasta.extend(format!(">{id}\n").bytes());
        for chunk in sequence.chunks(80) {
            fasta.extend(chunk);
            fasta.push(b'\n');
        }
        records.push(json!({"id":id,"kind":if equivalent {"sequence-equivalent-unphased-fragment"} else {"selected-source-block"},
            "sequence_fnv1a64":format!("{:016x}",checksum(&sequence)),"emitted_bp":sequence.len(),
            "scored_source_union_bp":directly_supported,"source_bridged_imputed_bp":bridge,"dp_block":dp_block,
            "original_source_span":if equivalent {Value::Null} else {json!(original_interval)},"suppressed_within_block_overlap_bp":suppressed,
            "trim_scope":"Source-coordinate trimming only; original axis memberships retained, NOT an exact source-to-axis base mapping.",
            "possible_source_extractions":if equivalent {alternatives.iter().map(|s|json!({"path":s.interval.path,"start":s.interval.start,"end":s.interval.end,"strand":s.strand})).collect::<Vec<_>>()} else {vec![json!({"path":interval.path,"start":interval.start,"end":interval.end,"strand":first.strand})]}, 
            "extraction":if equivalent {Value::Null} else {json!({"path":interval.path,"start":interval.start,"end":interval.end,"strand":first.strand})},
            "axis_component":threads.intervals[rows[0]].axis.component,
            "axis_span_start":rows.iter().map(|&i|threads.intervals[i].axis.start).min(),
            "axis_span_end":rows.iter().map(|&i|threads.intervals[i].axis.end).max(),
            "scored_axis_union_bp":union(rows.iter().map(|&i|{let a=&threads.intervals[i].axis;(a.start,a.end)}).collect()),
            "interval_indices":rows,"memberships":rows.iter().map(|&i|json!({"interval_index":i,"axis":threads.intervals[i].axis})).collect::<Vec<_>>(),
            "source_states":alternatives}));
        emitted.extend(rows.iter().copied());
        emitted_bp += sequence.len() as u64;
        supported_bp += directly_supported;
        bridged_bp += bridge;
        suppressed_bp += suppressed;
        Ok(())
    };
    for (block_id, block) in threads.blocks.iter().enumerate() {
        let mut start = 0;
        for end in 1..=block.interval_indices.len() {
            if end == block.interval_indices.len()
                || (!copy_source
                    && (block.source_gaps_bp[end - 1] > 0 || block.reference_gaps_bp[end - 1] > 0))
            {
                let rows = &block.interval_indices[start..end];
                let states = rows
                    .iter()
                    .map(|&i| {
                        let r = &threads.intervals[i];
                        &r.states[r.optimal_states[0]]
                    })
                    .collect();
                let previous = if start > 0 {
                    let r = &threads.intervals[block.interval_indices[start - 1]];
                    Some(&r.states[r.optimal_states[0]])
                } else {
                    None
                };
                add(rows, states, false, Some(block_id), previous)?;
                start = end;
            }
        }
    }
    let mut tie_reasons = BTreeMap::new();
    for (i, r) in threads
        .intervals
        .iter()
        .enumerate()
        .filter(|(_, r)| r.status == "unresolved-tied-paths")
    {
        let states: Vec<_> = r.optimal_states.iter().map(|&s| &r.states[s]).collect();
        let call = calls
            .calls
            .iter()
            .find(|c| c.group == r.axis.group)
            .unwrap();
        // A simultaneous multicopy bundle is not a set of interchangeable alleles.
        let multicopy = states.iter().any(|s| {
            call.bundles[s.bundle]
                .source_occurrences
                .iter()
                .map(|&id| {
                    let o = &calls.source_occurrences[id];
                    (&o.interval.path, o.interval.start, o.interval.end)
                })
                .collect::<BTreeSet<_>>()
                .len()
                > 1
        });
        let mut reason = "nonidentical-sequence-alternatives";
        if multicopy {
            reason = "multicopy-placement-tie-not-single-allele";
        } else if states
            .iter()
            .any(|s| s.strand.is_none() || s.mean_deviance > calls.parameters.max_mean_deviance)
        {
            reason = "unknown-orientation-or-poor-fit-alternative";
        } else {
            let strings = states
                .iter()
                .map(|s| extract(index, &s.interval, s.strand.as_deref().unwrap()))
                .collect::<io::Result<Vec<_>>>()?;
            if strings.windows(2).all(|p| p[0] == p[1]) {
                add(&[i], states, true, None, None)?;
                continue;
            }
        }
        tie_reasons.insert(i, reason);
    }
    drop(add);
    let mut reused_source = Vec::new();
    for (i, a) in records.iter().enumerate() {
        for b in &records[i + 1..] {
            if !a["dp_block"].is_null() && a["dp_block"] == b["dp_block"] {
                continue;
            }
            for x in a["possible_source_extractions"].as_array().unwrap() {
                for y in b["possible_source_extractions"].as_array().unwrap() {
                    let start = x["start"]
                        .as_u64()
                        .unwrap()
                        .max(y["start"].as_u64().unwrap());
                    let end = x["end"].as_u64().unwrap().min(y["end"].as_u64().unwrap());
                    if x["path"] == y["path"] && start < end {
                        reused_source.push(json!({"left_block":a["id"],"right_block":b["id"],"path":x["path"],"start":start,"end":end,"overlap_bp":end-start,
                            "status":"copy-or-placement-ambiguity-not-certified-independent-copies","policy":"Both literal block hypotheses retained in full output/evaluation denominators; no global donor-template deduplication."}));
                    }
                }
            }
        }
    }
    let mut unresolved = Vec::new();
    let mut bed = Vec::new();
    for (i, r) in threads.intervals.iter().enumerate() {
        if !emitted.contains(&i) {
            let reason = tie_reasons.get(&i).copied().unwrap_or(&r.status);
            unresolved.push(json!({"interval_index":i,"axis":r.axis,"reason":reason,"optimal_states":r.optimal_states}));
            bed.extend(
                format!(
                    "{}\t{}\t{}\t{}\t{}\n",
                    r.axis.component, r.axis.start, r.axis.end, r.axis.group, reason
                )
                .bytes(),
            );
        }
    }
    let mut gaps = Vec::new();
    for block in &threads.blocks {
        for (j, pair) in block.interval_indices.windows(2).enumerate() {
            let left = &threads.intervals[pair[0]];
            let right = &threads.intervals[pair[1]];
            let a = &left.states[left.optimal_states[0]];
            let b = &right.states[right.optimal_states[0]];
            if block.source_gaps_bp[j] > 0 || block.reference_gaps_bp[j] > 0 {
                gaps.push(json!({"left_interval":pair[0],"right_interval":pair[1],"path":block.source_path,"strand":block.strand,
                    "source_gap_bp":block.source_gaps_bp[j],"reference_gap_bp":block.reference_gaps_bp[j],
                    "source_gap_start":if block.source_gaps_bp[j]>0 {Some(if block.strand=="+" {a.interval.end} else {b.interval.end})} else {None},
                    "source_gap_end":if block.source_gaps_bp[j]>0 {Some(if block.strand=="+" {b.interval.start} else {a.interval.start})} else {None},
                    "source_overlap_start":if block.source_gaps_bp[j]<0 {Some(a.interval.start.max(b.interval.start))} else {None},
                    "source_overlap_end":if block.source_gaps_bp[j]<0 {Some(a.interval.end.min(b.interval.end))} else {None},
                    "policy":if copy_source {"same-donor-imputation-not-genotype-evidence"} else {"split-unemitted"}}));
            }
        }
    }
    for pair in threads.intervals.windows(2) {
        let a = &pair[0].axis;
        let b = &pair[1].axis;
        if a.component == b.component && b.start > a.end {
            unresolved.push(json!({"axis_component":a.component,"start":a.end,"end":b.start,"reason":"unscored-axis-gap-even-if-source-bridged"}));
            bed.extend(
                format!(
                    "{}\t{}\t{}\t.\tunscored-axis-gap-even-if-source-bridged\n",
                    a.component, a.end, b.start
                )
                .bytes(),
            );
        }
    }
    let unresolved_rows = threads.intervals.len() - emitted.len();
    let groups:Vec<_>=calls.calls.iter().map(|c|json!({"group":c.group,"genotype_status":c.status,"best_bundles":c.best_bundles,
        "bundles":c.bundles.iter().enumerate().map(|(b,x)|json!({"bundle":b,"identity":x.identity,"source_occurrences":x.source_occurrences})).collect::<Vec<_>>(),
        "emitted_interval_indices":emitted.iter().filter(|&&i|threads.intervals[i].axis.group==c.group).collect::<Vec<_>>(),
        "unscaffolded":threads.unscaffolded_groups.contains(&c.group),"repeated_axis":threads.repeated_axis_groups.contains(&c.group),
        "scope":"Only listed source_states were spelled; other simultaneous copies remain genotype-only, not reconstructed. Best bundles are local calls, not necessarily DP selections."})).collect();
    Ok(Reconstruction {
        fasta,
        provenance: json!({"version":FORMAT_VERSION,"model":MODEL,"gap_policy":if copy_source {"copy-source"} else {"split"},
        "scope":"Source-supported experimental blocks, not phased chromosomes or complete copy recovery. Scored spans copy donor sequence, not independently base-called alleles. Source bridges are imputed, not observed genotypes. Equivalent ties are independent unphased fragments, never concatenated marginals. Whole truth remains evaluation denominator.",
        "emitted_bp":emitted_bp,"scored_source_union_bp_sum":supported_bp,"source_bridged_imputed_bp":bridged_bp,"emitted_intervals":emitted.len(),
        "unresolved_intervals":unresolved_rows,"blocks":records,"gap_ledger":gaps,"group_ledger":groups,"suppressed_within_block_overlap_bp":suppressed_bp,"cross_block_reused_source_overlaps":reused_source}),
        unresolved: json!({"version":FORMAT_VERSION,"coordinate_system":threads.coordinate_system,"intervals":unresolved,
            "unscaffolded_groups":threads.unscaffolded_groups,"repeated_axis_groups":threads.repeated_axis_groups}),
        bed,
    })
}

pub fn run(
    calls_path: &Path,
    threads_path: &Path,
    names_path: &Path,
    sources: &[String],
    out: &Path,
    copy_source: bool,
) -> io::Result<Value> {
    let calls: Genotypes = read_json(calls_path)?;
    let threads: Threads = read_json(threads_path)?;
    let names_fp = fingerprint(names_path)?;
    if calls.panel.checksum_algorithm != "fnv1a64-content-v1"
        || !calls.panel.sidecars.iter().any(|(s, n, h)| {
            s == "names" && json!(n) == names_fp["bytes"] && json!(h) == names_fp["fnv1a64"]
        })
    {
        return Err(invalid("panel names fingerprint mismatch"));
    }
    let names = crate::syng::SyngNameMap::load(
        names_path
            .to_str()
            .ok_or_else(|| invalid("non-UTF8 names path"))?,
    )?;
    let lengths: BTreeMap<_, _> = names
        .path_to_name
        .into_iter()
        .zip(names.path_to_length)
        .collect();
    validate(&calls, &threads, &lengths)?;
    // Separate indexes disallow ambiguous same-name resolution across source files.
    enum SourceFile {
        Fasta(super::sequence_evaluation::Fasta),
        Agc(UnifiedSequenceIndex),
    }
    impl SequenceIndex for SourceFile {
        fn fetch_sequence(&self, n: &str, s: i32, e: i32) -> io::Result<Vec<u8>> {
            match self {
                Self::Agc(x) => x.fetch_sequence(n, s, e),
                Self::Fasta(f) => f
                    .get(n)
                    .and_then(|b| b.get(s as usize..e as usize))
                    .map(|b| b.to_vec())
                    .ok_or_else(|| invalid("FASTA extraction bounds/name")),
            }
        }
        fn get_sequence_length(&self, n: &str) -> io::Result<usize> {
            match self {
                Self::Agc(x) => x.get_sequence_length(n),
                Self::Fasta(f) => f
                    .get(n)
                    .map(Vec::len)
                    .ok_or_else(|| invalid("unknown FASTA source")),
            }
        }
    }
    let indexes = sources
        .iter()
        .map(|s| {
            if s.ends_with(".agc") {
                UnifiedSequenceIndex::from_files(&[s.clone()]).map(SourceFile::Agc)
            } else {
                super::sequence_evaluation::read_fasta(Path::new(s)).map(SourceFile::Fasta)
            }
        })
        .collect::<io::Result<Vec<_>>>()?;
    let mut chosen = BTreeMap::new();
    for o in &calls.source_occurrences {
        if chosen.contains_key(&o.interval.path) {
            continue;
        }
        let available: Vec<_> = indexes
            .iter()
            .enumerate()
            .filter_map(|(i, x)| x.get_sequence_length(&o.interval.path).ok().map(|n| (i, n)))
            .collect();
        if available.len() != 1 || available[0].1 as u64 != lengths[&o.interval.path] {
            return Err(invalid(format!(
                "source missing, ambiguous or wrong length: {}",
                o.interval.path
            )));
        }
        chosen.insert(o.interval.path.clone(), available[0].0);
    }
    struct Sources {
        indexes: Vec<SourceFile>,
        chosen: BTreeMap<String, usize>,
    }
    impl SequenceIndex for Sources {
        fn fetch_sequence(&self, n: &str, s: i32, e: i32) -> io::Result<Vec<u8>> {
            self.indexes[self.chosen[n]].fetch_sequence(n, s, e)
        }
        fn get_sequence_length(&self, n: &str) -> io::Result<usize> {
            self.indexes[self.chosen[n]].get_sequence_length(n)
        }
    }
    let source_files = sources
        .iter()
        .map(|p| fingerprint(&PathBuf::from(p)))
        .collect::<io::Result<Vec<_>>>()?;
    if let Some(p) = &calls.observations {
        let identities = |files: &[Value]| -> BTreeSet<(String, u64)> {
            files
                .iter()
                .map(|f| {
                    (
                        f["fnv1a64"].as_str().unwrap_or("").to_string(),
                        f["bytes"].as_u64().unwrap_or(0),
                    )
                })
                .collect()
        };
        if identities(&source_files) != identities(&p.source_files) {
            return Err(invalid(
                "reconstruction sources differ from compiled observation context",
            ));
        }
    }
    let source_bindings: Vec<_> = chosen
        .iter()
        .map(|(n, &i)| json!({"path":n,"length":lengths[n],"source_file_index":i}))
        .collect();
    let mut result = reconstruct(&calls, &threads, &Sources { indexes, chosen }, copy_source)?;
    result.provenance["inputs"] = json!({"calls":fingerprint(calls_path)?,"threads":fingerprint(threads_path)?,"panel_names":names_fp,"source_files":source_files,
        "panel":calls.panel,"sample_payload_checksum":calls.sample_payload_checksum,"catalog_payload_checksum":calls.catalog_payload_checksum,"source_bindings":source_bindings});
    if let Some(p) = &calls.observations {
        result.provenance["inputs"]["genotype_model"] = json!(calls.model);
        result.provenance["inputs"]["observations"] = json!(p);
    }
    atomic_write(&out.join("reconstruction.fa"), &result.fasta)?;
    write_json(&out.join("provenance.json"), &result.provenance)?;
    write_json(&out.join("unresolved.json"), &result.unresolved)?;
    atomic_write(&out.join("unresolved.bed"), &result.bed)?;
    Ok(
        json!({"emitted_bp":result.provenance["emitted_bp"],"blocks":result.provenance["blocks"].as_array().unwrap().len(),"unresolved_intervals":result.provenance["unresolved_intervals"]}),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    type Spec<'a> = Vec<Vec<(&'a str, u64, u64, &'a str, f64)>>;
    struct Sequences(BTreeMap<String, Vec<u8>>);
    impl SequenceIndex for Sequences {
        fn get_sequence_length(&self, n: &str) -> io::Result<usize> {
            Ok(self.0[n].len())
        }
        fn fetch_sequence(&self, n: &str, s: i32, e: i32) -> io::Result<Vec<u8>> {
            Ok(self.0[n][s as usize..e as usize].to_vec())
        }
    }
    fn fixture(spec: Spec<'_>, dna: &[(&str, &str)]) -> (Genotypes, Threads, Sequences) {
        let panel = PanelIdentity {
            checksum_algorithm: "fixture".into(),
            sidecars: vec![],
        };
        let mut catalog = Catalog {
            version: 1,
            panel: panel.clone(),
            catalog_accepted: false,
            validation: "fixture".into(),
            sources: Vec::new(),
            groups: Vec::new(),
            occurrences: Vec::new(),
            links: Vec::new(),
            features: Vec::new(),
        };
        let seq = Sequences(
            dna.iter()
                .map(|(n, s)| (n.to_string(), s.as_bytes().to_vec()))
                .collect(),
        );
        let mut calls_rows = Vec::new();
        let mut occurrences = Vec::new();
        let mut axis = Axis {
            version: 1,
            coordinate_system: "fixture".into(),
            intervals: Vec::new(),
        };
        for (g, rows) in spec.into_iter().enumerate() {
            let group = format!("g{g}");
            let reference = catalog.occurrences.len();
            let mut ids = Vec::new();
            let mut bundles: BTreeMap<String, (f64, Vec<usize>)> = BTreeMap::new();
            for (path, start, end, strand, score) in rows {
                let id = catalog.occurrences.len();
                let identity = source_identity(path).unwrap();
                ids.push(id);
                let source = if let Some(s) = catalog.sources.iter().position(|s| s.path == path) {
                    s
                } else {
                    let s = catalog.sources.len();
                    catalog.sources.push(Source {
                        id: s,
                        path: path.into(),
                        length: seq.0[path].len() as u64,
                    });
                    s
                };
                let interval = Interval {
                    path: path.into(),
                    start,
                    end,
                    strand: Some(strand.into()),
                };
                catalog.occurrences.push(Occurrence {
                    id,
                    group: g,
                    source,
                    interval: interval.clone(),
                    fully_contained_anchors: 0,
                    anchor_status: "fixture".into(),
                    feature_multiplicities: BTreeMap::new(),
                });
                occurrences.push(SourceOccurrence {
                    id,
                    group: group.clone(),
                    identity: identity.clone(),
                    interval,
                });
                bundles
                    .entry(identity)
                    .or_insert((score, Vec::new()))
                    .1
                    .push(id);
            }
            let bundles: Vec<_> = bundles
                .into_iter()
                .map(|(identity, (score, source_occurrences))| Bundle {
                    identity,
                    source_occurrences,
                    physical_multiplicities: BTreeMap::new(),
                    physical_feature_locations: 1,
                    modeled_features: 1,
                    positive_features: 1,
                    score,
                    mean_deviance: 0.0,
                    expected_total: 1.0,
                })
                .collect();
            let best = bundles
                .iter()
                .map(|b| b.score)
                .fold(f64::INFINITY, f64::min);
            let best_bundles: Vec<_> = bundles
                .iter()
                .enumerate()
                .filter(|(_, b)| b.score == best)
                .map(|(i, _)| i)
                .collect();
            calls_rows.push(GenotypeCall {
                group: group.clone(),
                status: if best_bundles.len() > 1 {
                    "tied"
                } else {
                    "informative"
                }
                .into(),
                factors: vec![],
                retained_features: 1,
                positive_features: 1,
                observed_total: 1.0,
                background_only_score: 0.0,
                bundles,
                best_bundles,
                score_gap: None,
            });
            catalog.groups.push(Group {
                id: group.clone(),
                scaffold: None,
                occurrences: ids,
            });
            let r = &catalog.occurrences[reference].interval;
            axis.intervals.push(AxisInterval {
                component: "axis".into(),
                start: r.start,
                end: r.end,
                group,
                reference_occurrence: reference,
                reference_strand: "+".into(),
                orientations: BTreeMap::new(),
            });
        }
        let calls = Genotypes {
            version: 1,
            model: super::super::genotype::MODEL.into(),
            panel,
            count_policy: COUNT_POLICY.into(),
            sample_payload_checksum: Some("1111111111111111".into()),
            catalog_payload_checksum: Some("2222222222222222".into()),
            observations: None,
            artifact_checksum_algorithm:
                "fnv1a64-payload-v1; sample=bincode-payload; catalog=compact-json-payload".into(),
            feature_details_included: false,
            experimental: true,
            catalog_accepted: false,
            ploidy: 1,
            parameters: Parameters {
                haploid_depth: 10.0,
                background: 0.1,
                max_mean_deviance: 10.0,
            },
            assumptions: "fixture".into(),
            read_lengths: BTreeMap::new(),
            global_factors_used: 0,
            excluded_features: 0,
            excluded_features_by_reason: BTreeMap::new(),
            source_occurrences: occurrences,
            calls: calls_rows,
        };
        let threads = super::super::threading::thread(&catalog, &calls, axis, 0.0).unwrap();
        validate(
            &calls,
            &threads,
            &seq.0
                .iter()
                .map(|(n, b)| (n.clone(), b.len() as u64))
                .collect(),
        )
        .unwrap();
        (calls, threads, seq)
    }
    fn strings(r: &Reconstruction) -> Vec<String> {
        String::from_utf8(r.fasta.clone())
            .unwrap()
            .split('>')
            .skip(1)
            .map(|s| s.lines().skip(1).collect())
            .collect()
    }
    #[test]
    fn exact_overlap_gap_policy_reverse_and_provenance() {
        let (calls, threads, seq) = fixture(
            vec![
                vec![("A#0#chr", 0, 5, "+", 0.0)],
                vec![("A#0#chr", 3, 8, "+", 0.0)],
                vec![("A#0#chr", 10, 12, "+", 0.0)],
            ],
            &[("A#0#chr", "ACGTACGTTTGA")],
        );
        let split = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&split), vec!["ACGTACGT", "GA"]);
        assert_eq!(split.provenance["emitted_bp"], 10);
        assert_eq!(split.provenance["source_bridged_imputed_bp"], 0);
        let copy = reconstruct(&calls, &threads, &seq, true).unwrap();
        assert_eq!(strings(&copy), vec!["ACGTACGTTTGA"]);
        assert_eq!(copy.provenance["source_bridged_imputed_bp"], 2);
        assert_eq!(copy.provenance["blocks"][0]["scored_source_union_bp"], 10);
        assert_eq!(copy.provenance["blocks"][0]["scored_axis_union_bp"], 10);
        assert_eq!(copy.provenance["blocks"][0]["axis_span_end"], 12);
        assert_eq!(copy.provenance["gap_ledger"][0]["source_gap_start"], 8);
        let (calls, threads, seq) = fixture(
            vec![
                vec![("A#0#chr", 0, 5, "+", 10.0), ("B#0#rev", 7, 12, "-", 0.0)],
                vec![("A#0#chr", 3, 8, "+", 10.0), ("B#0#rev", 4, 9, "-", 0.0)],
            ],
            &[("A#0#chr", "AAAAAAAAAAAA"), ("B#0#rev", "AAAACCGTACGA")],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&r), vec!["TCGTACGG"]);
        assert_eq!(r.provenance["emitted_bp"], 8);
        assert_eq!(
            r.provenance["blocks"][0]["extraction"],
            json!({"path":"B#0#rev","start":4,"end":12,"strand":"-"})
        );
    }
    #[test]
    fn sequence_equivalent_ties_not_ancestry_or_simultaneous_copy_deduplication() {
        let (calls, threads, seq) = fixture(
            vec![vec![
                ("A#0#chr", 0, 4, "+", 0.0),
                ("B#0#chr", 0, 4, "+", 0.0),
            ]],
            &[("A#0#chr", "ACGT"), ("B#0#chr", "ACGT")],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&r), vec!["ACGT"]);
        assert_eq!(
            r.provenance["blocks"][0]["kind"],
            "sequence-equivalent-unphased-fragment"
        );
        assert_eq!(
            r.provenance["blocks"][0]["source_states"]
                .as_array()
                .unwrap()
                .len(),
            2
        );
        let different = Sequences(BTreeMap::from([
            ("A#0#chr".into(), b"ACGT".to_vec()),
            ("B#0#chr".into(), b"ACGA".to_vec()),
        ]));
        let r = reconstruct(&calls, &threads, &different, false).unwrap();
        assert!(r.fasta.is_empty());
        assert_eq!(r.provenance["emitted_bp"], 0);
        let (calls, threads, seq) = fixture(
            vec![vec![
                ("A#0#chr", 0, 4, "+", 0.0),
                ("A#0#chr", 0, 4, "+", 0.0),
            ]],
            &[("A#0#chr", "ACGT")],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&r), vec!["ACGT"]);
        assert_eq!(
            r.provenance["blocks"][0]["source_states"][0]["source_occurrences"],
            json!([0, 1])
        );
        let (calls, threads, seq) = fixture(
            vec![vec![
                ("A#0#chr", 0, 4, "+", 0.0),
                ("A#0#copy", 0, 4, "+", 0.0),
            ]],
            &[("A#0#chr", "ACGT"), ("A#0#copy", "ACGT")],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert!(r.fasta.is_empty());
        assert_eq!(
            r.unresolved["intervals"][0]["reason"],
            "multicopy-placement-tie-not-single-allele"
        );
    }
    #[test]
    fn path_breaks_and_corrupt_linkage_states_blocks_rejected() {
        let (calls, threads, seq) = fixture(
            vec![
                vec![("A#0#chr", 0, 4, "+", 0.0), ("B#0#chr", 0, 4, "+", 10.0)],
                vec![("A#0#chr", 4, 8, "+", 10.0), ("B#0#chr", 4, 8, "+", 0.0)],
            ],
            &[("A#0#chr", "ACGTAAAA"), ("B#0#chr", "TTTTGGCC")],
        );
        assert_eq!(
            strings(&reconstruct(&calls, &threads, &seq, true).unwrap()),
            vec!["ACGT", "GGCC"]
        );
        let lengths = seq
            .0
            .iter()
            .map(|(n, b)| (n.clone(), b.len() as u64))
            .collect();
        let original = serde_json::to_value(&threads).unwrap();
        for (pointer, value) in [
            ("/sample_payload_checksum", json!("3333333333333333")),
            ("/version", json!(2)),
            ("/intervals/0/states/0/emission", json!(999)),
            ("/intervals/0/states/0/bundle", json!(999)),
            ("/intervals/0/states/0/interval/end", json!(99)),
            ("/intervals/0/states/0/strand", json!("-")),
            ("/intervals/0/optimal_states", json!([999])),
            ("/blocks/0/interval_indices", json!([0, 1])),
            ("/intervals/0/states", json!([])),
        ] {
            let mut bad = original.clone();
            *bad.pointer_mut(pointer).unwrap() = value;
            assert!(
                validate(&calls, &serde_json::from_value(bad).unwrap(), &lengths).is_err(),
                "{pointer}"
            );
        }
        let mut bad = serde_json::to_value(&calls).unwrap();
        bad["calls"][0]["bundles"][0]["source_occurrences"] = json!([999]);
        assert!(validate(&serde_json::from_value(bad).unwrap(), &threads, &lengths).is_err());
    }
    #[test]
    fn reference_gap_splits_deduplicate_physical_overlap_in_both_orientations() {
        for (strand, first, second, expected) in [
            ("+", (0, 4), (2, 6), vec!["AACC", "GG"]),
            ("-", (6, 10), (4, 8), vec!["TTAA", "CC"]),
        ] {
            let (calls, threads, seq) = fixture(
                vec![
                    vec![
                        ("A#0#chr", 0, 4, "+", 10.0),
                        ("B#0#chr", first.0, first.1, strand, 0.0),
                    ],
                    vec![
                        ("A#0#chr", 6, 10, "+", 10.0),
                        ("B#0#chr", second.0, second.1, strand, 0.0),
                    ],
                ],
                &[("A#0#chr", "AAAAAAAAAA"), ("B#0#chr", "AACCGGTTAA")],
            );
            assert_eq!(threads.blocks.len(), 1);
            let r = reconstruct(&calls, &threads, &seq, false).unwrap();
            assert_eq!(strings(&r), expected);
            assert_eq!(r.provenance["emitted_bp"], 6);
            assert_eq!(r.provenance["suppressed_within_block_overlap_bp"], 2);
            assert!(r.provenance["gap_ledger"][0]["source_gap_start"].is_null());
            assert!(r.provenance["gap_ledger"][0]["source_gap_end"].is_null());
            assert_eq!(
                r.provenance["gap_ledger"][0]["source_overlap_end"]
                    .as_u64()
                    .unwrap()
                    - r.provenance["gap_ledger"][0]["source_overlap_start"]
                        .as_u64()
                        .unwrap(),
                2
            );
            assert_eq!(
                r.provenance["blocks"][1]["suppressed_within_block_overlap_bp"],
                2
            );
            assert_eq!(r.provenance["blocks"][1]["axis_span_start"], 6);
            assert_eq!(
                r.provenance["cross_block_reused_source_overlaps"],
                json!([])
            );
        }
    }
    #[test]
    fn disconnected_template_reuse_stays_visible_but_distinct_physical_copies_survive() {
        use super::super::sequence_evaluation::{evaluate, Fasta};
        let (calls, threads, seq) = fixture(
            vec![
                vec![("A#0#chr", 0, 4, "+", 10.0), ("B#0#chr", 0, 4, "+", 0.0)],
                vec![("A#0#chr", 4, 8, "+", 10.0), ("B#0#chr", 0, 4, "+", 0.0)],
            ],
            &[("A#0#chr", "AAAAAAAA"), ("B#0#chr", "ACGT")],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&r), vec!["ACGT", "ACGT"]);
        assert_eq!(
            r.provenance["cross_block_reused_source_overlaps"][0]["overlap_bp"],
            4
        );
        let query: Fasta = BTreeMap::from([
            ("q1".into(), b"ACGT".to_vec()),
            ("q2".into(), b"ACGT".to_vec()),
        ]);
        let truth: Fasta = BTreeMap::from([("t1".into(), b"ACGT".to_vec())]);
        let paf="q1\t4\t0\t4\t+\tt1\t4\t0\t4\t4\t4\t60\tcg:Z:4M\nq2\t4\t0\t4\t+\tt1\t4\t0\t4\t4\t4\t60\tcg:Z:4M\n";
        let e = evaluate(&query, &truth, paf.as_bytes()).unwrap();
        assert_eq!(e["query_coverage"], 0.5);
        assert_eq!(e["counts"]["matches"], 4);
        let (calls, threads, seq) = fixture(
            vec![
                vec![("A#0#chr", 0, 4, "+", 10.0), ("B#0#copy1", 0, 4, "+", 0.0)],
                vec![("A#0#chr", 4, 8, "+", 10.0), ("B#0#copy2", 0, 4, "+", 0.0)],
            ],
            &[
                ("A#0#chr", "AAAAAAAA"),
                ("B#0#copy1", "ACGT"),
                ("B#0#copy2", "ACGT"),
            ],
        );
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert_eq!(strings(&r), vec!["ACGT", "ACGT"]);
        assert_eq!(
            r.provenance["cross_block_reused_source_overlaps"],
            json!([])
        );
        let mut truth = truth;
        truth.insert("t2".into(), b"ACGT".to_vec());
        let paf = paf
            .lines()
            .enumerate()
            .map(|(i, l)| {
                format!(
                    "{}\n",
                    if i == 1 {
                        l.replace("t1", "t2")
                    } else {
                        l.into()
                    }
                )
            })
            .collect::<String>();
        let e = evaluate(&query, &truth, paf.as_bytes()).unwrap();
        assert_eq!(e["query_coverage"], 1.0);
        assert_eq!(e["truth_coverage"], 1.0);
        assert_eq!(e["counts"]["matches"], 8);
    }
    #[test]
    fn unknown_or_poor_fit_tie_stays_unresolved_and_partial_failure_removes_results() {
        let (calls, mut threads, seq) = fixture(
            vec![vec![
                ("A#0#chr", 0, 4, "+", 0.0),
                ("B#0#chr", 0, 4, "+", 0.0),
            ]],
            &[("A#0#chr", "ACGT"), ("B#0#chr", "ACGT")],
        );
        threads.intervals[0].states[1].strand = None;
        let r = reconstruct(&calls, &threads, &seq, false).unwrap();
        assert!(r.fasta.is_empty());
        threads.intervals[0].states[1].strand = Some("+".into());
        threads.intervals[0].states[1].mean_deviance = 11.0;
        assert!(reconstruct(&calls, &threads, &seq, false)
            .unwrap()
            .fasta
            .is_empty());
        let dir = tempfile::tempdir().unwrap();
        let out = dir.path().join("failed");
        let names = [
            "reconstruction.fa",
            "provenance.json",
            "unresolved.json",
            "unresolved.bed",
        ];
        assert!(with_output_model(&out, MODEL, || {
            for name in names {
                atomic_write(&out.join(name), b"partial result")?;
            }
            std::fs::write(out.join("provenance.incomplete"), b"partial JSON")?;
            Err(invalid("injected after output writes"))
        })
        .is_err());
        for name in names {
            assert!(!out.join(name).exists());
        }
        assert!(!out.join("provenance.incomplete").exists());
        let manifest: Value = read_json(&out.join("manifest.json")).unwrap();
        assert_eq!(manifest["status"], "failed");
        assert!(with_output_model(&out, MODEL, || Ok(json!({}))).is_err());
    }
}
