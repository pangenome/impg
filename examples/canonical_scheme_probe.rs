//! One-off probe: record-extraction conventions on REAL reads, comparing the
//! syng best-orientation dual-view extraction (the production record
//! re-derivation) with canonical-scheme extractions, to measure:
//!  (a) how many TaggedWalks each read yields under each convention,
//!  (b) whether the canonical rc view adds anything over the canonical
//!      forward view (frame-freeness check),
//!  (c) the total record-instance count per convention over N reads.
//!
//! Usage: canonical_scheme_probe <panel.syng> <reads.fastq.gz> [N]

use impg::genome_inference::mem_records;
use impg::syng::{SyncmerParams, SyngIndex};
use std::io::{BufRead, BufReader};

fn main() -> std::io::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let panel = SyngIndex::load(&args[1], SyncmerParams::default())?;
    let n: usize = args
        .get(3)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200);
    let (reader, _) = niffler::get_reader(Box::new(std::fs::File::open(&args[2])?))
        .map_err(std::io::Error::other)?;
    let mut lines = BufReader::new(reader).lines();
    let k = panel.syncmer_length_bp() as usize;

    let mut syng_dual_total = 0usize;
    let mut canon_fwd_total = 0usize;
    let mut canon_dual_total = 0usize;
    let mut reads_with_double = 0usize;
    let mut reads_seen = 0usize;
    while reads_seen < n {
        let Some(_name) = lines.next().transpose()? else { break };
        let seq = lines.next().transpose()?.unwrap();
        lines.next().transpose()?;
        lines.next().transpose()?;
        reads_seen += 1;

        let syng = mem_records::tagged_mem_records(&panel, seq.as_bytes())?;
        syng_dual_total += syng.len();

        let canon_fwd = canonical_records(&panel, seq.as_bytes(), k)?;
        canon_fwd_total += canon_fwd.len();
        let rc = revcomp(seq.as_bytes());
        let canon_rc = canonical_records(&panel, &rc, k)?;
        // Mirrored form of a canon_fwd record:
        let mirror = |r: &[(i32, u64)]| -> Vec<(i32, u64)> {
            let span = seq.len() as u64;
            r.iter()
                .rev()
                .map(|&(node, pos)| (-node, span - k as u64 - pos))
                .collect()
        };
        let mut dual: Vec<Vec<(i32, u64)>> = canon_fwd.clone();
        let mut added_by_rc = 0usize;
        for r in &canon_rc {
            let mirrored = mirror(r);
            if !dual.iter().any(|d| d == &mirrored) {
                dual.push(mirrored);
                added_by_rc += 1;
            }
        }
        canon_dual_total += dual.len();
        if added_by_rc > 0 {
            reads_with_double += 1;
        }
    }
    println!(
        "reads {reads_seen} syng_dual {syng_dual_total} canon_fwd {canon_fwd_total} \
         canon_dual {canon_dual_total} reads_with_rc_additions {reads_with_double}"
    );
    Ok(())
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|&b| match b {
            b'A' | b'a' => b'T',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'T' | b't' => b'A',
            _ => b'N',
        })
        .collect()
}

/// Canonical-scheme records of one sequence, in its own frame: anchors =
/// canonical-qualified matched syncmer positions; MEMs = maximal GBWT-walks.
fn canonical_records(
    panel: &SyngIndex,
    seq: &[u8],
    k: usize,
) -> std::io::Result<Vec<Vec<(i32, u64)>>> {
    let anchors = canonical_anchor_walk(panel, seq, k)?;
    if anchors.is_empty() {
        return Ok(Vec::new());
    }
    let walk: Vec<impg::syng::SyngWalkStep> = anchors
        .iter()
        .map(|&(node, pos)| impg::syng::SyngWalkStep {
            signed_node: node,
            bp_pos: pos,
        })
        .collect();
    let mut records = Vec::new();
    for mem in panel.gbwt_mems_for_walk(&walk)? {
        records.push(
            walk[mem.step_start..mem.step_end]
                .iter()
                .map(|s| (s.signed_node, s.bp_pos))
                .collect(),
        );
    }
    Ok(records)
}

/// Canonical-qualified matched anchor positions of `seq` in its own frame.
fn canonical_anchor_walk(
    panel: &SyngIndex,
    seq: &[u8],
    k: usize,
) -> std::io::Result<Vec<(i32, u64)>> {
    if seq.len() < k {
        return Ok(Vec::new());
    }
    let length = seq.len() as u64;
    let fwd: std::collections::BTreeMap<u64, i32> = mem_records::raw_matched_syncmers(panel, seq)?
        .into_iter()
        .map(|(node, pos)| (pos, node))
        .collect();
    let rc = revcomp(seq);
    let rc_sel: std::collections::BTreeMap<u64, i32> =
        mem_records::raw_matched_syncmers(panel, &rc)?
            .into_iter()
            .map(|(node, q)| (length - k as u64 - q, -node))
            .collect();
    let mut positions: std::collections::BTreeSet<u64> = fwd.keys().copied().collect();
    positions.extend(rc_sel.keys().copied());
    let mut out = Vec::new();
    for p in positions {
        let start = p as usize;
        let Some(window) = seq.get(start..start + k) else { continue };
        if !window.iter().all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T')) {
            continue;
        }
        let mut rc_window = window.to_vec();
        rc_window.reverse();
        for b in &mut rc_window {
            *b = match *b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                _ => b'A',
            };
        }
        let canonical_form_is_forward = window <= rc_window.as_slice();
        let chosen = if canonical_form_is_forward {
            fwd.get(&p)
        } else {
            rc_sel.get(&p)
        };
        if let Some(&node) = chosen {
            out.push((node, p));
        }
    }
    Ok(out)
}
