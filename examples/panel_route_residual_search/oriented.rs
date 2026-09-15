//! Exact-cut oriented reciprocal construction. Never move a canonical boundary.
use super::*;

/// Search the actual source index in the requested traversal view, not the RC
/// word of a port whose anchor happens to be adjacent to the required cut.
fn exact_port(
    e: &mut routes::Evaluator<'_>,
    source: usize,
    cut: u64,
    reverse: bool,
    meter: &mut Meter,
) -> io::Result<Option<geometry::Port>> {
    let k = e.graph.k;
    let shift = if reverse { k - k / 2 } else { k / 2 };
    let Some(anchor) = cut.checked_sub(shift) else {
        return Ok(None);
    };
    let mut file = e.ports.source_file(e.graph, source)?;
    meter.source_index_opens += 1;
    let (mut lo, mut hi) = (0, e.graph.lanes[source].port_count);
    while lo < hi {
        let mid = lo + (hi - lo) / 2;
        let p = geometry::read_port(&mut file, k, mid, meter)?;
        let a = p.cut - if p.reverse { k - k / 2 } else { k / 2 };
        if (a, p.reverse) < (anchor, reverse) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    if lo == e.graph.lanes[source].port_count {
        return Ok(None);
    }
    let p = geometry::read_port(&mut file, k, lo, meter)?;
    Ok((p.cut == cut && p.reverse == reverse).then_some(p))
}

/// Checks ALL normalized seams, including those inside transported material.
/// None means joins supported, not complete assignment/capacity certification.
pub fn seam_rejection(
    e: &mut routes::Evaluator<'_>,
    assignment: &routes::Assignment,
    meter: &mut Meter,
) -> io::Result<Option<serde_json::Value>> {
    for (slot, route) in assignment.routes.iter().enumerate() {
        meter.work(route.segments.len() as u64)?;
        let normalized = route.normalize();
        for (seam, pair) in normalized.segments.windows(2).enumerate() {
            let exit = if pair[0].reverse {
                pair[0].start
            } else {
                pair[0].end
            };
            let entry = if pair[1].reverse {
                pair[1].end
            } else {
                pair[1].start
            };
            let left = exact_port(e, pair[0].source, exit, pair[0].reverse, meter)?;
            let right = exact_port(e, pair[1].source, entry, pair[1].reverse, meter)?;
            let reason = match (&left, &right) {
                (None, _) => Some("undeclared outgoing switch"),
                (_, None) => Some("undeclared incoming switch"),
                (Some(a), Some(b)) if a.word != b.word => {
                    Some("switch words/oriented cuts disagree")
                }
                _ => None,
            };
            if let Some(reason) = reason {
                return Ok(Some(
                    json!({"slot":slot,"seam":seam,"pair":pair,"left":left,"right":right,"reason":reason,"canonical_cuts_shifted":false}),
                ));
            }
            for port in [left.unwrap(), right.unwrap()] {
                meter.work(e.graph.k * 2)?;
                let anchor = port.cut
                    - if port.reverse {
                        e.graph.k - e.graph.k / 2
                    } else {
                        e.graph.k / 2
                    };
                let mut dna = e.sources.fetch(port.source, anchor, anchor + e.graph.k)?;
                if port.reverse {
                    dna = impg::graph::reverse_complement(&dna);
                }
                // Source/index provenance corruption is fatal, not geometry pruning.
                ensure(dna == port.word, "switch DNA provenance mismatch")?;
            }
        }
    }
    Ok(None)
}

/// Raw simultaneous transport at EXACT occupied canonical spans. This is not
/// accepted geometry until seam_rejection and public validation both succeed.
pub fn reciprocal_proposal(
    base: &routes::Assignment,
    slot: usize,
    lo: u64,
    hi: u64,
    donor: &routes::Segment,
    other: usize,
) -> io::Result<Option<routes::Assignment>> {
    if slot == other {
        return Ok(None);
    }
    let mut offset = 0;
    for occupied in &base.routes[other].segments {
        if occupied.source == donor.source
            && occupied.start <= donor.start
            && donor.end <= occupied.end
        {
            let start = offset
                + if occupied.reverse {
                    occupied.end - donor.end
                } else {
                    donor.start - occupied.start
                };
            let end = start + donor.end - donor.start;
            let mut removed = geometry::crop(&base.routes[slot], lo, hi);
            if occupied.reverse != donor.reverse {
                removed.reverse();
                for segment in &mut removed {
                    segment.reverse = !segment.reverse;
                }
            }
            let mut result = base.clone();
            result.routes[slot] =
                geometry::splice(&base.routes[slot], lo, hi, vec![donor.clone()])?;
            result.routes[other] = geometry::splice(&base.routes[other], start, end, removed)?;
            return Ok(Some(result));
        }
        offset += occupied.end - occupied.start;
    }
    Ok(None)
}
