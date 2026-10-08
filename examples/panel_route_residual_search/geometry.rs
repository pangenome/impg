//! Single-donor two-ended geometry. No all-pairs port or patch matrix.
use super::*;
use std::io::{Read, Seek, SeekFrom};
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Port {
    pub word: Vec<u8>,
    pub source: usize,
    pub cut: u64,
    pub reverse: bool,
}
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct Region {
    pub slot: usize,
    pub left: u64,
    pub right: u64,
    pub level: u8,
}
impl Region {
    pub fn refinements(&self, length: u64, max_level: u8) -> Vec<Self> {
        if self.level >= max_level {
            return vec![];
        }
        let shift = ((self.right - self.left) / 2).max(1);
        let mut out = Vec::new();
        for (left, right) in [
            (self.left.saturating_sub(shift), self.right),
            (self.left.saturating_add(shift), self.right),
            (self.left, self.right.saturating_sub(shift)),
            (self.left, (self.right + shift).min(length)),
        ] {
            if left < right && right <= length {
                let r = Self {
                    slot: self.slot,
                    left,
                    right,
                    level: self.level + 1,
                };
                if !out.contains(&r) {
                    out.push(r);
                }
            }
        }
        out
    }
}
/// Crop in molecule coordinates; source coordinates stay canonical forward.
pub fn crop(route: &routes::Route, lo: u64, hi: u64) -> Vec<routes::Segment> {
    let mut offset = 0;
    let mut out = Vec::new();
    for p in &route.segments {
        let end = offset + p.end - p.start;
        let a = lo.max(offset);
        let b = hi.min(end);
        if a < b {
            out.push(routes::Segment {
                source: p.source,
                reverse: p.reverse,
                start: if p.reverse {
                    p.end - (b - offset)
                } else {
                    p.start + a - offset
                },
                end: if p.reverse {
                    p.end - (a - offset)
                } else {
                    p.start + b - offset
                },
            });
        }
        offset = end;
    }
    out
}
pub fn splice(
    route: &routes::Route,
    lo: u64,
    hi: u64,
    pieces: Vec<routes::Segment>,
) -> io::Result<routes::Route> {
    let mut segments = crop(route, 0, lo);
    segments.extend(pieces);
    segments.extend(crop(route, hi, route.length()?));
    Ok(routes::Route { segments }.normalize())
}
#[derive(Debug)]
pub enum Reciprocal {
    Supported(routes::Assignment),
    UnsupportedOppositeStrand,
    NotContained,
}
pub fn reciprocal(
    base: &routes::Assignment,
    slot: usize,
    lo: u64,
    hi: u64,
    donor: &routes::Segment,
    other: usize,
) -> io::Result<Reciprocal> {
    if slot == other {
        return Ok(Reciprocal::NotContained);
    }
    let mut offset = 0;
    for p in &base.routes[other].segments {
        if p.source == donor.source && p.start <= donor.start && donor.end <= p.end {
            // At odd k, opposite port cuts differ by one canonical base. Reversing
            // a cut subwalk is not a valid reciprocal construction in general.
            if p.reverse != donor.reverse {
                return Ok(Reciprocal::UnsupportedOppositeStrand);
            }
            let a = offset
                + if p.reverse {
                    p.end - donor.end
                } else {
                    donor.start - p.start
                };
            let b = a + donor.end - donor.start;
            let removed = crop(&base.routes[slot], lo, hi);
            let mut result = base.clone();
            result.routes[slot] = splice(&base.routes[slot], lo, hi, vec![donor.clone()])?;
            result.routes[other] = splice(&base.routes[other], a, b, removed)?;
            return Ok(Reciprocal::Supported(result));
        }
        offset += p.end - p.start;
    }
    Ok(Reciprocal::NotContained)
}
pub fn read_port(file: &mut File, k: u64, index: u64, meter: &mut Meter) -> io::Result<Port> {
    meter.work(1)?;
    meter.port_reads += 1;
    file.seek(SeekFrom::Start(
        index
            .checked_mul(k + 17)
            .ok_or_else(|| invalid("port offset overflow"))?,
    ))?;
    let mut word = vec![0; k as usize];
    file.read_exact(&mut word)?;
    let mut s = [0; 8];
    let mut a = [0; 8];
    let mut r = [0; 1];
    file.read_exact(&mut s)?;
    file.read_exact(&mut a)?;
    file.read_exact(&mut r)?;
    ensure(
        r[0] <= 1 && word.iter().all(|x| b"ACGT".contains(x)),
        "invalid verified port record",
    )?;
    let reverse = r[0] == 1;
    Ok(Port {
        word,
        source: usize::try_from(u64::from_le_bytes(s)).map_err(io::Error::other)?,
        cut: u64::from_le_bytes(a) + if reverse { k - k / 2 } else { k / 2 },
        reverse,
    })
}
pub fn nearest(
    e: &mut routes::Evaluator<'_>,
    route: &routes::Route,
    position: u64,
    meter: &mut Meter,
) -> io::Result<Option<(u64, Port)>> {
    let mut offset = 0;
    for segment in &route.segments {
        meter.work(1)?;
        let length = segment.end - segment.start;
        if position <= offset + length {
            let local = (position - offset).min(length);
            let coordinate = if segment.reverse {
                segment.end - local
            } else {
                segment.start + local
            };
            let mut file = e.ports.source_file(e.graph, segment.source)?;
            meter.source_index_opens += 1;
            let count = e.graph.lanes[segment.source].port_count;
            let k = e.graph.k;
            let (mut lo, mut hi) = (0, count);
            // Source records are anchor-ordered, not cut-ordered across strands.
            let anchor = coordinate.saturating_sub(if segment.reverse { k - k / 2 } else { k / 2 });
            while lo < hi {
                let mid = lo + (hi - lo) / 2;
                let p = read_port(&mut file, k, mid, meter)?;
                let a = p.cut - if p.reverse { k - k / 2 } else { k / 2 };
                if a < anchor {
                    lo = mid + 1
                } else {
                    hi = mid
                }
            }
            let mut best: Option<(u64, Port)> = None;
            for index in lo.saturating_sub(3)..(lo + 4).min(count) {
                let p = read_port(&mut file, k, index, meter)?;
                if p.reverse != segment.reverse || p.cut < segment.start || p.cut > segment.end {
                    continue;
                }
                let at = offset
                    + if segment.reverse {
                        segment.end - p.cut
                    } else {
                        p.cut - segment.start
                    };
                if at == 0 || at == route.length()? {
                    continue;
                }
                if best.as_ref().is_none_or(|(old, _)| {
                    (at.abs_diff(position), at) < (old.abs_diff(position), *old)
                }) {
                    best = Some((at, p));
                }
            }
            return Ok(best);
        }
        offset += length;
    }
    Ok(None)
}
pub fn bucket(
    e: &mut routes::Evaluator<'_>,
    word: &[u8],
    meter: &mut Meter,
) -> io::Result<(u64, u64)> {
    let mut bounds = [0; 2];
    for (upper, result) in bounds.iter_mut().enumerate() {
        let (mut lo, mut hi) = (0, e.graph.port_count);
        while lo < hi {
            let mid = lo + (hi - lo) / 2;
            let p = read_port(&mut e.ports.global, e.graph.k, mid, meter)?;
            if p.word.as_slice() < word || (upper == 1 && p.word == word) {
                lo = mid + 1
            } else {
                hi = mid
            }
        }
        *result = lo;
    }
    Ok((bounds[0], bounds[1]))
}
