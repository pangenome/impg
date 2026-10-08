//! B2 lazy left-seeded two-donor continuation; return hubs discover the right end.
use super::*;
#[derive(Clone, Debug, Serialize)]
pub struct Cursor {
    pub base: usize,
    pub slot: usize,
    level: u8,
    cell: u64,
    stage: u8,
    lo: u64,
    entries: (u64, u64),
    i: u64,
    entry: Option<geometry::Port>,
    exit1: u64,
    first: Option<routes::Segment>,
    middles: (u64, u64),
    j: u64,
    middle: Option<geometry::Port>,
    exit2: u64,
    second: Option<routes::Segment>,
    returns: (u64, u64),
    r: u64,
    returning: Option<geometry::Port>,
    recipient_segment: usize,
}
impl Cursor {
    pub fn new(base: usize, slot: usize) -> Self {
        Self {
            base,
            slot,
            level: 0,
            cell: 0,
            stage: 0,
            lo: 0,
            entries: (0, 0),
            i: 0,
            entry: None,
            exit1: 0,
            first: None,
            middles: (0, 0),
            j: 0,
            middle: None,
            exit2: 0,
            second: None,
            returns: (0, 0),
            r: 0,
            returning: None,
            recipient_segment: 0,
        }
    }
    fn next_seed(&mut self, max_level: u8) -> bool {
        self.cell += 1;
        if self.cell >= 2 * (1u64 << self.level) - 1 {
            self.cell = 0;
            self.level += 1;
        }
        self.stage = 0;
        self.level > max_level
    }
    /// One bounded state-machine transition. No partial route gets validated,
    /// scored, or pruned by donor labels. Even contiguous seams are retained.
    pub fn advance(
        &mut self,
        e: &mut routes::Evaluator<'_>,
        base: &routes::Assignment,
        meter: &mut Meter,
    ) -> io::Result<(bool, Option<(Region, routes::Assignment)>)> {
        meter.work(1)?;
        let route = &base.routes[self.slot];
        match self.stage {
            0 => {
                let position = route.length()? * self.cell / (2 * (1u64 << self.level));
                if let Some((lo, port)) = geometry::nearest(e, route, position, meter)? {
                    self.lo = lo;
                    self.entries = geometry::bucket(e, &port.word, meter)?;
                    self.i = self.entries.0;
                    self.stage = 1;
                } else {
                    return Ok((self.next_seed(meter.limits.max_level), None));
                }
            }
            1 => {
                if self.i == self.entries.1 {
                    return Ok((self.next_seed(meter.limits.max_level), None));
                }
                self.entry = Some(geometry::read_port(
                    &mut e.ports.global,
                    e.graph.k,
                    self.i,
                    meter,
                )?);
                self.i += 1;
                self.exit1 = 0;
                self.stage = 2;
            }
            2 => {
                let entry = self.entry.as_ref().unwrap();
                if self.exit1 == e.graph.lanes[entry.source].port_count {
                    self.stage = 1;
                    return Ok((false, None));
                }
                let mut file = e.ports.source_file(e.graph, entry.source)?;
                meter.source_index_opens += 1;
                let exit = geometry::read_port(&mut file, e.graph.k, self.exit1, meter)?;
                self.exit1 += 1;
                if let Some(piece) = positive(entry, &exit) {
                    self.first = Some(piece);
                    self.middles = geometry::bucket(e, &exit.word, meter)?;
                    self.j = self.middles.0;
                    self.stage = 3;
                }
            }
            3 => {
                if self.j == self.middles.1 {
                    self.stage = 2;
                    return Ok((false, None));
                }
                self.middle = Some(geometry::read_port(
                    &mut e.ports.global,
                    e.graph.k,
                    self.j,
                    meter,
                )?);
                self.j += 1;
                self.exit2 = 0;
                self.stage = 4;
            }
            4 => {
                let entry = self.middle.as_ref().unwrap();
                if self.exit2 == e.graph.lanes[entry.source].port_count {
                    self.stage = 3;
                    return Ok((false, None));
                }
                let mut file = e.ports.source_file(e.graph, entry.source)?;
                meter.source_index_opens += 1;
                let exit = geometry::read_port(&mut file, e.graph.k, self.exit2, meter)?;
                self.exit2 += 1;
                if let Some(piece) = positive(entry, &exit) {
                    self.second = Some(piece);
                    self.returns = geometry::bucket(e, &exit.word, meter)?;
                    self.r = self.returns.0;
                    self.stage = 5;
                }
            }
            5 => {
                if self.r == self.returns.1 {
                    self.stage = 4;
                    return Ok((false, None));
                }
                self.returning = Some(geometry::read_port(
                    &mut e.ports.global,
                    e.graph.k,
                    self.r,
                    meter,
                )?);
                self.r += 1;
                self.recipient_segment = 0;
                self.stage = 6;
            }
            6 => {
                if self.recipient_segment == route.segments.len() {
                    self.stage = 5;
                    return Ok((false, None));
                }
                let index = self.recipient_segment;
                self.recipient_segment += 1;
                meter.work(route.segments.len() as u64)?;
                let segment = &route.segments[index];
                let p = self.returning.as_ref().unwrap();
                if segment.source == p.source
                    && segment.reverse == p.reverse
                    && segment.start <= p.cut
                    && p.cut <= segment.end
                {
                    let offset = route.segments[..index]
                        .iter()
                        .map(|s| s.end - s.start)
                        .sum::<u64>();
                    let hi = offset
                        + if segment.reverse {
                            segment.end - p.cut
                        } else {
                            p.cut - segment.start
                        };
                    if self.lo < hi && hi < route.length()? {
                        meter.work(
                            4 * base
                                .routes
                                .iter()
                                .map(|r| r.segments.len() as u64)
                                .sum::<u64>()
                                + 8,
                        )?;
                        let mut candidate = base.clone();
                        candidate.routes[self.slot] = geometry::splice(
                            route,
                            self.lo,
                            hi,
                            vec![self.first.clone().unwrap(), self.second.clone().unwrap()],
                        )?;
                        return Ok((
                            false,
                            Some((
                                Region {
                                    slot: self.slot,
                                    left: self.lo,
                                    right: hi,
                                    level: self.level,
                                },
                                candidate,
                            )),
                        ));
                    }
                }
            }
            _ => return Err(invalid("invalid chain cursor stage")),
        }
        Ok((false, None))
    }
}
fn positive(a: &geometry::Port, b: &geometry::Port) -> Option<routes::Segment> {
    (a.source == b.source
        && a.reverse == b.reverse
        && if a.reverse {
            a.cut > b.cut
        } else {
            a.cut < b.cut
        })
    .then(|| routes::Segment {
        source: a.source,
        start: a.cut.min(b.cut),
        end: a.cut.max(b.cut),
        reverse: a.reverse,
    })
}
