//! Shared geometry-only continuation. Scoring/feedback remain in each host.
use super::*;
#[derive(Clone, Debug, Serialize)]
pub enum Task {
    Coarse {
        base: usize,
        slot: usize,
        level: u8,
        cell: u64,
    },
    Region {
        base: usize,
        region: Region,
        guided: bool,
    },
    Pairs {
        base: usize,
        region: Region,
        guided: bool,
        lo: u64,
        hi: u64,
        left: (u64, u64),
        right: (u64, u64),
        i: u64,
        j: u64,
    },
    Candidate {
        base: usize,
        region: Region,
        guided: bool,
        lo: u64,
        hi: u64,
        donor: routes::Segment,
        other: usize,
    },
}

impl Task {
    pub fn base(&self) -> usize {
        match self {
            Self::Coarse { base, .. }
            | Self::Region { base, .. }
            | Self::Pairs { base, .. }
            | Self::Candidate { base, .. } => *base,
        }
    }
}
pub trait TaskQueue {
    fn len(&self) -> usize;
    fn push_back(&mut self, task: Task);
}
impl TaskQueue for VecDeque<Task> {
    fn len(&self) -> usize {
        VecDeque::len(self)
    }
    fn push_back(&mut self, task: Task) {
        VecDeque::push_back(self, task)
    }
}
pub struct Proposal {
    pub base: usize,
    pub region: Region,
    pub guided: bool,
    pub lo: u64,
    pub hi: u64,
    pub other: usize,
    pub candidate: routes::Assignment,
}
pub struct Context<'a, 'b> {
    pub assignment: &'a routes::Assignment,
    pub e: &'a mut routes::Evaluator<'b>,
    pub meter: &'a mut Meter,
    pub tasks: &'a mut dyn TaskQueue,
    pub events: &'a mut File,
    pub unsupported_opposite_strand: &'a mut usize,
}
impl Context<'_, '_> {
    fn push(&mut self, t: Task) -> io::Result<()> {
        self.meter.work(1)?;
        ensure(
            self.tasks.len() < self.meter.limits.max_tasks,
            "budget: pending_tasks",
        )?;
        self.tasks.push_back(t);
        Ok(())
    }
    pub fn step(&mut self, t: Task) -> io::Result<Option<Proposal>> {
        self.meter.work(1)?;
        match t {
            Task::Coarse {
                base,
                slot,
                level,
                cell,
            } => {
                let lengths = self.assignment.routes.len();
                if slot >= lengths {
                    return Ok(None);
                }
                let n = 1u64 << level;
                let length = self.assignment.routes[slot].length()?;
                let region = Region {
                    slot,
                    left: length * cell / (2 * n),
                    right: length * (cell + 2) / (2 * n),
                    level,
                };
                let (next_slot, next_level, next_cell) = if cell + 1 < 2 * n - 1 {
                    (slot, level, cell + 1)
                } else if level < self.meter.limits.max_level {
                    (slot, level + 1, 0)
                } else {
                    (slot + 1, 0, 0)
                };
                if region.left < region.right {
                    self.push(Task::Region {
                        base,
                        region,
                        guided: false,
                    })?;
                }
                self.push(Task::Coarse {
                    base,
                    slot: next_slot,
                    level: next_level,
                    cell: next_cell,
                })?;
            }
            Task::Region {
                base,
                region,
                guided,
            } => {
                let route = &self.assignment.routes[region.slot];
                let left = geometry::nearest(self.e, route, region.left, &mut self.meter)?;
                let right = geometry::nearest(self.e, route, region.right, &mut self.meter)?;
                line(
                    &mut self.events,
                    &json!({"event":"geometry_attempt","baseline":base,"region":region,"guided":guided,"left":left,"right":right}),
                )?;
                if let (Some((lo, p)), Some((hi, q))) = (left, right) {
                    if lo < hi {
                        let left = geometry::bucket(self.e, &p.word, &mut self.meter)?;
                        let right = geometry::bucket(self.e, &q.word, &mut self.meter)?;
                        self.push(Task::Pairs {
                            base,
                            region,
                            guided,
                            lo,
                            hi,
                            left,
                            right,
                            i: left.0,
                            j: right.0,
                        })?;
                    }
                }
            }
            Task::Pairs {
                base,
                region,
                guided,
                lo,
                hi,
                left,
                right,
                i,
                j,
            } => {
                if i >= left.1 || j >= right.1 {
                    return Ok(None);
                }
                let a = geometry::read_port(
                    &mut self.e.ports.global,
                    self.e.graph.k,
                    i,
                    &mut self.meter,
                )?;
                let b = geometry::read_port(
                    &mut self.e.ports.global,
                    self.e.graph.k,
                    j,
                    &mut self.meter,
                )?;
                let (ni, nj) = if j + 1 < right.1 {
                    (i, j + 1)
                } else {
                    (i + 1, right.0)
                };
                self.push(Task::Pairs {
                    base,
                    region: region.clone(),
                    guided,
                    lo,
                    hi,
                    left,
                    right,
                    i: ni,
                    j: nj,
                })?;
                if a.source == b.source
                    && a.reverse == b.reverse
                    && if a.reverse {
                        a.cut > b.cut
                    } else {
                        a.cut < b.cut
                    }
                {
                    let donor = routes::Segment {
                        source: a.source,
                        start: a.cut.min(b.cut),
                        end: a.cut.max(b.cut),
                        reverse: a.reverse,
                    };
                    self.push(Task::Candidate {
                        base,
                        region,
                        guided,
                        lo,
                        hi,
                        donor,
                        other: 0,
                    })?;
                }
            }
            Task::Candidate {
                base,
                region,
                guided,
                lo,
                hi,
                donor,
                other,
            } => {
                let assignment = self.assignment;
                let slots = assignment.routes.len();
                self.meter.work(
                    4 * assignment
                        .routes
                        .iter()
                        .map(|r| r.segments.len() as u64)
                        .sum::<u64>()
                        + 8,
                )?;
                let candidate = if other == 0 {
                    let mut a = assignment.clone();
                    a.routes[region.slot] =
                        geometry::splice(&a.routes[region.slot], lo, hi, vec![donor.clone()])?;
                    Some(a)
                } else if self.meter.limits.geometry_b2 {
                    oriented::reciprocal_proposal(
                        assignment,
                        region.slot,
                        lo,
                        hi,
                        &donor,
                        other - 1,
                    )?
                } else {
                    match geometry::reciprocal(assignment, region.slot, lo, hi, &donor, other - 1)?
                    {
                        geometry::Reciprocal::Supported(a) => Some(a),
                        geometry::Reciprocal::NotContained => None,
                        geometry::Reciprocal::UnsupportedOppositeStrand => {
                            *self.unsupported_opposite_strand += 1;
                            line(
                                &mut self.events,
                                &json!({"event":"unsupported_opposite_strand_reciprocal","baseline":base,"region":region,"donor":donor,"other_slot":other-1,"biological_infeasibility":false}),
                            )?;
                            None
                        }
                    }
                };
                // The reciprocal successor is retained regardless of hypothetical single validity.
                if other < slots {
                    self.push(Task::Candidate {
                        base,
                        region: region.clone(),
                        guided,
                        lo,
                        hi,
                        donor: donor.clone(),
                        other: other + 1,
                    })?;
                }
                let Some(candidate) = candidate else {
                    return Ok(None);
                };
                return Ok(Some(Proposal {
                    base,
                    region,
                    guided,
                    lo,
                    hi,
                    other,
                    candidate,
                }));
            }
        }
        Ok(None)
    }
}
