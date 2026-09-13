//! Read-only fixed-width adapter over the evaluator's verified handles.
use super::*;
use std::fs::File;
use std::io::{Read, Seek, SeekFrom};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub(super) struct Port {
    pub word: Vec<u8>,
    pub source: usize,
    pub anchor: u64,
    pub reverse: bool,
}
impl Port {
    pub fn cut(&self, k: u64) -> io::Result<u64> {
        self.anchor
            .checked_add(if self.reverse { k - k / 2 } else { k / 2 })
            .ok_or_else(|| invalid("port cut overflow"))
    }
}
pub(super) fn offset(k: u64, index: u64) -> io::Result<u64> {
    ensure(k > 0, "zero port width")?;
    index
        .checked_mul(
            k.checked_add(17)
                .ok_or_else(|| invalid("port width overflow"))?,
        )
        .ok_or_else(|| invalid("port offset overflow"))
}
pub(super) fn read(file: &mut File, k: u64, count: u64, index: u64) -> io::Result<Port> {
    ensure(index < count, "port index outside declared interval")?;
    let start = offset(k, index)?;
    let end = start
        .checked_add(k)
        .and_then(|n| n.checked_add(17))
        .ok_or_else(|| invalid("port end overflow"))?;
    ensure(end <= file.metadata()?.len(), "truncated port index")?;
    // Always absolute: evaluator operations and cloned handles may share cursors.
    file.seek(SeekFrom::Start(start))?;
    let mut word = vec![0; usize::try_from(k).map_err(io::Error::other)?];
    file.read_exact(&mut word)?;
    let mut src = [0; 8];
    let mut anchor = [0; 8];
    let mut rev = [0];
    file.read_exact(&mut src)?;
    file.read_exact(&mut anchor)?;
    file.read_exact(&mut rev)?;
    ensure(
        rev[0] <= 1 && word.iter().all(|b| b"ACGT".contains(b)),
        "invalid port record",
    )?;
    Ok(Port {
        word,
        source: usize::try_from(u64::from_le_bytes(src)).map_err(io::Error::other)?,
        anchor: u64::from_le_bytes(anchor),
        reverse: rev[0] == 1,
    })
}
pub(super) fn verified(p: &Port, g: &routes::Graph) -> io::Result<()> {
    ensure(
        p.source < g.lanes.len()
            && p.anchor
                .checked_add(g.k)
                .is_some_and(|end| end <= g.lanes[p.source].length),
        "port source/span outside graph",
    )
}

/// Virtual positions count as work even when the reversed index is outside len.
/// u128 represents the 2^64 virtual domain without a shift-width overflow.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Permutation {
    len: u64,
    ordinal: u128,
}
impl Permutation {
    pub fn new(len: u64) -> Self {
        Self { len, ordinal: 0 }
    }
    pub(super) fn validate(&self, base: u64, count: u64) -> io::Result<()> {
        let domain = if self.len == 0 {
            0
        } else {
            (self.len as u128).next_power_of_two()
        };
        ensure(
            base.checked_add(self.len).is_some_and(|n| n <= count) && self.ordinal <= domain,
            "invalid persisted permutation",
        )
    }
    pub fn next(&mut self) -> Option<Option<u64>> {
        let bits = if self.len <= 1 {
            0
        } else {
            64 - (self.len - 1).leading_zeros()
        };
        let domain = if self.len == 0 { 0 } else { 1u128 << bits };
        if self.ordinal >= domain {
            return None;
        }
        let index = if bits == 0 {
            0
        } else {
            (self.ordinal as u64).reverse_bits() >> (64 - bits)
        };
        self.ordinal += 1;
        Some((index < self.len).then_some(index))
    }
}

/// Source positions alone interleave adjacent offsets. None reads a legacy v2
/// cursor; Some(prefix) retains that old visited prefix while scanning paired order.
/// Skipped old and virtual positions each consume one ordinary machine operation.
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct SourcePermutation {
    len: u64,
    ordinal: u128,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    old_prefix: Option<u128>,
}
impl SourcePermutation {
    pub fn new(len: u64) -> Self {
        Self {
            len,
            ordinal: 0,
            old_prefix: Some(0),
        }
    }
    pub(super) fn update(&mut self) -> io::Result<()> {
        ensure(
            self.old_prefix.is_none(),
            "source scan already uses paired order",
        )?;
        self.old_prefix = Some(self.ordinal);
        self.ordinal = 0;
        Ok(())
    }
    pub(super) fn paired(&self) -> bool {
        self.old_prefix.is_some()
    }
    pub(super) fn validate(&self, base: u64, count: u64) -> io::Result<()> {
        Permutation {
            len: self.len,
            ordinal: self.ordinal,
        }
        .validate(base, count)?;
        if let Some(prefix) = self.old_prefix {
            Permutation {
                len: self.len,
                ordinal: prefix,
            }
            .validate(base, count)?;
        }
        Ok(())
    }
    pub fn next(&mut self) -> Option<Option<u64>> {
        let Some(prefix) = self.old_prefix else {
            let mut old = Permutation {
                len: self.len,
                ordinal: self.ordinal,
            };
            let next = old.next();
            self.ordinal = old.ordinal;
            return next;
        };
        let bits = if self.len <= 1 {
            0
        } else {
            64 - (self.len - 1).leading_zeros()
        };
        let domain = if self.len == 0 { 0 } else { 1u128 << bits };
        if self.ordinal >= domain {
            return None;
        }
        let ordinal = self.ordinal as u64;
        let index = if bits <= 1 {
            ordinal
        } else {
            2 * ((ordinal / 2).reverse_bits() >> (65 - bits)) + ordinal % 2
        };
        let old_ordinal = if bits == 0 {
            0
        } else {
            index.reverse_bits() >> (64 - bits)
        };
        self.ordinal += 1;
        Some((index < self.len && old_ordinal as u128 >= prefix).then_some(index))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Write;
    #[test]
    fn permutation_exact_coverage_virtual_work_and_overflow() {
        for len in 0..260 {
            let mut p = Permutation::new(len);
            let mut found = BTreeSet::new();
            let mut work = 0;
            while let Some(index) = p.next() {
                work += 1;
                if let Some(i) = index {
                    assert!(found.insert(i));
                }
            }
            assert_eq!(found, (0..len).collect());
            assert_eq!(work, if len == 0 { 0 } else { len.next_power_of_two() });
        }
        let mut p = Permutation::new(u64::MAX);
        assert_eq!(p.next(), Some(Some(0)));
        assert_eq!(p.next(), Some(Some(1 << 63)));
        p.ordinal = u64::MAX as u128;
        assert_eq!(p.next(), Some(None));
        assert_eq!(p.next(), None);
        assert!(offset(0, 0).is_err());
        assert!(offset(u64::MAX, 0).is_err());
        assert!(offset(31, u64::MAX).is_err());
    }
    #[test]
    fn source_pair_order_breaks_parity_stall_with_fixed_width_records() {
        for len in [4276u64, 18204] {
            let mut file = tempfile::tempfile().unwrap();
            for index in 0..len {
                file.write_all(b"ACG").unwrap();
                file.write_all(&0u64.to_le_bytes()).unwrap();
                file.write_all(&(index / 2).to_le_bytes()).unwrap();
                file.write_all(&[(index % 2) as u8]).unwrap();
            }
            let mut old = Permutation::new(len);
            for _ in 0..len.next_power_of_two() / 2 {
                if let Some(index) = old.next().unwrap() {
                    assert!(!read(&mut file, 3, len, index).unwrap().reverse);
                }
            }
            assert!(
                read(&mut file, 3, len, old.next().unwrap().unwrap())
                    .unwrap()
                    .reverse
            );
            let mut paired = SourcePermutation::new(len);
            for reverse in [false, true] {
                assert_eq!(
                    read(&mut file, 3, len, paired.next().unwrap().unwrap())
                        .unwrap()
                        .reverse,
                    reverse
                );
            }
        }
    }
    #[test]
    fn source_pair_order_prefix_conservation_and_charged_skips() {
        for len in 0u64..130 {
            let domain = if len == 0 { 0 } else { len.next_power_of_two() };
            for boundary in (0..=domain).filter(|n| *n <= 7 || *n % 17 == 0 || *n == domain) {
                let mut old = Permutation::new(len);
                let mut visited = BTreeSet::new();
                for _ in 0..boundary {
                    if let Some(i) = old.next().unwrap() {
                        visited.insert(i);
                    }
                }
                let json = serde_json::to_value(&old).unwrap();
                let mut paired: SourcePermutation = serde_json::from_value(json).unwrap();
                paired.update().unwrap();
                for base in [0, 1, 3, u64::MAX - len] {
                    paired.validate(base, base + len).unwrap();
                    let mut p = paired.clone();
                    let mut remaining = BTreeSet::new();
                    let mut work = 0;
                    while let Some(index) = p.next() {
                        work += 1;
                        if let Some(i) = index {
                            assert!(!visited.contains(&i));
                            assert!(remaining.insert(base + i));
                        }
                        let mut resumed: SourcePermutation =
                            serde_json::from_slice(&serde_json::to_vec(&p).unwrap()).unwrap();
                        assert_eq!(p.clone().next(), resumed.next());
                    }
                    assert_eq!(
                        remaining,
                        (0..len)
                            .filter(|i| !visited.contains(i))
                            .map(|i| base + i)
                            .collect()
                    );
                    assert_eq!(
                        work, domain,
                        "old/virtual skips must be charged, not looped away"
                    );
                }
            }
        }
        for len in [1u64 << 63, (1u64 << 63) + 1, u64::MAX] {
            let domain = (len as u128).next_power_of_two();
            for boundary in [0, 1, 3, domain / 2 + 1, domain] {
                let mut p = SourcePermutation {
                    len,
                    ordinal: boundary,
                    old_prefix: None,
                };
                p.update().unwrap();
                p.validate(0, len).unwrap();
                for ordinal in [0, 1, 2, 3, domain / 2, domain - 2, domain - 1] {
                    p.ordinal = ordinal;
                    let mut q: SourcePermutation =
                        serde_json::from_slice(&serde_json::to_vec(&p).unwrap()).unwrap();
                    assert_eq!(p.next(), q.next());
                    p.validate(0, len).unwrap();
                }
                p.ordinal = domain;
                assert_eq!(p.next(), None);
            }
        }
        let mut p = SourcePermutation::new(u64::MAX);
        assert_eq!(p.next(), Some(Some(0)));
        assert_eq!(p.next(), Some(Some(1)));
        assert_eq!(p.next(), Some(Some(1 << 63)));
        assert_eq!(p.next(), Some(Some((1 << 63) + 1)));
        p.ordinal = u64::MAX as u128;
        assert_eq!(p.next(), Some(None));
        assert_eq!(p.next(), None);
        assert!(SourcePermutation::new(2)
            .validate(u64::MAX, u64::MAX)
            .is_err());
    }
    #[test]
    fn fixed_width_invalid_truncation_cuts_and_shared_seeks() {
        let mut f = tempfile::tempfile().unwrap();
        for reverse in [0, 1] {
            f.write_all(b"ACG").unwrap();
            f.write_all(&2u64.to_le_bytes()).unwrap();
            f.write_all(&7u64.to_le_bytes()).unwrap();
            f.write_all(&[reverse]).unwrap();
        }
        let mut clone = f.try_clone().unwrap();
        let a = read(&mut f, 3, 2, 0).unwrap();
        assert_eq!(a.cut(3).unwrap(), 8);
        let b = read(&mut clone, 3, 2, 1).unwrap();
        assert_eq!(b.cut(3).unwrap(), 9);
        assert_eq!(read(&mut f, 3, 2, 0).unwrap(), a);
        assert!(read(&mut f, 3, 2, 2).is_err());
        f.set_len(39).unwrap();
        assert!(read(&mut f, 3, 2, 1).is_err());
        f.seek(SeekFrom::Start(19)).unwrap();
        f.write_all(&[2]).unwrap();
        assert!(read(&mut f, 3, 2, 0).is_err());
        f.seek(SeekFrom::Start(19)).unwrap();
        f.write_all(&[0]).unwrap();
        f.seek(SeekFrom::Start(0)).unwrap();
        f.write_all(b"N").unwrap();
        assert!(read(&mut f, 3, 2, 0).is_err());
        assert!(Port {
            anchor: u64::MAX,
            ..a
        }
        .cut(3)
        .is_err());
    }
}
