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
