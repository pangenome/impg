//! Count-only weighted collection FM index. No suffix locations survive construction.
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::io;

pub(crate) fn invalid(message: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, message.into())
}

/// Alternating even node tokens and odd internal gap tokens. 0 terminates the
/// tape; 1 separates MEMs. There is deliberately no gap before the first node.
pub fn encode_walk(walk: &[(i32, u64)]) -> io::Result<Vec<u64>> {
    let mut result = Vec::new();
    for (i, &(node, pos)) in walk.iter().enumerate() {
        if node == 0 || node == i32::MIN {
            return Err(invalid("invalid signed node"));
        }
        if i > 0 {
            let gap = pos
                .checked_sub(walk[i - 1].1)
                .filter(|g| *g > 0 && *g <= u32::MAX as u64)
                .ok_or_else(|| invalid("invalid walk spacing"))?;
            result.push(2 * gap + 1);
        }
        let zigzag = ((node as i64) << 1) ^ ((node as i64) >> 63);
        result.push(2 + 2 * zigzag as u64);
    }
    validate_pattern(&result)?;
    Ok(result)
}

pub(crate) fn validate_pattern(pattern: &[u64]) -> io::Result<()> {
    if pattern.is_empty() || pattern.len() % 2 == 0 {
        return Err(invalid("feature must begin and end on a node"));
    }
    for (i, &token) in pattern.iter().enumerate() {
        if i % 2 == 0 {
            let z = token.saturating_sub(2) / 2;
            if token < 4 || token % 2 != 0 || z > 2 * i32::MAX as u64 {
                return Err(invalid("invalid node token"));
            }
        } else if token < 3 || token % 2 != 1 || token / 2 > u32::MAX as u64 {
            return Err(invalid("invalid gap token"));
        }
    }
    Ok(())
}

pub fn reverse_complement(pattern: &[u64]) -> Vec<u64> {
    pattern
        .iter()
        .rev()
        .enumerate()
        .map(|(i, &token)| {
            if i % 2 == 1 {
                token
            } else {
                let z = (token - 2) / 2;
                let node = ((z >> 1) as i64) ^ -((z & 1) as i64);
                let node = -node;
                2 + 2 * (((node << 1) ^ (node >> 63)) as u64)
            }
        })
        .collect()
}

pub fn canonical(pattern: &[u64]) -> Vec<u64> {
    let rc = reverse_complement(pattern);
    if rc.as_slice() < pattern {
        rc
    } else {
        pattern.to_vec()
    }
}

#[derive(Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct WeightedBwt {
    bwt: Vec<u64>,
    /// Cumulative weight in suffix order, NOT a locating/sample-read table.
    weights: Vec<u64>,
    #[serde(skip)]
    ranks: BTreeMap<u64, (usize, Vec<usize>)>,
}

impl WeightedBwt {
    /// Prefix-doubling suffix construction: O(N log² N), O(N) temporary space.
    /// The tape contains each distinct record once, never all its substrings.
    pub fn build(records: &BTreeMap<Vec<u64>, u64>) -> io::Result<Self> {
        let mut text = Vec::new();
        let mut weight = Vec::new();
        for (record, &w) in records {
            validate_pattern(record)?;
            if w == 0 {
                return Err(invalid("zero MEM multiplicity"));
            }
            text.extend(record);
            weight.extend(std::iter::repeat_n(w, record.len()));
            text.push(1);
            weight.push(0);
        }
        text.push(0);
        weight.push(0);
        let n = text.len();
        let mut sa: Vec<usize> = (0..n).collect();
        let mut rank = text.clone();
        let mut next = vec![0; n];
        let mut width = 1;
        loop {
            let key = |i: usize| {
                (
                    rank[i] + 1,
                    if i + width < n {
                        rank[i + width] + 1
                    } else {
                        0
                    },
                )
            };
            sa.sort_unstable_by_key(|&i| key(i));
            next[sa[0]] = 0;
            for j in 1..n {
                next[sa[j]] = next[sa[j - 1]] + u64::from(key(sa[j - 1]) != key(sa[j]));
            }
            std::mem::swap(&mut rank, &mut next);
            if rank[sa[n - 1]] == (n - 1) as u64 {
                break;
            }
            width = width
                .checked_mul(2)
                .ok_or_else(|| invalid("suffix width overflow"))?;
        }
        let bwt = sa
            .iter()
            .map(|&i| text[if i == 0 { n - 1 } else { i - 1 }])
            .collect();
        let mut weights = vec![0u64];
        for &i in &sa {
            weights.push(
                weights
                    .last()
                    .unwrap()
                    .checked_add(weight[i])
                    .ok_or_else(|| invalid("count overflow"))?,
            );
        }
        let mut result = Self {
            bwt,
            weights,
            ranks: BTreeMap::new(),
        };
        result.rebuild()?;
        Ok(result)
    }

    pub fn rebuild(&mut self) -> io::Result<()> {
        if self.bwt.is_empty()
            || self.weights.len() != self.bwt.len() + 1
            || self.weights[0] != 0
            || self.weights.windows(2).any(|w| w[0] > w[1])
            || self.bwt.iter().filter(|&&t| t == 0).count() != 1
        {
            return Err(invalid("invalid weighted BWT structure"));
        }
        self.ranks.clear();
        for (i, &token) in self.bwt.iter().enumerate() {
            self.ranks.entry(token).or_default().1.push(i);
        }
        let mut total = 0;
        for (start, positions) in self.ranks.values_mut() {
            *start = total;
            total += positions.len();
        }
        Ok(())
    }

    /// Oriented count on the actual stored strings (canonical strings for a sample).
    pub fn count_oriented(&self, pattern: &[u64]) -> io::Result<u64> {
        validate_pattern(pattern)?;
        let (mut lo, mut hi) = (0, self.bwt.len());
        for token in pattern.iter().rev() {
            let Some((start, positions)) = self.ranks.get(token) else {
                return Ok(0);
            };
            lo = start + positions.partition_point(|&p| p < lo);
            hi = start + positions.partition_point(|&p| p < hi);
            if lo == hi {
                return Ok(0);
            }
        }
        Ok(self.weights[hi] - self.weights[lo])
    }

    /// RC-orbit count: C_oriented(f) + C_oriented(rc(f)), palindromes once.
    pub fn count(&self, pattern: &[u64]) -> io::Result<u64> {
        let count = self.count_oriented(pattern)?;
        let rc = reverse_complement(pattern);
        if rc == pattern {
            Ok(count)
        } else {
            count
                .checked_add(self.count_oriented(&rc)?)
                .ok_or_else(|| invalid("count overflow"))
        }
    }

    /// Enumerate positive adjacent-node feature counts directly from the FM index.
    /// Each suffix row supplies one possible three-token occurrence. Two LF
    /// steps recover its start row, whose weight is the same weight used by
    /// `count_oriented`. No read identities, suffix locations or record tape
    /// are reconstructed. Separators cannot form a valid node-gap-node triple.
    pub fn observed_pairs(&self) -> io::Result<BTreeMap<[u64; 3], u64>> {
        if self.bwt.is_empty()
            || self.ranks.values().map(|(_, p)| p.len()).sum::<usize>() != self.bwt.len()
        {
            return Err(invalid("weighted BWT must be rebuilt before enumeration"));
        }
        let mut lf = vec![0usize; self.bwt.len()];
        for (start, positions) in self.ranks.values() {
            for (rank, &row) in positions.iter().enumerate() {
                lf[row] = start + rank;
            }
        }
        let mut counts = BTreeMap::<[u64; 3], u64>::new();
        for (&last, (start, positions)) in &self.ranks {
            if last < 4 || last % 2 != 0 {
                continue;
            }
            for row in *start..start + positions.len() {
                let gap = self.bwt[row];
                let first = self.bwt[lf[row]];
                if gap < 3 || gap % 2 != 1 || first < 4 || first % 2 != 0 {
                    continue;
                }
                let pattern = [first, gap, last];
                validate_pattern(&pattern)?;
                let first_row = lf[lf[row]];
                let weight = self.weights[first_row + 1] - self.weights[first_row];
                if weight == 0 {
                    continue;
                }
                let key: [u64; 3] = canonical(&pattern)
                    .try_into()
                    .expect("canonical pair retains three tokens");
                let total = counts.entry(key).or_default();
                *total = total
                    .checked_add(weight)
                    .ok_or_else(|| invalid("count overflow"))?;
            }
        }
        Ok(counts)
    }

    pub fn symbols(&self) -> usize {
        self.bwt.len()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    fn walk(nodes: &[i32], gap: u64) -> Vec<u64> {
        encode_walk(
            &nodes
                .iter()
                .enumerate()
                .map(|(i, &n)| (n, i as u64 * gap))
                .collect::<Vec<_>>(),
        )
        .unwrap()
    }
    fn oracle(records: &BTreeMap<Vec<u64>, u64>, f: &[u64]) -> u64 {
        records
            .iter()
            .map(|(r, w)| {
                if r.len() < f.len() {
                    0
                } else {
                    r.windows(f.len()).filter(|s| *s == f).count() as u64 * w
                }
            })
            .sum()
    }
    #[test]
    fn observed_pairs_match_weighted_records_and_public_orbit_counts() {
        let records = BTreeMap::from([
            (walk(&[1, 2, 1, 2], 7), 3),
            (walk(&[-2, -1], 7), 5),
            (walk(&[1, 2], 8), 2),
            (walk(&[1, -1], 7), 4),
            (walk(&[4], 7), 11),
            (walk(&[5], 7), 13),
        ]);
        let mut expected = BTreeMap::<[u64; 3], u64>::new();
        for (record, &weight) in &records {
            for pair in record.windows(3).step_by(2) {
                *expected
                    .entry(canonical(pair).try_into().unwrap())
                    .or_default() += weight;
            }
        }
        let bwt = WeightedBwt::build(&records).unwrap();
        let bytes_before =
            bincode::serde::encode_to_vec(&bwt, bincode::config::standard()).unwrap();
        assert_eq!(bwt.observed_pairs().unwrap(), expected);
        let mut canonical_records = BTreeMap::new();
        for (record, &weight) in &records {
            *canonical_records.entry(canonical(record)).or_default() += weight;
        }
        assert_eq!(
            WeightedBwt::build(&canonical_records)
                .unwrap()
                .observed_pairs()
                .unwrap(),
            expected
        );
        for (pair, &count) in &expected {
            assert_eq!(bwt.count(pair).unwrap(), count);
        }
        assert_eq!(
            bytes_before,
            bincode::serde::encode_to_vec(&bwt, bincode::config::standard()).unwrap()
        );
        assert_eq!(bwt.count(&walk(&[1, -1], 7)).unwrap(), 4);
        let (mut restored, used): (WeightedBwt, usize) =
            bincode::serde::decode_from_slice(&bytes_before, bincode::config::standard()).unwrap();
        assert_eq!(used, bytes_before.len());
        assert!(restored.observed_pairs().is_err());
        restored.rebuild().unwrap();
        assert_eq!(restored.observed_pairs().unwrap(), expected);
        assert!(WeightedBwt::build(&BTreeMap::new())
            .unwrap()
            .observed_pairs()
            .unwrap()
            .is_empty());
    }

    #[test]
    fn exhaustive_weighted_substrings_spacing_repeats_boundaries_orientations() {
        let records = BTreeMap::from([
            (walk(&[1, 2, 1, 2], 7), 3),
            (walk(&[-2, -1], 7), 5),
            (walk(&[1, 2], 8), 2),
            (walk(&[1, -1], 7), 4),
        ]);
        let bwt = WeightedBwt::build(&records).unwrap();
        let canonical_records = records.iter().fold(BTreeMap::new(), |mut m, (r, w)| {
            *m.entry(canonical(r)).or_default() += w;
            m
        });
        let canonical_bwt = WeightedBwt::build(&canonical_records).unwrap();
        let expanded: Vec<_> = records
            .iter()
            .flat_map(|(r, &w)| std::iter::repeat_n(r, w as usize))
            .collect();
        for length in 1..=4u32 {
            for encoded in 0..4usize.pow(length) {
                let mut x = encoded;
                let nodes: Vec<i32> = (0..length)
                    .map(|_| {
                        let n = [1, 2, -1, -2][x % 4];
                        x /= 4;
                        n
                    })
                    .collect();
                for gap in [7, 8, 9] {
                    let f = walk(&nodes, gap);
                    let rc = reverse_complement(&f);
                    assert_eq!(bwt.count_oriented(&f).unwrap(), oracle(&records, &f));
                    let uncollapsed = expanded
                        .iter()
                        .map(|r| {
                            if r.len() < f.len() {
                                0
                            } else {
                                r.windows(f.len()).filter(|w| *w == f).count() as u64
                            }
                        })
                        .sum::<u64>();
                    assert_eq!(bwt.count_oriented(&f).unwrap(), uncollapsed);
                    assert_eq!(bwt.count(&f).unwrap(), canonical_bwt.count(&f).unwrap());
                    assert_eq!(
                        bwt.count(&f).unwrap(),
                        oracle(&records, &f) + if rc == f { 0 } else { oracle(&records, &rc) }
                    );
                }
            }
        }
        assert_eq!(bwt.count(&walk(&[1, 2], 7)).unwrap(), 11);
        assert_eq!(bwt.count(&walk(&[1, -1], 7)).unwrap(), 4); // palindrome
        assert!(bwt.count(&[1]).is_err());
        let spaced = encode_walk(&[(1, 100), (2, 105), (3, 112)]).unwrap();
        let interior = encode_walk(&[(2, 0), (3, 7)]).unwrap();
        let rc = reverse_complement(&spaced);
        assert_eq!(rc, encode_walk(&[(-3, 0), (-2, 7), (-1, 12)]).unwrap());
        let index = WeightedBwt::build(&BTreeMap::from([(spaced, 2), (rc, 3)])).unwrap();
        assert_eq!(index.count(&interior).unwrap(), 5);
        assert_eq!(
            index
                .count(&encode_walk(&[(2, 0), (3, 5)]).unwrap())
                .unwrap(),
            0
        );
    }
    #[test]
    fn empty_serialization_and_invalid_structure() {
        let bwt = WeightedBwt::build(&BTreeMap::new()).unwrap();
        assert_eq!(bwt.count(&walk(&[1], 1)).unwrap(), 0);
        let bytes = bincode::serde::encode_to_vec(&bwt, bincode::config::standard()).unwrap();
        let (mut loaded, _): (WeightedBwt, usize) =
            bincode::serde::decode_from_slice(&bytes, bincode::config::standard()).unwrap();
        loaded.rebuild().unwrap();
        loaded.weights[0] = 1;
        assert!(loaded.rebuild().is_err());
    }
}
