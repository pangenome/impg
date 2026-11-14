//! Parallel construction helpers for syng indices.
//!
//! The first implemented stage builds a deterministic global syncmer
//! dictionary in parallel. Existing serial GBWT path construction can then
//! replay sequences through that frozen dictionary.

use rayon::prelude::*;
use std::io;

use crate::syng::{
    encode_query_base, packed_syncmer_word_len, PackedSyncmer, SyncmerParams, SyncmerWindowValidity,
};
use crate::syng_ffi;

fn base_code(base: u8) -> u8 {
    encode_query_base(base).0
}

fn complement_code(code: u8) -> u8 {
    3 - code
}

fn forward_is_canonical(syncmer: &[u8]) -> bool {
    if syncmer.is_empty() {
        return true;
    }
    let mut left = 0usize;
    let mut right = syncmer.len() - 1;
    while left <= right {
        let fwd = base_code(syncmer[left]);
        let rev = complement_code(base_code(syncmer[right]));
        if fwd != rev {
            return fwd < rev;
        }
        if left == right {
            break;
        }
        left += 1;
        right -= 1;
    }
    true
}

/// Convert to C-safe numeric DNA. Invalid bytes become temporary zeroes; callers
/// extracting seeds must separately validate the full original window.
pub fn numeric_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| base_code(b)).collect()
}

/// Pack a full syncmer in the canonical orientation expected by KmerHash.
pub fn pack_canonical_syncmer(syncmer: &[u8]) -> PackedSyncmer {
    let word_len = syncmer.len().div_ceil(32);
    let mut words = vec![0u64; word_len];
    let len = syncmer.len();
    let use_forward = forward_is_canonical(syncmer);

    for i in 0..len {
        let code = if use_forward {
            base_code(syncmer[i])
        } else {
            complement_code(base_code(syncmer[len - 1 - i]))
        };
        words[i >> 5] |= (code as u64) << (2 * (i & 31));
    }

    words
}

struct SeqhashGuard(*mut syng_ffi::Seqhash);

impl SeqhashGuard {
    fn new(params: SyncmerParams) -> io::Result<Self> {
        let seqhash = unsafe {
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32)
        };
        if seqhash.is_null() {
            return Err(io::Error::other("seqhashCreate returned null"));
        }
        Ok(Self(seqhash))
    }
}

impl Drop for SeqhashGuard {
    fn drop(&mut self) {
        if !self.0.is_null() {
            unsafe {
                syng_ffi::impg_seqhashDestroy(self.0);
            }
        }
    }
}

fn extract_packed_syncmers_with(
    seqhash: *mut syng_ffi::Seqhash,
    syncmer_len: usize,
    seq: &[u8],
    seq_buf: &mut Vec<u8>,
) -> io::Result<Vec<PackedSyncmer>> {
    if seq.len() < syncmer_len {
        return Ok(Vec::new());
    }
    if seq.len() > std::os::raw::c_int::MAX as usize {
        return Err(io::Error::new(
            io::ErrorKind::InvalidInput,
            format!(
                "sequence length {} exceeds C int maximum ({})",
                seq.len(),
                std::os::raw::c_int::MAX
            ),
        ));
    }

    seq_buf.clear();
    seq_buf.reserve(seq.len() + 1);
    seq_buf.extend(seq.iter().map(|&b| base_code(b)));
    seq_buf.push(0);

    let mut out = Vec::new();
    let mut validity = SyncmerWindowValidity::new(seq);
    unsafe {
        let sit = syng_ffi::syncmerIterator(
            seqhash,
            seq_buf.as_mut_ptr() as *mut std::os::raw::c_char,
            seq.len() as std::os::raw::c_int,
        );
        if sit.is_null() {
            return Err(io::Error::other("syncmerIterator returned null"));
        }

        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(sit, std::ptr::null_mut(), &mut pos, std::ptr::null_mut()) {
            if pos < 0 {
                syng_ffi::impg_seqhashIteratorDestroy(sit);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!("syncmer iterator returned negative position {}", pos),
                ));
            }
            let start = pos as usize;
            if start + syncmer_len > seq.len() {
                syng_ffi::impg_seqhashIteratorDestroy(sit);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "syncmer iterator returned out-of-range position {} for sequence length {}",
                        start,
                        seq.len()
                    ),
                ));
            }
            if !validity.is_valid(start, syncmer_len) {
                continue;
            }
            out.push(pack_canonical_syncmer(&seq_buf[start..start + syncmer_len]));
        }

        syng_ffi::impg_seqhashIteratorDestroy(sit);
    }

    Ok(out)
}

/// Extract packed canonical full-syncmers from a sequence.
pub fn extract_packed_syncmers(
    params: SyncmerParams,
    seq: &[u8],
) -> io::Result<Vec<PackedSyncmer>> {
    let syncmer_len = (params.w + params.k) as usize;
    if seq.len() < syncmer_len {
        return Ok(Vec::new());
    }
    let guard = SeqhashGuard::new(params)?;
    let mut seq_buf = Vec::new();
    extract_packed_syncmers_with(guard.0, syncmer_len, seq, &mut seq_buf)
}

pub fn sort_dedup_packed_syncmers(mut syncmers: Vec<PackedSyncmer>) -> Vec<PackedSyncmer> {
    syncmers.sort_unstable();
    syncmers.dedup();
    syncmers
}

/// Build a deterministic global packed-syncmer dictionary in parallel.
pub fn build_packed_syncmer_dictionary(
    params: SyncmerParams,
    sequences: &[(String, Vec<u8>)],
) -> io::Result<Vec<PackedSyncmer>> {
    let syncmer_len = (params.w + params.k) as usize;
    let expected_words = packed_syncmer_word_len(params);
    let syncmers = sequences
        .par_iter()
        .map_init(
            || (SeqhashGuard::new(params), Vec::new()),
            |(guard_res, seq_buf), (_, seq)| match guard_res {
                Ok(guard) => extract_packed_syncmers_with(guard.0, syncmer_len, seq, seq_buf),
                Err(err) => Err(io::Error::new(err.kind(), err.to_string())),
            },
        )
        .try_reduce(Vec::new, |mut acc, mut part| {
            acc.append(&mut part);
            Ok(acc)
        })?;

    let dictionary = sort_dedup_packed_syncmers(syncmers);
    if let Some((i, packed)) = dictionary
        .iter()
        .enumerate()
        .find(|(_, packed)| packed.len() != expected_words)
    {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!(
                "dictionary entry {} has {} words, expected {}",
                i,
                packed.len(),
                expected_words
            ),
        ));
    }
    Ok(dictionary)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn canonical_packing_is_orientation_invariant() {
        let seq = b"ACGTACGTACGTA";
        let rc = b"TACGTACGTACGT";
        assert_eq!(
            pack_canonical_syncmer(seq),
            pack_canonical_syncmer(rc),
            "forward and reverse-complement syncmers should share packed canonical form"
        );
        // Multi-word (>32 bp) canonical packing invariant across word boundary.
        let long_seq = b"ACGGTTCAGATTACAGGCATCGATCGGATTACCATGACCAGTA";
        let mut long_rc = Vec::with_capacity(long_seq.len());
        for &b in long_seq.iter().rev() {
            long_rc.push(match b {
                b'A' => b'T',
                b'C' => b'G',
                b'G' => b'C',
                b'T' => b'A',
                other => other,
            });
        }
        assert_ne!(long_seq.as_slice(), long_rc.as_slice());
        assert_eq!(
            pack_canonical_syncmer(long_seq),
            pack_canonical_syncmer(&long_rc),
        );
        // Forward is canonical: A < T.
        let mut expected = vec![0u64; 2];
        for (i, &b) in long_seq.iter().enumerate() {
            let code = match b {
                b'A' => 0u64,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            expected[i / 32] |= code << (2 * (i % 32));
        }
        assert_ne!(expected[1], 0);
        assert_eq!(pack_canonical_syncmer(long_seq), expected);
    }

    #[test]
    fn parallel_dictionary_matches_sequential_extraction() {
        let params = SyncmerParams {
            w: 11,
            k: 5,
            seed: 7,
        };
        let seqs = vec![
            (
                "s1".to_string(),
                b"ACGTACGTACGTACGTTGCATGCATGCATGCA".to_vec(),
            ),
            (
                "s2".to_string(),
                b"TGCATGCATGCATGCAACGTACGTACGTACGT".to_vec(),
            ),
        ];
        let parallel = build_packed_syncmer_dictionary(params, &seqs).unwrap();
        let mut manual = Vec::new();
        for (_, seq) in &seqs {
            manual.extend(extract_packed_syncmers(params, seq).unwrap());
        }
        let manual = sort_dedup_packed_syncmers(manual);
        assert_eq!(parallel, manual);
    }
}
