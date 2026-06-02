//! Parallel construction helpers for syng indices.
//!
//! The first implemented stage builds a deterministic global syncmer
//! dictionary in parallel. Existing serial GBWT path construction can then
//! replay sequences through that frozen dictionary.

use rayon::prelude::*;
use std::io;

use crate::syng::{packed_syncmer_word_len, PackedSyncmer, SyncmerParams};
use crate::syng_ffi;

fn base_code(base: u8) -> u8 {
    match base {
        b'a' | b'A' | 0 => 0,
        b'c' | b'C' | 1 => 1,
        b'g' | b'G' | 2 => 2,
        b't' | b'T' | 3 => 3,
        _ => 0,
    }
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

/// Convert a sequence to syng's numeric DNA encoding.
pub fn numeric_sequence(seq: &[u8]) -> Vec<u8> {
    seq.iter().map(|&b| base_code(b)).collect()
}

/// Pack a full syncmer in the canonical orientation expected by KmerHash.
pub fn pack_canonical_syncmer(syncmer: &[u8]) -> PackedSyncmer {
    let word_len = (syncmer.len() + 31) / 32;
    let mut bytes = vec![0u8; word_len * 8];
    let use_forward = forward_is_canonical(syncmer);

    for i in 0..syncmer.len() {
        let code = if use_forward {
            base_code(syncmer[i])
        } else {
            complement_code(base_code(syncmer[syncmer.len() - 1 - i]))
        };
        bytes[i / 4] |= code << (2 * (i % 4));
    }

    bytes
        .chunks_exact(8)
        .map(|chunk| {
            let mut word = [0u8; 8];
            word.copy_from_slice(chunk);
            u64::from_le_bytes(word)
        })
        .collect()
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

    let mut seq_buf = numeric_sequence(seq);
    seq_buf.push(0);

    let mut out = Vec::new();
    unsafe {
        let seqhash =
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32);
        if seqhash.is_null() {
            return Err(io::Error::other("seqhashCreate returned null"));
        }

        let sit =
            syng_ffi::syncmerIterator(seqhash, seq_buf.as_mut_ptr() as *mut i8, seq.len() as i32);
        if sit.is_null() {
            syng_ffi::impg_seqhashDestroy(seqhash);
            return Err(io::Error::other("syncmerIterator returned null"));
        }

        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(sit, std::ptr::null_mut(), &mut pos, std::ptr::null_mut()) {
            let start = pos as usize;
            if start + syncmer_len > seq.len() {
                syng_ffi::impg_seqhashIteratorDestroy(sit);
                syng_ffi::impg_seqhashDestroy(seqhash);
                return Err(io::Error::new(
                    io::ErrorKind::InvalidData,
                    format!(
                        "syncmer iterator returned out-of-range position {} for sequence length {}",
                        start,
                        seq.len()
                    ),
                ));
            }
            out.push(pack_canonical_syncmer(&seq_buf[start..start + syncmer_len]));
        }

        syng_ffi::impg_seqhashIteratorDestroy(sit);
        syng_ffi::impg_seqhashDestroy(seqhash);
    }

    Ok(out)
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
    let expected_words = packed_syncmer_word_len(params);
    let syncmers = sequences
        .par_iter()
        .map(|(_, seq)| extract_packed_syncmers(params, seq))
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
    }
}
