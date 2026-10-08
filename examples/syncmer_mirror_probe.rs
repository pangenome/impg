//! One-off probe: is syng's closed-syncmer selection mirror-symmetric?
//!
//! Extracts the raw selected 63-mer start positions of a synthetic sequence
//! S and of rc(S), mirrors the rc positions (q -> L - 63 - q), and reports
//! the symmetric difference. Also reports the content-only prediction: a
//! 63-mer K (and its rc) should be selected iff the minimal canonical 8-mer
//! hash over K's 8-mer offsets 0..55 is at offset 0 or 54 — computed
//! independently in Rust from the same syng hash function.
//!
//! Usage: syncmer_mirror_probe [len] [seed]

use impg::syng::SyncmerParams;
use impg::syng_ffi;

fn make_seq(len: usize, seed: u64) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut s = Vec::with_capacity(len);
    let mut state: u64 = seed;
    for _ in 0..len {
        state = state
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        s.push(bases[((state >> 33) & 3) as usize]);
    }
    s
}

fn reverse_complement(seq: &[u8]) -> Vec<u8> {
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

/// Absolute positions of every syncmer selected by syng's C iterator.
fn raw_syncmer_positions(seq: &[u8], params: SyncmerParams) -> Vec<i32> {
    let seq_len = seq.len();
    let mut seq_buf: Vec<u8> = seq
        .iter()
        .map(|&b| match b {
            b'a' | b'A' | 0 => 0u8,
            b'c' | b'C' | 1 => 1u8,
            b'g' | b'G' | 2 => 2u8,
            b't' | b'T' | 3 => 3u8,
            _ => 0u8,
        })
        .collect();
    seq_buf.push(0);

    let mut positions = Vec::new();
    unsafe {
        let sh =
            syng_ffi::impg_seqhashCreateSafe(params.k as i32, params.w as i32, params.seed as i32);
        let sit = syng_ffi::syncmerIterator(
            sh,
            seq_buf.as_mut_ptr() as *mut std::os::raw::c_char,
            seq_len as i32,
        );
        let mut pos: i32 = 0;
        while syng_ffi::syncmerNext(sit, std::ptr::null_mut(), &mut pos, std::ptr::null_mut()) {
            positions.push(pos);
        }
        syng_ffi::impg_seqhashIteratorDestroy(sit);
        syng_ffi::impg_seqhashDestroy(sh);
    }
    positions
}

fn main() {
    let len: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(100_000);
    let seed: u64 = std::env::args()
        .nth(2)
        .and_then(|s| s.parse().ok())
        .unwrap_or(20260925);
    let params = SyncmerParams {
        k: 8,
        w: 55,
        seed: 7,
    };
    let syncmer_len = (params.k + params.w) as usize;

    for round in 0..4u64 {
        let seq = make_seq(len, seed + round);
        let rc = reverse_complement(&seq);
        let fwd = raw_syncmer_positions(&seq, params);
        let rc_sel = raw_syncmer_positions(&rc, params);
        let l = seq.len() as i64;
        // Mirror rc position q to forward position L - 63 - q.
        let mirrored: Vec<i64> = rc_sel
            .iter()
            .map(|&q| l - syncmer_len as i64 - q as i64)
            .collect();
        let fwd_set: std::collections::BTreeSet<i64> = fwd.iter().map(|&p| p as i64).collect();
        let rc_set: std::collections::BTreeSet<i64> = mirrored.into_iter().collect();
        let only_fwd: Vec<i64> = fwd_set.difference(&rc_set).copied().collect();
        let only_rc: Vec<i64> = rc_set.difference(&fwd_set).copied().collect();
        println!(
            "[round {round}] len {len} fwd {} rc {} only_fwd {} only_rc {}",
            fwd.len(),
            rc_sel.len(),
            only_fwd.len(),
            only_rc.len()
        );
        if !only_fwd.is_empty() || !only_rc.is_empty() {
            println!("  only_fwd (first 20): {:?}", &only_fwd[..only_fwd.len().min(20)]);
            println!("  only_rc  (first 20): {:?}", &only_rc[..only_rc.len().min(20)]);
        }
        // Edge-adjacency of the disagreements (boundary vs interior).
        let interior_only_fwd = only_fwd
            .iter()
            .filter(|&&p| p > 70 && p < l - 133)
            .count();
        let interior_only_rc = only_rc
            .iter()
            .filter(|&&p| p > 70 && p < l - 133)
            .count();
        println!("  interior-only fwd {} rc {}", interior_only_fwd, interior_only_rc);
    }

    // Context-sensitivity: which bases OUTSIDE the k-mer [p, p+k) influence
    // the selection of position p? (Determines whether read-edge anchors have
    // truncated-context decisions.)
    {
        let params8 = SyncmerParams { k: 8, w: 55, seed: 7 };
        let seq = make_seq(3000, 99);
        let base_positions: std::collections::BTreeSet<i32> =
            raw_syncmer_positions(&seq, params8).into_iter().collect();
        let mut influential_offsets: std::collections::BTreeSet<i32> = Default::default();
        for &p in base_positions.iter() {
            if p < 80 || p + 63 > 2920 {
                continue;
            }
            for offset in -8i32..=70i32 {
                if offset >= 0 && offset < 63 {
                    continue; // inside the k-mer: content itself
                }
                let index = (p + offset) as usize;
                let original = seq[index];
                for &replacement in [b'A', b'C', b'G', b'T'].iter() {
                    if replacement == original {
                        continue;
                    }
                    let mut mutated = seq.clone();
                    mutated[index] = replacement;
                    let mutated_positions: std::collections::BTreeSet<i32> =
                        raw_syncmer_positions(&mutated, params8).into_iter().collect();
                    if !mutated_positions.contains(&p) {
                        influential_offsets.insert(offset);
                    }
                }
            }
        }
        println!("influential offsets outside the k-mer: {:?}", influential_offsets);
    }

    // Edge behavior: selection of a 150bp window must equal the selection of
    // the same bases in the full sequence (content-only must hold at read
    // edges, or read anchors would disagree with panel steps).
    {
        let params8 = SyncmerParams { k: 8, w: 55, seed: 7 };
        let seq = make_seq(5000, 4242);
        let full: std::collections::BTreeSet<i32> =
            raw_syncmer_positions(&seq, params8).into_iter().collect();
        let mut windows_checked = 0usize;
        let mut windows_agree = 0usize;
        let mut disagreement_examples: Vec<(usize, Vec<i32>, Vec<i32>)> = Vec::new();
        for start in (0..5000 - 150).step_by(37) {
            let window = &seq[start..start + 150];
            let window_sel: Vec<i32> = raw_syncmer_positions(window, params8);
            let expected: Vec<i32> = full
                .range(start as i32..(start + 150) as i32)
                .map(|&p| p - start as i32)
                .filter(|&p| p >= 0 && (p as usize) + 63 <= 150)
                .collect();
            windows_checked += 1;
            if window_sel == expected {
                windows_agree += 1;
            } else if disagreement_examples.len() < 8 {
                disagreement_examples.push((start, window_sel.clone(), expected.clone()));
            }
        }
        println!(
            "150bp window selection vs full-sequence selection: {windows_agree}/{windows_checked} agree"
        );
        for (start, got, want) in disagreement_examples {
            println!("  window @{start}: got {got:?} expected {want:?}");
        }
    }
}
