// Probe: compare absolute syncmer positions from syng's C iterator against
// the walk_path positions in query_region_with_anchors. For a sequence where
// the first syncmer is NOT at position 0, the delta is exposed.
use impg::syng::{SyncmerParams, SyngIndex};
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

/// Return the absolute positions of every syncmer in `seq`, as reported
/// directly by syng's C iterator.
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
        let sit = syng_ffi::syncmerIterator(sh, seq_buf.as_mut_ptr() as *mut i8, seq_len as i32);
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
    let params = SyncmerParams {
        k: 8,
        w: 55,
        seed: 7,
    };

    // Try many seeds to find one whose first syncmer is NOT at position 0.
    // Syng's first syncmer lands anywhere in [0, w+k) depending on the seq.
    for seed in 0..50u64 {
        let seq = make_seq(2000, seed);
        let raw = raw_syncmer_positions(&seq, params);
        if raw.is_empty() {
            continue;
        }
        let first = raw[0];
        if first == 0 {
            // Skip the trivial case where the bug happens to cancel.
            continue;
        }

        // Build an index with this sequence and compare what
        // query_region_with_anchors reports for the shared syncmers of the
        // sequence against itself.
        let idx = SyngIndex::build(params, vec![("seq".to_string(), seq.clone())].into_iter());

        let anchored = idx
            .query_region_with_anchors("seq", 0, seq.len() as u64, 0)
            .unwrap();
        // Find the forward-strand self-hit: query_pos == target_pos.
        let self_hit = anchored
            .iter()
            .find(|h| h.genome == "seq" && h.strand == '+');
        if let Some(h) = self_hit {
            let first_anchor = h.anchors[0];
            println!(
                "seed={} raw_first_pos={} walk_path_first_q_pos={} walk_path_first_t_pos={}",
                seed, first, first_anchor.query_pos, first_anchor.target_pos,
            );
            println!(
                "  delta (raw - walk_path) = {}",
                first as i64 - first_anchor.query_pos as i64,
            );

            // The walk_path answer should equal the raw answer. Under the
            // bug, walk_path returns 0 for the first anchor but raw says `first`.
            let walk_first = first_anchor.query_pos as i32;
            let status = if walk_first == first { "OK" } else { "BUG" };
            println!("  status: {}", status);

            // Also emit the ORIGINAL-coord offsets for every shared syncmer
            // under the two conventions, to make the systematic shift obvious.
            let mut rows: Vec<(u64, u64)> = h
                .anchors
                .iter()
                .map(|a| (a.query_pos, a.target_pos))
                .collect();
            rows.sort();
            for (i, &(q, t)) in rows.iter().enumerate().take(6) {
                let raw_pos = raw.get(i).copied().unwrap_or(-1);
                println!("    syncmer #{i}: raw={raw_pos} walk_q={q} walk_t={t}");
            }

            // Only print the first offending seed to keep output short.
            break;
        }
    }
}
