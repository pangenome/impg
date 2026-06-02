//! Isolate startCount increment behavior: call syngBWTpathStartNew twice
//! at the same start node and verify j_last increments from 0 to 1.

use impg::syng::{SyncmerParams, SyngIndex};
use impg::syng_ffi;
use std::ffi::CString;
use std::path::Path;

const SYNG_SIMPLE_SIDE_COUNT_LIMIT: u32 = 65_000;

static SYNG_LOCK: std::sync::LazyLock<std::sync::Mutex<()>> =
    std::sync::LazyLock::new(|| std::sync::Mutex::new(()));

unsafe fn walk_first_step(
    gbwt: *mut syng_ffi::SyngBWT,
    start_node: i32,
    start_count: u32,
) -> (i32, u32) {
    let sbp = syng_ffi::syngBWTpathStartOld(gbwt, start_node, start_count);
    let mut next_node = 0i32;
    let mut offset = 0u32;
    assert!(syng_ffi::syngBWTpathNext(sbp, &mut next_node, &mut offset));
    syng_ffi::syngBWTpathDestroy(sbp);
    (next_node, offset)
}

unsafe fn roundtrip_gbwt(gbwt: *mut syng_ffi::SyngBWT, path: &Path) -> *mut syng_ffi::SyngBWT {
    let schema_text = syng_ffi::syng_schema_text();
    let schema = syng_ffi::oneSchemaCreateFromText(schema_text.as_ptr());
    assert!(!schema.is_null());

    let path_c = CString::new(path.to_str().unwrap()).unwrap();
    let type_c = CString::new("gbwt").unwrap();

    let out = syng_ffi::oneFileOpenWriteNew(path_c.as_ptr(), schema, type_c.as_ptr(), true, 1);
    assert!(!out.is_null());
    syng_ffi::syngBWTwrite(out, gbwt);
    syng_ffi::oneFileClose(out);

    let input = syng_ffi::oneFileOpenRead(path_c.as_ptr(), schema, type_c.as_ptr(), 1);
    assert!(!input.is_null());
    let loaded = syng_ffi::syngBWTread(input);
    syng_ffi::oneFileClose(input);
    syng_ffi::oneSchemaDestroy(schema);

    assert!(!loaded.is_null());
    loaded
}

#[test]
fn test_path_offset_above_u16_is_preserved() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    unsafe {
        syng_ffi::impg_syng_suppress_debug();
        let gbwt = syng_ffi::syngBWTcreate(63, 0);
        let sbp = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        syng_ffi::syngBWTpathAdd(sbp, 2, 70_000);
        syng_ffi::syngBWTpathFinish(sbp);

        assert_eq!(walk_first_step(gbwt, 1, 0), (2, 70_000));

        let dir = tempfile::Builder::new()
            .prefix("impg_test_syng_large_offset_")
            .tempdir()
            .unwrap();
        let loaded = roundtrip_gbwt(gbwt, &dir.path().join("idx.1gbwt"));
        syng_ffi::syngBWTdestroy(gbwt);

        assert_eq!(walk_first_step(loaded, 1, 0), (2, 70_000));
        syng_ffi::syngBWTdestroy(loaded);
    }
}

#[test]
fn test_one_edge_rskip_side_survives_load() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    unsafe {
        syng_ffi::impg_syng_suppress_debug();
        let gbwt = syng_ffi::syngBWTcreate(63, 0);
        let n_paths = SYNG_SIMPLE_SIDE_COUNT_LIMIT + 1;
        for _ in 0..n_paths {
            let sbp = syng_ffi::syngBWTpathStartNew(gbwt, 1);
            syng_ffi::syngBWTpathAdd(sbp, 2, 10);
            syng_ffi::syngBWTpathFinish(sbp);
        }

        let dir = tempfile::Builder::new()
            .prefix("impg_test_syng_one_edge_rskip_load_")
            .tempdir()
            .unwrap();
        let loaded = roundtrip_gbwt(gbwt, &dir.path().join("idx.1gbwt"));
        syng_ffi::syngBWTdestroy(gbwt);

        assert_eq!(walk_first_step(loaded, 1, 0), (2, 10));
        assert_eq!(walk_first_step(loaded, 1, n_paths - 1), (2, 10));
        syng_ffi::syngBWTdestroy(loaded);
    }
}

#[test]
fn test_start_count_increments_on_second_path() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    unsafe {
        syng_ffi::impg_syng_suppress_debug();

        let syncmer_len = 63; // default
        let gbwt = syng_ffi::syngBWTcreate(syncmer_len, 0);
        assert!(!gbwt.is_null());

        // First path starting at node 5
        let sbp1 = syng_ffi::syngBWTpathStartNew(gbwt, 5);
        let j1 = (*sbp1).j_last;
        eprintln!("First path at node 5: j_last = {}", j1);

        // Add a single follow node and finish (making it a complete path)
        syng_ffi::syngBWTpathAdd(sbp1, 6, 10);
        syng_ffi::syngBWTpathFinish(sbp1);

        // Second path at the SAME start node
        let sbp2 = syng_ffi::syngBWTpathStartNew(gbwt, 5);
        let j2 = (*sbp2).j_last;
        eprintln!("Second path at node 5: j_last = {}", j2);

        syng_ffi::syngBWTpathAdd(sbp2, 6, 10);
        syng_ffi::syngBWTpathFinish(sbp2);

        // Third at same node
        let sbp3 = syng_ffi::syngBWTpathStartNew(gbwt, 5);
        let j3 = (*sbp3).j_last;
        eprintln!("Third path at node 5: j_last = {}", j3);
        syng_ffi::syngBWTpathAdd(sbp3, 6, 10);
        syng_ffi::syngBWTpathFinish(sbp3);

        syng_ffi::syngBWTdestroy(gbwt);

        assert_eq!(j1, 0, "First path should have j_last=0");
        assert_eq!(
            j2, 1,
            "Second path should have j_last=1 (startCount incrementing)"
        );
        assert_eq!(j3, 2, "Third path should have j_last=2");
    }
}

/// Start at node 5 (not 1), longer paths.
#[test]
fn test_start_count_node5_long() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    unsafe {
        syng_ffi::impg_syng_suppress_debug();
        let gbwt = syng_ffi::syngBWTcreate(63, 0);
        let sbp_a = syng_ffi::syngBWTpathStartNew(gbwt, 5);
        let j_a = (*sbp_a).j_last;
        for k in 6..=55i32 {
            syng_ffi::syngBWTpathAdd(sbp_a, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_a);

        let sbp_b = syng_ffi::syngBWTpathStartNew(gbwt, 5);
        let j_b = (*sbp_b).j_last;
        for k in 6..=55i32 {
            syng_ffi::syngBWTpathAdd(sbp_b, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_b);
        syng_ffi::syngBWTdestroy(gbwt);

        eprintln!("node5 long: j_a={} j_b={}", j_a, j_b);
        assert_eq!(j_b, 1, "should increment");
    }
}

/// Start at node 1, SHORT paths.
#[test]
fn test_start_count_node1_short() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());
    unsafe {
        syng_ffi::impg_syng_suppress_debug();
        let gbwt = syng_ffi::syngBWTcreate(63, 0);
        let sbp_a = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_a = (*sbp_a).j_last;
        syng_ffi::syngBWTpathAdd(sbp_a, 2, 10);
        syng_ffi::syngBWTpathFinish(sbp_a);

        let sbp_b = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_b = (*sbp_b).j_last;
        syng_ffi::syngBWTpathAdd(sbp_b, 2, 10);
        syng_ffi::syngBWTpathFinish(sbp_b);
        syng_ffi::syngBWTdestroy(gbwt);

        eprintln!("node1 short: j_a={} j_b={}", j_a, j_b);
        assert_eq!(j_b, 1, "should increment");
    }
}

/// Two long identical forward paths, NO RC path in between.
#[test]
fn test_start_count_long_paths_no_rc() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    unsafe {
        syng_ffi::impg_syng_suppress_debug();
        let gbwt = syng_ffi::syngBWTcreate(63, 0);

        let sbp_a = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_a = (*sbp_a).j_last;
        for k in 2..=50i32 {
            syng_ffi::syngBWTpathAdd(sbp_a, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_a);

        let sbp_b = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_b = (*sbp_b).j_last;
        for k in 2..=50i32 {
            syng_ffi::syngBWTpathAdd(sbp_b, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_b);

        syng_ffi::syngBWTdestroy(gbwt);

        eprintln!("j_a={} j_b={}", j_a, j_b);
        assert_eq!(j_a, 0);
        assert_eq!(j_b, 1, "without RC: is j_b 1 or 0?");
    }
}

/// Like the basic test, but with many pathAdd calls per path (mimicking a long sequence)
/// and an RC path between the two forwards (mimicking the build loop).
#[test]
fn test_start_count_with_long_paths_and_rc() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    unsafe {
        syng_ffi::impg_syng_suppress_debug();

        let gbwt = syng_ffi::syngBWTcreate(63, 0);

        // Path A forward: nodes 1 -> 2 -> 3 -> ... -> 50
        let sbp_a = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_a = (*sbp_a).j_last;
        for k in 2..=50i32 {
            syng_ffi::syngBWTpathAdd(sbp_a, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_a);

        // Path A reverse complement: negative nodes -50 -> -49 -> ... -> -1
        let sbp_a_rc = syng_ffi::syngBWTpathStartNew(gbwt, -50);
        for k in (1..50i32).rev() {
            syng_ffi::syngBWTpathAdd(sbp_a_rc, -k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_a_rc);

        // Path B forward (identical to A): 1 -> 2 -> ... -> 50
        let sbp_b = syng_ffi::syngBWTpathStartNew(gbwt, 1);
        let j_b = (*sbp_b).j_last;
        eprintln!("Path A j_last = {}, Path B j_last = {}", j_a, j_b);
        for k in 2..=50i32 {
            syng_ffi::syngBWTpathAdd(sbp_b, k, 10);
        }
        syng_ffi::syngBWTpathFinish(sbp_b);

        syng_ffi::syngBWTdestroy(gbwt);

        assert_eq!(j_a, 0);
        assert_eq!(
            j_b, 1,
            "Path B should have j_last=1 (not 0) — bug repro if this fails"
        );
    }
}

fn make_seq(len: usize, seed: u8) -> Vec<u8> {
    let bases = [b'A', b'C', b'G', b'T'];
    let mut seq = Vec::with_capacity(len);
    let mut state: u32 = seed as u32;
    for _ in 0..len {
        state = state.wrapping_mul(1103515245).wrapping_add(12345);
        seq.push(bases[((state >> 16) % 4) as usize]);
    }
    seq
}

#[test]
fn test_identical_sequences_get_distinct_start_counts() {
    let _guard = SYNG_LOCK.lock().unwrap_or_else(|e| e.into_inner());

    let params = SyncmerParams::default();
    let seq = make_seq(1500, 42);

    // Two sequences with identical content but different names — mimics
    // what happens in yeast235 where AAA#0#chrIII and SGDref#0#chrIII are
    // byte-identical but named differently.
    let sequences = vec![
        ("seqA".to_string(), seq.clone()),
        ("seqB".to_string(), seq.clone()),
    ];

    let mut index = SyngIndex::build(params, sequences.into_iter());

    let start_a = index.name_map.path_starts[0]
        .as_ref()
        .expect("seqA missing");
    let start_b = index.name_map.path_starts[1]
        .as_ref()
        .expect("seqB missing");

    eprintln!(
        "seqA: start_node={} start_count={}",
        start_a.start_node, start_a.start_count
    );
    eprintln!(
        "seqB: start_node={} start_count={}",
        start_b.start_node, start_b.start_count
    );

    assert_eq!(
        start_a.start_node, start_b.start_node,
        "Identical sequences should share a start_node"
    );

    // The second sequence's start_count must be different from the first —
    // otherwise both paths collide at index 0 in the GBWT startCount table
    // and subsequent queries via syngBWTpathStartOld will fail.
    assert_ne!(
        start_a.start_count, start_b.start_count,
        "Two paths starting at the same node must have distinct start_count indices"
    );

    // Also round-trip: save, load, then try to walk each path via query_region.
    let dir = tempfile::Builder::new()
        .prefix("impg_test_ident_seqs_")
        .tempdir()
        .unwrap();
    let prefix = dir.path().join("idx");
    let prefix_str = prefix.to_str().unwrap();
    index.save(prefix_str).expect("save failed");
    let loaded = SyngIndex::load(prefix_str, SyncmerParams::default()).expect("load failed");

    // Query seqB — if the start_count was wrong, this will crash with
    // "syngBWTpathStartOld ... >= startCount" in the C layer.
    let intervals = loaded
        .query_region("seqB", 0, 1000, 0)
        .expect("query_region for seqB should succeed");
    eprintln!("query_region(seqB) returned {} intervals", intervals.len());
    assert!(!intervals.is_empty(), "Should find at least self-hit");

    drop(loaded);
    drop(index);
}
