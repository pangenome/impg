//! Raw FFI bindings to syng's C library.
//!
//! These declarations match the public API in syng's C headers.
//! Only consumed by `crate::syng` — do not use directly.

#![allow(non_camel_case_types, dead_code)]

use std::os::raw::{c_char, c_int, c_void};

// --- Syng type aliases matching utils.h ---
pub type I32 = i32;
pub type I64 = i64;
pub type U32 = u32;
pub type U64 = u64;
pub type U8 = u8;

// --- Opaque pointer types ---

/// Opaque: syng Array (dynamic array from array.h)
pub type Array = *mut c_void;

/// Opaque: syng Hash (hash table from hash.h)
pub type Hash = *mut c_void;

/// Opaque: syng Rskip (run-length skip list from rskip.h)
pub type Rskip = *mut c_void;

/// Opaque: ONElib OneFile handle
pub type OneFile = *mut c_void;

/// Opaque: ONElib OneSchema handle
pub type OneSchema = *mut c_void;

// --- Structs ---

/// Path metadata stored in the GBWT (syng.h: SyngPath).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct SyngPath {
    pub file: U32,
    pub path: U32,
    pub length: U64,
}

/// GBWT structure (syng.h: SyngBWT).
/// Fields are laid out to match C, but most are accessed only through C functions.
#[repr(C)]
pub struct SyngBWT {
    pub fixed_len: c_int,
    pub node: Array,
    pub status: Array,
    pub length: Array,
    pub path: Array,
    pub loc: *mut c_void,
    pub start_hash: Hash,
    pub start_hash_count: Array,
}

/// Path traversal state (syng.h: SyngBWTpath).
#[repr(C)]
pub struct SyngBWTpath {
    pub sb: *mut SyngBWT,
    pub last_node: I32,
    pub this_node: I32,
    pub last_off: U32,
    pub j_last: U32,
    pub j_max: U32,
}

/// Seqhash parameters (seqhash.h: Seqhash).
#[repr(C)]
pub struct Seqhash {
    pub seed: c_int,
    pub k: c_int,
    pub w: c_int,
    pub mask: U64,
    pub shift1: c_int,
    pub shift2: c_int,
    pub factor1: U64,
    pub factor2: U64,
    pub pattern_rc: [U64; 4],
}

/// Seqhash iterator state (seqhash.h: SeqhashIterator).
#[repr(C)]
pub struct SeqhashIterator {
    pub sh: *mut Seqhash,
    pub s: *mut c_char,
    pub s_end: *mut c_char,
    pub h: U64,
    pub h_rc: U64,
    pub hash: *mut U64,
    pub is_forward: *mut bool,
    pub base: c_int,
    pub i_start: c_int,
    pub i_min: c_int,
    pub min: U64,
    pub is_done: bool,
}

/// SeqPack for 2-bit DNA packing (seqio.h).
#[repr(C)]
pub struct SeqPack {
    pub unconv: [c_char; 5],
    pub unconvc: [c_char; 5],
    pub byte_expand: [U32; 256],
    pub byte_expandc: [U32; 256],
}

/// KmerHash table (kmerhash.h).
#[repr(C)]
pub struct KmerHash {
    pub len: c_int,
    pub dim: c_int,
    pub table: *mut I64,
    pub mask: I64,
    pub max: I64,
    pub plen: U64,
    pub psize: U64,
    pub pack: *mut U64,
    pub seqbuf: *mut c_char,
    pub finds: U64,
    pub deltas: U64,
    pub seq_pack: *mut SeqPack,
}

/// Syncmer parameters (syncmerset.h: SyncmerParams).
#[repr(C)]
#[derive(Debug, Clone, Copy)]
pub struct CSyncmerParams {
    pub w: c_int,
    pub k: c_int,
    pub seed: c_int,
}

/// SyncmerSet (syncmerset.h).
#[repr(C)]
pub struct SyncmerSet {
    pub params: CSyncmerParams,
    pub kh: *mut KmerHash,
    pub count: Array,
    pub tot_count: I64,
    pub this_count: Array,
    pub max_count: Array,
    pub loc: Array,
}

extern "C" {
    // --- SyngBWT (syng.h / syngbwt3.c) ---

    pub fn syngBWTcreate(fixed_len: c_int, max: I64) -> *mut SyngBWT;
    pub fn syngBWTdestroy(sb: *mut SyngBWT);
    pub fn syngBWTwrite(of: *mut c_void, sb: *mut SyngBWT);
    pub fn syngBWTread(of: *mut c_void) -> *mut SyngBWT;

    pub fn syngBWTpathStartNew(sb: *mut SyngBWT, start_node: I32) -> *mut SyngBWTpath;
    pub fn syngBWTpathAdd(sbp: *mut SyngBWTpath, next_node: I32, offset: U32);
    pub fn syngBWTpathFinish(sbp: *mut SyngBWTpath);

    pub fn syngBWTpathStartOld(sb: *mut SyngBWT, start_node: I32, count: U32) -> *mut SyngBWTpath;
    pub fn syngBWTpathNext(sbp: *mut SyngBWTpath, next_node: *mut I32, next_pos: *mut U32) -> bool;

    pub fn syngBWTmatchStart(sb: *mut SyngBWT, start_node: I32, high: *mut U32)
        -> *mut SyngBWTpath;
    pub fn syngBWTmatchNext(
        sbp: *mut SyngBWTpath,
        next_node: I32,
        next_off: U32,
        low: *mut U32,
        high: *mut U32,
    ) -> bool;
    pub fn syngBWTincomingRank(
        sb: *mut SyngBWT,
        node: I32,
        prev_node: I32,
        prev_off: U32,
        rank: U32,
        abs_rank: *mut U32,
    ) -> bool;
    pub fn syngBWTadvanceRank(
        sb: *mut SyngBWT,
        node: I32,
        abs_rank: U32,
        next_node: *mut I32,
        next_off: *mut U32,
        next_abs_rank: *mut U32,
    ) -> bool;

    pub fn syngBWTpathDestroy(sbp: *mut SyngBWTpath);
    pub fn syngBWTstat(sb: *mut SyngBWT);

    pub fn syngBWTlocFind(
        sb: *mut SyngBWT,
        loc: I64,
        file: *mut I64,
        path: *mut I64,
        offset: *mut I64,
    ) -> bool;

    // --- KmerHash (kmerhash.h) ---

    pub fn kmerHashCreate(initial_size: U64, len: c_int) -> *mut KmerHash;
    pub fn kmerHashDestroy(kh: *mut KmerHash);

    pub fn kmerHashAdd(kh: *mut KmerHash, dna: *mut c_char, index: *mut I64) -> bool;
    pub fn kmerHashFind(kh: *mut KmerHash, dna: *mut c_char, index: *mut I64) -> bool;
    pub fn kmerHashFindThreadSafe(
        kh: *mut KmerHash,
        dna: *mut c_char,
        index: *mut I64,
        buf: *mut U64,
    ) -> bool;

    pub fn kmerHashAddPacked(kh: *mut KmerHash, u: *mut U64, index: *mut I64) -> bool;
    pub fn kmerHashFindPacked(kh: *mut KmerHash, u: *mut U64, index: *mut I64) -> bool;

    pub fn kmerHashSeq(kh: *mut KmerHash, i: I64, buf: *mut c_char) -> *mut c_char;

    pub fn kmerHashWriteOneFile(kh: *mut KmerHash, of: *mut c_void) -> bool;
    pub fn kmerHashReadOneFile(of: *mut c_void) -> *mut KmerHash;

    // --- Seqhash (seqhash.h) ---

    pub fn seqhashCreate(k: c_int, w: c_int, seed: c_int) -> *mut Seqhash;
    /// Thread-safe variant — locks around srandom/random calls.
    pub fn impg_seqhashCreateSafe(k: c_int, w: c_int, seed: c_int) -> *mut Seqhash;

    pub fn seqhashWrite(sh: *mut Seqhash, f: *mut c_void);
    pub fn seqhashRead(f: *mut c_void) -> *mut Seqhash;

    pub fn seqhashIterator(sh: *mut Seqhash, s: *mut c_char, len: c_int) -> *mut SeqhashIterator;
    pub fn seqhashNext(
        si: *mut SeqhashIterator,
        kmer: *mut U64,
        pos: *mut c_int,
        is_f: *mut bool,
    ) -> bool;

    pub fn syncmerIterator(sh: *mut Seqhash, s: *mut c_char, len: c_int) -> *mut SeqhashIterator;
    pub fn syncmerNext(
        si: *mut SeqhashIterator,
        kmer: *mut U64,
        pos: *mut c_int,
        is_f: *mut bool,
    ) -> bool;

    // --- Rskip (rskip.h) ---

    pub fn rsNsym(rs: Rskip) -> c_int;
    pub fn rsLength(rs: Rskip) -> c_int;
    pub fn rsCount(rs: Rskip, symbol: I32) -> c_int;
    pub fn rsFind(rs: Rskip, k: U32, symbol: *mut I32) -> c_int;
    pub fn rsFindSyng(rs: Rskip, k: U32, symbol: *mut I32, offset: *mut U32) -> c_int;
    pub fn rsDestroy(rs: Rskip) -> c_int;

    // --- SyncmerSet (syncmerset.h) ---

    pub fn syncmerSetCreate(params: CSyncmerParams, initial_size: U64) -> *mut SyncmerSet;
    pub fn syncmerSetDestroy(sms: *mut SyncmerSet);
    pub fn syncmerSetWrite(sms: *mut SyncmerSet, of: *mut c_void) -> bool;
    pub fn syncmerSetRead(filename: *mut c_char) -> *mut SyncmerSet;

    pub fn syncmerParamsDefault() -> CSyncmerParams;

    // --- ONElib (ONElib.h) ---

    pub fn oneSchemaCreateFromText(text: *const c_char) -> *mut c_void;
    pub fn oneSchemaDestroy(schema: *mut c_void);

    pub fn oneFileOpenRead(
        path: *const c_char,
        schema: *mut c_void,
        file_type: *const c_char,
        nthreads: c_int,
    ) -> *mut c_void;
    pub fn oneFileOpenWriteNew(
        path: *const c_char,
        schema: *mut c_void,
        file_type: *const c_char,
        is_binary: bool,
        nthreads: c_int,
    ) -> *mut c_void;
    pub fn oneFileClose(of: *mut c_void);

    // --- utils.h ---

    pub fn timeUpdate(f: *mut c_void);
    pub fn timeTotal(f: *mut c_void);

    // --- impg_syng_helpers.c (wrappers for static-inline functions) ---

    pub fn impg_seqhashIteratorDestroy(si: *mut SeqhashIterator);
    pub fn impg_seqhashDestroy(sh: *mut Seqhash);

    /// Suppress syngBWT debug output by setting pathCount to a non-zero value.
    pub fn impg_syng_suppress_debug();
}

/// The syng ONEcode schema text, replicated from syng.h.
/// Used to open .1gbwt and .1khash files via ONElib.
pub static SYNG_SCHEMA_TEXT: &[u8] = include_bytes!("../vendor/syng/syng.h");

// We can't directly include the static char* from the header. Instead, we'll
// reconstruct it in Rust when needed or pass it from the C side. For now, we
// provide a function that returns the schema text as a CString.
pub fn syng_schema_text() -> std::ffi::CString {
    // This is the schema text from syng.h syngSchemaText
    std::ffi::CString::new(
        "1 3 def 1 0               schema for syng\n\
         .\n\
         P 3 seq                   SEQUENCE\n\
         S 4 path                  contains P (path) objects = syncmer sequences\n\
         S 3 gfa                   sequence graph - contains V (vertex) objects, probably with E lines\n\
         S 4 gbwt                  gbwt: a gfa with B, C, Z lines\n\
         .\n\
         D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n\
         .\n\
         O S 1 3 DNA               sequence of the node\n\
         .\n\
         O V 1 3 INT               graph node (vertex): length\n\
         D K 1 3 INT               coverage count of node\n\
         D E 3 3 INT 3 INT 3 INT   edge +: adjacent node (- if reversed), offset, count\n\
         D e 3 3 INT 3 INT 3 INT   edge -: adjacent node (- if reversed), offset, count\n\
         D B 1 8 INT_LIST          GBWT +: list of node indices from opposite-signed E lines\n\
         D b 1 8 INT_LIST          GBWT -: list of node indices from opposite-signed E lines\n\
         D C 1 8 INT_LIST          GBWT +: list of run-length counts\n\
         D c 1 8 INT_LIST          GBWT -: list of run-length counts\n\
         .\n\
         O P 3 3 INT 3 INT 3 INT   path: length in bp, source file number, sequence number in file\n\
         D I 1 6 STRING            identifier of path - used to return mapping information\n\
         D Z 4 3 INT 3 INT 3 INT 3 INT   GBWT path: starting node, pos, count, then length in nodes\n\
         D z 1 8 INT_LIST          alternative explicit list of node ids (-ve if reversed)\n\
         D o 1 8 INT_LIST          if z, then offsets of the nodes from start of sequence, 1:1 with z\n\
         D X 1 3 DNA               prefix before first node - required to fully reconstruct\n\
         D Y 1 3 DNA               suffix after last node - required to fully reconstruct\n\
         .\n\
         P 5 khash                 KMER HASH\n\
         S 7 syncset               SYNCMER SET\n\
         D h 3 3 INT 3 INT 3 INT   k, w, seed for the seqhash: for syncs k = |smer|, w+k = |syncmer|\n\
         O t 3 3 INT 3 INT 3 INT   max, len, dim for KmerHash table\n\
         D S 1 3 DNA               packed sequences aligned to 64-bit boundaries\n\
         D L 1 8 INT_LIST          locations in the table\n\
         D C 1 8 INT_LIST          kmer counts\n\
         D M 1 6 STRING            maximum count in any input - (1..127)\n\
         .\n\
         P 3 map                   MAP\n\
         O S 2 3 INT 3 INT         query sequence: index in source file (1-based) length\n\
         D I 1 6 STRING            identifier from source file (if requested)\n\
         D F 1 4 CHAR              filter: Z zero-length, Q quality, G poly-G (Illumina bad read)\n\
         D M 3 3 INT 3 INT 3 INT   mem: start, end (0-based), count\n\
         D X 2 3 INT 3 DNA         missing syncmer not found in graph: start coordinate, sequence\n\
         D U 3 3 INT 3 INT 3 INT   unique mapping: file, path, offset (negative offset = reverse)"
    ).expect("CString::new failed for syng schema text")
}
