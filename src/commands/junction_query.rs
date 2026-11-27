//! Source-backed junction queries on immutable assembly paths. An exact word and a uniquely
//! anchored same-path adjacency are distinct results. SYNG steps annotate, never censor, the
//! exhaustive search: absence of a syncmer (in particular on short arms) is not evidence of absence.
use crate::agc_index::AgcIndex;
use crate::syng::{existing_syng_sidecar_path, SyncmerParams, SyngIndex};
use aho_corasick::{AhoCorasick, AhoCorasickKind};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use serde_json::{json, Value};
use sha2::{Digest, Sha256};
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::fs::{self, File, OpenOptions};
use std::io::{self, BufRead, BufReader, BufWriter, Read, Write};
use std::os::unix::fs::{DirBuilderExt, MetadataExt, OpenOptionsExt};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};
use tempfile::NamedTempFile;

type Result<T> = io::Result<T>;
fn invalid(s: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, s.into())
}
fn ensure(ok: bool, msg: &str) -> Result<()> {
    if ok {
        Ok(())
    } else {
        Err(invalid(msg))
    }
}
const ROLE: [&str; 8] = [
    "flank1_40",
    "flank2_40",
    "join22",
    "join33",
    "joined40",
    "flank1_96",
    "flank2_96",
    "joined96",
];
const W: usize = 96;

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct Model {
    id: String,
    end1: (String, u64, char),
    end2: (String, u64, char),
    insert: Option<String>,
    #[serde(default)]
    coupled: Vec<[u64; 2]>,
}
struct Arm {
    cut: u64,
    reverse: bool,
    left: bool,
    lo: u64,
    seq: Vec<u8>,
    walk: Option<BTreeMap<i64, i32>>,
}
impl Arm {
    fn at(&self, lo: u64, hi: u64) -> Option<&[u8]> {
        self.seq
            .get(lo.checked_sub(self.lo)? as usize..hi.checked_sub(self.lo)? as usize)
    }
    fn width(&self, w: usize) -> Option<Vec<u8>> {
        let lo = if self.left != self.reverse {
            self.cut.checked_sub(w as u64)?
        } else {
            self.cut
        };
        let s = self.at(lo, lo + w as u64)?;
        clean(s).then(|| if self.reverse { rc(s) } else { s.to_vec() })
    }
}
struct Query {
    m: Model,
    a: Arm,
    b: Arm,
    refset: String,
    words: [Option<(usize, bool)>; 8],
}
fn clean(s: &[u8]) -> bool {
    s.iter().all(|c| matches!(c, b'A' | b'C' | b'G' | b'T'))
}
fn rc(s: &[u8]) -> Vec<u8> {
    s.iter()
        .rev()
        .map(|c| match c {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            _ => b'N',
        })
        .collect()
}
fn axis(name: &str) -> Result<&str> {
    let (a, c) = name
        .rsplit_once('#')
        .ok_or_else(|| invalid("expected PanSN path names"))?;
    ensure(
        !c.is_empty() && a.split('#').count() == 2 && a.split('#').all(|s| !s.is_empty()),
        "invalid PanSN path",
    )?;
    Ok(a)
}
fn private_dir(p: &Path) -> Result<()> {
    if !p.exists() {
        fs::DirBuilder::new().mode(0o700).create(p)?;
    }
    let m = fs::symlink_metadata(p)?;
    ensure(
        m.is_dir() && m.mode() & 0o077 == 0,
        "output directory must be private (0700)",
    )
}
fn digest(p: &Path) -> Result<String> {
    let mut f = File::open(p)?;
    let mut h = Sha256::new();
    let mut b = vec![0u8; 1 << 20];
    loop {
        let n = f.read(&mut b)?;
        if n == 0 {
            break;
        }
        h.update(&b[..n]);
    }
    Ok(format!("{:x}", h.finalize()))
}
fn stamp(p: &Path) -> Result<(u64, u64, u64, i64, i64)> {
    let m = fs::metadata(p)?;
    Ok((m.dev(), m.ino(), m.len(), m.mtime(), m.mtime_nsec()))
}
fn commit(mut t: NamedTempFile, p: &Path) -> Result<()> {
    t.flush()?;
    t.as_file().sync_all()?;
    t.persist(p)?;
    Ok(())
}
fn atom(p: &Path, b: &[u8]) -> Result<()> {
    let mut t = NamedTempFile::new_in(p.parent().unwrap())?;
    t.write_all(b)?;
    commit(t, p)
}
fn hex(b: &[u8]) -> String {
    format!("{:x}", Sha256::digest(b))
}

/// Inspect fully contained signed stored steps in a literal source interval. Unknown legacy
/// starts are explicitly unavailable; ambiguous source bases are never certified as matches.
fn steps(
    ix: &SyngIndex,
    p: usize,
    lo: u64,
    s: &[u8],
    rev: bool,
    cut: u64,
) -> Result<Option<BTreeMap<i64, i32>>> {
    if ix.name_map.path_starts[p].is_none() {
        return Ok(None);
    }
    let k = ix.syncmer_length_bp() as u64;
    let hi = lo + s.len() as u64;
    let mut result = BTreeMap::new();
    if s.len() >= k as usize {
        for (node, pos) in ix.walk_path_range(p, lo, hi)? {
            if pos < lo || pos + k > hi {
                continue;
            }
            let a = (pos - lo) as usize;
            let seq = &s[a..a + k as usize];
            if !clean(seq) {
                continue;
            }
            ensure(
                ix.syncmer_seq(node).eq_ignore_ascii_case(seq),
                "stored signed step disagrees with source assembly bases",
            )?;
            let (offset, id) = if rev {
                (cut as i64 - (pos + k) as i64, -node)
            } else {
                (pos as i64 - cut as i64, node)
            };
            ensure(
                result.insert(offset, id).is_none(),
                "duplicate step coordinate",
            )?;
        }
    }
    Ok(Some(result))
}
struct Source {
    agc: AgcIndex,
    names: Vec<String>,
}
impl Source {
    fn fetch(&self, p: usize, a: u64, b: u64) -> Result<Vec<u8>> {
        let s = self.agc.fetch_sequence(
            &self.names[p],
            a.try_into().map_err(invalid_int)?,
            b.try_into().map_err(invalid_int)?,
        )?;
        ensure(s.len() as u64 == b - a, "short assembly fetch")?;
        Ok(s)
    }
}
fn invalid_int(_: std::num::TryFromIntError) -> io::Error {
    invalid("coordinate exceeds AGC range API")
}
fn arm(
    ix: &SyngIndex,
    src: &Source,
    e: &(String, u64, char),
    left: bool,
    cuts: &[[u64; 2]],
) -> Result<Arm> {
    let p = *ix
        .name_map
        .name_to_path
        .get(&e.0)
        .ok_or_else(|| invalid("model reporting path missing from index"))? as usize;
    let len = ix.name_map.path_to_length[p];
    ensure(
        e.1 <= len && matches!(e.2, 'L' | 'R'),
        "invalid cut or retained side",
    )?;
    let rev = (e.2 == 'R') == left;
    let mut lo = e.1.saturating_sub(W as u64);
    let mut hi = (e.1 + W as u64).min(len);
    for pair in cuts {
        let t = pair[usize::from(!left)];
        ensure(t <= len && t.abs_diff(e.1) <= 10000, "invalid coupled cut")?;
        lo = lo.min(t);
        hi = hi.max(t);
    }
    let seq = src.fetch(p, lo, hi)?;
    let walk = steps(ix, p, lo, &seq, rev, e.1)?;
    Ok(Arm {
        cut: e.1,
        reverse: rev,
        left,
        lo,
        seq,
        walk,
    })
}
fn prepare(
    ix: &SyngIndex,
    src: &Source,
    models: &str,
    sets: &BTreeSet<String>,
) -> Result<(Vec<Query>, Vec<Vec<u8>>)> {
    let mut seen = HashSet::new();
    let mut queries = Vec::new();
    let mut ids = HashMap::new();
    let mut canon = Vec::new();
    for (i, l) in BufReader::new(File::open(models)?).lines().enumerate() {
        let mut m: Model = serde_json::from_str(&l?)
            .map_err(|_| invalid(format!("invalid model JSON at line {}", i + 1)))?;
        ensure(
            !m.id.is_empty()
                && m.id
                    .bytes()
                    .all(|v| v.is_ascii_alphanumeric() || b"_.:-".contains(&v))
                && seen.insert(m.id.clone()),
            "missing, invalid or duplicate model ID",
        )?;
        let group = axis(&m.end1.0)?.to_string();
        ensure(
            axis(&m.end2.0)? == group && sets.contains(&group),
            "select both reporting arms' reference assembly",
        )?;
        if let Some(v) = m.insert.as_mut() {
            v.make_ascii_uppercase();
            ensure(
                clean(v.as_bytes()),
                "insert must be ACGT, empty, or null (unknown)",
            )?;
        }
        let rep = [m.end1.1, m.end2.1];
        if m.coupled.is_empty() {
            m.coupled.push(rep)
        }
        ensure(
            m.coupled.contains(&rep),
            "coupled list lacks representative cut pair",
        )?;
        let a = arm(ix, src, &m.end1, true, &m.coupled)?;
        let b = arm(ix, src, &m.end2, false, &m.coupled)?;
        let mut distinct = HashSet::new();
        for &[x, y] in &m.coupled {
            ensure(distinct.insert([x, y]), "duplicate equivalent pair")?;
            if [x, y] == rep {
                continue;
            }
            ensure(
                m.insert.as_deref() == Some(""),
                "nontrivial equivalent cuts require a known empty insert",
            )?;
            let d = |p: u64, q: u64, r: bool| (p as i64 - q as i64) * if r { -1 } else { 1 };
            ensure(
                d(x, a.cut, a.reverse) == d(y, b.cut, b.reverse),
                "uncoupled shifts",
            )?;
            let aa = a.at(x.min(a.cut), x.max(a.cut)).unwrap();
            let bb = b.at(y.min(b.cut), y.max(b.cut)).unwrap();
            ensure(
                clean(aa)
                    && if a.reverse == b.reverse {
                        aa == bb
                    } else {
                        aa == rc(bb)
                    },
                "microhomology shift changes literal allele",
            )?;
        }
        let join = |w| -> Option<Vec<u8>> {
            let mut s = a.width(w)?;
            s.extend_from_slice(m.insert.as_ref()?.as_bytes());
            s.extend(b.width(w)?);
            Some(s)
        };
        let ws = [
            a.width(40),
            b.width(40),
            join(22),
            join(33),
            join(40),
            a.width(W),
            b.width(W),
            join(W),
        ];
        let words = ws.map(|v| {
            v.map(|s| {
                let rev = rc(&s);
                let flipped = rev < s;
                let key = if flipped { rev } else { s };
                let id = *ids.entry(key.clone()).or_insert_with(|| {
                    canon.push(key);
                    canon.len() - 1
                });
                (id, flipped)
            })
        });
        queries.push(Query {
            m,
            a,
            b,
            words,
            refset: group,
        });
    }
    ensure(!queries.is_empty(), "empty model ledger")?;
    Ok((queries, canon))
}

#[derive(Default, Serialize, Deserialize)]
struct Obs {
    n: u64,
    paths: u64,
    examples: Vec<(usize, u64, bool)>,
}
#[derive(Serialize, Deserialize)]
struct Unit {
    bind: String,
    group: String,
    paths: usize,
    bases: u64,
    checked_steps: usize,
    paths_without_checked_step: usize,
    counts: BTreeMap<usize, Obs>,
}
struct Patterns {
    ac: Option<AhoCorasick>,
    mapping: Vec<(usize, bool)>,
    collect: HashSet<usize>,
}
impl Patterns {
    fn new(canon: &[Vec<u8>], qs: &[Query]) -> Result<Self> {
        let mut strings = Vec::new();
        let mut mapping = Vec::new();
        for (i, w) in canon.iter().enumerate() {
            strings.push(w.clone());
            mapping.push((i, false));
            let rev = rc(w);
            if rev != *w {
                strings.push(rev);
                mapping.push((i, true));
            }
        }
        let ac = if strings.is_empty() {
            None
        } else {
            Some(
                AhoCorasick::builder()
                    .kind(Some(AhoCorasickKind::DFA))
                    .build(strings)
                    .map_err(io::Error::other)?,
            )
        };
        Ok(Self {
            ac,
            mapping,
            collect: qs
                .iter()
                .flat_map(|q| [q.words[4], q.words[7]])
                .flatten()
                .map(|v| v.0)
                .collect(),
        })
    }
}
fn assess(
    ix: &SyngIndex,
    src: &Source,
    pat: &Patterns,
    group: &str,
    paths: &[usize],
    limit: usize,
    bind: &str,
) -> Result<Unit> {
    let mut u = Unit {
        bind: bind.into(),
        group: group.into(),
        paths: 0,
        bases: 0,
        checked_steps: 0,
        paths_without_checked_step: 0,
        counts: BTreeMap::new(),
    };
    for &p in paths {
        let len = ix.name_map.path_to_length[p];
        ensure(
            src.agc.get_sequence_length(&src.names[p])? as u64 == len,
            "SYNG/AGC length mismatch",
        )?;
        let seq = src.fetch(p, 0, len)?;
        let width = seq.len().min(512);
        let lo = (seq.len() - width) / 2;
        let n =
            steps(ix, p, lo as u64, &seq[lo..lo + width], false, lo as u64)?.map_or(0, |v| v.len());
        u.checked_steps += n;
        u.paths_without_checked_step += usize::from(n == 0);
        let mut seen = HashSet::new();
        if let Some(ac) = &pat.ac {
            for m in ac.find_overlapping_iter(&seq) {
                let (id, rev) = pat.mapping[m.pattern().as_usize()];
                let o = u.counts.entry(id).or_default();
                o.n += 1;
                if seen.insert(id) {
                    o.paths += 1
                }
                if o.examples.len() < limit && pat.collect.contains(&id) {
                    o.examples.push((p, m.start() as u64, rev));
                }
            }
        }
        u.paths += 1;
        u.bases += len;
    }
    Ok(u)
}
fn obs<'a>(q: &Query, u: &'a Unit, i: usize) -> Option<&'a Obs> {
    q.words[i].and_then(|(k, _)| u.counts.get(&k))
}
fn write_result(
    ix: &SyngIndex,
    src: &Source,
    qs: &[Query],
    units: &[Unit],
    limit: usize,
    out: &Path,
) -> Result<()> {
    let mut summary = NamedTempFile::new_in(out)?;
    let mut hits = NamedTempFile::new_in(out)?;
    {
        let mut sw = BufWriter::new(summary.as_file_mut());
        let mut hw = BufWriter::new(hits.as_file_mut());
        for q in qs {
            let mut words = BTreeMap::new();
            let mut verified = Vec::new();
            let mut ambiguous = Vec::new();
            let mut shown = 0usize;
            let mut placements = 0u64;
            let count = |u: &Unit, i: usize| obs(q, u, i).map_or(0, |o| o.n);
            let reference = units
                .iter()
                .find(|u| u.group == q.refset)
                .expect("reporting assembly selected");
            for (i, label) in ROLE.iter().enumerate() {
                let mut sums = [0u64; 2];
                let mut by = BTreeMap::new();
                for u in units {
                    if let Some(o) = obs(q, u, i) {
                        sums[0] += o.n;
                        sums[1] += o.paths;
                        by.insert(&u.group, o.n);
                    }
                }
                let status = if q.words[i].is_none() {
                    "not_assessed"
                } else if sums[0] == 0 {
                    "absent_after_exhaustive_search"
                } else {
                    "present"
                };
                let reason = if q.words[i].is_none() {
                    Some(if q.m.insert.is_none() && matches!(i, 2 | 3 | 4 | 7) {
                        "unknown_insert"
                    } else {
                        "arm_edge_or_ambiguous_source_base"
                    })
                } else {
                    None
                };
                words.insert(*label,json!({"status":status,"reason":reason,"occurrences":if reason.is_none(){Some(sums[0])}else{None},"paths_with_hit":sums[1],"sets_with_hit":by.len(),"by_set":by}));
            }
            for u in units {
                let n = count(u, 4);
                if n == 0 {
                    continue;
                }
                placements += n;
                let unique = |i: usize, j: usize| {
                    [u, reference]
                        .iter()
                        .all(|g| count(g, i) == 1 && count(g, j) == 1)
                };
                let (j, width, placed) = if count(u, 7) > 0 && unique(5, 6) {
                    (7, W, true)
                } else {
                    (4, 40, unique(0, 1))
                };
                if placed {
                    verified.push(&u.group)
                } else {
                    ambiguous.push(&u.group)
                }
                let ins = q.m.insert.as_ref().unwrap().len() as u64;
                let width = width as u64;
                for &(pid, lo, canon_rev) in &obs(q, u, j).unwrap().examples {
                    if shown == limit {
                        break;
                    }
                    let rev = q.words[j].unwrap().1 ^ canon_rev;
                    let hi = lo + width * 2 + ins;
                    let (r1, r2, c1, c2) = if rev {
                        ([hi - width, hi], [lo, lo + width], hi - width, lo + width)
                    } else {
                        ([lo, lo + width], [hi - width, hi], lo + width, hi - width)
                    };
                    let target = src.fetch(pid, lo, hi)?;
                    let mut spell = q.a.width(width as usize).unwrap();
                    spell.extend_from_slice(q.m.insert.as_ref().unwrap().as_bytes());
                    spell.extend(q.b.width(width as usize).unwrap());
                    ensure(
                        target == if rev { rc(&spell) } else { spell.clone() },
                        "assembly range does not reproduce scanned join",
                    )?;
                    let strand = if spell == rc(&spell) {
                        "+/-"
                    } else if rev {
                        "-"
                    } else {
                        "+"
                    };
                    let mut graph = Vec::new();
                    for (r, c, arm) in [(r1, c1, &q.a), (r2, c2, &q.b)] {
                        let seq = &target[(r[0] - lo) as usize..(r[1] - lo) as usize];
                        let w = steps(ix, pid, r[0], seq, rev, c)?;
                        let shared = w.as_ref().zip(arm.walk.as_ref()).map(|(obs, rf)| {
                            obs.iter().filter(|(p, id)| rf.get(p) == Some(*id)).count()
                        });
                        graph.push(json!({"literal_checked_stored_steps":w.as_ref().map(|v|v.len()),"signed_steps_at_same_arm_offset_as_reference":shared}));
                    }
                    let coupled: Vec<[i64; 2]> =
                        q.m.coupled
                            .iter()
                            .filter(|r| r[0].abs_diff(q.a.cut) <= width)
                            .map(|r| {
                                let d = (r[0] as i64 - q.a.cut as i64)
                                    * if rev == q.a.reverse { 1 } else { -1 };
                                [c1 as i64 + d, c2 as i64 + d]
                            })
                            .collect();
                    serde_json::to_writer(
                        &mut hw,
                        &json!({"id":q.m.id,"set":u.group,"path":ix.name_map.path_to_name[pid],"strand":strand,"equivalents_outside_witness":q.m.coupled.len()-coupled.len(),"arm1":r1,"arm2":r2,"coupled_target_cuts0":coupled,"verified_arm_bases":width,"unique_flanks_in_reporting_and_target_set":placed,"arm_graph_evidence":graph}),
                    )?;
                    hw.write_all(b"\n")?;
                    shown += 1;
                }
            }
            let status = if !verified.is_empty() {
                "verified_same_path_join"
            } else if !ambiguous.is_empty() {
                "ambiguous_locus"
            } else if q.words[4].is_none() {
                "unassessable"
            } else {
                "no_verified_join"
            };
            let line = json!({"id":q.m.id,"literal":words,"placement":status,"unique_anchor_sets":verified,"ambiguous_only_sets":ambiguous,"coupled_reference_placements":q.m.coupled.len(),"same_path_40_occurrences":placements,"examples_shown":shown,"occurrences_not_shown":placements-shown as u64});
            serde_json::to_writer(&mut sw, &line)?;
            sw.write_all(b"\n")?;
        }
        sw.flush()?;
        hw.flush()?;
    }
    commit(summary, &out.join("summary.jsonl"))?;
    commit(hits, &out.join("hits.jsonl"))?;
    Ok(())
}

pub fn run(
    prefix: &str,
    agc: &str,
    models: &str,
    out: &Path,
    selected: &[String],
    limit: usize,
) -> Result<()> {
    private_dir(out)?;
    let guard = OpenOptions::new()
        .create(true)
        .truncate(false)
        .write(true)
        .mode(0o600)
        .open(out.join("owner.lock"))?;
    guard.try_lock().map_err(io::Error::other)?;
    let mut files = BTreeMap::from([
        ("agc", PathBuf::from(agc)),
        ("models", PathBuf::from(models)),
        ("executable", std::env::current_exe()?),
    ]);
    for k in ["meta", "names", "pstep", "spos", "1gbwt", "1khash"] {
        files.insert(
            k,
            if k.starts_with('1') {
                format!("{prefix}.{k}").into()
            } else {
                existing_syng_sidecar_path(prefix, k).into()
            },
        );
    }
    let stamps = files
        .values()
        .map(|p| stamp(p))
        .collect::<Result<Vec<_>>>()?;
    let mut hashes = BTreeMap::new();
    for (k, p) in &files {
        hashes.insert(k, digest(p)?);
    }
    for (i, line) in BufReader::new(File::open(&files["names"])?)
        .lines()
        .enumerate()
    {
        ensure(
            line?
                .split('\t')
                .next()
                .and_then(|v| v.parse::<usize>().ok())
                == Some(i),
            "noncontiguous path numbers in SYNG names",
        )?;
    }
    let ix = SyngIndex::load(prefix, SyncmerParams::default())?;
    let source = AgcIndex::build_from_files(&[agc.into()])?;
    let mut agc_paths = HashMap::new();
    for name in source.sequence_names() {
        let path = name
            .split_whitespace()
            .next()
            .ok_or_else(|| invalid("empty archive path"))?
            .to_owned();
        ensure(
            agc_paths.insert(path, name).is_none(),
            "nonunique archive path name",
        )?;
    }
    let mut groups = BTreeMap::<String, Vec<usize>>::new();
    let mut fetch = Vec::new();
    let mut indexed = HashSet::new();
    for (i, n) in ix.name_map.path_to_name.iter().enumerate() {
        ensure(indexed.insert(n.as_str()), "duplicate index path")?;
        let group = axis(n)?;
        let source = agc_paths
            .get(n)
            .ok_or_else(|| invalid("index path absent from AGC"))?;
        fetch.push(source.clone());
        if selected.is_empty() || selected.iter().any(|s| s == group) {
            groups.entry(group.to_string()).or_default().push(i);
        }
    }
    ensure(
        !groups.is_empty() && selected.iter().all(|s| groups.contains_key(s)),
        "unknown assembly set",
    )?;
    for n in agc_paths.keys() {
        if axis(n).is_ok_and(|a| groups.contains_key(a)) {
            ensure(
                indexed.contains(n.as_str()),
                "selected assembly has paths missing from graph",
            )?;
        }
    }
    let src = Source {
        agc: source,
        names: fetch,
    };
    let sets = groups.keys().cloned().collect();
    let (qs, canon) = prepare(&ix, &src, models, &sets)?;
    let pat = Patterns::new(&canon, &qs)?;
    let expected_paths = groups.values().map(Vec::len).sum::<usize>();
    let total_bases: u64 = groups
        .values()
        .flatten()
        .map(|&p| ix.name_map.path_to_length[p])
        .sum();
    let cfg = json!({"schema":"impg.junction-query.v1","sha256":hashes,"sets":sets,"syng_metadata_version":ix.metadata_version(),"syncmer_length":ix.syncmer_length_bp(),"smer_length":ix.params.k,"seed":ix.params.seed,"models":qs.len(),"distinct_canonical_words":canon.len(),"paths_expected":expected_paths,"bases_expected":total_bases,"examples_limit":limit});
    let cb = serde_json::to_vec_pretty(&cfg)?;
    let bind = hex(&cb);
    let cfg_path = out.join("config.json");
    if cfg_path.exists() {
        ensure(
            fs::read(cfg_path)? == cb,
            "configuration changed; use a new output directory",
        )?;
    } else {
        ensure(
            fs::read_dir(out)?.count() == 1,
            "nonempty unbound output directory",
        )?;
        atom(&cfg_path, &cb)?;
    }
    let d = out.join("units");
    private_dir(&d)?;
    let done = AtomicUsize::new(0);
    let all: Vec<_> = groups.iter().collect();
    eprintln!(
        "junction-query: {} models, {} words, {} complete sets, {} paths",
        qs.len(),
        canon.len(),
        all.len(),
        expected_paths
    );
    let units = all
        .par_iter()
        .enumerate()
        .map(|(i, &(group, paths))| -> Result<(Unit, String)> {
            let unit = d.join(format!("set-{i:04}.json"));
            let seal = unit.with_extension("sha256");
            let (u, sha) = if seal.exists() {
                let bytes = fs::read(&unit)?;
                let sha = hex(&bytes);
                ensure(fs::read_to_string(seal)? == sha, "checkpoint hash mismatch")?;
                (serde_json::from_slice::<Unit>(&bytes)?, sha)
            } else {
                let u = assess(&ix, &src, &pat, group, paths, limit, &bind)?;
                let bytes = serde_json::to_vec(&u)?;
                let sha = hex(&bytes);
                atom(&unit, &bytes)?;
                atom(&seal, sha.as_bytes())?;
                (u, sha)
            };
            ensure(
                u.bind == bind
                    && u.group == *group
                    && u.paths == paths.len()
                    && u.bases
                        == paths
                            .iter()
                            .map(|&p| ix.name_map.path_to_length[p])
                            .sum::<u64>(),
                "incompatible set checkpoint",
            )?;
            let n = done.fetch_add(1, Ordering::Relaxed) + 1;
            if n == 1 || n % 20 == 0 || n == all.len() {
                eprintln!("junction-query: {n}/{} sets sealed", all.len());
            }
            Ok((u, sha))
        })
        .collect::<Result<Vec<_>>>()?;
    let (units, seals): (Vec<_>, Vec<_>) = units.into_iter().unzip();
    write_result(&ix, &src, &qs, &units, limit, out)?;
    ensure(
        stamps
            == files
                .values()
                .map(|p| stamp(p))
                .collect::<Result<Vec<_>>>()?,
        "input changed during scan",
    )?;
    let receipt = json!({"config_sha256":bind,"set_sha256_in_config_order":seals,"summary_sha256":digest(&out.join("summary.jsonl"))?,"hits_sha256":digest(&out.join("hits.jsonl"))?,"paths_assessed":units.iter().map(|u|u.paths).sum::<usize>(),"bases_assessed":units.iter().map(|u|u.bases).sum::<u64>(),"models":qs.len(),"signed_midpath_steps_compared_to_source":units.iter().map(|u|u.checked_steps).sum::<usize>(),"paths_without_informative_midpoint_step":units.iter().map(|u|u.paths_without_checked_step).sum::<usize>(),"semantics":"All literal positions/strands in the selected AGC paths, not all human alleles. Uniqueness is conditional on arm length and the reporting and target assemblies. Counts are never capped. No patient genotype inferred."});
    atom(
        &out.join("complete.json"),
        &serde_json::to_vec_pretty(&receipt)?,
    )?;
    eprintln!("junction-query: complete");
    Ok(())
}
