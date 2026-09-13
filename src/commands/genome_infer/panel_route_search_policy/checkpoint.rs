//! Versioned, immutable-parent checkpoints; checksums detect corruption, not malicious forgery.
use super::*;
use std::fs::{self, File, OpenOptions};
use std::io::{Read, Write};

#[derive(Clone, Debug, Serialize, Deserialize, PartialEq, Eq)]
#[serde(deny_unknown_fields)]
pub(super) struct Bindings {
    pub policy: BTreeMap<String, String>,
    pub backend: String,
    pub graph: String,
    pub sample: String,
    pub depth_bits: u64,
    pub background_bits: u64,
    pub tie_bits: u64,
    pub max_feature_terms: usize,
    pub cache_terms: usize,
    pub count_policy: String,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Parent {
    pub(super) directory: PathBuf,
    pub(super) checkpoint: String,
    pub(super) ledger: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Seal {
    pub(super) version: u32,
    pub(super) checkpoint: String,
    pub(super) ledger: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub(super) struct Checkpoint {
    pub(super) version: u32,
    pub(super) bindings: Bindings,
    pub(super) parent: Option<Parent>,
    pub(super) budgets: Vec<Limits>,
    pub(super) state: State,
    pub(super) transition: Option<transition::Receipt>,
}
pub(super) fn hash_file(path: &Path) -> io::Result<String> {
    let mut f = File::open(path)?;
    let mut hash = 0xcbf29ce484222325u64;
    let mut buf = [0; 65536];
    let mut size = 0u64;
    loop {
        let n = f.read(&mut buf)?;
        if n == 0 {
            break;
        }
        size = size
            .checked_add(n as u64)
            .ok_or_else(|| invalid("file size overflow"))?;
        for b in &buf[..n] {
            hash = (hash ^ *b as u64).wrapping_mul(0x100000001b3);
        }
    }
    Ok(format!("fnv1a64:{size}:{hash:016x}"))
}
// A resumed ledger may append records, but must preserve its parent's bytes.
fn verify_ledger_prefix(child: &Path, parent: &Path) -> io::Result<()> {
    let mut child = File::open(child)?;
    let mut parent = File::open(parent)?;
    let mut expected = [0; 65536];
    let mut actual = [0; 65536];
    loop {
        let n = parent.read(&mut expected)?;
        if n == 0 {
            return Ok(());
        }
        child.read_exact(&mut actual[..n])?;
        ensure(
            actual[..n] == expected[..n],
            "inherited ledger prefix mismatch",
        )?;
    }
}
fn verify(directory: &Path) -> io::Result<(Checkpoint, Seal)> {
    let seal: Seal = genome::read_json(&directory.join("checkpoint-seal.json"))?;
    ensure(
        seal.version == 1
            && seal.checkpoint == hash_file(&directory.join("checkpoint.json"))?
            && seal.ledger == hash_file(&directory.join("evaluations.jsonl"))?,
        "checkpoint/ledger seal mismatch",
    )?;
    let checkpoint: Checkpoint = serde_json::from_reader(std::io::BufReader::new(File::open(
        directory.join("checkpoint.json"),
    )?))
    .map_err(io::Error::other)?;
    ensure(checkpoint.version == 2, "unsupported checkpoint version")?;
    validation::state(&checkpoint.state)?;
    validation::ledger(&checkpoint.state, &directory.join("evaluations.jsonl"))?;
    if checkpoint.transition.is_some() {
        transition::verify_boundary(&checkpoint, &seal.ledger)?;
    }
    Ok((checkpoint, seal))
}
pub(super) fn load(
    directory: &Path,
    bindings: &Bindings,
    limits: &Limits,
    extend: bool,
) -> io::Result<(State, Vec<Limits>, Parent)> {
    let directory = fs::canonicalize(directory)?;
    let (checkpoint, seal) = verify(&directory)?;
    ensure(
        &checkpoint.bindings == bindings,
        "resume policy/backend/graph/sample/parameter mismatch",
    )?;
    let old = checkpoint
        .budgets
        .last()
        .ok_or_else(|| invalid("missing cumulative budgets"))?;
    ensure(
        limits.max_work >= old.max_work
            && limits.max_evaluations >= old.max_evaluations
            && limits.max_state_bytes >= old.max_state_bytes
            && limits.max_optima >= old.max_optima,
        "resume budgets cannot decrease",
    )?;
    ensure(
        limits == old || extend,
        "budget extension requires --extend-budgets",
    )?;
    let mut ancestor = if checkpoint.transition.is_some() {
        None
    } else {
        checkpoint.parent.clone()
    };
    let mut seen = BTreeSet::from([directory.clone()]);
    let mut child_ledger = directory.join("evaluations.jsonl");
    while let Some(parent) = ancestor {
        ensure(
            seen.insert(fs::canonicalize(&parent.directory)?),
            "checkpoint ancestry cycle",
        )?;
        let (cp, s) = verify(&parent.directory)?;
        ensure(
            s.checkpoint == parent.checkpoint
                && s.ledger == parent.ledger
                && cp.bindings == *bindings,
            "immutable parent checkpoint/ledger mismatch",
        )?;
        let parent_ledger = parent.directory.join("evaluations.jsonl");
        verify_ledger_prefix(&child_ledger, &parent_ledger)?;
        child_ledger = parent_ledger;
        ancestor = if cp.transition.is_some() {
            None
        } else {
            cp.parent
        };
    }
    let extended = limits != old;
    let mut budgets = checkpoint.budgets;
    if extended {
        budgets.push(limits.clone());
    }
    Ok((
        checkpoint.state,
        budgets,
        Parent {
            directory,
            checkpoint: seal.checkpoint,
            ledger: seal.ledger,
        },
    ))
}
pub(super) fn ledger(out: &Path, parent: Option<&Parent>) -> io::Result<File> {
    let mut ledger = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(out.join("evaluations.jsonl"))?;
    if let Some(p) = parent {
        io::copy(
            &mut File::open(p.directory.join("evaluations.jsonl"))?,
            &mut ledger,
        )?;
        // Detect changed parent bytes during the copy before any new operation.
        ledger.flush()?;
        ensure(
            hash_file(&out.join("evaluations.jsonl"))? == p.ledger,
            "parent ledger changed during copy",
        )?;
    }
    Ok(ledger)
}
pub(super) fn save(
    out: &Path,
    bindings: Bindings,
    parent: Option<Parent>,
    budgets: Vec<Limits>,
    state: State,
    ledger: &mut File,
    transition: Option<transition::Receipt>,
) -> io::Result<()> {
    state.validate_links()?;
    ledger.flush()?;
    ledger.sync_all()?;
    write_checkpoint_json(
        &out.join("checkpoint.json"),
        &Checkpoint {
            version: 2,
            bindings,
            parent,
            budgets,
            state,
            transition,
        },
    )?;
    genome::write_json(
        &out.join("checkpoint-seal.json"),
        &Seal {
            version: 1,
            checkpoint: hash_file(&out.join("checkpoint.json"))?,
            ledger: hash_file(&out.join("evaluations.jsonl"))?,
        },
    )?;
    #[cfg(unix)]
    File::open(out)?.sync_all()?;
    Ok(())
}

// Same pretty-JSON bytes (no trailing newline) and atomic installation as the
// backend helper, without retaining a second full checkpoint-sized byte buffer.
fn write_checkpoint_json(path: &Path, value: &impl Serialize) -> io::Result<()> {
    let incomplete = path.with_extension("incomplete");
    let file = OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&incomplete)?;
    let mut writer = std::io::BufWriter::new(file);
    serde_json::to_writer_pretty(&mut writer, value).map_err(io::Error::other)?;
    writer.flush()?;
    writer.get_ref().sync_all()?;
    fs::rename(incomplete, path)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cumulative_ledger_preserves_exact_parent_prefix() {
        let d = tempfile::tempdir().unwrap();
        let parent = d.path().join("parent.jsonl");
        let child = d.path().join("child.jsonl");
        let original = vec![b'x'; 65537];
        fs::write(&parent, &original).unwrap();
        fs::write(&child, &original).unwrap();
        verify_ledger_prefix(&child, &parent).unwrap();
        let mut appended = original.clone();
        appended.extend_from_slice(b"\nnext record\n");
        fs::write(&child, &appended).unwrap();
        verify_ledger_prefix(&child, &parent).unwrap();
        appended[65536] = b'y';
        fs::write(&child, &appended).unwrap();
        assert!(verify_ledger_prefix(&child, &parent).is_err());
        fs::write(&child, &original[..65536]).unwrap();
        assert!(verify_ledger_prefix(&child, &parent).is_err());
        fs::write(&parent, b"").unwrap();
        verify_ledger_prefix(&child, &parent).unwrap();
    }
    #[test]
    fn streamed_checkpoint_bytes_match_existing_writer_without_newline() {
        let d = tempfile::tempdir().unwrap();
        let mut s = State::new(2);
        s.spawn(
            Context {
                family: 0,
                slot: 0,
                completed: vec![],
                segments: vec![],
                source: 0,
                cut: 0,
                reverse: false,
            },
            Op::SourceScan {
                base: 0,
                permutation: Permutation::new(u64::MAX),
            },
        )
        .unwrap();
        s.seed_focus(0);
        let value = (&s, u64::MAX, "escaped\ncheckpoint", Some(123u64));
        let old = d.path().join("old.json");
        let new = d.path().join("new.json");
        genome::write_json(&old, &value).unwrap();
        write_checkpoint_json(&new, &value).unwrap();
        let before = fs::read(old).unwrap();
        assert_eq!(before, fs::read(&new).unwrap());
        assert_ne!(before.last(), Some(&b'\n'));
        assert!(!new.with_extension("incomplete").exists());
    }
    #[test]
    fn interrupted_or_failed_stream_never_installs_a_checkpoint() {
        struct Failed;
        impl Serialize for Failed {
            fn serialize<S: serde::Serializer>(&self, serializer: S) -> Result<S::Ok, S::Error> {
                use serde::ser::SerializeStruct;
                let mut s = serializer.serialize_struct("Failed", 2)?;
                s.serialize_field("partial", &"x".repeat(16384))?;
                Err(serde::ser::Error::custom(
                    "injected checkpoint serialization failure",
                ))
            }
        }
        let d = tempfile::tempdir().unwrap();
        let path = d.path().join("checkpoint.json");
        assert!(write_checkpoint_json(&path, &Failed).is_err());
        assert!(!path.exists());
        assert!(!d.path().join("checkpoint-seal.json").exists());
        let incomplete = path.with_extension("incomplete");
        assert!(incomplete.exists());
        let bytes = fs::read(&incomplete).unwrap();
        assert!(!bytes.is_empty());
        assert!(write_checkpoint_json(&path, &State::new(1)).is_err());
        assert_eq!(bytes, fs::read(&incomplete).unwrap());
        assert!(!path.exists());
        // A failed attempt must not touch an already preserved destination either.
        fs::write(&path, b"preserved ancestor").unwrap();
        assert!(write_checkpoint_json(&path, &State::new(1)).is_err());
        assert_eq!(fs::read(path).unwrap(), b"preserved ancestor");
    }
}
