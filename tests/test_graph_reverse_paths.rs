//! A `-` interval is reverse-complemented before the graph is built, but its
//! path keeps the forward name `q:0-1000`. Every engine must emit that path
//! so it spells the forward sequence of `q`. Before the fix only POA reversed
//! the P line; seqwish and pggb spelled the reverse complement.

use std::collections::HashMap;
use std::process::Command;

fn dna(len: usize, seed: u64) -> Vec<u8> {
    let mut x = seed;
    (0..len)
        .map(|_| {
            x = x.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            b"ACGT"[(x >> 62) as usize]
        })
        .collect()
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'A',
        })
        .collect()
}

/// Path name to spelled sequence, from P lines.
fn spell_paths(gfa: &str) -> HashMap<String, Vec<u8>> {
    let segs: HashMap<&str, &[u8]> = gfa
        .lines()
        .filter_map(|l| l.strip_prefix("S\t"))
        .map(|l| {
            let mut f = l.split('\t');
            (f.next().unwrap(), f.next().unwrap().as_bytes())
        })
        .collect();
    gfa.lines()
        .filter_map(|l| l.strip_prefix("P\t"))
        .map(|l| {
            let mut f = l.split('\t');
            let name = f.next().unwrap().to_string();
            let mut seq = Vec::new();
            for step in f.next().unwrap().split(',') {
                let (id, orient) = step.split_at(step.len() - 1);
                let s = segs[id];
                if orient == "-" {
                    seq.extend(revcomp(s));
                } else {
                    seq.extend_from_slice(s);
                }
            }
            (name, seq)
        })
        .collect()
}

#[test]
fn reverse_interval_paths_spell_the_forward_sequence() {
    let dir = tempfile::tempdir().unwrap();
    let t = dna(2000, 5);
    let mut q = revcomp(&t[500..1500]);
    q[300] = if q[300] == b'A' { b'C' } else { b'A' };
    std::fs::write(
        dir.path().join("seqs.fa"),
        format!(
            ">t\n{}\n>q\n{}\n",
            String::from_utf8_lossy(&t),
            String::from_utf8_lossy(&q)
        ),
    )
    .unwrap();
    // q[300] is t[1199] on the reverse strand
    std::fs::write(
        dir.path().join("a.paf"),
        "q\t1000\t0\t1000\t-\tt\t2000\t500\t1500\t999\t1000\t60\tcg:Z:699=1X300=\n",
    )
    .unwrap();

    for engine in ["poa", "seqwish", "pggb", "syng-local"] {
        let out = Command::new(env!("CARGO_BIN_EXE_impg"))
            .current_dir(dir.path())
            .args([
                "query", "-a", "a.paf", "-r", "t:500-1500", "-d", "0", "-o", "gfa",
                "--sequence-files", "seqs.fa", "--gfa-engine", engine, "-t", "1",
            ])
            .output()
            .expect("run impg query");
        let stderr = String::from_utf8_lossy(&out.stderr);
        assert!(out.status.success(), "{engine}: {stderr}");
        let paths = spell_paths(&String::from_utf8_lossy(&out.stdout));
        let spelled = |name: &str| {
            paths
                .get(name)
                .unwrap_or_else(|| panic!("{engine}: no path {name}, paths {:?}", paths.keys()))
        };
        assert_eq!(spelled("t:500-1500"), &t[500..1500].to_vec(), "{engine}: t path");
        assert_eq!(spelled("q:0-1000"), &q, "{engine}: q path must spell q forward");
    }
}
