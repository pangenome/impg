use std::collections::BTreeMap;
use std::process::Command;

fn path_spellings(gfa: &str) -> BTreeMap<String, String> {
    let mut seg = std::collections::HashMap::new();
    for l in gfa.lines().filter(|l| l.starts_with("S\t")) {
        let f: Vec<&str> = l.split('\t').collect();
        seg.insert(f[1].to_string(), f[2].to_string());
    }
    let rc = |s: &str| -> String {
        s.bytes()
            .rev()
            .map(|b| match b { b'A' => 'T', b'C' => 'G', b'G' => 'C', b'T' => 'A', x => x as char })
            .collect()
    };
    gfa.lines()
        .filter(|l| l.starts_with("P\t"))
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            let spell = f[2]
                .split(',')
                .map(|s| {
                    let q = &seg[&s[..s.len() - 1]];
                    if s.ends_with('+') { q.clone() } else { rc(q) }
                })
                .collect();
            (f[1].to_string(), spell)
        })
        .collect()
}

/// impg crush skipped the Ygs sort that the :crush stage runs.
#[test]
fn crush_sorts_by_default_and_keeps_paths() {
    let input = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/test_data/crush/nested_bubbles_real.gfa");
    let dir = tempfile::tempdir().unwrap();
    let mut outs = Vec::new();
    for extra in [&[][..], &["--no-sort"][..]] {
        let out = dir.path().join(format!("o{}.gfa", extra.len()));
        let run = Command::new(env!("CARGO_BIN_EXE_impg"))
            .args(["crush", "-g", input, "-o"])
            .arg(&out)
            .args(extra)
            .args(["-v", "1"])
            .output()
            .expect("run impg crush");
        let stderr = String::from_utf8_lossy(&run.stderr).to_string();
        assert!(run.status.success(), "{stderr}");
        assert_eq!(stderr.contains("crush: sorting with Ygs"), extra.is_empty(), "{stderr}");
        outs.push(std::fs::read_to_string(&out).unwrap());
    }
    let before = path_spellings(&std::fs::read_to_string(input).unwrap());
    assert_eq!(path_spellings(&outs[0]), before);
    assert_eq!(path_spellings(&outs[1]), before);
}
