use std::collections::HashSet;
use std::process::Command;

/// A hairpin edge (2+ -> 2-) on a node that gfaffix splits for a shared
/// prefix: gfaffix 0.2.1 dropped the edge, so path p1 used a missing link.
#[test]
fn gfaffix_keeps_hairpin_edges_of_split_nodes() {
    let dir = tempfile::tempdir().unwrap();
    let input = dir.path().join("hairpin.gfa");
    let output = dir.path().join("out.gfa");
    std::fs::write(
        &input,
        "H\tVN:Z:1.0\nS\t1\tACG\nS\t2\tTTA\nS\t3\tTTC\n\
         L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t2\t-\t0M\n\
         P\tp1\t1+,2+,2-\t*\nP\tp2\t1+,3+\t*\n",
    )
    .unwrap();
    let status = Command::new(env!("CARGO_BIN_EXE_gfaffix"))
        .arg(&input)
        .arg("-o")
        .arg(&output)
        .stdout(std::process::Stdio::null())
        .status()
        .expect("run gfaffix");
    assert!(status.success());
    let gfa = std::fs::read_to_string(&output).unwrap();
    let flip = |o: &str| if o == "+" { "-" } else { "+" };
    let mut edges = HashSet::new();
    for l in gfa.lines().filter(|l| l.starts_with("L\t")) {
        let f: Vec<&str> = l.split('\t').collect();
        edges.insert((f[1].to_string(), f[2].to_string(), f[3].to_string(), f[4].to_string()));
    }
    for l in gfa.lines().filter(|l| l.starts_with("P\t")) {
        let steps: Vec<(String, String)> = l.split('\t').nth(2).unwrap().split(',')
            .map(|s| (s[..s.len() - 1].to_string(), s[s.len() - 1..].to_string()))
            .collect();
        for w in steps.windows(2) {
            let ((a, oa), (b, ob)) = (&w[0], &w[1]);
            let fwd = (a.clone(), oa.clone(), b.clone(), ob.clone());
            let rev = (b.clone(), flip(ob).to_string(), a.clone(), flip(oa).to_string());
            assert!(edges.contains(&fwd) || edges.contains(&rev), "path step {a}{oa} -> {b}{ob} has no edge:\n{gfa}");
        }
    }
}
