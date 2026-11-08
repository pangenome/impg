use std::process::Command;

fn impg(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_impg"))
        .args(args)
        .output()
        .expect("run impg")
}

/// sweepga options that impg never reads were listed in help and
/// silently ignored; now they are hidden and rejected.
#[test]
fn unused_sweepga_flags_are_hidden_and_rejected() {
    let help = String::from_utf8_lossy(&impg(&["graph", "--help"]).stdout).to_string();
    for flag in [
        "--self ",
        "--zstd ",
        "--scaffolds-only",
        "--max-pairs",
        "--group-by",
    ] {
        assert!(!help.contains(flag), "{flag} still in impg graph --help");
    }
    assert!(
        help.contains("--num-mappings"),
        "used options must stay in help"
    );
    let out = impg(&[
        "graph",
        "--sequence-files",
        "x.fa",
        "-g",
        "out.gfa",
        "--self",
    ]);
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert_eq!(out.status.code(), Some(2), "{stderr}");
    assert!(
        stderr.contains("--self is a sweepga option that impg does not use"),
        "{stderr}"
    );
}

/// --wfmash and --fastga were ignored too (impg always used --aligner);
/// they now set the aligner, as in sweepga.
#[test]
fn aligner_shorthands_set_the_aligner() {
    let dir = tempfile::tempdir().unwrap();
    let fasta = dir.path().join("x.fa");
    let seq: String = (0..2000)
        .map(|i| b"ACGT"[(i * 7 + i / 3) % 4] as char)
        .collect();
    std::fs::write(&fasta, format!(">a#1#c\n{seq}\n>b#1#c\n{seq}\n")).unwrap();
    let jobs = |flag: &str| {
        let out_dir = dir.path().join(flag.trim_start_matches('-'));
        let fasta = fasta.to_str().unwrap();
        let out = impg(
            &[
                "align",
                "--format",
                "joblist",
                "--sequence-files",
                fasta,
                "-o",
            ]
            .iter()
            .copied()
            .chain([out_dir.to_str().unwrap(), flag])
            .collect::<Vec<_>>(),
        );
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        std::fs::read_to_string(out_dir.join("align_jobs.txt")).unwrap()
    };
    assert!(
        jobs("--wfmash").starts_with("wfmash "),
        "{}",
        jobs("--wfmash")
    );
    assert!(
        jobs("--fastga").starts_with("FastGA "),
        "{}",
        jobs("--fastga")
    );
    let both = impg(&[
        "align",
        "--sequence-files",
        "x.fa",
        "-o",
        "o",
        "--wfmash",
        "--fastga",
    ]);
    assert_ne!(both.status.code(), Some(0));
    assert!(
        String::from_utf8_lossy(&both.stderr).contains("Cannot specify both --wfmash and --fastga")
    );
}
