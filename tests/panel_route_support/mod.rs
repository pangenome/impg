//! Only explicitly requested ignored gates may consult private oracle configuration.
pub fn configured_f7() -> String {
    assert!(
        cfg!(target_os = "linux"),
        "configured frozen-f7 gate requires Linux"
    );
    let path = std::env::var("IMPG_TEST_F7_FINITE")
        .expect("explicit frozen-f7 gate requires IMPG_TEST_F7_FINITE; no fallback path");
    let output = std::process::Command::new("sha256sum")
        .arg(&path)
        .output()
        .expect("configured Linux oracle gate requires sha256sum");
    assert!(
        output.status.success(),
        "configured frozen-f7 executable is missing/unreadable: {path}"
    );
    assert!(
        String::from_utf8_lossy(&output.stdout)
            .starts_with("a79200bd1ea1c4b3eacbe71effbe10adb6d61ec9904b381250c0d4052f8d8be4"),
        "configured frozen-f7 executable SHA256 mismatch"
    );
    path
}
