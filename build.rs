use std::fs;
use std::path::{Path, PathBuf};

fn main() {
    let out_dir = std::env::var("OUT_DIR").expect("OUT_DIR not set");

    // Derive the profile directory (e.g. target/release/) from OUT_DIR.
    // OUT_DIR is typically: target/release/build/<crate>-<hash>/out
    let profile_dir = PathBuf::from(&out_dir)
        .ancestors()
        .nth(3)
        .expect("cannot derive profile dir from OUT_DIR")
        .to_path_buf();

    let build_dir = profile_dir.join("build");

    // Copy wfmash binary
    copy_dep_binary(&build_dir, "wfmash-rs-", "wfmash", &profile_dir);

    // Copy wfmash's WFA2 shared libraries (needed at runtime)
    copy_dep_shared_libs(&build_dir, "wfmash-rs-", &profile_dir);

    // Copy FastGA binaries
    let fastga_binaries = [
        "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN", "ONEview",
    ];
    for name in &fastga_binaries {
        copy_dep_binary(&build_dir, "fastga-rs-", name, &profile_dir);
    }
}

/// Scan `build_dir` for a subdirectory matching `prefix*/out/<binary_name>`,
/// and copy it to `dest_dir`. Also copies to `$CARGO_HOME/bin/` if it exists
/// (handles `cargo install` case).
fn copy_dep_binary(build_dir: &Path, prefix: &str, binary_name: &str, dest_dir: &Path) {
    let Some(src) = find_binary_in_build(build_dir, prefix, binary_name) else {
        return; // silently skip if not found
    };

    // Copy to profile dir (e.g. target/release/)
    let dest = dest_dir.join(binary_name);
    if !dest.exists() {
        if let Err(e) = fs::copy(&src, &dest) {
            eprintln!("cargo:warning=Failed to copy {binary_name} to {}: {e}", dest.display());
        } else {
            set_executable(&dest);
        }
    }

    // Copy to CARGO_HOME/bin/ if it exists (for `cargo install`)
    if let Ok(cargo_home) = std::env::var("CARGO_HOME") {
        let cargo_bin = PathBuf::from(cargo_home).join("bin");
        if cargo_bin.is_dir() {
            let cargo_dest = cargo_bin.join(binary_name);
            if let Err(e) = fs::copy(&src, &cargo_dest) {
                eprintln!(
                    "cargo:warning=Failed to copy {binary_name} to {}: {e}",
                    cargo_dest.display()
                );
            } else {
                set_executable(&cargo_dest);
            }
        }
    }
}

/// Search `build_dir` for `<prefix><hash>/out/<binary_name>`.
fn find_binary_in_build(build_dir: &Path, prefix: &str, binary_name: &str) -> Option<PathBuf> {
    let entries = fs::read_dir(build_dir).ok()?;
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if name_str.starts_with(prefix) {
            let candidate = entry.path().join("out").join(binary_name);
            if candidate.is_file() {
                return Some(candidate);
            }
        }
    }
    None
}

/// Copy shared libraries (*.so*, *.dylib) from `<prefix><hash>/out/` subdirectories
/// to `dest_dir` and `$CARGO_HOME/lib/`. This handles wfmash's WFA2 runtime deps.
fn copy_dep_shared_libs(build_dir: &Path, prefix: &str, dest_dir: &Path) {
    let entries = match fs::read_dir(build_dir) {
        Ok(e) => e,
        Err(_) => return,
    };

    for entry in entries.flatten() {
        let name = entry.file_name();
        let name_str = name.to_string_lossy();
        if !name_str.starts_with(prefix) {
            continue;
        }

        // Search recursively under out/ for shared libraries
        let out_dir = entry.path().join("out");
        if !out_dir.is_dir() {
            continue;
        }

        collect_shared_libs(&out_dir, dest_dir);
    }
}

/// Recursively find and copy shared libraries from `dir` to `dest_dir`.
fn collect_shared_libs(dir: &Path, dest_dir: &Path) {
    let entries = match fs::read_dir(dir) {
        Ok(e) => e,
        Err(_) => return,
    };

    for entry in entries.flatten() {
        let path = entry.path();
        if path.is_dir() {
            collect_shared_libs(&path, dest_dir);
            continue;
        }

        let fname = entry.file_name();
        let fname_str = fname.to_string_lossy();

        // Match *.so* or *.dylib
        let is_shared_lib = fname_str.contains(".so") || fname_str.ends_with(".dylib");
        if !is_shared_lib || !path.is_file() {
            continue;
        }

        // Copy to profile dir
        let dest = dest_dir.join(&fname);
        if !dest.exists() {
            if let Err(e) = fs::copy(&path, &dest) {
                eprintln!(
                    "cargo:warning=Failed to copy {} to {}: {e}",
                    fname_str,
                    dest.display()
                );
            }
        }

        // Copy to CARGO_HOME/lib/ if it exists
        if let Ok(cargo_home) = std::env::var("CARGO_HOME") {
            let cargo_lib = PathBuf::from(cargo_home).join("lib");
            let _ = fs::create_dir_all(&cargo_lib);
            let cargo_dest = cargo_lib.join(&fname);
            if let Err(e) = fs::copy(&path, &cargo_dest) {
                eprintln!(
                    "cargo:warning=Failed to copy {} to {}: {e}",
                    fname_str,
                    cargo_dest.display()
                );
            }
        }
    }
}

#[cfg(unix)]
fn set_executable(path: &Path) {
    use std::os::unix::fs::PermissionsExt;
    if let Ok(metadata) = fs::metadata(path) {
        let mut perms = metadata.permissions();
        perms.set_mode(perms.mode() | 0o111);
        let _ = fs::set_permissions(path, perms);
    }
}

#[cfg(not(unix))]
fn set_executable(_path: &Path) {}
