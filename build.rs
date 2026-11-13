use std::fs;
use std::path::{Path, PathBuf};

fn main() {
    // Compile syng C library
    compile_syng();

    let out_dir = std::env::var("OUT_DIR").expect("OUT_DIR not set");

    // OUT_DIR is <profile>/build/<crate>-<hash>/out (cargo < 1.100)
    // or <profile>/build/<crate>/<hash>/out (cargo >= 1.100).
    let build_dir = PathBuf::from(&out_dir)
        .ancestors()
        .find(|p| p.file_name().is_some_and(|n| n == "build"))
        .unwrap_or_else(|| panic!("no build directory above OUT_DIR {out_dir}"))
        .to_path_buf();
    let profile_dir = build_dir
        .parent()
        .unwrap_or_else(|| panic!("no profile directory above {}", build_dir.display()))
        .to_path_buf();

    // Copy wfmash binary (required at runtime)
    copy_dep_binary(&build_dir, "wfmash-rs", "wfmash", &profile_dir, true);

    // Copy wfmash's WFA2 shared libraries (needed at runtime)
    copy_dep_shared_libs(&build_dir, "wfmash-rs", &profile_dir);

    // Copy FastGA binaries
    let fastga_binaries = [
        "FastGA", "FAtoGDB", "GIXmake", "GIXrm", "ALNtoPAF", "PAFtoALN", "ONEview",
    ];
    for name in &fastga_binaries {
        copy_dep_binary(&build_dir, "fastga-rs", name, &profile_dir, false);
    }
}

/// List the `out/` directories of the build-script runs of `krate`, in both
/// cargo layouts: `<krate>-<hash>/out` and `<krate>/<hash>/out`.
fn dep_out_dirs(build_dir: &Path, krate: &str) -> Vec<PathBuf> {
    let mut dirs = Vec::new();
    let entries = fs::read_dir(build_dir)
        .unwrap_or_else(|e| panic!("cannot read build directory {}: {e}", build_dir.display()));
    for entry in entries.flatten() {
        let name = entry.file_name();
        let name = name.to_string_lossy();
        if name.strip_prefix(krate).is_some_and(|rest| rest.starts_with('-')) {
            dirs.push(entry.path().join("out"));
        } else if name == krate {
            let runs = fs::read_dir(entry.path())
                .unwrap_or_else(|e| panic!("cannot read {}: {e}", entry.path().display()));
            dirs.extend(runs.flatten().map(|run| run.path().join("out")));
        }
    }
    dirs.retain(|d| d.is_dir());
    dirs
}

/// Find `<binary_name>` in the build output of `krate` and copy it to
/// `dest_dir`. Also copies to `$CARGO_HOME/bin/` if it exists (handles
/// `cargo install`). A missing binary fails the build if `required`, else it
/// is reported as a cargo warning.
fn copy_dep_binary(
    build_dir: &Path,
    krate: &str,
    binary_name: &str,
    dest_dir: &Path,
    required: bool,
) {
    let Some(src) = find_binary_in_build(build_dir, krate, binary_name) else {
        let msg = format!(
            "{binary_name} not found in the {krate} build output under {}",
            build_dir.display()
        );
        if required {
            panic!("{msg}");
        }
        println!("cargo:warning={msg}");
        return;
    };

    // Copy to profile dir (e.g. target/release/). Always overwrite: a stale
    // copy from an older dependency revision must not survive a version bump.
    let dest = dest_dir.join(binary_name);
    if let Err(e) = fs::copy(&src, &dest) {
        println!(
            "cargo:warning=Failed to copy {binary_name} to {}: {e}",
            dest.display()
        );
    } else {
        set_executable(&dest);
    }

    // Copy to CARGO_HOME/bin/ if it exists (for `cargo install`)
    if let Ok(cargo_home) = std::env::var("CARGO_HOME") {
        let cargo_bin = PathBuf::from(cargo_home).join("bin");
        if cargo_bin.is_dir() {
            let cargo_dest = cargo_bin.join(binary_name);
            if let Err(e) = fs::copy(&src, &cargo_dest) {
                println!(
                    "cargo:warning=Failed to copy {binary_name} to {}: {e}",
                    cargo_dest.display()
                );
            } else {
                set_executable(&cargo_dest);
            }
        }
    }
}

/// Search the build output of `krate` for `<binary_name>`.
///
/// Several hash directories can coexist, one per dependency revision. Take the
/// most recently modified one, otherwise a bumped revision keeps losing to a
/// leftover build of the old one.
fn find_binary_in_build(build_dir: &Path, krate: &str, binary_name: &str) -> Option<PathBuf> {
    dep_out_dirs(build_dir, krate)
        .into_iter()
        .map(|out| out.join(binary_name))
        .filter(|candidate| candidate.is_file())
        .max_by_key(|candidate| {
            fs::metadata(candidate)
                .and_then(|m| m.modified())
                .unwrap_or(std::time::UNIX_EPOCH)
        })
}

/// Copy shared libraries (*.so*, *.dylib) from the build output of `krate`
/// to `dest_dir` and `$CARGO_HOME/lib/`. This handles wfmash's WFA2 runtime deps.
fn copy_dep_shared_libs(build_dir: &Path, krate: &str, dest_dir: &Path) {
    for out_dir in dep_out_dirs(build_dir, krate) {
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

        // Copy to profile dir, overwriting any stale copy
        let dest = dest_dir.join(&fname);
        if let Err(e) = fs::copy(&path, &dest) {
            println!(
                "cargo:warning=Failed to copy {} to {}: {e}",
                fname_str,
                dest.display()
            );
        }

        // Copy to CARGO_HOME/lib/ if it exists
        if let Ok(cargo_home) = std::env::var("CARGO_HOME") {
            let cargo_lib = PathBuf::from(cargo_home).join("lib");
            let _ = fs::create_dir_all(&cargo_lib);
            let cargo_dest = cargo_lib.join(&fname);
            if let Err(e) = fs::copy(&path, &cargo_dest) {
                println!(
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

/// Compile syng C files into libsyng.a and link it.
fn compile_syng() {
    let syng_dir = PathBuf::from("vendor/syng");

    let c_files = [
        "syngbwt3.c",
        "rskip.c",
        "kmerhash.c",
        "seqhash.c",
        "syncmerset.c",
        "seqio.c",
        "ONElib.c",
        "utils.c",
        "array.c",
        "hash.c",
        "dict.c",
        "impg_syng_helpers.c",
    ];

    let mut build = cc::Build::new();
    build
        .include(&syng_dir)
        .opt_level(3)
        .define("ONEIO", None)
        // syng was written where char is signed. It is unsigned on aarch64, and
        // there syngBWTadvanceRank walked to a negative node and aborted.
        .flag_if_supported("-fsigned-char")
        .warnings(false);

    for file in &c_files {
        build.file(syng_dir.join(file));
    }

    build.compile("syng");

    // Link zlib (required by seqio.c)
    println!("cargo:rustc-link-lib=z");
    // Rerun if any syng source changes
    println!("cargo:rerun-if-changed=vendor/syng");
}
