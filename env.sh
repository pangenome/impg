#!/usr/bin/env bash
# Enter or source the Guix environment used to build impg from this checkout.
#
# Start an isolated shell:
#   ./env.sh
#   cargo build --release --locked
#
# Or add the Guix profile to the current shell:
#   source ./env.sh
#   cargo build --release --locked
#
set -euo pipefail

_impg_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_guix="${GUIX:-${HOME}/.guix-profile/bin/guix}"
_manifest="$_impg_dir/manifest.scm"

if [[ ! -x "$_guix" ]]; then
    echo "error: Guix command not found or not executable: $_guix" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    exec bash --noprofile --norc -c 'source "$1"; exec bash --noprofile --norc -i' bash "$_impg_dir/env.sh"
fi

set +u
_search_paths="$("$_guix" shell -m "$_manifest" --search-paths)"
if [[ -z "$_search_paths" ]]; then
    echo "error: Guix did not return a usable environment" >&2
    return 1 2>/dev/null || exit 1
fi
eval "$_search_paths"
set -u
export CC=gcc
export CXX=g++
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc

_profile_bin="${PATH%%:*}"
_profile="${_profile_bin%/bin}"
export LIBCLANG_PATH="$_profile/lib"
export SSL_CERT_FILE="$_profile/etc/ssl/certs/ca-certificates.crt"

_ld_so="$(gcc -print-file-name=ld-linux-x86-64.so.2)"
_gcc_lib="$(dirname "$(gcc -print-file-name=libgcc_s.so.1)")"
_shim_dir="${TMPDIR:-/tmp}/impg-guix-rust-${USER:-$(id -u)}-${_profile##*/}"

if command -v rustup >/dev/null 2>&1; then
    _cargo_real="$(rustup which cargo)"
    _rustc_real="$(rustup which rustc)"
    _rustdoc_real="$(rustup which rustdoc)"
else
    _cargo_real="$(command -v cargo)"
    _rustc_real="$(command -v rustc)"
    _rustdoc_real="$(command -v rustdoc)"
fi
_rust_lib="$(dirname "$_rustc_real")/../lib"

mkdir -p "$_shim_dir"
_write_rust_wrapper() {
    local _name="$1"
    local _real="$2"
    {
        printf '%s\n' '#!/usr/bin/env bash'
        printf 'exec %q --library-path %q %q "$@"\n' \
            "$_ld_so" "$_profile/lib:$_gcc_lib:$_rust_lib" "$_real"
    } >"$_shim_dir/$_name"
    chmod +x "$_shim_dir/$_name"
}
_write_rust_wrapper cargo "$_cargo_real"
_write_rust_wrapper rustc "$_rustc_real"
_write_rust_wrapper rustdoc "$_rustdoc_real"

export PATH="$_shim_dir:$PATH"
export CARGO="$_shim_dir/cargo"
export RUSTC="$_shim_dir/rustc"
export RUSTDOC="$_shim_dir/rustdoc"

echo "impg Guix build env ready: cargo $(cargo --version | cut -d' ' -f2), rustc $(rustc --version | cut -d' ' -f2)"

unset _impg_dir _guix _manifest _search_paths _profile_bin _profile
unset _ld_so _gcc_lib _shim_dir _cargo_real _rustc_real _rustdoc_real _rust_lib
unset -f _write_rust_wrapper
