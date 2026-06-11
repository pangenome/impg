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
# The guix-bioinformatics channel checkout is expected at
# /tmp/guix-bioinformatics-impg by default.  Override it with:
#   export GUIX_BIOINFORMATICS=/path/to/guix-bioinformatics

set -euo pipefail

_impg_dir="$(cd "$(dirname "${BASH_SOURCE[0]:-$0}")" && pwd)"
_guix="${GUIX:-/usr/local/guix-profiles/guix-pull/bin/guix}"
_guix_bio="${GUIX_BIOINFORMATICS:-/tmp/guix-bioinformatics-impg}"
_manifest="$_impg_dir/manifest.scm"

if [[ ! -x "$_guix" ]]; then
    echo "error: Guix command not found or not executable: $_guix" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ ! -d "$_guix_bio/gn/packages" ]]; then
    echo "error: guix-bioinformatics checkout not found: $_guix_bio" >&2
    echo "set GUIX_BIOINFORMATICS=/path/to/guix-bioinformatics" >&2
    return 1 2>/dev/null || exit 1
fi

if [[ "${BASH_SOURCE[0]}" == "$0" ]]; then
    exec "$_guix" shell --pure -L "$_guix_bio" -m "$_manifest" -- \
        bash --noprofile --norc -c '
            export CC=gcc
            export CXX=g++
            export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc
            export LIBCLANG_PATH="${GUIX_ENVIRONMENT}/lib"
            export SSL_CERT_FILE="${GUIX_ENVIRONMENT}/etc/ssl/certs/ca-certificates.crt"
            echo "impg Guix build shell ready: cargo $(cargo --version | cut -d" " -f2), rustc $(rustc --version | cut -d" " -f2)"
            exec bash --noprofile --norc -i
        '
fi

set +u
eval "$("$_guix" shell -L "$_guix_bio" -m "$_manifest" --search-paths)"
set -u
export CC=gcc
export CXX=g++
export CARGO_TARGET_X86_64_UNKNOWN_LINUX_GNU_LINKER=gcc

_profile_bin="${PATH%%:*}"
_profile="${_profile_bin%/bin}"
export LIBCLANG_PATH="$_profile/lib"
export SSL_CERT_FILE="$_profile/etc/ssl/certs/ca-certificates.crt"

echo "impg Guix build env ready: cargo $(cargo --version | cut -d' ' -f2), rustc $(rustc --version | cut -d' ' -f2)"

unset _impg_dir _guix _guix_bio _manifest _profile_bin _profile
