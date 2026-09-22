#!/usr/bin/env bash
#
# Build a source tarball that contains the git submodules.
#
# GitHub's auto-generated tag tarballs omit submodules, so they cannot build
# impg (vendor/gfaffix is a binary target, vendor/syng is compiled by build.rs).
# Upload the tarball this script writes as a release asset, and point the
# bioconda recipe at it.
#
# Usage: scripts/make-release-tarball.sh v0.5.0

set -euxo pipefail

TAG="${1:?usage: $0 <tag>   e.g. $0 v0.5.0}"

REPO_ROOT="$(git rev-parse --show-toplevel)"
cd "$REPO_ROOT"

PREFIX="impg-${TAG}"
OUT="${REPO_ROOT}/${PREFIX}.tar.gz"
STAGE="$(mktemp -d)"
trap 'rm -rf "$STAGE"' EXIT

echo "[tarball] staging ${PREFIX} from $(git rev-parse --short HEAD)"
git archive --prefix="${PREFIX}/" HEAD | tar -x -C "$STAGE"

echo "[tarball] adding submodules"
git submodule update --init --recursive
while read -r sm_path; do
  echo "[tarball]   ${sm_path}"
  git -C "$sm_path" archive --prefix="${PREFIX}/${sm_path}/" HEAD | tar -x -C "$STAGE"
done < <(git submodule status --recursive | awk '{print $2}')

echo "[tarball] verifying submodule contents landed"
test -f "${STAGE}/${PREFIX}/vendor/gfaffix/src/main.rs"
test -f "${STAGE}/${PREFIX}/vendor/syng/syngbwt3.c"
test -f "${STAGE}/${PREFIX}/Cargo.lock"

echo "[tarball] writing ${OUT}"
rm -f "$OUT"
tar -czf "$OUT" -C "$STAGE" "$PREFIX"

set +x
echo
echo "wrote  $OUT"
echo "size   $(du -h "$OUT" | cut -f1)"
echo "sha256 $(sha256sum "$OUT" | cut -d' ' -f1)"
