;;; Guix development environment for building impg from this checkout.
;;;
;;; Usage from this directory:
;;;   $HOME/.guix-profile/bin/guix shell \
;;;     -m manifest.scm
;;;
;;; Then build normally:
;;;   cargo build --release

(use-modules (guix profiles))

(specifications->manifest
 '("bash-minimal"
   "coreutils"
   "findutils"
   "git"
   "make"
   "nss-certs"
   "sed"
   "gcc-toolchain"
   "pkg-config"
   "cmake"
   "clang"
   "htslib"
   "glibc:static"
   "libdeflate"
   "zlib"
   "zstd"
   "bzip2"
   "xz"
   "openssl"))
