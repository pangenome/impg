;;; Guix development environment for building impg from this checkout.
;;;
;;; Usage from this directory:
;;;   /usr/local/guix-profiles/guix-pull/bin/guix shell \
;;;     -L /tmp/guix-bioinformatics-impg -m manifest.scm
;;;
;;; Then build normally:
;;;   cargo build --release

(use-modules (guix profiles)
             (gnu packages bash)
             (gnu packages base)
             (gnu packages bioinformatics)
             (gnu packages commencement)
             (gnu packages cmake)
             (gnu packages compression)
             (gnu packages jemalloc)
             (gnu packages llvm)
             (gnu packages maths)
             (gnu packages nss)
             (gnu packages pkg-config)
             (gnu packages rust)
             (gnu packages tls)
             (gnu packages version-control)
             (gn packages pangenome)
             (gn packages pangenome-rust))

(packages->manifest
 (list bash-minimal
       coreutils
       findutils
       git
       gnu-make
       nss-certs
       sed
       rust
       gcc-toolchain
       pkg-config
       cmake-minimal
       clang
       htslib
       libdeflate
       zlib
       zstd
       bzip2
       xz
       openssl
       jemalloc
       gsl
       spoa
       wfa2-lib/our
       wfmash
       fastga-rs))
