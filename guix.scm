;; To use this file to build a version of impg using git HEAD:
;;
;;   guix build -L . -f guix.scm                  # default build
;;
;; To get a development container using a recent guix (see `guix pull`)
;;
;;   guix shell --share=$HOME/.cargo -C -L . impg-shell-git
;;
;; and inside the container
;;
;;   rm -rf target/  # may be necessary when check-for-pregenerated-files fails
;;   CC=gcc cargo build --release
;;
;; note that we don't expect cargo to download packages, so ignore resolve errors
;;
;; List other packages in this guix.scm file
;;
;;   guix package -L . -A impg
;;
;; Installing guix (note that Debian comes with guix). Once installed update as a normal user with:
;;
;;   mkdir ~/opt
;;   guix pull -p ~/opt/guix # update guix takes a while - don't do this often!
;;
;; Use the update guix to build impg:
;;
;;   ~/opt/guix/bin/guix build -L . -f guix.scm
;;
;; Or get a shell (see above)
;;
;; If things do not work you may also have to update the guix-daemon in systemd. Guix mostly downloads binary
;; substitutes. If it wants to build a lot of software you probably have substitutes misconfigured.

;; by Pjotr Prins (c) 2025

(define-module (guix)
  #:use-module ((guix licenses) #:prefix license:)
  ;; #:use-module (guix build-system cmake)
  #:use-module (guix build-system cargo)
  #:use-module (guix download)
  #:use-module (guix gexp)
  #:use-module (guix git-download)
  #:use-module (guix packages)
  #:use-module (guix utils)
  #:use-module (gnu packages algebra)
  #:use-module (gnu packages base)
  #:use-module (gnu packages bash)
  #:use-module (gnu packages bioinformatics)
  #:use-module (gnu packages build-tools)
  #:use-module (gnu packages certs)
  #:use-module (gnu packages cmake)
  #:use-module (gnu packages commencement)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages cpp)
  #:use-module (gnu packages crates-io)
  #:use-module (gnu packages curl)
  #:use-module (gnu packages gcc)
  #:use-module (gnu packages jemalloc)
  #:use-module (gnu packages linux) ; for util-linux column
  #:use-module (gnu packages llvm)
  #:use-module (gnu packages maths)
  #:use-module (gnu packages multiprecision)
  #:use-module (gnu packages perl)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages python)
  #:use-module (gnu packages rust)
  #:use-module (gnu packages rust-apps) ; for cargo
  #:use-module (gnu packages tls)
  #:use-module (gnu packages version-control)
  #:use-module (srfi srfi-1)
  #:use-module (ice-9 popen)
  #:use-module (ice-9 rdelim)
  #:use-module (deps))

(define %source-dir (dirname (current-filename)))

(define %version
  (read-string (open-pipe "git describe --always --tags --long|tr -d $'\n'" OPEN_READ)))

(define-public rust-coitrees-0.4
  (package
    (name "rust-coitrees")
    (version "0.4.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "coitrees" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1qwb4c5gx30gl1kyi85rbq6z23l2f9lm0q02ym160n0fvc89c3r4"))))
    (build-system cargo-build-system)
    (arguments
     `(#:cargo-development-inputs (("rust-clap" ,rust-clap-4)
                                   ("rust-fnv" ,rust-fnv-1)
                                   ("rust-libc" ,rust-libc-0.2)
                                   ("rust-rand" ,rust-rand-0.8))))
    (home-page "https://github.com/dcjones/coitrees")
    (synopsis
     "very fast data structure for overlap queries on sets of intervals.")
    (description
     "This package provides a very fast data structure for overlap queries on sets of
intervals.")
    (license license:expat)))


(define-public impg-base-git
  (package
    (name "impg-base-git")
    (version %version)
    (source (local-file %source-dir #:recursive? #t))
    (build-system cargo-build-system)
    (inputs (list curl gnutls lzip openssl pkg-config zlib xz)) ;; mostly for htslib
    (arguments
     `(#:cargo-inputs (("rust-bincode" ,rust-bincode-1)
                       ("rust-coitrees" ,rust-coitrees-0.4)
                       ("rust-env-logger" ,rust-env-logger-0.11)
                       ("rust-natord" ,rust-natord-1)
                       ("rust-noodles" ,rust-noodles-0.99)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-spoa-rs" ,rust-spoa-rs-0.1)
                       )
       ;; #:cargo-development-inputs ()))
       #:cargo-package-flags '("--no-metadata" "--no-verify" "--allow-dirty")
     ))
    (synopsis "impg")
    (description
     "Tool for implicit graphs - i.e. building smaller 'local' pangenomes from positioned sequences.")
    (home-page "https://github.com/pangenome/impg")
    (license license:expat)))

(define-public impg-shell-git
  "Shell version to use 'cargo build'"
  (package
    (inherit impg-base-git)
    (name "impg-shell-git")
    (inputs
     (modify-inputs (package-inputs impg-base-git)
         (append binutils coreutils-minimal ;; for the shell
                 )))
    (propagated-inputs (list cmake rust rust-cargo nss-certs openssl perl gnu-make-4.2
                             coreutils-minimal which perl binutils gcc-toolchain pkg-config zlib
                             )) ;; to run cargo build in the shell
    (arguments
     `(
       #:cargo-development-inputs ()
       #:cargo-package-flags '("--no-metadata" "--no-verify" "--allow-dirty")
       #:phases (modify-phases %standard-phases
                               (delete 'check-for-pregenerated-files)
                               (delete 'configure)
                               (delete 'build)
                               (delete 'package)
                               (delete 'check)
                               (delete 'install)
                               )))))

impg-base-git ;; default deployment build with debug info
