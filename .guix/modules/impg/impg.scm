;;; IMplicit Pangenome Graph (IMPG)
;;; Copyright © 2025 Frederick M. Muriithi <fredmanglis@gmail.com>

(define-module (impg impg)
  #:use-module (guix gexp)
  #:use-module (guix utils)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system cargo)
  #:use-module ((guix licenses) #:prefix license:)

  #:use-module (gnu packages cmake)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages compression)

  #:use-module (impg bioinformatics))

(define %source-dir (dirname (dirname (dirname (current-source-directory)))))

(define vcs-file?
  (or (git-predicate %source-dir)
      (const #t)))

(define-public impg
  (package
    (name "impg")
    (version "0.3.1")
    (source
     (local-file "../../.."
                 "impg-checkout"
                 #:recursive? #t
                 #:select? vcs-file?))
    (build-system cargo-build-system)
    (arguments
     (list
      #:phases
      #~(modify-phases %standard-phases
          (add-before 'build 'remove-agc-as-core-dependency
            (lambda _
              (substitute* "Cargo.toml"
                (("^.*agc-rs.*$") "")
                (("^default =.*") ""))))
          (add-before 'build 'fix-dependency-sources
            (lambda _
              (substitute* "Cargo.toml"
                (("git = \"https://github.com/chfi/rs-handlegraph\", rev = \"3ac575e4216ce16a16667503a8875e469a40a97a\"")
                 "path = \"guix-vendor/rust-handlegraph-0.7.0-alpha.9.3ac575e-checkout\", version = \"0.7.0-alpha.9\"")
                (("git = \"https://github.com/AndreaGuarracino/spoa-rs.git\"")
                 "path = \"guix-vendor/rust-spoa-rs-0.1.0.6f4f102-checkout\", version = \"0.1.0\""))))
          (add-before 'build 'patch-include-paths-to-spoa
            (lambda _
              (substitute* "guix-vendor/rust-spoa-rs-0.1.0.6f4f102-checkout/build.rs"
                (("spoa/include") #$(file-append spoa "/include"))
                (("^ *out_dir\\.display.*$") "")
                (("\\{\\}/build/lib\",") #$(file-append spoa "/lib\""))))))))
    (inputs (cons* spoa
                   pkg-config
                   cmake-minimal
                   (list zstd "lib")
                   (cargo-inputs 'impg #:module '(impg rust-crates))))
    (synopsis "Tool to extract and compute graphs of sequences and alignments")
    (description "IMplicit Pangenome Graph (impg) is a tool that takes in
sequences and their relative alignments and extracts sections of the sequences
and alignments for analysis. It can also be used to compute the graphs for the
sequences.")
    (license license:expat)
    (home-page "https://github.com/pangenome/impg")))

impg
