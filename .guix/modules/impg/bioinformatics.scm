;;; IMplicit Pangenome Graph (IMPG)
;;; Copyright © 2025 Frederick M. Muriithi <fredmanglis@gmail.com>

(define-module (impg bioinformatics)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system meson)
  #:use-module ((guix licenses) #:prefix license:)

  #:use-module (gnu packages cmake)
  #:use-module (gnu packages check)
  #:use-module (gnu packages assembly)
  #:use-module (gnu packages pkg-config)
  #:use-module (gnu packages compression)
  #:use-module (gnu packages bioinformatics))

(define biosoup-0.11.0
  (package
    (inherit biosoup)
    (name "biosoup")
    (version "0.11.0")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/rvaser/biosoup")
             ;; Corresponds to version 0.11.0
             (commit "3e31aa1d9039a2689241aebd18c45933b2d0f5e3")))
       (file-name (git-file-name name version))
       (sha256
        (base32
         "0vn1hj3h152iwahnrzghqll34qaphchi07klb3j70vgc248micbz"))))))

(define bioparser-3.1.0
  (package
    (inherit bioparser)
    (name "bioparser")
    (version "3.1.0")
    (source
     (origin
       (method git-fetch)
       (uri (git-reference
             (url "https://github.com/rvaser/bioparser")
             ;; Corresponds to tag 3.1.0
             (commit "4fa7126293d2a0eb90125b58fb704f0eed33ffe0")))
       (file-name (git-file-name name version))
       (sha256
        (base32
         "057zb3g8qyvbkbfzgkygrf0iphz3v4icm13pchxzrafiy7zkbmrq"))))
    (inputs
     (modify-inputs (package-inputs bioparser)
       (replace "biosoup" biosoup-0.11.0)))))

(define-public spoa
  (let ((commit "08957f6b87ce4262358a88c6b2c3c7860cf60239"))
    (package
      (name "spoa")
      (version "4.0.8")
      (source (origin
                (method git-fetch)
                (uri (git-reference (url "https://github.com/rvaser/spoa")
                                    (commit commit)))
                (file-name (git-file-name "spoa" (string-append version "-" (string-take commit 7))))
                (sha256 (base32 "0vafy9ry3cdrymxshcfmiv4schb0va3yxb6g3p20l54wl8alhxfj"))))
      (build-system meson-build-system)
      (inputs (list zlib
                    simde
                    biosoup-0.11.0
                    bioparser-3.1.0
                    pkg-config
                    googletest
                    cpu-features
                    cmake-minimal))
      (synopsis "C++ implementation of the Partial Order Alignment algorithm")
      (description "Spoa (SIMD POA) is a c++ implementation of the partial order
alignment (POA) algorithm (as described in 10.1093/bioinformatics/18.3.452)
which is used to generate consensus sequences (as described in
10.1093/bioinformatics/btg109).  It supports three alignment modes: local
(Smith-Waterman), global (Needleman-Wunsch) and semi-global alignment (overlap),
and three gap modes: linear, affine and convex (piecewise affine). It also
supports Intel SSE4.1+ and AVX2 vectorization (marginally faster due to high
latency shifts), SIMDe and dispatching.")
      (license license:expat)
      (home-page "https://github.com/rvaser/spoa"))))
