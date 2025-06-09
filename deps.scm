(define-module (deps)
  #:use-module ((guix licenses) #:prefix license:)
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
  #:use-module (gnu packages crates-check)
  #:use-module (gnu packages crates-compression)
  #:use-module (gnu packages crates-crypto)
  #:use-module (gnu packages crates-io)
  #:use-module (gnu packages crates-web)
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
  )

(define-public rust-spoa-rs-0.1
  (package
    (name "rust-spoa-rs")
    (version "0.1.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "spoa_rs" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "014q4d9y4imm6s1njxnw06wmax5kq77irv1nsbc3n4swkpnl7a5h"))))
    (build-system cargo-build-system)
    (arguments
     `(#:cargo-inputs (("rust-cmake" ,rust-cmake-0.1)
                       ("rust-cxx" ,rust-cxx-1)
                       ("rust-cxx-build" ,rust-cxx-build-1))))
    (home-page "https://github.com/AndreaGuarracino/spoa-rs")
    (synopsis "Rust bindings for SPOA (SIMD POA)")
    (description "This package provides Rust bindings for SPOA (SIMD POA).")
    (license license:bsd-3)))

(define-public rust-noodles-refget-0.7
  (package
    (name "rust-noodles-refget")
    (version "0.7.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-refget" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0v4i51423gmsnxw1l99vk3rb96jj7y428s4k50ygy8rwpjn8l6fy"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bytes" ,rust-bytes-1)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-reqwest" ,rust-reqwest-0.12)
                       ("rust-serde" ,rust-serde-1)
                       ("rust-url" ,rust-url-2))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "refget client")
    (description "This package provides a refget client.")
    (license license:expat)))

(define-public rust-noodles-htsget-0.8
  (package
    (name "rust-noodles-htsget")
    (version "0.8.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-htsget" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "10c9abvslwfpvh7v8asp061j3230l0djrhvslzvmrnnxby0bcbsg"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-base64" ,rust-base64-0.22)
                       ("rust-bytes" ,rust-bytes-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-reqwest" ,rust-reqwest-0.12)
                       ("rust-serde" ,rust-serde-1)
                       ("rust-url" ,rust-url-2))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "An htsget client")
    (description "This package provides An htsget client.")
    (license license:expat)))

(define-public rust-noodles-gtf-0.45
  (package
    (name "rust-noodles-gtf")
    (version "0.45.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-gtf" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1g4hc9ws86yy187kdhss1j21p01j5n1hhx8sxrd5zraf6dy76g1g"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-lexical-core" ,rust-lexical-core-1)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-gff" ,rust-noodles-gff-0.50))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Gene Transfer Format (GTF) reader and writer")
    (description
     "This package provides Gene Transfer Format (GTF) reader and writer.")
    (license license:expat)))

(define-public rust-noodles-gff-0.50
  (package
    (name "rust-noodles-gff")
    (version "0.50.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-gff" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1a2jc118k321zxs4q71z3yw11nal8y1g3blmnqfllinhvndpd2f6"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-lexical-core" ,rust-lexical-core-1)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-percent-encoding" ,rust-percent-encoding-2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Generic Feature Format (GFF) reader and writer")
    (description
     "This package provides Generic Feature Format (GFF) reader and writer.")
    (license license:expat)))

(define-public rust-noodles-fastq-0.19
  (package
    (name "rust-noodles-fastq")
    (version "0.19.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-fastq" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "08s20w0dy9nrnj4nd93m5pwfkp1bpmi4czhkxf0m3jha5mivz6jh"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "FASTQ format reader and writer")
    (description "This package provides FASTQ format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-fasta-0.54
  (package
    (name "rust-noodles-fasta")
    (version "0.54.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-fasta" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1bmkw5cd7whdn3d6ygm78yh23x8bjbj9lmgsj4pgfvbwsjjiq62v"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-bytes" ,rust-bytes-1)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "FASTA format reader and writer")
    (description "This package provides FASTA format reader and writer.")
    (license license:expat)))

(define-public rust-libbz2-rs-sys-0.1
  (package
    (name "rust-libbz2-rs-sys")
    (version "0.1.3")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "libbz2-rs-sys" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1xdwpvw49v2llmdwzlxdcvcpnfxqa37c9hk9dchkd7h1il6a0r08"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-libc" ,rust-libc-0.2))))
    (home-page "https://github.com/trifectatechfoundation/libbzip2-rs")
    (synopsis "a drop-in compatible rust bzip2 implementation")
    (description
     "This package provides a drop-in compatible rust bzip2 implementation.")
    (license license:expat)))

(define-public rust-bzip2-sys-0.1
  (package
    (name "rust-bzip2-sys")
    (version "0.1.13+1.0.8")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "bzip2-sys" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "056c39pgjh4272bdslv445f5ry64xvb0f7nph3z7860ln8rzynr2"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-cc" ,rust-cc-1)
                       ("rust-pkg-config" ,rust-pkg-config-0.3))))
    (home-page "https://github.com/alexcrichton/bzip2-rs")
    (synopsis
     "Bindings to libbzip2 for bzip2 compression and decompression exposed as
Reader/Writer streams.")
    (description
     "This package provides Bindings to libbzip2 for bzip2 compression and decompression exposed as
Reader/Writer streams.")
    (license (list license:expat license:asl2.0))))

(define-public rust-bzip2-0.5
  (package
    (name "rust-bzip2")
    (version "0.5.2")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "bzip2" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0iya6nbj0p2y8jss0z05yncc5hadry164fw3zva01y06v4igpv29"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bzip2-sys" ,rust-bzip2-sys-0.1)
                       ("rust-libbz2-rs-sys" ,rust-libbz2-rs-sys-0.1))))
    (home-page "https://github.com/trifectatechfoundation/bzip2-rs")
    (synopsis
     "Bindings to libbzip2 for bzip2 compression and decompression exposed as
Reader/Writer streams.")
    (description
     "This package provides Bindings to libbzip2 for bzip2 compression and decompression exposed as
Reader/Writer streams.")
    (license (list license:expat license:asl2.0))))

(define-public rust-noodles-cram-0.84
  (package
    (name "rust-noodles-cram")
    (version "0.84.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-cram" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "17dr07nbyzj8h38352k93bpjsjlxshkdbf4a1wa9z9aw99lq069s"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-async-compression" ,rust-async-compression-0.4)
                       ("rust-bitflags" ,rust-bitflags-2)
                       ("rust-bstr" ,rust-bstr-1)
                       ("rust-byteorder" ,rust-byteorder-1)
                       ("rust-bzip2" ,rust-bzip2-0.5)
                       ("rust-flate2" ,rust-flate2-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-libdeflater" ,rust-libdeflater-1)
                       ("rust-md-5" ,rust-md-5-0.10)
                       ("rust-noodles-bam" ,rust-noodles-bam-0.81)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-fasta" ,rust-noodles-fasta-0.54)
                       ("rust-noodles-sam" ,rust-noodles-sam-0.77)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1)
                       ("rust-xz2" ,rust-xz2-0.1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "CRAM format reader and writer")
    (description "This package provides CRAM format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-bed-0.26
  (package
    (name "rust-noodles-bed")
    (version "0.26.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-bed" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1iqzy6a2qf7rn05rn3xyvhjlhc1qkcw94if8j660dhx4ga02h5iv"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-lexical-core" ,rust-lexical-core-1)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-tabix" ,rust-noodles-tabix-0.55))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "BED (Browser Extensible Data) reader")
    (description "This package provides BED (Browser Extensible Data) reader.")
    (license license:expat)))

(define-public rust-noodles-tabix-0.55
  (package
    (name "rust-noodles-tabix")
    (version "0.55.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-tabix" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0lf779jp9b7j9wn2pnx5srhpp31xr0qnq86fvdyqvhrfxrrn26s6"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-byteorder" ,rust-byteorder-1)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Tabix (TBI) format reader and writer")
    (description "This package provides Tabix (TBI) format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-vcf-0.79
  (package
    (name "rust-noodles-vcf")
    (version "0.79.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-vcf" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1jcxll1r51v2kxcbdjmabrrnlfyh78fb3dvic2wx7byyrhw578k6"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-tabix" ,rust-noodles-tabix-0.55)
                       ("rust-percent-encoding" ,rust-percent-encoding-2)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Variant Call Format (VCF) reader and writer")
    (description
     "This package provides Variant Call Format (VCF) reader and writer.")
    (license license:expat)))

(define-public rust-noodles-bcf-0.76
  (package
    (name "rust-noodles-bcf")
    (version "0.76.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-bcf" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1zljwpjgqha81hq23z09ps09yzdllv8kqvrqyfwjfhsld9whk8fk"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-byteorder" ,rust-byteorder-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-vcf" ,rust-noodles-vcf-0.79)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Binary Call Format (BCF) reader and writer")
    (description
     "This package provides Binary Call Format (BCF) reader and writer.")
    (license license:expat)))

(define-public rust-lexical-write-integer-1
  (package
    (name "rust-lexical-write-integer")
    (version "1.0.5")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-write-integer" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0y7rl3pkac4lhcfiwxzsb2p3m432if4af5jn4kxkda0lm7qxz7b2"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-lexical-util" ,rust-lexical-util-1)
                       ("rust-static-assertions" ,rust-static-assertions-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Efficient formatting of integers to strings")
    (description
     "This package provides Efficient formatting of integers to strings.")
    (license (list license:expat license:asl2.0))))

(define-public rust-lexical-write-float-1
  (package
    (name "rust-lexical-write-float")
    (version "1.0.5")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-write-float" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1gb2ip3r9wmbsgfa5d8cwgbw4hrgpyv5g9w1bas0yikzl9lcdby5"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-lexical-util" ,rust-lexical-util-1)
                       ("rust-lexical-write-integer" ,rust-lexical-write-integer-1)
                       ("rust-static-assertions" ,rust-static-assertions-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Efficient formatting of floats to strings")
    (description
     "This package provides Efficient formatting of floats to strings.")
    (license (list license:expat license:asl2.0))))

(define-public rust-lexical-util-1
  (package
    (name "rust-lexical-util")
    (version "1.0.6")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-util" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1cx0974y9x63ikra6l2vqcr4gmf8pipdrfzzfz0j9z9pym5y50js"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-static-assertions" ,rust-static-assertions-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Shared utilities for lexical creates")
    (description "This package provides Shared utilities for lexical creates.")
    (license (list license:expat license:asl2.0))))

(define-public rust-lexical-parse-integer-1
  (package
    (name "rust-lexical-parse-integer")
    (version "1.0.5")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-parse-integer" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0bpidpb0viqj9wx3z727y6d59smz5kj7km5nlwdi42pw4ap7l83j"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-lexical-util" ,rust-lexical-util-1)
                       ("rust-static-assertions" ,rust-static-assertions-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Efficient parsing of integers from strings")
    (description
     "This package provides Efficient parsing of integers from strings.")
    (license (list license:expat license:asl2.0))))

(define-public rust-lexical-parse-float-1
  (package
    (name "rust-lexical-parse-float")
    (version "1.0.5")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-parse-float" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1wmhcndf7gvfqvmd5v61nhbqgaybiw27q1cs41h81c5h3yq9qvyy"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-lexical-parse-integer" ,rust-lexical-parse-integer-1)
                       ("rust-lexical-util" ,rust-lexical-util-1)
                       ("rust-static-assertions" ,rust-static-assertions-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Efficient parsing of floats from strings")
    (description
     "This package provides Efficient parsing of floats from strings.")
    (license (list license:expat license:asl2.0))))

(define-public rust-lexical-core-1
  (package
    (name "rust-lexical-core")
    (version "1.0.5")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "lexical-core" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0n49vqm7njn1ia0z9jkyvap864i808abgd3hb9b7b43014cc6rdp"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-lexical-parse-float" ,rust-lexical-parse-float-1)
                       ("rust-lexical-parse-integer" ,rust-lexical-parse-integer-1)
                       ("rust-lexical-util" ,rust-lexical-util-1)
                       ("rust-lexical-write-float" ,rust-lexical-write-float-1)
                       ("rust-lexical-write-integer" ,rust-lexical-write-integer-1))))
    (home-page "https://github.com/Alexhuszagh/rust-lexical")
    (synopsis "Lexical, to- and from-string conversion routines")
    (description
     "This package provides Lexical, to- and from-string conversion routines.")
    (license (list license:expat license:asl2.0))))

(define-public rust-noodles-sam-0.77
  (package
    (name "rust-noodles-sam")
    (version "0.77.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-sam" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0a5w8brjpfx3wr306kvw1ldfw54nkxmyy44fvynbb4x8v4xr2a2x"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bitflags" ,rust-bitflags-2)
                       ("rust-bstr" ,rust-bstr-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-lexical-core" ,rust-lexical-core-1)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Sequence Alignment/Map (SAM) format reader and writer")
    (description
     "This package provides Sequence Alignment/Map (SAM) format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-csi-0.49
  (package
    (name "rust-noodles-csi")
    (version "0.49.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-csi" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0hlaw76c3b3r1w5zxn4yjqjdyzx53mzv8mwv3a980gap71fmi4f2"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bit-vec" ,rust-bit-vec-0.8)
                       ("rust-bstr" ,rust-bstr-1)
                       ("rust-byteorder" ,rust-byteorder-1)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Coordinate-sorted index (CSI) format reader and writer")
    (description
     "This package provides Coordinate-sorted index (CSI) format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-core-0.17
  (package
    (name "rust-noodles-core")
    (version "0.17.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-core" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0jqb1d1llg8nc9khhwgbjjsknb8f37cyjqah37bx2zb5b0hpb2ff"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Shared utilities when working with noodles")
    (description
     "This package provides Shared utilities when working with noodles.")
    (license license:expat)))

(define-public rust-miniz-oxide-0.8
  (package
    (name "rust-miniz-oxide")
    (version "0.8.8")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "miniz_oxide" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0al9iy33flfgxawj789w2c8xxwg1n2r5vv6m6p5hl2fvd2vlgriv"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-adler2" ,rust-adler2-2)
                       ("rust-compiler-builtins" ,rust-compiler-builtins-0.1)
                       ("rust-rustc-std-workspace-alloc" ,rust-rustc-std-workspace-alloc-1)
                       ("rust-rustc-std-workspace-core" ,rust-rustc-std-workspace-core-1)
                       ("rust-serde" ,rust-serde-1)
                       ("rust-simd-adler32" ,rust-simd-adler32-0.3))))
    (home-page "https://github.com/Frommi/miniz_oxide/tree/master/miniz_oxide")
    (synopsis
     "DEFLATE compression and decompression library rewritten in Rust based on miniz")
    (description
     "This package provides DEFLATE compression and decompression library rewritten in Rust based on miniz.")
    (license (list license:expat license:zlib license:asl2.0))))

(define-public rust-zlib-rs-0.5
  (package
    (name "rust-zlib-rs")
    (version "0.5.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "zlib-rs" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1a1vssif5m2hwsy574l1gb668q4k04ggqv88yvr9mq29g66r52w6"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-arbitrary" ,rust-arbitrary-1)
                       ("rust-quickcheck" ,rust-quickcheck-1))))
    (home-page "https://github.com/trifectatechfoundation/zlib-rs")
    (synopsis "memory-safe zlib implementation written in rust")
    (description
     "This package provides a memory-safe zlib implementation written in rust.")
    (license license:zlib)))

(define-public rust-libz-rs-sys-0.5
  (package
    (name "rust-libz-rs-sys")
    (version "0.5.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "libz-rs-sys" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0nlc7cdcrh8mxqb08yw5i7ghgpcs1ixq4kk4sx19dzk0sydwm2b4"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-zlib-rs" ,rust-zlib-rs-0.5))))
    (home-page "https://github.com/trifectatechfoundation/zlib-rs")
    (synopsis "memory-safe zlib implementation written in rust")
    (description
     "This package provides a memory-safe zlib implementation written in rust.")
    (license license:zlib)))

(define-public rust-cloudflare-zlib-sys-0.3
  (package
    (name "rust-cloudflare-zlib-sys")
    (version "0.3.6")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "cloudflare-zlib-sys" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "10sr9nl83aripyn8hviacn5qifyl09g7rdxzxc9fv6pdfvb8bdrj"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-cc" ,rust-cc-1))))
    (home-page "https://lib.rs/crates/cloudflare-zlib-sys")
    (synopsis "Cloudflare fork of zlib with performance improvements")
    (description
     "This package provides Cloudflare fork of zlib with performance improvements.")
    (license license:zlib)))

(define-public rust-flate2-1
  (package
    (name "rust-flate2")
    (version "1.1.1")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "flate2" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1kpycx57dqpkr3vp53b4nq75p9mflh0smxy8hkys4v4ndvkr5vbw"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-cloudflare-zlib-sys" ,rust-cloudflare-zlib-sys-0.3)
                       ("rust-crc32fast" ,rust-crc32fast-1)
                       ("rust-libz-ng-sys" ,rust-libz-ng-sys-1)
                       ("rust-libz-rs-sys" ,rust-libz-rs-sys-0.5)
                       ("rust-libz-sys" ,rust-libz-sys-1)
                       ("rust-miniz-oxide" ,rust-miniz-oxide-0.8))))
    (home-page "https://github.com/rust-lang/flate2-rs")
    (synopsis
     "DEFLATE compression and decompression exposed as Read/BufRead/Write streams.
Supports miniz_oxide and multiple zlib implementations. Supports zlib, gzip,
and raw deflate streams.")
    (description
     "This package provides DEFLATE compression and decompression exposed as Read/@code{BufRead/Write}
streams.  Supports miniz_oxide and multiple zlib implementations.  Supports
zlib, gzip, and raw deflate streams.")
    (license (list license:expat license:asl2.0))))

(define-public rust-bytes-1
  (package
    (name "rust-bytes")
    (version "1.10.1")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "bytes" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "0smd4wi2yrhp5pmq571yiaqx84bjqlm1ixqhnvfwzzc6pqkn26yp"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-portable-atomic" ,rust-portable-atomic-1)
                       ("rust-serde" ,rust-serde-1))))
    (home-page "https://github.com/tokio-rs/bytes")
    (synopsis "Types and traits for working with bytes")
    (description
     "This package provides Types and traits for working with bytes.")
    (license license:expat)))

(define-public rust-noodles-bgzf-0.41
  (package
    (name "rust-noodles-bgzf")
    (version "0.41.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-bgzf" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1r5p1l6dv5v7djhridjsl142zshaxafz7kha2kzc203yxdsmszph"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-byteorder" ,rust-byteorder-1)
                       ("rust-bytes" ,rust-bytes-1)
                       ("rust-crossbeam-channel" ,rust-crossbeam-channel-0.5)
                       ("rust-flate2" ,rust-flate2-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-libdeflater" ,rust-libdeflater-1)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1)
                       ("rust-tokio-util" ,rust-tokio-util-0.7))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Blocked gzip format (BGZF) reader and writer")
    (description
     "This package provides Blocked gzip format (BGZF) reader and writer.")
    (license license:expat)))

(define-public rust-noodles-bam-0.81
  (package
    (name "rust-noodles-bam")
    (version "0.81.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles-bam" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "04a86a5yl65rvig9scr1y0m3ja1dd1qhcw0j1g62y35il9dqgpq5"))))
    (build-system cargo-build-system)
    (arguments
     `(#:skip-build? #t
       #:cargo-inputs (("rust-bstr" ,rust-bstr-1)
                       ("rust-byteorder" ,rust-byteorder-1)
                       ("rust-futures" ,rust-futures-0.3)
                       ("rust-indexmap" ,rust-indexmap-2)
                       ("rust-memchr" ,rust-memchr-2)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-sam" ,rust-noodles-sam-0.77)
                       ("rust-pin-project-lite" ,rust-pin-project-lite-0.2)
                       ("rust-tokio" ,rust-tokio-1))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Binary Alignment/Map (BAM) format reader and writer")
    (description
     "This package provides Binary Alignment/Map (BAM) format reader and writer.")
    (license license:expat)))

(define-public rust-noodles-0.99
  (package
    (name "rust-noodles")
    (version "0.99.0")
    (source
     (origin
       (method url-fetch)
       (uri (crate-uri "noodles" version))
       (file-name (string-append name "-" version ".tar.gz"))
       (sha256
        (base32 "1lpacp75pw22sl5mam83yari9hvd2x9gzis7194acl8mzj34l3hi"))))
    (build-system cargo-build-system)
    (arguments
     `(#:cargo-inputs (("rust-noodles-bam" ,rust-noodles-bam-0.81)
                       ("rust-noodles-bcf" ,rust-noodles-bcf-0.76)
                       ("rust-noodles-bed" ,rust-noodles-bed-0.26)
                       ("rust-noodles-bgzf" ,rust-noodles-bgzf-0.41)
                       ("rust-noodles-core" ,rust-noodles-core-0.17)
                       ("rust-noodles-cram" ,rust-noodles-cram-0.84)
                       ("rust-noodles-csi" ,rust-noodles-csi-0.49)
                       ("rust-noodles-fasta" ,rust-noodles-fasta-0.54)
                       ("rust-noodles-fastq" ,rust-noodles-fastq-0.19)
                       ("rust-noodles-gff" ,rust-noodles-gff-0.50)
                       ("rust-noodles-gtf" ,rust-noodles-gtf-0.45)
                       ("rust-noodles-htsget" ,rust-noodles-htsget-0.8)
                       ("rust-noodles-refget" ,rust-noodles-refget-0.7)
                       ("rust-noodles-sam" ,rust-noodles-sam-0.77)
                       ("rust-noodles-tabix" ,rust-noodles-tabix-0.55)
                       ("rust-noodles-vcf" ,rust-noodles-vcf-0.79))))
    (home-page "https://github.com/zaeleus/noodles")
    (synopsis "Bioinformatics I/O libraries")
    (description "This package provides Bioinformatics I/O libraries.")
    (license license:expat)))
