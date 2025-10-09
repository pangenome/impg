;;; IMplicit Pangenome Graph (IMPG)
;;; Copyright © 2025 Frederick M. Muriithi <fredmanglis@gmail.com>

(define-module (impg rust-crates)
  #:use-module (guix packages)
  #:use-module (guix git-download)
  #:use-module (guix build-system cargo)

  #:use-module (gnu packages)

  #:export (lookup-cargo-inputs))

(define rust-adler2-2.0.1
  (crate-source "adler2" "2.0.1"
                "1ymy18s9hs7ya1pjc9864l30wk8p2qfqdi7mhhcc5nfakxbij09j"))

;; (define rust-agc-rs-0.1.0.aba33a3
;;   ;; TODO: Define standalone package if this is a workspace.
;;   (origin
;;     (method git-fetch)
;;     (uri (git-reference (url "https://github.com/pangenome/agc-rs.git")
;;                         (commit "aba33a36a5001dfac764749aacf6af0a770eb69b")))
;;     (file-name (git-file-name "rust-agc-rs" "0.1.0.aba33a3"))
;;     (sha256 (base32 "10br4cky2slnfida1wghha6ba58ghs41qmpqail9fg8hgxhsl6m3"))))

(define rust-aho-corasick-1.1.3
  (crate-source "aho-corasick" "1.1.3"
                "05mrpkvdgp5d20y2p989f187ry9diliijgwrs254fs9s1m1x6q4f"))

(define rust-anstream-0.6.20
  (crate-source "anstream" "0.6.20"
                "14k1iqdf3dx7hdjllmql0j9sjxkwr1lfdddi3adzff0r7mjn7r9s"))

(define rust-anstyle-1.0.11
  (crate-source "anstyle" "1.0.11"
                "1gbbzi0zbgff405q14v8hhpi1kz2drzl9a75r3qhks47lindjbl6"))

(define rust-anstyle-parse-0.2.7
  (crate-source "anstyle-parse" "0.2.7"
                "1hhmkkfr95d462b3zf6yl2vfzdqfy5726ya572wwg8ha9y148xjf"))

(define rust-anstyle-query-1.1.4
  (crate-source "anstyle-query" "1.1.4"
                "1qir6d6fl5a4y2gmmw9a5w93ckwx6xn51aryd83p26zn6ihiy8wy"))

(define rust-anstyle-wincon-3.0.10
  (crate-source "anstyle-wincon" "3.0.10"
                "0ajz9wsf46a2l3pds7v62xbhq2cffj7wrilamkx2z8r28m0k61iy"))

(define rust-anyhow-1.0.99
  (crate-source "anyhow" "1.0.99"
                "001icqvkfl28rxxmk99rm4gvdzxqngj5v50yg2bh3dzcvqfllrxh"))

(define rust-approx-0.5.1
  (crate-source "approx" "0.5.1"
                "1ilpv3dgd58rasslss0labarq7jawxmivk17wsh8wmkdm3q15cfa"))

(define rust-arrayvec-0.5.2
  (crate-source "arrayvec" "0.5.2"
                "12q6hn01x5435bprwlb7w9m7817dyfq55yrl4psygr78bp32zdi3"))

(define rust-autocfg-1.5.0
  (crate-source "autocfg" "1.5.0"
                "1s77f98id9l4af4alklmzq46f21c980v13z2r1pcxx6bqgw0d1n0"))

(define rust-bgzip-0.3.1
  (crate-source "bgzip" "0.3.1"
                "16zr2nclis3kgz0jxi7ayyk510ar5dvyfpf03fazajmn1ycdhkxn"))

(define rust-bincode-2.0.1
  (crate-source "bincode" "2.0.1"
                "0h5pxp3dqkigjwy926a03sl69n9wv7aq4142a20kw9lhn3bzbsin"))

(define rust-bincode-derive-2.0.1
  (crate-source "bincode_derive" "2.0.1"
                "029wmh26hq3hhs1gq629y0frn2pkl7ld061rk23fji8g8jd715dz"))

(define rust-bio-types-1.0.4
  (crate-source "bio-types" "1.0.4"
                "0zmdcvj44a088larkahcic5z61cwn2x80iym0w14albzid7zbp7l"))

(define rust-bitflags-1.3.2
  (crate-source "bitflags" "1.3.2"
                "12ki6w8gn1ldq7yz9y680llwk5gmrhrzszaa17g1sbrw2r2qvwxy"))

(define rust-bitflags-2.9.4
  (crate-source "bitflags" "2.9.4"
                "157kkcv8s7vk6d17dar1pa5cqcz4c8pdrn16wm1ld7jnr86d2q92"))

(define rust-bitvec-1.0.1
  (crate-source "bitvec" "1.0.1"
                "173ydyj2q5vwj88k6xgjnfsshs4x9wbvjjv7sm0h36r34hn87hhv"))

(define rust-boomphf-0.5.9
  (crate-source "boomphf" "0.5.9"
                "0braniw72g9yq5006sfgc1g8d4317bb524c694jw6nggizrvg3sf"))

(define rust-bstr-0.2.17
  (crate-source "bstr" "0.2.17"
                "08rjbhysy6gg27db2h3pnhvr2mlr5vkj797i9625kwg8hgrnjdds"))

(define rust-bumpalo-3.19.0
  (crate-source "bumpalo" "3.19.0"
                "0hsdndvcpqbjb85ghrhska2qxvp9i75q2vb70hma9fxqawdy9ia6"))

(define rust-bytemuck-1.23.2
  (crate-source "bytemuck" "1.23.2"
                "0xs637lsr9p73ackjkmbjw80dp1dfdw0ydhdk0gzjcnzpkpfm59r"))

(define rust-bytemuck-derive-1.10.1
  (crate-source "bytemuck_derive" "1.10.1"
                "0a9dczfzn2c1lgg7afhqrh2drmg34w49hxhipni6pjri49blw5ag"))

(define rust-byteorder-1.5.0
  (crate-source "byteorder" "1.5.0"
                "0jzncxyf404mwqdbspihyzpkndfgda450l0893pz5xj685cg5l0z"))

(define rust-bytes-1.10.1
  (crate-source "bytes" "1.10.1"
                "0smd4wi2yrhp5pmq571yiaqx84bjqlm1ixqhnvfwzzc6pqkn26yp"))

(define rust-bzip2-0.5.2
  (crate-source "bzip2" "0.5.2"
                "0iya6nbj0p2y8jss0z05yncc5hadry164fw3zva01y06v4igpv29"))

(define rust-bzip2-sys-0.1.13+1.0.8
  ;; TODO: Check bundled sources.
  (crate-source "bzip2-sys" "0.1.13+1.0.8"
                "056c39pgjh4272bdslv445f5ry64xvb0f7nph3z7860ln8rzynr2"))

(define rust-cc-1.2.37
  (crate-source "cc" "1.2.37"
                "0i5xlxsgd7jif1ry9k3ysnpsmbrckapqwq8d8l5vhkj0qs4ka6b5"))

(define rust-cfg-if-1.0.3
  (crate-source "cfg-if" "1.0.3"
                "1afg7146gbxjvkbjx7i5sdrpqp9q5akmk9004fr8rsm90jf2il9g"))

(define rust-clap-4.5.47
  (crate-source "clap" "4.5.47"
                "0c99f6m4a7d4ffgahid49h0ci2pv4ccdf417f76nl4wx5n801b3y"))

(define rust-clap-builder-4.5.47
  (crate-source "clap_builder" "4.5.47"
                "1mp1f0fz6wp9v87jc9372lg6r4514ja1l8cazf25hfz7a3vvpn9a"))

(define rust-clap-derive-4.5.47
  (crate-source "clap_derive" "4.5.47"
                "174z9g13s85la2nmi8gv8ssjwz77im3rqg5isiinw6hg1fp7xzdv"))

(define rust-clap-lex-0.7.5
  (crate-source "clap_lex" "0.7.5"
                "0xb6pjza43irrl99axbhs12pxq4sr8x7xd36p703j57f5i3n2kxr"))

(define rust-cmake-0.1.54
  (crate-source "cmake" "0.1.54"
                "1w41ma28yzad9x757s9sfq1wigjs9j902hbzc0nbxpc9vvws7jp7"))

(define rust-codespan-reporting-0.12.0
  (crate-source "codespan-reporting" "0.12.0"
                "108g41xqzhr8fx8hlpy5qzmqq8ylldbj37wndkaqm34yy1d2wvgy"))

(define rust-coitrees-0.4.0
  (crate-source "coitrees" "0.4.0"
                "1qwb4c5gx30gl1kyi85rbq6z23l2f9lm0q02ym160n0fvc89c3r4"))

(define rust-colorchoice-1.0.4
  (crate-source "colorchoice" "1.0.4"
                "0x8ymkz1xr77rcj1cfanhf416pc4v681gmkc9dzb3jqja7f62nxh"))

(define rust-console-0.16.1
  (crate-source "console" "0.16.1"
                "1x4x6vfi1s55nbr4i77b9r87s213h46lq396sij9fkmidqx78c5l"))

(define rust-core-affinity-0.8.3
  (crate-source "core_affinity" "0.8.3"
                "0hhkjybngi5n2ayjmbba2n2gh9fc8xbqgpzm2dp6q094nskv6d50"))

(define rust-crc32fast-1.5.0
  (crate-source "crc32fast" "1.5.0"
                "04d51liy8rbssra92p0qnwjw8i9rm9c4m3bwy19wjamz1k4w30cl"))

(define rust-crossbeam-channel-0.5.15
  (crate-source "crossbeam-channel" "0.5.15"
                "1cicd9ins0fkpfgvz9vhz3m9rpkh6n8d3437c3wnfsdkd3wgif42"))

(define rust-crossbeam-deque-0.8.6
  (crate-source "crossbeam-deque" "0.8.6"
                "0l9f1saqp1gn5qy0rxvkmz4m6n7fc0b3dbm6q1r5pmgpnyvi3lcx"))

(define rust-crossbeam-epoch-0.9.18
  (crate-source "crossbeam-epoch" "0.9.18"
                "03j2np8llwf376m3fxqx859mgp9f83hj1w34153c7a9c7i5ar0jv"))

(define rust-crossbeam-utils-0.8.21
  (crate-source "crossbeam-utils" "0.8.21"
                "0a3aa2bmc8q35fb67432w16wvi54sfmb69rk9h5bhd18vw0c99fh"))

(define rust-custom-derive-0.1.7
  (crate-source "custom_derive" "0.1.7"
                "1f81bavw1wnykwh21hh4yyzigs6zl6f6pkk9p3car8kq95yfb2pg"))

(define rust-cxx-1.0.184
  (crate-source "cxx" "1.0.184"
                "0s7ksxsbb6qigcsvlqccglfs7jg44gjlw257dbgd084x6vmhnjmy"))

(define rust-cxx-build-1.0.184
  (crate-source "cxx-build" "1.0.184"
                "02m4w407qlkq5hi72i8sah0svdfq43lhmprlpiff9a2n7swmbn97"))

(define rust-cxxbridge-cmd-1.0.184
  (crate-source "cxxbridge-cmd" "1.0.184"
                "0l01xis2b1d5dcn4x77vyfhbxcj3g0nwp2r2icmxrawxim36qbq5"))

(define rust-cxxbridge-flags-1.0.184
  (crate-source "cxxbridge-flags" "1.0.184"
                "0i31rrp2sbgns7f3lcpsspf2yy1c7g81fcn60a0cp1h933x4bl8g"))

(define rust-cxxbridge-macro-1.0.184
  (crate-source "cxxbridge-macro" "1.0.184"
                "0q1wpklc8yq8qw1gx8l2jsz24rfjidccj8c41am2sjj8qhxlmb02"))

(define rust-derive-new-0.5.9
  (crate-source "derive-new" "0.5.9"
                "0d9m5kcj1rdmdjqfgj7rxxhdzx0as7p4rp1mjx5j6w5dl2f3461l"))

(define rust-derive-new-0.6.0
  (crate-source "derive-new" "0.6.0"
                "1b8jv6jx0b8jgkz9kmz0ciqmnf74xkk0mmvkb5z1c87932kdwl6i"))

(define rust-displaydoc-0.2.5
  (crate-source "displaydoc" "0.2.5"
                "1q0alair462j21iiqwrr21iabkfnb13d6x5w95lkdg21q2xrqdlp"))

(define rust-either-1.15.0
  (crate-source "either" "1.15.0"
                "069p1fknsmzn9llaizh77kip0pqmcwpdsykv2x30xpjyija5gis8"))

(define rust-encode-unicode-1.0.0
  (crate-source "encode_unicode" "1.0.0"
                "1h5j7j7byi289by63s3w4a8b3g6l5ccdrws7a67nn07vdxj77ail"))

(define rust-env-filter-0.1.3
  (crate-source "env_filter" "0.1.3"
                "1l4p6f845cylripc3zkxa0lklk8rn2q86fqm522p6l2cknjhavhq"))

(define rust-env-logger-0.11.8
  (crate-source "env_logger" "0.11.8"
                "17q6zbjam4wq75fa3m4gvvmv3rj3ch25abwbm84b28a0j3q67j0k"))

(define rust-equivalent-1.0.2
  (crate-source "equivalent" "1.0.2"
                "03swzqznragy8n0x31lqc78g2af054jwivp7lkrbrc0khz74lyl7"))

(define rust-errno-0.3.14
  (crate-source "errno" "0.3.14"
                "1szgccmh8vgryqyadg8xd58mnwwicf39zmin3bsn63df2wbbgjir"))

(define rust-fastrand-2.3.0
  (crate-source "fastrand" "2.3.0"
                "1ghiahsw1jd68df895cy5h3gzwk30hndidn3b682zmshpgmrx41p"))

(define rust-find-msvc-tools-0.1.1
  (crate-source "find-msvc-tools" "0.1.1"
                "0b8rhghgjssjw9q8a3gg7f9kl8zhy9d7nqsc4s4nc52dyqq9knbz"))

(define rust-flate2-1.1.2
  (crate-source "flate2" "1.1.2"
                "07abz7v50lkdr5fjw8zaw2v8gm2vbppc0f7nqm8x3v3gb6wpsgaa"))

(define rust-flume-0.11.1
  (crate-source "flume" "0.11.1"
                "15ch0slxa8sqsi6c73a0ky6vdnh48q8cxjf7rksa3243m394s3ns"))

(define rust-fnv-1.0.7
  (crate-source "fnv" "1.0.7"
                "1hc2mcqha06aibcaza94vbi81j6pr9a1bbxrxjfhc91zin8yr7iz"))

(define rust-foldhash-0.2.0
  (crate-source "foldhash" "0.2.0"
                "1nvgylb099s11xpfm1kn2wcsql080nqmnhj1l25bp3r2b35j9kkp"))

(define rust-form-urlencoded-1.2.2
  (crate-source "form_urlencoded" "1.2.2"
                "1kqzb2qn608rxl3dws04zahcklpplkd5r1vpabwga5l50d2v4k6b"))

(define rust-fs-utils-1.1.4
  (crate-source "fs-utils" "1.1.4"
                "14r5wl14mz227v0lpy89lvjzfnxgdxigvrrmm6c4r52w03fakivg"))

(define rust-funty-2.0.0
  (crate-source "funty" "2.0.0"
                "177w048bm0046qlzvp33ag3ghqkqw4ncpzcm5lq36gxf2lla7mg6"))

(define rust-futures-core-0.3.31
  (crate-source "futures-core" "0.3.31"
                "0gk6yrxgi5ihfanm2y431jadrll00n5ifhnpx090c2f2q1cr1wh5"))

(define rust-futures-sink-0.3.31
  (crate-source "futures-sink" "0.3.31"
                "1xyly6naq6aqm52d5rh236snm08kw8zadydwqz8bip70s6vzlxg5"))

(define rust-getrandom-0.2.16
  (crate-source "getrandom" "0.2.16"
                "14l5aaia20cc6cc08xdlhrzmfcylmrnprwnna20lqf746pqzjprk"))

(define rust-getrandom-0.3.3
  (crate-source "getrandom" "0.3.3"
                "1x6jl875zp6b2b6qp9ghc84b0l76bvng2lvm8zfcmwjl7rb5w516"))

(define rust-gfa-0.10.1
  (crate-source "gfa" "0.10.1"
                "1x996rpfnflgi2j4dgaj5sdxdbf24zfm9d2ha0zy8aid0cd60cln"))

(define rust-glob-0.3.3
  (crate-source "glob" "0.3.3"
                "106jpd3syfzjfj2k70mwm0v436qbx96wig98m4q8x071yrq35hhc"))

(define rust-gzp-1.0.1
  (crate-source "gzp" "1.0.1"
                "0k9qhky0vm4kyqqqi8i8h99128mlfmvl9w53v9kgm9nql3lq18gc"))

(define rust-handlegraph-0.7.0-alpha.9.3ac575e
  ;; TODO: Define standalone package if this is a workspace.
  (origin
    (method git-fetch)
    (uri (git-reference (url "https://github.com/chfi/rs-handlegraph")
                        (commit "3ac575e4216ce16a16667503a8875e469a40a97a")))
    (file-name (git-file-name "rust-handlegraph" "0.7.0-alpha.9.3ac575e"))
    (sha256 (base32 "1x9lhc4hjyfixvhdxr6z0lanfcynnqsmx3dqaf6xw4dpx0i4mcgg"))))

(define rust-hashbrown-0.15.5
  (crate-source "hashbrown" "0.15.5"
                "189qaczmjxnikm9db748xyhiw04kpmhm9xj9k9hg0sgx7pjwyacj"))

(define rust-heck-0.5.0
  (crate-source "heck" "0.5.0"
                "1sjmpsdl8czyh9ywl3qcsfsq9a307dg4ni2vnlwgnzzqhc4y0113"))

(define rust-hermit-abi-0.5.2
  (crate-source "hermit-abi" "0.5.2"
                "1744vaqkczpwncfy960j2hxrbjl1q01csm84jpd9dajbdr2yy3zw"))

(define rust-hts-sys-2.2.0
  ;; TODO: Check bundled sources.
  (crate-source "hts-sys" "2.2.0"
                "1cmvdwssd6xjk6w1iigaj5rl9ibx4zaaskfb2ji2mlhw28f7z3g3"))

(define rust-icu-collections-2.0.0
  (crate-source "icu_collections" "2.0.0"
                "0izfgypv1hsxlz1h8fc2aak641iyvkak16aaz5b4aqg3s3sp4010"))

(define rust-icu-locale-core-2.0.0
  (crate-source "icu_locale_core" "2.0.0"
                "02phv7vwhyx6vmaqgwkh2p4kc2kciykv2px6g4h8glxfrh02gphc"))

(define rust-icu-normalizer-2.0.0
  (crate-source "icu_normalizer" "2.0.0"
                "0ybrnfnxx4sf09gsrxri8p48qifn54il6n3dq2xxgx4dw7l80s23"))

(define rust-icu-normalizer-data-2.0.0
  (crate-source "icu_normalizer_data" "2.0.0"
                "1lvjpzxndyhhjyzd1f6vi961gvzhj244nribfpdqxjdgjdl0s880"))

(define rust-icu-properties-2.0.1
  (crate-source "icu_properties" "2.0.1"
                "0az349pjg8f18lrjbdmxcpg676a7iz2ibc09d2wfz57b3sf62v01"))

(define rust-icu-properties-data-2.0.1
  (crate-source "icu_properties_data" "2.0.1"
                "0cnn3fkq6k88w7p86w7hsd1254s4sl783rpz4p6hlccq74a5k119"))

(define rust-icu-provider-2.0.0
  (crate-source "icu_provider" "2.0.0"
                "1bz5v02gxv1i06yhdhs2kbwxkw3ny9r2vvj9j288fhazgfi0vj03"))

(define rust-idna-1.1.0
  (crate-source "idna" "1.1.0"
                "1pp4n7hppm480zcx411dsv9wfibai00wbpgnjj4qj0xa7kr7a21v"))

(define rust-idna-adapter-1.2.1
  (crate-source "idna_adapter" "1.2.1"
                "0i0339pxig6mv786nkqcxnwqa87v4m94b2653f6k3aj0jmhfkjis"))

(define rust-ieee754-0.2.6
  (crate-source "ieee754" "0.2.6"
                "1771d2kvw1wga65yrg9m7maky0fzsaq9hvhkv91n6gmxmjfdl1wh"))

(define rust-indexmap-2.11.3
  (crate-source "indexmap" "2.11.3"
                "1hqs931f1sd3r92zj77ji9bs75f20amnj0s3aqas9zqkym29h4cj"))

(define rust-indicatif-0.18.0
  (crate-source "indicatif" "0.18.0"
                "1kg1wi3x9x15f22q99spfzcg5fzlmhcc5i6aqjxyssyh8vcld9kh"))

(define rust-is-terminal-polyfill-1.70.1
  (crate-source "is_terminal_polyfill" "1.70.1"
                "1kwfgglh91z33kl0w5i338mfpa3zs0hidq5j4ny4rmjwrikchhvr"))

(define rust-jiff-0.2.15
  (crate-source "jiff" "0.2.15"
                "0jby6kbs2ra33ji0rx4swcp66jzmcvgszc5v4izwfsgbn6w967xy"))

(define rust-jiff-static-0.2.15
  (crate-source "jiff-static" "0.2.15"
                "1d4l4pvlhz3w487gyhnzvagpbparspv4c8f35qk6g5w9zx8k8d03"))

(define rust-jobserver-0.1.34
  (crate-source "jobserver" "0.1.34"
                "0cwx0fllqzdycqn4d6nb277qx5qwnmjdxdl0lxkkwssx77j3vyws"))

(define rust-js-sys-0.3.79
  ;; TODO: Check bundled sources.
  (crate-source "js-sys" "0.3.79"
                "1l4qhb85znvz29ln5znng7wjjapwqpy4grw6l5rlxbaqhs5xliv2"))

(define rust-lazy-static-1.5.0
  (crate-source "lazy_static" "1.5.0"
                "1zk6dqqni0193xg6iijh7i3i44sryglwgvx20spdvwk3r6sbrlmv"))

(define rust-lexical-core-0.7.6
  (crate-source "lexical-core" "0.7.6"
                "1zjzab1fnaw4kj6ixyrskp4dyz761gdcab07m4bkvlk1l4mcc1v6"))

(define rust-libc-0.2.175
  (crate-source "libc" "0.2.175"
                "0hw5sb3gjr0ivah7s3fmavlpvspjpd4mr009abmam2sr7r4sx0ka"))

(define rust-libdeflate-sys-1.24.0
  ;; TODO: Check bundled sources.
  (crate-source "libdeflate-sys" "1.24.0"
                "0p41fqvblzdwhmkv96jzsziv7my7wc7qaqpbvyclbi36acr28n40"))

(define rust-libdeflater-1.24.0
  (crate-source "libdeflater" "1.24.0"
                "1frhpcm4kh0vaddn1s44j6m6dya40hdmp984lmkykp6nx73vqw5j"))

(define rust-liblzma-0.3.6
  (crate-source "liblzma" "0.3.6"
                "0r6pkykpajdypdyyij90d8s2ihhsz9m9ly7pm1dpfsg29frd4cd6"))

(define rust-liblzma-sys-0.3.13
  ;; TODO: Check bundled sources.
  (crate-source "liblzma-sys" "0.3.13"
                "0x9lni7a3x1rwdsribj311zpxb5n99kn256yad2z7vxck4ddznpg"))

(define rust-libz-ng-sys-1.1.22
  ;; TODO: Check bundled sources.
  (crate-source "libz-ng-sys" "1.1.22"
                "096qkwzy596zf88nfppr2vbw9fbjfr81k2ws4zf6wyrw58n8q4d7"))

(define rust-libz-rs-sys-0.5.2
  ;; TODO: Check bundled sources.
  (crate-source "libz-rs-sys" "0.5.2"
                "1kdy093bhxfkgx7li3raxigcc3qdqjn3hvrpjkblvv6r777vh3c4"))

(define rust-libz-sys-1.1.22
  ;; TODO: Check bundled sources.
  (crate-source "libz-sys" "1.1.22"
                "07b5wxh0ska996kc0g2hanjhmb4di7ksm6ndljhr4pi0vykyfw4b"))

(define rust-linear-map-1.2.0
  (crate-source "linear-map" "1.2.0"
                "1vh3sczl4xb5asdlpafdf3y4g9bp63fgs8y2a2sjgmcsn7v21bmz"))

(define rust-link-cplusplus-1.0.12
  (crate-source "link-cplusplus" "1.0.12"
                "10lcgfp9pnxpihp21s86xnq57vpr97m2k419d8rvkl57m8qcfy3z"))

(define rust-linux-raw-sys-0.11.0
  ;; TODO: Check bundled sources.
  (crate-source "linux-raw-sys" "0.11.0"
                "0fghx0nn8nvbz5yzgizfcwd6ap2pislp68j8c1bwyr6sacxkq7fz"))

(define rust-litemap-0.8.0
  (crate-source "litemap" "0.8.0"
                "0mlrlskwwhirxk3wsz9psh6nxcy491n0dh8zl02qgj0jzpssw7i4"))

(define rust-lock-api-0.4.13
  (crate-source "lock_api" "0.4.13"
                "0rd73p4299mjwl4hhlfj9qr88v3r0kc8s1nszkfmnq2ky43nb4wn"))

(define rust-log-0.4.28
  (crate-source "log" "0.4.28"
                "0cklpzrpxafbaq1nyxarhnmcw9z3xcjrad3ch55mmr58xw2ha21l"))

(define rust-matrixmultiply-0.3.10
  (crate-source "matrixmultiply" "0.3.10"
                "020sqwg3cvprfasbszqbnis9zx6c3w9vlkfidyimgblzdq0y6vd0"))

(define rust-memchr-2.7.5
  (crate-source "memchr" "2.7.5"
                "1h2bh2jajkizz04fh047lpid5wgw2cr9igpkdhl3ibzscpd858ij"))

(define rust-memmap-0.7.0
  (crate-source "memmap" "0.7.0"
                "0ns7kkd1h4pijdkwfvw4qlbbmqmlmzwlq3g2676dcl5vwyazv1b5"))

(define rust-miniz-oxide-0.8.9
  (crate-source "miniz_oxide" "0.8.9"
                "05k3pdg8bjjzayq3rf0qhpirq9k37pxnasfn4arbs17phqn6m9qz"))

(define rust-nalgebra-0.33.2
  (crate-source "nalgebra" "0.33.2"
                "0fvayv2fa6x4mfm4cq3m2cfcc2jwkiq4sm73209zszkh9gvcvbi6"))

(define rust-nalgebra-macros-0.2.2
  (crate-source "nalgebra-macros" "0.2.2"
                "1z6v9phhr1hwzyyblf792128lxfv1hy1sxl4cvikihcgmxr56ji5"))

(define rust-nanorand-0.7.0
  (crate-source "nanorand" "0.7.0"
                "1hr60b8zlfy7mxjcwx2wfmhpkx7vfr3v9x12shmv1c10b0y32lba"))

(define rust-natord-1.0.9
  (crate-source "natord" "1.0.9"
                "0z75spwag3ch20841pvfwhh3892i2z2sli4pzp1jgizbipdrd39h"))

(define rust-newtype-derive-0.1.6
  (crate-source "newtype_derive" "0.1.6"
                "1v3170xscs65gjx5vl1zjnqp86wngbzw3n2q74ibfnqqkx6x535c"))

(define rust-niffler-3.0.0
  (crate-source "niffler" "3.0.0"
                "0x1mzgfhpxr0mwwpsrmlkyalmbaiv97pspyjvymrzb1xr5f13lv2"))

(define rust-nom-5.1.3
  (crate-source "nom" "5.1.3"
                "0jyxc4d3pih60pp8hvzpg5ajh16s273cpnsdpzp04qv7g8w9m588"))

(define rust-noodles-0.100.0
  (crate-source "noodles" "0.100.0"
                "17lnhmzbp94g383sxxwmqd0aag9daa1vmgsqx7p9gzdkwwbrqpca"))

(define rust-noodles-bgzf-0.42.0
  (crate-source "noodles-bgzf" "0.42.0"
                "0fdllcmsdyqg6zays6y0s3lls1qjjdm5jkhh832x1by434zkw73w"))

(define rust-num-bigint-0.4.6
  (crate-source "num-bigint" "0.4.6"
                "1f903zd33i6hkjpsgwhqwi2wffnvkxbn6rv4mkgcjcqi7xr4zr55"))

(define rust-num-complex-0.4.6
  (crate-source "num-complex" "0.4.6"
                "15cla16mnw12xzf5g041nxbjjm9m85hdgadd5dl5d0b30w9qmy3k"))

(define rust-num-cpus-1.17.0
  (crate-source "num_cpus" "1.17.0"
                "0fxjazlng4z8cgbmsvbzv411wrg7x3hyxdq8nxixgzjswyylppwi"))

(define rust-num-integer-0.1.46
  (crate-source "num-integer" "0.1.46"
                "13w5g54a9184cqlbsq80rnxw4jj4s0d8wv75jsq5r2lms8gncsbr"))

(define rust-num-rational-0.4.2
  (crate-source "num-rational" "0.4.2"
                "093qndy02817vpgcqjnj139im3jl7vkq4h68kykdqqh577d18ggq"))

(define rust-num-traits-0.2.19
  (crate-source "num-traits" "0.2.19"
                "0h984rhdkkqd4ny9cif7y2azl3xdfb7768hb9irhpsch4q3gq787"))

(define rust-once-cell-1.21.3
  (crate-source "once_cell" "1.21.3"
                "0b9x77lb9f1j6nqgf5aka4s2qj0nly176bpbrv6f9iakk5ff3xa2"))

(define rust-once-cell-polyfill-1.70.1
  (crate-source "once_cell_polyfill" "1.70.1"
                "1bg0w99srq8h4mkl68l1mza2n2f2hvrg0n8vfa3izjr5nism32d4"))

(define rust-paste-1.0.15
  (crate-source "paste" "1.0.15"
                "02pxffpdqkapy292harq6asfjvadgp1s005fip9ljfsn9fvxgh2p"))

(define rust-percent-encoding-2.3.2
  (crate-source "percent-encoding" "2.3.2"
                "083jv1ai930azvawz2khv7w73xh8mnylk7i578cifndjn5y64kwv"))

(define rust-pkg-config-0.3.32
  (crate-source "pkg-config" "0.3.32"
                "0k4h3gnzs94sjb2ix6jyksacs52cf1fanpwsmlhjnwrdnp8dppby"))

(define rust-portable-atomic-1.11.1
  (crate-source "portable-atomic" "1.11.1"
                "10s4cx9y3jvw0idip09ar52s2kymq8rq9a668f793shn1ar6fhpq"))

(define rust-portable-atomic-util-0.2.4
  (crate-source "portable-atomic-util" "0.2.4"
                "01rmx1li07ixsx3sqg2bxqrkzk7b5n8pibwwf2589ms0s3cg18nq"))

(define rust-potential-utf-0.1.3
  (crate-source "potential_utf" "0.1.3"
                "12mhwvhpvvim6xqp6ifgkh1sniv9j2cmid6axn10fnjvpsnikpw4"))

(define rust-proc-macro2-1.0.101
  (crate-source "proc-macro2" "1.0.101"
                "1pijhychkpl7rcyf1h7mfk6gjfii1ywf5n0snmnqs5g4hvyl7bl9"))

(define rust-quick-error-1.2.3
  (crate-source "quick-error" "1.2.3"
                "1q6za3v78hsspisc197bg3g7rpc989qycy8ypr8ap8igv10ikl51"))

(define rust-quote-1.0.40
  (crate-source "quote" "1.0.40"
                "1394cxjg6nwld82pzp2d4fp6pmaz32gai1zh9z5hvh0dawww118q"))

(define rust-r-efi-5.3.0
  (crate-source "r-efi" "5.3.0"
                "03sbfm3g7myvzyylff6qaxk4z6fy76yv860yy66jiswc2m6b7kb9"))

(define rust-radium-0.7.0
  (crate-source "radium" "0.7.0"
                "02cxfi3ky3c4yhyqx9axqwhyaca804ws46nn4gc1imbk94nzycyw"))

(define rust-rand-core-0.6.4
  (crate-source "rand_core" "0.6.4"
                "0b4j2v4cb5krak1pv6kakv4sz6xcwbrmy2zckc32hsigbrwy82zc"))

(define rust-rawpointer-0.2.1
  (crate-source "rawpointer" "0.2.1"
                "1qy1qvj17yh957vhffnq6agq0brvylw27xgks171qrah75wmg8v0"))

(define rust-rayon-1.11.0
  (crate-source "rayon" "1.11.0"
                "13x5fxb7rn4j2yw0cr26n7782jkc7rjzmdkg42qxk3xz0p8033rn"))

(define rust-rayon-core-1.13.0
  (crate-source "rayon-core" "1.13.0"
                "14dbr0sq83a6lf1rfjq5xdpk5r6zgzvmzs5j6110vlv2007qpq92"))

(define rust-regex-1.11.2
  (crate-source "regex" "1.11.2"
                "04k9rzxd11hcahpyihlswy6f1zqw7lspirv4imm4h0lcdl8gvmr3"))

(define rust-regex-automata-0.1.10
  (crate-source "regex-automata" "0.1.10"
                "0ci1hvbzhrfby5fdpf4ganhf7kla58acad9i1ff1p34dzdrhs8vc"))

(define rust-regex-automata-0.4.10
  (crate-source "regex-automata" "0.4.10"
                "1mllcfmgjcl6d52d5k09lwwq9wj5mwxccix4bhmw5spy1gx5i53b"))

(define rust-regex-syntax-0.8.6
  (crate-source "regex-syntax" "0.8.6"
                "00chjpglclfskmc919fj5aq308ffbrmcn7kzbkz92k231xdsmx6a"))

(define rust-rust-htslib-0.46.0
  (crate-source "rust-htslib" "0.46.0"
                "13lnc33n5699v9cg1bgnfflfapfi9la9k37zfnpb9gh18v5gkimf"))

(define rust-rustc-hash-2.1.1
  (crate-source "rustc-hash" "2.1.1"
                "03gz5lvd9ghcwsal022cgkq67dmimcgdjghfb5yb5d352ga06xrm"))

(define rust-rustc-version-0.1.7
  (crate-source "rustc_version" "0.1.7"
                "1160jjsqhqr25cvhr48hmpp8v61bjvjcnxzb0cyf4373lmp3gxf5"))

(define rust-rustix-1.1.2
  (crate-source "rustix" "1.1.2"
                "0gpz343xfzx16x82s1x336n0kr49j02cvhgxdvaq86jmqnigh5fd"))

(define rust-rustversion-1.0.22
  (crate-source "rustversion" "1.0.22"
                "0vfl70jhv72scd9rfqgr2n11m5i9l1acnk684m2w83w0zbqdx75k"))

(define rust-ryu-1.0.20
  (crate-source "ryu" "1.0.20"
                "07s855l8sb333h6bpn24pka5sp7hjk2w667xy6a0khkf6sqv5lr8"))

(define rust-safe-arch-0.7.4
  (crate-source "safe_arch" "0.7.4"
                "08sk47n1kcm5w2di6bpgi2hsw8r2caz2230pwqvbdqfv5pl2vc4n"))

(define rust-scopeguard-1.2.0
  (crate-source "scopeguard" "1.2.0"
                "0jcz9sd47zlsgcnm1hdw0664krxwb5gczlif4qngj2aif8vky54l"))

(define rust-scratch-1.0.9
  (crate-source "scratch" "1.0.9"
                "1cj826qggwn482wbfnzij5g9p411qszai0dnfld4qzh93g2jx3yn"))

(define rust-semver-0.1.20
  (crate-source "semver" "0.1.20"
                "1b10m0hxrr947gp41lj9vnmgl5ddwx3d41vnblsg06ppvkz11x6l"))

(define rust-serde-1.0.225
  (crate-source "serde" "1.0.225"
                "07dxpjh0g1mq3md9yvn7jbgssgcizcircf23f04xml1mwbg28v7x"))

(define rust-serde-core-1.0.225
  (crate-source "serde_core" "1.0.225"
                "10v3z58j5k6xhdxh90xgrv20wlnz5fnl67n04jdm47nbl3wmd4v5"))

(define rust-serde-derive-1.0.225
  (crate-source "serde_derive" "1.0.225"
                "05j5zj2jdba3jnm7kh3fpljmhngmsa8pp5x495lpc7wbyynkda8f"))

(define rust-shlex-1.3.0
  (crate-source "shlex" "1.3.0"
                "0r1y6bv26c1scpxvhg2cabimrmwgbp4p3wy6syj9n0c4s3q2znhg"))

(define rust-simba-0.9.1
  (crate-source "simba" "0.9.1"
                "15gxgwcm6vs2wbbc5z4x8zsi1rhjl3nvqnxpl95hjrhnnaz894n9"))

(define rust-smallvec-1.15.1
  (crate-source "smallvec" "1.15.1"
                "00xxdxxpgyq5vjnpljvkmy99xij5rxgh913ii1v16kzynnivgcb7"))

(define rust-spin-0.9.8
  (crate-source "spin" "0.9.8"
                "0rvam5r0p3a6qhc18scqpvpgb3ckzyqxpgdfyjnghh8ja7byi039"))

(define rust-spoa-rs-0.1.0.6f4f102
  ;; TODO: Define standalone package if this is a workspace.
  (origin
    (method git-fetch)
    (uri (git-reference (url "https://github.com/AndreaGuarracino/spoa-rs.git")
                        (commit "6f4f1024cc09959b926d4616991672693013ed3e")))
    (file-name (git-file-name "rust-spoa-rs" "0.1.0.6f4f102"))
    (sha256 (base32 "13hmxqk70c0728v8l44qap25qla4w86im3kwfmay7jwwg9lk7lfg"))
    (patches
     (search-patches "impg/patches/do-not-build-cplusplus-spoa.patch"))))

(define rust-stable-deref-trait-1.2.0
  (crate-source "stable_deref_trait" "1.2.0"
                "1lxjr8q2n534b2lhkxd6l6wcddzjvnksi58zv11f9y0jjmr15wd8"))

(define rust-static-assertions-1.1.0
  (crate-source "static_assertions" "1.1.0"
                "0gsl6xmw10gvn3zs1rv99laj5ig7ylffnh71f9l34js4nr4r7sx2"))

(define rust-strsim-0.11.1
  (crate-source "strsim" "0.11.1"
                "0kzvqlw8hxqb7y598w1s0hxlnmi84sg5vsipp3yg5na5d1rvba3x"))

(define rust-strum-macros-0.26.4
  (crate-source "strum_macros" "0.26.4"
                "1gl1wmq24b8md527cpyd5bw9rkbqldd7k1h38kf5ajd2ln2ywssc"))

(define rust-succinct-0.5.2
  (crate-source "succinct" "0.5.2"
                "0654c9gq50x7djyf25zbzz3d2pc4x3z21wmjj3qbr6d9h4hbd63p"))

(define rust-syn-1.0.109
  (crate-source "syn" "1.0.109"
                "0ds2if4600bd59wsv7jjgfkayfzy3hnazs394kz6zdkmna8l3dkj"))

(define rust-syn-2.0.106
  (crate-source "syn" "2.0.106"
                "19mddxp1ia00hfdzimygqmr1jqdvyl86k48427bkci4d08wc9rzd"))

(define rust-synstructure-0.13.2
  (crate-source "synstructure" "0.13.2"
                "1lh9lx3r3jb18f8sbj29am5hm9jymvbwh6jb1izsnnxgvgrp12kj"))

(define rust-tap-1.0.1
  (crate-source "tap" "1.0.1"
                "0sc3gl4nldqpvyhqi3bbd0l9k7fngrcl4zs47n314nqqk4bpx4sm"))

(define rust-tempfile-3.22.0
  (crate-source "tempfile" "3.22.0"
                "0lza9r7dzm4k9fghw24yql6iz59wq8xgs46a7i29ir6xz88lvyl4"))

(define rust-termcolor-1.4.1
  (crate-source "termcolor" "1.4.1"
                "0mappjh3fj3p2nmrg4y7qv94rchwi9mzmgmfflr8p2awdj7lyy86"))

(define rust-thiserror-1.0.69
  (crate-source "thiserror" "1.0.69"
                "0lizjay08agcr5hs9yfzzj6axs53a2rgx070a1dsi3jpkcrzbamn"))

(define rust-thiserror-2.0.16
  (crate-source "thiserror" "2.0.16"
                "1h30bqyjn5s9ypm668yd9849371rzwk185klwgjg503k2hadcrrl"))

(define rust-thiserror-impl-1.0.69
  (crate-source "thiserror-impl" "1.0.69"
                "1h84fmn2nai41cxbhk6pqf46bxqq1b344v8yz089w1chzi76rvjg"))

(define rust-thiserror-impl-2.0.16
  (crate-source "thiserror-impl" "2.0.16"
                "0q3r1ipr1rhff6cgrcvc0njffw17rpcqz9hdc7p754cbqkhinpkc"))

(define rust-tinystr-0.8.1
  (crate-source "tinystr" "0.8.1"
                "12sc6h3hnn6x78iycm5v6wrs2xhxph0ydm43yyn7gdfw8l8nsksx"))

(define rust-typenum-1.18.0
  (crate-source "typenum" "1.18.0"
                "0gwgz8n91pv40gabrr1lzji0b0hsmg0817njpy397bq7rvizzk0x"))

(define rust-unicode-ident-1.0.19
  (crate-source "unicode-ident" "1.0.19"
                "17bx1j1zf6b9j3kpyf74mraary7ava3984km0n8kh499h5a58fpn"))

(define rust-unicode-width-0.2.1
  (crate-source "unicode-width" "0.2.1"
                "0k0mlq7xy1y1kq6cgv1r2rs2knn6rln3g3af50rhi0dkgp60f6ja"))

(define rust-unit-prefix-0.5.1
  (crate-source "unit-prefix" "0.5.1"
                "05rq0asf2f1q5vrcv4bwf0c3y6q20asqkiqpr8wqyrfxyb7h4d1j"))

(define rust-unty-0.0.4
  (crate-source "unty" "0.0.4"
                "1blhyv01qiv5sb72sal3xa1l8nzck3answawxkkiw3fd2x1phjbd"))

(define rust-url-2.5.7
  (crate-source "url" "2.5.7"
                "0nzghdv0kcksyvri0npxbjzyx2ihprks5k590y77bld355m17g08"))

(define rust-utf8-iter-1.0.4
  (crate-source "utf8_iter" "1.0.4"
                "1gmna9flnj8dbyd8ba17zigrp9c4c3zclngf5lnb5yvz1ri41hdn"))

(define rust-utf8parse-0.2.2
  (crate-source "utf8parse" "0.2.2"
                "088807qwjq46azicqwbhlmzwrbkz7l4hpw43sdkdyyk524vdxaq6"))

(define rust-vcpkg-0.2.15
  (crate-source "vcpkg" "0.2.15"
                "09i4nf5y8lig6xgj3f7fyrvzd3nlaw4znrihw8psidvv5yk4xkdc"))

(define rust-version-check-0.9.5
  (crate-source "version_check" "0.9.5"
                "0nhhi4i5x89gm911azqbn7avs9mdacw2i3vcz3cnmz3mv4rqz4hb"))

(define rust-virtue-0.0.18
  (crate-source "virtue" "0.0.18"
                "1cgp79pzzs117kjlc3jnnkixbyaqri12j40mx2an41qhrymv27h5"))

(define rust-wasi-0.11.1+wasi-snapshot-preview1
  (crate-source "wasi" "0.11.1+wasi-snapshot-preview1"
                "0jx49r7nbkbhyfrfyhz0bm4817yrnxgd3jiwwwfv0zl439jyrwyc"))

(define rust-wasi-0.14.7+wasi-0.2.4
  (crate-source "wasi" "0.14.7+wasi-0.2.4"
                "133fq3mq7h65mzrsphcm7bbbx1gsz7srrbwh01624zin43g7hd48"))

(define rust-wasip2-1.0.1+wasi-0.2.4
  (crate-source "wasip2" "1.0.1+wasi-0.2.4"
                "1rsqmpspwy0zja82xx7kbkbg9fv34a4a2if3sbd76dy64a244qh5"))

(define rust-wasm-bindgen-0.2.102
  (crate-source "wasm-bindgen" "0.2.102"
                "1lvs1swrikd06271123slv6c2m3f2zw284a7yjscyjb6fz929lja"))

(define rust-wasm-bindgen-backend-0.2.102
  (crate-source "wasm-bindgen-backend" "0.2.102"
                "0cp7jgjj3c9dig1rpw79ymli9290jqf3nsmi48zw1lyw9c8684rs"))

(define rust-wasm-bindgen-macro-0.2.102
  (crate-source "wasm-bindgen-macro" "0.2.102"
                "0xkcwca048bqjlilgl3y6gw1bqa1dv687gfxhknv2yrn7v5b8yhd"))

(define rust-wasm-bindgen-macro-support-0.2.102
  (crate-source "wasm-bindgen-macro-support" "0.2.102"
                "1c2jqla3wfl000m1krhldf4gn8wvp6vdzm8fppdy469shha80laa"))

(define rust-wasm-bindgen-shared-0.2.102
  (crate-source "wasm-bindgen-shared" "0.2.102"
                "1vnw2mzgvb454gbgrq6mc5rip4ckdvbbphh059ka7vn29jmb0pi5"))

(define rust-web-time-1.1.0
  (crate-source "web-time" "1.1.0"
                "1fx05yqx83dhx628wb70fyy10yjfq1jpl20qfqhdkymi13rq0ras"))

(define rust-wide-0.7.33
  (crate-source "wide" "0.7.33"
                "00yd2sg83xvfrjjlwndyk49fjx8jlmlrz8byigndig32rf7dmr8c"))

(define rust-winapi-0.3.9
  (crate-source "winapi" "0.3.9"
                "06gl025x418lchw1wxj64ycr7gha83m44cjr5sarhynd9xkrm0sw"))

(define rust-winapi-i686-pc-windows-gnu-0.4.0
  (crate-source "winapi-i686-pc-windows-gnu" "0.4.0"
                "1dmpa6mvcvzz16zg6d5vrfy4bxgg541wxrcip7cnshi06v38ffxc"))

(define rust-winapi-util-0.1.11
  (crate-source "winapi-util" "0.1.11"
                "08hdl7mkll7pz8whg869h58c1r9y7in0w0pk8fm24qc77k0b39y2"))

(define rust-winapi-x86-64-pc-windows-gnu-0.4.0
  (crate-source "winapi-x86_64-pc-windows-gnu" "0.4.0"
                "0gqq64czqb64kskjryj8isp62m2sgvx25yyj3kpc2myh85w24bki"))

(define rust-windows-aarch64-gnullvm-0.53.0
  (crate-source "windows_aarch64_gnullvm" "0.53.0"
                "0r77pbpbcf8bq4yfwpz2hpq3vns8m0yacpvs2i5cn6fx1pwxbf46"))

(define rust-windows-aarch64-msvc-0.53.0
  (crate-source "windows_aarch64_msvc" "0.53.0"
                "0v766yqw51pzxxwp203yqy39ijgjamp54hhdbsyqq6x1c8gilrf7"))

(define rust-windows-i686-gnu-0.53.0
  (crate-source "windows_i686_gnu" "0.53.0"
                "1hvjc8nv95sx5vdd79fivn8bpm7i517dqyf4yvsqgwrmkmjngp61"))

(define rust-windows-i686-gnullvm-0.53.0
  (crate-source "windows_i686_gnullvm" "0.53.0"
                "04df1in2k91qyf1wzizvh560bvyzq20yf68k8xa66vdzxnywrrlw"))

(define rust-windows-i686-msvc-0.53.0
  (crate-source "windows_i686_msvc" "0.53.0"
                "0pcvb25fkvqnp91z25qr5x61wyya12lx8p7nsa137cbb82ayw7sq"))

(define rust-windows-link-0.1.3
  (crate-source "windows-link" "0.1.3"
                "12kr1p46dbhpijr4zbwr2spfgq8i8c5x55mvvfmyl96m01cx4sjy"))

(define rust-windows-link-0.2.0
  (crate-source "windows-link" "0.2.0"
                "0r9w2z96d5phmm185aq92z54jp9h2nqisa4wgc71idxbc436rr25"))

(define rust-windows-sys-0.60.2
  ;; TODO: Check bundled sources.
  (crate-source "windows-sys" "0.60.2"
                "1jrbc615ihqnhjhxplr2kw7rasrskv9wj3lr80hgfd42sbj01xgj"))

(define rust-windows-sys-0.61.0
  ;; TODO: Check bundled sources.
  (crate-source "windows-sys" "0.61.0"
                "1ajpwsmzfcsa1r7i0dxzvfn24dp3525rcd7aq95ydvdj8171h0g2"))

(define rust-windows-targets-0.53.3
  (crate-source "windows-targets" "0.53.3"
                "14fwwm136dhs3i1impqrrip7nvkra3bdxa4nqkblj604qhqn1znm"))

(define rust-windows-x86-64-gnu-0.53.0
  (crate-source "windows_x86_64_gnu" "0.53.0"
                "1flh84xkssn1n6m1riddipydcksp2pdl45vdf70jygx3ksnbam9f"))

(define rust-windows-x86-64-gnullvm-0.53.0
  (crate-source "windows_x86_64_gnullvm" "0.53.0"
                "0mvc8119xpbi3q2m6mrjcdzl6afx4wffacp13v76g4jrs1fh6vha"))

(define rust-windows-x86-64-msvc-0.53.0
  (crate-source "windows_x86_64_msvc" "0.53.0"
                "11h4i28hq0zlnjcaqi2xdxr7ibnpa8djfggch9rki1zzb8qi8517"))

(define rust-wit-bindgen-0.46.0
  (crate-source "wit-bindgen" "0.46.0"
                "0ngysw50gp2wrrfxbwgp6dhw1g6sckknsn3wm7l00vaf7n48aypi"))

(define rust-writeable-0.6.1
  (crate-source "writeable" "0.6.1"
                "1fx29zncvbrqzgz7li88vzdm8zvgwgwy2r9bnjqxya09pfwi0bza"))

(define rust-wyhash-0.5.0
  (crate-source "wyhash" "0.5.0"
                "15f26hvx6nyp4d6iswha7rm3psidxa2k2iab1f1aqgsyq9iy3xms"))

(define rust-wyz-0.5.1
  (crate-source "wyz" "0.5.1"
                "1vdrfy7i2bznnzjdl9vvrzljvs4s3qm8bnlgqwln6a941gy61wq5"))

(define rust-yoke-0.8.0
  (crate-source "yoke" "0.8.0"
                "1k4mfr48vgi7wh066y11b7v1ilakghlnlhw9snzz8vi2p00vnhaz"))

(define rust-yoke-derive-0.8.0
  (crate-source "yoke-derive" "0.8.0"
                "1dha5jrjz9jaq8kmxq1aag86b98zbnm9lyjrihy5sv716sbkrniq"))

(define rust-zerofrom-0.1.6
  (crate-source "zerofrom" "0.1.6"
                "19dyky67zkjichsb7ykhv0aqws3q0jfvzww76l66c19y6gh45k2h"))

(define rust-zerofrom-derive-0.1.6
  (crate-source "zerofrom-derive" "0.1.6"
                "00l5niw7c1b0lf1vhvajpjmcnbdp2vn96jg4nmkhq2db0rp5s7np"))

(define rust-zerotrie-0.2.2
  (crate-source "zerotrie" "0.2.2"
                "15gmka7vw5k0d24s0vxgymr2j6zn2iwl12wpmpnpjgsqg3abpw1n"))

(define rust-zerovec-0.11.4
  (crate-source "zerovec" "0.11.4"
                "0fz7j1ns8d86m2fqg2a4bzi5gnh5892bxv4kcr9apwc6a3ajpap7"))

(define rust-zerovec-derive-0.11.1
  (crate-source "zerovec-derive" "0.11.1"
                "13zms8hj7vzpfswypwggyfr4ckmyc7v3di49pmj8r1qcz9z275jv"))

(define rust-zlib-rs-0.5.2
  (crate-source "zlib-rs" "0.5.2"
                "1wh0brb3cci6ifdwwz6xasznlrgb8pr99l1z8i15qpigyj9aw1ig"))

(define rust-zstd-0.13.3
  (crate-source "zstd" "0.13.3"
                "12n0h4w9l526li7jl972rxpyf012jw3nwmji2qbjghv9ll8y67p9"))

(define rust-zstd-safe-7.2.4
  (crate-source "zstd-safe" "7.2.4"
                "179vxmkzhpz6cq6mfzvgwc99bpgllkr6lwxq7ylh5dmby3aw8jcg"))

(define rust-zstd-sys-2.0.16+zstd.1.5.7
  ;; TODO: Check bundled sources.
  (crate-source "zstd-sys" "2.0.16+zstd.1.5.7"
                "0j1pd2iaqpvaxlgqmmijj68wma7xwdv9grrr63j873yw5ay9xqci"))


(define-cargo-inputs lookup-cargo-inputs
                     (impg =>
                           (list rust-adler2-2.0.1
                                 ;; rust-agc-rs-0.1.0.aba33a3
                                 rust-aho-corasick-1.1.3
                                 rust-anstream-0.6.20
                                 rust-anstyle-1.0.11
                                 rust-anstyle-parse-0.2.7
                                 rust-anstyle-query-1.1.4
                                 rust-anstyle-wincon-3.0.10
                                 rust-anyhow-1.0.99
                                 rust-approx-0.5.1
                                 rust-arrayvec-0.5.2
                                 rust-autocfg-1.5.0
                                 rust-bgzip-0.3.1
                                 rust-bincode-2.0.1
                                 rust-bincode-derive-2.0.1
                                 rust-bio-types-1.0.4
                                 rust-bitflags-1.3.2
                                 rust-bitflags-2.9.4
                                 rust-bitvec-1.0.1
                                 rust-boomphf-0.5.9
                                 rust-bstr-0.2.17
                                 rust-bumpalo-3.19.0
                                 rust-bytemuck-1.23.2
                                 rust-bytemuck-derive-1.10.1
                                 rust-byteorder-1.5.0
                                 rust-bytes-1.10.1
                                 rust-bzip2-0.5.2
                                 rust-bzip2-sys-0.1.13+1.0.8
                                 rust-cc-1.2.37
                                 rust-cfg-if-1.0.3
                                 rust-clap-4.5.47
                                 rust-clap-builder-4.5.47
                                 rust-clap-derive-4.5.47
                                 rust-clap-lex-0.7.5
                                 rust-cmake-0.1.54
                                 rust-codespan-reporting-0.12.0
                                 rust-coitrees-0.4.0
                                 rust-colorchoice-1.0.4
                                 rust-console-0.16.1
                                 rust-core-affinity-0.8.3
                                 rust-crc32fast-1.5.0
                                 rust-crossbeam-channel-0.5.15
                                 rust-crossbeam-deque-0.8.6
                                 rust-crossbeam-epoch-0.9.18
                                 rust-crossbeam-utils-0.8.21
                                 rust-custom-derive-0.1.7
                                 rust-cxx-1.0.184
                                 rust-cxx-build-1.0.184
                                 rust-cxxbridge-cmd-1.0.184
                                 rust-cxxbridge-flags-1.0.184
                                 rust-cxxbridge-macro-1.0.184
                                 rust-derive-new-0.5.9
                                 rust-derive-new-0.6.0
                                 rust-displaydoc-0.2.5
                                 rust-either-1.15.0
                                 rust-encode-unicode-1.0.0
                                 rust-env-filter-0.1.3
                                 rust-env-logger-0.11.8
                                 rust-equivalent-1.0.2
                                 rust-errno-0.3.14
                                 rust-fastrand-2.3.0
                                 rust-find-msvc-tools-0.1.1
                                 rust-flate2-1.1.2
                                 rust-flume-0.11.1
                                 rust-fnv-1.0.7
                                 rust-foldhash-0.2.0
                                 rust-form-urlencoded-1.2.2
                                 rust-fs-utils-1.1.4
                                 rust-funty-2.0.0
                                 rust-futures-core-0.3.31
                                 rust-futures-sink-0.3.31
                                 rust-getrandom-0.2.16
                                 rust-getrandom-0.3.3
                                 rust-gfa-0.10.1
                                 rust-glob-0.3.3
                                 rust-gzp-1.0.1
                                 rust-handlegraph-0.7.0-alpha.9.3ac575e
                                 rust-hashbrown-0.15.5
                                 rust-heck-0.5.0
                                 rust-hermit-abi-0.5.2
                                 rust-hts-sys-2.2.0
                                 rust-icu-collections-2.0.0
                                 rust-icu-locale-core-2.0.0
                                 rust-icu-normalizer-2.0.0
                                 rust-icu-normalizer-data-2.0.0
                                 rust-icu-properties-2.0.1
                                 rust-icu-properties-data-2.0.1
                                 rust-icu-provider-2.0.0
                                 rust-idna-1.1.0
                                 rust-idna-adapter-1.2.1
                                 rust-ieee754-0.2.6
                                 rust-indexmap-2.11.3
                                 rust-indicatif-0.18.0
                                 rust-is-terminal-polyfill-1.70.1
                                 rust-jiff-0.2.15
                                 rust-jiff-static-0.2.15
                                 rust-jobserver-0.1.34
                                 rust-js-sys-0.3.79
                                 rust-lazy-static-1.5.0
                                 rust-lexical-core-0.7.6
                                 rust-libc-0.2.175
                                 rust-libdeflate-sys-1.24.0
                                 rust-libdeflater-1.24.0
                                 rust-liblzma-0.3.6
                                 rust-liblzma-sys-0.3.13
                                 rust-libz-ng-sys-1.1.22
                                 rust-libz-rs-sys-0.5.2
                                 rust-libz-sys-1.1.22
                                 rust-linear-map-1.2.0
                                 rust-link-cplusplus-1.0.12
                                 rust-linux-raw-sys-0.11.0
                                 rust-litemap-0.8.0
                                 rust-lock-api-0.4.13
                                 rust-log-0.4.28
                                 rust-matrixmultiply-0.3.10
                                 rust-memchr-2.7.5
                                 rust-memmap-0.7.0
                                 rust-miniz-oxide-0.8.9
                                 rust-nalgebra-0.33.2
                                 rust-nalgebra-macros-0.2.2
                                 rust-nanorand-0.7.0
                                 rust-natord-1.0.9
                                 rust-newtype-derive-0.1.6
                                 rust-niffler-3.0.0
                                 rust-nom-5.1.3
                                 rust-noodles-0.100.0
                                 rust-noodles-bgzf-0.42.0
                                 rust-num-bigint-0.4.6
                                 rust-num-complex-0.4.6
                                 rust-num-integer-0.1.46
                                 rust-num-rational-0.4.2
                                 rust-num-traits-0.2.19
                                 rust-num-cpus-1.17.0
                                 rust-once-cell-1.21.3
                                 rust-once-cell-polyfill-1.70.1
                                 rust-paste-1.0.15
                                 rust-percent-encoding-2.3.2
                                 rust-pkg-config-0.3.32
                                 rust-portable-atomic-1.11.1
                                 rust-portable-atomic-util-0.2.4
                                 rust-potential-utf-0.1.3
                                 rust-proc-macro2-1.0.101
                                 rust-quick-error-1.2.3
                                 rust-quote-1.0.40
                                 rust-r-efi-5.3.0
                                 rust-radium-0.7.0
                                 rust-rand-core-0.6.4
                                 rust-rawpointer-0.2.1
                                 rust-rayon-1.11.0
                                 rust-rayon-core-1.13.0
                                 rust-regex-1.11.2
                                 rust-regex-automata-0.1.10
                                 rust-regex-automata-0.4.10
                                 rust-regex-syntax-0.8.6
                                 rust-rust-htslib-0.46.0
                                 rust-rustc-hash-2.1.1
                                 rust-rustc-version-0.1.7
                                 rust-rustix-1.1.2
                                 rust-rustversion-1.0.22
                                 rust-ryu-1.0.20
                                 rust-safe-arch-0.7.4
                                 rust-scopeguard-1.2.0
                                 rust-scratch-1.0.9
                                 rust-semver-0.1.20
                                 rust-serde-1.0.225
                                 rust-serde-core-1.0.225
                                 rust-serde-derive-1.0.225
                                 rust-shlex-1.3.0
                                 rust-simba-0.9.1
                                 rust-smallvec-1.15.1
                                 rust-spin-0.9.8
                                 rust-spoa-rs-0.1.0.6f4f102
                                 rust-stable-deref-trait-1.2.0
                                 rust-static-assertions-1.1.0
                                 rust-strsim-0.11.1
                                 rust-strum-macros-0.26.4
                                 rust-succinct-0.5.2
                                 rust-syn-1.0.109
                                 rust-syn-2.0.106
                                 rust-synstructure-0.13.2
                                 rust-tap-1.0.1
                                 rust-tempfile-3.22.0
                                 rust-termcolor-1.4.1
                                 rust-thiserror-1.0.69
                                 rust-thiserror-2.0.16
                                 rust-thiserror-impl-1.0.69
                                 rust-thiserror-impl-2.0.16
                                 rust-tinystr-0.8.1
                                 rust-typenum-1.18.0
                                 rust-unicode-ident-1.0.19
                                 rust-unicode-width-0.2.1
                                 rust-unit-prefix-0.5.1
                                 rust-unty-0.0.4
                                 rust-url-2.5.7
                                 rust-utf8-iter-1.0.4
                                 rust-utf8parse-0.2.2
                                 rust-vcpkg-0.2.15
                                 rust-version-check-0.9.5
                                 rust-virtue-0.0.18
                                 rust-wasi-0.11.1+wasi-snapshot-preview1
                                 rust-wasi-0.14.7+wasi-0.2.4
                                 rust-wasip2-1.0.1+wasi-0.2.4
                                 rust-wasm-bindgen-0.2.102
                                 rust-wasm-bindgen-backend-0.2.102
                                 rust-wasm-bindgen-macro-0.2.102
                                 rust-wasm-bindgen-macro-support-0.2.102
                                 rust-wasm-bindgen-shared-0.2.102
                                 rust-web-time-1.1.0
                                 rust-wide-0.7.33
                                 rust-winapi-0.3.9
                                 rust-winapi-i686-pc-windows-gnu-0.4.0
                                 rust-winapi-util-0.1.11
                                 rust-winapi-x86-64-pc-windows-gnu-0.4.0
                                 rust-windows-link-0.1.3
                                 rust-windows-link-0.2.0
                                 rust-windows-sys-0.60.2
                                 rust-windows-sys-0.61.0
                                 rust-windows-targets-0.53.3
                                 rust-windows-aarch64-gnullvm-0.53.0
                                 rust-windows-aarch64-msvc-0.53.0
                                 rust-windows-i686-gnu-0.53.0
                                 rust-windows-i686-gnullvm-0.53.0
                                 rust-windows-i686-msvc-0.53.0
                                 rust-windows-x86-64-gnu-0.53.0
                                 rust-windows-x86-64-gnullvm-0.53.0
                                 rust-windows-x86-64-msvc-0.53.0
                                 rust-wit-bindgen-0.46.0
                                 rust-writeable-0.6.1
                                 rust-wyhash-0.5.0
                                 rust-wyz-0.5.1
                                 rust-yoke-0.8.0
                                 rust-yoke-derive-0.8.0
                                 rust-zerofrom-0.1.6
                                 rust-zerofrom-derive-0.1.6
                                 rust-zerotrie-0.2.2
                                 rust-zerovec-0.11.4
                                 rust-zerovec-derive-0.11.1
                                 rust-zlib-rs-0.5.2
                                 rust-zstd-0.13.3
                                 rust-zstd-safe-7.2.4
                                 rust-zstd-sys-2.0.16+zstd.1.5.7)))
