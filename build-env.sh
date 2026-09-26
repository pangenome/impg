#!/usr/bin/env bash
# Reconstructed build environment (guix daemon offline): store-path profile
# matching repo manifest.scm (gcc-toolchain, htslib, libdeflate, zlib, zstd,
# bzip2, xz, openssl). Used only for local cargo builds in this session.
export HTS=/gnu/store/4ld0dqv4dv4xh4g5aw21pxdf3hd547k1-htslib-1.16
export DEFLATE=/gnu/store/2m6s21pnp7c24x89d77kmqi52ysws3qs-libdeflate-1.19
export ZLIB=/gnu/store/2jzflb96n50lknzn7znfq9ri4qlbifjc-zlib-1.3.1
export ZSTD=/gnu/store/1v5akbxgf2kg4yyx4vbn5i2m1z9pav4l-zstd-1.5.5
export BZ2=/gnu/store/0ngy7fmrp9g3q35wwccnjxzaaqq9f8wi-bzip2-1.0.8
export XZ=/gnu/store/0i4anagv296gx14v1rxkli1n2n17p5p9-xz-5.4.5
export SSL=/gnu/store/2kv5mhiv0w362p3dn06lyqlk51xqffgl-openssl-3.5.5
export JEM=/gnu/store/1vn63xy1cq4kc74w3l93b70wb385jbbz-jemalloc-5.3.0
export CPATH="$HTS/include:$DEFLATE/include:$ZLIB/include:$ZSTD/include:$BZ2/include:$XZ/include:$SSL/include:$JEM/include${CPATH:+:$CPATH}"
export LIBRARY_PATH="$HTS/lib:$DEFLATE/lib:$ZLIB/lib:$ZSTD/lib:$BZ2/lib:$XZ/lib:$SSL/lib:$JEM/lib${LIBRARY_PATH:+:$LIBRARY_PATH}"
export PKG_CONFIG_PATH="$HTS/lib/pkgconfig:$DEFLATE/lib/pkgconfig:$ZLIB/lib/pkgconfig:$ZSTD/lib/pkgconfig:$BZ2/lib/pkgconfig:$XZ/lib/pkgconfig:$SSL/lib/pkgconfig${PKG_CONFIG_PATH:+:$PKG_CONFIG_PATH}"
export CARGO_TARGET_DIR="${CARGO_TARGET_DIR:-target}"
