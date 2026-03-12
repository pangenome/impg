FROM debian:bullseye-slim AS binary

LABEL description="impg: implicit pangenome graph"
LABEL base_image="debian:bullseye-slim"
LABEL software="impg"
LABEL about.home="https://github.com/pangenome/impg"
LABEL about.license="SPDX:MIT"
LABEL maintainer="Andrea Guarracino <aguarracino@tgen.org>"

# System dependencies
RUN apt-get update \
    && apt-get install -y \
                       git \
                       bash \
                       curl \
                       ca-certificates \
                       build-essential \
                       cmake \
                       pkg-config \
                       zlib1g-dev \
                       libzstd-dev \
                       libbz2-dev \
                       liblzma-dev \
                       libcurl4-gnutls-dev \
                       libhts-dev \
                       libgsl-dev \
                       libjemalloc-dev \
                       libclang-dev \
    && rm -rf /var/lib/apt/lists/*

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"

# Portable build: no -march=native
ENV PORTABLE=1

RUN git clone --recursive https://github.com/pangenome/impg.git

RUN cd impg \
    && cargo build --release \
    && cp target/release/impg /usr/local/bin/impg \
    && cp target/release/gfaffix /usr/local/bin/gfaffix \
    && cp target/release/wfmash /usr/local/bin/wfmash 2>/dev/null || true \
    && cp target/release/libwfa2*.so* /usr/local/lib/ 2>/dev/null || true \
    && for bin in FastGA FAtoGDB GIXmake GIXrm ALNtoPAF PAFtoALN ONEview; do \
         cp target/release/$bin /usr/local/bin/$bin 2>/dev/null || true; \
       done \
    && rm -rf target .git \
    && rustup self uninstall -y \
    && apt-get clean \
    && apt-get purge \
    && rm -rf /var/lib/apt/lists/* /tmp/* /root/.cargo

RUN ldconfig \
    && chmod 777 /usr/local/bin/impg /usr/local/bin/gfaffix \
    && chmod 777 /usr/local/bin/wfmash 2>/dev/null || true \
    && for bin in FastGA FAtoGDB GIXmake GIXrm ALNtoPAF PAFtoALN ONEview; do \
         chmod 777 /usr/local/bin/$bin 2>/dev/null || true; \
       done

ENTRYPOINT ["impg"]
