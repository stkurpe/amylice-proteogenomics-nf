FROM rocker/r-ver:4.3.3

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        ca-certificates \
        curl \
        gcc \
        g++ \
        gfortran \
        make \
        pkg-config \
        zlib1g-dev \
        libbz2-dev \
        liblzma-dev \
        libcurl4-openssl-dev \
        libssl-dev \
        libxml2-dev \
        libuv1-dev \
        libgmp3-dev \
        libfontconfig1-dev \
        libfreetype6-dev \
        libfribidi-dev \
        libharfbuzz-dev \
        libjpeg-dev \
        libpng-dev \
        libtiff5-dev \
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'options(Ncpus=max(1, parallel::detectCores()-1)); install.packages(c("AmyloGram", "biogram", "seqinr", "jsonlite"), repos="http://cran.r-project.org"); stopifnot(requireNamespace("AmyloGram", quietly=TRUE)); stopifnot(requireNamespace("seqinr", quietly=TRUE))'

WORKDIR /work
