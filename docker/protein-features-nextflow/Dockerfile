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
    && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("Peptides", "protr", "zoo"), repos="http://cran.r-project.org"); stopifnot(requireNamespace("Peptides", quietly=TRUE)); stopifnot(requireNamespace("protr", quietly=TRUE)); stopifnot(requireNamespace("zoo", quietly=TRUE))'

WORKDIR /work
