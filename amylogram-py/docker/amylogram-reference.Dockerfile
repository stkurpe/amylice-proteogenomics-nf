FROM rocker/r-ver:4.3.3

RUN apt-get update \
  && apt-get install -y --no-install-recommends \
    libgmp3-dev \
    libuv1-dev \
    pkg-config \
    zlib1g-dev \
  && rm -rf /var/lib/apt/lists/*

RUN Rscript -e 'install.packages(c("AmyloGram", "biogram", "seqinr", "jsonlite"), repos="http://cran.r-project.org")' \
  && Rscript -e 'library(AmyloGram); library(biogram); library(seqinr); library(jsonlite)'

WORKDIR /work
