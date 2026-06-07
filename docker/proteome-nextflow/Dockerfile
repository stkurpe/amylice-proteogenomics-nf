FROM python:3.11-slim

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
        bcftools \
        gffread \
        procps \
        samtools \
    && rm -rf /var/lib/apt/lists/*
