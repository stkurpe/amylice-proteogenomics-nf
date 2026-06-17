FROM python:3.8-slim

ENV PIP_CACHE_DIR=/tmp/pip-cache

RUN python -m pip install --no-cache-dir \
    numpy \
    scipy \
    scikit-learn \
    pandas \
    xgboost \
    tabulate

WORKDIR /data
