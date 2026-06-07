PYTHON ?= python3
VENV ?= .venv
IMAGE ?= amylogram-py:local
FIXTURE_JSON ?= tests/fixtures/amylogram_reference.json
LOOKUP_BIN ?= tests/fixtures/amylogram_sixmer_probabilities.bin
BENCH_COUNT ?= 2000
BENCH_LENGTH ?= 300

.PHONY: help install-dev test parity benchmark demo docker-build docker-smoke

help:
	@printf "%s\n" "AmyloGram-Py commands:"
	@printf "%s\n" "  make install-dev   Install package in editable mode"
	@printf "%s\n" "  make test          Run unit tests"
	@printf "%s\n" "  make parity        Regenerate parity evidence report"
	@printf "%s\n" "  make benchmark     Regenerate lookup benchmark report"
	@printf "%s\n" "  make demo          Run local demo prediction"
	@printf "%s\n" "  make docker-build  Build production Docker image"
	@printf "%s\n" "  make docker-smoke  Run Docker smoke test with mounted inputs"

install-dev:
	$(PYTHON) -m venv $(VENV)
	$(VENV)/bin/python -m pip install setuptools wheel
	$(VENV)/bin/python -m pip install --no-build-isolation -e .

test:
	PYTHONPATH=src $(PYTHON) -m unittest discover -s tests -v

parity:
	$(PYTHON) evidence/run_parity_checks.py --report evidence/amylogram_parity_report.md

benchmark:
	$(PYTHON) evidence/benchmark_lookup_prediction.py \
		--count $(BENCH_COUNT) \
		--length $(BENCH_LENGTH) \
		--report evidence/lookup_benchmark.md

demo:
	PYTHONPATH=src $(PYTHON) -m amylogram_py.cli \
		evidence/demo_input.fasta \
		evidence/demo_predictions.csv \
		--sixmer-table $(LOOKUP_BIN) \
		--report-md evidence/demo_prediction_report.md \
		--report-json evidence/demo_prediction_report.json \
		--skipped-tsv evidence/demo_skipped.tsv \
		--top-k 2 \
		--top-tsv evidence/demo_top_hits.tsv

docker-build:
	docker build -f docker/Dockerfile -t $(IMAGE) .

docker-smoke:
	docker run --rm \
		-v "$(PWD)/evidence:/data/evidence:rw" \
		-v "$(PWD)/tests/fixtures:/data/tests/fixtures:ro" \
		$(IMAGE) \
		evidence/demo_input.fasta \
		evidence/docker_smoke_predictions.csv \
		--sixmer-table tests/fixtures/amylogram_sixmer_probabilities.bin \
		--report-md evidence/docker_smoke_prediction_report.md \
		--report-json evidence/docker_smoke_prediction_report.json \
		--skipped-tsv evidence/docker_smoke_skipped.tsv \
		--top-k 2 \
		--top-tsv evidence/docker_smoke_top_hits.tsv
