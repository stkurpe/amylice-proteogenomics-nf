#!/bin/bash
set -Eeuo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
cd "$PROJECT_DIR"
python3 tests/test_pipeline_static.py
python3 tests/test_pipeline_logging.py
PYTHONPATH="$PROJECT_DIR" python3 -m unittest discover -s tests -p 'test_proteome_pipeline.py' -v
