#!/bin/bash
set -Eeuo pipefail

cd /home/codex
python3 tests_amyloid/test_amyloid_predictor_contracts.py
