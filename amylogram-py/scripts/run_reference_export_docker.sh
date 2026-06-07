#!/usr/bin/env bash
set -Eeuo pipefail

image="${AMYLOGRAM_R_IMAGE:-amylogram-py-reference:4.3.3}"
out="${1:-tests/fixtures/amylogram_reference.json}"

mkdir -p "$(dirname "$out")"

if ! docker image inspect "$image" >/dev/null 2>&1; then
  docker build \
    -f docker/amylogram-reference.Dockerfile \
    -t "$image" \
    .
fi

docker run --rm \
  -v "$PWD:/work" \
  -w /work \
  "$image" \
  Rscript scripts/export_reference_fixtures.R "$out"
