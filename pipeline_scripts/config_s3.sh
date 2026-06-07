#!/bin/bash

CONFIG_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${CONFIG_DIR}/config.sh"

export S3_BUCKET="${S3_BUCKET:-s3://bioinfo-data-amylice-2026}"
