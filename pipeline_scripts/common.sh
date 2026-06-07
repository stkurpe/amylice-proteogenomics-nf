#!/bin/bash

# Shared helpers for the bioinformatics pipeline.

pipeline_die() {
    echo "ERROR: $*" >&2
    exit 1
}

require_sample_id() {
    local sample_id="${1:-}"

    if [[ -z "$sample_id" ]]; then
        pipeline_die "Sample id is empty"
    fi

    if [[ ! "$sample_id" =~ ^[A-Za-z0-9][A-Za-z0-9._-]{0,127}$ ]]; then
        pipeline_die "Unsafe sample id: $sample_id"
    fi

    case "$sample_id" in
        .*|*-|*/*|*'..'* ) pipeline_die "Unsafe sample id: $sample_id" ;;
    esac

    printf '%s' "$sample_id"
}

require_file() {
    local path="${1:-}"
    local label="${2:-file}"

    if [[ ! -f "$path" ]]; then
        pipeline_die "Missing ${label}: ${path}"
    fi
}
