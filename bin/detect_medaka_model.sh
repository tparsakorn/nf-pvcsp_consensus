#!/usr/bin/env bash
# Usage: detect_medaka_model.sh <reads.fastq> [override_model]
# Prints the medaka model name derived from the Dorado RG tag in the first read header,
# e.g. ...dna_r10.4.1_e8.2_400bps_hac@v5.2.0  ->  r1041_e82_400bps_hac_v5.2.0
set -euo pipefail
fastq="$1"
override="${2:-}"
if [ -n "$override" ]; then printf '%s' "$override"; exit 0; fi

rg=$(head -1 "$fastq" | grep -o 'RG:Z:[^[:space:]]*' || true)
model=$(printf '%s' "$rg" | sed -nE 's/.*dna_(r[0-9]+)\.([0-9]+)\.([0-9]+)_(e[0-9]+)\.([0-9]+)_([0-9]+bps)_([a-z]+)@(v[0-9.]+).*/\1\2\3_\4\5_\6_\7_\8/p')
if [ -z "$model" ]; then
    echo "ERROR: cannot derive medaka model from first header of $fastq (RG tag: '${rg}')." >&2
    echo "       Re-run with --medaka_model <model> (see: medaka tools list_models)." >&2
    exit 1
fi
printf '%s' "$model"
