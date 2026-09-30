#!/usr/bin/env bash
# Run nf-pvcsp on one or more barcodes.
#
# Usage (run from the folder that holds your FASTQ files):
#   /path/to/nf-pvcsp/run_pvcsp.sh [-d FASTQ_DIR] [-o OUTDIR] [-c MAX_CPUS] [-n] barcode07 barcode08
#
#   -d  folder containing <barcode>.fastq   (default: current folder)
#   -o  output folder                       (default: ./results)
#   -c  maximum CPUs to use                 (default: all cores)
#   -n  dry run: build the samplesheet and check the workflow, run nothing
#
# Barcode names are given WITHOUT ".fastq".
# Nextflow's work/ and .nextflow/ folders are created in the folder you run this from.
# Set NF_ENV to the name of your conda environment that contains Nextflow (default: nextflow);
# it is only activated if `nextflow` is not already on your PATH.
set -euo pipefail

PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
FASTQ_DIR="$PWD"
OUTDIR="$PWD/results"
MAX_CPUS=""
DRY=""

while getopts "d:o:c:nh" opt; do
    case "$opt" in
        d) FASTQ_DIR="$(cd "$OPTARG" && pwd)" ;;
        o) mkdir -p "$OPTARG"; OUTDIR="$(cd "$OPTARG" && pwd)" ;;
        c) MAX_CPUS="$OPTARG" ;;
        n) DRY="-preview" ;;
        h|*) sed -n '2,15p' "${BASH_SOURCE[0]}"; exit 0 ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -eq 0 ]; then
    echo "ERROR: give at least one barcode name, e.g.  run_pvcsp.sh barcode07 barcode08" >&2
    exit 1
fi

# 1) samplesheet (one file per run, kept next to the results)
mkdir -p "$OUTDIR/samplesheets"
SHEET="$OUTDIR/samplesheets/samplesheet_$(date +%Y%m%d_%H%M%S).csv"
echo "barcode_id,fastq" > "$SHEET"
for bc in "$@"; do
    fq="$FASTQ_DIR/${bc}.fastq"
    if [ ! -f "$fq" ]; then
        echo "ERROR: not found: $fq" >&2
        exit 1
    fi
    echo "${bc},${fq}" >> "$SHEET"
done
echo "Samplesheet: $SHEET"
cat "$SHEET"

# 2) make sure Nextflow is available
if ! command -v nextflow >/dev/null 2>&1; then
    CONDA_BASE="$(conda info --base 2>/dev/null || true)"
    [ -n "$CONDA_BASE" ] || CONDA_BASE="$HOME/miniconda3"
    CONDA_SH="$CONDA_BASE/etc/profile.d/conda.sh"
    if [ ! -f "$CONDA_SH" ]; then
        echo "ERROR: 'nextflow' not found and conda.sh not found at $CONDA_SH" >&2
        echo "       Install Nextflow, or set NF_ENV / activate your Nextflow environment first." >&2
        exit 1
    fi
    # conda's activate scripts are not "set -u" safe
    set +u
    # shellcheck disable=SC1090
    source "$CONDA_SH"
    conda activate "${NF_ENV:-nextflow}"
    set -u
fi

# 3) run (-resume reuses finished steps if you re-run after a crash)
EXTRA=()
[ -n "$MAX_CPUS" ] && EXTRA+=(--max_cpus "$MAX_CPUS")
nextflow run "$PIPE_DIR" --samplesheet "$SHEET" --outdir "$OUTDIR" -resume ${EXTRA[@]+"${EXTRA[@]}"} $DRY

echo
echo "Done. Results: $OUTDIR/<barcode>/"
