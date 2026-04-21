#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# NCBI AMRFinderPlus annotation of strainbench cluster representatives.
# Identifies AMR / stress / virulence genes by HMM + BLAST against curated DB.
#
# Prerequisites (once per machine):
#   conda create -n amrfinder -c bioconda ncbi-amrfinderplus
#   conda activate amrfinder
#   amrfinder -u    # download/update the database (once, then periodically)
#
# Running:
#   conda activate amrfinder && bash amrfinder.sh
#
# Customize via env vars:
#   AMRFINDER_THREADS   number of CPU threads      (default: 8)
#   AMRFINDER_OPTS      extra flags passed through (default: empty)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

AMRFINDER_THREADS="${AMRFINDER_THREADS:-8}"
AMRFINDER_OPTS="${AMRFINDER_OPTS:-}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="amrfinder_out"
OUTPUT_FILE="$OUTPUT_DIR/amrfinder.tsv"

if ! command -v amrfinder >/dev/null 2>&1; then
    echo "ERROR: amrfinder not found on PATH. Did you 'conda activate amrfinder'?" >&2
    exit 2
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running AMRFinderPlus..."
echo "  input:    $REPS_FASTA"
echo "  threads:  $AMRFINDER_THREADS"
echo "  output:   $OUTPUT_FILE"
echo

# -p = protein input; --plus enables stress/virulence detection in addition to AMR.
amrfinder \
    -p "$REPS_FASTA" \
    -o "$OUTPUT_FILE" \
    --plus \
    --threads "$AMRFINDER_THREADS" \
    $AMRFINDER_OPTS

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
