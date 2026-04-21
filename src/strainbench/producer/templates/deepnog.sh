#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# DeepNOG annotation of strainbench cluster representatives.
# Assigns each rep to an orthology group via a pre-trained neural network.
#
# Prerequisites (once per machine):
#   conda create -n deepnog -c conda-forge -y python=3.11 numpy biopython pytorch
#   conda activate deepnog
#   pip install deepnog
#   # First run will auto-download the model (~150 MB) to ~/deepnog_data/
#
# Running:
#   conda activate deepnog
#   bash deepnog.sh
#
# Customize via env vars:
#   DEEPNOG_DB        database to classify against    (default: eggNOG5)
#   DEEPNOG_TAX       taxon level (2 = Bacteria)      (default: 2)
#   DEEPNOG_THREADS   number of CPU workers           (default: 4)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

DEEPNOG_DB="${DEEPNOG_DB:-eggNOG5}"
DEEPNOG_TAX="${DEEPNOG_TAX:-2}"
DEEPNOG_THREADS="${DEEPNOG_THREADS:-4}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="deepnog_out"
OUTPUT_FILE="$OUTPUT_DIR/deepnog_classify.csv"

if ! command -v deepnog >/dev/null 2>&1; then
    echo "ERROR: deepnog not found on PATH." >&2
    echo "Did you 'conda activate deepnog'?" >&2
    exit 2
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running DeepNOG..."
echo "  input:         $REPS_FASTA"
echo "  database:      $DEEPNOG_DB"
echo "  tax level:     $DEEPNOG_TAX"
echo "  workers:       $DEEPNOG_THREADS"
echo "  output:        $OUTPUT_FILE"
echo "(First run downloads the model to ~/deepnog_data/ — be patient.)"
echo

deepnog infer \
    --out "$OUTPUT_FILE" \
    --database "$DEEPNOG_DB" \
    --tax "$DEEPNOG_TAX" \
    --num-workers "$DEEPNOG_THREADS" \
    "$REPS_FASTA"

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
