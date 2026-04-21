#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# DefenseFinder annotation.
#
# ⚠️ ARCHITECTURAL NOTE: DefenseFinder is designed to detect MULTI-PROTEIN
# defense systems (CRISPR-Cas, R-M, Abi, BREX, …) by looking at gene synteny.
# Running it on cluster representatives (out-of-context proteins) only catches
# single-protein systems and misses everything else.
#
# Better workflow: run DefenseFinder on each per-strain genome's protein FASTA
# (with --gff for synteny), then map the detected systems back to clusters
# via the protein IDs. That workflow isn't built into strainbench yet.
#
# This script runs the cluster-rep version anyway, which is useful for
# spot-checks but should not be your primary defense-system annotation.
#
# Prerequisites (once per machine):
#   conda create -n defensefinder -c bioconda mdmparis-defense-finder
#   conda activate defensefinder
#   defense-finder update   # downloads MacSyFinder models (once)
#
# Running:
#   conda activate defensefinder && bash defensefinder.sh
#
# Customize via env vars:
#   DEFENSEFINDER_THREADS   workers (default: 4)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

DEFENSEFINDER_THREADS="${DEFENSEFINDER_THREADS:-4}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="defensefinder_out"

if ! command -v defense-finder >/dev/null 2>&1; then
    echo "ERROR: defense-finder not found on PATH. 'conda activate defensefinder'?" >&2
    exit 2
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running DefenseFinder..."
echo "  input:    $REPS_FASTA"
echo "  workers:  $DEFENSEFINDER_THREADS"
echo "  output:   $OUTPUT_DIR/defense_finder_genes.tsv"
echo

defense-finder run \
    --workers "$DEFENSEFINDER_THREADS" \
    --out-dir "$OUTPUT_DIR" \
    "$REPS_FASTA"

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
