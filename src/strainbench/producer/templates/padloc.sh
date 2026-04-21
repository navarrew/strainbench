#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# PADLOC annotation.
#
# ⚠️ ARCHITECTURAL NOTE: like DefenseFinder, PADLOC detects MULTI-PROTEIN
# defense systems via gene synteny. It REQUIRES a GFF alongside the FASTA so
# it can read gene positions on a contig. This means it doesn't run on
# strainbench's cluster representatives (which have no genomic context).
#
# To use PADLOC properly, run it on per-strain genome protein+GFF pairs
# instead, then map detected systems back to clusters. That per-genome
# workflow isn't built into strainbench yet — this template is here for
# completeness.
#
# Prerequisites (once per machine):
#   conda create -n padloc -c bioconda padloc
#   conda activate padloc
#   padloc --db-update   # downloads HMM library (once, then periodically)
#
# Running:
#   conda activate padloc && bash padloc.sh
#
# Customize via env vars:
#   PADLOC_CPUS   number of CPU threads (default: 4)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

PADLOC_CPUS="${PADLOC_CPUS:-4}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="padloc_out"

if ! command -v padloc >/dev/null 2>&1; then
    echo "ERROR: padloc not found on PATH. Did you 'conda activate padloc'?" >&2
    exit 2
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running PADLOC..."
echo "  input:    $REPS_FASTA"
echo "  threads:  $PADLOC_CPUS"
echo "  output:   $OUTPUT_DIR/<basename>_padloc.csv"
echo

# PADLOC writes <input_basename>_padloc.csv into --outdir.
padloc \
    --faa "$REPS_FASTA" \
    --outdir "$OUTPUT_DIR" \
    --cpu "$PADLOC_CPUS"

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
