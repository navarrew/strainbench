#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# InterProScan annotation of strainbench cluster representatives.
# Domain-level annotations from many databases (Pfam, TIGRFAM, PANTHER, …)
# plus InterPro2GO mappings.
#
# WARNING: InterProScan is heavyweight — ~50 GB of bundled data, requires
# a JDK 11+, and a single ~3,600-cluster run can take 1–2 hours.
# Worth it for the domain-level resolution it provides.
#
# Prerequisites (once per machine):
#   Download InterProScan from:
#     https://interproscan-docs.readthedocs.io/en/latest/HowToDownload.html
#   Extract somewhere stable, e.g. /opt/interproscan-5.x.x
#   Run `python3 setup.py interproscan.properties` once to initialize.
#   Make sure 'interproscan.sh' is on PATH (or set INTERPROSCAN_BIN).
#
# Running:
#   bash interproscan.sh
#
# Customize via env vars:
#   INTERPROSCAN_BIN     path to interproscan.sh   (default: 'interproscan.sh' on PATH)
#   INTERPROSCAN_THREADS CPU threads               (default: 8)
#   INTERPROSCAN_APPS    comma-separated analyses  (default: 'Pfam,TIGRFAM,SUPERFAMILY,SMART,PANTHER')
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

INTERPROSCAN_BIN="${INTERPROSCAN_BIN:-interproscan.sh}"
INTERPROSCAN_THREADS="${INTERPROSCAN_THREADS:-8}"
INTERPROSCAN_APPS="${INTERPROSCAN_APPS:-Pfam,TIGRFAM,SUPERFAMILY,SMART,PANTHER}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="interproscan_out"
OUTPUT_FILE="$OUTPUT_DIR/interproscan.tsv"

if ! command -v "$INTERPROSCAN_BIN" >/dev/null 2>&1; then
    echo "ERROR: '$INTERPROSCAN_BIN' not found on PATH." >&2
    echo "Set INTERPROSCAN_BIN to the full path of interproscan.sh, or add" >&2
    echo "the InterProScan install dir to your PATH." >&2
    exit 2
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running InterProScan (this may take 1-2 hours)..."
echo "  input:        $REPS_FASTA"
echo "  analyses:     $INTERPROSCAN_APPS"
echo "  threads:      $INTERPROSCAN_THREADS"
echo "  output:       $OUTPUT_FILE"
echo

# -goterms enables InterPro2GO; -pa enables pathway annotations; -dp disables
# pre-calculated lookup (we want signature hits, not just the InterPro mapping).
"$INTERPROSCAN_BIN" \
    --input "$REPS_FASTA" \
    --output-file-base "$OUTPUT_DIR/interproscan" \
    --formats TSV \
    --applications "$INTERPROSCAN_APPS" \
    --cpu "$INTERPROSCAN_THREADS" \
    --goterms \
    --pathways

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
