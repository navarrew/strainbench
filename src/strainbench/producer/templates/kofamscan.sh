#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────────────
# kofamscan annotation of strainbench cluster representatives.
#
# Prerequisites (once per machine):
#   conda create -n kofam -c bioconda kofamscan
#   conda activate kofam
#   mkdir -p $HOME/databases/kofam && cd $HOME/databases/kofam
#   wget https://www.genome.jp/ftp/db/kofam/profiles.tar.gz && tar xzf profiles.tar.gz
#   wget https://www.genome.jp/ftp/db/kofam/ko_list.gz && gunzip ko_list.gz
#
# Running:
#   conda activate kofam
#   bash kofamscan.sh
#
# Customize by exporting env vars BEFORE running:
#   KOFAM_PROFILES   HMM profiles dir (default: ~/databases/kofam/profiles)
#   KOFAM_KO_LIST    ko_list file     (default: ~/databases/kofam/ko_list)
#   KOFAM_THREADS    number of CPU threads (default: 8)
# ──────────────────────────────────────────────────────────────────────────────
set -euo pipefail

KOFAM_PROFILES="${KOFAM_PROFILES:-$HOME/databases/kofam/profiles}"
KOFAM_KO_LIST="${KOFAM_KO_LIST:-$HOME/databases/kofam/ko_list}"
KOFAM_THREADS="${KOFAM_THREADS:-8}"

REPS_FASTA="${STRAINBENCH_REPS:-reps.faa}"
OUTPUT_DIR="kofamscan_out"
OUTPUT_FILE="$OUTPUT_DIR/kofamscan_detail.tsv"

# Pre-flight checks
if ! command -v exec_annotation >/dev/null 2>&1; then
    echo "ERROR: exec_annotation not found on PATH." >&2
    echo "Did you 'conda activate kofam'?" >&2
    exit 2
fi
if [ ! -d "$KOFAM_PROFILES" ]; then
    echo "ERROR: kofam profiles dir not found: $KOFAM_PROFILES" >&2
    echo "Download the KEGG profiles once (see header of this script)." >&2
    exit 3
fi
if [ ! -s "$KOFAM_KO_LIST" ]; then
    echo "ERROR: kofam ko_list file not found: $KOFAM_KO_LIST" >&2
    exit 3
fi
if [ ! -s "$REPS_FASTA" ]; then
    echo "ERROR: $REPS_FASTA not found. Run 'strainbench export-reps' first." >&2
    exit 4
fi

mkdir -p "$OUTPUT_DIR"
echo "Running kofamscan..."
echo "  input:     $REPS_FASTA"
echo "  profiles:  $KOFAM_PROFILES"
echo "  ko_list:   $KOFAM_KO_LIST"
echo "  threads:   $KOFAM_THREADS"
echo "  output:    $OUTPUT_FILE"
echo

exec_annotation \
    -p "$KOFAM_PROFILES" \
    -k "$KOFAM_KO_LIST" \
    --cpu "$KOFAM_THREADS" \
    -f detail-tsv \
    -o "$OUTPUT_FILE" \
    "$REPS_FASTA"

echo
echo "Done. Import with:"
echo "  strainbench import-annotation --workdir ."
