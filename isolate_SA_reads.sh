#!/usr/bin/env bash
# Usage: ./extract_sa_reads.sh <input.bam> <region> [output.bam]
#   <region> format: chr:start-end  (e.g. chr1:1000000-2000000)
#   [output.bam] is optional; defaults to <input>_sa_<region>.bam

module load samtools 2>/dev/null || true

set -euo pipefail

# ── Argument parsing ──────────────────────────────────────────────────────────
if [[ $# -lt 2 ]]; then
    echo "Usage: $0 <input.bam> <region> [output.bam]"
    echo "  region format: chr:start-end  (e.g. chr1:1000000-2000000)"
    exit 1
fi

INPUT_BAM="$1"
REGION="$2"
SAFE_REGION="${REGION//[^a-zA-Z0-9_]/_}"
DEFAULT_OUT="${INPUT_BAM%.bam}_sa_${SAFE_REGION}.bam"
OUTPUT_BAM="${3:-$DEFAULT_OUT}"

# ── Dependency check ──────────────────────────────────────────────────────────
if ! command -v samtools &>/dev/null; then
    echo "Error: samtools is not installed or not in PATH." >&2
    exit 1
fi

# ── Validate input ────────────────────────────────────────────────────────────
if [[ ! -f "$INPUT_BAM" ]]; then
    echo "Error: Input BAM file not found: $INPUT_BAM" >&2
    exit 1
fi

# Index if missing
if [[ ! -f "${INPUT_BAM}.bai" && ! -f "${INPUT_BAM%.bam}.bai" ]]; then
    echo "Index not found. Indexing ${INPUT_BAM}..."
    samtools index "$INPUT_BAM"
fi

echo "Extracting SA / supplementary reads in region: ${REGION}"

TMP_LOCUS=$(mktemp /tmp/sa_locus_XXXXXX.bam)
TMP_QNAMES=$(mktemp /tmp/sa_qnames_XXXXXX.txt)
TMP_FINAL=$(mktemp /tmp/sa_final_XXXXXX.bam)

cleanup() { rm -f "$TMP_LOCUS" "${TMP_LOCUS}.bai" "$TMP_QNAMES" "$TMP_FINAL" "${TMP_FINAL}.bai"; }
trap cleanup EXIT

# ── Step 1: Extract locus BAM (all reads overlapping the window) ──────────────
samtools view -b -h "$INPUT_BAM" "$REGION" -o "$TMP_LOCUS"
samtools index "$TMP_LOCUS"

# ── Step 2: Filter for SA-tagged or supplementary reads within the locus ─────
# -e expression matches:
#   flag & 2048  → supplementary alignment (0x800)
#   exists([SA]) → read carries an SA auxiliary tag (split-read primary)
samtools view -b -h \
    -e 'flag & 2048 || exists([SA])' \
    "$TMP_LOCUS" \
    | samtools view -b -o "$TMP_FINAL"

samtools index "$TMP_FINAL"

SA_COUNT=$(samtools view -c "$TMP_FINAL")
echo "  → Found ${SA_COUNT} SA / supplementary records in window."

if [[ "$SA_COUNT" -eq 0 ]]; then
    echo "No SA / supplementary reads found in region ${REGION}. Exiting."
    exit 0
fi


# ── Step 4: Sort and write final output ───────────────────────────────────────
samtools sort -o "$OUTPUT_BAM" "$TMP_FINAL"
samtools index "$OUTPUT_BAM"

echo ""
echo "Done! Output written to: ${OUTPUT_BAM}"
samtools flagstat "$OUTPUT_BAM"