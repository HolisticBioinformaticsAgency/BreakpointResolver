#!/usr/bin/env bash
# submit_mh_array.sh
# Submits one Slurm array job per barcode × platform combination.
# Usage: bash submit_mh_array.sh

set -euo pipefail

SCRIPT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/automated_precise_breakpoint_finding/batch_scripts/ont_pacbio_batch_scripts"
TASK_LIST="$SCRIPT_DIR/mh_tasks.txt"

cd "$SCRIPT_DIR"

# Barcodes 1001–1126 with both ONT and PacBio data
BOTH_PLATFORMS=(
  1001
  1002
  1003
  1004
  1005
  1006
  1007
  1008
  1097
  1098
  1099
  1100
  1101
  1102
  1103
  1104
)

# These samples have ONT data only
ONT_ONLY=(
  1113
  1114
  1115
  1116
  1117
  1118
  1119
  1125
  1126
)

# Regenerate task list every time you submit, in case barcodes change.
: > "$TASK_LIST"

for barcode in "${BOTH_PLATFORMS[@]}"; do
  printf "%s ONT\n" "$barcode" >> "$TASK_LIST"
  printf "%s PacBio\n" "$barcode" >> "$TASK_LIST"
done

for barcode in "${ONT_ONLY[@]}"; do
  printf "%s ONT\n" "$barcode" >> "$TASK_LIST"
done

N=$(wc -l < "$TASK_LIST")

echo "Generated task list: $TASK_LIST"
echo "Submitting $N tasks as a Slurm array..."
echo
cat "$TASK_LIST"
echo

sbatch \
  --job-name=mh_array \
  --array=1-"${N}" \
  --output="$SCRIPT_DIR/logs/mh_%A_%a.out" \
  --error="$SCRIPT_DIR/logs/mh_%A_%a.err" \
  --time=4:00:00 \
  --mem=16G \
  --cpus-per-task=2 \
  "$SCRIPT_DIR/run_mh_array_worker.sh"

echo "Submitted. Monitor with: squeue -u $USER"