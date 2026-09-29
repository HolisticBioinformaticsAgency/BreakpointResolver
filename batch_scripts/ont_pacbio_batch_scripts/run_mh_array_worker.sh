#!/usr/bin/env bash
# run_mh_array_worker.sh
# Slurm runs this script for each array task.
# Reads its own line from mh_tasks.txt using $SLURM_ARRAY_TASK_ID.
# All analysis parameters are read from params_batch.yaml.

SCRIPT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/automated_precise_breakpoint_finding/batch_scripts/ont_pacbio_batch_scripts"
TASK_LIST="$SCRIPT_DIR/mh_tasks.txt"
PARAMS="$SCRIPT_DIR/params_batch.yaml"

# Read barcode and platform for this task index
TASK_LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$TASK_LIST")
BARCODE=$(echo "$TASK_LINE" | awk '{print $1}')
PLATFORM=$(echo "$TASK_LINE" | awk '{print $2}')

echo "=== Task $SLURM_ARRAY_TASK_ID: $BARCODE $PLATFORM ==="
echo "Node: $SLURMD_NODENAME"
echo "Params: $PARAMS"
echo "Started: $(date)"

# ── Activate conda ────────────────────────────────────────────────────────────
source "/usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh"
conda activate ont
echo "Python: $(which python3)"

# ── Run worker ────────────────────────────────────────────────────────────────
mkdir -p "$SCRIPT_DIR/logs"
python3 "$SCRIPT_DIR/run_mh_worker.py" "$BARCODE" "$PLATFORM" "$PARAMS"

echo "Finished: $(date)"
