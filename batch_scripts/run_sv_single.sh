#!/usr/bin/env bash
#SBATCH --job-name=run_sv_single
#SBATCH --output=run_sv_single%j.log
#SBATCH --time=16:00:00
#SBATCH --mem=16G
#SBATCH --cpus-per-task=2

SCRIPT_DIR="/home/zlaw0001/vh83_scratch/projects/temp_dnascreen_copy/dnascreen/ONT_PacBio_CNV_calling/automated_precise_breakpoint_finding/batch_scripts"

# Optional: pass a custom params file as first argument, otherwise uses params_single.yaml
PARAMS="${1:-$SCRIPT_DIR/params_single.yaml}"

cd "$SCRIPT_DIR"

# ── Activate conda ─────────────────────────────────────────────────────────────
source /usr/local/anaconda/5.1.0-Python3.6-gcc5/etc/profile.d/conda.sh
conda activate ont
echo "Python: $(which python3)"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Starting single-BAM MH run"
echo "Params: $PARAMS"

python3 "$SCRIPT_DIR/run_sv_single.py" "$PARAMS"

echo "[$(date '+%Y-%m-%d %H:%M:%S')] Done."
