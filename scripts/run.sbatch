#!/usr/bin/bash

#SBATCH --time=0-04:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --job-name=qtmerge
#SBATCH --array=1-20
#SBATCH --mem=64GB
#SBATCH --partition=secondary
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --mail-user=syt3@illinois.edu
#SBATCH --mail-type=BEGIN,FAIL,END
#SBATCH --begin=now

source ~/.bashrc
micromamba activate treemerge

# ==========================================
# Parameters
# ==========================================
N=2000

# max_subset_size = min(10, floor(sqrt(N)))
SQRT_N=$(python3 -c "import math; print(int(math.sqrt(${N})))")
if [ "$SQRT_N" -gt 10 ]; then
    MAX_SUBSET_SIZE=$SQRT_N
else
    MAX_SUBSET_SIZE=10
fi

TASK_ID=$(printf "%02d" "${SLURM_ARRAY_TASK_ID}")

GENE_TREES=~/scratch/treemerge-data/${N}/${TASK_ID}/estimatedgenetre
OUTDIR=~/scratch/treemerge-data/${N}/${TASK_ID}/treemerge_estimatedgenetre

mkdir -p "$OUTDIR"

echo "Running dataset N=${N}, replicate ${TASK_ID}"
echo "max_subset_size=${MAX_SUBSET_SIZE}"

# ==========================================
# Run pipeline and record runtime/memory
# ==========================================
/usr/bin/time -v \
    -o "${OUTDIR}/time_memory.txt" \
    python run_dtm_pipeline.py \
        --gene_trees "$GENE_TREES" \
        --outdir "$OUTDIR" \
        --max_subset_size "$MAX_SUBSET_SIZE" \
    > "${SLURM_JOB_NAME}_${SLURM_JOB_ID}.oe"
