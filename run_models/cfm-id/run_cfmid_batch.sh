#!/bin/bash
#SBATCH --job-name=cfm_predict
#SBATCH --output=logs/cfm_predict_%a.out
#SBATCH --error=logs/cfm_predict_%a.err
#SBATCH --array=2-600
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=01:00:00


# Define variables
FILE=/RG/compbio/michaelel/cfm-id/GNPSnew.txt           # Input file
OUTPUT_DIR=/RG/compbio/michaelel/cfm-id/cfm_outputs     # Directory for output files
TOTAL_LINES=24943          # Total lines in the file
LINES_PER_JOB=$((TOTAL_LINES / 600))  # Lines each job processes
EXTRA_LINES=$((TOTAL_LINES % 600))   # Remainder lines
START_LINE=0

# Calculate start and end lines for each job
if [ "$SLURM_ARRAY_TASK_ID" -lt "$EXTRA_LINES" ]; then
    START_LINE=$((SLURM_ARRAY_TASK_ID * (LINES_PER_JOB + 1) + 1))
    END_LINE=$((START_LINE + LINES_PER_JOB))
else
    START_LINE=$((EXTRA_LINES * (LINES_PER_JOB + 1) + (SLURM_ARRAY_TASK_ID - EXTRA_LINES) * LINES_PER_JOB + 1))
    END_LINE=$((START_LINE + LINES_PER_JOB - 1))
fi

# Ensure output directory exists
mkdir -p $OUTPUT_DIR

# Process assigned lines
sed -n "${START_LINE},${END_LINE}p" $FILE | while IFS= read -r SMILES; do
    # Determine the global line number for the current line
    LINE_NUM=$((START_LINE + LINE_INDEX))

    # Set output file name
    OUTPUT_FILE=${OUTPUT_DIR}/GNPSnew_cfmpredict_${LINE_NUM}.txt

    # Run cfm-predict for the current SMILES
    singularity exec --bind $(pwd):/cfmid/public cfmid_latest.sif sh -c \
    "cfm-predict '${SMILES}' 0.001 /trained_models_cfmid4.0/[M+H]+/param_output.log /trained_models_cfmid4.0/[M+H]+/param_config.txt 1 ${OUTPUT_FILE}"

    LINE_INDEX=$((LINE_INDEX + 1))
done
