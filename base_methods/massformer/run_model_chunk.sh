#!/bin/bash
#SBATCH --job-name=predict_smis
#SBATCH --output=logs/MF_predict_chunk_%a.out
#SBATCH --error=logs/MF_predict_chunk_%a.err
#SBATCH --array=2-600
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00

python scripts/run_inference.py -c config/demo/demo_eval.yml -s predictions/GNPSnew.csv -o predictions/GNPSnew_chunk${SLURM_ARRAY_TASK_ID}.csv --chunk-id ${SLURM_ARRAY_TASK_ID} --num-chunks 600 -d -1

