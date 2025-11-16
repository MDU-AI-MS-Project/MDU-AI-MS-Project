#!/bin/bash
#SBATCH --job-name=predict_smis
#SBATCH --output=logs/rassp_predict_chunk_%a.out
#SBATCH --error=logs/rassp_predict_chunk_%a.err
#SBATCH --array=2-600
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00

USE_CUDA=0 python rassp/run_rassp.py rassp/sample_data/GNPSnew.txt rassp/sample_data/GNPSnew_rassp.json --file-type smiles --model-name "FormulaNet" --no-gpu --output-file-type json --chunk-id ${SLURM_ARRAY_TASK_ID} --num-chunks 600

