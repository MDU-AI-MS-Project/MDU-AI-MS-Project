#!/bin/bash
#SBATCH --job-name=predict_smis
#SBATCH --output=logs/predict_chunk_%a.out
#SBATCH --error=logs/predict_chunk_%a.err
#SBATCH --array=0-600
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00

gen_model="quickstart/iceberg/models/canopus_iceberg_generate.ckpt"
score_model="quickstart/iceberg/models/canopus_iceberg_score.ckpt"

labels="data/spec_datasets/GNPSnew_ion.tsv"
python src/ms_pred/dag_pred/predict_smis.py \
    --batch-size 16 \
    --max-nodes 100 \
    --gen-checkpoint $gen_model \
    --inten-checkpoint $score_model \
    --save-dir quickstart/iceberg/out_ion \
    --dataset-labels $labels \
    --num-workers 0 \
    --chunk-id ${SLURM_ARRAY_TASK_ID} \
    --num-chunks 600

    # --gpu
    # --binned-out
