#!/bin/bash
#SBATCH --job-name=predict_smis
#SBATCH --output=logs/scarf_predict_chunk_%a.out
#SBATCH --error=logs/scarf_predict_chunk_%a.err
#SBATCH --array=0-600
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=02:00:00

thread_model="quickstart/scarf/models/canopus_thread_model.ckpt"
weave_model="quickstart/scarf/models/canopus_weave_model.ckpt"

labels="data/spec_datasets/GNPSnew_ion.tsv"
python src/ms_pred/scarf_pred/predict_smis.py \
    --batch-size 32 \
    --sparse-out \
    --max-nodes 300 \
    --gen-checkpoint $thread_model \
    --inten-checkpoint $weave_model \
    --save-dir quickstart/scarf/out_ion \
    --dataset-labels $labels \
    --num-workers 0 \
    --chunk-id ${SLURM_ARRAY_TASK_ID} \
    --num-chunks 600

    # --gpu
    # --binned-out
