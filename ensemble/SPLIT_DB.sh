#!/bin/bash
#SBATCH --partition=work
#SBATCH --qos=basic
#SBATCH --array=2-10
#SBATCH --output="/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split%a.out"
#SBATCH --error="/RG/compbio/michaelel/data/splits/GNPSnew_canonical_output_split%a.err"
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --job-name="split_gen"

# Load conda environment (adjust if needed)
source ~/.bashrc
conda activate lsagnev6

# Execute the Python script with current split ID
srun python /RG/compbio/michaelel/data/splits/split_orig_data_distrib_v1.py $SLURM_ARRAY_TASK_ID
