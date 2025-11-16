#!/bin/bash
#SBATCH --partition=work
#SBATCH --qos=basic  # Added QoS as basic
#SBATCH --array=108-2106
#SBATCH --output="/RG/compbio/michaelel/data/run_sbatch_files/outerr/mslearnMol_seed_%a.out"
#SBATCH --error="/RG/compbio/michaelel/data/run_sbatch_files/outerr/mslearnMol_seed_%a.err"
#SBATCH --ntasks=1                   # Run a single task
#SBATCH --cpus-per-task=16            # Number of CPU cores per task
#SBATCH --mem=64GB
#SBATCH --job-name="mslearnMol_seed_%a"

### Modules
echo "activate environment"
conda activate lsagnev6

srun python /RG/compbio/michaelel/data/run_sbatch_files/mslearn_distrib_hypsrch_v2.py -w /RG/compbio/michaelel/data/run_sbatch_files -i $SLURM_ARRAY_TASK_ID -j /RG/compbio/michaelel/data/run_sbatch_files/mslearn_input_extended_v2_cleaned_valid.json