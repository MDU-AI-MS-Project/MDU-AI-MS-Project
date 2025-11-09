head -n 1 /RG/compbio/michaelel/data/run_sbatch_files/*outputMol.csv | head -n 1 > merged_output.csv && tail -n +2 -q /RG/compbio/michaelel/data/run_sbatch_files/*outputMol.csv >> merged_output.csv
