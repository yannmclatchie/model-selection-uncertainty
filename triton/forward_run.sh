#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1

# define dataset file
dataset="/data/datasets/forward_search_dataset.RDS"

# run models
for current_dataset in {1..2}
do
    singularity run -B /scratch,/m,/l,/share docker://andrjohns/stan-triton Rscript ./R/many-irrelevant/fit_many_models.R $dataset $current_dataset
done
