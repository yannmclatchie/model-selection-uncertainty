#!/bin/bash
#SBATCH --time=05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1

# run models
for current_dataset in {1..100}
do
    singularity run -B /scratch,/m,/l,/share /scratch/cs/bayes_ave/stan-triton.sif Rscript ./R/many-irrelevant/fit_many_models.R $1 $current_dataset $2
done
