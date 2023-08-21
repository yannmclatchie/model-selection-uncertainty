#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1
#SBATCH --output=slurm/%A_%a.out

# run models
for current_dataset in {1..100}
do
    singularity run -B /scratch,/m,/l,/share /scratch/cs/bayes_ave/stan-triton.sif Rscript ./R/many-irrelevant/high_risk_experiment.R $1 $current_dataset $2
done
