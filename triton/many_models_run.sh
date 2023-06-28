#!/bin/bash
#SBATCH --time=06:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1
#SBATCH --output=slurm/many-models/%A_%a.out

# load the required modules
# module load r
# module load gcc/11.2.0

# run models
for current_dataset in {1..50}
do
    singularity run -B /scratch,/m,/l,/share /scratch/cs/bayes_ave/stan-triton.sif Rscript ./R/many-irrelevant/fit_many_models.R $1 $current_dataset $2
    # srun Rscript ./R/many-irrelevant/fit_many_models.R $1 $current_dataset $2
done
