#!/bin/bash
#SBATCH --time=02:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=1G
#SBATCH --array=1

# run models
for current_dataset in {1..100}
do
    singularity run -B /scratch,/m,/l,/share docker://andrjohns/stan-triton Rscript ./R/fit_many_models.R $1 $current_dataset $2
done