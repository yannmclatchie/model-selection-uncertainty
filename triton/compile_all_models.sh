#!/bin/bash
#SBATCH --cpus-per-task=4
#SBATCH --mem=4G
#SBATCH --array=1
#SBATCH --time=00:10:00

singularity run -B /scratch,/m,/l,/share /scratch/cs/bayes_ave/stan-triton.sif Rscript ./R/compile_models.R "/scratch/work/mclatcy1/model-selection-uncertainty/stan"
