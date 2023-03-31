#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=($(seq 0 0.5 1))

# run models

for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  dataset="/scratch/cs/bayes_ave/yann/model-selection-uncertainty/data/datasets/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
