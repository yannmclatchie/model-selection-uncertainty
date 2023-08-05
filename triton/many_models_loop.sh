#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=("0.001" "0.002" "0.003" "0.004" "0.007" "0.012" "0.019" "0.032" "0.052" "0.085" "0.139" "0.227" "0.372" "0.609" "0.998")

# run models
for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  	dataset="data/datasets/many-irrelevant/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
