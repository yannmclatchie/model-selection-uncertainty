#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=("0.008" "0.010" "0.014" "0.019" "0.025" "0.034" "0.045" "0.060" "0.081" "0.108" "0.145" "0.194" "0.259" "0.347" "0.465" "0.622" "0.833" "1.116" "1.494" "2.000")

# run models
for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  	dataset="data/datasets/many-irrelevant/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
