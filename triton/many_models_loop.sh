#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=("0.062" "0.075" "0.090" "0.108" "0.130" "0.156" "0.187" "0.224" "0.269" "0.323" "0.387" "0.465" "0.558" "0.669" "0.803" "0.964" "1.157" "1.389" "1.667" "2.000")

# run models
for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  	dataset="data/datasets/many-irrelevant/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
