#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=('0.001' '0.002' '0.003' '0.006' '0.01' '0.018' '0.032' '0.056' '0.1' '0.178' '0.316' '0.562' '1.00' '1.778' '3.162' '5.623' '10.00')

# run models
for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  	dataset="data/datasets/many-irrelevant/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
