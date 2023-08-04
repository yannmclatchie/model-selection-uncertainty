#!/bin/bash

exec='stan/K_model_bias'
Ks=('2' '10' '100')
betas=('0.001' '0.00178' '0.00316' '0.00562' '0.01' '0.01778' '0.03162' '0.05623' '0.1' '0.17783' '0.31623' '0.56234' '1.00')

# run models
for K in "${Ks[@]}"
do
	for beta_delta in "${betas[@]}"
	do
	  	dataset="data/datasets/many-irrelevant/many_models_K${K}_beta${beta_delta}_datasets.RDS"
		sbatch triton/many_models_run.sh $dataset $exec
	done
done
