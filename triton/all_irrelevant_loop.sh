#!/bin/bash

exec='stan/K_model_bias'
Ks=({2..106..8})

# run models
for K in "${Ks[@]}"
do
	dataset="data/datasets/all-irrelevant/all_irrelevant_K${K}_datasets.RDS"
	sbatch triton/all_irrelevant_run.sh $dataset $exec
done
