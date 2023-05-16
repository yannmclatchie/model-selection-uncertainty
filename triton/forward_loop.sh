#!/bin/bash

iters=({1..20})
ns=('100' '200' '400')
rhos=('0' '0.5' '0.9')
priors=('normal' 'r2d2')

# run models
for prior_name in "${priors[@]}"
do
	for nobs in "${ns[@]}"
	do
		for rho in "${rhos[@]}"
		do
			for iter in "${iters[@]}"
			do
				dataset="data/datasets/forward/forward_dataset_${nobs}_${rho}.RDS"
				sbatch triton/forward_run.sh $dataset $nobs $rho $prior_name $iter
			done
		done
	done
done
