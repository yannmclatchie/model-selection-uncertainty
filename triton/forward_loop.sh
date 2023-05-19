#!/bin/bash

iters=({1..2})
# ns=('100' '200' '400')
ns=('100')
# rhos=('0' '0.5' '0.9')
rhos=('0')
priors=('normal' 'r2d2' 'projpred')

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
