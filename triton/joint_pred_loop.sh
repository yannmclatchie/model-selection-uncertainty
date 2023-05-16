#!/bin/bash

exec='stan/joint_predictive_lin_reg'
n='100'
batches=($(seq 1 100))

for batch in "${batches[@]}"
do
    # run models
    dataset="data/datasets/jointpred_n${n}_batch${batch}_datasets.RDS"
    sbatch triton/joint_pred_run.sh $dataset $batch $exec
done