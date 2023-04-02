#!/bin/bash

exec='stan/joint_predictive_lin_reg'
n='100'

# run models
dataset="data/datasets/jointpred_n{$n}_datasets.RDS"
sbatch triton/joint_pred_run.sh $dataset $exec
