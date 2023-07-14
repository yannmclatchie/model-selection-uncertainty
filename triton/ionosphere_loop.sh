# define dataset
dataset="/scratch/work/mclatcy1/model-selection-uncertainty/data/datasets/real/ionosphere.RDS"

# run models
for fold in {1..10}
do
	sbatch triton/ionosphere.sh $dataset $fold
done
