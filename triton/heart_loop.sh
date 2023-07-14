# define dataset
dataset="/scratch/work/mclatcy1/model-selection-uncertainty/data/datasets/real/heart.RDS"

# run models
for fold in {1..10}
do
	sbatch triton/heart.sh $dataset $fold
done
