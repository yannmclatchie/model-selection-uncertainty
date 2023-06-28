# define dataset
dataset="/scratch/work/mclatcy1/model-selection-uncertainty/data/datasets/real/crime.RDS"

# run models
for fold in {1..10}
do
	sbatch triton/crime.sh $dataset $fold
done
