# define dataset
dataset="/scratch/work/mclatcy1/model-selection-uncertainty/data/datasets/real/sonar.RDS"

# run models
for fold in {1..10}
do
	sbatch triton/sonar.sh $dataset $fold
done
