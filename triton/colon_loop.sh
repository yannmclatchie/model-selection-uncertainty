# define dataset
dataset="/scratch/work/mclatcy1/model-selection-uncertainty/data/datasets/real/colon.RDS"

# run models
for fold in {1..62}
do
	sbatch triton/colon.sh $dataset $fold
done
