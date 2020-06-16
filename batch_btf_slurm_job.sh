#!/bin/env bash

#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=24000
#SBATCH -A kg98
#SBATCH --array=1-108
#SBATCH --mail-user=kevin.aquino@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END


# Note: G_ind is passed on through another script.
Run=${SLURM_ARRAY_TASK_ID}

echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} / 108 For $G_ind ----- "
echo -e "\t\t\t --------------------------- \n"



matlab -nodisplay -r "setup_model_comparisons;batch_BTF_run('${G_ind}','${Run}'); exit"
