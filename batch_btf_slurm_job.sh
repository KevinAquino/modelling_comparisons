#!/bin/env bash

#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=24000
#SBATCH -A kg98
#SBATCH --array=1-2160
#SBATCH --job-name BTF

FILE='BTF_batch_list.txt'

BATCH_STRING=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${FILE})

G_ind=$(echo $BATCH_STRING | awk -F, '{print $1}')
Run=$(echo $BATCH_STRING | awk -F, '{print $2}')


echo -e "\t\t\t --------------------------- "
echo -e "\t\t\t ----- ${SLURM_ARRAY_TASK_ID} ${subject} ----- "
echo -e "\t\t\t --------------------------- \n"



matlab -nodisplay -r "setup_model_comparisons;batch_BTF_run('${G_ind}','${Run}'); exit"