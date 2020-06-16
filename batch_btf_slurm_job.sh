#!/bin/env bash

#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=1:00:00
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=24000
#SBATCH -A kg98

echo $G_ind $Run

matlab -nodisplay -r "setup_model_comparisons;batch_BTF_run('${G_ind}','${Run}'); exit"