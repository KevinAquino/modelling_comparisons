#!/bin/env bash

#SBATCH --job-name=BOLD
#SBATCH --account=kg98
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=4:00:00
#SBATCH --mail-user=kevin.aquino@monash.edu
#SBATCH --mail-type=FAIL
#SBATCH --mail-type=END
#SBATCH --export=ALL
#SBATCH --mem-per-cpu=24000
#SBATCH -A kg98


module load matlab/r2016a
echo $cm
matlab -nodisplay -r "run_multiple_BTF('${cm}'); exit"