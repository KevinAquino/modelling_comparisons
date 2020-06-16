# !/bin/bash

# This script is to set it all up (running on matlab 2018b)
module load matlab/r2018a
for g in `seq 1 20`;
	do
		for r in `seq 1 108`;
			do
				export G_ind=$g
				export Run=$r
				echo "G index: " $G_ind "Run index:" $Run
				sbatch --job-name="BTF_${G_ind}" batch_btf_slurm.job.sh
			done
    done
