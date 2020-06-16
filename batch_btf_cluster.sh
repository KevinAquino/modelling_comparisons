# !/bin/bash

# This script is to set it all up (running on matlab 2018b)
module load matlab/r2018a

# Now run this for each G value


for g in `seq 1 20`;
	do
		export G_ind=$g		
		echo "Sending G_ind: "$G_ind "To Massive"
		sbatch --job-name="BTF_"$G_ind batch_btf_slurm_job.sh
    done


