#!/bin/bash

# Number of subjects
SUBJECTS=146; 

for i in `seq 1 $SUBJECTS`;
	do
		export subjID=$i
		sbatch --job-name=job$i hopf_job.script
	done

