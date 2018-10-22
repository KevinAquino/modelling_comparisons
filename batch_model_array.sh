# !/bin/bash

ModelFile='c_vals.txt'

while IFS=$ModelFile read -r line || [[ -n "$line" ]]; do
    #cd $FMRIPREP_DIR
    # Now look at the subject:
    cm=$line
	export cm
	sbatch --job-name="$line" batch_models.sh
done < "$ModelFile"