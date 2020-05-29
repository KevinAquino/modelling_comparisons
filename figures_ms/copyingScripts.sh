# # !bin/bash
# # Get a range of subjects
# # No real GS
# subjects[0]=sub-10347
# subjects[1]=sub-10668
# subjects[2]=sub-10746

# # High Global signal
# subjects[3]=sub-10206
# # subjects[4]=sub-10356
# subjects[4]=sub-10958
# subjects[5]=sub-10912

# # Mixed High and low signal
# subjects[6]=sub-10189
# subjects[7]=sub-10376
# subjects[8]=sub-10934


# subjects[0]=sub-10570
# subjects[1]=sub-10356
# subjects[2]=sub-10575

# subjects[0]=sub-10217
# subjects[0]=sub-10987

# subjects[0]=sub-10958
# subjects[1]=sub-10492


# subjects[0]=sub-10575

subjects[0]=sub-10487

# for i in `seq 0 8`;
for i in `seq 0 0`;
do 
	save_folder=${subjects[i]}"/"
	subject=${subjects[i]}
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_dbscan.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_dbscan_liberal_regressors.tsv" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P+GMR_detrended_hpf.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_bold_space-MNI152NLin2009cAsym_dtissue_masked.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_bold_space-MNI152NLin2009cAsym_dtissue_masked_clusterorder.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/dbscan/$subject"_bold_space-MNI152NLin2009cAsym_gm_mask_gsordering_tissue.nii.gz" $save_folder
	scp kaqu0001@m3.massive.org.au:/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/fmriprep/$subject/func/$subject"_task-rest_bold_confounds.tsv" $save_folder	
	fslmaths $subject"_bold_space-MNI152NLin2009cAsym_dtissue_masked_clusterorder.nii.gz" $subject"_bold_space-MNI152NLin2009cAsym_dtissue_masked_clusterorder.nii.gz" -odt float
done

