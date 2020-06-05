% Grab the T1w/T2w Ration for the parcellation - to see if this, instead of degree, can also be used for the noisy model.


% Add CBIG tools, add freesurfer and fsl to the environment:
addpath(genpath('~/Documents/CBIG/stable_projects/registration/Wu2017_RegistrationFusion/'));
setenv('FREESURFER_HOME',['/Applications/freesurfer/']);
setenv('PATH',['/usr/local/fsl/bin/:',getenv('PATH')]);
setenv('FSLOUTPUTTYPE','NIFTI_GZ');

left_label  = '/Applications/freesurfer/subjects/fsaverage/label/lh.aparc.annot';
right_label = '/Applications/freesurfer/subjects/fsaverage/label/rh.aparc.annot';


t1t2='/Users/kevinaquino/projects/modelling_gustavo/empirical_data/Q1-Q6_RelatedParcellation210_AverageT1wDividedByT2w.nii.gz';

% Grab cortical t1t2ratio:
[data_lh,data_rh] = CBIG_RF_projectMNI2fsaverage(t1t2);
data_lh = parcelTimeSeries(data_lh.',left_label);
data_rh = parcelTimeSeries(data_rh.',right_label);     


% Grab subcortical t1t2 ratio:
aparc_aseg='empirical_data/aseg_mni.nii.gz';

timeSeries_text = ['empirical_data/t1t2values.txt'];
unix_command = ['fslmeants -i ','empirical_data/t1t2_in_aseg.nii.gz',' --label=',aparc_aseg,' -o ',timeSeries_text];
system(unix_command);

scTimeSeries = dlmread(timeSeries_text);
left_subcortex=scTimeSeries(:,[10,11,12,13,17,18,26],:)';
right_subcortex=scTimeSeries(:,[49,50,51,52,53,54,58],:)';



allInds = setdiff(1:35,4);
data_lh = [data_lh(allInds);left_subcortex];
data_rh = [data_rh(allInds);right_subcortex];

data_total = [data_lh;data_rh];  