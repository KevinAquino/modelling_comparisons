
t1t2_nifti = './empirical_data/Q1-Q6_RelatedParcellation210_AverageT1wDividedByT2w.nii.gz';
left_label='/Applications/freesurfer/subjects/fsaverage/label/lh.aparc.annot';
right_label='/Applications/freesurfer/subjects/fsaverage/label/rh.aparc.annot';

aparc_aseg_epi = '/Users/kevinaquino/projects/modelling_gustavo/empirical_data/fsa_aparcaseg.nii.gz';

setenv('FREESURFER_HOME',['/Applications/freesurfer/']);
[data_ts_lh,data_ts_rh] = CBIG_RF_projectMNI2fsaverage(t1t2_nifti);
% can use the CBIG tools here to transform the data into the surface space, not perfect but okay for our purposes - do we also want sub-cortex though? not sure 
data_ts_lh = parcelTimeSeries(data_ts_lh.',left_label);
data_ts_rh = parcelTimeSeries(data_ts_rh.',right_label);        


% Also now look at sub-cortex.
% Grab the map, do a transform then map it onto here.
% Will have to look at some stuff first        
t1t2_text = ['~/t1t2_on_fsa.txt'];
% Bottom works outside
% unix_command = ['fslmeants -i ','t1t2_aparc_aseg.nii.gz',' --label=',aparc_aseg_epi,' -o ',t1t2_text];
% system(unix_command);   

% Now to do this 


scTimeSeries = dlmread(t1t2_text);
left_subcortex=scTimeSeries(:,[10,11,12,13,17,18,26],:)';
right_subcortex=scTimeSeries(:,[49,50,51,52,53,54,58],:)';

allInds = setdiff(1:35,4);
data_ts_lh = [data_ts_lh(allInds,:);left_subcortex];
data_ts_rh = [data_ts_rh(allInds,:);right_subcortex];

data_total = [data_ts_lh;data_ts_rh];        

% Just cortex
allInds = setdiff(1:35,4);
data_ts_lh = [data_ts_lh(allInds,:)];
data_ts_rh = [data_ts_rh(allInds,:)];
t1t2Cortex = [data_ts_lh;data_ts_rh];        