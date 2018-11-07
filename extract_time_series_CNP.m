% load in the data from massive
% do the transformation here to aparc
% save the time series

% do this for all levels of processing
% do first across all subjects

base_folder_string = '/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/';

% addpath('/usr/local/freesurfer/devel/matlab/');
addpath('/home/kaqu0001/projects/eigen_decomposition');


fid = fopen('UCLA_data.txt');
tline = fgetl(fid);
counter = 1;
while ischar(tline)    
	subject_list{counter} = tline;	
    tline = fgetl(fid);
    counter = counter+1;    
end
fclose(fid);
% subject_list = subjectr_list(1);

% Look at data with ICA-AROMA GSR, and then look at ICA-AROMA DBSCAN
% Note no explicit physiological regression on either

analyses = {'ICA-AROMA+2P_GSR'};
tmpdir = '/home/kaqu0001/kg98_scratch/kevo/tmpdir';
system('mkdir -p /home/kaqu0001/kg98_scratch/kevo/tmpdir');
addpath([getenv('FREESURFER_HOME'),'/matlab']);
setenv('TMPDIR',tmpdir);
setenv('SUBJECTS_DIR',[base_folder_string,'freesurfer/']);

tic;
cd(tmpdir);

for analysis_type = 1:length(analyses)
    for subject=1:length(subject_list),
		subjectName = subject_list(subject);
        epi=[base_folder_string,'/fmriprep/',subjectName,'/func/',subjectName,'_task-rest_bold_space-T1w_variant-AROMAnonaggr+2P_preproc.nii.gz'];
        outputName=[base_folder_string,'/fmriprep/',subjectName,'/func/',subjectName,'_task-rest_bold_space-self-AROMAnonaggr+2P_preproc'];
        % Now perform GSR to this (or should we perform meanGMTR)
        % First just wrtite the code without
		unix_string = ['mri_vol2surf --mov ',epi,' --trgsubject ',subjectName,' --hemi lh --o ',outputName,'.lh.nii.gz'];
		unix(unix_string);
		unix_string = ['mri_vol2surf --mov ',epi,' --trgsubject ',subjectName,' --hemi rh --o ',outputName,'.rh.nii.gz'];
		unix(unix_string);		
        data_ts_lh = MRIread([outputName,'.lh.nii.gz']);data_ts_lh = squeeze(data_ts_lh.vol);
        data_ts_rh = MRIread([outputName,'.rh.nii.gz']);data_ts_rh = squeeze(data_ts_rh.vol);
        % data_ts_rh = MRIread(['timeseries' num2str(subject) '.rh.nii']);data_ts_rh = squeeze(data_ts_rh.vol);
        
        data_ts_lh = parcelTimeSeries(data_ts_lh,'~/kg98/kevo/fsaverage/label/lh.aparc.annot');
        data_ts_rh = parcelTimeSeries(data_ts_rh,'~/kg98/kevo/fsaverage/label/rh.aparc.annot');
         % data_ts_lh = parcelTimeSeries(data_ts_lh,'/home/kaqu0001/projects/eigen_decomposition/lh.HCP-MMP1.annot');
         % data_ts_rh = parcelTimeSeries(data_ts_rh,'/home/kaqu0001/projects/eigen_decomposition/rh.HCP-MMP1.annot');
        
        % Getting rid of corpus callosum in the parcellation
       % data_ts_lh = data_ts_lh(setdiff(1:35,4),:);
       % data_ts_rh = data_ts_rh(setdiff(1:35,4),:);      
        data_total = [data_ts_lh;data_ts_rh];        
        correlation_type(:,:,subject,analysis_type) = corr(data_total.');
        time_series(:,:,subject,analysis_type) = data_total;

	end
end
toc;

save('AROMA+2P+GSR');

