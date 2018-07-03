% load in the data from massive
% do the transformation here to aparc
% save the time series

% do this for all levels of processing
% do first across all subjects

base_folder_string = '/home/kaqu0001/kg98_scratch/Linden/ResProjects/rfMRI_denoise/UCLA/data/';
addpath('/usr/local/freesurfer/devel/matlab/');
addpath('/home/kaqu0001/projects/eigen_decomposition');


fid = fopen('UCLA_subset.txt');
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

analyses = {'ICA-AROMA_GSR','ICA-AROMA_DBSCAN'};
analysis_folders = {'','/home/kaqu0001/kg98_scratch/kevo/GSR_data/dibwrat'};
tmpdir = '/home/kaqu0001/kg98_scratch/kevo/GSR_data/dibwrat/tmpdir';
addpath([getenv('FREESURFER_HOME'),'/matlab']);
setenv('TMPDIR',tmpdir);

tic;
cd(tmpdir);

for analysis_type = 1:length(analyses)
    for subject=1:length(subject_list),

    	% first for the GSR data:
    	if(analysis_type == 1)
    		epi_tseries = [base_folder_string,subject_list{subject},'/func/prepro/ICA-AROMA_output/ica_dibwrat',subject_list{subject},'_task-rest_bold.nii'];
    	elseif(analysis_type == 2)
    		epi_tseries = [analysis_folders{analysis_type},'dbscan_corrected_ica_dibwrat',subject_list{subject},'_task-rest_bold.nii'];    		
    	end

		
		% first copy the data locally
		unix_string = ['cp ',epi_tseries,' ',tmpdir,'/epi_prepro.nii.gz'];
		unix(unix_string);

		unix_string = ['mri_vol2surf --mov epi_prepro' num2str(subject) '.nii.gz --mni152reg --trgsubject fsaverage --hemi lh --o timeseries' num2str(subject) '.lh.nii.gz'];
		unix(unix_string);
		unix_string = ['mri_vol2surf --mov epi_prepro' num2str(subject) '.nii.gz --mni152reg --trgsubject fsaverage --hemi rh --o timeseries' num2str(subject) '.rh.nii.gz'];
		unix(unix_string);		
        data_ts_lh = MRIread(['timeseries' num2str(subject) '.lh.nii']);data_ts_lh = squeeze(data_ts_lh.vol);
        data_ts_rh = MRIread(['timeseries' num2str(subject) '.rh.nii']);data_ts_rh = squeeze(data_ts_rh.vol);
        
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
		% then have to move it into freesurver
		% unix_string = 'mri_vol2vol --mov epi_prepro.nii.gz --targ $SUBJECTS_DIR/fsaverage/mri/orig.mgz --regheader --o xformed_tseries.nii';
		% unix(unix_string)
		% data = MRIread('xformed_tseries.nii');
% 		keyboard
	end
end
toc;

save saved_data_noBPF_HCP

% then do for different ACMs -- use Stuart's ones
% use distance ones


