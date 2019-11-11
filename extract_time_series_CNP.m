% load in the data from massive
% do the transformation here to aparc
% save the time series

% do this for all levels of processing
% do first across all subjects

base_folder_string = '/scratch/kg98/Linden/ResProjects/GSR/data/CNP/derivatives/';

% addpath('/usr/local/freesurfer/devel/matlab/');
addpath('/home/kaqu0001/projects/eigen_decomposition');


fid = fopen('UCLA_data_all.txt');
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

% analyses = {'ICA-AROMA+2P_GSR'};
tmpdir = '/home/kaqu0001/kg98_scratch/kevo/tmpdir';
system('mkdir -p /home/kaqu0001/kg98_scratch/kevo/tmpdir');
addpath([getenv('FREESURFER_HOME'),'/matlab']);
setenv('TMPDIR',tmpdir);
setenv('SUBJECTS_DIR',[base_folder_string,'freesurfer/']);
left_label='/home/kaqu0001/kg98/kevo/fsaverage/label/lh.aparc.annot';
right_label='/home/kaqu0001/kg98/kevo/fsaverage/label/rh.aparc.annot';
tic;
% cd(tmpdir);
badSub = [];
for subject=1:length(subject_list),
    disp('========================================================================================================================')
    disp(['=======================================SUBJECT',subject_list{subject},'================================================'])
    disp('========================================================================================================================')

    try 
		subjectName = subject_list{subject};
        epi{1} = [base_folder_string,'/fmriprep/',subjectName,'/dbscan/',subjectName,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf.nii.gz'];        
        epi{2} = [base_folder_string,'/fmriprep/',subjectName,'/dbscan/',subjectName,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P+GMR_detrended_hpf.nii.gz'];
        epi{3} = [base_folder_string,'/fmriprep/',subjectName,'/dbscan/',subjectName,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_dbscan.nii.gz'];

        aparc_aseg=[base_folder_string,'/fmriprep/',subjectName,'/anat/',subjectName,'_T1w_space-MNI152NLin2009cAsym_preproc_aparcaseg_roi.nii.gz'];
        aparc_aseg_epi=[base_folder_string,'/fmriprep/',subjectName,'/anat/',subjectName,'_bold_space-MNI152NLin2009cAsym_preproc_aparcaseg_roi.nii.gz'];
        % Here we are xforming the APARC_ASEG
        unix_command = ['mri_vol2vol --mov ',aparc_aseg,' --targ ',epi{1},' --regheader --o ',aparc_aseg_epi];
        system(unix_command);
        % First just wrtite the code without
        for analysis_type=1:3,
            clear data_ts_lh data_ts_rh
            [data_ts_lh,data_ts_rh] = CBIG_RF_projectMNI2fsaverage(epi{analysis_type});
    		% can use the CBIG tools here to transform the data into the surface space, not perfect but okay for our purposes - do we also want sub-cortex though? not sure 
            data_ts_lh = parcelTimeSeries(data_ts_lh.',left_label);
            data_ts_rh = parcelTimeSeries(data_ts_rh.',right_label);        


            % Also now look at sub-cortex.
            % Grab the map, do a transform then map it onto here.
            % Will have to look at some stuff first        
            timeSeries_text = [base_folder_string,'/fmriprep/',subjectName,'/dbscan/','tser.txt'];
            unix_command = ['fslmeants -i ',epi{analysis_type},' --label=',aparc_aseg_epi,' -o ',timeSeries_text];
            system(unix_command);

            scTimeSeries = dlmread(timeSeries_text);
            left_subcortex=scTimeSeries(:,[10,11,12,13,17,18,26],:)';
            right_subcortex=scTimeSeries(:,[49,50,51,52,53,54,58],:)';

            allInds = setdiff(1:35,4);
            data_ts_lh = [data_ts_lh(allInds,:);left_subcortex];
            data_ts_rh = [data_ts_rh(allInds,:);right_subcortex];
            
            data_total = [data_ts_lh;data_ts_rh];        

            % Getting rid of Corpus Callosum in the aparc parcellation

            % Saving the data.
            correlation_type(:,:,subject,analysis_type) = corr(data_total.');
            time_series(:,:,subject,analysis_type) = data_total;
        end
    catch
        % Don't do anything but here save the bad subjects
        badSub = [badSub,subject];            
    end
end
toc;

save('UCLA_ALL_SUBS_updated_DiCER');

