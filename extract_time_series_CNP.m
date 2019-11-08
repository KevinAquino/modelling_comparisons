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

analyses = {'ICA-AROMA+2P_GSR'};
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
        epi{1} = [base_folder_string,'/fmriprep/',subjectName,'/func/',subjectName,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc+2P.nii.gz'];        
        epi{2} = [base_folder_string,'/fmriprep/',subjectName,'/func/',subjectName,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-smoothAROMAnonaggr_preproc+2P+GMR.nii.gz'];
        epi{3} = [base_folder_string,'/fmriprep/',subjectName,'/dbscan/',subjectName,'task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_dbscan.nii.gz'];

        


        % First just wrtite the code without
        for analysis_type=1:3,
            [data_ts_lh,data_ts_rh] = CBIG_RF_ProjectMNI2fsaverage(epi{analysis_type});
    		% can use the CBIG tools here to transform the data into the surface space, not perfect but okay for our purposes - do we also want sub-cortex though? not sure 
            data_ts_lh = parcelTimeSeries(data_ts_lh.',left_label);
            data_ts_rh = parcelTimeSeries(data_ts_rh.',right_label);        
            data_total = [data_ts_lh;data_ts_rh];        

            % Getting rid of Corpus Callosum in the aparc parcellation
            allInds = setdiff(1:70,[4,35+4]);
            data_total = data_total(allInds,:);

            % Also now look at sub-cortex.
            % Grab the map, do a transform then map it onto here.
            

            % Saving the data.
            correlation_type(:,:,subject,analysis_type) = corr(data_total.');
            time_series(:,:,subject,analysis_type) = zscore(data_total(:,4:end),[],2);
        end
    catch
        % Don't do anything but here save the bad subjects
        badSub = [badSub,subject];            
    end
end
toc;

save('AROMA+2P+GSR');

