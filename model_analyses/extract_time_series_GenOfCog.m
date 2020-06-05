% load in the data from massive
% do the transformation here to aparc
% save the time series

% do this for all levels of processing
% do first across all subjects

base_folder_string = '/scratch/kg98/kristina/Projects/GenofCog/derivatives/dicer/';

% addpath('/usr/local/freesurfer/devel/matlab/');
addpath('/home/kaqu0001/projects/eigen_decomposition');
addpath(genpath('/home/kaqu0001/projects/CBIG/stable_projects/registration/Wu2017_RegistrationFusion/'));


fid = fopen('GenCog_data_all.txt');
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
% setenv('SUBJECTS_DIR',[base_folder_string,'freesurfer/']);
left_label='/home/kaqu0001/kg98/kevo/fsaverage/label/lh.aparc.annot';
right_label='/home/kaqu0001/kg98/kevo/fsaverage/label/rh.aparc.annot';
tic;
% cd(tmpdir);
badSub = [];
all_sample{1}=[1:73];
all_sample{2}=[74:147];
all_sample{3}=[148:220];
all_sample{4}=[221:293];
all_sample{5}=[294:366];
all_sample{6}=[367:440];

dicer_list=all_sample{listNum};
% for subject=1:length(subject_list),
for subject=dicer_list,
    disp('========================================================================================================================')
    disp(['=============================SUBJECT ',subject_list{subject},'  ',num2str(subject),'/440================================================'])
    disp('========================================================================================================================')
    disp(['Total number of bad subjects: ',num2str(length(badSub))]);
    try 
		subjectName = subject_list{subject};
        epi{1} = [base_folder_string,subjectName,'/',subjectName,'_filtered_func_data_clean_mni.nii.gz'];
        epi{2} = [tmpdir,'/',subjectName,'_filtered_func_data_clean_mni_GMR.nii.gz'];
        epi{3} = [base_folder_string,subjectName,'/',subjectName,'_filtered_func_data_clean_mni_dbscan.nii.gz'];

        % First do GMR
        % Extract GM signal
        brain_signal_txt=['/home/kaqu0001/kg98_scratch/kevo/tmpdir/',subjectName,'_gm_signal.txt'];
        unix_command=['fslmeants -i ',epi{1},' -o ',brain_signal_txt,' --label=',base_folder_string,subjectName,'/',subjectName,'_dtissue_func.nii.gz'];
        system(unix_command);
        % Now perform GMR
        unix_command=['fsl_regfilt -i ',epi{1},' -o ',epi{2},' -d ',brain_signal_txt,' -f 4 -a'];
        system(unix_command);
        % Here is the aseg that was derived externally: (from within fsaverage/mri.2mm/)
        % mri_vol2vol --mov aseg.mgz --targ /usr/local/fsl/data/standard/MNI152_T1_2mm.nii.gz --o ~/aseg_mni_152.nii.gz --reg reg.2mm.mni152.dat --nearest       
        aparc_aseg_epi=['aseg_mni_152.nii.gz'];
        % % Here we are xforming the APARC_ASEG
        % unix_command = ['mri_vol2vol --mov ',aparc_aseg,' --targ ',epi{1},' --regheader --o ',aparc_aseg_epi];

        % Can use default from atlas from freesurfer outputs (can do in before)

        
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
            timeSeries_text = ['/home/kaqu0001/kg98_scratch/kevo/tmpdir/',subjectName,'_tser.txt'];
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

save(['GenCog_ALL_DiCER',num2str(listNum)]);

