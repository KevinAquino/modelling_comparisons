
% folder='sub-10274/';
% subject='sub-10274';


folder='sub-10280/';
subject='sub-10280';

functype='_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc.nii.gz';
inputScan=[folder,subject,functype];
mask=[folder,subject,'_bold_space-MNI152NLin2009cAsym_dtissue_masked.nii.gz'];

functypeCell={'_task-rest_bold_space-MNI152NLin2009cAsym_preproc.nii.gz', ... 
	'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf.nii.gz', ...
	'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P+GMR_detrended_hpf.nii.gz', ...
	'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_dbscan.nii.gz'};

% functypeCell={'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf.nii.gz','_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P+GMR_detrended_hpf.nii.gz','_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_dbscan.nii.gz'};
funcLabels={'MPP\newline ','AROMA+2P\newline ','AROMA+2P\newline+GMR','AROMA+2P\newline+DiCER'}

% conf = dlmread([subject,'_task-rest_bold_confounds.tsv'],'',2,0);


subjects={subject};


fh=carpetSummary_modelling(subjects,1,functypeCell,funcLabels,folder);

figFolder='/Users/kevinaquino/projects/modelling_gustavo/manuscript/figs/';
saveeps(fh,[figFolder,'/UCLA_noisy_subject_biphasic'],[50 50]);

