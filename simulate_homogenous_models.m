% Generic script here to run homogenous models
% These models run independently of functional imaging data
% and just need to be run with an input of structual imaging matrix
% ============================================================

% Here we have the four processing streams that we want to model the responses to
fmri_dataset = 'UCLA'
test=0;

switch fmri_dataset
	case 'UCLA'
		% Load up the previously calcualted QC-FC metrics (From Aquino et al. 2019)
		% load('~/Documents/fMRIClusterCorrect/stats/CNP_eps_08.mat');
		load('/Users/aquino/projects/modelling_comparisons/figures_ms/CNP_eps_08.mat')
		subjects_restricted=metadata.ParticipantID;
		% Now do the subject pruning:
		subjects_restricted=subjects_restricted(setdiff(1:length(subjects_restricted),[5 7]));

		% Load in the time series: (with accompanying structural connectivity matrix)
		load('empirical_data/UCLA_time_series_four_groups.mat');
		C=ADJ_average;
		% Find the number of subjects that have passed QC-FC:
		subject_ids_restricted=find(ismember(metadata.participant_id,subjects_restricted));
		time_series = time_series(:,:,subject_ids_restricted,:);

		n_subjects=size(time_series,3);
		n_volumes=size(time_series,1);

		empirical_params.sc_matrix_name = 'HCP_APARC_ASEG';
		empirical_params.TR=2;
		% % Remove out the bad subjects (parsing of the subject data here)
		% time_series=time_series(:,:,setdiff(Nsubs,badSub),:);
		% time_series=time_series(:,:,1:100,:);
	case 'GenCog'
	case 'HCP'
end

% Empirical specification
empirical_params.normFactor = 0.2; % This normalization factor is what the MAX C is
empirical_params.fmri_dataset=fmri_dataset;
empirical_params.sc_matrix = C/max(C(:))*empirical_params.normFactor;
empirical_params.time_series=time_series;

% General patterns here - level of different Global coupling.
simulation_params.G = linspace(0,4.5,20);
% Specify the number of runs:
simulation_params.N_RUNS=n_subjects;
simulation_params.N_FRAMES=n_volumes;


% A structure here to capture all of the models. Now of course all of these can't really be 
% solved on one desktop so they will be sent to the cluster to be solved from matlab. Perhaps each
% sent out all at once from matlab.

% MODELS={'DECO+WANG+BALANCED','HOPF+GLOBAL','NOISY+DEGREE'}
% MODELS={'NOISY+DEGREE'};
MODELS={'BTF'};


if(test)
	% Test case variables
	simulation_params.G = [1 2 3];
	simulation_params.N_RUNS=5;
end

for nm=1:length(MODELS),
	% Run the models:
	model_time_series=permute(time_series(:,:,:,1),[2 1 3]);
	simulation_params.MODEL=MODELS{nm}
	% Note: the Hopf Model needs a description of the mean w_f, and the NDM needs a time series
	% They are not nescessary in detail - but calculated here based on time series.	
	run_network_model(empirical_params,simulation_params);
end

