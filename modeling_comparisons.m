% Generic script here to test the different models under different pre-processing, the first thing here is to only look at using ten subjects,
% Then extend this onto every single subject -- for sure.

% ============================================================

% Here we have the four processing streams that we want to model the responses to
fmri_dataset = 'UCLA'

switch fmri_dataset
	case 'UCLA'
		% Load in the time series: (with accompanying structural connectivity matrix)
		load('empirical_data/UCLA_time_series_four_groups.mat');
		C=ADJ_average;
		Nsubs=1:size(time_series,3);
		% Remove out the bad subjects (parsing of the subject data here)
		time_series=time_series(:,:,setdiff(Nsubs,badSub),:);
		time_series=time_series(:,:,1:100,:);
	case 'GenCog'
	case 'HCP'
end

preprocessing_stream = noiseOptions;


% Load in the structural matrix:
% load('empirical_data/exemplarSC.mat');


sc_matrix = C/max(C(:))*0.2;
G = linspace(0,10,20);

% A structure here to capture all of the models. Now of course all of these can't really be 
% solved on one desktop so they will be sent to the cluster to be solved from matlab. Perhaps each
% sent out all at once from matlab.
model_class = struct;

model_class.global = struct;
model_class.global.models = {'HOPF+GLOBAL','DECO+WANG','BTF','NOISY+DEGREE'};
% Possibly for the BTF model we can use the brain dynamics toolkit we can possibly look at the
% responses in neural space and then look at the BOLD activity with a BOLD forward model.

model_class.globalAndNode = struct;
model_class.globalAndNode.models = {'HOPF+HETEROGENOUS','DECO+WANG+BALANCED'};

model_class.globalAndEdge = struct;
model_class.globalAndEdge.models = {'HOPF+ANEC'};




% Here set up the global coupling value for each model, this is done 
% ============================================================

% Here now the general script get it all working
for model_class_type = fields(model_class).',
	disp([' Working with: the following classes of models:', model_class_type{1}]);

	% Here just grab the models
	models = getfield(getfield(model_class,model_class_type{1}),'models');
	for model = models;

		disp(['Currently at: ',model]);
		% Now go through each preprocessing stream
		for prepro_num=1:length(preprocessing_stream);
			prepro = preprocessing_stream{prepro_num},
			disp(['Using model:',model{1},' with processing stream:',prepro]);
			model_time_series=permute(time_series(:,:,:,prepro_num),[2 1 3]);

			% Here now run the models
			run_network_model(sc_matrix,model_time_series,G,model{1},prepro,fmri_dataset);

			% Now here need to add the actual model with parameters that are equivalent
			% Have to work out the inputs really
			
			% We don't actually want to save outputs in a generic function because we want the capability to be general enough to
			% use the cluster.


			% The function has to have the following inputs: the time series of the empirical data, the SC matrix (Which will be global) and then of course the global coupling G.

		end
	end
end


% The idea here is to set it up for each type of model there will be 6 types of models, and run over 4 different preprocessing types.
% Cool stuff to consider:
% 	- Variation of the parameter estimates etc with preprocessing types


% First off for the global models it should be easy as it already has all the stuff we need