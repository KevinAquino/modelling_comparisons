% Generic script here to test the different models under different pre-processing, the first thing here is to only look at using ten subjects,
% Then extend this onto every single subject -- for sure.

% ============================================================

% Here we have the four processing streams that we want to model the responses to
% preprocessing_stream ={'MINIMAL','ICA-AROMA','ICA-AROMA_GSR','ICA-AROMA_DBSCAN'};
preprocessing_stream ={'ICA-AROMA','ICA-AROMA_GSR','ICA-AROMA_DBSCAN'};

% Load in the time series:
load('/Users/kevinaquino/projects/modelling_gusatvo/empirical_data/ten_subjects_4types.mat');
% Load in the structural matrix:
load('/Users/kevinaquino/projects/modelling_gusatvo/empirical_data/exemplarSC.mat');
sc_matrix = C/max(C(:))*0.2;
G = linspace(0,6,20);

% A structure here to capture all of the models. Now of course all of these can't really be 
% solved on one desktop so they will be sent to the cluster to be solved from matlab. Perhaps each
% sent out all at once from matlab.
model_class = struct;

model_class.global = struct;
model_class.global.models = {'HOPF+GLOBAL','DECO+WANG','BTF','NOISY+DEGREE'};
% Possibly for the BTF model we can use the brain dynamics toolkit (to avoid doing the wrong things :-O) we can possibly look at the
% responses in neural space and then look at the BOLD activity with a BOLD forward model.

model_class.globalAndNode = struct;
model_class.globalAndNode.models = {'HOPF+HETEROGENOUS','DECO+WANG+BALANCED'};

model_class.globalAndEdge = struct;
model_class.globalAndEdge.models = {'HOPF+ANEC'};


bigFig = figure('color','white');
counter = 1;
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
		for prepro = preprocessing_stream,
			disp(['Using model:',model{1},' with processing stream:',prepro{1}]);
			load(['results/',model{1},'/',prepro{1},'.mat']);
			if(isempty(FCD))
				G =[];
			
			else
				figure('color','white');plot(G,FCD,'k.-','MarkerSize',24);hold on;plot(G,grandFCcorr,'r.-','MarkerSize',24);
				title([model{1},' using data: ',prepro{1}],'Interpreter','none');
				xlabel('Global Coupling (G)');ylabel('FCorr (red), FCD (black)');
				ylim([0 1]);xlim([0 G(end)]);
				set(gca,'fontSize',18);


				figure(bigFig);
				subplot(4,3,counter);
				plot(G,FCD,'k.-','MarkerSize',24);hold on;plot(G,grandFCcorr,'r.-','MarkerSize',24);
				title([model{1},' using data: ',prepro{1}],'Interpreter','none');
				xlabel('Global Coupling (G)');ylabel('FCorr (red), FCD (black)');
				ylim([0 1]);xlim([0 G(end)]);
				set(gca,'fontSize',18);
				counter = counter+1;
			end
			


			% switch prepro{1}
			% 	case 'MINIMAL'
			% 		time_series = time_series_aparc_processed(:,:,:,4);
			% 	case 'ICA-AROMA'
			% 		time_series = time_series_aparc_processed(:,:,:,1);
			% 	case 'ICA-AROMA_GSR'
			% 		time_series = time_series_aparc_processed(:,:,:,2);
			% 	case 'ICA-AROMA_DBSCAN'
			% 		time_series = time_series_aparc_processed(:,:,:,3);
			% end
				
			
			% % Here now run the models
			% run_network_model(sc_matrix,time_series,G,model{1},prepro{1});

			% Now here need to add the actual model with parameters that are equivalent
			% Have to work out the inputs really
			
			% We don't actually want to save outputs in a generic function because we want the capability to be general enough to
			% use the cluster.


			% The function has to have the following inputs: the time series of the empirical data, the SC matrix (Which will be global) and then of course the global coupling G.

		end
	end
end

set(bigFig,'PaperPosition',[0.25 0.25 100 40],'InvertHardCopy','off');
print(bigFig,['Comparison_figure.png'],'-dpng','-opengl');


% The idea here is to set it up for each type of model there will be 6 types of models, and run over 4 different preprocessing types.
% Cool stuff to consider:
% 	- Variation of the parameter estimates etc with preprocessing types


% First off for the global models it should be easy as it already has all the stuff we need