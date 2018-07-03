% Generic script here to test the different models under different pre-processing 

% ============================================================

% Here we have the four processing streams that we want to model the responses to
preprocessing_stream ={'MINIMAL','ICA-AROMA','ICA-AROMA+GSR','ICA-AROMA+DBSCAN'};

% A structure here to capture all of the models. Now of course all of these can't really be 
% solved on one desktop so they will be sent to the cluster to be solved from matlab. Perhaps each
% sent out all at once from matlab.
model_class = struct;

model_class.global = struct;
model_class.global.models = {'HOPF','DECO+WANG','BTF'};
% Possibly for the BTF model we can use the brain dynamics toolkit (to avoid doing the wrong things :-O) we can possibly look at the
% responses in neural space and then look at the BOLD activity with a BOLD forward model.

model_class.globalAndNode = struct;
model_class.globalAndNode.models = {'HOPF','DECO+WANG'};

model_class.globalAndEdge = struct;
model_class.globalAndEdge.models = {'HOPF'};


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
			% Now here need to add the actual model with parameters that are equivalent
			% Have to work out the inputs really
		end
	end
end


% The idea here is to set it up for each type of model there will be 6 types of models, and run over 4 different preprocessing types.
% Cool stuff to consider:
% 	- Variation of the parameter estimates etc with preprocessing types

% Lets see!