% Bespoke method for the heterogenous model, this is because its not as general as the first bits of code.

% Using GSR and GSR'd model.
load('figures_ms/UCLA_data_params.mat')
MODELS={'HOPF+HETEROGENOUS'};

simulation_params.preproMethod=2;
simulation_params.UseGSR=1;

for nm=1:length(MODELS),
	% Run the models:	
	simulation_params.MODEL=MODELS{nm}
	% Note: the Hopf Model needs a description of the mean w_f, and the NDM needs a time series
	% They are not nescessary in detail - but calculated here based on time series.	
	run_network_model(empirical_params,simulation_params);
end

% Copy these results into another directory

!cp -r UCLA/results/HOPF+HETEROGENOUS UCLA/results/HOPF+HETEROGENOUS+GSR

% Now try do the same analysis with DiCER and not GSR'ing the model as well

simulation_params.preproMethod=3;
simulation_params.UseGSR=0;

for nm=1:length(MODELS),
	% Run the models:	
	simulation_params.MODEL=MODELS{nm}
	% Note: the Hopf Model needs a description of the mean w_f, and the NDM needs a time series
	% They are not nescessary in detail - but calculated here based on time series.	
	run_network_model(empirical_params,simulation_params);
end

!cp -r UCLA/results/HOPF+HETEROGENOUS UCLA/results/HOPF+HETEROGENOUS+DICER


% Now try do the same analysis with just AROMA (default) and not GSR'ing the model as well

simulation_params.preproMethod=1;
simulation_params.UseGSR=0;

for nm=1:length(MODELS),
	% Run the models:	
	simulation_params.MODEL=MODELS{nm}
	% Note: the Hopf Model needs a description of the mean w_f, and the NDM needs a time series
	% They are not nescessary in detail - but calculated here based on time series.	
	run_network_model(empirical_params,simulation_params);
end

!cp -r UCLA/results/HOPF+HETEROGENOUS UCLA/results/HOPF+HETEROGENOUS+AROMA



% Here do the first steps of Hopf ANEC Have to train on 80% Then of course test it.

MODELS={'HOPF+ANEC'};

simulation_params.preproMethod=2;
simulation_params.UseGSR=1;
simulation_params.MODEL = MODELS{1};
% Run ANEC at least 20 times (could run more times for better cross validation)
load('UCLA/results/HOPF+ANEC/subjects_for_anec.mat');

for n_perm=2:10;%length(rand_perms),
	display(['CV step... RUN ',num2str(n_perm),'.']);
	% Run the models:	
	simulation_params.rand_perms=rand_perms;
	simulation_params.n_perm=n_perm;
	% Now reduce it to only fit ANEC to 80% of the data, that will later be tested on the remaining 20% of the data.
	empirical_params_80pc=empirical_params;
	empirical_params_80pc.time_series = empirical_params.time_series(:,:,rand_perms(n_perm,:),:);
	run_network_model(empirical_params_80pc,simulation_params);
	unix_command=['cp -r UCLA/results/HOPF+ANEC/simulation.mat ','UCLA/results/HOPF+ANEC/simulation_GSR_',num2str(n_perm),'.mat'];
	system(unix_command)
end

% DICER


MODELS={'HOPF+ANEC'};

simulation_params.preproMethod=3;
simulation_params.UseGSR=0;
simulation_params.MODEL = MODELS{1};
% Run ANEC at least 20 times (could run more times for better cross validation)
load('UCLA/results/HOPF+ANEC/subjects_for_anec.mat');

for n_perm=2:10;%length(rand_perms),
	display(['CV step... RUN ',num2str(n_perm),'.']);
	% Run the models:	
	simulation_params.rand_perms=rand_perms;
	simulation_params.n_perm=n_perm;
	% Now reduce it to only fit ANEC to 80% of the data, that will later be tested on the remaining 20% of the data.
	empirical_params_80pc=empirical_params;
	empirical_params_80pc.time_series = empirical_params.time_series(:,:,rand_perms(n_perm,:),:);
	run_network_model(empirical_params_80pc,simulation_params);
	unix_command=['cp -r UCLA/results/HOPF+ANEC/simulation.mat ','UCLA/results/HOPF+ANEC/simulation_DICER_',num2str(n_perm),'.mat'];
	system(unix_command)
end


