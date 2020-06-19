% Bespoke method for the heterogenous model, this is because its not as general as the first bits of code.

% Using GSR and GSR'd model.

MODELS={'HOPF+HETEROGENOUS'};

simulation_params.preproMethod=2;
simulation_params.UseGSR=1;

for nm=1:length(MODELS),
	% Run the models:
	model_time_series=permute(time_series(:,:,:,:),[2 1 3 4]);
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
	model_time_series=permute(time_series(:,:,:,:),[2 1 3 4]);
	simulation_params.MODEL=MODELS{nm}
	% Note: the Hopf Model needs a description of the mean w_f, and the NDM needs a time series
	% They are not nescessary in detail - but calculated here based on time series.	
	run_network_model(empirical_params,simulation_params);
end

!cp -r UCLA/results/HOPF+HETEROGENOUS UCLA/results/HOPF+HETEROGENOUS+DICER