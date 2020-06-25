% Function here simulates the model with the given sc matrix and then provides some outputs that are then saved
% in a file that is saved for different types


function run_network_model(empirical_params,simulation_params)
% Now when this is solved for each element of G 
% it is then save seperately, this is to make sure 
% it is consisent once we go to the cluster level. 


% now the first thing one has to do is to check the model
% Another thing -- we have to consider actually adding GSR to the models as well. 
fmri_dataset=empirical_params.fmri_dataset;

switch simulation_params.MODEL
	case 'HOPF+GLOBAL',
		% Here it is the hopf global model i.e. setting a=-0.01.
		folder = [fmri_dataset,'/results/HOPF+GLOBAL'];
		[ts_simulated] = run_hopf_homogenous_MODEL(empirical_params,simulation_params);

	case 'BTF',
		% Here will have to grab the BTF model and then look at stuff
		folder = [fmri_dataset,'/results/BTF'];		
		[ts_simulated] = run_BTF_model(empirical_params,simulation_params);

	case 'HOPF+HETEROGENOUS',
		% Here look at the heterogenous 
		folder = [fmri_dataset,'/results/HOPF+HETEROGENOUS'];
		[ts_simulated,simulation_params] = run_hopf_heterogenous_MODEL(empirical_params,simulation_params);		
		% Will have to instead have something more uniform, but for now keep this structure.

	case 'DECO+WANG+BALANCED',
		% Here it is the balanced EI model which is essentially the noise degree model
		%  -- This is esentiall the noise-degree model but for now keep it seperate
		folder = [fmri_dataset,'/results/DECO+WANG+BALANCED'];
		[ts_simulated] = run_BEI_MODEL(empirical_params,simulation_params);
		
	case 'HOPF+ANEC',
		% Here it is the HOPF model with ANEC 
		% Will have to work through this as well. 
		% This here works now -- I think however, will have to implement this version of the FCD as well and add it all together up to the final FCD.
		folder = [fmri_dataset,'/results/HOPF+ANEC'];
		[ts_simulated,simulation_params] = run_hopf_ANEC_MODEL(empirical_params,simulation_params);	

	case 'NOISY+DEGREE',
		% Here is simply the noise and degree model
		folder = [fmri_dataset,'/results/NOISY+DEGREE'];
		% Maybe here run it then save it appropriately.
		ts_simulated = run_noisy_degree_model(empirical_params,simulation_params);
end


system(['mkdir -p ',folder]);
% If no batching was induced/setup make sure it doesnt get tested
if(~isfield(simulation_params,'batch'))
	simulation_params.batch.on=0;
end


if(simulation_params.batch.on)
	G_index = simulation_params.batch.g_ind;
	N_run_model = simulation_params.batch.run;
	save([folder,'/','G_ind_',num2str(G_index),'_RUN_',num2str(N_run_model),'simulation','.mat'],'simulation_params','ts_simulated');
else
	save([folder,'/','simulation','.mat'],'simulation_params','ts_simulated');
end



% There is a slight issue -- FCD can be calculated in different ways -- either using the phases (as Gustavo is using now) or by standard FCD measures, will have to work out later.