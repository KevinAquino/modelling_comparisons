% Function here simulates the model with the given sc matrix and then provides some outputs that are then saved
% in a file that is saved for different types


% Outputs that have to be saved per G (would be good to keep it quite standard)
% The FC matrix 
% The FCD KS-estimate
% The estimated time series as well (which is needed)

function run_model(sc_matrix,time_series,G,model,preproType,fmri_dataset)
% Now when this is solved for each element of G 
% it is then save seperately, this is to make sure 
% it is consisent once we go to the cluster level. 


% now the first thing one has to do is to check the model
% Another thing -- we have to consider actually adding GSR to the models as well. 

switch model	
	case 'HOPF+GLOBAL',
		% Here it is the hopf global model i.e. setting a=-0.01.
		folder = [fmri_dataset,'/results/HOPF+GLOBAL'];
		[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,time_series,G,folder);
	case 'DECO+WANG',
		% Here look at the model what does this use? not the balanced EI model (will have to look)
		folder = [fmri_dataset,'/results/DECO+WANG'];
		FCD = [];ts_simulated_all = [];grandFCcorr = [];
	case 'BTF',
		% Here will have to grab the BTF model and then look at stuff
		folder = [fmri_dataset,'/results/BTF'];
		FCD = [];ts_simulated_all = [];grandFCcorr = [];
	case 'HOPF+HETEROGENOUS',
		% Here look at the heterogenous 
		folder = [fmri_dataset,'/results/HOPF+HETEROGENOUS'];
		[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_heterogenous_bif(sc_matrix,time_series,G,folder);
		% Will have to instead have something more uniform, but for now keep this structure.
	case 'DECO+WANG+BALANCED',
		% Here it is the balanced EI model which is essentially the noise degree model
		%  -- This is esentiall the noise-degree model but for now keep it seperate
		folder = [fmri_dataset,'/results/DECO+WANG+BALANCED'];
		FCD = [];ts_simulated_all = [];grandFCcorr = [];
	case 'HOPF+ANEC',
		% Here it is the HOPF model with ANEC 
		% Will have to work through this as well. 
		% This here works now -- I think however, will have to implement this version of the FCD as well and add it all together up to the final FCD.
		folder = [fmri_dataset,'/results/HOPF+ANEC'];
		[ts_simulated_all,grandFCcorr,FCD] = run_hopf_model_heterogenous_edge_bif(sc_matrix,time_series,G,folder);

	case 'NOISY+DEGREE',
		% Here is simply the noise and degree model
		folder = [fmri_dataset,'/results/NOISY+DEGREE'];
		% Maybe here run it then save it appropriately.
		[ts_simulated_all,FCcorr,grandFCcorr,FCD] =run_noisy_degree_model(sc_matrix,time_series,G,folder);
end


if(~isdir(folder))
	system(['mkdir -p ',folder]);
end

save([folder,'/',preproType,'.mat'],'G','ts_simulated_all','grandFCcorr','FCD');

% There is a slight issue -- FCD can be calculated in different ways -- either using the phases (as Gustavo is using now) or by standard FCD measures, will have to work out later.