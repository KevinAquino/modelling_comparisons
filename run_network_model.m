% Function here simulates the model with the given sc matrix and then provides some outputs that are then saved
% in a file that is saved for different types


% Outputs that have to be saved per G (would be good to keep it quite standard)
% The FC matrix 
% The FCD matrix
% The estimated time series as well (which is needed)

function run_model(sc_matrix,time_series,G,model)
% Now when this is solved for each element of G 
% it is then save seperately, this is to make sure 
% it is consisent once we go to the cluster level. 


% now the first thing one has to do is to check the model

switch model	
	case 'HOPF+GLOBAL',
		% Here it is the hopf global model i.e. setting a=0.
		
	case 'DECO+WANG',
		% Here look at the model what does this use? not the balanced EI model (will have to look)
	case 'BTF',
		% Here will have to grab the BTF model and then look at stuff
	case 'HOPF+HETEROGENOUS',
		% Here look at the heterogenous 
	case 'DECO+WANG+BALANCED',
		% Here it is the balanced EI model which is essentially the noise degree model
	case 'HOPF+ANEC',
		% Here it is the HOPF model with ANEC
	case 'NOISY+DEGREE',
		% Here is simply the noise and degree model
		folder = 'results/NOISY+DEGREE';
		% Maybe here run it then save it appropriately.
		[ts_simulated_all,FCcorr,grandFCcorr] =run_noisy_degree_model(sc_matrix,time_series,G,folder);
end