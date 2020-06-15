% This here sets up the modelling for the BEI model
function ts_simulated_all = run_BEI_MODEL(empirical_params,simulation_params)


	% Connectivity parameters
	G=simulation_params.G;
	C=empirical_params.sc_matrix;
	% Here work out how long the simulation needs to be (in seconds)
	total_time=empirical_params.TR*simulation_params.N_FRAMES;
	% First need to work out if the balanced parameters  have been calculated

	precalculated_Ji_file=['models/BEI_model/',empirical_params.sc_matrix_name,'_precalculated_Ji.mat'];
	% First find if there any previous calculations
	if exist([pwd,'/',precalculated_Ji_file], 'file') == 2
    	load(precalculated_Ji_file);
    	if(size(J_precalculated,1)~=length(G))
    		disp('Your Ji file is of different length - please rename it to avoid losses!');    		
    		return
    	end
	else
	     J_precalculated = precalculate_Balanced_weights(G,C);
	     save(precalculated_Ji_file,'J_precalculated','G');
	end

	% After these have been calculated (takes a long time) the simulation is made for multiple runs
	for g_ind=1:length(G),		
		for nr=1:simulation_params.N_RUNS,
			disp(['Run ',num2str(nr),'/',num2str(simulation_params.N_RUNS),'...']);
			ts_simulated_all(:,:,g_ind,nr) = balanced_EI_model(empirical_params.sc_matrix,G(g_ind),J_precalculated(g_ind,:)',total_time,empirical_params.TR);
		end 
	end
	disp(['Finished BEI simulation.']);
end