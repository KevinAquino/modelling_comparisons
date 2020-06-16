% This here runs accesses the BTF execution commands
%

function [ts_simulated_all] = run_BTF_model(empirical_params,simulation_params)
	
	% This is the parameter within one segment of the BTF run - its run in 2s chunks
	SIM_time=2;
	% Extra iterations in the BTF code to remove transients, i.e. run another 4*TR times, and remove the first 4TRs worth
	EXTRA_ITERATIONS=4;
	% Connectivity parameters	
	% NOTE: In the bdtoolkit and in the BTF there is an extra normaliztion, so we do it here explicitly to remove confusion
	C=empirical_params.sc_matrix;
	C=C/max(C(:));
	C = C - diag(diag(C));
	G=simulation_params.G;
	% Also we need a translation factor that is used in the model since there is a normalization in the BTF model
	BTF_normFactor = empirical_params.normFactor;

	TR=empirical_params.TR;
	% Here work out how long the simulation needs to be (in seconds)
	total_time=empirical_params.TR*simulation_params.N_FRAMES;
	% First need to work out if the balanced parameters  have been calculated
	N_iterations=ceil(simulation_params.N_FRAMES*empirical_params.TR/SIM_time)+EXTRA_ITERATIONS;
	% Testing parameter 
	if(simulation_params.test)
		N_iterations=4;
	end

	if(simulation_params.batch.on)
		% Run the BTF for 1 complete run
		G_index = simulation_params.batch.g_ind;
		ts_simulated_all = BTF_model(C,G(G_index)*BTF_normFactor,N_iterations,TR);
	else
		% If there is no batching then will have to run this for all global coupling values. this will take a very long time ~40 days for 100 subjects X 20 G values, 300 seconds.
		for gind=1:length(G),
			for nr=1:simulation_params.N_RUNS,
				ts_simulated_all(:,:,g_ind,nr) = BTF_model(C,G(gind)*BTF_normFactor,N_iterations,TR);
			end
		end
	end

	

	