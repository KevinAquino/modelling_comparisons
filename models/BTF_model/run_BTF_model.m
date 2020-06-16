% This here runs accesses the BTF execution commands
%

function [ts_simulated_all] = run_BTF_model(empirical_params,simulation_params)

	% This is the parameter within one segment of the BTF run - its run in 2s chunks
	SIM_time=2;
	% Extra runs to remove transients, i.e. run another 4*TR times, and remove the first 4TRs worth
	EXTRA_RUNS=4;
	% Connectivity parameters
	G=simulation_params.G;
	C=empirical_params.sc_matrix;
	% Here work out how long the simulation needs to be (in seconds)
	total_time=empirical_params.TR*simulation_params.N_FRAMES;
	% First need to work out if the balanced parameters  have been calculated
	N_iterations=ceil(simulation_params.N_FRAMES*empirical_params.TR/SIM_time)+EXTRA_RUNS;

	% Run the BTF model here
	[Vin,BOLD] = BTF_model(C,G,N_iterations);

	% Will 