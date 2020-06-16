% This is for the batch runs
function batch_BTF_run(G_ind,Run)
	% Here will have the interface for the models (make getenvs)
	simulation_params.batch.on = 1;

	% Make getenvs for this
	simulation_params.batch.g_ind = str2num(G_ind);
	simulation_params.batch.run = str2num(Run);

	% Run this for the particular instance
	simulate_homogenous_models;

