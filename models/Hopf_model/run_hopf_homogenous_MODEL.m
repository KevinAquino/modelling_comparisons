% Code here to run the homogenous model of the hopf model.
function ts_simulated_all = run_hopf_homogenous_MODEL(empirical_params,simulation_params)


	% Connectivity parameters
	G=simulation_params.G;
	C=empirical_params.sc_matrix;
	% Here work out how long the simulation needs to be (in seconds)
	total_time=empirical_params.TR*simulation_params.N_FRAMES;
	
	% For the hopf model we need to calculate the frequency peaks, this is not very variable, however we will calculate it here:
	all_subjects = empirical_params.time_series;
	for j=1:size(all_subjects,3)
		time_series=all_subjects(:,:,j,1);
		f_diff(:,j) = calculate_peak_frequency(squeeze(time_series)',empirical_params.TR);
	end
	% After this, just use the mean frequency for each node
	f_peak=mean(f_diff,2);

	% Store this for future use
	simulation_params.f_peak=f_peak;
	
	% After these have been calculated (takes a long time) the simulation is made for multiple runs
	N=size(C,1);
	a=0*ones(N,2);
	for g_ind=1:length(G),		
		for nr=1:simulation_params.N_RUNS,
			disp(['Run ',num2str(nr),'/',num2str(simulation_params.N_RUNS),'...']);			
			ts_simulated_all(:,:,g_ind,nr) = run_hopf_no_ts(C,G(g_ind),f_peak,empirical_params.TR,total_time,a);
		end 
	end
	disp(['Finished Hopf homogenous simulation.']);
end