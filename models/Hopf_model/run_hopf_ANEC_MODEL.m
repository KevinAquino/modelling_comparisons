% Code here to run the heterogenous model of the hopf model.
% Note we run this only on GSR data.

function [ts_simulated_all,simulation_params] = run_hopf_ANEC_MODEL(empirical_params,simulation_params)


	% Special parameters in the simulation_params for this model
	% simulation_params.UseGSR -- whether or not to use GSR on the mode
	% simulation_params.preproMethod -- whether or not to use 

	% Connectivity parameters
	G=simulation_params.G;
	C=empirical_params.sc_matrix;
	% Here work out how long the simulation needs to be (in seconds)
	total_time=empirical_params.TR*simulation_params.N_FRAMES;
	
	% For the hopf model we need to calculate the frequency peaks, this is not very variable, however we will calculate it here:
	% Note we calculate different frequencies dependending on what the prepro method is.
	all_subjects = empirical_params.time_series;
	for j=1:size(all_subjects,3)
		time_series=all_subjects(:,:,j,simulation_params.preproMethod);
		f_diff(:,j) = calculate_peak_frequency(squeeze(time_series)',empirical_params.TR);
	end
	% After this, just use the mean frequency for each node
	f_peak=mean(f_diff,2);

	% Store this for future use
	simulation_params.f_peak=f_peak;
	
	% Here now is the Hopf ANEC optimization, this optimization can take a while (not too long though!)
	[~,Coptim]=hopf_optimization_ANEC(C,G,f_peak,empirical_params.time_series(:,:,:,simulation_params.preproMethod),empirical_params.TR,total_time,simulation_params.UseGSR);
	simulation_params.ANEC_C = Coptim;
	% Note that this model is a homogenous Hopf bifurcation.

	% After these have been calculated (takes a long time) the simulation is made for multiple runs
	for g_ind=1:length(G),		
		for nr=1:simulation_params.N_RUNS,
			disp(['Run ',num2str(nr),'/',num2str(simulation_params.N_RUNS),'...']);
			% Running the default Hopf homogenous model now with the ANEC derived matrix.
			xs = run_hopf_no_ts(squeeze(Coptim(g_ind,:,:)),G(g_ind),f_peak,empirical_params.TR,total_time);		
			if(simulation_params.UseGSR)
				% Here using the parameters we got and do GSR as the data was fit with GSR
				xs = RegressNoiseSignal(xs.',mean(xs.')).';
			end
			ts_simulated_all(:,:,g_ind,nr) = xs;
		end 
	end

	% disp(['Finished Hopf homogenous simulation.']);
end