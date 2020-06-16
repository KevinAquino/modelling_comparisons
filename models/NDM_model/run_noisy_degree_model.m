% This function runs the simple noisy degree model
function ts_simulated_all = run_noisy_degree_model(empirical_params,simulation_params)



	% Connectivity parameters
	G=simulation_params.G;
	C=empirical_params.sc_matrix;

	% This runs the model using the same global signal for each time point - it need not be that way, it can just be an autocorrelated
	% version of random noise (just matching the autocorrelation)

	% Work out the degree
	D=sum(C);

	% Calculating the NDM for each subject for each G
	all_subjects = empirical_params.time_series;
	for g_ind=1:length(G),
		for subject=1:size(all_subjects,3)
			time_series=all_subjects(:,:,subject,1);
			gs=mean(time_series,2);
			for j=1:length(D),
				ts_simulated_all(:,j,g_ind,subject) = G(g_ind)*gs*D(j) + 0.5*randn(length(gs),1);
			end
		end

	end