% HOPF Model without data:
% 
% SC Matrix normalization
function [G_all,metric] = coarse_grained_Hopf(C,FC_empirical,comparison)
	% First step here: normalization: ( to get point of reference )
	C = C/max(C(:))*0.2;
	N=size(C,1);

	% Frequency specification for each node (its relatively constant for all time series)
	f_diff=ones(N,1)*0.05;
	TR=0.72;

	% Parameters to tune the number of sims for the coarse and the fine resolutions
	N_coarse_sims=30;
	N_refined_sims=100;

	% Parameter chosen for how far to go in the initial sweep
	MAX_COARSE_G_SWEEP=6;

	% COARSE SIMULATIONS
	% ==================
	% coarse parameter tuning to first get an overal idea on which G-values we should use, 
	% we are looking to find G such we do not have NaNs
	

	G_all=linspace(0,MAX_COARSE_G_SWEEP,N_coarse_sims);
	cor_vals=run_hopf_with_g_values(C,G_all,f_diff,TR);
	cor_xformed=reshape(cor_vals,[numel(C),length(G_all)]);

	figure('color','white');
	title('Coarse resolution for simulations')	
	plot(G_all,mean(cor_xformed));xlabel('G');ylabel('<FC>');set(gca,'fontSize',14)

	

	% FINE SIMULATIONS
	% ==================
	% fine parameter simulations, here we aim to use the max G_value that we dont get divergence
	% this will change PER network so look for the NaNs of the FC, then grab the highest G value

	non_nan_G_values=~isnan(mean(cor_xformed));
	max_G_value=max(G_all(non_nan_G_values));

	% Now refine the search we did above:
	G_all=linspace(0,max_G_value,N_refined_sims);
	cor_vals=run_hopf_with_g_values(C,G_all,f_diff,TR);
	cor_xformed=reshape(cor_vals,[numel(C),length(G_all)]);

	figure('color','white');
	title('Fine resolution for simulations')	
	plot(G_all,mean(cor_xformed));xlabel('G');ylabel('<FC>');set(gca,'fontSize',14)	

	% Now comes the comparisons
	% 1. Get a FC summary measure
	% Run through this metric - either similarity or EVs?
	switch comparison
		case 'eigenvalues'
			eig_data = eig(FC_empirical);
			for j=1:N_refined_sims,
				simulation_eigs=eig(cor_vals(:,:,j));
				metric(j) = corr(simulation_eigs,eig_data);
			end
		otherwise
			disp('Comparison metric not found!');
	end
		
	figure('color','white');plot(G_all,metric);xlabel('G');ylabel(['Metric: ',comparison]);set(gca,'fontSize',14)
end

function cor_vals = run_hopf_with_g_values(C,G_all,f_diff,TR)
	for Gn = 1:length(G_all);
		G = G_all(Gn);
		% Frequency of each node:
		xs = run_hopf_no_ts(C,G,f_diff,TR);
		cor_vals(:,:,Gn)=corr(xs);	
	end
end