function [ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,time_series,G,folder)

	% Set up the structure for the Hopf model
	for sub=1:size(time_series,3),
		subject_ts{sub} = time_series(:,:,sub);
		corrFC(:,:,sub) = corr(squeeze(time_series(:,:,sub)).');		
	end

	meancorrFC = mean(corrFC,3);
	C = sc_matrix;
	wG = G;	
	TR = 2;
	Tmax = 147;%Check whats your tmax.
	subjects = size(time_series,3);
	N = size(sc_matrix,1);

	phfcddata = [];
	% Calculation of the total FCD
	TR = 2;
	for subject=1:10,
		ts = time_series(:,:,subject);
	phfcddata = [phfcddata,phase_fcd(ts,TR)];	
	end

	% ====== Calculation of the peak freq. for each node (should be made into a function if it is used multiple times) ==============

	% Filtering parameters (to discover the peak of the power spectrum)
	% This sets up the parameters to perform a band-pass filter of the data.

	f_diff = calculate_peak_frequency(time_series,TR);
	
	% ============================== Parameters for the Hopf model ==============================
	% 
	% Setting the global Hopf parameter
	a=-0.01*ones(N,2);	
	% With this, set up the omega vector:
	omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
	% Setting up the dt, and sig - strength of noise
	dt = 0.1;
	sig = 0.04;

	% Now here we can run the Global model for each global coupling vector
	upperTriangle = find(triu(ones(size(sc_matrix)),1));

	for coupling_index = 1:length(G),
		% Set up the weighted coupling matrix:
		wC = G(coupling_index)*C;		
		xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);		
		% Now after this, we can look at the fitting corr matrix.
		ts_simulated_all(:,:,coupling_index) = xs;
		corrFC_sim = corr(xs);
		grandFCcorr(coupling_index) = corr(corrFC_sim(upperTriangle),meancorrFC(upperTriangle));
		% Calculat ethe instantaneous phase from the dynamic FC measure
		phfcd = phase_fcd(xs,TR);

		% Now calculate the FCD
		[~,~,FCD(coupling_index)] = kstest2(phfcd,phfcddata);
	end

	% This is default setting the bifurcation parameter as a;
	bifpar = a;


end