function run_hopf_model_heterogenous_bif(sc_matrix,time_series,G,folder)

	% Set up the structure for the Hopf model
	for sub=1:size(time_series,3),
		subject_ts{sub} = time_series(:,:,sub);
	end

	C = sc_matrix;
	wG = G;
	ldata = 2;%Number of subjects.
	Tmax = 147;%Check whats your tmax.

	Cfg.opt_a.nIters = 200;%Number of iterations for initial optimization
	Cfg.opt_a.updateStrength = 0.1;%Update strength for a optimization
	Cfg.opt_a.abortCrit = 0.1; %maximally 2.5% error
	Cfg.opt_a.gref = 3;%This is reference G, for the initial optimization only. Ignore if not using newOpt

	%NOTE: totalH is the time series for all subjects (nodes x timepoints)
	%inside a 1 x numberOfSubjects array. sc is the C matrix.


	[FC_simul, FC_emp, fitting, meta, ksP, PhasesD, Phases, bifpar]=...
	    hbif_NewOpt(C, subject_ts, Tmax, wG, ldata, Cfg);