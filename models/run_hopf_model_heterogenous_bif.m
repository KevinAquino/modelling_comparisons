function [ts_simulated_all,fitting,bifpar,FCD] = run_hopf_model_heterogenous_bif(sc_matrix,time_series,G,folder)


	[N,Tmax,subjects] = size(time_series);
	% Set up the structure for the Hopf model
	for sub=1:subjects,
		subject_ts{sub} = time_series(:,:,sub);
	end
	% keyboard
	C = sc_matrix;
	wG = G;
	ldata = subjects;%Number of subjects.	
	TR = 2;

	Cfg.opt_a.nIters = 200;%Number of iterations for initial optimization
	Cfg.opt_a.updateStrength = 0.1;%Update strength for a optimization
	Cfg.opt_a.abortCrit = 0.1; %maximally 2.5% error
	Cfg.opt_a.gref = 3;%This is reference G, for the initial optimization only. Ignore if not using newOpt

	%NOTE: totalH is the time series for all subjects (nodes x timepoints)
	%inside a 1 x numberOfSubjects array. sc is the C matrix.


	[FC_simul, FC_emp, fitting, meta, ksP, PhasesD, Phases, bifpar, ts_simulated_all]=...
	    hbif_NewOpt(C, subject_ts, Tmax, wG, ldata, Cfg);


    % Using the same way to calculate FCD as been done in all the other scripts
	phfcddata = [];
	for subject=1:subjects,
		ts = time_series(:,:,subject);
		phfcddata = [phfcddata,phase_fcd(ts,TR)];	
	end

	for coupling_index = 1:length(G)
		% This processing of the time series is TEMPRORARY will change when the code above is cleaned up!
		phfsim = phase_fcd(ts_simulated_all(1:10:end,:,coupling_index),TR);
		[~,~,FCD(coupling_index)] = kstest2(phfsim,phfcddata);
	end