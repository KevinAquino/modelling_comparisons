function run_hopf_model(sc_matrix,time_series,G,folder)

	for sub=1:10,
		subject_ts{sub} = time_series(:,:,sub);
	end
	C = sc_matrix;
	wG = G;
	ldata = 10;%Number of subjects.
	run_hbif;
	

	keyboard
	% Note:here instead we can send it off to the cluster