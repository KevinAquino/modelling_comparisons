function run_multiple_HBTF(cm)
	disp('USING..........')
	disp(cm)
	% Load in the empirical SC matrix
	load('empirical_data/exemplarSC.mat')
	addpath(genpath('~/projects/bdtoolkit')); % Add the bdtoolkit
	addpath(genpath(pwd))

	% Do some post-processing of the SCM:
	C = C/max(C(:));
	C = C - diag(diag(C));

	% After doing this run the BTF model for parameter

	% A really simple way to do this..
	[Vin,BOLD] = run_BTF_model(C,cm);
	save([num2str(cm),'.mat'],'Vin','BOLD');

end