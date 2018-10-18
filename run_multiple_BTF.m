% Load in the empirical SC matrix
load('empirical_data/exemplarSC.mat')
addpath(genpath('~/projects/bdtoolkit')); % Add the bdtoolkit
addpath(genpath(pwd))

% Do some post-processing of the SCM:
C = C/max(C(:));
C = C - diag(diag(C));

% After doing this run the BTF model for parameter


for cm=[0.1,0.2];
	% A really simple way to do this..
	[Vin,BOLD_interp] = run_BTF_model(C,cm);
	save([num2str(cm),'.mat'],'Vin','BOLD_interp');
end