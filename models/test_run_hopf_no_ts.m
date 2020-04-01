% Structural connectivity test:

% First find the optimal parameters



load('empirical_data/SC_HCP_cortexOnly.mat')
% SC Matrix normalization
C=AdjDens;
C = C/max(C(:))*0.2;

% Frequency specification for each node (its realtively )
f_diff=ones(N,1)*0.05;
TR=0.72;

% Global coupling parameter (tune this)
G = 5;
% Frequency of each node:

xs = run_hopf_no_ts(C,G,f_diff,TR);



imagesc(corr(xs));

G_vec = linspace(0,10,20);
[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(C,ts_all,G_vec,'',TR);