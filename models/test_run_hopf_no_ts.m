% HOPF Model without data:
% 
% SC Matrix normalization
C=AdjDens;
C = C/max(C(:))*0.2;

% Frequency specification for each node (its relatively constant for all time series)
f_diff=ones(N,1)*0.05;
TR=0.72;

% Global coupling parameter (tune this)
G = 5;
% Frequency of each node:
xs = run_hopf_no_ts(C,G,f_diff,TR);


% Calculating Hopf model for the given data, to work out 
G_vec = linspace(0,10,40);
[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(C,all_time_series,G_vec,'',TR);

figure('color','white');
plot(G_vec,grandFCcorr,'r.-','lineWidth',2,'MarkerSize',30);
hold on;
plot(G_vec,FCD,'k.-','lineWidth',2,'MarkerSize',30);
title('Looking for G to max. FCorr, min. FCD');
xlabel('G')
ylabel('FCD (black),FCorr (red)')
set(gca,'fontSize',18);