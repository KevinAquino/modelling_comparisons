% A function here to look at the data for UCLA

% Calculate FCD for UCLA
time_series=empirical_params.time_series;

% Calculation of FCD"
for k=1:3,
	FCD_d = [];
	for ns=1:size(time_series,3),
		t1 = zscore(time_series(:,:,ns,k));
		% FCD_d = [FCD_d,fcd_calculator(t1',10,5)];
		FCD_d = [FCD_d phase_fcd(t1,2)];
	end
	FCD_data_all(:,k)= FCD_d;
end

save('figures_ms/FCD_data.mat','FCD_data_all');

% Calculation of FC means:
for k=1:3
	for j=1:size(time_series,3),
		tt=time_series(:,:,j,k);
		FC(:,:,j,k) = corr(tt);
	end;
	FC_mean(:,:,k) = nanmean(FC(:,:,:,k),3);
end


% Now to look at FCD and FC for Hopf-GSR
load('UCLA/results/HOPF+HETEROGENOUS+GSR/simulation.mat');
G=simulation_params.G;

all_FCDS=calc_fcd(G,ts_simulated,FCD_data_all);
all_fits = calc_fit_all_FC(FC_mean,ts_simulated,G);




load('UCLA/results/HOPF+HETEROGENOUS+DICER/simulation.mat');
G=simulation_params.G;

all_FCDS=calc_fcd(G,ts_simulated,FCD_data_all);
all_fits = calc_fit_all_FC(FC_mean,ts_simulated,G);


figure('color','white');
title('DiCER Optimization')
subplot(1,3,1);
plot(G,all_fits(4,:),'k.-','MarkerSize',15);
hold on
plot(G,all_FCDS(:,4),'r.-','MarkerSize',15);
xlabel('G');
ylabel('FCD/FC')
set(gca,'fontSize',18);

IND=2;
subplot(1,3,2);
imagesc(ts_simulated(:,:,IND,1)');
subplot(1,3,3);
for j=1:size(ts_simulated,4),
	FC_sim(:,:,j) = corr(ts_simulated(:,:,IND,j));
end
FC_sim= mean(FC_sim,3);
Combo=tril(FC_mean(:,:,3),-1) + triu(FC_sim,1);
nice_aparc_plotter(Combo,0.3*[-1 1]);
colormap(redwhitebluemap(100));