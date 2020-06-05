
% Now something for FCD - maybe show non GSR and GSR in a big histogram
% load('/Users/kevinaquino/projects/modelling_gustavo/results/HOPF+GLOBAL/ICA-AROMA.mat')

% hold on;
for k=1:3,
	figure('color','white');
	[h,bins]=hist(FCD_data_all(:,k),50);
	bar(bins,h/sum(h),'barWidth',1,'FaceAlpha',0.5,'EdgeColor','none');
	set(gca,'fontSize',18);
	xlabel('Time-lagged Phase Coherence');
	ylabel('Normalized count')
	ylim([0 0.2]);
end
% legend({'ICA-AROMA','ICA-AROMA+GSR','ICA-AROMA+DiCER'});





clear FC FC_mean
for k=1:3
	for j=1:100,
		FC(:,:,j,k) = corr(time_series(:,:,j,k));
	end;
	FC_mean(:,:,k) = nanmean(FC(:,:,:,k),3);
end


% Fitting corrs for BE
inds=find(triu(ones(size(FC_mean(:,:,1))),1));
for j=1:length(G),
	% Model vs ICA-AROMA data
	FC_model=corr(ts_simulated_all(:,:,j));
	FC_emp = FC_mean(:,:,1);
	sim(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));	

	% Model vs GSR data
	FC_emp = FC_mean(:,:,2);
	sim2(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	
	

	% Model GSR vs Data GSR
	sim_t=ts_simulated_all(:,:,j)';
	sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	FC_model=corr(sim_t');
	sim3(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));

	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer_gsm(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
end


[~,ind]=max(sim);
[~,indGSR]=max(sim3);
[~,indDiCER]=max(sim_dicer_gsm)

ind_all=[ind,indGSR,indDiCER];

for k=1:3,
	sim_t=ts_simulated_all(:,:,ind)';
	if(k>1)
		sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	end
	[phase_model pmat] = phase_fcd(sim_t',2);
	figure('color','white');
	[h,bins]=hist(phase_model,50);
	bar(bins,h/sum(h),'barWidth',1,'FaceAlpha',0.5,'EdgeColor','none');
	set(gca,'fontSize',18);
	xlabel('Time-lagged Phase Coherence');
	ylabel('Normalized count')
	ylim([0 0.2]);

	figure('color','white');
	imagesc(pmat)
	axis image;
	xlabel('$\tau_i$','Interpreter','LaTeX');
	set(gca,'fontSize',18);
	colorbar
	colormap(hot);
	caxis([0 1]);
end



load('/Users/kevinaquino/projects/modelling_gustavo/results/BEI/BEI_simulation_88.mat')
G = G(1:18);

all_FCDS=calc_fcd(G,ts_simulated_all,FCD_data_all);
figure('color','white');
plot(G,all_FCDS);
legend({'MODEL vs ICA-AROMA','MODEL vs ICA-AROMA+GSR','MODEL+GSR vs ICA-AROMA+GSR','MODEL vs ICA-AROMA+DiCER','MODEL + GSR vs ICA-AROMA+DiCER'});
set(gca,'fontSize',18);


all_fits = calc_fit_all_FC(FC_mean,ts_simulated_all,G);
figure;plot(G,(all_fits' + 1 - all_FCDS )/2);



for k=1:3
	for j=1:100,
		tt=time_series(:,[1:34,42:41+34],j,k);
		FC(:,:,j,k) = corr(tt);
	end;
	FC_mean(:,:,k) = nanmean(FC(:,:,:,k),3);
end


load('/Users/kevinaquino/projects/modelling_gustavo/results/BTF/ICA-AROMA.mat');
all_FCDS=calc_fcd(G,ts_simulated_all,FCD_data_all);
figure('color','white');
plot(G,all_FCDS);
legend({'MODEL vs ICA-AROMA','MODEL vs ICA-AROMA+GSR','MODEL+GSR vs ICA-AROMA+GSR','MODEL vs ICA-AROMA+DiCER','MODEL + GSR vs ICA-AROMA+DiCER'});
set(gca,'fontSize',18);

load('/Users/kevinaquino/projects/modelling_gustavo/results/HOPF+GLOBAL/ICA-AROMA.mat')
all_FCDS=calc_fcd(G,ts_simulated_all,FCD_data_all);
figure('color','white');
plot(G,all_FCDS);
legend({'MODEL vs ICA-AROMA','MODEL vs ICA-AROMA+GSR','MODEL+GSR vs ICA-AROMA+GSR','MODEL vs ICA-AROMA+DiCER','MODEL + GSR vs ICA-AROMA+DiCER'});
set(gca,'fontSize',18);


t1 = zscore(time_series(:,:,14,1));
gs=mean(t1');
D=sum(C);
for gg=1:length(G),
	for j=1:length(D),
		ts_simulated_all(:,j,gg) = G(gg)*gs*D(j) + 0.5*randn(1,length(gs));
	end
end
all_FCDS=calc_fcd(G,ts_simulated_all,FCD_data_all);
figure('color','white');
plot(G,all_FCDS);
legend({'MODEL vs ICA-AROMA','MODEL vs ICA-AROMA+GSR','MODEL+GSR vs ICA-AROMA+GSR','MODEL vs ICA-AROMA+DiCER','MODEL + GSR vs ICA-AROMA+DiCER'});
set(gca,'fontSize',18);



g_use=[1:4:20];
for j=1:5,
	
	subplot(2,6,j);
	sim_t=ts_simulated_all(:,:,g_use(j))';
	% FCD_1 = phase_fcd(sim_t,2);
	FCD_1 = fcd_calculator(sim_t,10,5);

	h=hist(FCD_1(find(triu(ones(size(FCD_1)),1))),50);
	[~,~,FCD1]=kstest2(FCD_d(find(triu(ones(size(FCD_d)),1))),FCD_1(find(triu(ones(size(FCD_1)),1))));
	bar(h/sum(h),'barWidth',1,'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none');
	title(['G=',num2str(G(g_use(j))),' FCD = ',num2str(round(FCD1*100)/100)]);


	sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	% FCD_2 = phase_fcd(sim_t,2);

	FCD_2 = fcd_calculator(sim_t,10,5);
	hold on
	subplot(2,6,j+6);
	h=hist(FCD_2(find(triu(ones(size(FCD_2)),1))),50);
	bar(h/sum(h),'barWidth',1,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none');
	[~,~,FCD2]=kstest2(FCD_d(find(triu(ones(size(FCD_d)),1))),FCD_2(find(triu(ones(size(FCD_2)),1))));
	title(['G=',num2str(G(g_use(j))),' FCD GSR = ',num2str(round(FCD2*100)/100)]);
end
subplot(2,6,[6 12]);
% FCD_d = phase_fcd(t1',2);
hist(FCD_d(find(triu(ones(size(FCD_d)),1))),50);
title('Data');