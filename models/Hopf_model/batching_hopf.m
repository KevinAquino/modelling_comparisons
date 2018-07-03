% default HOPF model for each of these.

for tn=1:5,
	for j=1:24,
		subject_ts_all{j,tn} = squeeze(time_series(:,:,j,tn));
	end
end


for tn = 1:5,	
	disp(['RUNNING HOPF MODEL THROUGH DATA PREPROCESSED WITH ...............' analyses{tn}]);
	subject_ts = subject_ts_all(:,tn)
	if(tn==5)
		run_hbif_GSR;
	else
		run_hbif;
	end
	fits{tn} = fitting;
	fits_sim{tn} = FC_simul;
	fits_emp{tn} = FC_emp;
	fits_bifar{tn} = bifpar;
end

%

hh = ones(size(C));
for j=1:7,
	FC_emp = fits_emp{1};
	FC_sim = fits_sim{1}(:,:,j);
	inds = setdiff(find(triu(hh,1)),find(isnan(FC_emp(:))));
	fitting(j) = corr(FC_emp(inds),FC_sim(inds));
end

fits{1} = fitting;

figure('Color','white');
hold on
for tn=1:5,
	plot(wG,fits{tn},'.-','lineWidth',2,'MarkerSize',30);
end
xlabel('G');ylabel('Fitting Cor');
legend(analyses,'fontSize',18,'location','eastOutSide','Interpreter','none')
set(gca,'fontSize',18);



figure('Color','white');
for tn=1:5,
	subplot(2,5,tn)
	imagesc(fits_emp{tn}(:,:,end));
	title(['FCemp ' analyses{tn} ],'Interpreter','none','fontSize',18)
	axis image
	axis off;
	caxis([0 1]);
	subplot(2,5,tn+5)
	imagesc(fits_sim{tn}(:,:,end));
	title(['FCsimul ' analyses{tn} ],'Interpreter','none','fontSize',18)
	axis image
	axis off;
	caxis([0 1]);
end

figure('Color','white');
hold on
M = [];
for tn=1:5,
	h = plot(wG,fits{tn},'.-','lineWidth',2,'MarkerSize',30);
	[fitting] = noiseModelPlotting(fits_emp{tn},0);
	plot(wG,max(fitting)*ones(size(wG)),'--','lineWidth',2,'Color',h.Color);
	M{tn,1} = analyses{tn};
	M{tn,2} = ['noise model: ' analyses{tn}];
end
M = reshape(M.',numel(M),1);
xlabel('G');ylabel('Fitting Cor');
legend(M,'fontSize',18,'location','eastOutSide','Interpreter','none')
set(gca,'fontSize',18);

for j=1:5;bif_all(:,j) = fits_bifar{j}(end,:);end

figure('Color','white');
hold on;
for tn=1:5,
	plot(bif_all(:,tn),'lineWidth',3);
end
xlabel('region');ylabel('a');
legend(analyses,'fontSize',18,'location','eastOutSide','Interpreter','none')
set(gca,'fontSize',18);



D = sum(C,1);
all_plusd(:,1:5) = bif_all;
all_plusd(:,6) = D;
figure('Color','white');
imagesc(corr(all_plusd))
lbs = cell(6,1);lbs(1:5) = analyses;lbs{6} = 'Degree';
set(gca,'XTickLabel',lbs,'YTickLabel',lbs,'TickLabelInterpreter','none','fontSize',18,'XTickLabelRotation',45);
caxis([0 1]);
title('Correlation of a')
colorbar
axis image
colormap jet




for subject=1:24,
	for j=1:68,
		pows = abs(fft(time_series(j,:,subject,5)));
		pows = pows/sum(pows);
		ratio(j,subject) = sum(pows(1:37))/sum(pows(38:74));
	end
end
avg_ratio = mean(ratio,2);

clear lbs
D = sum(C,1);
all_plusd(:,1:5) = bif_all;
all_plusd(:,6) = D;
% all_plusd(:,7) = avg_ratio;
figure('Color','white');
imagesc(abs(corr(all_plusd)))
lbs = cell(6,1);lbs(1:5) = analyses;lbs{6} = 'Degree';%lbs{7} = 'low to high P';
set(gca,'XTickLabel',lbs,'YTickLabel',lbs,'TickLabelInterpreter','none','fontSize',18,'XTickLabelRotation',45);
caxis([0 1]);
title('Correlation of a')
colorbar
axis image
colormap jet