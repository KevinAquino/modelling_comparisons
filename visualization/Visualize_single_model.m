% Here we look at visualizing a single model

% Load up the initial time series
load('UCLA_time_series_four_groups.mat');

% Now look at the individual empirical time series:
clear FC mean_FC
for j=1:3,
	for subject=1:110,
		FC(:,:,subject) = corr(time_series(:,:,subject,j));
	end
	mean_FC(:,:,j) = nanmean(FC,3);
end

% Define colormap:
cmap = [flipud(BF_getcmap('blues',9,0));[1,1,1],;BF_getcmap('reds',9,0)].';

for j=1:3,
	figure('color','white');	
	[square_mat,inds,total_order,all_regions] = nice_aparc_plotter(mean_FC(:,:,j),[-0.3 0.3],'black');
	if(j==1)
		caxis([-0.7 0.7]);
	end
	colormap(cmap');
end
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% Now look at the BEI
load('results/BEI/BEI_simulation_88.mat');
% Maybe look at the regions
gs=mean(ts_simulated_all(:,:,13),2);
tss = ts_simulated_all(:,:,13);
tss = tss - ((pinv(gs)*(tss))'*gs')';
% FC_sim=corr(ts_simulated_all(:,:,13));
FC_sim=corr(tss);
figure;
nice_aparc_plotter(FC_sim,[-1 1],'black');
colormap(cmap');
for j=1:3,
	FC_emp=mean_FC(:,:,j);
	for rg=1:7,
		region_inds=cell2mat(inds(rg,:));
		emp_fcs=FC_emp(region_inds,:);
		sim_fcs=FC_sim(region_inds,:);
		cor_region(rg,j) = corr(sim_fcs(:),emp_fcs(:));
	end
end
figure('color','white');
bar(cor_region,'BarWidth',1);
set(gca,'XTick',[1:7],'XTickLabel',all_regions);
ylabel('Corr of Network');
set(gca,'fontSize',18);
legend(noiseOptions);
colormap(cell2mat(theColors(1:3)))

C=ADJ_average;
C=C/max(C(:))*0.2;
figure;
nice_aparc_plotter(log10(C),[-7 -2],'white');

