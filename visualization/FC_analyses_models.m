data='UCLA';% Here we show the FC analyses for each of the models

% First load the FC matrices that were previously saved.

load('figures_ms/meanFC_UCLA.mat');

% Here to extract the inds for the upper triangle.
inds=find(triu(ones(size(mFC(:,:,1))),1));
% Number of regions
N=size(mFC,1);

% Taking the empirical matrix and making it flat.
FC_emp_flat=reshape(mFC,N*N,size(mFC,3));

% Color specification
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% Where to save it all:
savePrefix='figures_ms/';

% Now look for every model.
MODELS={'DECO+WANG+BALANCED','HOPF+GLOBAL','NOISY+DEGREE','BTF'};%,'HOPF+HETEROGENOUS+GSR','HOPF+HETEROGENOUS+DICER'};
% MODELS={'HOPF+HETEROGENOUS+GSR','HOPF+HETEROGENOUS+DICER'};
% MODELS={'HOPF+ANEC'};

for nm=1:length(MODELS),
	MODEL=MODELS{nm};
	% Load the model simulation results:
	load([data,'/results/',MODEL,'/simulation.mat']);

	% Now grab the models
	% Look at MODEL FC
	for nrun=1:size(ts_simulated,4)
		for g_ind=1:length(simulation_params.G),
			FC_model_all(:,:,g_ind,nrun)=corr(ts_simulated(:,:,g_ind,nrun));
		end
	end
	% Look at MODEL+GSR FC
	for nrun=1:size(ts_simulated,4)
		for g_ind=1:length(simulation_params.G),
			xs=ts_simulated(:,:,g_ind,nrun);
			xs=RegressNoiseSignal(xs',mean(xs'))';
			FC_model_GSR_all(:,:,g_ind,nrun)=corr(xs);
		end
	end	
	FC_model = mean(FC_model_all,4);
	FC_model_GSR = mean(FC_model_GSR_all,4);

	FC_model_flat=reshape(FC_model,N*N,length(simulation_params.G));
	FC_model_GSR_flat=reshape(FC_model_GSR,N*N,length(simulation_params.G));

	% Now have a look for every G value - different versions
	for g_ind=1:length(simulation_params.G),

		whole_brain_corr(g_ind,1)=corr(atanh(FC_emp_flat(inds,1)),atanh(FC_model_flat(inds,g_ind)));
		whole_brain_corr(g_ind,2)=corr(atanh(FC_emp_flat(inds,2)),atanh(FC_model_flat(inds,g_ind)));
		whole_brain_corr(g_ind,3)=corr(atanh(FC_emp_flat(inds,2)),atanh(FC_model_GSR_flat(inds,g_ind)));
		whole_brain_corr(g_ind,4)=corr(atanh(FC_emp_flat(inds,3)),atanh(FC_model_flat(inds,g_ind)));
	end

	figure('Color','white');
	hold on;
	plot(simulation_params.G,whole_brain_corr(:,1),'Color',theColors{1},'LineWidth',3,'Marker','.','MarkerSize',32);
	plot(simulation_params.G,whole_brain_corr(:,2),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle',':');
	plot(simulation_params.G,whole_brain_corr(:,3),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32);
	plot(simulation_params.G,whole_brain_corr(:,4),'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32);
	set(gca,'fontSize',18);
	xlabel('G');
	ylabel('$R_s$','Interpreter','LaTeX');
	xlim([0 4.5]);
	ylim([-.1 .7]);
	saveeps(gcf,[savePrefix,'UCLA_',MODEL,'_FC'],[20 20]);
	% keyboard
	% Here now show the actual FC matrices as a column vector?
	switch MODEL
	case 'DECO+WANG+BALANCED'
		ind_matrices=[11 14 9];
	case 'HOPF+GLOBAL'
		ind_matrices=[5 10 3];
	case 'NOISY+DEGREE'
		ind_matrices=[6 20 3];
	case 'BTF'
		ind_matrices=[15 14 10];
	end

	figure('Color','white');
	for j=[1,3],
		subplot(3,1,j)
		nice_aparc_plotter(FC_model(:,:,ind_matrices(j)),0.7*[-1 1],'black',1);
		if(j==3)
			caxis(0.3*[-1 1]);
		end
	end
	subplot(3,1,2)
	[~,inds_reg,total_order,all_regions] = nice_aparc_plotter(FC_model_GSR(:,:,ind_matrices(2)),0.3*[-1 1],'black',1);
	colormap(redwhitebluemap(100));
	% keyboard

	savePng(gcf,[savePrefix,'UCLA_',MODEL,'_FC_MATS'],[20 60]);

	% % Also show the region by region one.
	% for j=1:3,
	% FC_emp=mFC(:,:,j);
	% FC_mod=FC_model(:,:,ind_matrices(j));
	% if(j==2)
	% 	FC_mod=FC_model_GSR(:,:,ind_matrices(2));
	% 	% FC_emp=mFC(:,:,2);
	% end
	% 	for rg=1:7,
	% 		region_inds=cell2mat(inds_reg(rg,:));
	% 		emp_fcs=FC_emp(region_inds,:);
	% 		sim_fcs=FC_mod(region_inds,:);
	% 		cor_region(rg,j) = corr(sim_fcs(:),emp_fcs(:));
	% 	end
	% end
	% figure('color','white');
	% bar(cor_region,'BarWidth',1);
	% set(gca,'XTick',[1:7],'XTickLabel',all_regions);
	% ylabel('Corr of Network');
	% set(gca,'fontSize',18);
	% legend(noiseOptions);
	% colormap(cell2mat(theColors(1:3)))
end