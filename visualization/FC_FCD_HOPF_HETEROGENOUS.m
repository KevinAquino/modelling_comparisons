data='UCLA';% Here we show the FC analyses for each of the models

% First load the FC matrices that were previously saved.
load('figures_ms/meanFC_UCLA.mat');
% LOAD FCD of data
load('figures_ms/FCD_data.mat');


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
% MODELS={'HOPF+HETEROGENOUS+GSR','HOPF+HETEROGENOUS+DICER'};
MODELS={'HOPF+HETEROGENOUS+AROMA'}

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

	% Now also calculate FCDs
	all_FCDS=calc_fcd(simulation_params.G,ts_simulated,FCD_data_all);	

	% Now have a look for every G value - different versions
	figure('Color','white');
	hold on;
	switch MODEL
	case 'HOPF+HETEROGENOUS+GSR'	
		for g_ind=1:length(simulation_params.G)				
			whole_brain_corr(g_ind)=corr(atanh(FC_emp_flat(inds,2)),atanh(FC_model_GSR_flat(inds,g_ind)));		
		end
		plot(simulation_params.G,whole_brain_corr,'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32);
		plot(simulation_params.G,all_FCDS(:,3),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');
	case 'HOPF+HETEROGENOUS+DICER'	
		for g_ind=1:length(simulation_params.G)		
			whole_brain_corr(g_ind)=corr(atanh(FC_emp_flat(inds,3)),atanh(FC_model_flat(inds,g_ind)));		
		end
		plot(simulation_params.G,whole_brain_corr,'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32);
		plot(simulation_params.G,all_FCDS(:,4),'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');
	
	case 'HOPF+HETEROGENOUS+AROMA'	
		for g_ind=1:length(simulation_params.G)		
			whole_brain_corr(g_ind)=corr(atanh(FC_emp_flat(inds,1)),atanh(FC_model_flat(inds,g_ind)));		
		end
		plot(simulation_params.G,whole_brain_corr,'Color',theColors{1},'LineWidth',3,'Marker','.','MarkerSize',32);
		plot(simulation_params.G,all_FCDS(:,1),'Color',theColors{1},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');
	end
	


	set(gca,'fontSize',18);
	xlabel('G');
	ylabel('$R_s$,$FCD$','Interpreter','LaTeX');
	xlim([0 4.5]);
	ylim([0 1]);
	saveeps(gcf,[savePrefix,'UCLA_',MODEL,'_FC_FCD'],[20 20]);
	% % keyboard
	% % Here now show the actual FC matrices as a column vector?
	% switch MODEL
	% case 'DECO+WANG+BALANCED'
	% 	ind_matrices=[11 14 9];
	% case 'HOPF+GLOBAL'
	% 	ind_matrices=[5 10 3];
	% case 'NOISY+DEGREE'
	% 	ind_matrices=[6 20 3];
	% case 'BTF'
	% 	ind_matrices=[15 14 10];
	% end

	% figure('Color','white');
	% for j=[1,3],
	% 	subplot(3,1,j)
	% 	nice_aparc_plotter(FC_model(:,:,ind_matrices(j)),0.7*[-1 1],'black',1);
	% 	if(j==3)
	% 		caxis(0.3*[-1 1]);
	% 	end
	% end
	% subplot(3,1,2)
	% [~,inds_reg,total_order,all_regions] = nice_aparc_plotter(FC_model_GSR(:,:,ind_matrices(2)),0.3*[-1 1],'black',1);
	% colormap(redwhitebluemap(100));
	% % keyboard

	% savePng(gcf,[savePrefix,'UCLA_',MODEL,'_FC_MATS'],[20 60]);

end