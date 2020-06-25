% This here is to show the figures for the ANEC under different Cross validation procedures.

data='UCLA';% Here we show the FC analyses for each of the models
Nsubjects=108; % Number of subjects in UCLA
% First load the FC matrices that were previously saved. (we dont use them its just to find variables)
load('figures_ms/meanFC_UCLA.mat');
load('figures_ms/UCLA_data_params.mat')

% Here to extract the inds for the upper triangle.
inds=find(triu(ones(size(mFC(:,:,1))),1));
% Number of regions
N=size(mFC,1);


% Color specification
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% Where to save it all:
savePrefix='figures_ms/';

figure;
hold on;
% This here is where we show ANEC over multiple iterations of the cross validation.
MODEL='HOPF+ANEC';
for iteration=1:10,
	disp(['Iteration:',num2str(iteration)]);
	clear FC_data
	% Here load it up the first iteration:
	load([data,'/results/',MODEL,'/simulation_DICER_',num2str(iteration),'.mat']);
	% Calculate Mean for FC data, and FCD for that data
	% First find left over subjects:
	test_subjects=setdiff([1:Nsubjects],simulation_params.rand_perms(iteration,:));
	% Now generate the mean for GSR'd data
	test_data=empirical_params.time_series(:,:,test_subjects,3);

	% Now calculate mean FC and FCD
	for n=1:size(test_data,3)
		FC_data(:,:,n) = corr(test_data(:,:,n));
	end
	FC_emp=mean(FC_data,3);
	% Taking the empirical matrix and making it flat.
	FC_emp_flat=reshape(FC_emp,N*N,1);
	ts_simulated=ts_simulated(:,:,:,1:size(test_data,3));
	% Calculate Mean for Simulation and FC for simulation - make sure its done with the same number of subjects
	for nrun=1:size(test_data,3)
		for g_ind=1:length(simulation_params.G),
			% xs=ts_simulated(:,:,g_ind,nrun);
			% xs=RegressNoiseSignal(xs',mean(xs'))';
			% FC_model_all(:,:,g_ind,nrun)=corr(xs);
			FC_model_all(:,:,g_ind,nrun)=corr(ts_simulated(:,:,g_ind,nrun));
		end
	end

	FC_model = mean(FC_model_all,4);
	FC_model_flat=reshape(FC_model,N*N,length(simulation_params.G));

	% Note data already has GSR
	% Now look at the data
	for g_ind=1:length(simulation_params.G)				
		whole_brain_corr(g_ind)=corr(atanh(FC_emp_flat(inds)),atanh(FC_model_flat(inds,g_ind)));		
	end
	plot(simulation_params.G,whole_brain_corr,'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32);

	% Now calculate the FCD
	% Calculation of FCD"
	FCD_d = [];
	for ns=1:size(test_data,3),
		t1 = zscore(test_data(:,:,ns));
		% FCD_d = [FCD_d,fcd_calculator(t1',10,5)];
		FCD_d = [FCD_d phase_fcd(t1,2)];
	end
	FCD_ANEC = calc_fcd_single(simulation_params.G,ts_simulated,FCD_d);
	plot(simulation_params.G,FCD_ANEC,'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');
	drawnow
end