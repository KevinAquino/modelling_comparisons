% Here this code will show all the carpet plots and associated FC matrices and VE1s for each model.
data='UCLA';

% Load in the VE1 data
load('figures_ms/VE1_data_UCLA.mat');
VE1_data=cell2mat(VE1_UCLA);
mean_VE1=mean(VE1_data);
% First the BEI model.
MODELS={'DECO+WANG+BALANCED','HOPF+GLOBAL','NOISY+DEGREE','BTF'};
for nm=1:length(MODELS),
	MODEL=MODELS{nm};
	% Load the model simulation results:
	load([data,'/results/',MODEL,'/simulation.mat']);
	% After this then load in the script.
	show_g=[1 6 11 16 20];
	if(strcmp(MODEL,'BTF'))
		show_g=[1 6 11 15 16];
		show_sim=4;
	else
		show_sim=1;
	end
	summary_plot(ts_simulated,simulation_params.G,show_g,mean_VE1,0,['figures_ms/',data,'_',MODEL,'_NO_GSR'],show_sim);
end


% Now GSR'd

MODELS={'DECO+WANG+BALANCED','HOPF+GLOBAL','BTF'};
for nm=1:length(MODELS),
	MODEL=MODELS{nm};
	% Load the model simulation results:
	load([data,'/results/',MODEL,'/simulation.mat']);
	% After this then load in the script.
	show_g=[1 6 11 16 20];
	if(strcmp(MODEL,'BTF'))
		show_g=[1 6 11 15 16];
	end
	summary_plot(ts_simulated,simulation_params.G,show_g,mean_VE1,1,['figures_ms/',data,'_',MODEL,'_GSR'],show_sim);
end