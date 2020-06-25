% Here we show the FCD analyses for each of the models

% First load the FC matrices that were previously saved.

load('figures_ms/FCD_data.mat');


% Color specification
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% Where to save it all:
savePrefix='figures_ms/';

% Now look for every model.
MODELS={'DECO+WANG+BALANCED','HOPF+GLOBAL','NOISY+DEGREE','BTF'};

for nm=1:length(MODELS),
	MODEL=MODELS{nm};
	% Load the model simulation results:
	load([data,'/results/',MODEL,'/simulation.mat']);

	all_FCDS=calc_fcd(simulation_params.G,ts_simulated,FCD_data_all);	
	

	figure('Color','white');
	hold on;
	plot(simulation_params.G,all_FCDS(:,1),'Color',theColors{1},'LineWidth',3,'Marker','.','MarkerSize',32);
	plot(simulation_params.G,all_FCDS(:,2),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle',':');
	plot(simulation_params.G,all_FCDS(:,3),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32);
	plot(simulation_params.G,all_FCDS(:,4),'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32);
	set(gca,'fontSize',18);
	xlabel('G');
	ylabel('FCD','Interpreter','LaTeX');
	xlim([0 4.5]);
	ylim([0 1]);
	saveeps(gcf,[savePrefix,'UCLA_',MODEL,'_FCD'],[20 20]);
	% keyboard
end