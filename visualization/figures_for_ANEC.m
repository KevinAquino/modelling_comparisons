% This here is to show the figures for the ANEC under different Cross validation procedures.

data='UCLA';% Here we show the FC analyses for each of the models
Nsubjects=108; % Number of subjects in UCLA
% First load the FC matrices that were previously saved. (we dont use them its just to find variables)
% load('figures_ms/meanFC_UCLA.mat');
load('figures_ms/UCLA_data_params.mat')
simulation_params.G = linspace(0,4.5,20);
% Here to extract the inds for the upper triangle.
% inds=find(triu(ones(size(mFC(:,:,1))),1));
% Number of regions
% N=size(mFC,1);


% Color specification
tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% Where to save it all:
savePrefix='figures_ms/';

figure('color','white');
hold on;
load('/Users/kevinaquino/projects/modelling_gustavo/figures_ms/all_FC_FD_data_ANEC_DICER.mat');
plot(simulation_params.G,squeeze(all_FC_FD_data(:,1,:)),'Color',[0.9 0.9 0.9]);
plot(simulation_params.G,squeeze(all_FC_FD_data(:,2,:)),'Color',[0.7 0.7 0.7],'LineStyle','-.');
plot(simulation_params.G,mean(squeeze(all_FC_FD_data(:,1,:)),2),'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32);
plot(simulation_params.G,mean(squeeze(all_FC_FD_data(:,2,:)),2),'Color',theColors{3},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');

xlabel('G');
ylabel('$R_s,FCD$','Interpreter','LaTeX');
ylim([-0.1 1]);
xlim([0 4.5]);
set(gca,'fontSize',36)

saveeps(gcf,[savePrefix,'UCLA_ANEC_DiCER'],[20 20]);




figure('color','white');
hold on;
load('/Users/kevinaquino/projects/modelling_gustavo/figures_ms/all_FC_FD_data_ANEC_GSR.mat');
plot(simulation_params.G,squeeze(all_FC_FD_data(:,1,:)),'Color',[0.9 0.9 0.9]);
plot(simulation_params.G,squeeze(all_FC_FD_data(:,2,:)),'Color',[0.7 0.7 0.7],'LineStyle','-.');
plot(simulation_params.G,mean(squeeze(all_FC_FD_data(:,1,:)),2),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32);
plot(simulation_params.G,mean(squeeze(all_FC_FD_data(:,2,:)),2),'Color',theColors{2},'LineWidth',3,'Marker','.','MarkerSize',32,'LineStyle','-.');

xlabel('G');
ylabel('$R_s,FCD$','Interpreter','LaTeX');
ylim([-0.1 1]);
xlim([0 4.5]);
set(gca,'fontSize',36)

saveeps(gcf,[savePrefix,'UCLA_ANEC_GSR'],[20 20]);