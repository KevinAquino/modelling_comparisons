data='UCLA';% Here we show the FC analyses for each of the models

% First load the FC matrices that were previously saved.
load('figures_ms/meanFC_UCLA.mat');
% LOAD FCD of data
load('figures_ms/FCD_data.mat');


% Here to extract the inds for the upper triangle.
% inds=find(triu(ones(size(mFC(:,:,1))),1));
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
MODELS={'HOPF+HETEROGENOUS+GSR','HOPF+HETEROGENOUS+DICER','HOPF+HETEROGENOUS+AROMA'};

for nm=1:length(MODELS),
	MODEL=MODELS{nm};
	% Load the model simulation results:
	load([data,'/results/',MODEL,'/simulation.mat']);
	bif_1(:,nm) = simulation_params.bifpar(:,2);
	if(nm==3)
		bif_1(:,nm) = simulation_params.bifpar(:,3);
	end
end
figure;
[~,inds,total_order,~] = nice_aparc_plotter(ones(82,82),[-0.3 0.3],'black');

figure('color','white');
subplot(311)
bar(bif_1(total_order,3),'FaceColor',theColors{1},'BarWidth',1);%'Marker','.','MarkerSize',32);
for j=1:7,
	line([inds{j,1}(end)+0.5 inds{j,1}(end)+0.5],[min(bif_1(:,3)) max(bif_1(:,3))],'LineWidth',2,'Color','black')
	line([inds{j,2}(end)+0.5 inds{j,2}(end)+0.5],[min(bif_1(:,3)) max(bif_1(:,3))],'LineWidth',2,'Color','black')
end
axis tight
set(gca,'XTick',[],'FontSize',18,'box','on');ylabel('$a_i$','Interpreter','LaTeX')
subplot(312)
hold on;
bar(bif_1(total_order,1),'FaceColor',theColors{2},'BarWidth',1);%'Marker','.','MarkerSize',32);
for j=1:7,
	line([inds{j,1}(end)+0.5 inds{j,1}(end)+0.5],[min(bif_1(:,1)) max(bif_1(:,1))],'LineWidth',2,'Color','black')
	line([inds{j,2}(end)+0.5 inds{j,2}(end)+0.5],[min(bif_1(:,1)) max(bif_1(:,1))],'LineWidth',2,'Color','black')
end
axis tight
set(gca,'XTick',[],'FontSize',18,'box','on');ylabel('$a_i$','Interpreter','LaTeX')
subplot(313)
hold on;
bar(bif_1(total_order,2),'FaceColor',theColors{3},'BarWidth',1);%'Marker','.','MarkerSize',32);
for j=1:7,
	line([inds{j,1}(end)+0.5 inds{j,1}(end)+0.5],[min(bif_1(:,2)) max(bif_1(:,2))],'LineWidth',2,'Color','black')
	line([inds{j,2}(end)+0.5 inds{j,2}(end)+0.5],[min(bif_1(:,2)) max(bif_1(:,2))],'LineWidth',2,'Color','black')
end
axis tight
set(gca,'XTick',[],'FontSize',18,'box','on');ylabel('$a_i$','Interpreter','LaTeX')


saveeps(gcf,[savePrefix,'UCLA_Bifurfaction'],[20 20]);