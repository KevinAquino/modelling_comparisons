figFolder = '/Users/kevinaquino/Documents/Journal_publications/InPrep/ClusterCorrectfMRIPaper/figs';
mFDVE = figure('color','white');
files={'CNP_eps_08.mat','BZ_eps_08.mat','cambridge_me.mat'}
datasites={'UCLA','BZ','CME'};

	% just loading here for odering for all pipelines
	id_labels=tdfread('Gordon/CommunityModified.txt');
	id_all = id_labels.DA;
	id_all = ['DA ';id_all];
	all_types = unique(id_all,'rows');
	for j=1:size(all_types,1),
		all_types_cell{j} = all_types(j,:);
	end
	% Can change the order if nesc.
	ids_cell = cell(13,1);
	for j=1:size(id_all,1)
		id_find = find(strcmp(id_all(j,:),all_types_cell));
		ids_cell{id_find,1} = [ids_cell{id_find,1},j];
	end
	% For plotting:
	all_idz=cell2mat(ids_cell.');
	for j=1:length(ids_cell),
		nums(j) = numel(ids_cell{j});
	end
	line_breaks=cumsum(nums);
	line_breaks = line_breaks(1:end-1);
	pos_labels = cumsum(nums) - nums/2;



for opendataset=1:3,		
	disp(['DATASET: ',datasites{opendataset},'.............']);
	load(files{opendataset})

	% load('CNP_data_2P_GSR_postAROMA.mat')
	% QC-FC distance dependence:
	N_methods = 3;
	numPrePro= 3;
	plotting=1;
	% allQC = allQC(1:4);
	allQC = allQC([1 2 4]);
	% allQC = allQC(1:N_methods);
	first_pc = first_pc(:,[1 2 4]);
	FC=FC(:,:,:,[1 2 4]);

	% noiseOptions = noiseOptions(1:4);
	% noiseOptions = noiseOptions(1:N_methods);
	noiseOptions = noiseOptions([1 2 4]);
	noiseOptions{3} = 'AROMA+2P+DiCER';
	FSize=18;
	fh = figure('color','white');
	hold on;
	tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
	col_vec=[1 4 2 3 5 6];
	theColors = tempColors(col_vec);

	for j=1:N_methods,
		% subplot(ceil(N_methods/2),2,j)
		subplot(1,3,j)
		BF_PlotQuantiles(ROIDistVec(allQC(j).NaNFilter),allQC(j).QCFC,11,0,0,tempColors{col_vec(j)})
		hold on;
		plot([0:200],zeros(1,201),'--','Color','k','LineWidth',2)
		xlabel('distance (mm)')
		ylabel('QC-FC')
		set(gca,'fontSize',16)
		title(noiseOptions{j},'fontSize',18,'Interpreter','none')
		% text('FontSize',24,'BoxOpacity',0,'Font','Arial Bold')
		ylim([-0.1 0.3])
	end


	if(plotting)
		saveeps(fh,[figFolder,'/QC_FC_metrics',datasites{opendataset}],[30 10]);
	end

	% Creating a bar chart for QC

	x = noiseOptions;
	% y = num2cell(100 - [allData.PercExcluded]);
	xy = cell(x);
	for i = 1:length(x)
		% if y{i} ~= 100
			% xy{i} = strcat(x{i});
		% elseif y{i} == 100
			xy{i} = x{i};
	end
	% end



	Fig_QCFC_Dist = figure('color','w', 'units', 'centimeters', 'pos', [0 0 16 9], 'name',['Fig_QCFC_Dist']); box('on'); movegui(Fig_QCFC_Dist,'center');
	sp = subplot(2,3,[2:3]);
	if(opendataset~=2)
		text(0,-25,'A','fontSize',36,'fontWeight','bold')
	else
		text(0,-25,'C','fontSize',36,'fontWeight','bold')
	end
	pos = get(sp,'Position');
	% set(gca,'Position',[pos(1)*2.75, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
	% set(gca,'Position',[pos(1)*3.25, pos(2)*1.2, pos(3)*1.2, pos(4)*1]); % [left bottom width height]

	% Create data
	data = {[allQC(:).QCFC_PropSig_unc]'};
	data_std = cell(1,length(data)); [data_std{:}] = deal(zeros(size(data{1})));

	% Create table
	T = table(data{1},'RowNames',noiseOptions,'VariableNames',{'QCFC_PropSig'})

	% Create bar chart
	clear extraParams
	extraParams.xTickLabels = xy;
	extraParams.xLabel = ''; % 'Pipeline'
	% extraParams.yLabel = 'QC-FC (%)';
	extraParams.yLabel = 'QC-FC uncorrected (%)';
	extraParams.theColors = theColors;
	extraParams.theLines = {'-','-','-','-','-','-'};
	extraParams.yLimits = [0 35];
	extraParams.FSize = 18;

	TheBarChart(data,data_std,false,extraParams)
	set(gca,'fontSize',18)
	ax = gca;

	% Set axis stuff
	ax.FontSize = FSize;
	% ax.Interpreter = 'none'
	% ------------------------------------------------------------------------------
	% QCFC distributions
	% ------------------------------------------------------------------------------
	sp = subplot(2,3,[4:6]);
	data2 = {allQC(:).QCFC};
	rain_test;
	yy=get(gca,'Ylim');
	xlim([-0.6 0.6]);

	if(opendataset~=2)
		text(-.7,yy(2),'B','fontSize',36,'fontWeight','bold');	
	else
		text(-.7,yy(2),'D','fontSize',36,'fontWeight','bold');
	end
	

	% pos = get(sp,'Position');
	% set(gca,'Position',[pos(1)*1, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
	% data2 = {allQC(:).QCFC};
	% clear extraParams
	% % extraParams.theLabels = {allData(:).noiseOptionsNames};
	% extraParams.customSpot = '';
	% extraParams.add0Line = true;
	% extraParams.theColors = theColors;
	% BF_JitteredParallelScatter(data2,1,1,0,extraParams);
	% ax = gca;

	% % Set axis stuff
	% ax.FontSize = FSize;
	% % ax.XTick = [1:size(data{1},1)];
	% ax.XTick = [];
	% ax.XTickLabel = [];
	% ax.XLim = ([0 numPrePro+1]);

	% 	ax.YLim = ([-0.6 1]);
	% % end

	xlabel('QC-FC (Pearson''s r)')

	% % add text
	% TextRotation = 0;
	% strprec = '%0.2f';
	% data3 = {allQC(:).QCFC_AbsMed};

	% text(1:size(data3,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data3,2)),num2str([data3{1,:}]',strprec),... 
	% 'HorizontalAlignment','right',... 
	% 'VerticalAlignment','middle',...
	% 'Color','black',...
	% 'FontSize', FSize,...
	% 'Rotation',TextRotation)

	% keyboard
	% view(90,90)

	if(plotting)
		% saveeps(Fig_QCFC_Dist,[figFolder,'/QC_FC_metrics_uncorr',datasites{opendataset}],[30 40]);
		savePng(Fig_QCFC_Dist,[figFolder,'/QC_FC_metrics_uncorr',datasites{opendataset}],[30 30]);
	end


	lbs{1,1} = {'A'};lbs{1,2} = {'D'};
	lbs{2,1} = {'B'};lbs{2,2} = {'E'};
	lbs{3,1} = {'C'};lbs{3,2} = {'F'};

	% Take average of FC matrices
	m2 = squeeze(mean(tanh(FC),3));
	fc_mat = figure('color','white');
	for j=1:N_methods,
		subplot(2,N_methods,j)
		m3 = m2(:,:,j);
		% imagesc(m3(all_idz,all_idz));		
		imagesc(m2(:,:,j));
		axis image
		if(j==1)
			caxis([-0.3,0.3]);
		else
			caxis([-0.1,0.1]);
		end


		set(gca,'YDir','normal','fontSize',18);
		lb_col=[0.5 0.5 0.5];
		text(-100,380,lbs{j,1},'fontSize',36,'fontWeight','bold')
		% for g_id=1:length(line_breaks)
		% 	line([1 333],[line_breaks(g_id) line_breaks(g_id)],'Color',lb_col,'LineWidth',1);
		% 	line([line_breaks(g_id) line_breaks(g_id)],[1 333],'Color',lb_col,'LineWidth',1);
		% end
		% for g_id=1:length(pos_labels)
		% 	if(g_id==9)
		% 		spacer = -15;
		% 	else
		% 		spacer=0;
		% 	end
		% 	text(-20 - spacer,pos_labels(g_id),deblank(all_types_cell{g_id}),'HorizontalAlignment','right','fontSize',18);
		% 	text(pos_labels(g_id),-20 - spacer,deblank(all_types_cell{g_id}),'HorizontalAlignment','right','fontSize',18,'Rotation',90);
		% end

		title(noiseOptions{j},'fontSize',18,'Interpreter','none')
		colormap([flipud(BF_getcmap('blues',9,0));[1,1,1],;BF_getcmap('reds',9,0)])
		p=colorbar;
		pp=get(p,'Limits');
		set(p,'Ticks',[pp(1) 0 pp(2)]);
		axis off	
		subplot(2,N_methods,j+N_methods);
		fc_nm = m2(:,:,j);
		mat = find(triu(ones(333,333),1));
		bb=linspace(-0.5,0.5,50);

		raincloud_plot(fc_nm(mat),'color',theColors{j},'box_on',1);
		% disp(mean(fc_nm(mat)));
		% [histograms_fc,bins] = hist(fc_nm(mat),bb);
		% histograms_fc=histograms_fc/numel(fc_nm(mat));
		% p = bar(bins,histograms_fc,'BarWidth',1);
		% p.FaceColor = theColors{j};
		% p.LineWidth=2;
		xlabel('Pearson''s r');
		% ylabel('Density');
		set(gca,'fontSize',18,'box','on','YTick',[],'XTick',[-0.4 0 0.5 1])
		
		Ypos=get(gca,'YLim')
		% if(j==1)
			xlim([-0.4 1]);
		% else
		% end
		text(-0.7,Ypos(2),lbs{j,2},'fontSize',36,'fontWeight','bold')
		% keyboard

		% text(-.5,0.45,lbs{j,2},'fontSize',36,'fontWeight','bold')
		% xlim([-0.25 0.75])
		% ylim([0 0.4]);
	end

	% keyboard
	if(plotting)
		saveeps(fc_mat,[figFolder,'/FC_average_mat',datasites{opendataset}],[40 20]);
	end



	VE_violin = figure('color','white');
	% sp = subplot(1,2,1);
	% pos = get(sp,'Position');
	% set(gca,'Position',[pos(1)*1, pos(2)*1.2, pos(3)*1.4, pos(4)*1]); % [left bottom width height]
	for j=1:N_methods,
		data2{j} = first_pc(:,j);
	end
	% % data2 = {first_pc};
	% clear extraParams
	% % extraParams.theLabels = noiseOptions(1:4);
	% extraParams.theLabels = noiseOptions(1:3);
	% extraParams.customSpot = '';
	% extraParams.add0Line = true;
	% extraParams.theColors = theColors;
	% BF_JitteredParallelScatter(data2,1,1,0,extraParams);
	% ax = gca;

	% % Set axis stuff
	% ax.FontSize = FSize;
	% % ax.XTick = [1:size(data2,1)];
	% ax.XTick = [];
	% ax.XTickLabel = [];
	% ax.XLim = ([0 4]);

	% 	ax.YLim = ([0 50]);
	% % end
	rain_test;
	xlabel('VE by 1st PC (%)')

	% add text
	TextRotation = 0;
	strprec = '%0.2f';
	for j=1:N_methods,
		data3{j} = median(first_pc(:,j))
	end

	% Max for VE1 for AROMA+2P
	[vv,ii] = max(h1{2}.XData)
	sub_l=plot([h1{2}.XData(ii),h2{2}.XData(ii),h3{2}.XData(ii)],[h1{2}.YData(ii),h2{2}.YData(ii),h3{2}.YData(ii)],'k*','MarkerSize',10);
	legend([h1{1} h2{1} h3{1} sub_l], [noiseOptions,metadata.ParticipantID(ii)]);

	% Follow through each time

	% text(1:size(data3,2),repmat(ax.YLim(2) - ax.YLim(2)*.05,1,size(data3,2)),num2str([data3{1,:}]',strprec),... 
	% 'HorizontalAlignment','right',... 
	% 'VerticalAlignment','middle',...
	% 'Color','black',...
	% 'FontSize', FSize,...
	% 'Rotation',TextRotation)

	% view(90,90)
	% set(gca,'XLim',[0.5 3.5])
	% for j=1:N_methods
	% 	text(j,ax.YLim(2)+5,noiseOptions{j},'FontSize',FSize,'FontWeight','bold')
	% end
	figure;plot(cell2mat(data2));set(gca,'XTickLabel',metadata.ParticipantID,'XTickLabelRotation',90,'XTick',[1:length(data2{1})])

	if(plotting)
		% saveeps(VE_violin,[figFolder,'/VE_violin',datasites{opendataset}],[25 15])
		savePng(VE_violin,[figFolder,'/VE_violin',datasites{opendataset}],[25 15]);
	end


	% figure('color','white')
	% hold on;
	% for j=1:N_methods,
	% 	fc_nm = m2(:,:,j);
	% 	[histograms_fc,bins] = hist(fc_nm(:),1000);
	% 	p = stairs(bins,histograms_fc);
	% 	p.Color = theColors{j};
	% 	p.LineWidth=2;
	% 	xlabel('Corr');
	% 	ylabel('count')	
	% 	set(gca,'fontSize',18,'box','on')
	% 	xlim([-1 1])
	% end
	% legend(noiseOptions)



	% fdGSSE=figure('color','white');
	% plot(metadata.fdJenk_m,metadata.GSSD,'.','MarkerSize',18);
	% r = corr(metadata.fdJenk_m,metadata.GSSD);
	% text([max(metadata.fdJenk_m)/2],max(metadata.GSSD)*5/6,['R^2=',num2str(r)],'fontSize',24);
	% xlabel('mFD (mm)');ylabel('GSSD (arb)');
	% set(gca,'fontSize',24);
	% if(plotting)
	% 	saveeps(fdGSSE,[figFolder,'/fdGSSE'],[30 20]);
	

	% gssdVsVE=figure('color','white')
	% for j=1:N_methods,
	% 	subplot(ceil(N_methods/2),2,j)
	% 	r = corr(metadata.GSSD,first_pc(:,j));
	% 	plot(metadata.GSSD,first_pc(:,j),'.','MarkerSize',18,'Color',theColors{j});
	% 	text([max(metadata.GSSD)/2],35,['R^2=',num2str(r)],'fontSize',18);
	% 	title(noiseOptions{j},'fontSize',18,'Interpreter','none')
	% 	xlabel('GSSD');ylabel('1st PC VE (%)');
	% 	ylim([0 40])
	% 	set(gca,'fontSize',18);
	% end
	% if(plotting)
	% 	saveeps(gssdVsVE,[figFolder,'/gssdVsVE'],[30 20]);
	% end

	
	for nm=1:N_methods,
		figure(mFDVE);
		subplot(3,3,nm + (opendataset-1)*3)
		r = corr(metadata.fdJenk_m,first_pc(:,nm));
		plot(metadata.fdJenk_m,first_pc(:,nm),'.','MarkerSize',18,'Color',theColors{nm});
		text([max(metadata.fdJenk_m)/2],35,['R^2=',num2str(r)],'fontSize',18);
		title(noiseOptions{nm},'fontSize',18,'Interpreter','none')
		xlabel('mFD (mm)');ylabel('1st PC VE (%)');
		ylim([0 40])
		set(gca,'fontSize',18);
	end


	% keyboard



	% figure('color','white')
	% hold on;
	% for j=1:N_methods,
	% 	fc_nm = FC(:,:,:,j);
	% 	[histograms_fc,bins] = hist(fc_nm(:),1000);
	% 	p = stairs(bins,histograms_fc);
	% 	p.Color = theColors{j};
	% 	p.LineWidth=2;
	% 	xlabel('Corr');
	% 	ylabel('count')
	% 	legend(noiseOptions)
	% 	set(gca,'fontSize',18,'box','on')
	% 	xlim([-1 1])
	% end



	% meanFCvsVE=figure('color','white')
	% for j=1:N_methods,
	% 	mn_fc = (reshape(FC(:,:,:,j),[333*333,size(FC,3)]));

	% 	mn_fc(isinf(mn_fc)) = nan;
	% 	mn_fc =nanmean(abs(mn_fc));
	% 	subplot(ceil(N_methods/2),2,j)
	% 	r = corr(mn_fc.',first_pc(:,j));
	% 	plot(mn_fc,first_pc(:,j),'.','MarkerSize',18,'Color',theColors{j});
	% 	text([max(mn_fc)/2],35,['R^2=',num2str(r)],'fontSize',18);
	% 	title(noiseOptions{j},'fontSize',18,'Interpreter','none')
	% 	xlabel('<|FC|>');ylabel('1st PC VE (%)');
	% 	ylim([0 40])
	% 	set(gca,'fontSize',18);
	% end
	% if(plotting)
	% 	saveeps(meanFCvsVE,[figFolder,'/meanFCvsVE'],[40 30]);
	% end



	% figure;
	% fdGSSE=figure('color','white');
	% subplot(1,3,1)
	% % fdGSSE=figure('color','white');
	% plot(metadata.fdJenk_m,metadata.GSSD,'.','MarkerSize',18);
	% ylim([0 7]);xlim([0.05 0.3])
	% r = corr(metadata.fdJenk_m,metadata.GSSD);
	% text([max(metadata.fdJenk_m)/2],max(metadata.GSSD)*5/6,['R^2=',num2str(r)],'fontSize',24);
	% xlabel('mFD (mm)');ylabel('GSSD (arb)');
	% title('Preproc','fontSize',18,'Interpreter','none')
	% set(gca,'fontSize',24);

	% subplot(1,3,2)
	% plot(metadata.fdJenk_m,(metadata.GSSD_prepro(:,1)),'.','MarkerSize',18);
	% ylim([0 7]);xlim([0.05 0.3])
	% r = corr(metadata.fdJenk_m,(metadata.GSSD_prepro(:,1)));
	% text([max(metadata.fdJenk_m)/2],max(metadata.GSSD)*5/6,['R^2=',num2str(r)],'fontSize',24);
	% title('AROMA+2P','fontSize',18,'Interpreter','none')
	% xlabel('mFD (mm)');ylabel('GSSD (arb)');
	% set(gca,'fontSize',24);
	% if(plotting)
	% 	% saveeps(fdGSSE,[figFolder,'/fdGSSE'],[30 20]);
	% end
	% subplot(1,3,3)
	% plot(metadata.fdJenk_m,(metadata.GSSD_prepro(:,2)),'.','MarkerSize',18);
	% r = corr(metadata.fdJenk_m,(metadata.GSSD_prepro(:,2)));
	% text([max(metadata.fdJenk_m)/2],max(metadata.GSSD)*5/6,['R^2=',num2str(r)],'fontSize',24);
	% xlabel('mFD (mm)');ylabel('GSSD (arb)');
	% ylim([0 7]);xlim([0.05 0.3])
	% title('AROMA+2P+DBSCAN','fontSize',18,'Interpreter','none')
	% set(gca,'fontSize',24);
	% if(plotting)
	% 	saveeps(fdGSSE,[figFolder,'/fdGSSE_prepros'],[50 20]);
	% end

	% fc_nm_gmr = FC(:,:,:,2);
	% fc_nm_dbs = FC(:,:,:,4);

	% figure('color','white');plot(tanh(fc_nm_gmr(:)),tanh(fc_nm_dbs(:)),'k.');xlim([-1 1]);ylim([-1 1]);
	% hold on;line([-1,1],[-1,1],'Color','r')
	% hold on;line([-1,1],[0,0],'Color','g')
	% hold on;line([0,0],[-1,1],'Color','g')
	% % xlabel('GS r (z-transformed)');ylabel('DBSCAN r (z-transformed)');set(gca,'fontSize',18);
	% xlabel('GS r');ylabel('DBSCAN r');set(gca,'fontSize',18);
	% axis image



	% fc_nm_gmr = mean(tanh(FC(:,:,:,2)),3);
	% fc_nm_dbs = mean(tanh(FC(:,:,:,4)),3);

	% figure('color','white');plot((fc_nm_gmr(:)),(fc_nm_dbs(:)),'k.');xlim([-1 1]);ylim([-1 1]);
	% hold on;line([-1,1],[-1,1],'Color','r')
	% hold on;line([-1,1],[0,0],'Color','g')
	% hold on;line([0,0],[-1,1],'Color','g')
	% % xlabel('GS r (z-transformed)');ylabel('DBSCAN r (z-transformed)');set(gca,'fontSize',18);
	% xlabel('GS r');ylabel('DBSCAN r');set(gca,'fontSize',18);
	% axis image
end

	% if(plotting)
	% 	saveeps(mFDVE,[figFolder,'/meanFdVsVE'],[50 40]);
	% end