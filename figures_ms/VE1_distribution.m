tempColors = num2cell([255,105,97;97,168,255;178,223,138;117,112,179;255,179,71;255,231,152;]./255,2);  
col_vec=[1 4 2 3 5 6];
theColors = tempColors(col_vec);

% load('/Users/aquino/projects/modelling_comparisons/figures_ms/CNP_eps_08.mat')
load('~/Documents/fMRIClusterCorrect/stats/CNP_eps_08.mat');
% load('/Users/aquino/projects/modelling_comparisons/figures_ms/CNP_eps_08.mat')
subjects_restricted=metadata.ParticipantID;

% Now do the subject pruning:
subjects_restricted=subjects_restricted(setdiff(1:length(subjects_restricted),[5 7]));
% remaining_subs=1:length()
load('UCLA_time_series_four_groups.mat')
% load('/Users/aquino/projects/modelling_comparisons/empirical_data/UCLA_time_series_four_groups.mat')
load('/Users/kevinaquino/projects/modelling_gustavo/empirical_data/UCLA_time_series_four_groups.mat')
subject_ids_restricted=find(ismember(metadata.participant_id,subjects_restricted));
C=ADJ_average;C=C + C';C=C/max(C(:))*0.2;
% Have to leave out another subject too

% subject_ids_restricted=setdiff(subject_ids_restricted,[5 7]);

time_series = time_series(:,:,subject_ids_restricted,:);


[nt,N,~,~]=size(time_series);

for j=1:length(noiseOptions),
	for n=1:size(subject_ids_restricted),
		[coef,score,~,~,explained] = pca(zscore(time_series(:,:,n,j)).');			
		first_pc(n,j) = explained(1);
	end

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
figure;
rain_test;
xlabel('VE by 1st PC (%)')
set(gcf,'Color','white');


% QCFC and mean FC for the three methods.

% First flatten and get UT
uT=find(triu(ones(N,N),1));
for j=1:length(noiseOptions)
	for n=1:size(subject_ids_restricted),
		FC_all(:,:,n,j) = corr(time_series(:,:,n,j));
	end
end

for j=1:length(noiseOptions),
% Get corr as a flat matrix
	[flattened_cor]= reshape(squeeze(FC_all(:,:,:,j)),[N*N,length(subject_ids_restricted)]);
	% disp(['Found: ',num2str(sum(isnan(flattened_cor(:)))/numel(flattened_cor)*100),'% of NaNs ']);
	% nonNans=find(~isnan(flattened_cor));	
	% Do QC-FC
	[QC_FC(:,j) p(:,j)] = corr(flattened_cor',metadata.FD(subject_ids_restricted));
	data2{j} = abs(QC_FC(:,j));
	prop_j(j) = sum(p(uT,j)<0.01)

end
figure;
rain_test;
xlabel('QC-FC');




% Now do some mean FC as well.
linearVector = linspace(0,1,100).';
reverseLinearVector = linearVector(end:-1:1);
cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];cmapB = [linearVector,linearVector,ones(size(linearVector))];
cmap = [cmapB;cmapA];
for j=1:length(noiseOptions),
	mFC=mean(FC_all(:,:,:,j),3);
	figure('color','white');
	[square_mat,inds,total_order,all_regions] = nice_aparc_plotter(mFC,[-0.3 0.3],'black');
	colormap(cmap);
	if(j==1)
		caxis([-0.7 0.7]);
	end
end



figure('color','white');
imagesc(log10(C))
[square_mat,inds,total_order,all_regions] = nice_aparc_plotter(log10(C),[-4 -1],'white');


for n=1:7,
	reg_ind=cell2mat(inds(n,1:2));
	miniMat_C=C(:,total_order(reg_ind));
	miniMat_FC=mFC(:,total_order(reg_ind));
	figure;plot(log10(miniMat_C(:)),miniMat_FC(:),'.')
end

