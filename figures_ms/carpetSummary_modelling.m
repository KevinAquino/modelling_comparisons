function fh=carpetSummary_modelling(subjects,numSubs,functypeCell,funcLabelCell,folder)

fh = figure('color','white');
[ha, pos] = tight_subplot(2+length(functypeCell), numSubs, [0.01 0.01]);

ha = reshape(ha,numSubs,2+length(functypeCell)).';
pos = reshape(pos,numSubs,2+length(functypeCell)).';
cols={[.64 .08 .18],[0 .45 .74],[.93 .69 .13],[0.47 .67 .19],[.65 .65 .65]};
lastPlot=[];
% functypeCell={'','',''};

for subID=1:numSubs,
	subject=subjects{subID};
	% The mask nifti
	mask=MRIread([folder,subject,'_bold_space-MNI152NLin2009cAsym_dtissue_masked.nii.gz']);
	% The ordering for GS
	% gsOrdering = MRIread([subject,'_bold_space-MNI152NLin2009cAsym_gm_mask_gsordering_tissue.nii.gz']);
	% The clustered order plots	
	% sub-10274_task-rest
	% keyboard
	clustOrder = MRIread([folder,subject,'_task-rest_bold_space-MNI152NLin2009cAsym_variant-AROMAnonaggr_preproc+2P_detrended_hpf_clusterorder.nii']);
	% The estimates form data
	conf = dlmread([folder,subject,'_task-rest_bold_confounds.tsv'],'',2,0);
	gm 		= find(mask.vol==4);	
	
	% Load up the fMRI time series
	for nf=1:length(functypeCell),
		func=MRIread([folder,subject,functypeCell{nf}]);
		ts 	= zscore(reshape(func.vol,prod(func.volsize),func.nframes),[],2);			
		
		newOrder_gm = fromPythonOrderingToMatlab(gm,1+clustOrder.vol(gm),mask.vol);
		tsCLUST = ts([gm(newOrder_gm)],:);
		if(nf==1)
			GMsignal=zscore(mean(ts(gm,:)));
		end

		axes(ha(nf+2,subID));
		imagesc([tsCLUST]);
		caxis(1.2*[-1 1]);		
		colormap gray
		set(ha(nf+2,subID),'YTick',[])
		if(subID==1);		
			setLabelPosition(funcLabelCell{nf});
		end
	end

	axes(ha(1,subID));	
	plot(GMsignal,'Color',cols{1},'lineWidth',2);axis tight;ylim([-2 2]);
	if(subID==1);
		% ylabel('FD','FontSize',18);
		setLabelPosition('GM')
		title(subject,'fontSize',18)	
	else
		set(ha(4,subID),'YTickLabel',[])		
	end
	set(ha(4,subID),'XTickLabel',[],'fontSize',18)


	axes(ha(2,subID));
	plot(conf(:,7),'Color',cols{3},'lineWidth',2);axis tight;ylim([0 0.5]);
	if(subID==1);
		% ylabel('FD','FontSize',18);
		setLabelPosition('FD')
	else
		set(ha(4,subID),'YTickLabel',[])
	end
	set(ha(4,subID),'XTickLabel',[],'fontSize',18)

	drawnow
	

end



% keyboard

end

function setLabelPosition(label)
		yl = ylabel(label,'FontSize',18,'FontWeight','bold');
		position = get(yl,'Position');
		set(yl,'Position',[-10 position(2) -1])
end