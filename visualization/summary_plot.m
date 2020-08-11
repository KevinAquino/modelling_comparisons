function summary_plot(ts_simulated_all,G,show_g,ve_data,GS,savePrefix,simulation_iteration)
		
	if(nargin<7)
		simulation_iteration=1;	
	end

	% Calculate VE1 for all simulations
	for ns=1:size(ts_simulated_all,4),
		for j=1:length(G),
			t1=zscore(ts_simulated_all(:,:,j,ns))';
			if(GS);t1=gsr_model(t1);end % GSR MODEL			
			[~,~,~,~,expl] = pca(t1');
			VE1_all(j,ns) = expl(1);
		end
	end
	

	VE1=mean(VE1_all,2);
	SEM_V1=std(VE1_all,[],2)/sqrt(size(ts_simulated_all,4));
	STD_V1=std(VE1_all,[],2);	
	% Plots of the carpet time series.
	figure('color','white')	
	for j=1:5,
		subplot(1,5,j);
		t1=ts_simulated_all(:,:,show_g(j),simulation_iteration)';
		cor_all=corr(t1',mean(t1)');
		[~,inds]=sort(cor_all);

		% GS if needed:
		% keyboard
		if(GS);t1=gsr_model(t1);end % GSR MODEL
		imagesc(zscore(t1(inds,:)')');
		axis off
	end
	colormap gray;
	savePng(gcf,[savePrefix,'_carpets'],[60 5]);


	% Plots here of the current FC matrices (not averaged)
	figure('color','white')

	for j=1:5,
		subplot(1,5,j);
		t1=ts_simulated_all(:,:,show_g(j),simulation_iteration)';
		% GSR IF NEEDED
		if(GS);t1=gsr_model(t1);end % GSR MODEL
		nice_aparc_plotter(corr(t1'),[-1 1],'black',0.5);
		colormap(redwhitebluemap)
		caxis([-1 1])
		% axis off
	end
	% colormap(cmap)

	savePng(gcf,[savePrefix,'_FCs'],[60 20]);

	cols={'--',':','-.'};

	figure('color','white')	
	plot(G,(VE1),'k.','lineWidth',1,'MarkerSize',32);	
	hold on
	
	for k=1:3
		line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
	end
	xlabel('Global Coupling (G)','fontSize',18);
	ylabel('VE1','fontSize',18)
	set(gca,'fontSize',18)
	legend({'Model','AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');	
	h=errorbar(G,VE1,STD_V1);	
	set(h,'LineWidth',2,'Color','black');
	xlim([0 4.5])
	saveeps(gcf,[savePrefix,'_VE1'],[20 20]);
end

function t1=gsr_model(t1)
	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);	
	% t1=t1'-(pinv(mean(t1'))'*(t1))'*mean(t1');
	% t1=t1';
end