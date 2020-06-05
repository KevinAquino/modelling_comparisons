function summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,GS)
		

	for j=1:length(G),
		t1=zscore(ts_simulated_all(:,:,j))';
		if(GS);t1=gsr_model(t1);end % GSR MODEL
		cor_all=corr(t1');
		GFC(j) = mean(cor_all(:));
		[~,~,~,~,expl] = pca(t1');
		VE1(j) = expl(1);
	end

	% First plots
	figure('color','white')	
	for j=1:5,
		subplot(2,5,j);
		t1=ts_simulated_all(:,:,show_g(j))';
		cor_all=corr(t1',mean(t1)');
		[~,inds]=sort(cor_all);

		% GS if needed:
		% keyboard
		if(GS);t1=gsr_model(t1);end % GSR MODEL
		imagesc(zscore(t1(inds,:)')');

		axis off
		caxis(4*[-1 1]);
		colormap gray;
		subplot(2,5,j+5);

		imagesc(corr(t1'));
		axis image
		caxis([0 1]);
	end
	linearVector = linspace(0,1,100).';
	reverseLinearVector = linearVector(end:-1:1);
	cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
	cmapB = [linearVector,linearVector,ones(size(linearVector))];
	cmap = [cmapB;cmapA];

	% CMAT PLOTs plots
	figure('color','white')

	for j=1:5,
		subplot(1,5,j);
		t1=ts_simulated_all(:,:,show_g(j))';
		% GSR IF NEEDED
		if(GS);t1=gsr_model(t1);end % GSR MODEL
		imagesc(corr(t1'));
		axis image
		caxis([-1 1])
		axis off
	end
	colormap(cmap)


	cols={'--',':','-.'};

	figure('color','white')
	subplot(211)
	plot(G,(VE1),'k','lineWidth',3)
	hold on
	for k=1:3
		line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
	end
	xlabel('Global Coupling (G)','fontSize',18);
	ylabel('VE1','fontSize',18)
	set(gca,'fontSize',18)
	legend({'Model','AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
	subplot(212)
	plot(G,(GFC),'k','lineWidth',3)
	hold on
	for k=1:3
		line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
	end
	legend({'Model','AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
	xlabel('Global Coupling (G)','fontSize',18);
	ylabel('GFC','fontSize',18)
	set(gca,'fontSize',18)
end
function t1=gsr_model(t1)
	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);	
	% t1=t1'-(pinv(mean(t1'))'*(t1))'*mean(t1');
	% t1=t1';
end