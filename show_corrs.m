% % figure('color','white')
% % show_g=[1:4:20];
% % for j=1:5,
% % 	subplot(2,5,j);
% % 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% % 	cor_all=corr(t1',mean(t1)');
% % 	[~,inds]=sort(cor_all);
% % 	imagesc(zscore(t1(inds,:)')');
% % 	axis off
% % 	caxis(4*[-1 1]);
% % 	colormap gray;
% % 	subplot(2,5,j+5);

% % 	imagesc(corr(t1'));
% % 	axis image
% % 	caxis([0 1]);
% % end


% % for j=1:length(G),
% % 	t1=zscore(ts_simulated_all(:,:,j))';
% % 	cor_all=corr(t1');
% % 	GFC(j) = mean(cor_all(:));
% % end


% % linearVector = linspace(0,1,100).';
% % reverseLinearVector = linearVector(end:-1:1);
% % cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
% % cmapB = [linearVector,linearVector,ones(size(linearVector))];
% % cmap = [cmapB;cmapA];

% % figure('color','white')
% % show_g=[1:4:20];
% % for j=1:5,
% % 	subplot(1,5,j);
% % 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% % 	imagesc(corr(t1'));
% % 	axis image
% % 	caxis([-1 1])
% % 	axis off
% % end
% % colormap(cmap)


% % for hh=1:20,
% % 	[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,squeeze(permute(time_series(:,:,1:10,1),[2 1 4 3])),G,'');
% % 	for j=1:length(G),
% % 		t1=zscore(ts_simulated_all(:,:,j))';
% % 		cor_all=corr(t1');
% % 		GFC(hh,j) = mean(cor_all(:));
% % 	end
% % end





% % for hh=1:20,
% % 	[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,squeeze(permute(time_series(:,:,1:10,1),[2 1 4 3])),G,'');
% % 	for j=1:length(G),
% % 		t1=zscore(ts_simulated_all(:,:,j))';
% % 		cor_all=corr(t1');
% % 		GFC(hh,j) = mean(cor_all(:));		
% % 		[~,~,~,~,expl] = pca(t1');
% % 		VE1(hh,j) = expl(1);
% % 	end
% % end





% % for k=1:3,
% % 	t1 = zscore(time_series(:,:,14,k));
% % 	[~,~,~,~,expl] = pca(t1);
% % 	co = corr(t1);
% % 	GFC_exp(k) = mean(abs(co(:)));
% % 	ve_data(k) = expl(1);
% % end


% % cols={'--',':','-.'};

% % figure('color','white')
% % subplot(211)
% % plot(G,mean(VE1),'k','lineWidth',3)
% % hold on
% % for k=1:3
% % 	line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% % end
% % xlabel('Global Coupling (G)','fontSize',18);
% % ylabel('VE1','fontSize',18)
% % set(gca,'fontSize',18)
% % legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% % subplot(212)
% % plot(G,mean(GFC),'k','lineWidth',3)
% % hold on
% % for k=1:3
% % 	line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% % end
% % legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% % xlabel('Global Coupling (G)','fontSize',18);
% % ylabel('GFC','fontSize',18)
% % set(gca,'fontSize',18)




	
% % for j=1:length(G),
% % 	t1=zscore(ts_simulated_all(:,:,j))';
% % 	cor_all=corr(t1');
% % 	GFC(j) = mean(cor_all(:));
% % 	[~,~,~,~,expl] = pca(t1');
% % 	VE1(j) = expl(1);
% % end


% show_g=[1 9 13 14 15];



% figure('color','white')
% % show_g=[1:4:20];
% for j=1:5,
% 	subplot(2,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	cor_all=corr(t1',mean(t1)');
% 	[~,inds]=sort(cor_all);
% 	imagesc(zscore(t1(inds,:)')');
% 	axis off
% 	caxis(4*[-1 1]);
% 	colormap gray;
% 	subplot(2,5,j+5);

% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([0 1]);
% end
% linearVector = linspace(0,1,100).';
% reverseLinearVector = linearVector(end:-1:1);
% cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
% cmapB = [linearVector,linearVector,ones(size(linearVector))];
% cmap = [cmapB;cmapA];

% figure('color','white')

% for j=1:5,
% 	subplot(1,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([-1 1])
% 	axis off
% end
% colormap(cmap)


% cols={'--',':','-.'};

% figure('color','white')
% subplot(211)
% plot(G,(VE1),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('VE1','fontSize',18)
% set(gca,'fontSize',18)
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% subplot(212)
% plot(G,(GFC),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('GFC','fontSize',18)
% set(gca,'fontSize',18)



% clear ts_simulated_all VE1;
% vals = {'G_0.05.mat','G_0.1.mat','G_0.15.mat','G_0.2.mat','G_0.25.mat','G_0.3.mat','G_0.35.mat','G_0.4.mat','G_0.45.mat','G_0.5.mat','G_0.55.mat','G_0.6.mat','G_0.65.mat','G_0.7.mat','G_0.75.mat','G_0.8.mat','G_0.85.mat','G_0.9.mat','G_0.95.mat','G_1.05.mat','G_1.1.mat','G_1.15.mat','G_1.2.mat','G_1.25.mat','G_1.3.mat','G_1.35.mat','G_1.4.mat','G_1.45.mat','G_1.5.mat','G_1.55.mat','G_1.6.mat','G_1.65.mat','G_1.7.mat','G_1.75.mat','G_1.8.mat','G_1.85.mat','G_1.9.mat','G_1.95.mat','G_1.mat','G_2.05.mat','G_2.1.mat','G_2.15.mat','G_2.2.mat','G_2.25.mat','G_2.3.mat','G_2.35.mat','G_2.4.mat','G_2.45.mat','G_2.5.mat','G_2.55.mat','G_2.6.mat','G_2.65.mat','G_2.7.mat','G_2.75.mat','G_2.8.mat','G_2.85.mat','G_2.9.mat','G_2.95.mat','G_2.mat','G_3.05.mat','G_3.1.mat','G_3.15.mat','G_3.2.mat','G_3.25.mat','G_3.3.mat','G_3.35.mat','G_3.4.mat','G_3.45.mat','G_3.5.mat','G_3.55.mat','G_3.6.mat','G_3.65.mat','G_3.7.mat','G_3.75.mat','G_3.8.mat','G_3.85.mat','G_3.9.mat','G_3.95.mat','G_3.mat','G_4.05.mat','G_4.1.mat','G_4.15.mat','G_4.2.mat','G_4.25.mat','G_4.3.mat','G_4.35.mat','G_4.4.mat','G_4.45.mat','G_4.5.mat','G_4.55.mat','G_4.6.mat','G_4.65.mat','G_4.7.mat','G_4.75.mat','G_4.8.mat','G_4.85.mat','G_4.9.mat','G_4.95.mat','G_4.mat','G_5.mat'};
% for j=1:100;
% 	G(j) = we;
% 	load(vals{j});
% 	ts_simulated_all(:,:,j) = bds;	
% 	grand_fc(j) = fitting;
% end

% ts_simulated_all = ts_simulated_all(:,:,1:5:end);
% G = G(1:5:end);
% grand_fc = grand_fc(1:5:end);


	
% for j=1:length(G),
% 	t1=zscore(ts_simulated_all(:,:,j))';
% 	cor_all=corr(t1');
% 	GFC(j) = mean(cor_all(:));
% 	[~,~,~,~,expl] = pca(t1');
% 	VE1(j) = expl(1);
% end




% figure('color','white')
% show_g=[1:4:20];
% for j=1:5,
% 	subplot(2,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	cor_all=corr(t1',mean(t1)');
% 	[~,inds]=sort(cor_all);
% 	imagesc(zscore(t1(inds,:)')');
% 	axis off
% 	caxis(4*[-1 1]);
% 	colormap gray;
% 	subplot(2,5,j+5);

% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([0 1]);
% end
% linearVector = linspace(0,1,100).';
% reverseLinearVector = linearVector(end:-1:1);
% cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
% cmapB = [linearVector,linearVector,ones(size(linearVector))];
% cmap = [cmapB;cmapA];

% figure('color','white')

% for j=1:5,
% 	subplot(1,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([-1 1])
% 	axis off
% end
% colormap(cmap)


% cols={'--',':','-.'};

% figure('color','white')
% subplot(211)
% plot(G,(VE1),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('VE1','fontSize',18)
% set(gca,'fontSize',18)
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% subplot(212)
% plot(G,(GFC),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('GFC','fontSize',18)
% set(gca,'fontSize',18)






% %% GSR


% figure('color','white')
% show_g=[1:4:20];
% for j=1:5,
% 	subplot(2,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';	
% 	cor_all=corr(t1',mean(t1)');
% 	[~,inds]=sort(cor_all);
% 	% GMR
% 	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);
% 	imagesc(zscore(t1(inds,:)')');
% 	axis off
% 	caxis(4*[-1 1]);
% 	colormap gray;
% 	subplot(2,5,j+5);

% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([0 1]);
% end


% for j=1:length(G),
% 	t1=zscore(ts_simulated_all(:,:,j))';
% 	cor_all=corr(t1');
% 	GFC(j) = mean(cor_all(:));
% end


% linearVector = linspace(0,1,100).';
% reverseLinearVector = linearVector(end:-1:1);
% cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
% cmapB = [linearVector,linearVector,ones(size(linearVector))];
% cmap = [cmapB;cmapA];

% figure('color','white')
% show_g=[1:4:20];
% for j=1:5,
% 	subplot(1,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([-1 1])
% 	axis off
% end
% colormap(cmap)


% for hh=1:20,
% 	[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,squeeze(permute(time_series(:,:,1:10,1),[2 1 4 3])),G,'');
% 	for j=1:length(G),
% 		t1=zscore(ts_simulated_all(:,:,j))';
% 		cor_all=corr(t1');
% 		GFC(hh,j) = mean(cor_all(:));
% 	end
% end





% for hh=1:20,
% 	[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,squeeze(permute(time_series(:,:,1:10,1),[2 1 4 3])),G,'');
% 	for j=1:length(G),
% 		t1=zscore(ts_simulated_all(:,:,j))';
% 		cor_all=corr(t1');
% 		GFC(hh,j) = mean(cor_all(:));		
% 		[~,~,~,~,expl] = pca(t1');
% 		VE1(hh,j) = expl(1);
% 	end
% end




% cols={'--',':','-.'};

% figure('color','white')
% subplot(211)
% plot(G,mean(VE1),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('VE1','fontSize',18)
% set(gca,'fontSize',18)
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% subplot(212)
% plot(G,mean(GFC),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('GFC','fontSize',18)
% set(gca,'fontSize',18)




	
% for j=1:length(G),
% 	t1=zscore(ts_simulated_all(:,:,j))';
% 	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);
% 	cor_all=corr(t1');
% 	GFC(j) = mean(abs(cor_all(:)));
% 	[~,~,~,~,expl] = pca(t1');
% 	VE1(j) = expl(1);
% end


% show_g=[1 9 13 14 15];



% figure('color','white')
% % show_g=[1:4:20];
% for j=1:5,
% 	subplot(2,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	cor_all=corr(t1',mean(t1)');
% 	[~,inds]=sort(cor_all);

% 	% GMR
% 	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);
% 	imagesc(zscore(t1(inds,:)')');
% 	axis off
% 	caxis(4*[-1 1]);
% 	colormap gray;
% 	subplot(2,5,j+5);

% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([0 1]);
% end
% linearVector = linspace(0,1,100).';
% reverseLinearVector = linearVector(end:-1:1);
% cmapA = [ones(size(linearVector)),reverseLinearVector,reverseLinearVector];
% cmapB = [linearVector,linearVector,ones(size(linearVector))];
% cmap = [cmapB;cmapA];

% figure('color','white')

% for j=1:5,
% 	subplot(1,5,j);
% 	t1=zscore(ts_simulated_all(:,:,show_g(j)))';
% 	% GMR
% 	t1=t1-(pinv(mean(t1))'*(t1'))'*mean(t1);
% 	imagesc(corr(t1'));
% 	axis image
% 	caxis([-1 1])
% 	axis off
% end
% colormap(cmap)


% cols={'--',':','-.'};

% figure('color','white')
% subplot(211)
% plot(G,(VE1),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[ve_data(k) ve_data(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('VE1','fontSize',18)
% set(gca,'fontSize',18)
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% subplot(212)
% plot(G,(GFC),'k','lineWidth',3)
% hold on
% for k=1:3
% 	line([G(1) G(end)],[GFC_exp(k) GFC_exp(k)],'LineStyle',cols{k},'Color','black','lineWidth',2)
% end
% legend({'AROMA+2P','AROMA+2P+GMR','AROMA+2P+DiCER'},'Location','northwest');
% xlabel('Global Coupling (G)','fontSize',18);
% ylabel('GFC','fontSize',18)
% set(gca,'fontSize',18)



% clear ts_simulated_all VE1;
% vals = {'G_0.05.mat','G_0.1.mat','G_0.15.mat','G_0.2.mat','G_0.25.mat','G_0.3.mat','G_0.35.mat','G_0.4.mat','G_0.45.mat','G_0.5.mat','G_0.55.mat','G_0.6.mat','G_0.65.mat','G_0.7.mat','G_0.75.mat','G_0.8.mat','G_0.85.mat','G_0.9.mat','G_0.95.mat','G_1.05.mat','G_1.1.mat','G_1.15.mat','G_1.2.mat','G_1.25.mat','G_1.3.mat','G_1.35.mat','G_1.4.mat','G_1.45.mat','G_1.5.mat','G_1.55.mat','G_1.6.mat','G_1.65.mat','G_1.7.mat','G_1.75.mat','G_1.8.mat','G_1.85.mat','G_1.9.mat','G_1.95.mat','G_1.mat','G_2.05.mat','G_2.1.mat','G_2.15.mat','G_2.2.mat','G_2.25.mat','G_2.3.mat','G_2.35.mat','G_2.4.mat','G_2.45.mat','G_2.5.mat','G_2.55.mat','G_2.6.mat','G_2.65.mat','G_2.7.mat','G_2.75.mat','G_2.8.mat','G_2.85.mat','G_2.9.mat','G_2.95.mat','G_2.mat','G_3.05.mat','G_3.1.mat','G_3.15.mat','G_3.2.mat','G_3.25.mat','G_3.3.mat','G_3.35.mat','G_3.4.mat','G_3.45.mat','G_3.5.mat','G_3.55.mat','G_3.6.mat','G_3.65.mat','G_3.7.mat','G_3.75.mat','G_3.8.mat','G_3.85.mat','G_3.9.mat','G_3.95.mat','G_3.mat','G_4.05.mat','G_4.1.mat','G_4.15.mat','G_4.2.mat','G_4.25.mat','G_4.3.mat','G_4.35.mat','G_4.4.mat','G_4.45.mat','G_4.5.mat','G_4.55.mat','G_4.6.mat','G_4.65.mat','G_4.7.mat','G_4.75.mat','G_4.8.mat','G_4.85.mat','G_4.9.mat','G_4.95.mat','G_4.mat','G_5.mat'};
% for j=1:100;
% 	G(j) = we;
% 	load(vals{j});
% 	ts_simulated_all(:,:,j) = bds;	
% 	grand_fc(j) = fitting;
% end

% ts_simulated_all = ts_simulated_all(:,:,1:5:end);
% G = G(1:5:end);
% grand_fc = grand_fc(1:5:end);


	
% for j=1:length(G),
% 	t1=zscore(ts_simulated_all(:,:,j))';
% 	cor_all=corr(t1');
% 	GFC(j) = mean(cor_all(:));
% 	[~,~,~,~,expl] = pca(t1');
% 	VE1(j) = expl(1);
% end



% NEW CODE:

% Prepare data:
for k=1:3,
	t1 = zscore(time_series(:,:,14,k));
	[~,~,~,~,expl] = pca(t1);
	co = corr(t1);
	GFC_exp(k) = mean(abs(co(:)));
	ve_data(k) = expl(1);
end

% BTF:
load('/Users/kevinaquino/projects/modelling_gustavo/results/BTF/ICA-AROMA.mat');
show_g=[1 9 13 14 15];
summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,1);

% B E-I:
cd('/Users/kevinaquino/projects/original_EI_code/');
clear ts_simulated_all VE1;
vals = {'G_0.05.mat','G_0.1.mat','G_0.15.mat','G_0.2.mat','G_0.25.mat','G_0.3.mat','G_0.35.mat','G_0.4.mat','G_0.45.mat','G_0.5.mat','G_0.55.mat','G_0.6.mat','G_0.65.mat','G_0.7.mat','G_0.75.mat','G_0.8.mat','G_0.85.mat','G_0.9.mat','G_0.95.mat','G_1.05.mat','G_1.1.mat','G_1.15.mat','G_1.2.mat','G_1.25.mat','G_1.3.mat','G_1.35.mat','G_1.4.mat','G_1.45.mat','G_1.5.mat','G_1.55.mat','G_1.6.mat','G_1.65.mat','G_1.7.mat','G_1.75.mat','G_1.8.mat','G_1.85.mat','G_1.9.mat','G_1.95.mat','G_1.mat','G_2.05.mat','G_2.1.mat','G_2.15.mat','G_2.2.mat','G_2.25.mat','G_2.3.mat','G_2.35.mat','G_2.4.mat','G_2.45.mat','G_2.5.mat','G_2.55.mat','G_2.6.mat','G_2.65.mat','G_2.7.mat','G_2.75.mat','G_2.8.mat','G_2.85.mat','G_2.9.mat','G_2.95.mat','G_2.mat','G_3.05.mat','G_3.1.mat','G_3.15.mat','G_3.2.mat','G_3.25.mat','G_3.3.mat','G_3.35.mat','G_3.4.mat','G_3.45.mat','G_3.5.mat','G_3.55.mat','G_3.6.mat','G_3.65.mat','G_3.7.mat','G_3.75.mat','G_3.8.mat','G_3.85.mat','G_3.9.mat','G_3.95.mat','G_3.mat','G_4.05.mat','G_4.1.mat','G_4.15.mat','G_4.2.mat','G_4.25.mat','G_4.3.mat','G_4.35.mat','G_4.4.mat','G_4.45.mat','G_4.5.mat','G_4.55.mat','G_4.6.mat','G_4.65.mat','G_4.7.mat','G_4.75.mat','G_4.8.mat','G_4.85.mat','G_4.9.mat','G_4.95.mat','G_4.mat','G_5.mat'};
for j=1:100;
	load(vals{j});
	G(j) = we;	
	ts_simulated_all(:,:,j) = bds;	
	grand_fc(j) = fitting;
end

ts_simulated_all = ts_simulated_all(:,:,1:5:end);
G = G(1:5:end);
grand_fc = grand_fc(1:5:end);
show_g=[1:4:20];

summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,1);


% HOPF
load('/Users/kevinaquino/projects/modelling_gustavo/results/HOPF+GLOBAL/ICA-AROMA.mat')
show_g=[1:4:20];
summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,1);

% for hh=1:20,
% 	[ts_simulated_all,grandFCcorr,bifpar,FCD] = run_hopf_model_homogenous_bif(sc_matrix,squeeze(permute(time_series(:,:,1:10,1),[2 1 4 3])),G,'');
% 	for j=1:length(G),
% 		t1=zscore(ts_simulated_all(:,:,j))';
% 		cor_all=corr(t1');
% 		GFC(hh,j) = mean(cor_all(:));		
% 		[~,~,~,~,expl] = pca(t1');
% 		VE1(hh,j) = expl(1);
% 	end
% end

% NDM:
clear ts_simulated_all;
t1 = zscore(time_series(:,:,14,1));
gs=mean(t1');
D=sum(C);
for gg=1:length(G),
	for j=1:length(D),
		ts_simulated_all(:,j,gg) = G(gg)*gs*D(j) + 0.5*randn(1,length(gs));
	end
end

summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,0);

load('/Users/kevinaquino/projects/modelling_gustavo/results/BEI/BEI_simulation_88.mat')
show_g=[1:4:20];
summary_plot(ts_simulated_all,G,show_g,ve_data,GFC_exp,0);

% Now something for FCD - maybe show non GSR and GSR in a big histogram
% load('/Users/kevinaquino/projects/modelling_gustavo/results/HOPF+GLOBAL/ICA-AROMA.mat')
figure;
t1 = zscore(time_series(:,:,14,1));
FCD_d = fcd_calculator(t1',10,5);
g_use=[1:4:20];
for j=1:5,
	
	subplot(2,6,j);
	sim_t=ts_simulated_all(:,:,g_use(j))';
	% FCD_1 = phase_fcd(sim_t,2);
	FCD_1 = fcd_calculator(sim_t,10,5);

	h=hist(FCD_1(find(triu(ones(size(FCD_1)),1))),50);
	[~,~,FCD1]=kstest2(FCD_d(find(triu(ones(size(FCD_d)),1))),FCD_1(find(triu(ones(size(FCD_1)),1))));
	bar(h/sum(h),'barWidth',1,'FaceColor','blue','FaceAlpha',0.5,'EdgeColor','none');
	title(['G=',num2str(G(g_use(j))),' FCD = ',num2str(round(FCD1*100)/100)]);


	sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	% FCD_2 = phase_fcd(sim_t,2);

	FCD_2 = fcd_calculator(sim_t,10,5);
	hold on
	subplot(2,6,j+6);
	h=hist(FCD_2(find(triu(ones(size(FCD_2)),1))),50);
	bar(h/sum(h),'barWidth',1,'FaceColor','red','FaceAlpha',0.5,'EdgeColor','none');
	[~,~,FCD2]=kstest2(FCD_d(find(triu(ones(size(FCD_d)),1))),FCD_2(find(triu(ones(size(FCD_2)),1))));
	title(['G=',num2str(G(g_use(j))),' FCD GSR = ',num2str(round(FCD2*100)/100)]);
end
subplot(2,6,[6 12]);
% FCD_d = phase_fcd(t1',2);
hist(FCD_d(find(triu(ones(size(FCD_d)),1))),50);
title('Data');

for k=1:3
	for j=1:100,
		FC(:,:,j,k) = corr(time_series(:,:,j,k));
	end;
	FC_mean(:,:,k) = nanmean(FC(:,:,:,k),3);
end


% Fitting corrs for BE
inds=find(triu(ones(size(FC_mean(:,:,1))),1));
for j=1:length(G),
	% Model vs ICA-AROMA data
	FC_model=corr(ts_simulated_all(:,:,j));
	FC_emp = FC_mean(:,:,1);
	sim(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));	

	% Model vs GSR data
	FC_emp = FC_mean(:,:,2);
	sim2(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	
	

	% Model GSR vs Data GSR
	sim_t=ts_simulated_all(:,:,j);
	sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	FC_model=corr(sim_t);
	sim3(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));

	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer_gsm(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
end

figure;plot(G,sim);
hold on;
plot(G,sim2);
plot(G,sim3);
plot(G,sim_dicer);
plot(G,sim_dicer_gsm);

legend({'BEI vs ICA-AROMA','BEI vs ICA-AROMA+GSR','BEI+GSR vs ICA-AROMA+GSR','BEI vs ICA-AROMA+DiCER','BEI + GSR vs ICA-AROMA+DiCER'});



% Fitting corrs for BTF



for k=1:3
	for j=1:100,
		tt=time_series(:,[1:34,42:41+34],j,k);
		FC(:,:,j,k) = corr(tt);
	end;
	FC_mean(:,:,k) = nanmean(FC(:,:,:,k),3);
end


inds=find(triu(ones(size(FC_mean(:,:,1))),1));
for j=1:length(G),
	% Model vs ICA-AROMA data
	FC_model=corr(ts_simulated_all(:,:,j));
	FC_emp = FC_mean(:,:,1);
	sim(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));	

	% Model vs GSR data
	FC_emp = FC_mean(:,:,2);
	sim2(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	
	

	% Model GSR vs Data GSR
	sim_t=ts_simulated_all(:,:,j);
	sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
	FC_model=corr(sim_t);
	sim3(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));

	% Dicer estimates
	FC_emp = FC_mean(:,:,3);
	sim_dicer_gsm(j)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
end

figure;plot(G,sim);
hold on;
plot(G,sim2);
plot(G,sim3);
plot(G,sim_dicer);
plot(G,sim_dicer_gsm);

legend({'BEI vs ICA-AROMA','BEI vs ICA-AROMA+GSR','BEI+GSR vs ICA-AROMA+GSR','BEI vs ICA-AROMA+DiCER','BEI + GSR vs ICA-AROMA+DiCER'});

legend({'MODEL vs ICA-AROMA','MODEL vs ICA-AROMA+GSR','MODEL+GSR vs ICA-AROMA+GSR','MODEL vs ICA-AROMA+DiCER','MODEL + GSR vs ICA-AROMA+DiCER'});

%%% Here show the fits summarized

% Here maybe show a summary?


% Be awesome FITS fs mFD?


% for sub=1:100,
% 	% FC_emp = FC(:,:,sub,2);
% 	for j=1:length(G),
% 		sim_t=ts_simulated_all(:,:,j);
% 		sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
% 		FC_model=corr(sim_t);		
% 		sim_fit(j,sub)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));		
% 	end
% 	% fit_max(sub) = max(sim_fit(:,sub));
% end



% % Prepare data:
% for sub=1:100,
% 	t1 = zscore(time_series(:,:,sub,1));
% 	[~,~,~,~,expl] = pca(t1);
% 	co = corr(t1);
% 	GFC_exp(sub) = mean(abs(co(:)));
% 	ve_data(sub) = expl(1);
% end

% nonNans=find(~isnan(fit_max));
% corr(ve_data(nonNans)',fit_max(nonNans)')


% Findthe best fit
[~,ind]=max(sim_dicer_gsm);
sim_t=ts_simulated_all(:,:,ind)';
sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
FC_model=corr(sim_t);
figure;imagesc(triu(FC_model,1) + tril(FC_mean(:,:,2),-1))
caxis([-1 1]);colormap(cmap);axis image;colorbar
