% just a results figure:

figure('color','white');hold on;
plot(G,grandFCcorr_test,'r*-','lineWidth',2);plot(G,FCD_test,'r+--');
plot(G,grandFCcorr_retest,'k*-','lineWidth',2);plot(G,FCD_retest,'k+--');
legend({'No ANEC corr','No ANEC FCD','ANEC Retest corr','ANEC Retest FCD'});
xlabel('Global coupling');ylabel('FC similarity or FCD');set(gca,'fontSize',18)

figure('color','white');plot(G,FC_fit,'r*-');hold on;plot(G,FCD,'k*-');xlabel('Global coupling');ylabel('FC similarity');
title('ANEC optimization');legend('FC','FCD');
set(gca,'fontSize',18);

% Show Differences
uT=find(triu(ones(82,82),1));

cmap = [flipud(BF_getcmap('blues',9,0));[1,1,1],;BF_getcmap('reds',9,0)].';
for ns=1:size(ts2,3),meanCor_rt(:,:,ns)=corr(ts2(:,:,ns));end
meanCor_rt=mean(meanCor_rt,3);	
meanCor_rt(ii) = NaN;

corr_sim_test=corr(ts_simulated_all_test(:,:,2)');
corr_sim_test(ii) = NaN;
compMat = triu(meanCor_rt,1) + tril(corr_sim_test,-1);

figure('color','white');imagesc(atanh(compMat));axis image;caxis([-0.3 0.3]);colormap(cmap');colorbar;title('FCs no ANEC vs retestData');set(gca,'fontSize',18);
hf=figure('color','white');plot(corr_sim_test(uT),meanCor_rt(uT),'r.');hold on;

% Now with Anec
corr_sim_test=corr(ts_simulated_all_retest(:,:,2)');
corr_sim_test(ii) = NaN;
compMat = triu(meanCor_rt,1) + tril(corr_sim_test,-1);
figure('color','white');imagesc(atanh(compMat));axis image;caxis([-0.3 0.3]);colormap(cmap');colorbar;title('FCs retest ANEC vs retestData');set(gca,'fontSize',18);
figure(hf);plot(corr_sim_test(uT),meanCor_rt(uT),'k.');
legend({'No ANEC','ANEC retest'});xlabel('atanh(FCsim)');ylabel('atanh(FCemp)');set(gca,'fontSize',18);
% Original with ANEC
corr_sim_test=corr(ts_simulated_all(:,:,2)');
compMat = triu(meanCor,1) + tril(corr_sim_test,-1);
compMat(ii) = NaN;
figure('color','white');imagesc(atanh(compMat));axis image;caxis([-0.3 0.3]);colormap(cmap');colorbar;title('FC ANEC FC training');set(gca,'fontSize',18);


cc=[[1,1,1];BF_getcmap('reds',9,0)];
sc_comp = triu(sc_matrix,1) + tril(all_anecs{2},-1);
figure('color','white');imagesc(sc_comp)
colormap(cc);caxis([0 0.01]);
colorbar;
title('SC vs ANEC');
set(gca,'fontSize',18);axis image;