global_signal=movmean(randn(1,500),10);

nVoxels=1000;
allSignal=randn(nVoxels,500);

beta_weights=randn(nVoxels,1)*5 + 2;

newSignal=(global_signal'*(beta_weights'))' + allSignal;
figure;
% imagesc(zscore(newSignal',[],2));
imagesc(zscore(newSignal,[],2));
% caxis(2*[-1 1]);
colormap gray

% Maybe have to show this again using cluster ordering, GSO etc -- would be neat to see this!

time_pt=linspace(0,360,500);
% delay=(rand(nVoxels,1)-0.5)*2;
for j=1:nVoxels,
	% sig=interp1(time_pt,global_signal,time_pt+delay(j));
	newSignal2(j,:)=beta_weights(j)*(randn(1,500)*1 + global_signal)+allSignal(j,:);
end
figure;imagesc(zscore(newSignal2,[],2));
colormap gray
caxis(2*[-1 1]);

corr(nanmean(newSignal2)',global_signal')


% Might have to check again regarding the distributions (to see what kind of patterns are useful!)