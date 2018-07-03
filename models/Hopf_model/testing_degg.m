nt = 660;
rn = 0.5;
D = sum(C,1);
signal = movmean(randn(nt,1),5);
signal = randn(nt,1);
clear bold_simple
noiseWeight = 0.1;	

for n=1:length(D);
		bold_simple(:,n) = (D(n))*(signal) + noiseWeight*randn(size(signal));
end
cor_test = corr(bold_simple + rn*randn(size(bold_simple)));
figure;
imagesc(cor_test)
caxis([0 1]);colormap jet

fcd_model_simple = fcd_calculator(bold_simple.' +  rn*randn(size(bold_simple.')),30,20);

xlb = linspace(0,660*2,size(fcd_model_simple,1));

iss = find(tril(ones(size(fcd_model_simple,1)),-1));
figure('color','white');subplot(211);imagesc(xlb,xlb,fcd_model_simple);caxis([0 1]);colormap jet
axis image
xlabel('time shift (s)');ylabel('time shift (s)');
title(['$\beta = ' num2str(rn) '$'],'Interpreter','LaTeX');
set(gca,'fontSize',18)
subplot(212)
hist(fcd_model_simple(iss),100)
xlabel('FCD');ylabel('count');
set(gca,'fontSize',18)
xlim([0 1])

% sig1 = randn(nt,1);
% for j=2:6;sigg{j} = movmean(sig1,j);end
% sigg{1} = sig1;
% figure('Color','white');for j=1:6;subplot(6,1,j);plot(sigg{j});title(['Moving average window = ' num2str(j-1)]);set(gca,'fontSize',24);end