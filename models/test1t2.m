
upperTriangle = find(triu(ones(size(sc_matrix)),1));

normD=D/max(D);
t1t2_ctx = t1t2_ratio([1:34 42:(42+33)]);
t1t2_ctx=t1t2_ctx-min(t1t2_ctx);
t1t2_ctx=t1t2_ctx/max(t1t2_ctx);
% t1t2_ctx = t1t2_ctx/max(t1t2_ctx)*max(D);
sig=movmean(randn(1,Tmax),4);

Nregions=size(t1t2_ctx,1);
ts = time_s([1:34 42:(42+33)],:,subject);
clear ts_simulated fcorr_1 fcorr_2;
allk=linspace(-1,1,50);

for k=1:length(allk)
	for coupling_index=1:length(G),
		% Coupling index here:
		if(G(coupling_index) == 0)
			coupling_factor = 0;
		else
			coupling_factor = 1/G(coupling_index);
		end

		for j=1:Nregions;
					ts_simulated(j,:) = zscore(mean(ts,1))*1*(normD(j)+allk(k)*t1t2_ctx(j)) + coupling_factor*randn(1,Tmax);% + allk(k)*t1t2_ctx(j)*zscore(sig);
					ts_simulated2(j,:) = mean(ts,1)*1*normD(j) + coupling_factor*randn(1,Tmax);
					% ts_simulated(j,:) = mean(ts,1)*D(j) + coupling_factor*randn(1,Tmax);
		end
		% figure;subplot(131);imagesc(corr(ts_simulated'));title('ND+T1T2');axis image;subplot(132);imagesc(corr(ts_simulated2'));title('ND');axis image;subplot(133);imagesc(corr(ts'));title('EMP');axis image;
		% caxis([-0.5 0.5]);
		cemp=corr(ts');
		c1=corr(ts_simulated');
		c2=corr(ts_simulated2');	
		fcorr_1(coupling_index,k)=corr(c1(upperTriangle),cemp(upperTriangle));
		fcorr_2(coupling_index,k)=corr(c2(upperTriangle),cemp(upperTriangle));
	end
end

% figure;plot(G,fcorr_1,'r.-');hold on;plot(G,fcorr_2,'k.-');legend({'ND+T1T2','ND'});

figure;subplot(211);imagesc(fcorr_1);colorbar;subplot(212);imagesc(fcorr_2);colorbar;

max(fcorr_1(:))
max(fcorr_2(:))