% FC calcualtor

for j=1:61;
	window{j}=[1:30]+10*(j-1);
end
N = 68;
Isubdiag = find(tril(ones(N),-1));

for subject=1:24,
	data = subject_ts{subject};
	for w_no=1:length(window),
		ts1 = zscore(data(:,window{w_no}),[],2);
		% ts1 = RegressNoiseSignal(ts1,mean(ts1));

		FC1 = corr(ts1.');
		for w_no2 = w_no:length(window),
			ts2 = zscore(data(:,window{w_no2}),[],2);
			% ts2 = RegressNoiseSignal(ts2,mean(ts2));
			FC2 = corr(ts2.');
			cor_all = corrcoef(FC1(Isubdiag),FC2(Isubdiag));
			fcd_est(w_no,w_no2,subject) = cor_all(2);
		end
	end
end

fcd_all = mean(fcd_est,3);
% fcd_all = 0.5*(fcd_all + fcd_all.');

figure;imagesc(fcd_all);
caxis([0 1]);
colormap jet