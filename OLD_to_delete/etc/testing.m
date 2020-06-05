
bb = BOLD_interp(:,14:end);
val2= bb(:,1e4:end);
gs = mean(val2);
hh = RegressNoiseSignal(val2,gs);
figure;imagesc(hh)
d2 = (zscore(hh,[],2));
figure;
imagesc(corr(d2.'))