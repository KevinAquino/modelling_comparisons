%sillytest.m

for subject=1:24,
	for j=1:68,
		pows = abs(fft(time_series(j,:,subject,5)));
		ratio(j,subject) = sum(pows(1:37))/sum(pows(38:74));
	end
end
