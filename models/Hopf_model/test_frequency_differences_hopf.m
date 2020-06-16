% A script to check to see if this is unchanged with GSR or not.
all_subjects = empirical_params.time_series;
clear f_diff;
for pp=1:3,
	for j=1:size(all_subjects,3)
		time_series=all_subjects(:,:,j,pp);
		f_diff(:,j,pp) = calculate_peak_frequency(squeeze(time_series)',empirical_params.TR);
	end
end


figure;
hold on;
for pp=1:3,
	f_pp=f_diff(:,:,pp);
	mn(:,pp)=mean(f_pp,2);
	er=std(f_pp,[],2);
	errorbar([1:82],mn(:,pp),er/sqrt(108));	
end
legend(noiseOptions)
% Shows that they are overlapping really, no need to worry about this.