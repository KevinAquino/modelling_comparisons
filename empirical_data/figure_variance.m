

for ng=1:length(G,
	[ COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(ts_simulated_all(:,:,ng));
	ve1(ng) = EXPLAINED(1);
end

% Maybe also calculate the mean FC?

for subject=14,
	for np=1:3,
		[ COEFF, SCORE, LATENT, TSQUARED, EXPLAINED, MU] = pca(time_series(:,:,subject,np));
		ve1_emp(subject,np) = EXPLAINED(1);
	end
end