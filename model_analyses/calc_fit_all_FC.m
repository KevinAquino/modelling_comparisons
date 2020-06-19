function all_fits = calc_fit_all_FC(FC_mean,ts_simulated_all,G)

% Fitting corrs for BE
inds=find(triu(ones(size(FC_mean(:,:,1))),1));
for NR=1:size(ts_simulated_all,4),
	for j=1:length(G),
		% Model vs ICA-AROMA data
		FC_model=corr(ts_simulated_all(:,:,j,NR));
		FC_emp = FC_mean(:,:,1);
		all_fits(1,j,NR)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
		% Dicer estimates
		FC_emp = FC_mean(:,:,3);
		all_fits(2,j,NR)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));	

		% Model vs GSR data
		FC_emp = FC_mean(:,:,2);
		all_fits(3,j,NR)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
		
		

		% Model GSR vs Data GSR
		sim_t=ts_simulated_all(:,:,j,NR)';
		sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
		FC_model=corr(sim_t');
		all_fits(4,j,NR)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));

		% Dicer estimates
		FC_emp = FC_mean(:,:,3);
		all_fits(5,j,NR)=corr(atanh(FC_model(inds)),atanh(FC_emp(inds)));
	end
end

% keyboard
all_fits=mean(all_fits,3);