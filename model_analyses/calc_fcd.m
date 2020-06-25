function all_FCDS=calc_fcd(G,ts_simulated_all,FCD_data_all)

	for j=1:length(G),

		for NR=1:size(ts_simulated_all,4),
			sim_t=ts_simulated_all(:,:,j,NR)';
			FCD_model(:,NR) = phase_fcd(sim_t',2);			
			sim_t=sim_t-(pinv(mean(sim_t))'*(sim_t'))'*mean(sim_t);
			FCD_model_GSR(:,NR) = phase_fcd(sim_t',2);
		end
			% Model vs AROMA
			[~,~,ksVal1] = kstest2(FCD_model(:),FCD_data_all(:,1));
			[~,~,ksVal2] = kstest2(FCD_model(:),FCD_data_all(:,2));
			[~,~,ksVal3] = kstest2(FCD_model_GSR(:),FCD_data_all(:,2));
			[~,~,ksVal4] = kstest2(FCD_model(:),FCD_data_all(:,3));

			all_FCDS(j,1) = ksVal1;
			all_FCDS(j,2) = ksVal2;
			all_FCDS(j,3) = ksVal3;
			all_FCDS(j,4) = ksVal4;			
		end
end