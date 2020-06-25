ts_simulated = zeros(153,82);
inds_good=cell(20,1);
inds_bad=cell(20,1);

for g_ind=1:20,
	for RUN=1:108,
		file_name=['G_ind_',num2str(g_ind),'_RUN_',num2str(RUN),'simulation.mat'];
		if(isfile(file_name))
			load(file_name);
			ts_simulated_all_mat(:,:,g_ind,RUN) = ts_simulated;
			% n=n+1;
			inds_good{g_ind} = [inds_good{g_ind},RUN];
		else
			inds_bad{g_ind} = [inds_bad{g_ind},RUN];
			ts_simulated_all_mat(:,:,g_ind,RUN)=ts_simulated*0;
		end
		
	end
end


% Work out how many left
for g_ind=1:20,
	num_completed(g_ind)=length(inds_good{g_ind});
end

min_completed=min(num_completed);


for g_ind=1:20,
	ts_all(:,:,g_ind,:) = ts_simulated_all_mat(:,:,g_ind,inds_good{g_ind}(1:min_completed));
end

ts_simulated=ts_all;
save('./../simulation.mat','ts_simulated','simulation_params');