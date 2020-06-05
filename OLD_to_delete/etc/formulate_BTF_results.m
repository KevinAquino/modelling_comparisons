% This here is to just get the BTF into a form that can be put into a figure.



% Here we have the four processing streams that we want to model the responses to
% preprocessing_stream ={'MINIMAL','ICA-AROMA','ICA-AROMA_GSR','ICA-AROMA_DBSCAN'};
preprocessing_stream ={'ICA-AROMA','ICA-AROMA_GSR','ICA-AROMA_DBSCAN'};
% preprocessing_stream ={'ICA-AROMA'};

% Load in the time series:
load('/Users/kevinaquino/projects/modelling_gusatvo/empirical_data/ten_subjects_aparc_hcp.mat');
% Load in the structural matrix:
load('/Users/kevinaquino/projects/modelling_gusatvo/empirical_data/exemplarSC.mat');


TR = 2;
for prepro = preprocessing_stream,
	GSR_flag = 0;
	G=[0:0.05:1];
	switch prepro{1}
		case 'MINIMAL'
			time_series = time_series_aparc(:,:,:,4);
		case 'ICA-AROMA'
			time_series = time_series_aparc(:,:,:,1);
		case 'ICA-AROMA_GSR'
			time_series = time_series_aparc(:,:,:,2);
			GSR_flag = 1;
		case 'ICA-AROMA_DBSCAN'
			time_series = time_series_aparc(:,:,:,3);
	end
	subjects = size(time_series_aparc,3);
	% Using the same way to calculate FCD as been done in all the other scripts
	phfcddata = [];
	for subject=1:subjects,
		ts = time_series(:,:,subject);
		phfcddata = [phfcddata,phase_fcd(ts,TR)];	
		corrFC(:,:,subject) = corr(squeeze(time_series(:,:,subject)).');	
	end

	meancorrFC = mean(corrFC,3);

	upperTriangle = find(triu(ones(size(C)),1));

	BOLD = [];
	for coupling_index = 1:length(G)
		load(['/Users/kevinaquino/projects/modelling_gusatvo/results_testing/',num2str(100*G(coupling_index)),'_BTF.mat']);
		hh = zscore(BOLD(:,100:20:end),[],2);
		if(GSR_flag)
			gs=mean(hh);
			hh = RegressNoiseSignal(hh,gs);
		end
		% This processing of the time series is TEMPRORARY will change when the code above is cleaned up!
		phfsim = phase_fcd(hh,2);
		[~,~,FCD(coupling_index)] = kstest2(phfsim,phfcddata);
		corrFC_sim = corr(hh.');
		grandFCcorr(coupling_index) = corr(corrFC_sim(upperTriangle),meancorrFC(upperTriangle));
		ts_simulated_all(:,:,coupling_index) = hh.';

	end
	G = G/0.2;
	save(['/Users/kevinaquino/projects/modelling_gusatvo/results/BTF/',prepro{1}],'G','grandFCcorr','FCD','ts_simulated_all');
end 



