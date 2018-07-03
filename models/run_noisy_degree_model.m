% This here runs the noisy degree model, and might seem clunky how it is set up but it is needed here so that it is 
% consistent with the other models when they are being run on a cluster (or just killed because they have been running too long)
function [ts_simulated_all,FCcorr,grandFCcorr] =run_noisy_degree_model(sc_matrix,time_series,G,folder)


% We first grab the global signal (if there is one) form the estimated model, then find the optimal model for the G-value.

% Do it for the grand model, then look at the individual ones i think... 

% First do it for one subject, not sure what it corresponds to in terms of a G-value


% Now find the degree, as this is the most important part of the model
D = mean(sc_matrix,1);


% Here to make sure we exclude the main diagonal
upperTriangle = find(triu(ones(size(sc_matrix)),1));


for subject=1:10,
	for coupling_index=1:length(G),
		% Find the equivalent G
		ts = time_series(:,:,subject);
		corrFC(:,:,subject) = corr(ts.');
		coupling_factor = max(D)/G(coupling_index);
		for j=1:68;
			ts_simulated(j,:) = mean(ts,1)*D(j) + coupling_factor*randn(1,147);
		end
		ts_simulated_all(:,:,subject,coupling_index) = ts_simulated;

		% Look at the simulated data:
		corrFC_sim = corr(ts_simulated.');
		corrFC_emp = corrFC(:,:,subject);
		FCcorr(coupling_index,subject) = corr(corrFC_sim(upperTriangle),corrFC_emp(upperTriangle));
	end
end

% Here now we look at the average as well -- cool to look at both really!

meancorrFC = mean(corrFC,3);
globalSignal = movmean(randn(1,147),3);

for coupling_index=1:length(G),	
	coupling_factor = max(D)/G(coupling_index);
	for j=1:68;
		ts_simulated(j,:) = globalSignal*D(j) + coupling_factor*randn(1,147);
	end
	% Look at the simulated data:
	corrFC_sim = corr(ts_simulated.');	
	grandFCcorr(coupling_index) = corr(corrFC_sim(upperTriangle),meancorrFC(upperTriangle));
end