% This here runs the noisy degree model, and might seem clunky how it is set up but it is needed here so that it is 
% consistent with the other models when they are being run on a cluster (or just killed because they have been running too long)
function [ts_simulated_all,FCcorr,grandFCcorr,FCD] =run_noisy_degree_model(sc_matrix,time_series,G,folder)


% We first grab the global signal (if there is one) form the estimated model, then find the optimal model for the G-value.

% Do it for the grand model, then look at the individual ones i think... 

% First do it for one subject, not sure what it corresponds to in terms of a G-value


% Now find the degree, as this is the most important part of the model
D = mean(sc_matrix,1);
TR = 2;
Tmax = size(time_series,2);
Nregions=size(time_series,1);
Nsubs=size(time_series,3);
% Here to make sure we exclude the main diagonal
upperTriangle = find(triu(ones(size(sc_matrix)),1));
phfcddata = [];


for subject=1:Nsubs,
	ts = time_series(:,:,subject);
	% Find the equivalent G
	corrFC(:,:,subject) = corr(ts.');
	corrFC_emp = corrFC(:,:,subject);		
	for coupling_index=1:length(G),

		% Coupling index here:
		if(G(coupling_index) == 0)
			coupling_factor = 0;
		else
			coupling_factor = max(D)/G(coupling_index);
		end

		for j=1:Nregions;
			ts_simulated(j,:) = mean(ts,1)*D(j) + coupling_factor*randn(1,Tmax);
		end
		ts_simulated_all(:,:,subject,coupling_index) = ts_simulated;
		% Look at the simulated data:
		corrFC_sim = corr(ts_simulated.');
		FCcorr(coupling_index,subject) = corr(corrFC_sim(upperTriangle),corrFC_emp(upperTriangle));
		phfcd(coupling_index,:,subject) = phase_fcd(ts_simulated,TR);
	end
	phfcddata = [phfcddata,phase_fcd(ts,TR)];	
end


% Now work out the fcd for all subjects and put it together
for coupling_index=1:length(G),
	phfcd_all = phfcd(coupling_index,:,:);
	phfcd_all = phfcd_all(:);
	[~,~,FCD(coupling_index)] = kstest2(phfcd_all,phfcddata);
end


% Here now we look at the average as well -- cool to look at both really now it doesnt look as good maybe will have to make 
% a number of GS

meancorrFC = mean(corrFC,3);
globalSignal = movmean(randn(1,Tmax),3);

for coupling_index=1:length(G),	
	if(G(coupling_index) == 0)
		coupling_factor = 0;
	else
		coupling_factor = max(D)/G(coupling_index);
	end
	for j=1:Nregions;
		ts_simulated(j,:) = globalSignal*D(j) + coupling_factor*randn(1,Tmax);
	end
	% Look at the simulated data:
	corrFC_sim = corr(ts_simulated.');	
	grandFCcorr(coupling_index) = corr(corrFC_sim(upperTriangle),meancorrFC(upperTriangle));
	% phfcd = phase_fcd(ts_simulated,TR);
% 	keyboard
	% [~,~,FCD(coupling_index)] = kstest2(phfcd,phfcddata);
end

% Question should we be calculating the ks-test for all subjects or just the grand average?