% FC calcualtor

function fcd_est = fcd_calculator(data,windowWidth,overlap)

% Some parameters for the FCD calculation
if(nargin<3)
	overlap = 20;
end

if(nargin<2)
	windowWidth = 30;
end

nonOverlap = windowWidth - overlap;
% Calculate the number of windows:
N = size(data,1);
noWindows = floor((size(data,2)-windowWidth)/nonOverlap + 1);

for j=1:noWindows;
	window{j}=[1:windowWidth]+nonOverlap*(j-1);
end

% Lower Triangle
Isubdiag = find(tril(ones(N),-1));


for w_no=1:length(window),
	ts1 = zscore(data(:,window{w_no}),[],2);
	% ts1 = RegressNoiseSignal(ts1,mean(ts1));

	FC1 = corr(ts1.');
	for w_no2 = w_no:length(window),
		ts2 = zscore(data(:,window{w_no2}),[],2);
		% ts2 = RegressNoiseSignal(ts2,mean(ts2));
		FC2 = corr(ts2.');
		cor_all = corrcoef(FC1(Isubdiag),FC2(Isubdiag));
		fcd_est(w_no,w_no2) = cor_all(2);
	end
end

% Make it a symmetric matrix - remembering to get rid of the double diagonal.
fcd_est = fcd_est + fcd_est.' - eye(size(fcd_est,1));