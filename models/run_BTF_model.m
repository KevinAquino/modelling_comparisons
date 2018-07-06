function run_BTF_model(C);


% Function here to solve BTF, solve it over a long time
C = C/max(C(:))*0.2;
C = C - diag(diag(C));

% Initial set up with default parameters:
sys = BTF2003SDE(C*3);


% Solve the model for 1 second
sol = bdSolve(sys,[0 1000]);
Y = reshape(sol.y,68,3,10001);


lastValues = Y(:,:,end);

Vtest = squeeze(Y(:,1,:));


% TR = 0.1e-3;
% k=2;                  % 2nd order butterworth filter
% fnq=1/(2*TR);
% flp = 100;            % low frequency of the filter
% Wn=[flp/fnq]; 		 % butterworth bandpass non-dimensional frequency
% [bfilt2,afilt2]=butter(k,Wn);   % construct the filter
% % signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));


for n=1:68,
	% inputNeural(n,:) = movmean(abs(gradient(Vtest(n,:),0.1)),2500);
	inputNeural(n,:) = abs(gradient(Vtest(n,:),0.1));
	% inputNeural(n,:) = filtfilt(bfilt2,afilt2,(gradient(Vtest(n,:),0.1)));
end

% I have to low pass probably as a first pass take moving average of 250 ms

Vin = inputNeural(:,1:250:end);

% W = squeeze(Y(:,2,1:2500:end));
% Z = squeeze(Y(:,3,1:2500:end));

% Solve this for 147 more seconds;
for n=1:147*2,
	% Now set the last point as the intial value for the simulation to continue the simulation
	sys2 = sys;
	for j=1:3,
		sys2.vardef(j).value = lastValues(:,j);
	end
	% Solve again for 1 more second
	sol2 = bdSolve(sys2,[0 1000]);
	Y = reshape(sol2.y,68,3,10001);
	Vtest = squeeze(Y(:,1,:));

	for n=1:68,
		inputNeural(n,:) = movmean(abs(gradient(Vtest(n,:),0.1)),2500);
		% inputNeural(n,:) = filtfilt(bfilt2,afilt2,(gradient(Vtest(n,:),0.1)));
	end

	Vin = [Vin,inputNeural(:,1:250:end)];
	lastValues = Y(:,:,end);

	
	% W = [W, squeeze(Y(:,2,1:2500:end))];
	% Z = [Z, squeeze(Y(:,3,1:2500:end))];

end


dt = 250/10000;
time = 0:dt:(size(Vin,2)-1)*dt;



% Now model the HRF with a standard two-gamma distribution
modelHrf = gampdf(time, 6, 1) - gampdf(time, 16, 1)/6;		
hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';

TRs = 0:2:147*2;
% Convolve the BOLD response 
for n=1:68,
	BOLD(n,:) = conv(Vin(n,:),hrf,'same');
	BOLD_interp(n,:) = interp1(time,BOLD(n,:),TRs);
end


BOLD_interp(isnan(BOLD_interp)) = 0;

BOLD = BOLD_interp;

bb = BOLD_interp(:,35:end);
figure;imagesc(corr(bb.'));