function run_BTF_model(C);


% Function here to solve BTF, solve it over a long time
% C = C/max(C(:))*0.2;
% C = C - diag(diag(C));

% Initial set up with default parameters:
% sys = BTF2003SDE(C*3);
sys = BTF2003(C);
% keyboard

% Solve the model for 2 seconds
sol = bdSolve(sys,[0 2000]);
for j=1:size(C,1)*3,
	y_all(j,:) = interp1(sol.x,sol.y(j,:),[1:1:2000]);
end


Y = reshape(y_all,size(C,1),3,2000);

% Take the last values
lastValues = Y(:,:,end);
Vtest = squeeze(Y(:,1,:));


for n=1:size(C,1),
	% inputNeural(n,:) = movmean(abs(gradient(Vtest(n,:),0.1)),2500);
	inputNeural(n,:) = abs(gradient(Vtest(n,:),1));
	% inputNeural(n,:) = filtfilt(bfilt2,afilt2,(gradient(Vtest(n,:),0.1)));
end

% I have to low pass probably as a first pass take moving average of 250 ms

Vin = inputNeural;

% keyboard;

% W = squeeze(Y(:,2,1:2500:end));
% Z = squeeze(Y(:,3,1:2500:end));

% Solve this for 147*2 more seconds;
for n=1:10,	
	disp(['Iteration: ',num2str(n),'....']);
	% Now set the last point as the intial value for the simulation to continue the simulation
	sys2 = sys;
	for j=1:3,
		sys2.vardef(j).value = lastValues(:,j);
	end
	% Solve again for 1 more second
	sol2 = bdSolve(sys2,[0 2000]);
	% Y = reshape(sol2.y,68,3,20001);

	for j=1:size(C,1)*3,
		y_all(j,:) = interp1(sol.x,sol.y(j,:),[1:1:2000]);
	end


	Y = reshape(y_all,size(C,1),3,2000);


	Vtest = squeeze(Y(:,1,:));

	for n=1:68,
		% inputNeural(n,:) = movmean(abs(gradient(Vtest(n,:),0.1)),2500);
		inputNeural(n,:) = abs(gradient(Vtest(n,:),0.1));
		% inputNeural(n,:) = filtfilt(bfilt2,afilt2,(gradient(Vtest(n,:),0.1)));
	end

	Vin = [Vin,inputNeural];
	lastValues = Y(:,:,end);

	
	% W = [W, squeeze(Y(:,2,1:2500:end))];
	% Z = [Z, squeeze(Y(:,3,1:2500:end))];

end

keyboard

% dt = 2000/20000*250/1000;
% 1 ms
dt = 1;

time = 0:dt:(size(Vin,2)-1)*dt;
% time_sampled = time(1:250:end);
% time = 0:dt:(size(Vin,2)-1)*dt;


% % Now model the HRF with a standard two-gamma distribution
% modelHrf = gampdf(time/1e3, 6, 1) - gampdf(time/1e3, 16, 1)/6;		
% hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';
% % keyboard;
% TRs = 0:2:40*2;
% % Convolve the BOLD response 
% for n=1:68,
% 	BOLD(n,:) = conv(Vin(n,:),hrf,'same');
% 	BOLD_interp(n,:) = interp1(time/1e3,BOLD(n,:),TRs);
% end


% BOLD_interp(isnan(BOLD_interp)) = 0;

% BOLD = BOLD_interp;

% bb = BOLD_interp(:,14:end);
keyboard
for n=1:68,b(n,:) = BOLD(80,Vin(n,:));end
figure;
imagesc(corr(b.'));

keyboard;