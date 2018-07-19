function run_BTF_model(C,cparam)


% Function here to solve BTF, solve it over a long time
% C = C/max(C(:))*0.2;
% C = C - diag(diag(C));

% Initial set up with default parameters:
% sys = BTF2003SDE(C*3);
sys = BTF2003(C);
sys.pardef=bdSetValue(sys.pardef,'C',cparam);
% keyboard

% Solve the model for 40 seconds to remove transients
sol = bdSolve(sys,[0 40000]);
% for j=1:size(C,1)*3,
% 	y_all(j,:) = interp1(sol.x,sol.y(j,:),[1:1:40000]);
% end

% keyboard
Y = reshape(sol.y(:,end),size(C,1),3,1);
for j=1:3,	
	sys.vardef(j).value = Y(:,j);
end

% After this, let the system run.
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
Vtin = [];
% W = squeeze(Y(:,2,1:2500:end));
% Z = squeeze(Y(:,3,1:2500:end));

% Solve this for 147*2 more seconds;
for n=1:200,	
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
		inputNeural(n,:) = abs(gradient(Vtest(n,:),1));
		% inputNeural(n,:) = filtfilt(bfilt2,afilt2,(gradient(Vtest(n,:),0.1)));
	end

	Vin = [Vin,inputNeural];
	Vtin = [Vtin,Vtest];
	lastValues = Y(:,:,end);

	
	% W = [W, squeeze(Y(:,2,1:2500:end))];
	% Z = [Z, squeeze(Y(:,3,1:2500:end))];

end

% keyboard

% dt = 2000/20000*250/1000;
% 1 ms
dt = 1;

% time = 0:dt:(size(Vin,2)-1)*dt;
% time_sampled = time(1:250:end);
time = 0:dt:(size(Vin,2)-1)*dt;

% 
% % % Now model the HRF with a standard two-gamma distribution
% modelHrf = gampdf(time/1e3, 6, 1) - gampdf(time/1e3, 16, 1)/6;		
% hrf = circshift(modelHrf.',round(length(modelHrf)/2)).';
% % % keyboard;
% % TRs = 0:0.01:100*4;
% % % Convolve the BOLD response 
% for n=1:68,
% 	BOLD(n,:) = conv(Vin(n,:),hrf,'same');
% % 	BOLD_interp(n,:) = interp1(time(1e5:end)/1e3,zscore(BOLD(n,1e5:end)),TRs);
% end

for n=1:68,
	out=tobold(Vin(n,:),time/1e3);
	bold_region(n,:) = out.y;
end

disp('Here will have to test the convolution model vs the Balloon model this might be a missing ingredient. ')
keyboard;
TRs = 0:0.01:100*4;
% % Downsample the BOLD response 
% for n=1:68,	
% 	BOLD_interp(n,:) = interp1(time(1e5:end)/1e3,zscore(BOLD(n,1e5:end)),TRs);
% end

gs = mean(BOLD_interp);
% keyboard
noNans = find(~isnan(gs));

% keyboard;
% BOLD_interp(isnan(BOLD_interp)) = 0;

% BOLD = BOLD_interp;

% bb = BOLD_interp(:,14:end);
% 
% for n=1:68,b(n,:) = BOLD(200,Vin(n,:));end
% zt = b(:,1e5:end);
% zt = zscore(zt,[],2);
% figure;
% imagesc(corr(zt.'));
figure;imagesc(corr(BOLD_interp(:,noNans).'))
caxis([0 1]);
title(['G = ',num2str(cparam/0.2)]);
axis image;
drawnow;
