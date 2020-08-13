% time_scale_analyses.m
%
% function time_scales()

% Load up the time series
load('GenCog_ALL_DiCER.mat');
% Next load up the HCP
load('/Users/kevinaquino/projects/modelling_gustavo/empirical_data/t1t2Cortex_68_HCP_surface_data.mat');

% Now look up time scales (for each region)

% Will do now for just 100 subjects to just check if i can replicate the findings.
region_list = [1:34,42:76];
for ns=1:440,
	ts=squeeze(time_series(:,:,ns,1));
	% Now look for every subject
	for region=1:68,
		nr = region_list(region);
		out = CO_AutoCorrShape(ts(nr,:),'posDrown');
		tau1(region,ns) = out.decayTimescale; % Murray
		tau2(region,ns) = out.sumacf; % Watanabe
	end

end
% end

figure;plot(t1t2Cortex,mean(tau2,2),'.','MarkerSize',18);


% Now look at using Ito's method
s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.6 0 3],'Lower',[0.2 -1 1],'Upper',[1 0.1 10]);
f = fittype('A*(exp(-x/tau) + C)','options',s);

region_list = [1:34,42:76];
for ns=1:440,
	ts=squeeze(time_series(:,:,ns,1));
	% Now look for every subject
	for region=1:68,
		nr = region_list(region);
		ac=autocorr(ts(nr,:),100);
		% Look here at the fitting

	    [c, gof] = fit([1:100]',ac(2:end)',f);
	    tau3(region,ns) = c.tau;
		% out = CO_AutoCorrShape(ts(nr,:),'posDrown');
		% tau1(region,ns) = out.decayTimescale; % Murray
		% tau2(region,ns) = out.sumacf; % Watanabe
	end

end

figure;plot(t1t2Cortex,mean(tau3,2),'.','MarkerSize',18);



figure('color','white'); 
subplot(121);
plot(t1t2Cortex,mean(tau2,2),'.','MarkerSize',18);
ylabel('$\tau$','Interpreter','LaTeX','fontSize',18);
xlabel('$T_1/T_2$','Interpreter','LaTeX','fontSize',18);
title('HCTSA')

subplot(122);
plot(t1t2Cortex,mean(tau3,2),'.','MarkerSize',18);
ylabel('$\tau$','Interpreter','LaTeX','fontSize',18);
xlabel('$T_1/T_2$','Interpreter','LaTeX','fontSize',18);
title('Ito')


figure;hold on;
plot(t1t2Cortex,mean(tau2,2),'.','MarkerSize',18);
plot(t1t2Cortex,mean(tau3,2),'.','MarkerSize',18);
legend({'HCTSA','Ito'})



figure('color','white'); 
subplot(121);
hold on
plot(t1t2Cortex,mean(tau2,2),'.','MarkerSize',18);
errorbar(t1t2Cortex,mean(tau2,2),std(tau2,[],2),'LineStyle','none')
ylabel('$\tau$','Interpreter','LaTeX','fontSize',18);
xlabel('$T_1/T_2$','Interpreter','LaTeX','fontSize',18);
title('HCTSA')

subplot(122);
hold on
plot(t1t2Cortex,mean(tau3,2),'.','MarkerSize',18);
errorbar(t1t2Cortex,mean(tau3,2),std(tau3,[],2),'LineStyle','none')
ylabel('$\tau$','Interpreter','LaTeX','fontSize',18);
xlabel('$T_1/T_2$','Interpreter','LaTeX','fontSize',18);
title('Ito')