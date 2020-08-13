function [taus] = calc_tscale(time_series,t1t2Cortex)

	nTot = size(time_series,3);

	% Will do now for just 100 subjects to just check if i can replicate the findings.
	if(nTot>1),
		region_list = [1:34,42:76];
	else
		region_list = [1:68];
	end

	for ns=1:nTot,
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

	% figure;plot(t1t2Cortex,mean(tau2,2),'.','MarkerSize',18);


	% Now look at using Ito's method
	s = fitoptions('Method','NonlinearLeastSquares','StartPoint',[0.6 0 3],'Lower',[0.2 -1 1],'Upper',[1 0.1 10]);
	f = fittype('A*(exp(-x/tau) + C)','options',s);
	
	for ns=1:nTot,
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
	taus = [tau1 tau2 tau3];
	% figure;plot(t1t2Cortex,mean(tau3,2),'.','MarkerSize',18);
