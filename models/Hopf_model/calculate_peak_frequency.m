function f_diff = calculate_peak_frequency(time_series,TR)

	k=2;                  % 2nd order butterworth filter
	fnq=1/(2*TR);
	flp = .04;            % lowpass frequency of filter
	fhi = .07;            % highpass frequency of the filter
	Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
	[bfilt2,afilt2]=butter(k,Wn);   % construct the filter

	[N,Tmax,subjects] = size(time_series); % Some parameters here extracted here.


	% This step is to extract the frequency for each node, first we get the power spectrum for each subject 
	insub = 1;
	for nsub=1:subjects	    
		% Grab the signal data
	    signaldata=squeeze(time_series(:,:,nsub));	    

	    % Work out the time in s (from the volumes)
	    Ts = Tmax*TR;
	    % Determine the frequency vector from the fft

	    freq = (0:Tmax/2-1)/Ts;
	    nfreqs=length(freq);
	    
	    % Here determine the power spectrum for each node
	    for seed=1:N
	        x=detrend(demean(signaldata(seed,:)));
	        ts =zscore(filtfilt(bfilt2,afilt2,x));
	        pw = abs(fft(ts));
	        PowSpect(:,seed,insub) = pw(1:floor(Tmax/2)).^2/(Tmax/TR);
	    end
	    insub=insub+1;
	end


	% This here is to look at the Power distribution, and to average the values for each subject:
	Power_Areas=mean(PowSpect,3);
	for seed=1:N'
		% This gaussian filter ithink is to smooth out the actual fequencies, will have to check
	    Power_Areas(:,seed)=gaussfilt(freq,Power_Areas(:,seed)',0.01);
	end

	% Curiously, is this the same for each area? it would be interesting to see how much this actually makes a difference.
	% For data that has higher resolution it probably matters a little more. 

	[maxpowdata,index]=max(Power_Areas);
	f_diff = freq(index);