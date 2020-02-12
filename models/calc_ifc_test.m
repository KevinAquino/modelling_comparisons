function iFC = calc_ifc_test(time_series,TR)
    % keyboard
	delt = TR;            % sampling interval
	k=2;                  % 2nd order butterworth filter
	fnq=1/(2*delt);
	flp = .04;           % lowpass frequency of filter
	fhi = .07;           % highpass
	Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
	[bfilt2,afilt2]=butter(k,Wn);   % construct the filter
	[N,Tmax] = size(time_series);

	signaldata=squeeze(time_series);  %%seed, time, sub    
    Phase_BOLD_data=zeros(N,Tmax);
    % keyboard
    for seed=1:N
        signaldata(seed,:)=signaldata(seed,:)-mean(signaldata(seed,:));
        signal_filt_data =filtfilt(bfilt2,afilt2,signaldata(seed,:));
        Phase_BOLD_data(seed,:) = angle(hilbert(signal_filt_data));
    end
    iFC=zeros(Tmax,N,N);
    for t=1:Tmax
        for n=1:N
            for p=1:N
                iFC(t,n,p)=cos(Phase_BOLD_data(n,t)-Phase_BOLD_data(p,t));
            end
        end
    end

end