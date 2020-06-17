% A seperate function to calculate the power integral, to make the calculations robust

function vsig = calculatePowerRatio(time_series,total_time,TR)

    % Nodes = number of nodes in the network
    % Instances = number of instances of the time series, for empirical data it is subjects for simulation it is runs
    [~,Nodes,Instances] = size(time_series);

    % Here are the parameters for the calculation of the power ratio
    low_cutoff = 0.04;                          % 0.04 Hz low cutoff 
    high_cutoff = 0.07;                         % 0.07 Hz high cutoff
    delt = TR;                                  % sampling interval
    fnq = 1/(2*delt);                           % Nyquist frequency
    nyq_cutoff = fnq-0.001;%.249;           % highpass needs to be limited by Nyquist frequency, which in turn depends on TR
    k = 2;                                      % 2nd order butterworth filter


    % Calculation of variables needed for the filtering
    TT=total_time;
    Ts = TT*TR;
    freq = (0:TT/2-1)/Ts;
    [~, idxMinFreq] = min(abs(freq-low_cutoff));
    [~, idxMaxFreq] = min(abs(freq-high_cutoff));
    nFreqs = length(freq);

    % WIDE BANDPASS filtering conversion 
    Wn = [low_cutoff/fnq nyq_cutoff/fnq];        % butterworth bandpass non-dimensional frequency
    [bfilt_wide, afilt_wide] = butter(k,Wn);    % construct the filter


    PowSpect_filt_wide = zeros(nFreqs, Nodes, Instances);
    for seed=1:Nodes
        
        for t_instance=1:Instances

            signaldata = time_series(:,seed,t_instance);
            x=zscore(detrend(demean(signaldata)));
            
            ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
            pw_filt_wide = abs(fft(ts_filt_wide));
            PowSpect_filt_wide(:,seed,t_instance) = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
        end
        
    end

    
    Power_Areas_filt_wide_unsmoothed = mean(PowSpect_filt_wide,3);    
    Power_Areas_filt_wide_smoothed = zeros(nFreqs, Nodes);
    vsig = zeros(1, Nodes);
    for seed=1:Nodes        
        Power_Areas_filt_wide_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_wide_unsmoothed(:,seed)',0.01);
        %relative power in frequencies of interest (.04 - .07 Hz) with respect
        %to entire power of bandpass-filtered data (.04 - just_below_nyquist)
        vsig(seed) =...
            sum(Power_Areas_filt_wide_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_wide_smoothed(:,seed));
    end
