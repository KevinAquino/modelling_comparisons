function phfcd = phase_fcd(xs,TR)
    k=2;                  % 2nd order butterworth filter
    fnq=1/(2*TR);
    flp = .04;            % lowpass frequency of filter
    fhi = .07;            % highpass frequency of the filter
    Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
    [bfilt2,afilt2]=butter(k,Wn);   % construct the filter
    [nn,N] = size(xs);
    Tmax = nn;

    upperTriangle = find(tril(ones(N),-1));    
    % N = size(xs,2);


    BOLD=xs';
    % nn = size(xs,1)
    Phase_BOLD=zeros(N,nn);
    signal_filt=zeros(N,nn);
    for seed=1:N
        BOLD(seed,:)=demean(detrend(BOLD(seed,:)));
        signal_filt(seed,:) =filtfilt(bfilt2,afilt2,BOLD(seed,:));
        Phase_BOLD(seed,:) = angle(hilbert(signal_filt(seed,:)));
    end
    T=10:Tmax-10;
    for t=T
        for i=1:N
            for j=1:i-1
                patt(i,j)=cos(adif(Phase_BOLD(i,t),Phase_BOLD(j,t)));
            end
        end        
        pattern(t-9,:)=patt(upperTriangle);
    end
    
    kk=1;
    npattmax=size(pattern,1);
    for t=1:npattmax-2
        p1=mean(pattern(t:t+2,:));
        for t2=t+1:npattmax-2
            p2=mean(pattern(t2:t2+2,:));
            phfcd(kk)=dot(p1,p2)/norm(p1)/norm(p2);
            kk=kk+1;
        end
    end
    