
%Function to simulate a FCM based on the Hopf bifurcation. The time series
%needs to have a N x t dimension within an array {1xsamplesize}. wG is the
%coupling strength vector or scalar, ldata is the number of subjects in the
%sample. opt == 1 performs optimization of the bifurcation parameter "a",
%any other value fixes a=0.001. Outputs are the simulated and empirical FC
%as well as the fitting, metastability and ksd across wG. Also outputs the
%empirical (PhasesD) and simulated phases in each node across time and the
%bifurcation working parameter (bifpar) of each node.
%Victor Saenger, 2015. Adapated from Gustavo Deco.
%
%changes by Jens Schwarzbach
%pertaining to the optimization of a
%1. moved the initialization of a outside the loop over different wG in order
%to inherit the so far best a values to the next weight
%2. moved the initialization of minm outside the iter loop (because it was
%overwritten, which prevented finding the optimal a)
%3. moved the updating of a to the end of the iter loop to prevent
%selecting the a values for optimIter + 1 rather than optimIter
%4. number of iterations for a-optimization is now configurable
%5. update strength for a-optimization is now configurable
%6. used wide bandpass filters for computing relative power in frequencies
%   of interest for empirical and for simulated data
%general
%- optimization for speed
%- removed unnecessary loop over a scalar containing the current weight in the main body
%- return arguments FC_emp, Phases now also have the added dimension idx_g
%  to get one matrix per weight
%- added return argument trackminm1 for diagnostic purposes
%- added input argument simulID (used for saving diagnostic info. e.g.
%  movies, see below)
%- stores simulation movies
%- model fit parallelized (1 process per subject)
%- TR is configurable

%Changes by Victor Saenger
%New optimization method which linearizes the a parameter across G.
%Line 276 added the new funciton aLin
%Cuts execution time by more than a half.


function [FC_simul, FC_emp, fitting, meta, ksP, PhasesD, Phases, bifpar] = hbif_NewOpt_GSR(C,tseries,Tmax,wG,ldata,Cfg)

if ~isfield(Cfg, 'simulID'), Cfg.simulID = 'Unknown'; else end;

if ~isfield(Cfg, 'TRsec'), Cfg.TRsec = 2; else end;

if ~isfield(Cfg, 'opt_a'), Cfg.opt_a = []; else end;
if ~isfield(Cfg.opt_a, 'nIters'), Cfg.opt_a.nIters = 100; else end;
if ~isfield(Cfg.opt_a, 'updateStrength'), Cfg.opt_a.updateStrength = 0.1; else end
if ~isfield(Cfg.opt_a, 'abortCrit'), Cfg.opt_a.abortCrit = 0.1; else end

if ~isfield(Cfg, 'plots'), Cfg.plots = []; else end;
if ~isfield(Cfg.plots, 'showOptimization'), Cfg.plots.showOptimization = 0; else end;
if ~isfield(Cfg.plots, 'makeVideo'), Cfg.plots.makeVideo = 0; else end;


rng('shuffle');
nNodes = length(C);
nSubs = ldata; %JVS: I would really like to replace this later and do subject selection outside this function
si = 1:ldata; %same here; I just keep this for compatibility reasons such that Victor can run his datasets unchanged
nWeights = numel(wG);
fprintf(1, 'Fitting models for %d subjects and %d different weights\n', nSubs, nWeights);

FC_simul = zeros(nNodes, nNodes, nWeights);
fitting = zeros(1, nWeights);
meta = zeros(1, nWeights);
ksP = zeros(1, nWeights);
Phases = zeros(nNodes, Tmax, nSubs, nWeights);
bifpar = zeros(nWeights, nNodes);

%--------------------------------------------------------------------------
%CALCULATE FUNCTIONAL CONNECTIVITY MATRIX
%--------------------------------------------------------------------------
r = zeros(nNodes, nNodes, nSubs);
ts = zeros(nNodes, Tmax, nSubs);
% keyboard

for i = 1:nSubs;
    ts(:,:,i) = tseries{si(i)};
    r(:,:,i) = corrcoef(ts(:,:,i)');
end

FC_emp=mean(r,3);
C=C/max(max(C))*0.2;%


%--------------------------------------------------------------------------
%COMPUTE POWER SPECTRA FOR
%NARROWLY FILTERED DATA WITH LOW BANDPASS (0.04 to 0.07 Hz)
%WIDELY FILTERED DATA (0.04 Hz to justBelowNyquistFrequency)
%[justBelowNyquistFrequency depends on TR,
%for a TR of 2s this is 0.249 Hz]
%--------------------------------------------------------------------------
TT=Tmax;
Ts = TT*Cfg.TRsec;
freq = (0:TT/2-1)/Ts;
[~, idxMinFreq] = min(abs(freq-0.04));
[~, idxMaxFreq] = min(abs(freq-0.07));
nFreqs = length(freq);

delt = 2;                                   % sampling interval
fnq = 1/(2*delt);                           % Nyquist frequency
k = 2;                                      % 2nd order butterworth filter

%WIDE BANDPASS
flp = .04;                                  % lowpass frequency of filter
fhi = fnq-0.001;%.249;                      % highpass needs to be limited by Nyquist frequency, which in turn depends on TR
Wn = [flp/fnq fhi/fnq];                     % butterworth bandpass non-dimensional frequency
[bfilt_wide, afilt_wide] = butter(k,Wn);    % construct the filter

%NARROW LOW BANDPASS
flp = .04;                                  % lowpass frequency of filter
fhi = .07;                                  % highpass
Wn=[flp/fnq fhi/fnq];                       % butterworth bandpass non-dimensional frequency
[bfilt_narrow,afilt_narrow] = butter(k,Wn); % construct the filter


PowSpect_filt_narrow = zeros(nFreqs, nNodes, nSubs);
PowSpect_filt_wide = zeros(nFreqs, nNodes, nSubs);
for seed=1:nNodes
    
    for idxSub=1:nSubs
        signaldata = tseries{si(idxSub)};
        x=detrend(demean(signaldata(seed,:)));
        
        ts_filt_narrow =zscore(filtfilt(bfilt_narrow,afilt_narrow,x));
        pw_filt_narrow = abs(fft(ts_filt_narrow));
        PowSpect_filt_narrow(:,seed,idxSub) = pw_filt_narrow(1:floor(TT/2)).^2/(TT/2);
        
        ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));
        pw_filt_wide = abs(fft(ts_filt_wide));
        PowSpect_filt_wide(:,seed,idxSub) = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
    end
    
end


Power_Areas_filt_narrow_unsmoothed = mean(PowSpect_filt_narrow,3);
Power_Areas_filt_wide_unsmoothed = mean(PowSpect_filt_wide,3);
Power_Areas_filt_narrow_smoothed = zeros(nFreqs, nNodes);
Power_Areas_filt_wide_smoothed = zeros(nFreqs, nNodes);
vsig = zeros(1, nNodes);
for seed=1:nNodes
    Power_Areas_filt_narrow_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_narrow_unsmoothed(:,seed)',0.01);
    Power_Areas_filt_wide_smoothed(:,seed)=gaussfilt(freq,Power_Areas_filt_wide_unsmoothed(:,seed)',0.01);
    %relative power in frequencies of interest (.04 - .07 Hz) with respect
    %to entire power of bandpass-filtered data (.04 - just_below_nyquist)
    vsig(seed) =...
        sum(Power_Areas_filt_wide_smoothed(idxMinFreq:idxMaxFreq,seed))/sum(Power_Areas_filt_wide_smoothed(:,seed));
end
vmax=max(vsig); %consider computing this later where needed
vmin=min(vsig);%consider computing this later where needed

%a-minimization seems to only work if we use the indices for frequency of
%maximal power from the narrowband-smoothed data
[~, idxFreqOfMaxPwr]=max(Power_Areas_filt_narrow_smoothed);
f_diff = freq(idxFreqOfMaxPwr);

%FOR EACH AREA AND TIMEPOINT COMPUTE THE INSTANTANEOUS PHASE IN THE RANGE
%OF .04 TO .07 Hz
PhasesD = zeros(nNodes, Tmax, nSubs);
for idxSub = 1:nSubs
    signaldata=tseries{si(idxSub)};
    for seed=1:nNodes
        x = demean(detrend(signaldata(seed,:)));
        xFilt = filtfilt(bfilt_narrow,afilt_narrow,x);    % zero phase filter the data
        Xanalytic = hilbert(demean(xFilt));
        PhasesD(seed,:,idxSub) = angle(Xanalytic);
    end 
end

%f_diff  previously computed frequency with maximal power (of narrowly filtered data) by area
omega = repmat(2*pi*f_diff',1,2); %angular velocity
omega(:,1) = -omega(:,1);


%FROM HERE ON SIMULATIONS AND FITTING
dt = 0.1;
sig = 0.04; %was 0.04
dsig = sqrt(dt)*sig; % to avoid sqrt(dt) at each time step

a = repmat(-0.05*ones(nNodes,1),1,2);
wC = Cfg.opt_a.gref*C;
sumC = repmat(sum(wC,2),1,2);
trackminm1 = zeros(Cfg.opt_a.nIters, nWeights); %for tracking the minimization (good for debugging)
xs = zeros(3000/2,nNodes);

  %Optimize for G ref;
      minm=100;
        bestIter = 1; %JVS for tracking purposes
        for iter = 1:Cfg.opt_a.nIters
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            for t=1:dt:1000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
            end
            
            for t=1:dt:3000%3000 seconds worth of simulated data; why dont we limit it to Tmax*Cfg.TRsec?
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
                
                if mod(t,2)==0
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            
            xs = RegressNoiseSignal(xs.',mean(xs.')).';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            vsigs = zeros(1, nNodes);
            for seed=1:nNodes
                
                x=detrend(demean(xs(1:nn,seed)'));%was x
                ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));%zscore before!
                
                TT=length(x);
                Ts = TT*Cfg.TRsec;
                freq = (0:TT/2-1)/Ts;
                [~, idxMinFreqS]=min(abs(freq-0.04));
                [~, idxMaxFreqS]=min(abs(freq-0.07));
                
                pw_filt_wide = abs(fft(ts_filt_wide));
                Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
                Pow=gaussfilt(freq,Pow1,0.01);
                
                vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
                
            end
            
            vsmin=min(vsigs);%vsig y vsigs max min points not equal
            vsmax=max(vsigs);
            bb=(vmax-vmin)/(vsmax-vsmin);
            aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
            vsigs=aa+bb*vsigs;
            minm1=max(abs(vsig-vsigs)./vsig);
            trackminm1(iter, Cfg.opt_a.gref) = minm1;%JVS tracking of minm1
            
            if minm1<minm
                minm=minm1;
                a1=a;
                bestIter = iter; %JVS: tracking
                best_vsigs = vsigs; %JVS: tracking
            end
            
            %--------------------------------------------------------------
            %FEEDBACK
            %--------------------------------------------------------------
            if Cfg.plots.showOptimization
                showOptimPlot(h_track_opt, idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
            end
            fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
            %--------------------------------------------------------------
            
            %CRITERION REACHED?
            if minm<Cfg.opt_a.abortCrit %default is 0.1
                break;
            end
            
            %UPDATE a VALUES FOR NEXT ITER
            if ~any(isnan(vsigs))
                a(:,1)=a(:,1)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
                a(:,2)=a(:,2)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
            else
                %JVS: this should not happen according to Gustavo, but it
                %can when a values run crazy
                warning('There are NaNs in the power spectra. Probably a-values too strong.');
                a = a1;
            end
        end
        a=a1;
  %%%%%
  %Linearize for all Gs
   aLinear=aLin(a(:,1)',Cfg.opt_a.gref,wG);
   
%SIMULATE FOR EACH G USING LINEARIZED A VALUES:
for idx_g = 1:nWeights
    we = wG(idx_g);
    wC = we*C;
    xs = zeros(3000/2,nNodes);
    sumC = repmat(sum(wC,2),1,2);
    
    fprintf(1, '-----------------------------------------\n');
    fprintf(1, 'G(%d/%d) = %5.3f\n', idx_g, numel(wG), we);
    fprintf(1, '-----------------------------------------\n');
    
    %Choose current idx_g from aLinear and optimize for 30 iters only
    av=aLinear(idx_g,:);
    a=repmat(av',1,2);

        minm=100;
        bestIter = 1; %JVS for tracking purposes
        for iter = 1:100
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
            nn=0;
            for t=1:dt:1000
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
            end
            
            for t=1:dt:3000%3000 seconds worth of simulated data; why dont we limit it to Tmax*Cfg.TRsec?
                suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
                zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
                z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
                
                if mod(t,2)==0
                    nn=nn+1;
                    xs(nn,:)=z(:,1)';
                end
            end
            
            xs = RegressNoiseSignal(xs.',mean(xs.')).';
            
            %%%%%%%%%%%%%%%%%%%%%%%%%
            vsigs = zeros(1, nNodes);
            for seed=1:nNodes
                
                x=detrend(demean(xs(1:nn,seed)'));%was x
                ts_filt_wide =zscore(filtfilt(bfilt_wide,afilt_wide,x));%zscore before!
                
                
                TT=length(x);
                Ts = TT*Cfg.TRsec;
                freq = (0:TT/2-1)/Ts;
                [~, idxMinFreqS]=min(abs(freq-0.04));
                [~, idxMaxFreqS]=min(abs(freq-0.07));
                
                pw_filt_wide = abs(fft(ts_filt_wide));
                Pow1 = pw_filt_wide(1:floor(TT/2)).^2/(TT/2);
                Pow=gaussfilt(freq,Pow1,0.01);
                
                vsigs(seed)=sum(Pow(idxMinFreqS:idxMaxFreqS))/sum(Pow);
                
            end
            
            vsmin=min(vsigs);%vsig y vsigs max min points not equal
            vsmax=max(vsigs);
            bb=(vmax-vmin)/(vsmax-vsmin);
            aa=vmin-bb*vsmin;  %% adaptation of local bif parameters
            vsigs=aa+bb*vsigs;
            minm1=max(abs(vsig-vsigs)./vsig);
            trackminm1(iter, idx_g) = minm1;%JVS tracking of minm1
            
            if minm1<minm
                minm=minm1;
                a1=a;
                bestIter = iter; %JVS: tracking
                best_vsigs = vsigs; %JVS: tracking
            end
            
            %--------------------------------------------------------------
            %FEEDBACK
            %--------------------------------------------------------------
            if Cfg.plots.showOptimization
                showOptimPlot(h_track_opt, idx_g, we, iter, a, a1, vsig, vsigs, bestIter, best_vsigs, trackminm1, Cfg)
            end
            fprintf(1, 'iter: %03d, minm1: %5.3f\n', iter, minm1);
            %--------------------------------------------------------------
            
            %CRITERION REACHED?
            if minm<Cfg.opt_a.abortCrit %default is 0.1
                break;
            end
            
            %UPDATE a VALUES FOR NEXT ITER
            if ~any(isnan(vsigs))
                a(:,1)=a(:,1)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
                a(:,2)=a(:,2)+Cfg.opt_a.updateStrength*(1-vsigs./vsig)';
            else
                %JVS: this should not happen according to Gustavo, but it
                %can when a values run crazy
                warning('There are NaNs in the power spectra. Probably a-values too strong.');
                a = a1;
            end
        end
    
    a=a1;%use those avalues that have been found to be optimal
    bifpar(idx_g,:)=a(:,1)';
    %%%%%%%%%%%%%%
    %%% Final simul
    fprintf(1, 'SIMULATING OPTIMIZED MODEL.\n');
    
    xs=zeros(Tmax*nSubs,nNodes);
    z = 0.1*ones(nNodes,2); % --> x = z(:,1), y = z(:,2)
    nn=0;
    for t=1:dt:3000 %JVS is it really necessary to swing in for 3000secs?
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
    end
    
    for t=1:dt:Tmax*Cfg.TRsec*nSubs %JVS: was 15000, now faster
        suma = wC*z - sumC.*z; % sum(Cij*xi) - sum(Cij)*xj
        zz = z(:,end:-1:1); % flipped z, because (x.*x + y.*y)
        z = z + dt*(a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma) + dsig*randn(nNodes,2);
        
        if mod(t,Cfg.TRsec)==0
            nn=nn+1;
            xs(nn,:)=z(:,1)';
        end
    end
    
    fprintf(1, 'COMPUTING MODEL FIT.\n');
    
    xs = RegressNoiseSignal(xs.',mean(xs.')).';
    FC_simul(:, :, idx_g) = corrcoef(xs(1:nn,:)); %Now one FC_simul per G
%     keyboard
    cc=corrcoef(squareform(tril(FC_emp,-1)),squareform(tril(FC_simul(:, :, idx_g),-1)));%atanh(FC...
    fitting(idx_g)=cc(2);
    
    %%%%%%%%%%%%%%%%%%
    
    metastability22 = zeros(1, nSubs);
    for idxSub=1:nSubs
        tini=(idxSub-1)*Tmax;
        for seed=1:nNodes
            ts_simul = detrend(demean(xs(tini+1:tini+Tmax,seed)'));
            ts_simul_filt_narrow = filtfilt(bfilt_narrow,afilt_narrow,ts_simul);
            Xanalytic = hilbert(ts_simul_filt_narrow);
            Phases(seed,:,idxSub, idx_g) = angle(Xanalytic);
        end
        
        T=10:Tmax-10;
        sync = zeros(1, numel(T));
        for t=T
            ku=sum(complex(cos(Phases(:,t,idxSub, idx_g)),sin(Phases(:,t,idxSub, idx_g))))/nNodes;
            sync(t-9, idx_g)=abs(ku);
        end
        
        metastability22(idxSub)=std(sync(:, idx_g));
        
    end
    
    meta(idx_g)=mean(metastability22);
    
    %can replace with parfor if PARALLEL PROCESSING TOOL is installed:
    for i = 1:nSubs
        pcD(:,i)=patternCons(PhasesD(:,:,i),nNodes,Tmax);
        pcS(:,i)=patternCons(Phases(:,:,i, idx_g),nNodes,Tmax); 
    end
    
    [~, ~, ksP(idx_g)]=kstest2(pcS(:),pcD(:));
      
    fprintf(1, 'DONE.\n\n');
end


