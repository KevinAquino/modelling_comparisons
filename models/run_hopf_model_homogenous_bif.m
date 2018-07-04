function [ts_simulated_all,grandFCcorr,bifpar] = run_hopf_model_homogenous_bif(sc_matrix,time_series,G,folder)

	% Set up the structure for the Hopf model
	for sub=1:size(time_series,3),
		subject_ts{sub} = time_series(:,:,sub);
		corrFC(:,:,sub) = corr(squeeze(time_series(:,:,sub)).');		
	end

	meancorrFC = mean(corrFC,3);
	C = sc_matrix;
	wG = G;	
	TR = 2;
	Tmax = 147;%Check whats your tmax.
	subjects = size(time_series,3);
	N = size(sc_matrix,1);

	% ====== Calculation of the peak freq. for each node (should be made into a function if it is used multiple times) ==============

	% Filtering parameters (to discover the peak of the power spectrum)
	% This sets up the parameters to perform a band-pass filter of the data.

	k=2;                  % 2nd order butterworth filter
	fnq=1/(2*TR);
	flp = .04;            % lowpass frequency of filter
	fhi = .07;            % highpass frequency of the filter
	Wn=[flp/fnq fhi/fnq]; % butterworth bandpass non-dimensional frequency
	[bfilt2,afilt2]=butter(k,Wn);   % construct the filter


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


	
	% ============================== Parameters for the Hopf model ==============================
	% 
	% Setting the global Hopf parameter
	a=-0.01*ones(N,2);	
	% With this, set up the omega vector:
	omega = repmat(2*pi*f_diff',1,2); omega(:,1) = -omega(:,1);
	% Setting up the dt, and sig - strength of noise
	dt = 0.1;
	sig = 0.04;
	
	% Now here we can run the Global model for each global coupling vector
	upperTriangle = find(triu(ones(size(sc_matrix)),1));

	for coupling_index = 1:length(G),
		% Set up the weighted coupling matrix:
		wC = G(coupling_index)*C;
		xs = solve_hopf_ode(omega,a,wC,dt,Tmax,TR,sig);
		% Now after this, we can look at the fitting corr matrix.
		ts_simulated_all(:,:,coupling_index) = xs;
		corrFC_sim = corr(xs);
		grandFCcorr(coupling_index) = corr(corrFC_sim(upperTriangle),meancorrFC(upperTriangle));

	end

	% This is default setting the bifurcation parameter as a;
	bifpar = a;


end